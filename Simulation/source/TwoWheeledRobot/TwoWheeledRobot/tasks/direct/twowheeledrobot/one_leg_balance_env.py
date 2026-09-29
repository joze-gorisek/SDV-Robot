"""One-leg balance environment.

Left side of robot rests on the ground wheel; right (top) joints are locked at
full extension.  The policy controls:
  - left/right wheel current  (pitch + velocity regulation, same as PureNNBalanceEnv)
  - front_left / back_left joint angle  (roll regulation)

Observation (6D):
  roll, roll_rate, pitch, pitch_rate,
  x_rel (left-wheel odometry since reset), velocity (left wheel only)
  (action/current history is carried by the GRU hidden state, not observed)

Action (3D):
  wheel_L in [-1, 1]  (passed to CurrentActionProcessor; right wheel is lifted → always zero)
  joint_FL, joint_BL in [-1, 1]  (maps to center ± range in extension degrees)

Termination:
  - pitch fall (inherited: |pitch| > 40° for 5 consecutive steps)
  - roll window (|roll| < 20° OR |roll| > 72° for 5 consecutive steps)
  - physics broken / invalid state (inherited)
  - timeout (inherited)
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import torch

from .one_leg_balance_env_cfg import OneLegBalanceEnvCfg
from .pure_nn_components import pitch_from_projected_gravity, roll_from_projected_gravity
from .pure_nn_balance_env import PureNNBalanceEnv
from .residual_lqr_env import R_WHEEL

# CyberGear index within _cg_ids: [front_left, front_right, back_left, back_right]
_IDX_FL = 0
_IDX_FR = 1
_IDX_BL = 2
_IDX_BR = 3


class OneLegBalanceEnv(PureNNBalanceEnv):
    cfg: OneLegBalanceEnvCfg

    def __init__(self, cfg: OneLegBalanceEnvCfg, render_mode: str | None = None, **kwargs):
        # Right-wheel-touch termination needs a real PhysX contact sensor, not
        # a joint-extension proxy. Force it on for this task specifically
        # (other tasks keep enable_wheel_contacts=False by default) -- must be
        # set BEFORE super().__init__(), which calls _setup_scene() and
        # instantiates the ContactSensor from these cfg fields. Same pattern
        # already proven in scripts/lqr_control_one_leg_jump.py.
        cfg.enable_wheel_contacts = True
        cfg.robot_cfg.spawn.activate_contact_sensors = True
        super().__init__(cfg, render_mode, **kwargs)

        # Resolve the RIGHT wheel's contact-sensor body index (DDSM115_Simplified_01,
        # per the USD joint dump: DDSM115_Desni's body0 -- the right wheel).
        # find_bodies here is on the ContactSensor (bodies it was configured to
        # track via wheel_contacts' prim_path regex), not the articulation.
        right_ids, right_names = self.wheel_contacts.find_bodies("DDSM115_Simplified_01")
        if len(right_ids) != 1:
            raise RuntimeError(f"Expected exactly one right-wheel contact body, found {right_names}")
        self._right_wheel_contact_id = torch.tensor(right_ids, device=self.device, dtype=torch.long)

        # Right-wheel-touch termination counter (same counter-based pattern as
        # _roll_fall_counter / _wrong_direction_counter -- a few consecutive
        # steps above the force threshold, not an instant trigger, so a brief
        # settling bounce right after reset doesn't spuriously end the episode).
        self._right_wheel_touch_counter = torch.zeros(self.num_envs, device=self.device, dtype=torch.long)
        self._right_wheel_touched = torch.zeros(self.num_envs, device=self.device, dtype=torch.bool)

        # Pre-compute constants so they are not recomputed every step.
        self._bottom_joint_center_rad = math.radians(cfg.bottom_joint_center_deg)
        self._bottom_joint_range_rad = math.radians(cfg.bottom_joint_range_deg)
        self._top_joint_center_rad = math.radians(cfg.top_joint_center_deg)
        self._top_joint_range_rad = math.radians(cfg.top_joint_range_deg)
        self._reset_joint_ext_rad = math.radians(cfg.reset_joint_ext_deg)

        # Roll termination counter (separate from parent's pitch _fall_counter).
        self._roll_fall_counter = torch.zeros(self.num_envs, device=self.device, dtype=torch.long)

        # Joint action history for smoothness penalty. All four CyberGear joints
        # are now policy-controlled (5D action).
        # Shape: [num_envs, 4] — [FL, BL, FR, BR] actions in [-1, 1].
        self._prev_joint_actions = torch.zeros(self.num_envs, 4, device=self.device)
        self._cur_joint_actions = torch.zeros(self.num_envs, 4, device=self.device)

        # Station keeping: position since episode start, integrated from
        # left-wheel odometry (mirrors the STM32 DDSM115 encoder integration).
        self._x_rel = torch.zeros(self.num_envs, device=self.device)

        # Two-phase reward latch: True once roll has first crossed
        # phase2_roll_deg this episode. Gates both the reward phase (reach vs
        # balance) and the "fell back to two-leg stance" termination (which must
        # NOT fire at the flat spawn, where roll starts at ~0 < roll_low).
        self._has_reached = torch.zeros(self.num_envs, device=self.device, dtype=torch.bool)

        # Phase 1 wrong-direction termination: counter-based (like
        # _roll_fall_counter), flag consumed by _get_dones.
        self._wrong_direction_counter = torch.zeros(self.num_envs, device=self.device, dtype=torch.long)
        self._wrong_direction = torch.zeros(self.num_envs, device=self.device, dtype=torch.bool)

        # Cache ext_sign scalars for bottom and top joints (avoid repeated indexing).
        ext = self._cg_ext_sign  # shape [1, 4]
        self._ext_sign_fl = float(ext[0, _IDX_FL])   # +1
        self._ext_sign_fr = float(ext[0, _IDX_FR])   # -1
        self._ext_sign_bl = float(ext[0, _IDX_BL])   # -1
        self._ext_sign_br = float(ext[0, _IDX_BR])   # +1

        # Precomputed raw joint angles for the one-leg-stance spawn: left legs
        # (FL/BL) grounded at spawn_left_ext_deg, right legs (FR/BR) retracted
        # to spawn_right_ext_deg so the right wheel lifts off.
        self._spawn_left_ext_rad = math.radians(cfg.spawn_left_ext_deg)
        self._spawn_right_ext_rad = math.radians(cfg.spawn_right_ext_deg)
        self._spawn_raw_fl = self._ext_sign_fl * self._spawn_left_ext_rad
        self._spawn_raw_fr = self._ext_sign_fr * self._spawn_right_ext_rad
        self._spawn_raw_bl = self._ext_sign_bl * self._spawn_left_ext_rad
        self._spawn_raw_br = self._ext_sign_br * self._spawn_right_ext_rad

        # ── Passive parallelogram-bearing spawn angles (loop-closure FK) ──────
        # front_left/back_left/front_right/back_right are DRIVEN (position-
        # controlled) joints, but each sits in a closed parallelogram loop with
        # passive "bearing_joints" (stiffness=0, robot_cfg.py) that have NO
        # authored default pose anywhere in the USD -- confirmed by directly
        # inspecting World0.usd: every joint's state:angular:physics:position
        # is unauthored, so Isaac Lab's default_joint_pos falls back to 0 rad
        # for them. The parent reset (standup_env) writes default_joint_pos for
        # every non-CG joint, so these bearings get hard-reset to 0 rad every
        # episode regardless of what the driven CG joints are commanded to.
        # 0 rad is NOT the geometrically-consistent value for this loop at
        # EITHER the flat-spawn (45/-45) or one-leg-stance (45/-45 left,
        # 0/0 right) configuration -- verified by solving the loop-closure
        # constraint numerically from World0.usd's actual joint transforms
        # (Newton least-squares on the left/right sub-loops, residual ~1e-8
        # m/rad, i.e. essentially exact):
        #   left  loop (front_left=+45deg, back_left=-45deg, UNCHANGED from the
        #   original flat spawn -- this bug predates the one-leg-stance spawn,
        #   it was just invisible with both wheels grounded):
        #     Revolute_9   = -1.0004237149366428 rad (-57.32 deg)
        #     tochangeleft =  0.4300510778806665 rad ( 24.64 deg)  -- NOT a
        #       writable joint (see below), value kept here only for reference.
        #     Revolute_10  =  1.0004237046476674 rad ( 57.32 deg)
        #   right loop (front_right=0deg, back_right=0deg -- the new retracted
        #   spawn): solves to ~0 rad for all three, i.e. the 0-rad default
        #   already happens to be correct here, no override needed.
        # Without this, the robot spawned with its left-side wheel mount
        # folded into an inconsistent, overstressed configuration ("wheel
        # above the joint instead of below").
        #
        # tochangeleft/tochangeright are NOT queryable via find_joints() --
        # confirmed at runtime (ValueError: "tochangeleft: []", only 10 joints
        # exposed: front/back_left/right, Revolute_1/2/9/10, DDSM115_*).
        # PhysX collapsed them into an internal loop-closure (maximal-
        # coordinate) constraint rather than exposing them as regular
        # reduced-coordinate DOFs -- each parallelogram loop has one more
        # joint than PhysX's articulation tree can represent as independent
        # DOFs, so the "closing" joint is handled implicitly. Only Revolute_9
        # and Revolute_10 are set explicitly; PhysX resolves the implicit
        # tochangeleft constraint on its own once front_left/back_left/
        # Revolute_9/Revolute_10 are all consistent, since that fully
        # determines the loop.
        self._left_bearing_ids, _ = self.robot.find_joints("Revolute_9")
        self._left_bearing3_ids, _ = self.robot.find_joints("Revolute_10")
        self._left_bearing_targets_rad = (
            -1.0004237149366428,  # Revolute_9
            1.0004237046476674,   # Revolute_10
        )

        # ── Iteration-based spawn-roll curriculum ──────────────────────────────
        # roll = spawn_roll_center_deg ± band(iteration). The iteration is
        # derived from self.common_step_counter // curriculum_num_steps_per_env,
        # plus an offset set by train.py on resume (set_curriculum_iteration_offset).
        self._curriculum_iter_offset = 0
        self._curriculum_last_stage = -1
        self._spawn_roll_center_rad = math.radians(cfg.spawn_roll_center_deg)
        stages = sorted(cfg.spawn_roll_curriculum, key=lambda s: s[0])
        self._curriculum_starts = [int(s[0]) for s in stages]
        self._curriculum_bands_rad = [math.radians(float(s[1])) for s in stages]

        # Observation normalization tensor matching cfg.observation_scale (10D).
        self._obs_scale = torch.tensor(
            cfg.observation_scale, device=self.device, dtype=torch.float32
        ).view(1, -1)
        if self._obs_scale.shape[1] != cfg.observation_space:
            raise ValueError("observation_scale length must match observation_space (6)")

        # Disable per-step disturbance forces. The inherited PureNNBalanceEnv
        # calls set_external_force_and_torque every policy step (a deprecated
        # GPU sync point). One-leg training does not use disturbances, so clear
        # body_ids to skip the call entirely and recover a significant fraction
        # of the collection time.
        self._body_ids = None

        print(
            "[OneLegBalanceEnv] one-leg balance controller active, "
            f"dt={self.step_dt:.3f}s  obs={cfg.observation_space}  act={cfg.action_space}  "
            f"roll_window=[{cfg.roll_low_threshold_deg:.0f}°, {cfg.roll_high_threshold_deg:.0f}°]  "
            f"bottom_joint={cfg.bottom_joint_center_deg:.0f}°±{cfg.bottom_joint_range_deg:.0f}°  "
            f"top_joint={cfg.top_joint_center_deg:.0f}°±{cfg.top_joint_range_deg:.0f}° (unlocked)  "
            f"spawn=one-leg-stance(roll={cfg.spawn_roll_center_deg:.0f}°±band, pitch={cfg.reset_pitch_deg:.0f}°, "
            f"left_legs={cfg.spawn_left_ext_deg:.0f}°, right_legs={cfg.spawn_right_ext_deg:.0f}°)  "
            f"roll_curriculum={cfg.spawn_roll_curriculum}"
        )

    # ── spawn-roll curriculum ───────────────────────────────────────────────

    def set_curriculum_iteration_offset(self, offset: int) -> None:
        """Align the spawn-roll curriculum with the true training iteration.

        Called by train.py after a resume (the env's step counter restarts at 0
        on a fresh process, so without this the curriculum would fall back to
        stage 0). Fresh runs pass offset=0.
        """
        self._curriculum_iter_offset = int(offset)

    def set_spawn_roll_deg(self, center_deg: float, band_deg: float = 0.0) -> None:
        """Override the spawn-roll angle at runtime (e.g. an eval sweep in
        play.py). center_deg/spawn_roll_curriculum are cached into
        _spawn_roll_center_rad/_curriculum_bands_rad/_curriculum_starts at
        __init__ and never re-read from self.cfg afterward, so mutating
        env_cfg.spawn_roll_center_deg post-construction has NO effect --
        this is the actual runtime hook. Call env.reset() right after this to
        force every env to respawn at the new angle immediately, rather than
        waiting for natural episode termination.
        """
        self._spawn_roll_center_rad = math.radians(center_deg)
        self._curriculum_starts = [0]
        self._curriculum_bands_rad = [math.radians(band_deg)]
        self._curriculum_last_stage = -1
        self.cfg.spawn_roll_center_deg = center_deg  # keep cfg/logging consistent, not read by _reset_idx
        self.cfg.spawn_roll_curriculum = ((0, band_deg),)

    def _current_iteration(self) -> int:
        """Training iteration derived from the total env step count.

        rsl-rl collects curriculum_num_steps_per_env policy steps per env per
        iteration, and every env is stepped together, so one step() call ==
        one policy step for the whole batch. common_step_counter counts those.
        """
        steps_per_iter = max(1, int(self.cfg.curriculum_num_steps_per_env))
        return self._curriculum_iter_offset + int(self.common_step_counter) // steps_per_iter

    def _current_spawn_roll_band_rad(self) -> float:
        """Half-band (rad) of the spawn-roll uniform range for this iteration."""
        it = self._current_iteration()
        band = self._curriculum_bands_rad[0]
        stage = 0
        for i, start in enumerate(self._curriculum_starts):
            if it >= start:
                band = self._curriculum_bands_rad[i]
                stage = i
        # Announce a stage advance once, so training logs show the transition.
        if stage != self._curriculum_last_stage:
            self._curriculum_last_stage = stage
            print(
                f"[OneLegBalanceEnv] spawn-roll curriculum -> stage {stage + 1}"
                f"/{len(self._curriculum_starts)} at iteration ~{it}: "
                f"roll = {self.cfg.spawn_roll_center_deg:.0f}° ± {math.degrees(band):.0f}°",
                flush=True,
            )
        return band

    # ── helpers ───────────────────────────────────────────────────────────────

    def _compute_roll_rate(self) -> torch.Tensor:
        """Approximate roll rate from body-frame angular velocity around Y axis."""
        # Body Y axis is the roll axis (tilt left/right).
        # Sign: positive omega_y → robot tilts right (consistent with roll formula).
        return self.bno080.data.ang_vel_b[:, 1]

    def _compute_velocity(self) -> torch.Tensor:
        """Ground speed from the left (grounded) wheel only.

        The parent's two-wheel average underestimates by ~2× here because the
        right wheel is lifted and idle.
        """
        raw_wheel_vel = self.robot.data.joint_vel[:, self._wheel_ids] * self._wheel_sign
        return raw_wheel_vel[:, 0] * R_WHEEL

    # ── physics step ──────────────────────────────────────────────────────────

    def _pre_physics_step(self, actions: torch.Tensor) -> None:
        # Action layout (5D): [wheel_L, joint_FL, joint_BL, joint_FR, joint_BR]
        # Right wheel is lifted — insert zero so the parent's 4D wheel processor
        # receives [wheel_L, 0] without any changes to PureNNBalanceEnv.
        actions_4d = torch.zeros(self.num_envs, 4, device=self.device)
        actions_4d[:, 0] = actions[:, 0]   # wheel_L
        # actions_4d[:, 1] stays 0          # wheel_R (lifted, no current)

        # Parent handles: wheel current processing, disturbance forces, pitch bookkeeping.
        # It also sets all 4 CG joints to zero — we override all of them below.
        super()._pre_physics_step(actions_4d)

        # ── all four joints: map policy action [-1, 1] → physical extension ───
        # order stored in _cur_joint_actions: [FL, BL, FR, BR]
        self._prev_joint_actions = self._cur_joint_actions.clone()
        joint_actions = actions[:, 1:5].clamp(-1.0, 1.0)
        self._cur_joint_actions = joint_actions

        # Bottom (FL, BL): narrow range around bottom_center. Top (FR, BR): wide
        # range around top_center so the policy can extend to tip and fold to 0°.
        ext_fl = self._bottom_joint_center_rad + joint_actions[:, 0] * self._bottom_joint_range_rad
        ext_bl = self._bottom_joint_center_rad + joint_actions[:, 1] * self._bottom_joint_range_rad
        ext_fr = self._top_joint_center_rad + joint_actions[:, 2] * self._top_joint_range_rad
        ext_br = self._top_joint_center_rad + joint_actions[:, 3] * self._top_joint_range_rad

        raw_fl = self._ext_sign_fl * ext_fl
        raw_bl = self._ext_sign_bl * ext_bl
        raw_fr = self._ext_sign_fr * ext_fr
        raw_br = self._ext_sign_br * ext_br

        raw_fl.clamp_(float(self._cg_joint_lo[0, _IDX_FL]), float(self._cg_joint_hi[0, _IDX_FL]))
        raw_bl.clamp_(float(self._cg_joint_lo[0, _IDX_BL]), float(self._cg_joint_hi[0, _IDX_BL]))
        raw_fr.clamp_(float(self._cg_joint_lo[0, _IDX_FR]), float(self._cg_joint_hi[0, _IDX_FR]))
        raw_br.clamp_(float(self._cg_joint_lo[0, _IDX_BR]), float(self._cg_joint_hi[0, _IDX_BR]))

        # Write in _cg_ids order: [FL, FR, BL, BR].
        cg_targets = torch.stack([raw_fl, raw_fr, raw_bl, raw_br], dim=1)
        self.robot.set_joint_position_target(cg_targets, joint_ids=self._cg_ids)

    # ── observations ──────────────────────────────────────────────────────────

    def _get_observations(self) -> dict:
        self._enforce_cybergear_joint_state_limits()

        grav_b = self.bno080.data.projected_gravity_b
        roll = roll_from_projected_gravity(grav_b)
        roll_rate = self._compute_roll_rate()
        pitch = pitch_from_projected_gravity(grav_b)
        pitch_rate = -self.bno080.data.ang_vel_b[:, 0]

        # Add sensor noise (same style as parent).
        roll = roll + torch.randn_like(roll) * self.cfg.pitch_noise_std
        roll_rate = roll_rate + torch.randn_like(roll_rate) * self.cfg.pitch_rate_noise_std
        pitch = pitch + self._pitch_bias + torch.randn_like(pitch) * self.cfg.pitch_noise_std
        pitch_rate = pitch_rate + torch.randn_like(pitch_rate) * self.cfg.pitch_rate_noise_std

        velocity = self._compute_velocity()

        # 6-D observation: physical state only. Action/current history is left to
        # the GRU hidden state, so prev-current and prev-action are not fed in.
        obs_raw = torch.cat([
            roll.unsqueeze(1),           # [0] roll
            roll_rate.unsqueeze(1),      # [1] roll rate
            pitch.unsqueeze(1),          # [2] pitch
            pitch_rate.unsqueeze(1),     # [3] pitch rate
            self._x_rel.unsqueeze(1),    # [4] drift from start (m)
            velocity.unsqueeze(1),       # [5] ground speed, left wheel (m/s)
        ], dim=1)

        obs = torch.nan_to_num(
            obs_raw / self._obs_scale, nan=0.0, posinf=10.0, neginf=-10.0
        ).clamp(-10.0, 10.0)

        # Apply optional 1-step observation delay (inherited random flag).
        obs_out = torch.where(self._obs_delay_samples.view(-1, 1) > 0, self._obs_delay, obs)
        self._obs_delay = obs.clone()

        return {"policy": obs_out}

    # ── rewards ───────────────────────────────────────────────────────────────

    def _get_rewards(self) -> torch.Tensor:
        grav_b = self.bno080.data.projected_gravity_b
        roll = roll_from_projected_gravity(grav_b)
        roll_rate = self._compute_roll_rate()
        pitch = pitch_from_projected_gravity(grav_b)
        pitch_rate = -self.bno080.data.ang_vel_b[:, 0]

        velocity = self._compute_velocity()
        # Integrate left-wheel odometry (once per control step; rewards run
        # before resets, so freshly reset envs are re-zeroed in _reset_idx).
        self._x_rel = self._x_rel + velocity * self.step_dt

        # Update pitch-based termination flags (idempotent within one step).
        yaw_error = torch.zeros_like(pitch)
        yaw_rate = torch.zeros_like(pitch)
        self._update_termination_flags(pitch, pitch_rate, velocity, yaw_error, yaw_rate)

        # ── two-phase latch ───────────────────────────────────────────────────
        # Latch "has reached the one-leg stance" the first time roll crosses
        # phase2_roll_deg. Drives both the reward phase and the low-roll
        # termination gate below.
        phase2_roll = math.radians(self.cfg.phase2_roll_deg)
        self._has_reached = self._has_reached | (roll >= phase2_roll)

        # Update roll-based termination counter. Tipping fully over (> high) is
        # always fatal; falling back below the low threshold only counts as a
        # failure AFTER the robot has reached the stance — otherwise it would
        # terminate instantly at the flat (roll ≈ 0) spawn.
        roll_lo = math.radians(self.cfg.roll_low_threshold_deg)
        roll_hi = math.radians(self.cfg.roll_high_threshold_deg)
        too_high = roll.abs() > roll_hi
        too_low = (roll.abs() < roll_lo) & self._has_reached
        outside_window = too_high | too_low
        self._roll_fall_counter = torch.where(
            outside_window, self._roll_fall_counter + 1, torch.zeros_like(self._roll_fall_counter)
        )
        roll_fell = self._roll_fall_counter >= self.cfg.roll_fall_consecutive_steps

        # Phase 1 only: tipping the wrong way (negative roll) past this angle
        # for several consecutive steps is a failure. Counter-based (not
        # instant) — an instant trigger punished the transient negative dips a
        # noisy, exploring policy naturally produces while attempting the tip,
        # which taught it to never attempt the tip at all.
        wrong_direction_roll = math.radians(self.cfg.phase1_wrong_direction_roll_deg)
        past_wrong_direction = (roll < wrong_direction_roll) & ~self._has_reached
        self._wrong_direction_counter = torch.where(
            past_wrong_direction, self._wrong_direction_counter + 1, torch.zeros_like(self._wrong_direction_counter)
        )
        wrong_direction = self._wrong_direction_counter >= self.cfg.phase1_wrong_direction_consecutive_steps
        # Stored for _get_dones, which does not recompute roll itself (mirrors
        # the _roll_fall_counter pattern). NOTE: DirectRLEnv.step() calls
        # _get_dones() BEFORE _get_rewards(), so _get_dones consumes the flag
        # set here on the FOLLOWING step -- termination lags by one step.
        self._wrong_direction = wrong_direction

        terminal_penalty = self._last_terminal_penalty.clone()
        terminal_penalty = torch.where(
            roll_fell & ~self._last_fall,
            torch.full_like(terminal_penalty, self.cfg.fall_penalty),
            terminal_penalty,
        )
        terminal_penalty = torch.where(
            wrong_direction & ~self._last_fall,
            torch.full_like(terminal_penalty, self.cfg.fall_penalty),
            terminal_penalty,
        )

        # ── right-wheel-touch termination (real contact sensor) ────────────────
        # Net contact force on the right wheel body, from the PhysX ContactSensor
        # forced on for this task in __init__ (NOT a joint-extension proxy).
        right_wheel_force = torch.linalg.vector_norm(
            self.wheel_contacts.data.net_forces_w[:, self._right_wheel_contact_id[0]], dim=-1
        )
        right_wheel_touching = right_wheel_force > self.cfg.right_wheel_touch_force_n
        self._right_wheel_touch_counter = torch.where(
            right_wheel_touching, self._right_wheel_touch_counter + 1, torch.zeros_like(self._right_wheel_touch_counter)
        )
        right_wheel_touched = self._right_wheel_touch_counter >= self.cfg.right_wheel_touch_consecutive_steps
        # Stored for _get_dones (same cross-call pattern as _wrong_direction).
        self._right_wheel_touched = right_wheel_touched
        terminal_penalty = torch.where(
            right_wheel_touched & ~self._last_fall,
            torch.full_like(terminal_penalty, self.cfg.fall_penalty),
            terminal_penalty,
        )

        # ── reward components ─────────────────────────────────────────────────
        cfg = self.cfg
        alive = torch.ones_like(roll) * cfg.rew_alive
        # Pitch stability applies in BOTH phases — the robot must not fall
        # fore/aft while climbing or balancing.
        r_pitch = -cfg.rew_pitch * pitch.pow(2)
        r_pitch_rate = -cfg.rew_pitch_rate * pitch_rate.pow(2)
        r_current = -cfg.rew_current * self._action_processor.command_current.pow(2).sum(dim=1)

        # Joint action/smoothness penalties are split bottom vs top. Bottom
        # (FL, BL) are the fine-balance legs and stay penalized in both
        # phases. Top (FR, BR) are the tip-up actuators — they must be free to
        # move fast and hard to lever the body up in phase 1 (a smoothness
        # penalty there directly suppresses the "jerk" motion the tip-up
        # needs and rewards the policy for just standing still on both
        # wheels). The top penalty only switches on in phase 2, once the
        # stance is reached and joint motion should return to fine control.
        # _cur_joint_actions order: [FL, BL, FR, BR].
        bottom_actions = self._cur_joint_actions[:, 0:2]
        top_actions = self._cur_joint_actions[:, 2:4]
        bottom_delta = (self._cur_joint_actions - self._prev_joint_actions)[:, 0:2]
        top_delta = (self._cur_joint_actions - self._prev_joint_actions)[:, 2:4]

        r_joint_action_bottom = -cfg.rew_joint_action * bottom_actions.pow(2).sum(dim=1)
        r_joint_smooth_bottom = -cfg.rew_joint_smooth * bottom_delta.pow(2).sum(dim=1)
        r_joint_action_top = -cfg.rew_joint_action * top_actions.pow(2).sum(dim=1)
        r_joint_smooth_top = -cfg.rew_joint_smooth * top_delta.pow(2).sum(dim=1)

        # Top-joint (FR/BR) extension, in degrees, kept only for observability
        # (logged as top_joint_ext_deg). The soft shaping penalty and the hard
        # extension-limit termination that used to sit here are both removed —
        # neither is part of the specified reward table, and the "must fold to
        # lift the right wheel" constraint is no longer enforced.
        ext_fr_deg = cfg.top_joint_center_deg + top_actions[:, 0] * cfg.top_joint_range_deg
        ext_br_deg = cfg.top_joint_center_deg + top_actions[:, 1] * cfg.top_joint_range_deg
        top_ext_deg = torch.stack([ext_fr_deg, ext_br_deg], dim=1)

        # ── SINGLE-PHASE REWARD ─────────────────────────────────────────────────
        # The robot now spawns already in the one-leg stance (see _reset_idx),
        # so there is no "reach" phase to reward toward the stance — only
        # "balance" (hold it) applies, for the whole episode. The old two-phase
        # reach/balance split is commented out below rather than deleted, in
        # case a future curriculum stage spawns below phase2_roll_deg and needs
        # it back. self._has_reached / phase1_wrong_direction_* still gate
        # TERMINATION only (the "fell back below roll_low" check must not fire
        # before the robot has ever been in the stance) -- this section is
        # reward-only.
        #
        # reach_target = math.radians(cfg.reach_roll_target_deg)
        # r_reach = -cfg.rew_reach * (roll - reach_target).pow(2)
        # r_reach_rate = cfg.rew_reach_rate * roll_rate

        # Balance signals (roll_rate / velocity / station-keeping). No
        # roll-angle target here — the policy finds its own equilibrium around
        # ~47° via the roll-rate penalty.
        r_roll_rate = -cfg.rew_roll_rate * roll_rate.pow(2)
        r_velocity = -cfg.rew_velocity * velocity.pow(2)
        r_position = -cfg.rew_position * self._x_rel.pow(2)

        # Matches the specified reward table exactly (12 rows: Preživetje,
        # Naklon, Sprememba naklona, Ukaz levih/desnih sklepov, Hitra sprememba
        # levih/desnih sklepov, Tok koles, Predčasna prekinitev, Linearna
        # hitrost, Pozicija, Sprememba nagiba) -- no top-joint-extension term.
        common = alive + r_pitch + r_pitch_rate + r_joint_action_bottom + r_joint_smooth_bottom + r_current
        # reward_reach = common + r_reach + r_reach_rate
        reward_balance = (
            common + r_roll_rate + r_velocity + r_position
            + r_joint_action_top + r_joint_smooth_top
        )
        reward = reward_balance + terminal_penalty
        # For logging only: report the top-joint penalty that was actually applied this step.
        r_joint_action = r_joint_action_bottom + r_joint_action_top
        r_joint_smooth = r_joint_smooth_bottom + r_joint_smooth_top
        reward = torch.nan_to_num(reward, nan=0.0, posinf=1.0, neginf=-100.0).clamp(-100.0, 1.0)

        self._episode_reward += reward
        self.extras["log"] = {
            "reward": reward.mean(),
            "episode_reward": self._episode_reward.mean(),
            "reward_alive": alive.mean(),
            "reward_roll_rate_penalty": r_roll_rate.mean(),
            "reward_pitch_penalty": r_pitch.mean(),
            "reward_pitch_rate_penalty": r_pitch_rate.mean(),
            # "reward_reach": r_reach.mean(),            # single-phase reward: reach phase disabled
            # "reward_reach_rate": r_reach_rate.mean(),  # single-phase reward: reach phase disabled
            "reached_rate": self._has_reached.float().mean(),
            "reward_velocity_penalty": r_velocity.mean(),
            "reward_position_penalty": r_position.mean(),
            "reward_joint_action_penalty": r_joint_action.mean(),
            "reward_joint_smooth_penalty": r_joint_smooth.mean(),
            "top_joint_action_abs": top_actions.abs().mean(),
            "top_joint_ext_deg": top_ext_deg.mean(),
            "reward_current_penalty": r_current.mean(),
            "reward_terminal_penalty": terminal_penalty.mean(),
            "roll_abs_deg": roll.abs().mean() * (180.0 / math.pi),
            "roll_rate_abs": roll_rate.abs().mean(),
            "pitch_abs_deg": pitch.abs().mean() * (180.0 / math.pi),
            "pitch_rate_abs": pitch_rate.abs().mean(),
            "velocity_abs": velocity.abs().mean(),
            "x_rel_abs": self._x_rel.abs().mean(),
            "roll_fall_rate": roll_fell.float().mean(),
            "termination_wrong_direction": wrong_direction.float().mean(),
            "right_wheel_force_n": right_wheel_force.mean(),
            "termination_right_wheel_touch": right_wheel_touched.float().mean(),
            "termination_fall": self._last_fall.float().mean(),
            "termination_physics_broken": self._last_physics_broken.float().mean(),
            "termination_invalid_state": self._last_invalid_state.float().mean(),
            "termination_timeout": self._last_timeout.float().mean(),
            "current_rms": torch.sqrt(self._action_processor.command_current.pow(2).mean()),
            "disturbance_force_n": self._last_disturbance_force.norm(dim=1).mean(),
        }
        return reward

    # ── dones ─────────────────────────────────────────────────────────────────

    def _get_dones(self) -> tuple[torch.Tensor, torch.Tensor]:
        grav_b = self.bno080.data.projected_gravity_b
        pitch = pitch_from_projected_gravity(grav_b)
        pitch_rate = -self.bno080.data.ang_vel_b[:, 0]

        velocity = self._compute_velocity()
        yaw_error = torch.zeros_like(pitch)
        yaw_rate = torch.zeros_like(pitch)
        self._update_termination_flags(pitch, pitch_rate, velocity, yaw_error, yaw_rate)

        # _roll_fall_counter is already gated by the two-phase latch in
        # _get_rewards (low-roll only counts after the stance is reached), so a
        # flat-spawned env is never terminated for being "too upright" before it
        # has climbed up.
        roll_fell = self._roll_fall_counter >= self.cfg.roll_fall_consecutive_steps
        terminated = (
            self._last_fall | self._last_physics_broken | self._last_invalid_state
            | roll_fell | self._wrong_direction | self._right_wheel_touched
        )
        timeout = self._last_timeout
        return terminated, timeout

    # ── reset ─────────────────────────────────────────────────────────────────

    def _reset_idx(self, env_ids: Sequence[int] | None):
        # Parent resets motor randomization, action processor, disturbances,
        # pitch observation delay/bias, episode reward, and termination flags.
        # It also writes a pitch-only root pose and sets CG joints to 0 — we
        # override both of those below.
        super()._reset_idx(env_ids)

        if env_ids is None:
            env_ids = self.robot._ALL_INDICES
        env_ids_t = (
            env_ids
            if isinstance(env_ids, torch.Tensor)
            else torch.tensor(env_ids, device=self.device, dtype=torch.long)
        )
        n = len(env_ids_t)

        # ── one-leg-stance spawn: roll = center ± curriculum band, pitch = 0 ──
        # Sample a per-env roll uniformly in [center - band, center + band]; the
        # band widens with the training iteration (spawn_roll_curriculum). Roll
        # is toward the one-leg equilibrium (positive), so the robot spawns
        # already tilted onto the left wheel and must learn to hold it.
        band = self._current_spawn_roll_band_rad()
        roll = self._spawn_roll_center_rad + torch.empty(n, device=self.device).uniform_(-band, band)
        pitch = torch.full((n,), math.radians(self.cfg.reset_pitch_deg), device=self.device)

        # ── spawn height ──────────────────────────────────────────────────────
        # Fixed generous margin (spawn_upright_z + spawn_clearance +
        # spawn_lateral_h), the value the original one-leg reset and
        # scripts/lqr_control_one_leg_jump.py both use. It sits STRICTLY above
        # ground for the extended-leg pose so the robot never spawns clipped;
        # it settles the last few mm under gravity within the first control
        # step. A tilt-aware (cos/sin) formula sits ~6 cm lower at 47° and spawns
        # the robot penetrating the floor -> PhysX ejects it at up to
        # max_depenetration_velocity (1 m/s), which at this tilt is mostly
        # horizontal and makes the robot slide off on spawn. Keep the margin.
        spawn_z = self.cfg.spawn_upright_z + self.cfg.spawn_clearance + self.cfg.spawn_lateral_h

        # ── combined roll + pitch quaternion (wxyz) ───────────────────────────
        # Compose: q_pitch * q_roll (pitch applied after roll in world frame).
        # q_roll = [cos(r/2), 0, sin(r/2), 0], q_pitch = [cos(p/2), -sin(p/2), 0, 0]
        hp = pitch * 0.5
        hr = roll * 0.5
        qw = torch.cos(hp) * torch.cos(hr)
        qx = -torch.sin(hp) * torch.cos(hr)
        qy = torch.cos(hp) * torch.sin(hr)
        qz = -torch.sin(hp) * torch.sin(hr)

        root_state = self.robot.data.default_root_state[env_ids_t].clone()
        root_state[:, :3] += self.scene.env_origins[env_ids_t]
        root_state[:, 2] = self.scene.env_origins[env_ids_t, 2] + spawn_z
        root_state[:, 3] = qw
        root_state[:, 4] = qx
        root_state[:, 5] = qy
        root_state[:, 6] = qz
        root_state[:, 7:] = 0.0

        # Small initial pitch rate.
        pitch_rate_range = self.cfg.reset_pitch_rate_range_radps
        root_state[:, 10] = torch.empty(n, device=self.device).uniform_(
            -pitch_rate_range, pitch_rate_range
        )
        self.robot.write_root_pose_to_sim(root_state[:, :7], env_ids_t)
        self.robot.write_root_velocity_to_sim(root_state[:, 7:], env_ids_t)

        # ── CG joint positions: one-leg stance (asymmetric) ───────────────────
        # Left legs (FL/BL) grounded at spawn_left_ext_deg, right legs (FR/BR)
        # retracted to spawn_right_ext_deg so the right wheel lifts off. Raw
        # angles precomputed in __init__ (_spawn_raw_*).
        cg_pos = torch.zeros(n, 4, device=self.device)
        cg_pos[:, _IDX_FL] = self._spawn_raw_fl
        cg_pos[:, _IDX_FR] = self._spawn_raw_fr
        cg_pos[:, _IDX_BL] = self._spawn_raw_bl
        cg_pos[:, _IDX_BR] = self._spawn_raw_br
        cg_vel = torch.zeros(n, 4, device=self.device)
        self.robot.write_joint_state_to_sim(cg_pos, cg_vel, self._cg_ids, env_ids_t)
        self.robot.set_joint_position_target(cg_pos, env_ids=env_ids_t, joint_ids=self._cg_ids)

        # ── Passive left-side bearing joints: set to the loop-closure-consistent
        # angles (see __init__ for the FK derivation), not the USD's
        # unauthored 0 rad default. Only Revolute_9/Revolute_10 are directly
        # writable (tochangeleft is an implicit PhysX loop-closure constraint,
        # not a queryable joint -- see __init__). Right-side bearings are left
        # untouched -- their FK solution for the 0/0 retracted spawn is ~0 rad,
        # matching the default, so no override is needed there.
        for bearing_ids, target_rad in zip(
            (self._left_bearing_ids, self._left_bearing3_ids),
            self._left_bearing_targets_rad,
        ):
            pos = torch.full((n, 1), target_rad, device=self.device)
            vel = torch.zeros(n, 1, device=self.device)
            self.robot.write_joint_state_to_sim(pos, vel, bearing_ids, env_ids_t)

        # ── CG joint stiffness floor ──────────────────────────────────────────
        # The parent reset (standup_env) writes stiffness from its domain
        # randomization (cg_fixed_kp=60 by default) via
        # write_joint_stiffness_to_sim(). cybergear_joints is an
        # ImplicitActuatorCfg (robot_cfg.py), so robot.data.joint_stiffness
        # reflects the real PhysX DOF-drive stiffness and
        # write_joint_stiffness_to_sim() is the correct path. Clamp all four
        # CG joints to at least min_cg_stiffness_nm_rad so a low-stiffness
        # config cannot prevent them from holding their commanded positions.
        env_ids_cpu = env_ids_t.cpu()
        stiffness = self.robot.data.joint_stiffness[env_ids_t].clone()
        min_kp = self.cfg.min_cg_stiffness_nm_rad
        for jid in self._cg_ids:
            stiffness[:, jid] = stiffness[:, jid].clamp(min=min_kp)
        self.robot.write_joint_stiffness_to_sim(stiffness, env_ids=env_ids_cpu)

        # ── one-leg specific state ────────────────────────────────────────────
        self._roll_fall_counter[env_ids_t] = 0
        self._prev_joint_actions[env_ids_t] = 0.0
        self._cur_joint_actions[env_ids_t] = 0.0
        self._x_rel[env_ids_t] = 0.0
        # Spawn is already in the one-leg stance, so latch phase-2 (balance)
        # immediately for any env spawned at/above phase2_roll_deg — this
        # activates the balance-phase reward and the "fell back below roll_low"
        # termination from the first step. Envs sampled below the threshold
        # (only possible at the widest band) start in reach phase and latch
        # once roll climbs past it, exactly as before.
        phase2_roll = math.radians(self.cfg.phase2_roll_deg)
        self._has_reached[env_ids_t] = roll >= phase2_roll
        self._wrong_direction_counter[env_ids_t] = 0
        self._wrong_direction[env_ids_t] = False
        self._right_wheel_touch_counter[env_ids_t] = 0
        self._right_wheel_touched[env_ids_t] = False
