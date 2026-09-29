"""Configuration for the one-leg balance environment."""

import math

from isaaclab.sim import SimulationCfg
from isaaclab.utils import configclass

from .pure_nn_balance_env_cfg import PureNNBalanceEnvCfg


@configclass
class OneLegBalanceEnvCfg(PureNNBalanceEnvCfg):

    # ── Physics timing ────────────────────────────────────────────────────────
    # 1ms × 20 substeps = 20ms control step (50Hz) → 1000Hz physics.
    decimation: int = 20
    sim: SimulationCfg = SimulationCfg(dt=0.001, render_interval=20)

    # ── RL interface ──────────────────────────────────────────────────────────
    observation_space: int = 6
    action_space: int = 5       # [wheel_L, joint_FL, joint_BL, joint_FR, joint_BR]  (right wheel lifted — always zero)
    episode_length_s: float = 10.0

    # ── Roll termination window ───────────────────────────────────────────────
    # Robot must stay within this roll window to survive.
    # Too upright  (< low)  → fell back to two-leg stance.
    # Too tilted   (> high) → tipped completely over.
    roll_low_threshold_deg: float = 20.0
    roll_high_threshold_deg: float = 72.0
    # Consecutive steps outside the window before termination (matches parent fall logic).
    roll_fall_consecutive_steps: int = 5

    # ── Right-wheel-touch termination (real contact sensor, not joint extension) ──
    # The right wheel must stay lifted. Detected via a PhysX ContactSensor on
    # the DDSM115_Simplified_01 (right wheel) body -- see __init__, which forces
    # enable_wheel_contacts/activate_contact_sensors on for this task. Net
    # contact force above this threshold for this many consecutive control
    # steps ends the episode with fall_penalty (same "Predčasna prekinitev"
    # weight as every other termination cause, not a separate reward term).
    right_wheel_touch_force_n: float = 0.1   # N -- matches the "touching floor" threshold used elsewhere (lqr_control_one_leg_jump.py)
    right_wheel_touch_consecutive_steps: int = 5

    # Phase 1 only: the robot must tip toward +reach_roll_target_deg. If roll
    # instead goes negative (wrong direction) past this angle for several
    # consecutive steps, end the episode rather than let it wander/settle on
    # the wrong side. Widened from an earlier -15 deg with an instant trigger —
    # that combination punished the transient negative dips a noisy, exploring
    # policy naturally produces while attempting the tip, which discouraged it
    # from ever trying and left it stuck balancing flat in phase 1 forever.
    # A counter (not instant) plus a looser angle gives room for a real attempt
    # to recover before being killed. Not checked once has_reached is latched
    # (phase 2), since balance dynamics can legitimately dip roll rate through
    # small negative excursions.
    phase1_wrong_direction_roll_deg: float = -30.0
    phase1_wrong_direction_consecutive_steps: int = 5

    # Pitch termination is inherited from parent (fall_pitch_threshold_deg = 25°).
    # For the tilted one-leg stance, pitch dynamics are less severe so loosen it.
    fall_pitch_threshold_deg: float = 40.0

    # ── Joint configuration ───────────────────────────────────────────────────
    # All four CyberGear joints are policy-controlled (5D action). The top
    # joints (front_right, back_right) are the tip-up actuators: they need a
    # WIDE range so the policy can extend them (to lever the body up toward the
    # one-leg equilibrium) and fold them (toward 0°, the resting one-leg stance).
    # The bottom joints (front_left, back_left) are the fine balance legs and
    # keep a narrow range around their center.
    # Initial (reset) extension of all four joints for the flat symmetric spawn.
    reset_joint_ext_deg: float = 45.0

    # CyberGear position gain (kp) for all four leg joints during training.
    # Overrides the parent standup default (60) for the one-leg task.
    cg_fixed_kp: float = 30.0        # Nm/rad

    # Minimum CyberGear stiffness enforced for all four leg joints at each
    # episode reset. Guards against future low-stiffness tuning breaking the
    # position-hold for the active joints.
    min_cg_stiffness_nm_rad: float = 3.0

    # Bottom joints (front_left, back_left) controllable range.
    # Policy action [-1, 1] maps linearly to [center - range, center + range].
    # Expressed as extension degrees (always positive regardless of raw sign).
    # Widened to match top_joint_range_deg -- the previous ±15° cap is removed;
    # the hard physical joint limit (raw_fl/raw_bl.clamp_() in
    # _pre_physics_step, from the CyberGear's actual raw travel) is still the
    # real safety bound, this only removes the artificial policy-side cap.
    bottom_joint_center_deg: float = 45.0
    bottom_joint_range_deg: float = 45   # policy can command the full range, clamped by hardware limits

    # Top joints (front_right, back_right) controllable range — wide, for tipping.
    top_joint_center_deg: float = 45.0
    top_joint_range_deg: float = 45      # policy can command 0°..90° extension

    # ── Reset — spawn IN the one-leg stance, roll-band curriculum ─────────────
    # The robot spawns already at the one-leg equilibrium (~47° roll) balanced
    # on the LEFT wheel, right legs (FR/BR) retracted so the right wheel lifts,
    # and must learn (pure RL, no LQR) to hold the stance. Each episode samples
    # roll = spawn_roll_center_deg ± band, where band WIDENS with the training
    # iteration per spawn_roll_curriculum. pitch spawns at 0° with a small rate.
    spawn_roll_center_deg: float = 47.0
    # Iteration-based curriculum: (start_iteration, roll_half_band_deg). The
    # active stage is the last one whose start_iteration <= current iteration.
    # 0–350: ±5°, 350–700: ±8°, 700+: ±11°.
    spawn_roll_curriculum: tuple = ((0, 5.0), (350, 8.0), (700, 11.0))
    # Rollout length used to convert env step count -> training iteration for
    # the curriculum. MUST match num_steps_per_env in the rsl-rl agent cfg
    # (rsl_rl_one_leg_balance_cfg.py); a mismatch silently mis-times the stages.
    curriculum_num_steps_per_env: int = 64

    # Spawn leg pose (extension degrees, always positive; ext_sign applied in
    # the env). Right legs retract to lift the right wheel; left legs hold the
    # grounded balance stance.
    spawn_right_ext_deg: float = 0.0    # FR/BR retracted (right wheel lifts off)
    spawn_left_ext_deg: float = 45.0    # FL/BL grounded balance legs

    reset_pitch_deg: float = 0.0
    # Spawn at rest — zero initial pitch rate so the robot doesn't lurch/slide on
    # spawn. Raise this to reintroduce a small pitch perturbation for robustness.
    reset_pitch_rate_range_radps: float = 0.0
    reset_velocity_range_mps: float = 0.05

    # ── Two-phase reward ──────────────────────────────────────────────────────
    # Phase 1 ("reach"): from the flat spawn, pull roll UP toward the target.
    # Phase 2 ("balance"): once roll first crosses phase2_roll_deg (latched for
    # the rest of the episode), switch to the balance rewards below. The switch
    # is per-env and one-way — once balancing, a brief dip back below the
    # threshold does not revert to the reach reward.
    reach_roll_target_deg: float = 47.0   # roll angle the reach phase pulls toward
    phase2_roll_deg: float = 40.0         # balance reward kicks in at this roll
    # Steepened from 2.0. At 2.0 the policy camped at its ~15° static leg-
    # extension ceiling: alive(+1.0) − 2.0·(15°−47°)² ≈ +0.38/step, a safe,
    # profitable local optimum. At 3.5 that same camp is ≈ alive − 3.5·(0.558)²
    # ≈ −0.09/step — no longer profitable, so the policy is pushed to climb —
    # while still strictly better than falling (−100 ended immediately vs
    # ≈ −0.09/step ≈ −40 over the rest of the episode), so it should NOT learn
    # to suicide instead. Climbing pays off steeply: 30° → ≈ +0.69/step,
    # 47° → +1.0/step.
    rew_reach: float = 3.5                # weight on -(roll - target)^2 in phase 1
    # r_reach alone only rewards being CLOSE to the target, not moving toward
    # it — a weak gradient for an initially-random policy to stumble into a
    # large coordinated tipping motion. This adds a direct, linear bonus for
    # positive roll rate (climbing toward the target) — and an equal penalty
    # for negative roll rate (retreating) — so any motion in the right
    # direction pays off immediately, not just the end state.
    rew_reach_rate: float = 1.5           # weight on +roll_rate in phase 1

    # ── Balance reward weights (phase 2) ──────────────────────────────────────
    # Primary balance signal: minimize roll rate. The policy naturally discovers
    # the stable angle (~47°) without being told where to stand.
    rew_roll_rate: float = 2.0

    # Pitch and forward-motion regulation (inherited wheels + LQR style).
    rew_pitch: float = 5.0
    rew_pitch_rate: float = 0.3
    # Velocity/position are now observable (obs[11]/obs[10]) so the policy can
    # regulate them directly.  Keep these penalties GENTLE — an inverted pendulum
    # must drive the wheel under its CoM to catch a fall, and heavy velocity or
    # position penalties suppress exactly that recovery motion.  Priority order
    # is survival ≫ balance (pitch/roll_rate) ≫ station-keeping.
    rew_velocity: float = 0.4
    # Station keeping: penalize drift from the episode start position
    # (x_rel integrated from left-wheel odometry, same as STM32 deployment).
    rew_position: float = 0.15

    # Joint command penalties — keep them small so they shape but do not dominate.
    rew_joint_action: float = 0.1
    rew_joint_smooth: float = 2.0

    # Wheel current penalty (inherited).
    rew_current: float = 0.003

    # Alive reward — constant survival signal.
    rew_alive: float = 1.0

    # Terminal penalties.
    fall_penalty: float = -100.0
    physics_broken_penalty: float = -100.0
    invalid_state_penalty: float = -100.0

    # ── Observation normalization (6D) ────────────────────────────────────────
    # Order: roll, roll_rate, pitch, pitch_rate, x_rel, velocity.
    # The GRU hidden state carries action/current history internally, so the
    # previous-action and previous-current signals are no longer fed in.
    observation_scale: tuple = (
        math.radians(50.0),   # [0]  roll  — normalized by ~tipping angle
        4.0,                   # [1]  roll rate
        math.radians(25.0),   # [2]  pitch
        4.0,                   # [3]  pitch rate
        0.5,                   # [4]  x_rel — drift from start (m), left-wheel odometry
        0.5,                   # [5]  linear velocity (m/s), left wheel only
    )

    # Disable yaw rate from parent obs builder (obs is rebuilt from scratch here).
    include_yaw_rate: bool = False
