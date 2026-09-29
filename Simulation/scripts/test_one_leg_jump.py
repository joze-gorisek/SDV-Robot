"""
One-leg DYNAMIC TIP feasibility diagnostic — open-loop, no RL policy.

Answers a hard go/no-go question: can this robot reach the one-leg stance
(~47 deg roll) at all? A *static* CoM-shift (slowly extending the right legs)
caps at a ~15-20 deg roll ceiling — the trained policy camps exactly there.
Reaching 40 deg+ requires a *dynamic, inertial* leg thrust that overshoots
that static equilibrium and carries the CoM past the tipping point (where the
left wheel becomes the sole pivot). This script tests whether such a thrust is
physically possible.

A two-wheeled robot is PITCH-unstable, so it cannot simply be spawned and left
to sit upright — it topples fore/aft before any tip test can run. Each trial
therefore begins with an LQR pitch-stabilization warmup (the exact 4-state
controller from scripts/lqr_control.py) that holds the robot upright on both
wheels; the LQR stays active THROUGH the thrust so the robot can't fall fore/
aft, isolating the roll-tip question. Only pitch is stabilized — the LQR never
acts on roll.

How it works: from the stabilized upright state, the same SHIFT/FOLD state
machine extends the right legs (FR, BR), but now the ramp SPEED and joint
STIFFNESS are swept from slow/soft (quasi-static) to fast/stiff (a genuine
inertial thrust). A slow reference trial establishes the static ceiling; every
dynamic trial reports its peak roll relative to that ceiling. Optionally a brief
extra LEFT-wheel pulse (--wheel-assist-nm, on top of the LQR effort, capped at
the ~2.0 Nm hardware ceiling) tests whether forward acceleration couples into
the tip. Peak CyberGear torque is reported so a thrust that silently saturates
the 12 Nm limit is visible, not assumed.

Feasibility verdict per setting:
  FEASIBLE    peak roll crosses the overshoot bar (default 40 deg), repeatably,
              with leg torque <= 12 Nm (and assist <= 2 Nm).
  PARTIAL     dynamic overshoot over the static ceiling but short of the bar.
  STATIC ONLY no dynamic gain over the ceiling.

The script bypasses env.step() and drives the articulation + physics loop
directly, so no RL termination/reset logic can interfere mid-tip.

Usage (inside the Isaac Lab container, from /workspace):

    # Sweep thrust SPEED and STIFFNESS to find whether any dynamic thrust crosses 40 deg
    isaaclab/isaaclab.sh -p Project/scripts/test_one_leg_jump.py \\
        --task Template-Twowheeledrobot-Standup-v0 --right-ext-deg 90 \\
        --shift-duration-list 0.4,0.15,0.08,0.05,0.03 --thrust-kp-list 60,120,250 \\
        --overshoot-target-deg 40 --trials 5 --headless \\
        --csv /workspace/Project/logs/dynamic_tip_feasibility.csv

    # If pure-leg thrust falls short, add a timed left-wheel assist pulse
    ... --shift-duration-list 0.05 --thrust-kp-list 250 \\
        --wheel-assist-nm 2.0 --wheel-assist-window-s 0.05 --wheel-assist-sign 1
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

_EXTENSION_SOURCE_PATH = Path(__file__).resolve().parents[1] / "source" / "TwoWheeledRobot"
if _EXTENSION_SOURCE_PATH.is_dir():
    sys.path.insert(0, str(_EXTENSION_SOURCE_PATH))

from isaaclab.app import AppLauncher

parser = argparse.ArgumentParser(description="One-leg CoM-shift tipping diagnostic (no velocity, no wheel current).")
parser.add_argument("--task", type=str, default="Template-Twowheeledrobot-Standup-v0")
parser.add_argument("--target-roll-deg", type=float, default=47.0, help="Roll angle the tip should reach.")
parser.add_argument("--start-ext-deg", type=float, default=45.0, help="Initial (at-rest) extension of all four joints.")
parser.add_argument("--right-ext-deg", type=float, default=90.0, help="FR/BR extension target during the slow shift.")
parser.add_argument(
    "--right-ext-list", type=str, default=None,
    help="Comma-separated FR/BR shift targets to sweep, e.g. '60,70,80,90'. Overrides --right-ext-deg.",
)
parser.add_argument("--shift-duration-s", type=float, default=0.4, help="Time to ramp FR/BR from start to right-ext-deg (slow = quasi-static, not a push).")
parser.add_argument(
    "--shift-duration-list", type=str, default=None,
    help="Comma-separated shift durations (s) to sweep, e.g. '0.4,0.15,0.08,0.05,0.03'. "
         "Primary DYNAMIC-thrust knob — a 0.03 s ramp of 45->90 deg is a genuine inertial thrust. "
         "Overrides --shift-duration-s.",
)
parser.add_argument(
    "--thrust-kp-list", type=str, default=None,
    help="Comma-separated CyberGear stiffness (Nm/rad) values to sweep during the thrust, "
         "e.g. '60,120,250'. Higher kp = more torque authority to accelerate the legs "
         "(bounded by the 12 Nm CyberGear peak). Overrides --hold-kp for the trial.",
)
parser.add_argument(
    "--overshoot-target-deg", type=float, default=40.0,
    help="Feasibility bar: peak roll must exceed this to count as crossing the tipping point.",
)
parser.add_argument(
    "--wheel-assist-nm", type=float, default=0.0,
    help="Peak LEFT-wheel effort (Nm) during the thrust, to test whether a forward pulse couples "
         "into the tip. Clamped to the ~2.0 Nm hardware ceiling.",
)
parser.add_argument(
    "--wheel-assist-window-s", type=float, default=None,
    help="Duration of the left-wheel assist pulse from t=0 (default: the shift duration).",
)
parser.add_argument(
    "--wheel-assist-sign", type=int, default=1, choices=(-1, 1),
    help="Sign/direction of the forward wheel pulse; sweep both to see which couples into +roll.",
)
parser.add_argument(
    "--fold-roll-deg", type=float, default=90.0,
    help="Once roll reaches this angle, start retracting FR/BR to 0 deg. Default 90 deg effectively "
         "DISABLES the fold during a feasibility climb, so peak roll measures the thrust's unimpeded "
         "reach (folding at a low angle would cap the dynamic overshoot). Lower it only to test the "
         "full shift->fold stance maneuver.",
)
parser.add_argument("--fold-duration-s", type=float, default=0.3, help="Time to ramp FR/BR from right-ext-deg to 0 deg.")
parser.add_argument("--hold-kp", type=float, default=60.0, help="CyberGear stiffness for all four joints (settle + trials when --thrust-kp-list is unset).")
parser.add_argument(
    "--lqr-warmup-s", type=float, default=1.0,
    help="Seconds of LQR pitch-stabilization to hold the robot upright on both wheels BEFORE the "
         "thrust. LQR stays active through the thrust so the robot can't fall fore/aft, isolating "
         "the roll-tip question.",
)
parser.add_argument(
    "--warmup-clearance", type=float, default=0.0,
    help="Extra height (m) above the env's proven flat-spawn height (spawn_upright_z + clearance + "
         "lateral_h) to drop from; LQR + gravity settle it onto the wheels during the warmup. "
         "Default 0 = minimal drop. Raise slightly if it ever spawns clipping.",
)
parser.add_argument("--flight-s", type=float, default=3.0, help="Observation window after the shift begins.")
parser.add_argument("--liftoff-mm", type=float, default=5.0, help="Right-wheel rise that counts as liftoff.")
parser.add_argument("--tip-roll-deg", type=float, default=72.0, help="Roll beyond this counts as tipped over.")
parser.add_argument("--trials", type=int, default=3, help="Repeats per right-ext setting.")
parser.add_argument("--print-interval-s", type=float, default=0.0, help="Seconds between prints (0 = print every physics step).")
parser.add_argument("--csv", type=str, default=None, help="Optional CSV path for per-step logs.")
AppLauncher.add_app_launcher_args(parser)
args_cli, hydra_args = parser.parse_known_args()
args_cli.num_envs = 1
sys.argv = [sys.argv[0]] + hydra_args

# Left-wheel effort ceiling: i_max 2.0 A x DDSM115_KT 0.75 Nm/A ~= 1.5 Nm sustained,
# ~2.0 Nm peak. Reject unphysical assist so feasibility can't be claimed with torque
# the real motor can't deliver.
_WHEEL_ASSIST_MAX_NM = 2.0
_CG_PEAK_TORQUE_NM = 12.0
if abs(args_cli.wheel_assist_nm) > _WHEEL_ASSIST_MAX_NM:
    print(
        f"[jump-diag] WARNING: --wheel-assist-nm {args_cli.wheel_assist_nm:.2f} exceeds the "
        f"{_WHEEL_ASSIST_MAX_NM:.1f} Nm hardware ceiling — clamping.",
        flush=True,
    )
    args_cli.wheel_assist_nm = math.copysign(_WHEEL_ASSIST_MAX_NM, args_cli.wheel_assist_nm)

app_launcher = AppLauncher(args_cli)
simulation_app = app_launcher.app

# ─── everything below runs after Isaac Sim has launched ──────────────────────

import gymnasium as gym
import torch

import isaaclab_tasks  # noqa: F401
import TwoWheeledRobot.tasks  # noqa: F401
from isaaclab.envs import DirectMARLEnvCfg, DirectRLEnvCfg, ManagerBasedRLEnvCfg
from isaaclab_tasks.utils.hydra import hydra_task_config

# CyberGear index order within _cg_ids: [front_left, front_right, back_left, back_right]
IDX_FL, IDX_FR, IDX_BL, IDX_BR = 0, 1, 2, 3

# ── 4-state pitch LQR (replicated from scripts/lqr_control.py) ────────────────
# A two-wheeled robot is PITCH-unstable; without active balance it tips over
# fore/aft before any tip test can run. These are the exact split 4-state LQR
# per-wheel force gains and sign conventions from lqr_control.py, used to hold
# the robot upright on both wheels during a warmup and throughout the thrust.
LQR_K_POS = 23.8188            # force per wheel per m of wheel position
LQR_K_VEL = 16.7186            # force per wheel per m/s
LQR_K_PITCH = 77.3774          # force per wheel per rad
LQR_K_PITCH_RATE = 4.5646      # force per wheel per rad/s
LQR_WHEEL_FEEDBACK_SIGN = (-1.0, +1.0)   # raw joint encoder -> physical forward-positive
LQR_WHEEL_EFFORT_SIGN = (-1.0, +1.0)     # physical torque -> USD joint effort (mirror)
LQR_PITCH_SENSOR_SIGN = -1.0             # projected-gravity pitch -> LQR pitch convention
LQR_R_WHEEL = 0.05035                    # wheel radius (m)
LQR_TORQUE_LIMIT_NM = 2.0                # DDSM115 peak ~2.7 A * 0.75 Nm/A ~= 2.0 Nm


def roll_from_projected_gravity(g: torch.Tensor) -> torch.Tensor:
    return torch.atan2(g[:, 0], -g[:, 2])


@hydra_task_config(args_cli.task, "rsl_rl_cfg_entry_point")
def main(env_cfg: ManagerBasedRLEnvCfg | DirectRLEnvCfg | DirectMARLEnvCfg, agent_cfg=None):
    env_cfg.scene.num_envs = 1
    env_cfg.seed = 0

    env = gym.make(args_cli.task, cfg=env_cfg)
    env.reset()
    inner = env.unwrapped
    device = inner.device
    dt = float(inner.cfg.sim.dt)
    render = not args_cli.headless

    cg_ids = inner._cg_ids
    ext_sign = inner._cg_ext_sign[0].cpu()          # [+1, -1, -1, +1]
    cg_lo = inner._cg_joint_lo[0].cpu()
    cg_hi = inner._cg_joint_hi[0].cpu()
    env0 = torch.tensor([0], device=device, dtype=torch.long)

    def raw_angle(idx: int, ext_deg: float) -> float:
        """Extension degrees -> raw joint angle (rad), clamped to hard limits."""
        raw = float(ext_sign[idx]) * math.radians(ext_deg)
        return min(max(raw, float(cg_lo[idx])), float(cg_hi[idx]))

    def leg_targets(ext_fl_bl: float, ext_fr_br: float) -> torch.Tensor:
        t = torch.zeros(1, 4, device=device)
        t[0, IDX_FL] = raw_angle(IDX_FL, ext_fl_bl)
        t[0, IDX_BL] = raw_angle(IDX_BL, ext_fl_bl)
        t[0, IDX_FR] = raw_angle(IDX_FR, ext_fr_br)
        t[0, IDX_BR] = raw_angle(IDX_BR, ext_fr_br)
        return t

    # Right wheel body for liftoff detection. Try a few patterns; if none match,
    # print the available body names so the pattern can be fixed (a missing body
    # is why an earlier version reported lift = 0 always).
    right_wheel_body = None
    for pattern in (".*[Dd]esni.*", ".*DDSM115_Desni.*", ".*[Rr]ight.*[Ww]heel.*", ".*wheel.*[Rr].*"):
        try:
            ids, names = inner.robot.find_bodies(pattern)
        except Exception:
            ids, names = [], []
        if len(ids) > 0:
            right_wheel_body = ids[0]
            print(f"[jump-diag] right wheel body: '{names[0]}' (id {right_wheel_body}, pattern '{pattern}')", flush=True)
            break
    if right_wheel_body is None:
        try:
            all_names = inner.robot.body_names
        except Exception:
            all_names = "<unavailable>"
        print(
            "[jump-diag] right wheel body NOT found — liftoff falls back to a roll>5° proxy. "
            f"Available body names: {all_names}",
            flush=True,
        )

    def set_cg_stiffness(kp: float) -> None:
        # cybergear_joints is an explicit IdealPDActuatorCfg (robot_cfg.py): it
        # computes PD torque in Python from the actuator's own .stiffness
        # tensor, not from PhysX's DOF drive, so write_joint_stiffness_to_sim()
        # would be a silent no-op here. Mutate the actuator's tensor directly.
        inner.robot.actuators["cybergear_joints"].stiffness[0:1] = kp

    def place_robot(spawn_z: float, ext_fl_bl: float, ext_fr_br: float) -> None:
        """Teleport upright at roll = 0, legs at the given extensions, ZERO
        velocity (root and joints). No velocity is ever injected anywhere in
        this script."""
        root_state = inner.robot.data.default_root_state[0:1].clone()
        root_state[:, :3] += inner.scene.env_origins[0:1]
        root_state[:, 2] = inner.scene.env_origins[0, 2] + spawn_z
        root_state[:, 3:7] = torch.tensor([1.0, 0.0, 0.0, 0.0], device=device)
        root_state[:, 7:] = 0.0
        inner.robot.write_root_pose_to_sim(root_state[:, :7], env0)
        inner.robot.write_root_velocity_to_sim(root_state[:, 7:], env0)

        cg_pos = leg_targets(ext_fl_bl, ext_fr_br)
        cg_vel = torch.zeros(1, 4, device=device)
        inner.robot.write_joint_state_to_sim(cg_pos, cg_vel, cg_ids, env0)
        inner.robot.set_joint_position_target(cg_pos, joint_ids=cg_ids, env_ids=env0)

        # Wheels: zero state, zero velocity. No current/effort is ever applied
        # to them anywhere in this script — they are always passive.
        wheel_ids = inner._wheel_ids
        wheel_zero = torch.zeros(1, 2, device=device)
        inner.robot.write_joint_state_to_sim(wheel_zero, wheel_zero, wheel_ids, env0)

    def step_sim(targets: torch.Tensor, wheel_effort_L: float = 0.0, wheel_effort_R: float = 0.0) -> None:
        inner.robot.set_joint_position_target(targets, joint_ids=cg_ids, env_ids=env0)
        # Per-wheel effort (Nm, already in USD joint sign). The LQR pitch
        # stabilizer drives both wheels during warmup + thrust; the optional
        # tip assist adds to the left. Default 0.0 keeps a caller fully passive.
        wheel_effort = torch.zeros(1, 2, device=device)
        wheel_effort[0, 0] = wheel_effort_L
        wheel_effort[0, 1] = wheel_effort_R
        inner.robot.set_joint_effort_target(wheel_effort, joint_ids=inner._wheel_ids, env_ids=env0)
        inner.robot.write_data_to_sim()
        inner.sim.step(render=render)
        inner.scene.update(dt)

    def lqr_wheel_efforts() -> tuple[float, float]:
        """Per-wheel effort (Nm, USD joint sign) from the 4-state pitch LQR,
        replicating scripts/lqr_control.py's split controller. Holds the robot
        upright on both wheels; acts only on pitch, never on roll."""
        grav = inner.bno080.data.projected_gravity_b[0]
        pitch = math.atan2(float(grav[1]), -float(grav[2]))
        pitch_rate = -float(inner.bno080.data.ang_vel_b[0, 0])
        lqr_pitch = LQR_PITCH_SENSOR_SIGN * pitch
        lqr_pitch_rate = LQR_PITCH_SENSOR_SIGN * pitch_rate

        raw_pos = inner.robot.data.joint_pos[0, inner._wheel_ids]
        raw_vel = inner.robot.data.joint_vel[0, inner._wheel_ids]
        efforts = [0.0, 0.0]
        for i in (0, 1):
            pos_m = float(raw_pos[i]) * LQR_WHEEL_FEEDBACK_SIGN[i] * LQR_R_WHEEL
            vel_m = float(raw_vel[i]) * LQR_WHEEL_FEEDBACK_SIGN[i] * LQR_R_WHEEL
            force = -(
                LQR_K_POS * pos_m + LQR_K_VEL * vel_m
                + LQR_K_PITCH * lqr_pitch + LQR_K_PITCH_RATE * lqr_pitch_rate
            )
            torque = force * LQR_R_WHEEL   # ground force (N) -> wheel torque (Nm)
            torque = max(-LQR_TORQUE_LIMIT_NM, min(LQR_TORQUE_LIMIT_NM, torque))
            efforts[i] = LQR_WHEEL_EFFORT_SIGN[i] * torque
        return -efforts[0], -efforts[1]

    def lqr_warmup() -> dict:
        """Spawn upright and hold the robot on both wheels with the LQR pitch
        stabilizer for lqr_warmup_s. Replaces the old passive settle, which let
        the pitch-unstable robot topple before any test could run. Returns the
        state at the end of the warmup (checked for a clean upright start)."""
        # Reuse the env's proven flat-spawn height for the 45°-leg pose
        # (spawn_upright_z alone is for FOLDED legs and would clip). This sits a
        # few mm above true contact, so the drop the LQR must absorb is tiny.
        spawn_z = (
            inner.cfg.spawn_upright_z + inner.cfg.spawn_clearance
            + inner.cfg.spawn_lateral_h + args_cli.warmup_clearance
        )
        place_robot(spawn_z, args_cli.start_ext_deg, args_cli.start_ext_deg)
        set_cg_stiffness(args_cli.hold_kp)
        hold = leg_targets(args_cli.start_ext_deg, args_cli.start_ext_deg)
        for _ in range(round(args_cli.lqr_warmup_s / dt)):
            eL, eR = lqr_wheel_efforts()
            step_sim(hold, eL, eR)
        return read_state()

    def read_state() -> dict:
        grav = inner.bno080.data.projected_gravity_b
        roll = float(roll_from_projected_gravity(grav)[0])
        roll_rate = float(inner.bno080.data.ang_vel_b[0, 1])
        pitch = math.atan2(float(grav[0, 1]), -float(grav[0, 2]))
        pitch_rate = -float(inner.bno080.data.ang_vel_b[0, 0])
        body_z = float(inner.robot.data.root_pos_w[0, 2])
        wheel_z = (
            float(inner.robot.data.body_pos_w[0, right_wheel_body, 2])
            if right_wheel_body is not None
            else float("nan")
        )
        # Peak |applied torque| across the four CG joints, so a thrust that
        # silently saturates the 12 Nm CyberGear limit is visible, not assumed.
        cg_torque = inner.robot.data.applied_torque[0, cg_ids].abs()
        max_cg_torque = float(cg_torque.max())
        return {
            "roll": roll, "roll_rate": roll_rate, "pitch": pitch, "pitch_rate": pitch_rate,
            "body_z": body_z, "wheel_z": wheel_z, "max_cg_torque": max_cg_torque,
        }

    def run_trial(right_ext_deg: float, shift_dur: float, kp: float, trial: int,
                  static_ceiling: float | None, csv_rows: list) -> dict:
        target_roll = math.radians(args_cli.target_roll_deg)
        tip_roll = math.radians(args_cli.tip_roll_deg)
        fold_roll = math.radians(args_cli.fold_roll_deg)
        overshoot_target = args_cli.overshoot_target_deg
        assist_nm = args_cli.wheel_assist_nm * args_cli.wheel_assist_sign
        assist_win = args_cli.wheel_assist_window_s if args_cli.wheel_assist_window_s is not None else shift_dur
        print_every = 1 if args_cli.print_interval_s <= 0.0 else max(1, round(args_cli.print_interval_s / dt))

        # LQR pitch-stabilization warmup: spawn upright and hold the robot on
        # both wheels until it is stable, THEN thrust from that clean state.
        base = lqr_warmup()
        set_cg_stiffness(kp)
        wheel_z0 = base["wheel_z"]
        warmup_ok = abs(base["pitch"]) < math.radians(10.0) and abs(base["roll"]) < math.radians(10.0)

        m = {
            "right_ext_deg": right_ext_deg,
            "shift_dur": shift_dur,
            "kp": kp,
            "trial": trial,
            "warmup_pitch_deg": math.degrees(base["pitch"]),
            "warmup_ok": warmup_ok,
            "peak_roll_deg": base["roll"] * 180.0 / math.pi,
            "t_peak_s": 0.0,
            "rate_at_cross": None,       # rad/s at first crossing of target roll
            "peak_roll_rate": 0.0,       # max |roll_rate| during SHIFT (inertial energy)
            "lifted": False,
            "max_lift_mm": 0.0,
            "max_leg_torque_nm": 0.0,
            "tipped": False,
            "crossed_tipping_point": False,
            "dynamic_gain_deg": None,
            "final_roll_deg": 0.0,
            "t_fold_s": None,
        }
        if not warmup_ok:
            print(
                f"  [warn] warmup did not stabilize (pitch {m['warmup_pitch_deg']:+.1f}°) — "
                "raise --lqr-warmup-s or adjust --warmup-clearance.",
                flush=True,
            )
        phase = "SHIFT"
        prev_roll = base["roll"]
        n_steps = round(args_cli.flight_s / dt)
        for k in range(n_steps):
            t = (k + 1) * dt

            if phase == "SHIFT":
                frac = min(1.0, t / shift_dur)
                ext_fr_br = args_cli.start_ext_deg + frac * (right_ext_deg - args_cli.start_ext_deg)
                if prev_roll >= fold_roll:
                    phase = "FOLD"
                    m["t_fold_s"] = t
                    fold_start_ext = ext_fr_br
                    fold_start_t = t
            if phase == "FOLD":
                frac = min(1.0, (t - fold_start_t) / args_cli.fold_duration_s)
                ext_fr_br = fold_start_ext + frac * (0.0 - fold_start_ext)

            # LQR keeps holding pitch throughout the thrust so the robot can't
            # fall fore/aft — this isolates the roll-tip question. The optional
            # tip assist adds a forward pulse to the LEFT wheel during the shift.
            eL, eR = lqr_wheel_efforts()
            if phase == "SHIFT" and t <= assist_win:
                eL += assist_nm

            targets = leg_targets(args_cli.start_ext_deg, ext_fr_br)
            step_sim(targets, wheel_effort_L=eL, wheel_effort_R=eR)
            s = read_state()

            if not math.isnan(s["wheel_z"]) and not math.isnan(wheel_z0):
                lift_mm = (s["wheel_z"] - wheel_z0) * 1000.0
                m["max_lift_mm"] = max(m["max_lift_mm"], lift_mm)
                if lift_mm >= args_cli.liftoff_mm:
                    m["lifted"] = True
            elif s["roll"] > math.radians(5.0):
                m["lifted"] = True

            m["max_leg_torque_nm"] = max(m["max_leg_torque_nm"], s["max_cg_torque"])
            if phase == "SHIFT":
                m["peak_roll_rate"] = max(m["peak_roll_rate"], abs(s["roll_rate"]))

            roll_deg = math.degrees(s["roll"])
            pitch_deg = math.degrees(s["pitch"])
            if roll_deg > m["peak_roll_deg"]:
                m["peak_roll_deg"], m["t_peak_s"] = roll_deg, t
            if m["rate_at_cross"] is None and prev_roll < target_roll <= s["roll"]:
                m["rate_at_cross"] = s["roll_rate"]
            if s["roll"] >= tip_roll:
                m["tipped"] = True
            prev_roll = s["roll"]

            if csv_rows is not None:
                csv_rows.append(
                    f"{right_ext_deg},{shift_dur:.3f},{kp:.0f},{trial},{t:.4f},{phase},"
                    f"{roll_deg:.3f},{s['roll_rate']:.4f},{s['body_z']:.4f},{m['max_lift_mm']:.1f},"
                    f"{s['max_cg_torque']:.3f},{ext_fr_br:.2f}"
                )
            if (k + 1) % print_every == 0:
                print(
                    f"  t={t:5.3f}s {phase:5s} roll={roll_deg:+7.2f}° "
                    f" pitch={pitch_deg:+7.2f}° "
                    f"rate={s['roll_rate']:+6.2f}rad/s lift={m['max_lift_mm']:5.1f}mm "
                    f"tq={s['max_cg_torque']:5.1f}Nm FR/BR_ext={ext_fr_br:+5.1f}°",
                    flush=True,
                )
            if m["tipped"]:
                break

        m["final_roll_deg"] = math.degrees(prev_roll)
        m["crossed_tipping_point"] = m["peak_roll_deg"] >= overshoot_target
        if static_ceiling is not None:
            m["dynamic_gain_deg"] = m["peak_roll_deg"] - static_ceiling

        # ── controlled-arrival verdict (original goal: land near target gently) ─
        tgt = args_cli.target_roll_deg
        if not m["lifted"]:
            m["verdict"] = "NO LIFTOFF"
        elif m["tipped"]:
            m["verdict"] = "TIPPED"
        elif m["peak_roll_deg"] < tgt - 5.0:
            m["verdict"] = f"UNDERSHOOT (peak {m['peak_roll_deg']:.1f}°)"
        elif m["peak_roll_deg"] > tgt + 8.0:
            m["verdict"] = f"OVERSHOOT (peak {m['peak_roll_deg']:.1f}°)"
        else:
            rate = m["rate_at_cross"]
            note = f", rate@{tgt:.0f}°={rate:+.2f}" if rate is not None else ""
            m["verdict"] = "GOOD" + note

        # ── feasibility verdict (the go/no-go: can a dynamic thrust cross the bar?) ─
        sat = "  [TORQUE-SATURATED >12Nm — unphysical]" if m["max_leg_torque_nm"] > _CG_PEAK_TORQUE_NM else ""
        if m["crossed_tipping_point"]:
            gain = f" (gain +{m['dynamic_gain_deg']:.1f}°)" if m["dynamic_gain_deg"] is not None else ""
            m["feasibility_verdict"] = f"FEASIBLE — peak {m['peak_roll_deg']:.1f}° >= {overshoot_target:.0f}°{gain}{sat}"
        elif static_ceiling is not None and m["peak_roll_deg"] > static_ceiling + 3.0:
            m["feasibility_verdict"] = (
                f"PARTIAL — dynamic overshoot +{m['dynamic_gain_deg']:.1f}° but short of "
                f"{overshoot_target:.0f}°; try faster shift / higher kp / wheel assist{sat}"
            )
        else:
            m["feasibility_verdict"] = f"STATIC ONLY — no dynamic gain over the ~{static_ceiling:.0f}° ceiling{sat}" \
                if static_ceiling is not None else f"peak {m['peak_roll_deg']:.1f}° (no reference)"
        return m

    # ── sweep setup ───────────────────────────────────────────────────────────
    if args_cli.right_ext_list:
        ext_values = [float(v) for v in args_cli.right_ext_list.split(",")]
    else:
        ext_values = [args_cli.right_ext_deg]
    shift_durs = (
        [float(v) for v in args_cli.shift_duration_list.split(",")]
        if args_cli.shift_duration_list else [args_cli.shift_duration_s]
    )
    kps = (
        [float(v) for v in args_cli.thrust_kp_list.split(",")]
        if args_cli.thrust_kp_list else [args_cli.hold_kp]
    )

    csv_rows = [] if args_cli.csv else None
    STATIC_REF_SHIFT_S = 0.4   # slow quasi-static reference → the static ceiling

    # Static-ceiling reference: one slow, no-assist, hold_kp trial at the largest
    # extension (LQR-stabilized like every trial). Its peak roll is the static
    # ceiling every dynamic trial is compared against.
    _saved_assist = args_cli.wheel_assist_nm
    args_cli.wheel_assist_nm = 0.0
    with torch.inference_mode():
        print("\n── static-ceiling reference (LQR warmup, slow 0.40s shift, no assist) ──", flush=True)
        ref = run_trial(max(ext_values), STATIC_REF_SHIFT_S, args_cli.hold_kp, 0, None, None)
    args_cli.wheel_assist_nm = _saved_assist
    static_ceiling = ref["peak_roll_deg"]
    print(
        f"[jump-diag] STATIC CEILING ≈ {static_ceiling:.1f}° (peak roll of the slow reference; "
        f"warmup pitch {ref['warmup_pitch_deg']:+.1f}°)",
        flush=True,
    )

    print()
    print("=" * 92)
    print("ONE-LEG DYNAMIC TIP FEASIBILITY DIAGNOSTIC (LQR-stabilized warmup, no policy)")
    print(f"  target {args_cli.target_roll_deg:.0f}°  overshoot bar {args_cli.overshoot_target_deg:.0f}°  "
          f"start ext {args_cli.start_ext_deg:.0f}°  right-ext {ext_values}°  "
          f"LQR warmup {args_cli.lqr_warmup_s:.1f}s")
    print(f"  shift durations {shift_durs}s  thrust kp {kps} Nm/rad  "
          f"wheel-assist {args_cli.wheel_assist_nm * args_cli.wheel_assist_sign:+.2f}Nm  "
          f"static ceiling ≈ {static_ceiling:.1f}°")
    print("=" * 92)

    results = []
    with torch.inference_mode():
        for ext in ext_values:
            for dur in shift_durs:
                for kp in kps:
                    for trial in range(args_cli.trials):
                        print(
                            f"\n── ext {ext:.0f}°, shift {dur:.3f}s, kp {kp:.0f}, "
                            f"trial {trial + 1}/{args_cli.trials} ──",
                            flush=True,
                        )
                        results.append(run_trial(ext, dur, kp, trial, static_ceiling, csv_rows))

    # ── summary table ─────────────────────────────────────────────────────────
    print()
    print("=" * 92)
    print(f"{'ext°':>4} {'shift':>6} {'kp':>4} {'trl':>3} {'peak°':>6} {'gain°':>6} "
          f"{'pkRate':>6} {'lift':>5} {'tq':>5}  feasibility")
    print("-" * 92)
    for m in results:
        gain = f"{m['dynamic_gain_deg']:+.1f}" if m["dynamic_gain_deg"] is not None else " n/a"
        print(
            f"{m['right_ext_deg']:4.0f} {m['shift_dur']:6.3f} {m['kp']:4.0f} {m['trial'] + 1:3d} "
            f"{m['peak_roll_deg']:6.1f} {gain:>6} {m['peak_roll_rate']:6.2f} "
            f"{m['max_lift_mm']:5.1f} {m['max_leg_torque_nm']:5.1f}  {m['feasibility_verdict']}"
        )
    print("=" * 92)

    # ── overall go/no-go ──────────────────────────────────────────────────────
    # A setting is FEASIBLE only if it crosses the bar on ALL its trials with
    # physical torque (repeatable, not a fluke).
    from collections import defaultdict
    by_setting: dict = defaultdict(list)
    for m in results:
        by_setting[(m["right_ext_deg"], m["shift_dur"], m["kp"])].append(m)
    feasible_settings = [
        key for key, ms in by_setting.items()
        if all(x["crossed_tipping_point"] for x in ms)
        and all(x["max_leg_torque_nm"] <= _CG_PEAK_TORQUE_NM for x in ms)
    ]
    print()
    if feasible_settings:
        # minimal = slowest shift / lowest kp that still works (the "target aggressiveness")
        minimal = sorted(feasible_settings, key=lambda k: (-k[1], k[2]))[0]
        print(f"VERDICT: FEASIBLE — dynamic tip crosses {args_cli.overshoot_target_deg:.0f}° repeatably.")
        print(f"  {len(feasible_settings)} setting(s) work. Minimal (slowest/softest): "
              f"ext {minimal[0]:.0f}°, shift {minimal[1]:.3f}s, kp {minimal[2]:.0f} Nm/rad.")
        print("  → Proceed to Part B (learning). This is the target thrust aggressiveness.")
    else:
        best = max(results, key=lambda x: x["peak_roll_deg"])
        print(f"VERDICT: NOT FEASIBLE with the swept settings. Best peak roll = "
              f"{best['peak_roll_deg']:.1f}° (bar {args_cli.overshoot_target_deg:.0f}°, "
              f"static ceiling {static_ceiling:.1f}°).")
        print("  → Try faster --shift-duration-list, higher --thrust-kp-list, or --wheel-assist-nm 2.0 "
              "before concluding the tip is physically impossible.")
    print("=" * 92)

    if args_cli.csv and csv_rows is not None:
        header = ("right_ext_deg,shift_dur_s,kp,trial,t_s,phase,roll_deg,roll_rate_rad_s,"
                  "body_z_m,lift_mm,max_cg_torque_nm,fr_br_ext_deg")
        Path(args_cli.csv).parent.mkdir(parents=True, exist_ok=True)
        Path(args_cli.csv).write_text(header + "\n" + "\n".join(csv_rows) + "\n")
        print(f"[jump-diag] per-step log written to {args_cli.csv}")

    env.close()


if __name__ == "__main__":
    main()
    simulation_app.close()
