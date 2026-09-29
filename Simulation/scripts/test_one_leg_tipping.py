"""
One-leg static tipping angle diagnostic.

Spawns the robot at a fixed roll angle with bottom CyberGear joints held at a
fixed angle and top joints locked at their upper limit. No control is applied —
gravity does the work. Watch the GUI and read the terminal to see whether the
robot settles or tips over.

Usage (run from repo root):

    # Test roll=45 deg with bottom joints at 45 deg
    /isaac-lab/isaaclab.sh -p scripts/test_one_leg_tipping.py \\
        --task Template-Twowheeledrobot-Standup-v0 \\
        --roll-deg 45.0 \\
        --bottom-joint-deg 45.0

    # Binary search: try different roll angles
    /isaac-lab/isaaclab.sh -p scripts/test_one_leg_tipping.py \\
        --task Template-Twowheeledrobot-Standup-v0 \\
        --roll-deg 55.0 \\
        --bottom-joint-deg 45.0

Terminal output every second:
    t=1.00s  roll=+52.3°  body_z=0.082m  FL=45.0°  BL=45.0°  status=HOLDING
    t=2.00s  roll=+53.1°  body_z=0.081m  FL=45.0°  BL=45.0°  status=HOLDING
    ...

status=HOLDING  → robot is close to the spawn angle, gravity balanced
status=FALLING  → roll is increasing fast, robot is tipping over
status=SETTLING → roll is decreasing, robot found stable equilibrium

Stop with Ctrl+C. Record the roll angle at which HOLDING transitions to FALLING
— that is your tipping angle for this joint configuration.
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

parser = argparse.ArgumentParser(description="One-leg static tipping angle diagnostic.")
parser.add_argument("--task", type=str, default="Template-Twowheeledrobot-Standup-v0")
parser.add_argument(
    "--roll-deg",
    type=float,
    default=45.0,
    help="Initial roll angle in degrees. Robot is spawned at this angle and released.",
)
parser.add_argument(
    "--bottom-joint-deg",
    type=float,
    default=45.0,
    help="Target angle in degrees for both grounded-side (left) CyberGear joints.",
)
parser.add_argument(
    "--top-joint-deg",
    type=float,
    default=0.0,
    help="Lock angle in degrees for top-side (right) CyberGear joints. Default is near upper limit.",
)
parser.add_argument(
    "--grounded-side",
    choices=("left", "right"),
    default="left",
    help="Which side is on the ground. Left = front_left + back_left active, right = front_right + back_right active.",
)
parser.add_argument(
    "--max-time-s",
    type=float,
    default=10.0,
    help="Stop simulation after this many seconds.",
)
parser.add_argument(
    "--print-interval-s",
    type=float,
    default=0.5,
    help="Print state to terminal every this many seconds.",
)
parser.add_argument(
    "--fall-roll-threshold-deg",
    type=float,
    default=80.0,
    help="Declare TIPPED_OVER when |roll| exceeds this value.",
)
AppLauncher.add_app_launcher_args(parser)
args_cli, hydra_args = parser.parse_known_args()
args_cli.num_envs = 1
sys.argv = [sys.argv[0]] + hydra_args

app_launcher = AppLauncher(args_cli)
simulation_app = app_launcher.app

# ─── everything below runs after Isaac Sim has launched ──────────────────────

import gymnasium as gym
import torch
import TwoWheeledRobot.tasks  # noqa: F401
import isaaclab_tasks  # noqa: F401
from isaaclab_tasks.utils.hydra import hydra_task_config
from isaaclab.envs import DirectRLEnvCfg, ManagerBasedRLEnvCfg, DirectMARLEnvCfg


def roll_from_projected_gravity(projected_gravity_b: torch.Tensor) -> torch.Tensor:
    """Roll angle around the forward (Y) axis from body-frame projected gravity.

    Positive roll = robot tips to the right (positive X side).
    Sign convention mirrors pitch_from_projected_gravity in pure_nn_components.py
    which uses atan2(grav_y, -grav_z) for pitch.
    """
    return torch.atan2(projected_gravity_b[:, 0], -projected_gravity_b[:, 2])


def build_roll_quaternion(roll_rad: float) -> torch.Tensor:
    """Return a wxyz quaternion for a pure roll rotation around the body Y axis.

    This is applied on top of the default upright spawn orientation.
    Positive roll_rad tilts the robot to the right (left leg goes down).
    """
    half = roll_rad * 0.5
    # Rotation around world Y axis: q = [cos(θ/2), 0, sin(θ/2), 0] in wxyz
    return torch.tensor([math.cos(half), 0.0, math.sin(half), 0.0], dtype=torch.float32)


def joint_angle_to_action_fraction(angle_rad: float, lo_rad: float, hi_rad: float) -> float:
    """Convert a physical joint angle to the [-1, 1] action fraction.

    Matches the formula in lqr_control.py::action_to_env():
        action = 2 * (target - lo) / (hi - lo) - 1
    """
    span = hi_rad - lo_rad
    if abs(span) < 1e-9:
        return 0.0
    return float(2.0 * (angle_rad - lo_rad) / span - 1.0)


@hydra_task_config(args_cli.task, "rsl_rl_cfg_entry_point")
def main(env_cfg: ManagerBasedRLEnvCfg | DirectRLEnvCfg | DirectMARLEnvCfg, agent_cfg=None):
    env_cfg.scene.num_envs = 1
    env_cfg.seed = 0

    env = gym.make(args_cli.task, cfg=env_cfg)
    env.reset()

    # ── resolve joint index order ─────────────────────────────────────────────
    # StandupEnv stores _cg_ids as [front_left, front_right, back_left, back_right]
    # (the order they appear in robot_cfg.py joint_names_expr).
    # Within the 6D action: indices 0-3 are CyberGear, 4-5 are wheel currents.
    # CyberGear action order matches _cg_ids order.
    inner = env.unwrapped
    cg_lo = inner._cg_joint_lo[0].cpu()   # shape [4]
    cg_hi = inner._cg_joint_hi[0].cpu()   # shape [4]

    # Index mapping within the 4 CyberGear joints:
    # 0=front_left  1=front_right  2=back_left  3=back_right
    IDX_FRONT_LEFT  = 0
    IDX_FRONT_RIGHT = 1
    IDX_BACK_LEFT   = 2
    IDX_BACK_RIGHT  = 3

    if args_cli.grounded_side == "left":
        active_ids  = [IDX_FRONT_LEFT, IDX_BACK_LEFT]
        locked_ids  = [IDX_FRONT_RIGHT, IDX_BACK_RIGHT]
        bottom_angle_rad = math.radians(args_cli.bottom_joint_deg)
        top_angle_rad    = math.radians(args_cli.top_joint_deg)
    else:
        active_ids  = [IDX_FRONT_RIGHT, IDX_BACK_RIGHT]
        locked_ids  = [IDX_FRONT_LEFT, IDX_BACK_LEFT]
        # Right joints have inverted limits — clamp target to their valid range
        bottom_angle_rad = math.radians(-args_cli.bottom_joint_deg)
        top_angle_rad    = math.radians(-args_cli.top_joint_deg)

    # Build the fixed 6D action tensor [cg0, cg1, cg2, cg3, wheel_L, wheel_R]
    # _cg_ext_sign = [+1, -1, -1, +1] for [front_left, front_right, back_left, back_right].
    # Mirrored joints (back_left, front_right) have inverted raw limits so
    # positive extension maps to a negative raw angle. Apply the sign before
    # converting to an action fraction so "45 deg" means the same physical
    # extension direction on every joint.
    ext_sign = inner._cg_ext_sign[0].cpu()   # shape [4]

    cg_action = torch.zeros(1, 4, device=inner.device)
    for idx in active_ids:
        lo = float(cg_lo[idx])
        hi = float(cg_hi[idx])
        raw_angle = float(ext_sign[idx]) * bottom_angle_rad
        cg_action[0, idx] = joint_angle_to_action_fraction(raw_angle, lo, hi)
    for idx in locked_ids:
        lo = float(cg_lo[idx])
        hi = float(cg_hi[idx])
        raw_angle = float(ext_sign[idx]) * top_angle_rad
        cg_action[0, idx] = joint_angle_to_action_fraction(raw_angle, lo, hi)

    # Wheel currents are always zero — passive observation
    wheel_action = torch.zeros(1, 2, device=inner.device)
    action = torch.cat([cg_action, wheel_action], dim=1)

    # ── spawn robot at the requested roll angle ───────────────────────────────
    roll_rad = math.radians(args_cli.roll_deg)
    roll_quat = build_roll_quaternion(roll_rad).to(inner.device)

    # Compute spawn height matching standup_env_cfg.py formula:
    #   z = upright_z * |cos(roll)| + lateral_h * |sin(roll)| + clearance
    # Without this the robot spawns underground at non-zero roll and PhysX
    # ejects it upward.
    spawn_z = (
        inner.cfg.spawn_upright_z  * abs(math.cos(roll_rad))
        + inner.cfg.spawn_lateral_h * abs(math.sin(roll_rad))
        + inner.cfg.spawn_clearance
    )

    root_state = inner.robot.data.default_root_state[0:1].clone()
    root_state[:, :3] += inner.scene.env_origins[0:1]
    root_state[:, 2]   = inner.scene.env_origins[0, 2] + spawn_z
    root_state[:, 3:7] = roll_quat
    root_state[:, 7:]  = 0.0
    inner.robot.write_root_pose_to_sim(root_state[:, :7], torch.tensor([0], device=inner.device))
    inner.robot.write_root_velocity_to_sim(root_state[:, 7:], torch.tensor([0], device=inner.device))

    # ── print configuration summary ───────────────────────────────────────────
    print()
    print("=" * 60)
    print("ONE-LEG TIPPING DIAGNOSTIC")
    print("=" * 60)
    print(f"  Spawn roll angle  : {args_cli.roll_deg:+.1f} deg")
    print(f"  Bottom joints     : {args_cli.grounded_side} side at {args_cli.bottom_joint_deg:.1f} deg")
    print(f"  Top joints locked : opposite side at {args_cli.top_joint_deg:.1f} deg")
    print(f"  Wheel current     : 0 A (passive)")
    print(f"  Max test time     : {args_cli.max_time_s:.1f} s")
    print(f"  Tipped-over limit : |roll| > {args_cli.fall_roll_threshold_deg:.1f} deg")
    print()
    print(f"  CG action vector  : {action[0].cpu().tolist()}")
    print()
    print("  HOLDING  = roll stable near spawn angle")
    print("  SETTLING = roll moving toward a lower equilibrium")
    print("  FALLING  = roll accelerating away — tipping over")
    print("  TIPPED   = |roll| exceeded fall threshold")
    print("=" * 60)
    print()

    # ── main simulation loop ──────────────────────────────────────────────────
    control_dt     = inner.step_dt
    print_every_n  = max(1, round(args_cli.print_interval_s / control_dt))
    max_steps      = round(args_cli.max_time_s / control_dt)
    fall_threshold = math.radians(args_cli.fall_roll_threshold_deg)

    prev_roll_rad  = roll_rad
    tipped         = False
    step           = 0

    with torch.inference_mode():
        while simulation_app.is_running() and step < max_steps and not tipped:
            obs, _, terminated, truncated, _ = env.step(action)
            step += 1
            time_s = step * control_dt

            # ── read current state ────────────────────────────────────────────
            grav_b   = inner.bno080.data.projected_gravity_b
            roll_now = float(roll_from_projected_gravity(grav_b)[0])
            body_z   = float(inner.robot.data.root_pos_w[0, 2])
            roll_rate = (roll_now - prev_roll_rad) / control_dt
            prev_roll_rad = roll_now

            # ── read actual joint positions (converted to extension degrees) ───
            joint_pos = inner.robot.data.joint_pos[0, inner._cg_ids].cpu()
            fl_deg = math.degrees(float(ext_sign[IDX_FRONT_LEFT])  * float(joint_pos[IDX_FRONT_LEFT]))
            fr_deg = math.degrees(float(ext_sign[IDX_FRONT_RIGHT]) * float(joint_pos[IDX_FRONT_RIGHT]))
            bl_deg = math.degrees(float(ext_sign[IDX_BACK_LEFT])   * float(joint_pos[IDX_BACK_LEFT]))
            br_deg = math.degrees(float(ext_sign[IDX_BACK_RIGHT])  * float(joint_pos[IDX_BACK_RIGHT]))

            # ── classify status ───────────────────────────────────────────────
            roll_deg_now = math.degrees(roll_now)
            if abs(roll_now) >= fall_threshold:
                status = "TIPPED"
                tipped = True
            elif abs(roll_rate) > math.radians(30):       # fast angular motion
                if abs(roll_now) > abs(roll_rad) * 1.05:  # moving away from spawn
                    status = "FALLING"
                else:
                    status = "SETTLING"
            else:
                status = "HOLDING"

            # ── periodic terminal print ───────────────────────────────────────
            if step % print_every_n == 0 or tipped:
                print(
                    f"t={time_s:5.2f}s  "
                    f"roll={roll_deg_now:+6.2f}°  "
                    f"body_z={body_z:.4f}m  "
                    f"FL={fl_deg:+5.1f}° BL={bl_deg:+5.1f}°  "
                    f"FR={fr_deg:+5.1f}° BR={br_deg:+5.1f}°  "
                    f"status={status}  [ext deg]"
                )

    # ── final verdict ─────────────────────────────────────────────────────────
    print()
    print("=" * 60)
    if tipped:
        print(f"RESULT: TIPPED OVER at roll={math.degrees(prev_roll_rad):+.1f}°")
        print(f"        Try a smaller --roll-deg value.")
    else:
        print(f"RESULT: HELD for {args_cli.max_time_s:.1f}s, final roll={math.degrees(prev_roll_rad):+.1f}°")
        print(f"        Robot appears stable at this angle. Try a larger --roll-deg value.")
    print()
    print(f"  Spawn angle    : {args_cli.roll_deg:+.1f} deg")
    print(f"  Final roll     : {math.degrees(prev_roll_rad):+.1f} deg")
    print(f"  Bottom joints  : {args_cli.bottom_joint_deg:.1f} deg extension"
          f" (front_left raw={fl_deg:.1f}°  back_left raw={bl_deg:.1f}°)")
    print("=" * 60)

    env.close()


if __name__ == "__main__":
    main()
    simulation_app.close()
