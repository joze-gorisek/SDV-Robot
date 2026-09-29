"""
Robot mass/inertia/joint-property diagnostic.

Purpose:
    Dump per-body mass + inertia and per-joint PhysX limits (position, velocity,
    effort, friction, armature) plus the configured actuator gains (kp/kd), all
    parsed straight from the USD at spawn time. No control loop is run and no
    physics steps happen beyond env.reset() -- this only reads static/default
    properties.

    Written to chase down a "sim needs way more torque than real hardware"
    discrepancy (e.g. lqr_control_one_leg_jump.py's one-leg thrust saturating
    the 12 Nm CyberGear limit just to hold a pose the real robot holds fine).
    The usual causes are: a CAD-imported body carrying an incorrect/placeholder
    mass or inertia, a suspiciously low PhysX velocity limit silently braking a
    joint regardless of applied torque, or nonzero friction/armature values that
    were never set anywhere in robot_cfg.py (i.e. they came straight from the
    USD import, unreviewed). This script flags all three automatically.

Usage (inside the Isaac Lab container, from /workspace):

    isaaclab/isaaclab.sh -p Project/scripts/inspect_robot_mass_properties.py \\
        --task Template-Twowheeledrobot-Standup-v0 --headless

    # Compare sim total mass against a scale reading of the real robot:
    isaaclab/isaaclab.sh -p Project/scripts/inspect_robot_mass_properties.py \\
        --task Template-Twowheeledrobot-Standup-v0 --headless --real-mass-kg 3.9

Edit here when:
    You need to add another PhysX/USD-imported property to the dump.

Avoid changing here without also checking:
    robot_cfg.py for the actuator groups this cross-checks against, and
    lqr_control_one_leg_jump.py (LQR_PHYSICAL_PARAMS) for the design-time mass
    estimate this compares the sim total against.
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

parser = argparse.ArgumentParser(description="Dump per-body mass/inertia and per-joint PhysX properties.")
parser.add_argument("--task", type=str, default="Template-Twowheeledrobot-Standup-v0")
parser.add_argument(
    "--real-mass-kg", type=float, default=None,
    help="Measured real-robot total mass (kg), e.g. from a scale. If given, compares against the sim total.",
)
parser.add_argument(
    "--mass-flag-pct", type=float, default=30.0,
    help="Flag any single body carrying more than this percent of total sim mass.",
)
parser.add_argument(
    "--vel-limit-flag", type=float, default=10.0,
    help="Flag any joint whose PhysX velocity limit (rad/s) is below this -- low enough to brake the "
         "joint regardless of applied torque.",
)
AppLauncher.add_app_launcher_args(parser)
args_cli, hydra_args = parser.parse_known_args()
args_cli.num_envs = 1
sys.argv = [sys.argv[0]] + hydra_args

app_launcher = AppLauncher(args_cli)
simulation_app = app_launcher.app

# ─── everything below runs after Isaac Sim has launched ──────────────────────

import gymnasium as gym
import TwoWheeledRobot.tasks  # noqa: F401
import isaaclab_tasks  # noqa: F401
from isaaclab.envs import DirectMARLEnvCfg, DirectRLEnvCfg, ManagerBasedRLEnvCfg
from isaaclab_tasks.utils.hydra import hydra_task_config

# Design-time mass estimate the LQR gains in lqr_control_one_leg_jump.py were
# derived against -- NOT read from the USD, just a cross-reference. A large gap
# between this and the sim's actual USD-parsed total below is a real red flag.
LQR_REFERENCE_BODY_MASS_KG = 2.6
LQR_REFERENCE_WHEEL_CART_MASS_KG = 1.53


@hydra_task_config(args_cli.task, "rsl_rl_cfg_entry_point")
def main(env_cfg: ManagerBasedRLEnvCfg | DirectRLEnvCfg | DirectMARLEnvCfg, agent_cfg=None):
    env_cfg.scene.num_envs = 1
    env_cfg.seed = 0

    env = gym.make(args_cli.task, cfg=env_cfg)
    env.reset()
    robot = env.unwrapped.robot
    data = robot.data

    # ── body mass / inertia ─────────────────────────────────────────────────
    masses = data.default_mass[0]            # (num_bodies,) kg
    inertias = data.default_inertia[0]        # (num_bodies, 9) [Ixx,Iyx,Izx,Ixy,Iyy,Izy,Ixz,Iyz,Izz]
    total_mass = float(masses.sum())
    order = sorted(range(len(robot.body_names)), key=lambda i: float(masses[i]), reverse=True)

    print()
    print("=" * 92)
    print("BODY MASS / INERTIA  (sorted heaviest first, from the USD at spawn time)")
    print("=" * 92)
    print(f"{'body_name':30s} {'mass_kg':>9s} {'% total':>8s} {'Ixx':>10s} {'Iyy':>10s} {'Izz':>10s}")
    print("-" * 92)
    mass_flags = []
    for i in order:
        name = robot.body_names[i]
        m = float(masses[i])
        pct = 100.0 * m / total_mass if total_mass > 0 else 0.0
        ixx, iyy, izz = float(inertias[i, 0]), float(inertias[i, 4]), float(inertias[i, 8])
        flag = ""
        if pct >= args_cli.mass_flag_pct:
            flag = f"  <== {pct:.0f}% of total mass"
            mass_flags.append(name)
        print(f"{name:30s} {m:9.4f} {pct:7.1f}% {ixx:10.6f} {iyy:10.6f} {izz:10.6f}{flag}")
    print("-" * 92)
    print(f"TOTAL MASS (sim, all bodies): {total_mass:.4f} kg")
    print(
        f"Reference only (NOT read from the USD): lqr_control_one_leg_jump.py's LQR_PHYSICAL_PARAMS assumes "
        f"body_mass_kg={LQR_REFERENCE_BODY_MASS_KG:.2f} + wheel_cart_mass_kg={LQR_REFERENCE_WHEEL_CART_MASS_KG:.2f} "
        f"= {LQR_REFERENCE_BODY_MASS_KG + LQR_REFERENCE_WHEEL_CART_MASS_KG:.2f} kg (design-time estimate)."
    )
    if args_cli.real_mass_kg is not None:
        diff = total_mass - args_cli.real_mass_kg
        pct_diff = 100.0 * diff / args_cli.real_mass_kg if args_cli.real_mass_kg > 0 else 0.0
        print(
            f"Measured real robot mass: {args_cli.real_mass_kg:.4f} kg  |  sim - real = {diff:+.4f} kg "
            f"({pct_diff:+.1f}%)"
        )
        if abs(pct_diff) >= 15.0:
            print(
                f"  !! FLAG: sim total mass differs from the measured real mass by {pct_diff:+.1f}% -- "
                "this alone could explain a sim-only torque/holding deficit."
            )
    else:
        print("(pass --real-mass-kg <scale reading> to auto-compare against the real robot)")
    if mass_flags:
        print(f"  !! FLAG: {', '.join(mass_flags)} each carry >= {args_cli.mass_flag_pct:.0f}% of total sim mass.")

    # ── joint properties: PhysX limits + configured actuator gains ─────────
    joint_names = robot.joint_names
    num_joints = len(joint_names)
    kp_by_joint: list[float | None] = [None] * num_joints
    kd_by_joint: list[float | None] = [None] * num_joints
    actuator_by_joint: list[str] = ["?"] * num_joints
    for act_name, actuator in robot.actuators.items():
        for col, jid in enumerate(actuator.joint_indices.tolist()):
            actuator_by_joint[jid] = act_name
            if hasattr(actuator, "stiffness"):
                kp_by_joint[jid] = float(actuator.stiffness[0, col])
                kd_by_joint[jid] = float(actuator.damping[0, col])

    print()
    print("=" * 110)
    print("JOINT PROPERTIES  (PhysX solver-enforced limits + configured actuator gains)")
    print("=" * 110)
    print(
        f"{'joint_name':24s} {'actuator':17s} {'lo_deg':>8s} {'hi_deg':>8s} {'vel_lim':>9s} "
        f"{'eff_lim':>8s} {'friction':>9s} {'armature':>9s} {'kp':>7s} {'kd':>7s}"
    )
    print("-" * 110)
    vel_flags, friction_flags, armature_flags, unlimited_flags = [], [], [], []
    # PhysX reports an unconfigured joint limit as +/-FLT_MAX (~3.4028235e38 rad),
    # not zero or NaN -- converted to degrees that prints as a ~40-digit number.
    # Treat anything past this as "no limit configured" rather than a real angle.
    _UNLIMITED_RAD = 1.0e30
    for j in range(num_joints):
        name = joint_names[j]
        lo_rad = float(data.joint_pos_limits[0, j, 0])
        hi_rad = float(data.joint_pos_limits[0, j, 1])
        lo_deg = "unlimited" if lo_rad <= -_UNLIMITED_RAD else f"{math.degrees(lo_rad):.1f}"
        hi_deg = "unlimited" if hi_rad >= _UNLIMITED_RAD else f"{math.degrees(hi_rad):.1f}"
        vel_lim = float(data.joint_vel_limits[0, j])
        eff_lim = float(data.joint_effort_limits[0, j])
        friction = float(data.joint_friction_coeff[0, j])
        armature = float(data.joint_armature[0, j])
        kp = f"{kp_by_joint[j]:7.2f}" if kp_by_joint[j] is not None else "      -"
        kd = f"{kd_by_joint[j]:7.2f}" if kd_by_joint[j] is not None else "      -"
        notes = []
        if actuator_by_joint[j] in ("cybergear_joints", "bearing_joints") and lo_deg == "unlimited" and hi_deg == "unlimited":
            notes.append("NO POSITION LIMIT IN PHYSX (only the env's software target-clamp restrains it)")
            unlimited_flags.append(name)
        if vel_lim < args_cli.vel_limit_flag:
            notes.append("LOW VEL_LIM")
            vel_flags.append(name)
        if friction > 1.0e-6:
            notes.append("NONZERO FRICTION (not set in robot_cfg.py)")
            friction_flags.append(name)
        if armature > 1.0e-6:
            notes.append("NONZERO ARMATURE (not set in robot_cfg.py)")
            armature_flags.append(name)
        note_str = f"  <== {', '.join(notes)}" if notes else ""
        print(
            f"{name:24s} {actuator_by_joint[j]:17s} {lo_deg:>8s} {hi_deg:>8s} {vel_lim:9.3f} "
            f"{eff_lim:8.3f} {friction:9.4f} {armature:9.5f} {kp} {kd}{note_str}"
        )
    print("-" * 110)
    if unlimited_flags:
        print(
            f"  !! FLAG: {', '.join(unlimited_flags)} have NO PhysX position limit configured -- PhysX will let "
            "these rotate arbitrarily far under enough torque/momentum; only the env's software clamp on the "
            "TARGET (not the actual joint) keeps them near the intended range during normal position control."
        )
    if vel_flags:
        print(
            f"  !! FLAG: {', '.join(vel_flags)} have a PhysX velocity limit below {args_cli.vel_limit_flag:.1f} "
            "rad/s -- PhysX will brake the joint to this speed no matter how much torque is applied."
        )
    if friction_flags:
        print(
            f"  !! FLAG: {', '.join(friction_flags)} have nonzero joint friction coming from the USD import "
            "(robot_cfg.py never sets joint_friction_coeff for any actuator group)."
        )
    if armature_flags:
        print(
            f"  !! FLAG: {', '.join(armature_flags)} have nonzero armature coming from the USD import "
            "(robot_cfg.py never sets armature for any actuator group) -- armature adds artificial "
            "reflected inertia at the joint, which directly slows angular acceleration under a fixed torque."
        )
    if not (vel_flags or friction_flags or armature_flags or unlimited_flags or mass_flags):
        print("No automatic flags raised -- mass/inertia and joint properties look unremarkable by these heuristics.")
    print("=" * 110)

    env.close()


if __name__ == "__main__":
    main()  # type: ignore[reportCallIssue]
    simulation_app.close()
