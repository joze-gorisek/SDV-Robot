# Copyright (c) 2022-2026, The Isaac Lab Project Developers (https://github.com/isaac-sim/IsaacLab/blob/main/CONTRIBUTORS.md).
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Play and export RSL-RL checkpoints for registered Isaac Lab tasks.

Purpose:
    Launch Isaac Sim, load a trained checkpoint, step the selected environment,
    and export the policy to TorchScript and ONNX for evaluation/deployment.

Edit here when:
    You need to adjust checkpoint resolution, export behavior, or play-loop
    diagnostics. Most task-specific behavior belongs in environment/config files.

Avoid changing here without also checking:
    scripts/rsl_rl/train.py, scripts/rsl_rl/cli_args.py, UART deployment, and
    registered task IDs in twowheeledrobot/__init__.py.
"""

"""Launch Isaac Sim Simulator first."""

import argparse
import sys
from pathlib import Path

_EXTENSION_SOURCE_PATH = Path(__file__).resolve().parents[2] / "source" / "TwoWheeledRobot"
if _EXTENSION_SOURCE_PATH.is_dir():
    sys.path.insert(0, str(_EXTENSION_SOURCE_PATH))

from isaaclab.app import AppLauncher

# local imports
import cli_args  # isort: skip

# add argparse arguments
parser = argparse.ArgumentParser(description="Train an RL agent with RSL-RL.")
parser.add_argument("--video", action="store_true", default=False, help="Record videos during training.")
parser.add_argument("--video_length", type=int, default=200, help="Length of the recorded video (in steps).")
parser.add_argument(
    "--disable_fabric", action="store_true", default=False, help="Disable fabric and use USD I/O operations."
)
parser.add_argument(
    "--cinematic",
    action="store_true",
    default=False,
    help="Presentation-quality scene: robot-tracking camera, key light, matte ground. For video capture only.",
)
parser.add_argument(
    "--cinematic_follow_yaw",
    action="store_true",
    default=False,
    help=(
        "Camera orbits with the robot's heading so the viewing angle stays fixed relative to the robot."
        " Horizon stays level (yaw only) — the balance lean is preserved on screen. Requires --cinematic."
    ),
)
parser.add_argument(
    "--cinematic_orbit",
    type=float,
    default=None,
    metavar="TURNS",
    help=(
        "Sweep the camera TURNS full revolutions around the robot over the capture (e.g. 2)."
        " Overrides --cinematic_follow_yaw. Requires --cinematic."
    ),
)
parser.add_argument(
    "--disturbance",
    type=str,
    default=None,
    choices=["none", "human_push", "double_human_push", "payload", "payload_push", "slope", "sine_diagnostic"],
    help=(
        "Force an external disturbance every episode, bypassing the curriculum stage gate"
        " (disturbances are otherwise off below stage 3). 'double_human_push' gives two force pulses per episode."
    ),
)
parser.add_argument(
    "--episode_length",
    type=float,
    default=None,
    help=(
        "Override episode_length_s (seconds). Raise it for uninterrupted video takes — the default 10 s"
        " timeout respawns the robot mid-shot. Does not affect early terminations (falls)."
    ),
)
parser.add_argument("--num_envs", type=int, default=1, help="Number of environments to simulate.")
parser.add_argument("--task", type=str, default=None, help="Name of the task.")
parser.add_argument(
    "--agent", type=str, default="rsl_rl_cfg_entry_point", help="Name of the RL agent configuration entry point."
)
parser.add_argument("--seed", type=int, default=None, help="Seed used for the environment")
parser.add_argument(
    "--use_pretrained_checkpoint",
    action="store_true",
    help="Use the pre-trained checkpoint from Nucleus.",
)
parser.add_argument("--real-time", action="store_true", default=False, help="Run in real-time, if possible.")
parser.add_argument("--num_steps", type=int, default=0, help="Stop after this many steps (0 = run forever).")
parser.add_argument(
    "--spawn_roll_deg",
    type=float,
    default=None,
    help="Exact spawn roll angle in degrees (one-leg balance task only). "
    "Overrides the reset roll center and disables the random perturbation, "
    "so every env spawns at exactly this angle.",
)
parser.add_argument(
    "--measure_survival",
    action="store_true",
    default=False,
    help="Measure the per-episode survival rate at the spawn angle. Counts how "
    "many episodes reach the full episode length (survival) vs terminate early "
    "(fall), then prints the survival rate and mean episode length. Use together "
    "with --spawn_roll_deg, --num_envs and --num_steps.",
)
parser.add_argument(
    "--spawn_roll_deg_list",
    type=str,
    default=None,
    help="Comma-separated spawn roll angles (deg), e.g. '44,46,47,48'. Sweeps all of "
    "them in ONE Isaac Sim process/app launch instead of one process per angle -- "
    "forces every env to respawn at each angle in turn (env.reset()) and runs "
    "--num_steps control steps per angle, printing a [SURVIVAL] block for each "
    "(requires --measure_survival and --num_steps > 0; overrides --spawn_roll_deg).",
)
# append RSL-RL cli arguments
cli_args.add_rsl_rl_args(parser)
# append AppLauncher cli args
AppLauncher.add_app_launcher_args(parser)
# parse the arguments
args_cli, hydra_args = parser.parse_known_args()
# always enable cameras to record video
if args_cli.video:
    args_cli.enable_cameras = True

# clear out sys.argv for Hydra
sys.argv = [sys.argv[0]] + hydra_args

# launch omniverse app
app_launcher = AppLauncher(args_cli)
simulation_app = app_launcher.app

"""Rest everything follows."""

import math
import os
import time

import gymnasium as gym
import torch

# Work around a PyTorch/inspect incompatibility with namespace packages
# (isaaclab has no __init__.py). torch.library.register_fake's decorator --
# run at import time of torch.distributed.tensor._collective_utils, which
# rsl_rl.runners pulls in transitively -- inspects the calling frame purely
# to build a cosmetic label for the op registration. If the resolved frame's
# module happens to be the isaaclab namespace package, inspect.getfile()
# raises TypeError("... is a built-in module") because namespace packages
# have no conventional __file__. This has no effect on actual op
# registration/dispatch (it's diagnostic bookkeeping only), so short-circuit
# it rather than crash.
#
# Two separate call sites need covering: torch._library.utils.get_source()
# (via inspect.getframeinfo) and Library._register_fake's direct
# inspect.getmodule(frame) call. Best-effort: if the private API doesn't
# exist in the installed torch version, just proceed unpatched.
try:
    import torch._library.utils as _torch_library_utils

    _torch_library_utils.get_source = lambda *args, **kwargs: ""
except Exception:
    pass

# _register_fake tolerates caller_module being None, so swallowing the
# TypeError and reporting "no module" is safe.
import inspect as _inspect  # noqa: E402

_orig_getmodule = _inspect.getmodule


def _safe_getmodule(obj, *args, **kwargs):
    try:
        return _orig_getmodule(obj, *args, **kwargs)
    except TypeError:
        return None


_inspect.getmodule = _safe_getmodule

from rsl_rl.runners import OnPolicyRunner

from isaaclab.envs import (
    DirectMARLEnv,
    DirectMARLEnvCfg,
    DirectRLEnvCfg,
    ManagerBasedRLEnvCfg,
    ViewerCfg,
    multi_agent_to_single_agent,
)
from isaaclab.utils.assets import retrieve_file_path
from isaaclab.utils.dict import print_dict

from isaaclab_rl.rsl_rl import (
    RslRlBaseRunnerCfg,
    RslRlVecEnvWrapper,
    export_policy_as_jit,
    export_policy_as_onnx,
    handle_deprecated_rsl_rl_cfg,
)
from isaaclab_rl.utils.pretrained_checkpoint import get_published_pretrained_checkpoint

import isaaclab_tasks  # noqa: F401
from isaaclab_tasks.utils import get_checkpoint_path
from isaaclab_tasks.utils.hydra import hydra_task_config

import TwoWheeledRobot.tasks  # noqa: F401


class CinematicCamera:
    """Per-step camera driver for rotating shots.

    ViewerCfg's built-in tracking only reads positions (root_pos_w / body_pos_w)
    and adds the eye/lookat offsets along world axes, so the camera translates
    with the robot but never turns with it. This drives the camera directly
    instead, rotating the eye offset about the robot each step.

    Two modes:
      * follow-yaw — azimuth locked to the robot's heading, so the shot stays
        at a fixed angle relative to the robot as it drives around.
      * orbit — azimuth swept a fixed number of turns over the capture.

    Both rotate only about world Z. The camera is never rolled with the body:
    at the ~47 deg one-leg lean, rolling the camera would level the robot on
    screen and tilt the world instead, which reads as "not balancing".
    """

    def __init__(self, env, eye, lookat, turns: float | None, total_steps: int | None, track_robot: bool = True):
        self._env = env
        self._track_robot = track_robot
        self._lookat = tuple(lookat)
        # Polar form of the eye offset, so only the azimuth needs updating.
        dx = eye[0] - lookat[0]
        dy = eye[1] - lookat[1]
        self._radius = math.hypot(dx, dy)
        self._height = eye[2] - lookat[2]
        self._base_azimuth = math.atan2(dy, dx)
        self._turns = turns
        # Orbit pacing: complete `turns` revolutions over the capture length.
        # Falls back to 1000 steps (20 s) when not recording, so the motion is
        # still visible in a live viewport.
        self._period_steps = total_steps if total_steps else 1000

    def _heading(self) -> float:
        """World-frame heading of the robot, in radians, from its forward axis.

        Taken as the ground projection of the body x-axis rather than a
        quaternion-to-euler yaw: the one-leg stance tilts x ~47 deg out of
        horizontal, and a naive euler yaw couples to that lean and makes the
        camera swing as the robot rocks. The projection stays stable.
        """
        from isaaclab.utils.math import quat_apply

        robot = self._env.scene["robot"]
        quat = robot.data.root_quat_w[0:1]
        fwd = quat_apply(quat, torch.tensor([[1.0, 0.0, 0.0]], device=quat.device))[0]
        fx, fy = float(fwd[0]), float(fwd[1])
        if math.hypot(fx, fy) < 1e-4:
            # Forward axis is near-vertical; heading is undefined, hold the last.
            return self._base_azimuth
        return math.atan2(fy, fx)

    def update(self, timestep: int) -> None:
        if self._turns is not None:
            azimuth = self._base_azimuth + 2.0 * math.pi * self._turns * (timestep / self._period_steps)
        else:
            azimuth = self._base_azimuth + self._heading()

        if self._track_robot:
            origin = self._env.scene["robot"].data.root_pos_w[0]
            ox, oy, oz = float(origin[0]), float(origin[1]), float(origin[2])
        else:
            ox, oy, oz = 0.0, 0.0, 0.0

        target = (ox + self._lookat[0], oy + self._lookat[1], oz + self._lookat[2])
        eye = (
            target[0] + self._radius * math.cos(azimuth),
            target[1] + self._radius * math.sin(azimuth),
            target[2] + self._height,
        )
        self._env.sim.set_camera_view(eye=eye, target=target)


def _apply_cinematic_scene() -> None:
    """Re-light the scene and mute the ground grid for presentation captures.

    Runs after env creation, so it edits prims the env already spawned:
    the default scene is a single dome light at intensity 2000 plus the stock
    grid ground plane, which reads as a debug view. Here the dome is dimmed to
    an ambient fill and a distant key light is added for a directional contact
    shadow under the stance wheel -- that shadow is what makes the robot look
    like it is standing on the floor rather than floating above it.

    Purely cosmetic: touches no physics, no collision, no articulation state.
    Best-effort -- a failure here must never block a capture, and the prim
    paths involved ("/World/Light", the grid asset's shader) are conventions of
    the stock assets rather than guaranteed API.
    """
    import isaaclab.sim as sim_utils

    # Dim the existing dome to a fill light. Spawning a second dome would
    # double the ambient term and flatten the image further.
    try:
        from pxr import UsdLux

        import omni.usd

        stage = omni.usd.get_context().get_stage()
        dome = stage.GetPrimAtPath("/World/Light")
        if dome.IsValid():
            UsdLux.LightAPI(dome).GetIntensityAttr().Set(600.0)
    except Exception as exc:  # noqa: BLE001
        print(f"[WARNING]: cinematic — could not dim the dome light ({exc}); using it as-is.")

    # Key light: angled from above/front so the robot casts a shadow toward the
    # camera side. `angle` softens the shadow edge (0.53 deg is the real sun).
    try:
        import torch as _torch

        from isaaclab.utils.math import quat_from_euler_xyz

        quat = quat_from_euler_xyz(
            _torch.tensor(-0.85),  # tilt down from horizontal
            _torch.tensor(0.0),
            _torch.tensor(0.6),  # swing around to the camera side
        )
        key_cfg = sim_utils.DistantLightCfg(intensity=2500.0, color=(1.0, 0.97, 0.92), angle=2.0)
        key_cfg.func("/World/CinematicKey", key_cfg, orientation=tuple(quat.tolist()))
    except Exception as exc:  # noqa: BLE001
        print(f"[WARNING]: cinematic — key light failed ({exc}); scene keeps the default lighting.")

    # Mute the ground grid. The stock plane's grid lines are the single biggest
    # "this is a simulator" tell. diffuse_tint multiplies the grid texture, so
    # this softens the lines rather than removing them outright.
    try:
        from pxr import Gf, Sdf

        from isaaclab.sim.utils import change_prim_property

        change_prim_property(
            prop_path="/World/ground/Looks/theGrid/Shader.inputs:diffuse_tint",
            value=Gf.Vec3f(0.22, 0.22, 0.24),
            type_to_create_if_not_exist=Sdf.ValueTypeNames.Color3f,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"[WARNING]: cinematic — ground tint failed ({exc}); keeping the default grid.")

    print("[INFO]: Cinematic scene applied (tracking camera, key light, muted ground).")


@hydra_task_config(args_cli.task, args_cli.agent)
def main(env_cfg: ManagerBasedRLEnvCfg | DirectRLEnvCfg | DirectMARLEnvCfg, agent_cfg: RslRlBaseRunnerCfg):
    """Play with RSL-RL agent."""
    # grab task name for checkpoint path
    task_name = args_cli.task.split(":")[-1]
    train_task_name = task_name.replace("-Play", "")

    # override configurations with non-hydra CLI arguments
    agent_cfg: RslRlBaseRunnerCfg = cli_args.update_rsl_rl_cfg(agent_cfg, args_cli)
    # Convert deprecated policy: RslRlPpoActorCriticCfg → actor/critic: RslRlMLPModelCfg
    import importlib.metadata as _meta
    agent_cfg = handle_deprecated_rsl_rl_cfg(agent_cfg, _meta.version("rsl-rl-lib"))
    env_cfg.scene.num_envs = args_cli.num_envs if args_cli.num_envs is not None else env_cfg.scene.num_envs

    # optional: spawn all envs at an exact roll angle (one-leg balance task).
    # The one-leg task samples roll = spawn_roll_center_deg ± curriculum band.
    # Pin it to an exact angle by centring on the requested roll and collapsing
    # the curriculum to a single zero-width stage. --spawn_roll_deg_list sweeps
    # multiple angles within this SAME process (see the sweep block below,
    # after env/policy are ready) -- it only needs a valid single value here to
    # get the env constructed; the real per-angle work happens via
    # env.unwrapped.set_spawn_roll_deg() + env.reset() at runtime, since
    # spawn_roll_center_deg is cached at env __init__ and mutating env_cfg
    # after construction has no effect.
    _sweep_angles = (
        [float(v) for v in args_cli.spawn_roll_deg_list.split(",")] if args_cli.spawn_roll_deg_list else None
    )
    _initial_spawn_roll_deg = args_cli.spawn_roll_deg if args_cli.spawn_roll_deg is not None else (
        _sweep_angles[0] if _sweep_angles else None
    )
    if _initial_spawn_roll_deg is not None:
        if hasattr(env_cfg, "spawn_roll_center_deg"):
            env_cfg.spawn_roll_center_deg = _initial_spawn_roll_deg
            env_cfg.spawn_roll_curriculum = ((0, 0.0),)
            print(f"[INFO] Spawn roll override: all envs spawn at exactly {_initial_spawn_roll_deg:.1f} deg roll.")
        else:
            print("[WARNING] --spawn_roll_deg/--spawn_roll_deg_list given, but this task has no spawn_roll_center_deg — ignored.")

    # set the environment seed
    # note: certain randomizations occur in the environment initialization so we set the seed here
    env_cfg.seed = agent_cfg.seed
    env_cfg.sim.device = args_cli.device if args_cli.device is not None else env_cfg.sim.device
    # With fabric on (the default), PhysX writes transforms to the Fabric/USDRT
    # layer and the USD stage is never updated -- so anything reading the stage
    # (Stage Recorder, USD-based capture tooling) sees a frozen scene. Turning
    # fabric off makes physics write through to USD. Slower, but fine at the
    # handful of envs used for rendering/recording.
    if args_cli.disable_fabric:
        env_cfg.sim.use_fabric = False

    # specify directory for logging experiments
    # Resolve relative to the current working directory first (original
    # behavior), but if that does not exist, fall back to the repo-root-relative
    # location. Training and play may be launched from different directories
    # (e.g. /workspace vs /workspace/Project inside the container) while the
    # checkpoints always live under <repo>/logs — this makes play find them
    # regardless of where it is launched from.
    log_root_path = os.path.abspath(os.path.join("logs", "rsl_rl", agent_cfg.experiment_name))
    if not os.path.isdir(log_root_path):
        _repo_root = Path(__file__).resolve().parents[2]
        repo_log_root = str(_repo_root / "logs" / "rsl_rl" / agent_cfg.experiment_name)
        if os.path.isdir(repo_log_root):
            print(f"[INFO] '{log_root_path}' not found — using repo-root logs at '{repo_log_root}'.")
            log_root_path = repo_log_root
    print(f"[INFO] Loading experiment from directory: {log_root_path}")
    if args_cli.use_pretrained_checkpoint:
        resume_path = get_published_pretrained_checkpoint("rsl_rl", train_task_name)
        if not resume_path:
            print("[INFO] Unfortunately a pre-trained checkpoint is currently unavailable for this task.")
            return
    elif args_cli.checkpoint:
        resume_path = retrieve_file_path(args_cli.checkpoint)
    else:
        resume_path = get_checkpoint_path(log_root_path, agent_cfg.load_run, agent_cfg.load_checkpoint)

    log_dir = os.path.dirname(resume_path)

    # set the log directory for the environment (works for all environment types)
    env_cfg.log_dir = log_dir

    if args_cli.episode_length is not None:
        # Only moves the timeout reset; fall terminations are unchanged, so a
        # policy that loses balance still ends its episode when it should.
        print(f"[INFO]: episode_length_s {env_cfg.episode_length_s} → {args_cli.episode_length}")
        env_cfg.episode_length_s = args_cli.episode_length

    if args_cli.disturbance is not None:
        # benchmark_disturbance_kind forces the kind on every reset regardless of
        # curriculum_stage, which defaults to 1 — below the stage-3 gate where the
        # generator starts assigning disturbances at all.
        if not hasattr(env_cfg, "benchmark_disturbance_kind"):
            print(f"[WARNING]: --disturbance given, but {task_name} has no disturbance support — ignoring.")
        else:
            env_cfg.benchmark_disturbance_kind = args_cli.disturbance
            print(f"[INFO]: Forcing disturbance every episode: {args_cli.disturbance}")
            # World-frame pushes: the one-leg stance is rolled ~47 deg, so a
            # body-frame X push would spend most of its magnitude vertically.
            if hasattr(env_cfg, "disturbance_forces_global"):
                env_cfg.disturbance_forces_global = True
                print("[INFO]: Disturbance forces applied in world frame (horizontal push).")

    if args_cli.cinematic:
        if args_cli.num_envs > 1:
            # Fleet shot: hold a fixed world camera above the whole env grid
            # rather than tracking one robot. The grid is laid out on
            # env_spacing centres around the origin, so the extent grows as
            # sqrt(num_envs) * spacing; back the camera off to match and look
            # down the diagonal so the rows do not overlap into a single line.
            spacing = getattr(getattr(env_cfg, "scene", None), "env_spacing", 2.5)
            half_extent = 0.5 * math.sqrt(args_cli.num_envs) * spacing
            dist = max(4.0, 3.0 * half_extent)
            env_cfg.viewer = ViewerCfg(
                eye=(0.60 * dist, -0.60 * dist, 0.45 * dist),
                lookat=(0.0, 0.0, 0.0),
                origin_type="world",
                resolution=(1920, 1080),
            )
            print(f"[INFO]: Fleet camera: {args_cli.num_envs} envs, grid half-extent {half_extent:.1f} m.")
        else:
            # Single-robot shot: camera rides with the robot (asset_root), so it
            # stays framed as the robot drives instead of leaving the viewport.
            # Eye/lookat are relative to the robot root: 1.1 m out at azimuth
            # -90 deg (square side-on, along -Y), tilted 20 deg down onto the
            # root point. To re-derive an eye from angles:
            #   x = r*cos(az),  y = r*sin(az),  z = lookat_z + r*tan(elevation)
            env_cfg.viewer = ViewerCfg(
                eye=(0.0, -1.1, 0.4),
                lookat=(0.0, 0.0, 0.0),
                origin_type="asset_root",
                asset_name="robot",
                env_index=0,
                resolution=(1920, 1080),
            )
        if args_cli.cinematic_orbit is not None or args_cli.cinematic_follow_yaw:
            # CinematicCamera drives the camera itself. Leaving origin_type as
            # "asset_root" makes the built-in tracker re-apply the fixed offset
            # on every render step, fighting the driver: the camera flips
            # between the two positions each frame and the captured frame lands
            # on whichever wrote last. "world" is not updated per-step, so the
            # driver has sole control. It still follows the robot — the driver
            # reads root_pos_w itself.
            env_cfg.viewer.origin_type = "world"

    # create isaac environment
    env = gym.make(args_cli.task, cfg=env_cfg, render_mode="rgb_array" if args_cli.video else None)

    if args_cli.cinematic:
        _apply_cinematic_scene()

    if args_cli.disturbance is not None:
        # OneLegBalanceEnv clears _body_ids in __init__ to skip the per-step
        # set_external_force_and_torque call (a GPU sync point it does not need
        # during training). The parent's disturbance block is gated on that
        # attribute, so without restoring it the sampled pushes are never
        # applied and the flag silently does nothing.
        unwrapped_env = env.unwrapped
        if getattr(unwrapped_env, "_body_ids", None) is None and hasattr(unwrapped_env, "_resolve_push_body_ids"):
            unwrapped_env._body_ids = unwrapped_env._resolve_push_body_ids()
            if unwrapped_env._body_ids is None:
                print("[WARNING]: could not resolve a body to push — disturbances will not be applied.")
            else:
                print("[INFO]: Re-enabled external disturbance forces (disabled by default on this task).")

    # Optional per-step camera driver for rotating shots. Built from the same
    # eye/lookat as the cinematic ViewerCfg above, and takes over from the
    # built-in tracker once the loop starts.
    cinematic_camera = None
    if args_cli.cinematic_orbit is not None or args_cli.cinematic_follow_yaw:
        if not args_cli.cinematic:
            print("[WARNING]: --cinematic_orbit/--cinematic_follow_yaw need --cinematic — ignoring.")
        else:
            if args_cli.cinematic_orbit is not None and args_cli.cinematic_follow_yaw:
                print("[WARNING]: --cinematic_orbit overrides --cinematic_follow_yaw.")
            cinematic_camera = CinematicCamera(
                env.unwrapped,
                eye=env_cfg.viewer.eye,
                lookat=env_cfg.viewer.lookat,
                turns=args_cli.cinematic_orbit,
                total_steps=args_cli.video_length if args_cli.video else None,
                # Multi-env shots orbit the grid centre, not one robot.
                track_robot=args_cli.num_envs == 1,
            )
            mode = f"orbit {args_cli.cinematic_orbit} turns" if args_cli.cinematic_orbit else "follow-yaw"
            print(f"[INFO]: Cinematic camera: {mode}.")

    # convert to single-agent instance if required by the RL algorithm
    if isinstance(env.unwrapped, DirectMARLEnv):
        env = multi_agent_to_single_agent(env)

    # wrap for video recording
    if args_cli.video:
        # Timestamped prefix: RecordVideo's default name ("rl-video") plus the
        # step trigger yields the same filename every run, so consecutive takes
        # overwrite each other. Stamping the launch time keeps every take.
        video_name_prefix = time.strftime("play-%Y%m%d-%H%M%S")
        if args_cli.cinematic:
            video_name_prefix += "-cinematic"
        video_kwargs = {
            "video_folder": os.path.join(log_dir, "videos", "play"),
            "step_trigger": lambda step: step == 0,
            "name_prefix": video_name_prefix,
            "video_length": args_cli.video_length,
            "disable_logger": True,
        }
        print("[INFO] Recording videos during training.")
        print_dict(video_kwargs, nesting=4)
        env = gym.wrappers.RecordVideo(env, **video_kwargs)

    # wrap around environment for rsl-rl
    env = RslRlVecEnvWrapper(env, clip_actions=agent_cfg.clip_actions)

    print(f"[INFO]: Loading model checkpoint from: {resume_path}")
    # load previously trained model
    if agent_cfg.class_name == "OnPolicyRunner":
        runner = OnPolicyRunner(env, agent_cfg.to_dict(), log_dir=None, device=agent_cfg.device)
    elif agent_cfg.class_name == "DistillationRunner":
        from rsl_rl.runners import DistillationRunner
        runner = DistillationRunner(env, agent_cfg.to_dict(), log_dir=None, device=agent_cfg.device)
    else:
        raise ValueError(f"Unsupported runner class: {agent_cfg.class_name}")
    runner.load(resume_path)

    # obtain the trained policy for inference
    policy = runner.get_inference_policy(device=env.unwrapped.device)

    # extract the neural network module
    # we do this in a try-except to maintain backwards compatibility.
    try:
        # version 2.3 onwards
        policy_nn = runner.alg.policy
    except AttributeError:
        try:
            # version 2.2 and below
            policy_nn = runner.alg.actor_critic
        except AttributeError:
            # rsl-rl >= 4.0: separate actor / critic models
            policy_nn = runner.alg.actor

    # extract the normalizer
    if hasattr(policy_nn, "actor_obs_normalizer"):
        normalizer = policy_nn.actor_obs_normalizer
    elif hasattr(policy_nn, "student_obs_normalizer"):
        normalizer = policy_nn.student_obs_normalizer
    else:
        normalizer = None

    # export policy to jit/onnx — wrapped so a failed export never blocks play
    export_model_dir = os.path.join(os.path.dirname(resume_path), "exported")
    try:
        export_policy_as_jit(policy_nn, normalizer=normalizer, path=export_model_dir, filename="policy.pt")
        print(f"[INFO]: JIT policy exported to: {export_model_dir}/policy.pt")
    except Exception as e:
        print(f"[WARNING]: JIT export failed ({e}) — skipping.")
    try:
        export_policy_as_onnx(policy_nn, normalizer=normalizer, path=export_model_dir, filename="policy.onnx")
        print(f"[INFO]: ONNX policy exported to: {export_model_dir}/policy.onnx")
    except Exception as e:
        print(f"[WARNING]: ONNX export failed ({e}) — skipping.")

    dt = env.unwrapped.step_dt

    def run_segment(num_steps: int, angle_label: float, run_forever: bool) -> None:
        """Step the env, optionally measuring survival, for num_steps control
        steps (or indefinitely if run_forever). One call = one sweep stage
        when --spawn_roll_deg_list is given, or the whole run otherwise."""
        obs = env.get_observations()
        unwrapped = env.unwrapped
        max_ep_len = int(unwrapped.max_episode_length)
        step_count = torch.zeros(unwrapped.num_envs, dtype=torch.long, device=unwrapped.device)
        survived = 0   # episodes that reached the full episode length
        fell = 0       # episodes that terminated early
        ep_len_sum = 0  # total steps over all finished episodes (for the mean)

        timestep = 0
        push_was_active = False  # rising-edge latch for the push indicator
        while simulation_app.is_running():
            start_time = time.time()
            # Drive the camera *before* stepping: the frame RecordVideo captures
            # is rendered inside env.step(), so updating afterwards would apply
            # each camera position to the following frame.
            if cinematic_camera is not None:
                cinematic_camera.update(timestep)

            with torch.inference_mode():
                actions = policy(obs)
                obs, _, dones, extras = env.step(actions)
                # reset recurrent states for episodes that have terminated
                policy_nn.reset(dones)
            timestep += 1

            # Announce each push on its rising edge. Nothing in the viewport
            # renders external forces, so this is the only way to correlate a
            # stumble on screen with an actual disturbance.
            if args_cli.disturbance is not None and hasattr(unwrapped, "_last_disturbance_force"):
                force_n = float(unwrapped._last_disturbance_force[0].norm())
                torque_nm = float(unwrapped._last_disturbance_torque[0].norm())
                push_active = force_n > 1.0e-6 or torque_nm > 1.0e-6
                if push_active and not push_was_active:
                    print(f"[step {timestep:5d}]  ### PUSH  force={force_n:.2f} N  torque={torque_nm:.3f} Nm")
                push_was_active = push_active

            # classify each finished episode as survival (timed out) or fall (early termination)
            if args_cli.measure_survival:
                step_count += 1
                done_idx = torch.nonzero(dones, as_tuple=False).flatten()
                if len(done_idx) > 0:
                    # RslRlVecEnvWrapper reports timeouts (episode reached full length)
                    # in extras["time_outs"]; fall back to the env flag if absent.
                    time_outs = extras.get("time_outs") if isinstance(extras, dict) else None
                    if time_outs is None:
                        time_outs = getattr(unwrapped, "_last_timeout", None)
                    for i in done_idx.tolist():
                        n = int(step_count[i].item())
                        ep_len_sum += n
                        is_timeout = bool(time_outs[i]) if time_outs is not None else (n >= max_ep_len)
                        if is_timeout:
                            survived += 1
                        else:
                            fell += 1
                        step_count[i] = 0

            # periodic console output so the user can see it is running
            if timestep % 100 == 0:
                if hasattr(unwrapped, "_controller"):
                    # SMC / PETASMC env — print controller internals
                    ctrl = unwrapped._controller
                    print(
                        f"[step {timestep:5d}]  "
                        f"tilt={ctrl.last.get('tilt', 0.0):+6.2f}°  "
                        f"IET_bal={ctrl.last.get('iet_bal', 0.0)*1000:.1f}ms  "
                        f"T_min_bal={float(ctrl.et_min_iet_bal[0] if hasattr(ctrl.et_min_iet_bal, '__len__') else ctrl.et_min_iet_bal)*1000:.1f}ms  "
                        f"K_max_bal={float(ctrl.et_k_max_bal[0] if hasattr(ctrl.et_k_max_bal, '__len__') else ctrl.et_k_max_bal):.2f}A"
                    )
                else:
                    # Stand-up env — print episode stats
                    ep_len = unwrapped.episode_length_buf.float().mean().item()
                    print(f"[step {timestep:5d}]  spawn_roll={angle_label:.1f}deg  mean_ep_len={ep_len:.1f}")
                    # CyberGear joint torques and wheel current (env 0 only)
                    if hasattr(unwrapped, "_cg_ids"):
                        cg_torque = unwrapped.robot.data.applied_torque[:, unwrapped._cg_ids].abs().cpu()
                        wheel_cur = unwrapped._action_processor.command_current.abs().cpu()
                        print(
                            f"           CG max |torque| (Nm) FL={cg_torque[:, 0].max():5.3f}  FR={cg_torque[:, 1].max():5.3f}"
                            f"  BL={cg_torque[:, 2].max():5.3f}  BR={cg_torque[:, 3].max():5.3f}"
                            f"  |  wheel L max={wheel_cur[:, 0].max():5.3f}A  R max={wheel_cur[:, 1].max():5.3f}A"
                        )

            if args_cli.video and timestep == args_cli.video_length:
                break
            if not run_forever and timestep >= num_steps:
                print(f"[INFO]: Reached {num_steps} steps — stopping.")
                break

            # time delay for real-time evaluation
            sleep_time = dt - (time.time() - start_time)
            if args_cli.real_time and sleep_time > 0:
                time.sleep(sleep_time)

        # survival-rate summary
        if args_cli.measure_survival:
            n_ep = survived + fell
            if n_ep == 0:
                print(f"[SURVIVAL] spawn_roll = {angle_label:.1f} deg   No episode finished — increase --num_steps.")
            else:
                rate = 100.0 * survived / n_ep
                mean_len = ep_len_sum / n_ep
                print("\n" + "=" * 60)
                print(f"[SURVIVAL] spawn_roll = {angle_label:.1f} deg   (episode length = {max_ep_len} steps)")
                print(f"[SURVIVAL] episodes: {n_ep}   survived: {survived}   fell: {fell}")
                print(f"[SURVIVAL] survival rate: {rate:.1f} %   mean episode length: {mean_len:.1f} steps")
                print("=" * 60)

    if _sweep_angles is not None:
        if args_cli.num_steps <= 0:
            raise ValueError("--spawn_roll_deg_list requires --num_steps > 0 (fixed length per angle).")
        if not hasattr(env.unwrapped, "set_spawn_roll_deg"):
            raise RuntimeError("--spawn_roll_deg_list given, but this task has no set_spawn_roll_deg.")
        print(f"[INFO]: Sweeping {len(_sweep_angles)} spawn angles in one process: {_sweep_angles}")
        for angle in _sweep_angles:
            print(f"\n[INFO] === spawn_roll = {angle:.1f} deg ===")
            # Force every env to respawn at the new angle immediately, rather
            # than waiting out the previous episode's remaining length.
            # Wrapped in inference_mode() because run_segment()'s env.step()
            # calls are too -- some env state tensors (e.g. _prev_actions)
            # get reallocated (.clone()-style, not just mutated) during those
            # steps, which marks them as "inference tensors" that can only be
            # written in-place from inside an inference_mode() context. This
            # reset is the first time the program writes to them from outside
            # one (normally all resets happen automatically inside step(),
            # already inside the caller's inference_mode() block).
            with torch.inference_mode():
                env.unwrapped.set_spawn_roll_deg(angle, 0.0)
                env.reset()
                all_dones = torch.ones(env.unwrapped.num_envs, dtype=torch.bool, device=env.unwrapped.device)
                policy_nn.reset(all_dones)
            run_segment(args_cli.num_steps, angle, run_forever=False)
    else:
        print("[INFO]: Simulation running. Press Ctrl+C to stop.")
        angle_label = args_cli.spawn_roll_deg if args_cli.spawn_roll_deg is not None else float("nan")
        run_segment(args_cli.num_steps, angle_label, run_forever=(args_cli.num_steps <= 0))

    # close the simulator
    env.close()


if __name__ == "__main__":
    # run the main function
    main()
    # close sim app
    simulation_app.close()
