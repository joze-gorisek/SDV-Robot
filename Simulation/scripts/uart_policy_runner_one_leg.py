#!/usr/bin/env python3
"""
Run an exported one-leg balance GRU policy against hardware over UART.

Purpose:
    Convert STM32 JSON telemetry into the 6-D observation used by
    OneLegBalanceEnv, run the exported *recurrent* TorchScript policy, and send
    the left-wheel current plus the two left CyberGear joint targets back over
    serial.

This is the OneLegBalance counterpart of scripts/uart_policy_runner.py (which
serves the stateless 18-D / 6-D StandupEnv policy). Keep the two separate — the
observation layout, action layout, and statefulness all differ.

Stateful (GRU) contract — READ THIS:
    The exported TorchScript policy is RECURRENT. It carries a GRU hidden state
    (rnn_num_layers × 1 × rnn_hidden_dim = 1 × 1 × 32 floats) *inside* the module
    across successive forward() calls, so you do NOT pass it in or out here — but
    you MUST call ``policy.reset()`` whenever balancing (re)starts, otherwise the
    hidden state leaks stale history from the previous run into the new one.
    This runner calls reset() once at startup; call it again on every fresh
    entry into the one-leg stance.

    For a hand-written STM32 C implementation, use the ONNX export instead
    (policy.onnx): its signature is (obs[12], h_in[1,1,32]) -> (actions[3],
    h_out[1,1,32]). Persist h between 50 Hz steps, feed h_out back as the next
    h_in, and zero h on stance entry. See STM32_DEPLOYMENT.md.

Edit here when:
    Hardware packet fields, angle units, the 6-D observation normalization, the
    5-D action scaling, or safe-error behavior changes. Keep this in lockstep
    with one_leg_balance_env.py::_get_observations()/_pre_physics_step() and the
    observation_scale / joint / current constants in
    one_leg_balance_env_cfg.py + pure_nn_balance_env_cfg.py.

STM32 -> host, one JSON object per line:
    {
      "roll": 0.82, "pitch": 0.0, "yaw": 0.0,
      "gyro": [0.0, 0.0, 0.0],     # body-frame [wx, wy, wz] rad/s
      "wheel_vel_left": 0.0        # left DDSM115 wheel velocity, rad/s
    }

Host -> STM32, one JSON object per line:
    {
      "cg_target": [fl, fr, bl, br],   # radians; fr/br are the locked top joints (0.0)
      "wheel_current": [left, right],  # amps; right is always 0.0 (lifted)
      "action": [wheel_L, joint_FL, joint_BL, joint_FR, joint_BR]   # raw policy outputs in [-1, 1]
    }

Angles are radians by default. Pass --degrees if the STM32 packet uses degrees.
The exported policy is expected to be the JIT file produced by
scripts/rsl_rl/play.py for Template-Twowheeledrobot-OneLegBalance-v0.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass
from typing import Any

import serial
import torch


# ── Physical / control constants (mirror the Isaac env) ─────────────────────
CONTROL_DT = 0.020          # 50 Hz control step (must match the sim step_dt)
R_WHEEL = 0.05035           # wheel radius (m), residual_lqr_env.LQR_PHYSICAL_PARAMS
I_MAX_A = 2.0               # wheel tanh current ceiling, pure_nn_balance_env_cfg.i_max_a
WHEEL_SIGN_LEFT = -1.0      # _wheel_sign[0]; forward-positive ground speed

# CyberGear joint order everywhere: [front_left, front_right, back_left, back_right]
CG_EXT_SIGN = (1.0, -1.0, -1.0, 1.0)   # _cg_ext_sign
BOTTOM_JOINT_CENTER_RAD = math.radians(45.0)   # bottom_joint_center_deg
BOTTOM_JOINT_RANGE_RAD = math.radians(15.0)    # bottom_joint_range_deg
TOP_JOINT_CENTER_RAD = math.radians(45.0)      # top_joint_center_deg
TOP_JOINT_RANGE_RAD = math.radians(45.0)       # top_joint_range_deg (0..90 ext)

# Raw joint limits (rad): FL/BR = [-10deg, +90deg], FR/BL = [-90deg, +10deg]
LO_EXT = math.radians(10.0)
LIM = math.pi / 2.0
CG_JOINT_LO = (-LO_EXT, -LIM, -LIM, -LO_EXT)
CG_JOINT_HI = (LIM, LO_EXT, LO_EXT, LIM)

# 6-D observation normalization — MUST match one_leg_balance_env_cfg.observation_scale
OBS_SCALE = torch.tensor(
    [
        math.radians(50.0),  # 0  roll
        4.0,                 # 1  roll rate
        math.radians(25.0),  # 2  pitch
        4.0,                 # 3  pitch rate
        0.5,                 # 4  x_rel (m)
        0.5,                 # 5  velocity (m/s)
    ],
    dtype=torch.float32,
)


@dataclass
class RobotState:
    roll: float
    pitch: float
    yaw: float
    gyro: list[float]        # [wx, wy, wz] rad/s
    wheel_vel_left: float    # rad/s


def projected_gravity_from_roll_pitch(roll: float, pitch: float) -> tuple[float, float, float]:
    """Roll/pitch -> Isaac-style projected_gravity_b (same as the standup runner).

    Upright is [0, 0, -1]; roll shows up mostly in x, pitch mostly in y. If the
    hardware signs do not match simulation, fix the sign at this input layer
    before trusting the policy on the real robot.
    """
    sr, cr = math.sin(roll), math.cos(roll)
    sp, cp = math.sin(pitch), math.cos(pitch)
    return (-sr * cp, sp, -cr * cp)


def roll_from_projected_gravity(g: tuple[float, float, float]) -> float:
    """atan2(grav_x, -grav_z) — identical to pure_nn_components.roll_from_projected_gravity."""
    return math.atan2(g[0], -g[2])


def pitch_from_projected_gravity(g: tuple[float, float, float]) -> float:
    """atan2(grav_y, -grav_z) — identical to pure_nn_components.pitch_from_projected_gravity."""
    return math.atan2(g[1], -g[2])


def parse_state(line: bytes, use_degrees: bool) -> RobotState:
    msg: dict[str, Any] = json.loads(line.decode("utf-8"))
    roll = float(msg["roll"])
    pitch = float(msg["pitch"])
    yaw = float(msg.get("yaw", 0.0))
    gyro = [float(x) for x in msg["gyro"]]
    wheel_vel_left = float(msg.get("wheel_vel_left", 0.0))

    if use_degrees:
        roll = math.radians(roll)
        pitch = math.radians(pitch)
        yaw = math.radians(yaw)
        gyro = [math.radians(x) for x in gyro]

    if len(gyro) != 3:
        raise ValueError("gyro must contain [wx, wy, wz]")

    return RobotState(roll=roll, pitch=pitch, yaw=yaw, gyro=gyro, wheel_vel_left=wheel_vel_left)


def actions_to_commands(action: torch.Tensor) -> tuple[list[float], list[float]]:
    """Map the 5-D policy output [wheel_L, joint_FL, joint_BL, joint_FR, joint_BR]
    to hardware commands.

    Mirrors one_leg_balance_env._pre_physics_step (all four joints are policy-
    controlled) and the parent CurrentActionProcessor tanh scaling (wheel
    current). The right wheel is lifted (zero current).
    """
    action = torch.clamp(action, -1.0, 1.0)

    # Left wheel: tanh(raw) * i_max, same as CurrentActionProcessor.net_current.
    wheel_l = float(math.tanh(float(action[0])) * I_MAX_A)

    # Bottom joints (FL, BL): narrow range. Top joints (FR, BR): wide range.
    ext_fl = BOTTOM_JOINT_CENTER_RAD + float(action[1]) * BOTTOM_JOINT_RANGE_RAD
    ext_bl = BOTTOM_JOINT_CENTER_RAD + float(action[2]) * BOTTOM_JOINT_RANGE_RAD
    ext_fr = TOP_JOINT_CENTER_RAD + float(action[3]) * TOP_JOINT_RANGE_RAD
    ext_br = TOP_JOINT_CENTER_RAD + float(action[4]) * TOP_JOINT_RANGE_RAD
    raw_fl = CG_EXT_SIGN[0] * ext_fl   # front_left,  ext_sign = +1
    raw_bl = CG_EXT_SIGN[2] * ext_bl   # back_left,   ext_sign = -1
    raw_fr = CG_EXT_SIGN[1] * ext_fr   # front_right, ext_sign = -1
    raw_br = CG_EXT_SIGN[3] * ext_br   # back_right,  ext_sign = +1

    # Clamp every joint to its hard mechanical limit.
    cg_target = [
        min(max(raw_fl, CG_JOINT_LO[0]), CG_JOINT_HI[0]),
        min(max(raw_fr, CG_JOINT_LO[1]), CG_JOINT_HI[1]),
        min(max(raw_bl, CG_JOINT_LO[2]), CG_JOINT_HI[2]),
        min(max(raw_br, CG_JOINT_LO[3]), CG_JOINT_HI[3]),
    ]
    wheel_current = [wheel_l, 0.0]   # right wheel lifted
    return cg_target, wheel_current


class OneLegObservationBuilder:
    """Stateful 6-D observation builder matching OneLegBalanceEnv.

    The only state it carries is the running left-wheel odometry (x_rel); the
    action/current history is now left to the GRU hidden state, so it is neither
    observed nor tracked here. Call reset() together with policy.reset() whenever
    balancing restarts.
    """

    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        self.x_rel = 0.0

    def build(self, state: RobotState) -> torch.Tensor:
        g = projected_gravity_from_roll_pitch(state.roll, state.pitch)
        roll = roll_from_projected_gravity(g)
        pitch = pitch_from_projected_gravity(g)
        roll_rate = state.gyro[1]        # ang_vel_b[1]
        pitch_rate = -state.gyro[0]      # -ang_vel_b[0]

        # Left-wheel ground speed + odometry integration (mirrors _compute_velocity).
        velocity = state.wheel_vel_left * WHEEL_SIGN_LEFT * R_WHEEL
        self.x_rel += velocity * CONTROL_DT

        obs_raw = torch.tensor(
            [roll, roll_rate, pitch, pitch_rate, self.x_rel, velocity],
            dtype=torch.float32,
        )
        obs = obs_raw / OBS_SCALE
        obs = torch.nan_to_num(obs, nan=0.0, posinf=10.0, neginf=-10.0).clamp(-10.0, 10.0)
        return obs.unsqueeze(0)


def _maybe_reset_policy(policy: torch.jit.ScriptModule) -> None:
    """Zero the GRU hidden state if the exported module supports it."""
    reset = getattr(policy, "reset", None)
    if callable(reset):
        try:
            reset()
        except Exception as exc:  # noqa: BLE001
            print(f"[uart_policy_runner_one_leg] policy.reset() failed: {exc}", file=sys.stderr)


def run(args: argparse.Namespace) -> int:
    policy = torch.jit.load(args.policy, map_location="cpu")
    policy.eval()
    _maybe_reset_policy(policy)   # start each session from a clean hidden state

    obs_builder = OneLegObservationBuilder()
    next_tick = time.monotonic()

    with serial.Serial(args.port, args.baud, timeout=args.timeout) as uart:
        print(f"[uart_policy_runner_one_leg] opened {args.port} @ {args.baud}", file=sys.stderr)
        while True:
            line = uart.readline()
            if not line:
                if args.fail_silent:
                    continue
                print("[uart_policy_runner_one_leg] UART timeout", file=sys.stderr)
                continue

            # A "reset" packet lets the STM32 signal a fresh stance entry, so both
            # the GRU hidden state and the odometry/history restart together.
            stripped = line.strip()
            if stripped == b"reset" or stripped == b'{"reset":true}':
                _maybe_reset_policy(policy)
                obs_builder.reset()
                print("[uart_policy_runner_one_leg] hidden state + odometry reset", file=sys.stderr)
                continue

            try:
                state = parse_state(line, args.degrees)
                obs = obs_builder.build(state)
                with torch.inference_mode():
                    action = policy(obs).squeeze(0).detach().cpu().to(torch.float32)
                action = torch.clamp(action, -1.0, 1.0)
                cg_targets, wheel_current = actions_to_commands(action)

                response = {
                    "cg_target": cg_targets,
                    "wheel_current": wheel_current,
                    "action": action.tolist(),
                }
                if args.echo_state:
                    response["x_rel"] = obs_builder.x_rel
                uart.write((json.dumps(response, separators=(",", ":")) + "\n").encode("utf-8"))

                if args.rate_hz > 0.0:
                    next_tick += 1.0 / args.rate_hz
                    sleep_s = next_tick - time.monotonic()
                    if sleep_s > 0.0:
                        time.sleep(sleep_s)
                    else:
                        next_tick = time.monotonic()
            except Exception as exc:  # noqa: BLE001
                print(f"[uart_policy_runner_one_leg] bad packet: {exc}; raw={line!r}", file=sys.stderr)
                if args.safe_on_error:
                    safe = {
                        "cg_target": [0.0, 0.0, 0.0, 0.0],
                        "wheel_current": [0.0, 0.0],
                        "action": [0.0, 0.0, 0.0, 0.0, 0.0],
                    }
                    uart.write((json.dumps(safe, separators=(",", ":")) + "\n").encode("utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the one-leg balance GRU policy over UART.")
    parser.add_argument("--policy", required=True, help="Path to exported JIT policy.pt (recurrent)")
    parser.add_argument("--port", required=True, help="Serial port, e.g. /dev/ttyACM0 or COM5")
    parser.add_argument("--baud", type=int, default=115200)
    parser.add_argument("--timeout", type=float, default=0.1)
    parser.add_argument("--rate-hz", type=float, default=50.0, help="Optional host-side loop cap. Use 0 to disable.")
    parser.add_argument("--degrees", action="store_true", help="UART angles and gyro are degrees / deg/s.")
    parser.add_argument("--echo-state", action="store_true", help="Echo x_rel in the response for debugging.")
    parser.add_argument("--safe-on-error", action="store_true", help="Send zero-safe command after malformed packets.")
    parser.add_argument("--fail-silent", action="store_true", help="Do not print UART timeout messages.")
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
