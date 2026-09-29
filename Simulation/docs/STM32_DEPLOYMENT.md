# STM32 Deployment Notes

These notes describe the cleaned stand-up policy interface. The Isaac Lab task is:

```text
Template-Twowheeledrobot-Standup-v0
```

## Policy Contract

The trained actor takes 18 normalized `float32` observations and returns 6 actions in `[-1, 1]`.

Observations:

```text
0..2    projected gravity in body frame
3..5    body angular velocity / 10 rad/s
6..9    CyberGear joint extension fraction
10..11  DDSM115 wheel velocity / 20.94 rad/s
12..17  previous action
```

If the STM32 sends roll/pitch/yaw instead of projected gravity, convert roll and pitch to projected gravity before inference. Yaw does not affect gravity direction and is not used by this policy.

CyberGear extension fractions are normalized joint positions in this order:

```text
front_left, front_right, back_left, back_right
```

The two mirrored joints use opposite signs so that positive normalized value means "more extended" on every leg. The conversion used by the Python UART runner is:

```text
cg_norm = ((cg_angle * [1,-1,-1,1] + 10deg) / 100deg) * 2 - 1
```

Actions:

```text
0..3  CyberGear target angles
4     left DDSM115 current command
5     right DDSM115 current command
```

The simulation maps wheel actions to current with `wheel_current_max` from `standup_env_cfg.py` and then torque with `DDSM115_KT` from `sim_params.py`. The conservative policy envelope is `0.96 Nm`, or `1.28 A` at `0.75 Nm/A`. The motor model clamps current to `2.7 A`, peak torque to `2.0 Nm`, and reduces available torque linearly to zero at the `200 rpm` no-load speed.

## Sign Convention

The left wheel USD is mirrored, so simulation negates left wheel torque:

```c
torque_left  = -current_left  * DDSM115_KT;
torque_right =  current_right * DDSM115_KT;
```

Keep the firmware-side motor direction mapping consistent with the physical wiring, not blindly with the USD. The important external behavior is that positive wheel action should help the learned policy perform the same maneuver on hardware as in simulation.

## Files To Keep Aligned

- `source/TwoWheeledRobot/TwoWheeledRobot/tasks/direct/twowheeledrobot/standup_env.py`
- `source/TwoWheeledRobot/TwoWheeledRobot/tasks/direct/twowheeledrobot/standup_env_cfg.py`
- `source/TwoWheeledRobot/TwoWheeledRobot/tasks/direct/twowheeledrobot/sim_params.py`
- `source/TwoWheeledRobot/TwoWheeledRobot/tasks/direct/twowheeledrobot/agents/rsl_rl_standup_cfg.py`
- `RealImplementationCode/`

When changing observation normalization, action scaling, motor signs, or the actor network shape, update the STM32 inference code in the same commit.

## Export

Train:

```bash
python scripts/rsl_rl/train.py \
  --task Template-Twowheeledrobot-Standup-v0 \
  --headless --num_envs 4096
```

Play/export from a checkpoint:

```bash
python scripts/rsl_rl/play.py \
  --task Template-Twowheeledrobot-Standup-v0 \
  --num_envs 1
```

Check exported policies under `logs/rsl_rl/standup_two_wheel/<run>/exported/`.

## Pure NN Balance Controller

Task:

```text
Template-Twowheeledrobot-PureNNBalance-v0
```

Policy rate is 50 Hz (`dt = 0.02 s`). The deployed current model takes 8 normalized `float32` observations and returns two wheel current commands in amperes:

```text
0  x_rel / 1.0 m
1  linear_velocity / 1.0 m/s
2  pitch / 25 deg
3  pitch_rate / 4.0 rad/s
4  yaw_error / pi rad
5  yaw_rate / 4.0 rad/s
6  previous_left_current / I_max
7  previous_right_current / I_max
```

The inference-ready ONNX wrapper applies `tanh(actor(obs)) * I_max`; default `I_max = 2.0 A`. Default training does not apply additional command smoothing or hard slew limiting. Optional hardware-safety testing can enable:

```text
I_filtered = alpha * I_previous + (1 - alpha) * I_network
optional slew limit = 0.3 A per 20 ms sample
```

Train all curriculum stages and export the current-output ONNX model:

```bash
python scripts/train_pure_nn_curriculum.py --num_envs 4096 --headless
```

Manual export from an existing `policy.pt`:

```bash
python scripts/export_pure_nn_current_onnx.py \
  --policy logs/rsl_rl/pure_nn_balance_two_wheel/<run>/exported/policy.pt \
  --output logs/rsl_rl/pure_nn_balance_two_wheel/<run>/exported/policy_current.onnx \
  --i-max-a 2.0
```

The export script runs ONNX Runtime validation and fails if max error is `>= 1e-4`.

Benchmark scenarios:

```bash
python scripts/benchmark_pure_nn_balance.py \
  --policy logs/rsl_rl/pure_nn_balance_two_wheel/<run>/exported/policy.pt \
  --num_envs 64 --headless
```

## UART Runner

Host-side runner:

```bash
python scripts/uart_policy_runner.py \
  --policy logs/rsl_rl/standup_two_wheel/<run>/exported/policy.pt \
  --port /dev/ttyACM0 \
  --baud 115200
```

STM32 sends JSON lines:

```json
{"roll":0.0,"pitch":1.57,"yaw":0.0,"gyro":[0.0,0.0,0.0],"cg":[0.0,0.0,0.0,0.0],"ddsm":[0.0,0.0]}
```

Host sends JSON lines:

```json
{"cg_target":[0.0,0.0,0.0,0.0],"wheel_current":[0.0,0.0],"action":[0.0,0.0,0.0,0.0,0.0,0.0]}
```

DDSM115 velocity is part of the policy observation. Send `[left, right]` wheel angular velocity in rad/s using the same sign convention as simulation.

## One-Leg Balance Controller (recurrent GRU — STATEFUL)

Task:

```text
Template-Twowheeledrobot-OneLegBalance-v0
```

**This policy is recurrent — the deployment contract is different from every
other task here.** The actor is a GRU (`rnn_type="gru"`, `rnn_num_layers=1`,
`rnn_hidden_dim=32`) followed by a `[32, 32]` MLP head. Policy rate is 50 Hz
(`dt = 0.02 s`). The robot balances on the **left** wheel at ~47° roll; the
right wheel is lifted (zero current). All four CyberGear joints are
policy-controlled — the top joints (front_right, back_right) are the tip-up
actuators and are no longer locked.

### Observation (6 normalized `float32`)

Order and normalization scale (must match `one_leg_balance_env_cfg.observation_scale`):

```text
0   roll        / 50 deg     roll = atan2(grav_x, -grav_z)
1   roll_rate   / 4.0 rad/s  = gyro_y  (ang_vel_b[1])
2   pitch       / 25 deg     pitch = atan2(grav_y, -grav_z)
3   pitch_rate  / 4.0 rad/s  = -gyro_x (-ang_vel_b[0])
4   x_rel       / 0.5 m      left-wheel odometry since stance entry
5   velocity    / 0.5 m/s    left-wheel ground speed
```

After dividing by the scale, apply `nan_to_num` then clamp to `[-10, 10]`.

Only physical state is observed — the previous wheel currents and previous
joint actions are NOT fed in; the GRU hidden state carries that history
internally. (The firmware still needs the gyro to form roll_rate/pitch_rate and
the left wheel velocity to form velocity/x_rel.)

`x_rel` and `velocity` come from the **left wheel only** (`velocity =
wheel_vel_left * wheel_sign_left * R_WHEEL`, `R_WHEEL = 0.05035 m`,
`wheel_sign_left = -1`); `x_rel` is the running integral `x_rel += velocity *
0.02` since the stance began. Do **not** use a two-wheel average — the lifted
right wheel would halve the estimate.

### Action (5 raw outputs in `[-1, 1]`)

```text
0  wheel_L   -> left wheel current  = tanh(action0) * I_max   (I_max = 2.0 A)
1  joint_FL  -> front_left  target  = ext_sign_FL * (45 deg + action1 * 15 deg)   # bottom, narrow
2  joint_BL  -> back_left   target  = ext_sign_BL * (45 deg + action2 * 15 deg)   # bottom, narrow
3  joint_FR  -> front_right target  = ext_sign_FR * (45 deg + action3 * 45 deg)   # top, wide (0..90)
4  joint_BR  -> back_right  target  = ext_sign_BR * (45 deg + action4 * 45 deg)   # top, wide (0..90)
```

`ext_sign = [FL,FR,BL,BR] = [+1,-1,-1,+1]`. All four CyberGear joints are
policy-controlled (the top joints are the tip-up actuators — no longer locked).
Emit `cg_target = [FL, FR, BL, BR]`. Right wheel current is always `0.0`. Clamp
each joint to its hard limit (FL/BR `[-10°,+90°]`, FR/BL `[-90°,+10°]`).

### Stateful GRU — hidden-state handling (the critical part)

The GRU carries a hidden state `h` of shape `[rnn_num_layers, 1, rnn_hidden_dim]
= [1, 1, 32]` (32 `float32`) that must persist **between** 50 Hz control steps.

- **TorchScript export (`policy.pt`)** keeps `h` *inside* the module across
  `forward()` calls and exposes a `reset()` method. `scripts/uart_policy_runner_one_leg.py`
  uses this: it calls `policy.reset()` at startup and whenever the STM32 sends a
  `reset` line, and never handles `h` explicitly.
- **ONNX export (`policy.onnx`)** exposes `h` *externally* — this is the format
  to use for a hand-written STM32 C implementation:

  ```text
  inputs:  obs   [1, 6]
           h_in  [1, 1, 32]
  outputs: actions [1, 5]
           h_out   [1, 1, 32]
  ```

  Firmware loop: keep a 32-float `h` buffer, feed it as `h_in`, run inference,
  store `h_out` back into `h` for the next step.

- **Reset semantics:** there is no natural "episode reset" in the field. You
  MUST zero `h` (and restart the `x_rel` odometry and the prev-action / prev-
  current history) every time the robot (re)enters the one-leg stance. Failing
  to reset leaks stale history from the previous run and degrades balance. The
  runner exposes this via a `reset` UART line.

### UART Runner

Host-side runner (separate from the stateless standup runner):

```bash
python scripts/uart_policy_runner_one_leg.py \
  --policy logs/rsl_rl/one_leg_balance_two_wheel/<run>/exported/policy.pt \
  --port /dev/ttyACM0 \
  --baud 115200
```

STM32 sends JSON lines (gyro is body-frame `[wx, wy, wz]`):

```json
{"roll":0.82,"pitch":0.0,"yaw":0.0,"gyro":[0.0,0.0,0.0],"wheel_vel_left":0.0}
```

Host sends JSON lines (`cg_target` = `[FL, FR, BL, BR]`, all four controlled):

```json
{"cg_target":[0.78,-0.78,-0.78,0.78],"wheel_current":[0.0,0.0],"action":[0.0,0.0,0.0,0.0,0.0]}
```

To restart balancing, send a bare `reset` line (or `{"reset":true}`) — this
zeroes the GRU hidden state and the odometry/history together.

### Files To Keep Aligned (one-leg)

- `source/TwoWheeledRobot/TwoWheeledRobot/tasks/direct/twowheeledrobot/one_leg_balance_env.py`
- `source/TwoWheeledRobot/TwoWheeledRobot/tasks/direct/twowheeledrobot/one_leg_balance_env_cfg.py`
- `source/TwoWheeledRobot/TwoWheeledRobot/tasks/direct/twowheeledrobot/agents/rsl_rl_one_leg_balance_cfg.py`
- `scripts/uart_policy_runner_one_leg.py`

When changing the observation layout, action scaling, joint/current constants,
or the GRU shape, update `uart_policy_runner_one_leg.py` and the STM32 inference
code in the same commit.
