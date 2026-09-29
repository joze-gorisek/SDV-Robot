# Project Guide

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Isaac Lab RL project for a two-wheeled balancing robot (4 CyberGear leg joints + 2 DDSM115 wheel motors), trained with RSL-RL PPO and deployed to an STM32F446RE over UART. Four registered tasks (see `source/TwoWheeledRobot/TwoWheeledRobot/tasks/direct/twowheeledrobot/__init__.py`):

| Task ID | Env | Purpose |
|---|---|---|
| `Template-Twowheeledrobot-Standup-v0` | `StandupEnv` | Self-right from fallen poses |
| `Template-Twowheeledrobot-ResidualLQR-v0` | `ResidualLqrEnv` | RL residual on top of LQR balance |
| `Template-Twowheeledrobot-PureNNBalance-v0` | `PureNNBalanceEnv` | Pure NN two-wheel balance |
| `Template-Twowheeledrobot-OneLegBalance-v0` | `OneLegBalanceEnv` | Balance on left wheel at ~47° roll |

## Runtime Environment

The user runs Isaac Sim **inside a Docker container** where this repo is mounted at `/workspace/Project` and Isaac Lab at `/workspace/isaaclab`. On the host, Isaac Lab is at `/home/jaka/IsaacLab`. All train/play commands go through the Isaac Lab launcher:

```bash
# inside container (paths as user sees them):
isaaclab/isaaclab.sh -p Project/scripts/rsl_rl/train.py \
    --task Template-Twowheeledrobot-OneLegBalance-v0 --headless --num_envs 8000

isaaclab/isaaclab.sh -p Project/scripts/rsl_rl/play.py \
    --task Template-Twowheeledrobot-OneLegBalance-v0 --num_envs 4
```

- GPU is 8 GB — ~4000–8000 envs fits; 10000 OOMs if a stale run is still resident. Stale python3 processes from crashed runs must be killed before relaunching.
- Play without `--headless` needs `--/app/extensions/registryEnabled=false` if the container has no internet (extension registry sync fails otherwise).
- Checkpoints and TensorBoard logs land in `logs/rsl_rl/<experiment_name>/<timestamp>/`; `play.py` auto-exports `exported/policy.pt` (TorchScript) and `policy.onnx`. Monitor with `tensorboard --logdir logs/rsl_rl/<experiment_name>` — env metrics (`roll_abs_deg`, `x_rel_abs`, `current_rms`, per-term reward penalties…) are logged via `extras["log"]` under the `Episode/` prefix.
- **Logs die with the container**: `container.py stop` runs `docker compose down --volumes`, which deletes the container's writable layer (`/workspace/logs`). Symlink logs onto the host mount before training so they survive: `mkdir -p /workspace/Project/logs && ln -sfn /workspace/Project/logs /workspace/logs` (re-run after each new container; the data lands on the host SSD).
- **Pause/resume** (`train.py`): rsl-rl auto-saves `model_<iter>.pt` every `save_interval`. `train.py` also catches SIGINT/SIGTERM (Ctrl+C, `kill`, `docker stop`) and saves a checkpoint at the exact iteration before exiting. Resume with `--resume --load_run <run_folder>`: it continues in the SAME run folder, restores the iteration count, and treats `max_iterations` as the TOTAL target (not additional). `train.py` also calls `env.set_curriculum_iteration_offset()` if the env defines it (`hasattr`-guarded; `OneLegBalanceEnv` implements it — the env derives its iteration from `common_step_counter`, which restarts at 0 in a fresh process, so without the offset a resumed run would fall back to curriculum stage 0).
- `scripts/list_envs.py` prints registered task IDs; `scripts/test_one_leg_tipping.py` (static tipping angle) and `scripts/test_one_leg_jump.py` (open-loop jump-launch tuning, bypasses env.step) are the standalone one-leg physics diagnostics; the play loop prints per-100-step CyberGear torque maxima and wheel current (`robot.data.applied_torque` — there is no `joint_effort` attribute).

## Hard Constraints

- **Network size is [32, 32] — non-negotiable.** The policy must fit STM32F446RE flash/SRAM (~6 KB of float32 weights). Never widen the hidden layers.
- **No observation normalizers** (`actor_obs_normalization=False`). Observations are manually normalized in the env via fixed `observation_scale` tuples so the exact same arithmetic can be replicated in fixed-point firmware. Changing an observation layout breaks checkpoint compatibility AND the firmware contract — update `STM32_DEPLOYMENT.md` and `scripts/uart_policy_runner.py` in lockstep.
- **rsl-rl >= 4.0 API**: `train.py`/`play.py` call `handle_deprecated_rsl_rl_cfg()` to convert the deprecated `policy: RslRlPpoActorCriticCfg` style configs. Every agent cfg needs `obs_groups = {"actor": ["policy"], "critic": ["policy"]}`. In rsl-rl 4.x the trained module is `runner.alg.actor` (not `.policy` / `.actor_critic`).
- Do not import `DistillationRunner` at module level in play/train scripts — it pulls `tensordict`, which crashes against Isaac Sim's bundled torch. Import it lazily if ever needed.

## Environment Inheritance Chain

```
DirectRLEnv
  └── StandupEnv            (joint indexing, DDSM115 motor model, CG limits, resets)
        ├── ResidualLqrEnv  (LQR + residual actions)
        └── PureNNBalanceEnv (CurrentActionProcessor, disturbances, pitch obs/noise/delay)
              └── OneLegBalanceEnv (roll control, all 4 joints policy-controlled, x_rel odometry)
```

Subclass envs call `super()._pre_physics_step()` / `super()._reset_idx()` and then **override** the parts they change (e.g. `OneLegBalanceEnv` re-writes the root pose after the parent reset). Shared building blocks (`CurrentActionProcessor`, `DisturbanceGenerator`, `pitch/roll_from_projected_gravity`, `BalanceReward`) live in `pure_nn_components.py`.

## Robot Conventions (easy to get wrong)

- **CyberGear joint order** everywhere (`_cg_ids`, obs, actions): `[front_left, front_right, back_left, back_right]`.
- **Extension sign** `_cg_ext_sign = [+1, -1, -1, +1]` — mirrored joints; multiply by it so "positive = more extended" on every leg. Raw joint limits: FL/BR = [-10°, +90°], FR/BL = [-90°, +10°].
- **Wheel order** `_wheel_ids = [left, right]` (USD joints `DDSM115_Levi`, `DDSM115_Desni`), forward-positive sign `_wheel_sign = [-1, +1]`.
- **Pitch** = `atan2(grav_y, -grav_z)`, **roll** = `atan2(grav_x, -grav_z)` from body-frame projected gravity (`pure_nn_components.py`).
- Physics timing: `PHYSICS_DT = 1 ms × decimation 20 = 20 ms` control step (50 Hz) — defined in `sim_params.py`, overridable per-env cfg.
- CyberGear joints are position-controlled (kp/kd); wheels are current/effort-controlled through the DDSM115 motor model.

## One-Leg Balance Task Specifics

- 6D obs (`roll, roll_rate, pitch, pitch_rate, x_rel, velocity`) / 5D action (`[wheel_L, joint_FL, joint_BL, joint_FR, joint_BR]`); right wheel lifted (zero current). All four CG joints are policy-controlled — bottom (FL/BL) narrow range for balance, top (FR/BR) wide range as tip-up actuators. Two-phase reward (reach→balance) latched by `_has_reached` at `phase2_roll_deg`=40°. `_reset_idx` **seeds** the latch from the sampled spawn roll rather than zeroing it, so most envs start already in phase 2. Phase 1 is still reachable at the wider curriculum bands: spawn is 47° ± band, so ±8° reaches 39° and ±11° reaches 36° — both dip below the 40° latch. Phase-1-only logic (`_wrong_direction`, the low-roll termination gate) therefore still matters from iteration ~350 on. Action/current history is not observed — the GRU hidden state carries it.
- **Recurrent policy**: this task uses a GRU (`RslRlPpoActorCriticRecurrentCfg`, `rnn_type="gru"`, hidden 32, 1 layer) + [32,32] MLP head — the only task that isn't a plain MLP. The exported policy is **stateful** (GRU hidden state persists between 50 Hz steps and must be zeroed on stance entry). Deploy via `scripts/uart_policy_runner_one_leg.py` (TorchScript, internal state + `reset()`) or the ONNX export (external `h_in`/`h_out`); the standup `uart_policy_runner.py` is stateless and separate. The [32,32] hidden-layer cap still holds.
- Robot spawns **already in the one-leg stance**: roll = `spawn_roll_center_deg` (47°) ± a half-band that widens with the training iteration via `spawn_roll_curriculum` = `((0, 5°), (350, 8°), (700, 11°))`; pitch spawns at 0°. The iteration is derived as `common_step_counter // curriculum_num_steps_per_env` — that cfg value (64) **must match `num_steps_per_env` in `rsl_rl_one_leg_balance_cfg.py`** or the stages mis-time silently. `set_spawn_roll_deg(center, band)` is the runtime override for eval sweeps (the cached `_spawn_roll_*` fields are read, not `cfg`, so mutating the cfg post-construction does nothing); call `env.reset()` after it. The spawn height formula (`spawn_upright_z`) is calibrated for folded (raw-0) legs, not the extended `reset_joint_ext_deg` pose used here, so `_reset_idx` adds a fixed `spawn_lateral_h` safety margin to guarantee the robot never spawns clipped into the ground; it settles the last few mm onto the ground under gravity within the first control step.
- Termination: roll outside [20°, 72°] or |pitch| > 40°, both for 5 consecutive steps.
- `x_rel` (obs[10]) is integrated from **left-wheel odometry only** (`_compute_velocity()`), mirroring what the STM32 computes from the DDSM115 encoder. Ground speed must use the left wheel only — the two-wheel average from the parent env underestimates by ~2× here.
- Reward-weight lesson learned: velocity/position penalties must stay well below pitch/roll-rate weights, or the policy refuses the corrective wheel motion needed to catch falls (survival ≫ balance ≫ station-keeping).
- Typical training: ~500–1000 iterations at 4000–8000 envs; reward converges around iteration 300–400.

## Lint & Diagnostics

- Lint/format: `ruff check .` and `ruff format .` (config in `pyproject.toml`); `pre-commit run --all-files` runs ruff plus codespell, license-header, and whitespace hooks. Type checking: `pyright` (basic mode, targets `source` and `scripts`).
- There is no pytest suite — correctness is checked via standalone diagnostic scripts run inside the Isaac Lab container, e.g. `scripts/list_envs.py` (prints registered task IDs), `scripts/test_one_leg_tipping.py` (one-leg physics diagnostic), `scripts/test_policy_angle_sweep.py` (sanity-checks an exported policy's outputs across synthetic poses without hardware), `scripts/validate_onnx_policy.py` (ONNX vs. TorchScript parity).
- `scripts/lqr_control.py` is the LQR diagnostic/evaluation entry point (free-spin motor test, `lqr-model`, `lqr-floor` modes); `scripts/calculate_lqr_gains.py` prints reference LQR gains.
- `scripts/plot_training_curves.py` renders publication-quality training-curve figures from a TensorBoard event file (reads `Episode/` tags, no Isaac Sim needed, plain `python`); `scripts/lqr_control _one_leg_jump.py` is the one-leg-jump analogue of `lqr_control.py` (open-loop jump/LQR diagnostics for the one-leg task).

## Other Docs Worth Reading

- `DEVELOPER_GUIDE.md` — file-by-file "edit here first" table, but it predates `PureNNBalance`/`OneLegBalance` and only documents Standup/ResidualLQR; treat it as partial, not authoritative for the newer tasks.
- `STM32_DEPLOYMENT.md` — firmware policy contract (obs order, normalization, action scaling) for all four tasks, including the stateful GRU contract for `OneLegBalance`.
- `RealImplementationCode/` — STM32 firmware-side files.
