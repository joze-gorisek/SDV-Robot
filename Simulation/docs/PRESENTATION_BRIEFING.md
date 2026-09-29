# Project Briefing — Two-Wheeled Wheel-Legged Balancing Robot

**Purpose of this document:** a self-contained context dump to paste into a fresh chat
that will help build a **student paper competition presentation**. Everything below is
drawn from the actual repository state and the real training logs — no invented numbers.
Where something is *not* yet demonstrated, it is marked explicitly. Do not let a deck
claim more than the "Evidence status" section supports.

Repo: `TwoWheeledRobot` · working branch `JakaERK2026v2` (suggests **ERK 2026**, the
Slovenian Electrotechnical and Computer Science Conference — student paper track).

---

## 1. One-paragraph summary

A two-wheeled, wheel-legged balancing robot (four CyberGear leg joints + two DDSM115 hub
motors) is trained in Isaac Lab with RSL-RL PPO and deployed to an STM32F446RE
microcontroller over UART. The distinguishing constraint is that every policy must fit a
microcontroller's flash and SRAM budget and be reproducible in fixed-point firmware
arithmetic — so the network is capped at `[32, 32]` hidden units and **no learned
observation normalizer is allowed** (observations are normalized by hand-fixed constants
inside the environment so the firmware can replicate the exact arithmetic). The newest and
most novel result is a **recurrent (GRU) policy that balances the robot on a single wheel**
at roughly 47° of roll, holding the stance for full 10-second episodes with a measured
fall rate of zero in simulation.

---

## 2. Hardware and physical parameters

| Quantity | Value | Source |
|---|---|---|
| Leg actuators | 4 × CyberGear, position-controlled (kp/kd) | `robot_cfg.py` |
| Wheel actuators | 2 × DDSM115 hub motors, current/torque-controlled | `robot_cfg.py` |
| Microcontroller | STM32F446RE (512 KB flash, 128 KB SRAM) | `rsl_rl_one_leg_balance_cfg.py` |
| IMU | BNO080-class (projected gravity + body angular rate) | env code |
| Body mass | 2.60 kg | `residual_lqr_env.py:83` |
| Wheel/cart mass | 1.53 kg | `residual_lqr_env.py:84` |
| Wheel radius | 0.05035 m | `residual_lqr_env.py:85` |
| Body CoM height | 0.14 m | `residual_lqr_env.py:86` |
| Body pitch inertia | 0.01296 kg·m² | `residual_lqr_env.py:87` |
| Body yaw inertia | 0.03151 kg·m² | `residual_lqr_env.py:88` |
| Track width | 0.383 m | `residual_lqr_env.py:89` |
| Motor torque constant `K_t` | 0.75 Nm/A | `sim_params.py:83` |
| Continuous / peak current | 1.5 A / 2.7 A | `sim_params.py:84–85` |
| Rated / peak torque | 0.96 Nm / 2.0 Nm | `sim_params.py:86–87` |
| No-load speed | 200 rpm | `sim_params.py:89` |
| Ground friction (static / dynamic) | 0.6 / 0.4 | `sim_params.py:65–66` |
| Control rate | 50 Hz (1 ms physics × decimation 20) | `sim_params.py:53–54` |

**Robot conventions that are easy to get wrong** (worth a backup/appendix slide):

- CyberGear joint order everywhere: `[front_left, front_right, back_left, back_right]`.
- Extension sign `[+1, −1, −1, +1]` — the mirrored legs must be multiplied by this so
  "positive = more extended" holds on every leg. Raw limits: FL/BR `[−10°, +90°]`,
  FR/BL `[−90°, +10°]`.
- Wheel order `[left, right]`, forward-positive sign `[−1, +1]`.
- `pitch = atan2(g_y, −g_z)`, `roll = atan2(g_x, −g_z)` from body-frame projected gravity.

---

## 3. The four tasks (the project arc)

All are registered Gym tasks in
`source/TwoWheeledRobot/TwoWheeledRobot/tasks/direct/twowheeledrobot/__init__.py`.

| Task ID | Env class | Purpose |
|---|---|---|
| `Template-Twowheeledrobot-Standup-v0` | `StandupEnv` | Self-right from fallen poses |
| `Template-Twowheeledrobot-ResidualLQR-v0` | `ResidualLqrEnv` | RL residual added on top of an analytical LQR |
| `Template-Twowheeledrobot-PureNNBalance-v0` | `PureNNBalanceEnv` | Pure NN two-wheel balance (no LQR) |
| `Template-Twowheeledrobot-OneLegBalance-v0` | `OneLegBalanceEnv` | **Balance on the left wheel at ~47° roll** |

Inheritance chain (each subclass calls `super()` then overrides what it changes):

```
DirectRLEnv
  └── StandupEnv              joint indexing, DDSM115 motor model, CG limits, resets
        ├── ResidualLqrEnv    LQR + residual actions
        └── PureNNBalanceEnv  CurrentActionProcessor, disturbances, pitch obs/noise/delay
              └── OneLegBalanceEnv   roll control, all 4 joints, x_rel odometry, GRU
```

Shared building blocks (`CurrentActionProcessor`, `DisturbanceGenerator`,
`pitch/roll_from_projected_gravity`, `BalanceReward`) live in `pure_nn_components.py`.

**Narrative value:** this is a clean escalation — classical LQR baseline → RL correcting
LQR → RL replacing LQR → RL doing something LQR cannot do at all (the one-leg stance).
That is a strong spine for a talk.

---

## 4. The headline contribution: one-leg balance with a recurrent policy

### 4.1 Why it is interesting

The robot lifts its right wheel and balances on the **left wheel alone**, tilted ~47° in
roll. This is simultaneously an inverted pendulum in pitch and a tipping problem in roll,
with only one ground contact point. The right wheel must stay off the ground — enforced by
a real PhysX contact sensor, not a joint-angle proxy. The controller is pure RL; no LQR
assist.

### 4.2 Observation / action contract (6D → 5D)

Observation (normalized by fixed constants, then `nan_to_num`, then clamped to `[−10, 10]`):

| # | Signal | Normalization | Physical source |
|---|---|---|---|
| 0 | `roll` | / 50° | `atan2(g_x, −g_z)` |
| 1 | `roll_rate` | / 4.0 rad/s | `gyro_y` |
| 2 | `pitch` | / 25° | `atan2(g_y, −g_z)` |
| 3 | `pitch_rate` | / 4.0 rad/s | `−gyro_x` |
| 4 | `x_rel` | / 0.5 m | **left-wheel odometry only**, integrated since stance entry |
| 5 | `velocity` | / 0.5 m/s | **left wheel only** |

Action (5 raw outputs in `[−1, 1]`):

```
0  wheel_L   → left wheel current = tanh(a0) · I_max        (I_max = 2.0 A)
1  joint_FL  → ext_sign_FL · (45° + a1 · 45°)    bottom — fine balance leg
2  joint_BL  → ext_sign_BL · (45° + a2 · 45°)    bottom — fine balance leg
3  joint_FR  → ext_sign_FR · (45° + a3 · 45°)    top — tip-up actuator
4  joint_BR  → ext_sign_BR · (45° + a4 · 45°)    top — tip-up actuator
```

Right wheel current is always exactly `0.0` (it is in the air).

**Critical detail worth a slide:** `x_rel` and `velocity` must come from the left wheel
*only*. The parent environment's two-wheel average would underestimate ground speed by
roughly 2× here, because the lifted right wheel contributes zero. This mirrors exactly what
the STM32 computes from the single DDSM115 encoder that is actually turning.

### 4.3 Network architecture — and why it is recurrent

This is the only task in the project that is **not** a plain MLP:

```
obs(6) → GRU(hidden 32, 1 layer) → MLP [32, 32] → actions(5)
```

- `rnn_type="gru"`, `rnn_hidden_dim=32`, `rnn_num_layers=1`, `activation="relu"`.
- ≈ 4.4 K GRU parameters + ≈ 2.2 K MLP head ≈ **6.6 K parameters ≈ 26 KB as float32**.
- Fits the F446RE (512 KB flash / 128 KB SRAM) comfortably; per-step compute is trivial at 50 Hz.

**Why a GRU:** the previous wheel currents and previous joint commands are deliberately
*not* in the observation vector. The hidden state carries that history internally instead.
This is the substantive design claim — the recurrent state replaces the stacked
`prev_action` observations used by the other tasks, which keeps the observation vector down
to 6 physical quantities the firmware already has.

**The deployment consequence:** the exported policy is **stateful**. The 32-float hidden
state must persist across 50 Hz control steps and be fed back each step, and it must be
**zeroed every time the robot re-enters the one-leg stance** — there is no natural "episode
reset" in the field. Failing to reset leaks stale history from the previous run and degrades
balance. This is a genuinely interesting sim2real wrinkle and is a good discussion point.

### 4.4 Training setup

| Parameter | Value |
|---|---|
| Algorithm | PPO (RSL-RL ≥ 4.0) |
| `num_steps_per_env` | 64 (also the BPTT sequence length) |
| `max_iterations` | 1300 |
| Learning rate | 3.0e-4, adaptive schedule |
| `entropy_coef` | 0.005 |
| `gamma` / `lam` | 0.99 / 0.95 |
| `clip_param` | 0.2, `desired_kl` 0.01 |
| Epochs / minibatches | 4 / 4 |
| `init_noise_std` | 0.2 |
| Episode length | 10.0 s = **500 control steps at 50 Hz** |
| Envs | ~4000–8000 (8 GB GPU; 10000 OOMs) |
| Observation normalizer | **disabled** (`actor_obs_normalization=False`) |

**Spawn-roll curriculum.** The robot spawns *already in* the one-leg stance at
roll = 47° ± band, where the band widens with the training iteration:

```
iterations   0–350   →  ±5°
iterations 350–700   →  ±8°
iterations   700+    →  ±11°
```

The iteration is derived as `common_step_counter // curriculum_num_steps_per_env`, and that
config value (64) **must match `num_steps_per_env`** in the agent config or the stages
mis-time silently.

**Two-phase reward with a latch.** Phase 1 ("reach") pulls roll up toward the 47° target;
phase 2 ("balance") takes over once roll first crosses 40°, latched one-way per environment
for the rest of the episode. `_reset_idx` *seeds* the latch from the sampled spawn roll
rather than zeroing it, so most environments begin already in phase 2. Phase 1 still matters
at the wider curriculum bands (47° − 8° = 39° and 47° − 11° = 36° both dip below the 40°
latch), so the phase-1 logic remains live from iteration ~350 onward.

Reward weights (phase 2 unless noted):

| Term | Weight |
|---|---|
| Alive | +1.0 |
| Reach, `−(roll − target)²` (phase 1) | 3.5 |
| Reach rate, `+roll_rate` (phase 1) | 1.5 |
| Roll rate | 2.0 |
| Pitch | 5.0 |
| Pitch rate | 0.3 |
| Velocity | 0.4 |
| Position (`x_rel` drift) | 0.15 |
| Joint action | 0.1 |
| Joint smoothness | 2.0 |
| Wheel current | 0.003 |
| Fall / physics-broken / invalid-state penalty | −100.0 |

**Reward-tuning lesson that is genuinely worth presenting:** velocity and position penalties
must stay well below the pitch and roll-rate weights. An inverted pendulum has to drive its
wheel *under* its centre of mass to catch a fall; heavy velocity or position penalties
suppress exactly that recovery motion. The priority order is
**survival ≫ balance ≫ station-keeping**. There is a second documented instance of the same
class of problem: at `rew_reach = 2.0` the policy found a profitable local optimum by camping
at its ~15° static leg-extension ceiling (`alive(+1.0) − 2.0·(15°−47°)² ≈ +0.38/step`);
steepening to 3.5 made that camp cost ≈ −0.09/step and forced it to climb, while still being
far better than the −100 of falling, so it did not learn to self-terminate instead. These are
concrete, quantitative reward-shaping anecdotes — much more compelling than a generic
"we tuned the reward" line.

**Termination conditions:**

- Roll outside `[20°, 72°]` for 5 consecutive steps (the low bound only applies after the
  stance has been reached).
- `|pitch| > 40°` for 5 consecutive steps.
- Right-wheel contact force > 0.1 N for 5 consecutive steps (PhysX contact sensor).
- Phase 1 only: roll below −30° for 5 consecutive steps ("wrong direction").
- Timeout at 500 steps.

### 4.5 Results — measured, from the actual TensorBoard log

Run: `logs/rsl_rl/one_leg_balance_two_wheel/2026-07-16_13-34-22`, 1300 iterations.
(Extracted directly from the event file; the tags live under the `Episode/` prefix, written
via `extras["log"]`.)

**Final performance, mean over the last 50 iterations:**

| Metric | Value |
|---|---|
| Mean episode length | **496.2 / 500 steps** (≈ 99.2 % of the 10 s episode) |
| Fall rate (`roll_fall_rate`) | **0.000** |
| Terminations by timeout | 0.002 |
| Terminations by right-wheel touch | 0.000 |
| Mean total reward | **349.0** |
| Mean per-step reward | 0.705 (of a 1.0 alive maximum) |
| Held roll angle | 53.37° |
| Mean absolute pitch | **2.83°** |
| Mean absolute roll rate | 0.134 rad/s |
| Position drift `x_rel` | **0.082 m** |
| Wheel current RMS | 0.245 A (of a 2.0 A envelope) |
| "Reached stance" rate | 99.9 % |

**Learning curve (use this for the results figure):**

| Iteration | Mean reward | Episode length | Fall rate | Pitch (°) | Roll rate | `x_rel` (m) |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | −208.9 | 60.4 | 0.017 | 10.55 | 0.491 | 0.091 |
| 50 | −245.2 | 62.3 | 0.031 | 7.30 | 0.518 | 0.102 |
| 100 | −159.5 | 92.1 | 0.005 | 9.84 | 0.348 | 0.109 |
| 150 | +49.5 | 343.8 | 0.000 | 5.06 | 0.187 | 0.143 |
| 200 | 174.6 | 430.0 | 0.000 | 4.02 | 0.165 | 0.189 |
| 300 | 309.9 | 499.0 | 0.000 | 3.21 | 0.142 | 0.135 |
| 500 | 328.7 | 499.0 | 0.000 | 2.93 | 0.139 | 0.111 |
| 700 | 301.2 | 454.9 | 0.000 | 3.04 | 0.137 | 0.099 |
| 1000 | 320.5 | 487.7 | 0.000 | 2.86 | 0.150 | 0.086 |
| 1299 | 337.1 | 492.0 | 0.000 | 2.83 | 0.129 | 0.078 |

**How to read this narratively:** iterations 0–100 are the exploration phase (negative
reward, ~60-step episodes, the policy is falling). The transition happens sharply between
iteration 100 and 150 — episode length jumps from 92 to 344 and reward crosses zero. By
iteration 300 the policy already survives full episodes with zero falls. Everything after
that is refinement: pitch tightens from 3.2° to 2.8°, roll rate from 0.142 to 0.129,
position drift from 0.135 m to 0.078 m. The dip around iteration 700 lines up with the
curriculum widening to ±11°, and the policy recovers from it — that is a nice, honest detail
to point out rather than hide.

Note: `Policy/mean_std` falls from 0.199 to 0.032, i.e. the policy becomes nearly
deterministic — consistent with a converged, confident controller.

**Throughput:** ≈ 21 000–22 700 total FPS; ~11–12 s per iteration collection time;
1300 iterations ran in roughly 4 hours wall-clock (15:34 → 19:39 by checkpoint timestamps).

---

## 5. Supporting work (candidate appendix / context material)

### 5.1 LQR baseline

An analytical LQR was derived from the linearized pitch model and validated in stages.
Gain matrix: `K = [−27.83, −6.05, −1.00, −3.21]`.

- **Suspended-air test** (`LQR_MODEL_SUSPENDED_TEST_REPORT.md`): sign convention verified —
  `+1°` pitch → `+0.486 A`, `−1°` → `−0.486 A`; current routed correctly through the DDSM115
  torque/saturation/torque-speed model; wheels reached ≈ 199 rpm against a modeled 200 rpm
  no-load speed. All acceptance checks passed. No balancing performance was claimed, correctly,
  because the wheels had no ground contact.
- **First floor-contact test at 0.5° initial pitch:** the motors-disabled baseline hit the 10°
  stop threshold in 0.30 s; the LQR-enabled run completed the full 5.00 s window with pitch
  bounded to `[−0.49°, +2.29°]`. However it oscillates, saturates current in 30.1 % of samples
  (peak request 6.54 A against a 2.7 A clamp), and drifts 0.389 m. **No Q/R tuning has been
  performed.** Be careful not to present the LQR as a tuned, fair baseline — it is a validated
  but untuned reference.

### 5.2 Pure NN two-wheel balance

Documented in detail in `paper_section_3_4_implementation_notes.md` (this file is already
close to paper prose and is the best single source for methods text).

- Architecture `8 → 32 → 32 → 2`, ReLU, **1410 parameters**.
- Observations: `x_rel`, velocity, pitch, pitch rate, yaw error, yaw rate, previous left and
  right currents.
- Actions: two wheel currents via `tanh(a) · I_max`, `I_max = 2.0 A`.
- Reward: quadratic balance objective — alive +1.0, penalties on pitch (12.0), pitch rate (0.8),
  velocity (0.15), position (0.4), yaw error (0.1), current (0.01), current change (0.03).
- 5-stage curriculum (1000/1000/1500/1500/2000 iterations) introducing disturbances progressively.

### 5.3 Sim2real modeling — the domain randomization story

This is a genuine strength of the work and deserves a slide. Implemented and *active*:

| Effect | Range |
|---|---|
| Left/right motor gain asymmetry | independently sampled from (0.8, 1.2) |
| Motor current deadzone | (0.03, 0.20) A |
| Motor current bias | (−0.08, 0.08) A |
| Current-loop (electrical) lag | first-order, τ ∈ (0.005, 0.020) s |
| Per-motor current saturation | (1.6, 2.4) A |
| Observation delay | 0 or 1 sample |
| Action delay | 0 or 1 sample |
| Pitch sensor bias | ±0.5° |
| Pitch noise | 0.15° std |
| Pitch-rate noise | 0.02 rad/s std |
| Velocity noise | 0.01 m/s std |
| Torque-speed limiting | linear derate to zero at 200 rpm |

**Be honest about what is configured but NOT applied** (present in the config but with no
active implementation in the environment): body mass scale, CoM height scale, pitch inertia
scale, wheel radius scale, and ground friction randomization. If a reviewer asks "did you
randomize mass?", the answer is no — the ranges exist in config but are inert. Have this ready
as a backup slide rather than being caught by it.

Disturbances during pure-NN training (curriculum stages 3–5): short "human push" force
impulses (`Fx ∈ [−6, 6] N`, 0.05–0.20 s), persistent payload / CoM-shift equivalent pitch
torques (`Mx ∈ [−0.25, 0.25] Nm`), and combinations. A sine disturbance exists only as a
diagnostic benchmark, not in the training curriculum. **Note:** the one-leg task's
`disturbance_force_n` logged exactly 0.000 throughout the run above — external disturbances
were *not* active during one-leg training. Do not claim push-robustness for the one-leg policy.

### 5.4 STM32 deployment path

- Train → `play.py` auto-exports `exported/policy.pt` (TorchScript) and `policy.onnx`.
- TorchScript keeps the GRU hidden state *inside* the module and exposes `reset()`.
- ONNX exposes it *externally* — this is the form for hand-written STM32 C:
  ```
  inputs:  obs [1,6],  h_in  [1,1,32]
  outputs: actions [1,5],  h_out [1,1,32]
  ```
- Host-side runner `scripts/uart_policy_runner_one_leg.py` speaks JSON lines over UART at
  115200 baud; a bare `reset` line zeroes the GRU state, the odometry, and the history together.
- `scripts/validate_onnx_policy.py` checks ONNX vs. TorchScript parity, failing above 1e-4.
- Firmware-side source lives in `RealImplementationCode/` (`controler.c`, `cybergear.c`,
  `DDSM115.c`, `StateEstimator.c`, `state_machine.c`, `JumpStrategy.c`, `StartupStrategy.c`,
  `telemetry.c`, …).

---

## 6. Evidence status — read this before writing any claim

**Demonstrated:**
- One-leg balance policy trains to convergence in simulation with zero falls and 99.2 % episode
  survival. Numbers in §4.5 are real and reproducible from the committed event file.
- LQR sign convention and motor-model path validated (suspended + one 0.5° floor test).
- ONNX/TorchScript parity tooling exists and enforces 1e-4.

**NOT demonstrated — do not claim:**
- **No hardware results.** There are no measured real-robot one-leg balance results anywhere
  in the repo. The STM32 path is a fully specified *contract*, plus host-side runner code —
  not a validated deployment.
- **No exported policy artifacts currently exist.** A search for `*.onnx` and `policy*.pt`
  across `logs/rsl_rl/` returns nothing; every run's `exported/` directory is empty. Run
  `play.py` on `model_1299.pt` before claiming an export exists or quoting an on-device
  footprint as measured.
- **No disturbance rejection for the one-leg task.** `disturbance_force_n` was identically zero
  for the whole run.
- **No LQR-vs-RL comparison for the one-leg stance.** LQR was never tuned, and the one-leg
  stance has no LQR baseline at all. Any "RL beats LQR" framing would be unsupported.
- **No benchmark sweep for one-leg.** `benchmark_pure_nn_balance.py` targets the two-wheel
  pure-NN task only.
- There is no pytest suite; correctness is checked through standalone diagnostic scripts.

**Honest framing that still lands well:** "a compact recurrent policy learns a one-wheel
balancing stance that has no closed-form controller in this project, converges in ~300 PPO
iterations, and fits a 26 KB microcontroller budget — with a fully specified stateful
deployment contract; hardware validation is in progress." That is a legitimate and
competitive student-paper claim.

---

## 7. Figure and exhibit sources

- `scripts/plot_training_curves.py` — renders publication-quality curves straight from a
  TensorBoard event file. Plain `python`, no Isaac Sim needed. Reads `Episode/` tags, supports
  EMA smoothing, vertical phase markers, and PDF/PNG/SVG output. **Note:** its axis labels are
  currently Slovenian (`iteracija`). Use it for the reward, episode-length, pitch, and roll-rate
  curves. Curriculum transitions to mark: `--phases 350 700`.
  - Caveat: the host Python has no `tensorboard` module installed and `pip` is
    externally-managed (PEP 668), so this script cannot run on the host as-is. Either run it
    inside the container, or use a minimal event-file parser — the numbers in §4.5 were
    extracted that way.
- `logs/rsl_rl/one_leg_balance_two_wheel/2026-07-16_13-34-22/events.out.tfevents.*` — the
  source of truth for every metric in §4.5.
- `scripts/test_one_leg_tipping.py` — static tipping-angle diagnostic; good for a figure
  motivating *why* ~47° is the equilibrium.
- `scripts/test_one_leg_jump.py`, `scripts/lqr_control_one_leg_jump.py` — open-loop jump-launch
  diagnostics.
- `scripts/test_policy_angle_sweep.py` — sanity-checks an exported policy across synthetic poses
  without hardware; could produce a nice "policy response surface" figure.
- `logs/render/*.usd` — recorded simulation renders, a source for stance screenshots.
- `Integrated_Modeling_and_Control_Optimization_of_Biped_Wheel-Legged_Robot.pdf` — reference
  paper in the repo root (external work; cite properly if used).

**Figures that do not exist yet and would have to be built:** controller block diagram,
training-pipeline diagram, the one-leg stance geometry sketch, and the sim2real
domain-randomization diagram. `paper_section_3_4_implementation_notes.md` ends with an explicit
"Figures needed" and "Tables needed" list — reuse it.

---

## 8. Suggested deck framing (for the presentation chat to refine)

**Presentation type:** Structured Argument (academic, analytical).

**Recommended single argument:** *A 6.6 K-parameter recurrent policy learns a one-wheel
balancing stance that fits a microcontroller.* The MCU constraint is what makes the recurrence
interesting — the GRU hidden state substitutes for the action-history observations the
firmware would otherwise have to buffer, which keeps the observation vector at six physical
quantities the STM32 already measures. Everything else (standup, residual LQR, pure NN balance)
becomes one context slide plus appendix.

**Rejected alternatives and why:** "RL vs LQR" is not supportable (untuned LQR, no one-leg
baseline). "Full project narrative" is too much for a competition slot and would dilute the
novel result. "Sim2real on hardware" overclaims — there are no hardware numbers.

**Proposed spine (Situation / Complication / Resolution):**

1. Wheel-legged robots need balance controllers that run on the cheap MCU already on the robot.
2. The one-wheel stance has no closed-form controller here, and the MCU budget forbids the wide
   networks that usually absorb that difficulty.
3. A GRU(32) + MLP[32,32] policy, trained with PPO under fixed hand-normalized observations,
   holds the stance for full episodes with zero falls and 26 KB of weights.

**Ghost-deck action-title draft** (these are the argument; refine before building):

1. Wheel-legged robots must balance on the microcontroller they already carry
2. The one-wheel stance has no closed-form controller — and the MCU budget rules out a big network
3. Can a 6 K-parameter recurrent policy hold a one-wheel stance at 50 Hz?
4. A GRU hidden state replaces the action history, keeping the observation at six measured quantities
5. Two-phase reward and a widening spawn-roll curriculum make the stance learnable
6. The policy converges in ~300 iterations, from falling in 1.2 s to surviving the full 10 s episode
7. At convergence the stance holds with zero falls, 2.8° pitch error and 8 cm of drift
8. Sensor noise, motor asymmetry and one-step latency are randomized so the policy transfers
9. The stateful contract is the real deployment cost: 32 floats persisted and zeroed on stance entry
10. Conclusions

**Open decisions the presentation chat must ask the user about:**

- **Language** — Slovenian or English? (ERK accepts both; existing plot labels are Slovenian.)
- **Time slot** — determines slide budget (≤ 1 slide/minute).
- **Institution template or colors**, and the venue's aspect ratio (default 16:9).
- **Author name, affiliation, co-authors/supervisor, and the paper title.**
- **Contact details and any preprint/repo link or QR code** for the conclusions slide.
- Whether hardware results will exist by the presentation date — if yes, the argument can be
  strengthened; if no, keep the honest "in progress" framing.

**Design constraints that apply** (from the academic-pptx skill): white background, one
sans-serif face, ≤ 3 colors (default navy `1F4E79` / mid-blue `2E75B6`), action titles as
complete sentences, one exhibit per results slide with the key finding annotated directly on
the chart, ≥ 20 pt body text, ≤ ~40 words per slide, a References slide, and Conclusions as the
last non-appendix slide (never "Thank You").

---

## 9. Key file map

| File | What it holds |
|---|---|
| `PROJECT_GUIDE.md` | Authoritative project guide — task table, conventions, hard constraints |
| `paper_section_3_4_implementation_notes.md` | Near-paper-prose methods for the pure-NN controller |
| `STM32_DEPLOYMENT.md` | Firmware policy contract for all four tasks, incl. the stateful GRU |
| `one_leg_balance_env.py` / `_cfg.py` | The one-leg task: rewards, terminations, curriculum |
| `agents/rsl_rl_one_leg_balance_cfg.py` | PPO + GRU architecture config, with deployment notes |
| `pure_nn_components.py` | Shared blocks: action processor, disturbances, angle helpers, reward |
| `sim_params.py` | Physics timing, ground, DDSM115 and CyberGear constants |
| `residual_lqr_env.py` | LQR physical parameters, A/B matrices, gains |
| `LQR_MODEL_SUSPENDED_TEST_REPORT.md` | LQR validation results |
| `RealImplementationCode/` | STM32 firmware sources |
| `DEVELOPER_GUIDE.md` | Partial — predates PureNNBalance and OneLegBalance |

**Hard constraints to never violate in any claim or future work:** network capped at
`[32, 32]`; no observation normalizers (manual fixed-constant normalization only, for
firmware parity); rsl-rl ≥ 4.0 API where the trained module is `runner.alg.actor`; changing an
observation layout breaks both checkpoint compatibility and the firmware contract, so
`STM32_DEPLOYMENT.md` and the UART runners must be updated in the same commit.
