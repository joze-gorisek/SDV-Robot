"""RSL-RL PPO configuration for the one-leg balance controller.

Recurrent (GRU) policy variant
------------------------------
The actor/critic are now recurrent: a single-layer GRU (hidden dim 32) feeds a
[32, 32] MLP head.  The GRU carries a hidden state between control steps, so it
can integrate the observation history (motor lag, the randomized 1-step
observation delay, wheel/joint dynamics) instead of relying on the stacked
`prev_*` observations to expose history.  `handle_deprecated_rsl_rl_cfg()` in
train.py/play.py converts this deprecated `RslRlPpoActorCriticRecurrentCfg`
into the new `RslRlRNNModelCfg` (class_name "RNNModel") for rsl-rl >= 4.0.

STM32 deployment impact — the firmware contract changes:
  * The exported policy is now STATEFUL. `play.py` exports it with the GRU
    hidden state as an extra input/output (the exporter has a dedicated GRU
    path), so the STM32 must persist a 32-float hidden state across the 50 Hz
    control loop and feed it back each step.
  * Hidden state must be zeroed whenever balancing (re)starts on hardware —
    there is no natural "episode reset" in the field, so define one (e.g. on
    entering the one-leg stance).
  * `scripts/uart_policy_runner.py` and `STM32_DEPLOYMENT.md` describe a
    stateless MLP and must be updated in lockstep before hardware use.
  * Parameter/flash budget grows (GRU ≈ 4.4 K params + MLP head ≈ 2.2 K ≈
    6.6 K params, ~26 KB float32) but still fits the F446RE (512 KB flash,
    128 KB SRAM) comfortably; per-step compute is trivial at 50 Hz.
"""

from isaaclab.utils import configclass
from isaaclab_rl.rsl_rl import (
    RslRlOnPolicyRunnerCfg,
    RslRlPpoActorCriticRecurrentCfg,
    RslRlPpoAlgorithmCfg,
)


@configclass
class OneLegBalancePPORunnerCfg(RslRlOnPolicyRunnerCfg):
    # num_steps_per_env must stay in sync with the env cfg's
    # curriculum_num_steps_per_env (iteration counting for the spawn-roll
    # curriculum). For the recurrent policy it also sets the BPTT sequence length.
    num_steps_per_env = 64
    # Spawn-roll curriculum (one_leg_balance_env_cfg.spawn_roll_curriculum):
    # 0–350 iters ±5°, 350–700 ±8°, 700+ ±11° around the 47° stance -- the
    # widest stage has no upper bound in spawn_roll_curriculum, it just holds
    # ±11° for every iteration past 700, so extending training here gives that
    # final stage more convergence time (700-1300 = 600 iters at ±11°, vs 350
    # for the previous two stages).
    max_iterations = 1300
    save_interval = 50
    experiment_name = "one_leg_balance_two_wheel"
    obs_groups = {"actor": ["policy"], "critic": ["policy"]}

    # Must be exactly RslRlPpoActorCriticRecurrentCfg — handle_deprecated_rsl_rl_cfg
    # dispatches on `type(policy) is ...`, so a subclass would not convert.
    policy: RslRlPpoActorCriticRecurrentCfg = RslRlPpoActorCriticRecurrentCfg(
        init_noise_std=0.2,
        noise_std_type="log",
        actor_obs_normalization=False,
        critic_obs_normalization=False,
        rnn_type="gru",
        rnn_hidden_dim=32,            # hidden state carried between control steps
        rnn_num_layers=1,             # keep it single-layer for the STM32 budget
        actor_hidden_dims=[32, 32],   # MLP head after the GRU
        critic_hidden_dims=[32, 32],
        activation="relu",
    )

    algorithm: RslRlPpoAlgorithmCfg = RslRlPpoAlgorithmCfg(
        value_loss_coef=1.0,
        use_clipped_value_loss=True,
        clip_param=0.2,
        entropy_coef=0.005,           # slightly higher than balance — more exploration needed
        num_learning_epochs=4,
        num_mini_batches=4,
        learning_rate=3.0e-4,
        schedule="adaptive",
        gamma=0.99,
        lam=0.95,
        desired_kl=0.01,
        max_grad_norm=0.5,
        normalize_advantage_per_mini_batch=True,
    )
