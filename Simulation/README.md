# Simulation (Isaac Lab + RSL-RL)

Isaac Lab extension `TwoWheeledRobot` with four registered tasks (stand-up, residual LQR,
pure-NN two-wheel balance, one-leg balance), plus training, export and analysis scripts.

```
source/TwoWheeledRobot/   Isaac Lab extension (envs, configs, robot model, USD assets)
scripts/rsl_rl/           train.py / play.py (play.py also exports TorchScript + ONNX)
scripts/                  LQR tools, UART policy runners, ONNX validation, plots, benchmarks
pretrained/               final one-leg checkpoint (model_1299.pt), params and TensorBoard log
data/                     DDSM115 free-spin measurement used for the motor model
docs/                     developer guide, STM32 policy contract, paper notes, presentation briefing
```

Setup, following the Isaac Lab installation guide, then:

```bash
python -m pip install -e source/TwoWheeledRobot
python scripts/list_envs.py
python scripts/rsl_rl/train.py --task Template-Twowheeledrobot-OneLegBalance-v0 --headless --num_envs 4096
python scripts/rsl_rl/play.py  --task Template-Twowheeledrobot-OneLegBalance-v0 --num_envs 1 \
       --checkpoint pretrained/one_leg_balance_two_wheel/model_1299.pt
```

The STM32 policy interface (observation and action layout, sign conventions, GRU state handling)
is specified in [`docs/STM32_DEPLOYMENT.md`](docs/STM32_DEPLOYMENT.md).
