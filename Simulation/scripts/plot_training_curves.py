"""Plot publication-quality training curves from a TensorBoard event file.

Reads the scalars logged by the RSL-RL runner (tags prefixed with "Episode/")
and renders clean matplotlib figures. Optional vertical lines mark curriculum
phase transitions.

Requires only `tensorboard` and `matplotlib` (no Isaac Sim). Run with plain python:

    python scripts/plot_training_curves.py \
        --logdir logs/rsl_rl/one_leg_balance_two_wheel/2026-07-05_14-38-13 \
        --tags Episode/reward Episode/roll_fall_rate \
        --phases 500 800 1100 \
        --outdir figs

With no --tags, the script lists all available scalar tags and exits.
"""

import argparse
import os

import matplotlib.pyplot as plt
from tensorboard.backend.event_processing import event_accumulator


def load_scalars(logdir: str):
    ea = event_accumulator.EventAccumulator(
        logdir, size_guidance={event_accumulator.SCALARS: 0}
    )
    ea.Reload()
    return ea


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--logdir", required=True, help="Run directory containing the event file.")
    p.add_argument("--tags", nargs="*", default=None, help="Scalar tags to plot (e.g. Episode/reward).")
    p.add_argument("--phases", nargs="*", type=float, default=None, help="Iterations to mark with vertical lines.")
    p.add_argument("--smooth", type=float, default=0.0, help="EMA smoothing factor in [0,1). 0 = raw.")
    p.add_argument("--outdir", default="figs", help="Output directory for the figures.")
    p.add_argument("--fmt", default="pdf", choices=["pdf", "png", "svg"], help="Output file format.")
    args = p.parse_args()

    ea = load_scalars(args.logdir)
    available = ea.Tags().get("scalars", [])

    if not args.tags:
        print("Available scalar tags:")
        for t in available:
            print("  ", t)
        return

    os.makedirs(args.outdir, exist_ok=True)
    plt.rcParams.update({"font.size": 11, "figure.figsize": (5.0, 3.2), "axes.grid": True})

    for tag in args.tags:
        if tag not in available:
            print(f"[skip] tag not found: {tag}")
            continue
        events = ea.Scalars(tag)
        steps = [e.step for e in events]
        values = [e.value for e in events]

        # optional exponential-moving-average smoothing
        if args.smooth > 0.0:
            sm, prev = [], values[0]
            for v in values:
                prev = args.smooth * prev + (1.0 - args.smooth) * v
                sm.append(prev)
            values = sm

        fig, ax = plt.subplots()
        ax.plot(steps, values, linewidth=1.4)
        if args.phases:
            for x in args.phases:
                ax.axvline(x, color="0.5", linestyle="--", linewidth=0.8)
        ax.set_xlabel("iteracija")
        ax.set_ylabel(tag.split("/")[-1])
        fig.tight_layout()

        safe = tag.replace("/", "_")
        out = os.path.join(args.outdir, f"{safe}.{args.fmt}")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"[saved] {out}")


if __name__ == "__main__":
    main()
