#!/usr/bin/env python3

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEMES_DIR = ROOT / "demes"
FIG_DIR = ROOT / "figures"

MODELS = [
    ("constant", DEMES_DIR / "constant.yaml", "#1b9e77"),
    ("bottleneck", DEMES_DIR / "bottleneck.yaml", "#d95f02"),
    ("expansion", DEMES_DIR / "expansion.yaml", "#7570b3"),
]

FULL_X_MAX = 6000


def load_schedule(path: Path):
    with open(path, "r", encoding="utf-8") as fh:
        graph = yaml.safe_load(fh)

    demes_list = graph.get("demes", [])
    if len(demes_list) != 1:
        raise ValueError(f"{path} must contain exactly one deme.")

    epochs = demes_list[0].get("epochs", [])
    if not epochs:
        raise ValueError(f"{path} must define at least one epoch.")

    change_times = sorted(
        int(round(float(epoch["start_time"])))
        for epoch in epochs[1:]
        if "start_time" in epoch
    )
    max_time = max(change_times) if change_times else 0

    x = [max_time]
    y = [int(round(float(epochs[0]["start_size"])))]

    for idx, epoch in enumerate(epochs[1:], start=1):
        if "start_time" not in epoch or "start_size" not in epoch:
            raise ValueError(f"{path} epoch {idx} must define start_time and start_size.")
        x.append(int(round(float(epoch["start_time"]))))
        y.append(int(round(float(epoch["start_size"]))))

    x.append(0)
    y.append(y[-1])
    return x, y


def main():
    FIG_DIR.mkdir(exist_ok=True)
    out_path = FIG_DIR / "demes_models.png"
    schedules = [(label, *load_schedule(path), color) for label, path, color in MODELS]
    max_y = max(max(y) for _, _, y, _ in schedules)

    fig, axes = plt.subplots(3, 1, figsize=(8.5, 10.5), sharex=True, sharey=True)

    for ax, (label, x, y, color) in zip(axes, schedules):
        if x[0] < FULL_X_MAX:
            x = [FULL_X_MAX] + x
            y = [y[0]] + y

        ax.step(x, y, where="post", linewidth=3.0, color=color)
        ax.set_title(label.capitalize())
        ax.grid(alpha=0.18, linewidth=0.6)
        ax.set_xlim(FULL_X_MAX, 0)
        ax.set_ylim(0, max_y * 1.05)
        ax.set_xlabel("Generations Before Present")

    axes[0].set_ylabel("Population Size")
    axes[1].set_ylabel("Population Size")
    axes[2].set_ylabel("Population Size")
    fig.suptitle("ABCslim_trees Demographic Models", y=0.98)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    print(out_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
