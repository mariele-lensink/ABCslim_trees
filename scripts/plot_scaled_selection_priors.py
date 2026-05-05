#!/usr/bin/env python3
"""Plot selection-coefficient priors after SLiM scaling and truncation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PARAMS = ("gdfe", "idfe")


def read_selection_priors(path: Path) -> dict[str, np.ndarray]:
    values = {param: [] for param in PARAMS}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            for param in PARAMS:
                values[param].append(float(row[param]))
    return {param: np.asarray(vals, dtype=float) for param, vals in values.items()}


def summarize(param: str, q: float, biological: np.ndarray) -> dict[str, float | str]:
    scaled = biological * q
    capped = np.maximum(scaled, -1.0)
    capped_rows = scaled < -1.0
    abs_capped = np.abs(capped)
    return {
        "parameter": param,
        "Q": q,
        "n": biological.size,
        "n_mean_capped": int(capped_rows.sum()),
        "pct_mean_capped": 100.0 * capped_rows.mean(),
        "min_abs_after_cap": abs_capped.min(),
        "q05_abs_after_cap": np.quantile(abs_capped, 0.05),
        "median_abs_after_cap": np.quantile(abs_capped, 0.50),
        "q95_abs_after_cap": np.quantile(abs_capped, 0.95),
        "max_abs_after_cap": abs_capped.max(),
    }


def write_summary(path: Path, rows: list[dict[str, float | str]]) -> None:
    fieldnames = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_q(priors: dict[str, np.ndarray], q: float, out_svg: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey="row")
    fig.suptitle(f"Selection prior after SLiM scaling and -1 truncation guardrail (Q={q:g})")

    transforms = (
        ("Biological prior", lambda x: x),
        ("Scaled in SLiM", lambda x: x * q),
        ("After -1 floor", lambda x: np.maximum(x * q, -1.0)),
    )

    for row_idx, param in enumerate(PARAMS):
        biological = priors[param]
        transformed = [func(biological) for _, func in transforms]
        logs = [np.log10(np.abs(vals)) for vals in transformed]
        x_min = min(log_vals.min() for log_vals in logs)
        x_max = max(log_vals.max() for log_vals in logs)
        bins = np.linspace(x_min, x_max, 61)

        for col_idx, ((title, _), log_vals) in enumerate(zip(transforms, logs)):
            ax = axes[row_idx, col_idx]
            ax.hist(log_vals, bins=bins, color="#3b6ea8", alpha=0.82, edgecolor="white")
            ax.set_title(title)
            ax.set_xlabel("log10(|s|)")
            if col_idx == 0:
                ax.set_ylabel(f"{param} count")
            ax.axvline(0.0, color="#9a3412", linewidth=1.2, linestyle="--")
            ax.grid(axis="y", alpha=0.25)

        scaled = biological * q
        pct_capped = 100.0 * np.mean(scaled < -1.0)
        axes[row_idx, 2].text(
            0.97,
            0.92,
            f"{pct_capped:.1f}% of prior means floored at -1",
            ha="right",
            va="top",
            transform=axes[row_idx, 2].transAxes,
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.9},
        )

    fig.tight_layout()
    fig.savefig(out_svg)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Apply ABCtrees.slim selection scaling to gdfe/idfe and plot the resulting priors."
    )
    parser.add_argument(
        "--param-file",
        type=Path,
        default=Path("priors/priors_2.6.25.csv"),
        help="Input parameter CSV with gdfe and idfe columns.",
    )
    parser.add_argument(
        "--q",
        type=float,
        nargs="+",
        default=[50.0],
        help="One or more SLiM scaling factors to diagnose.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("figures"),
        help="Directory for SVG plots and summary CSV.",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    priors = read_selection_priors(args.param_file)

    summary_rows = []
    for q in args.q:
        q_label = f"{q:g}".replace(".", "p")
        plot_q(priors, q, args.outdir / f"scaled_selection_prior_Q{q_label}.svg")
        for param, values in priors.items():
            summary_rows.append(summarize(param, q, values))

    write_summary(args.outdir / "scaled_selection_prior_summary.csv", summary_rows)


if __name__ == "__main__":
    main()
