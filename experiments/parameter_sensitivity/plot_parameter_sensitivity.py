#!/usr/bin/env python3
"""Plot success and routing sensitivity from summary_metrics.csv."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PARAMETERS = [
    ("coverage_threshold", "Coverage switch", "Reads"),
    ("reuse_threshold", "Reuse threshold", "Threshold"),
    ("phase_thresholds", "Phase-load thresholds", "Low / high"),
    ("cross_weight_scale", "Cross-weight scale", "Scale"),
    ("route_min_fraction", "Minimum route fraction", "Fraction"),
    ("phase_band", "Phase band", "Band width"),
]
DEFAULT = {
    "coverage_threshold": (40.0, "40"),
    "reuse_threshold": (0.50, "0.50"),
    "phase_thresholds": (1.0, "3/6"),
    "cross_weight_scale": (1.0, "1.00"),
    "route_min_fraction": (0.10, "0.10"),
    "phase_band": (6.0, "auto"),
}
COLORS = {
    "Srinivas-110": "#0072B2",
    "Chandak-108": "#E69F00",
    "Erlich-152": "#009E73",
}
MARKERS = {"Srinivas-110": "o", "Chandak-108": "s", "Erlich-152": "D"}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def build_series(rows: list[dict[str, str]], dataset: str, parameter: str,
                 metric: str) -> tuple[list[float], list[float], list[str]]:
    baseline = next(row for row in rows
                    if row["dataset"] == dataset and row["parameter"] == "baseline")
    baseline_value = float(baseline[metric])
    x0, label0 = DEFAULT[parameter]
    points = [(x0, baseline_value, label0)]
    for row in rows:
        if row["dataset"] == dataset and row["parameter"] == parameter:
            points.append((float(row["plot_x"]), float(row[metric]),
                           row["parameter_value"]))
    points.sort(key=lambda item: item[0])
    return ([p[0] for p in points], [p[1] for p in points], [p[2] for p in points])


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D9DEE3", linewidth=0.65, alpha=0.75)
    ax.tick_params(labelsize=8.5, width=0.8, length=3)


def plot(rows: list[dict[str, str]], output: Path, metric: str,
         ylabel: str, delta: bool) -> None:
    datasets = [name for name in COLORS
                if any(row["dataset"] == name for row in rows)]
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.4))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.10, top=0.86,
                        wspace=0.27, hspace=0.46)
    for ax, (parameter, title, xlabel) in zip(axes.flat, PARAMETERS):
        for dataset in datasets:
            xs, ys, labels = build_series(rows, dataset, parameter, metric)
            if delta:
                default_y = ys[xs.index(DEFAULT[parameter][0])]
                ys = [value - default_y for value in ys]
            ax.plot(xs, ys, color=COLORS[dataset], marker=MARKERS[dataset],
                    markersize=5.0, linewidth=1.7, label=dataset)
            if parameter == "phase_thresholds":
                ax.set_xticks(xs, labels)
        ax.axvline(DEFAULT[parameter][0], color="#6C757D", linewidth=0.9,
                   linestyle="--", zorder=0)
        if delta:
            ax.axhline(0.0, color="#343A40", linewidth=0.8, zorder=0)
        ax.set_title(title, fontsize=10.0, fontweight="semibold", pad=7)
        ax.set_xlabel(xlabel, fontsize=9.2)
        style_axis(ax)
    axes[0, 0].set_ylabel(ylabel, fontsize=9.5)
    axes[1, 0].set_ylabel(ylabel, fontsize=9.5)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=max(1, len(labels)),
               frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 0.975))
    for suffix in ("svg", "pdf", "png"):
        fig.savefig(output.with_suffix(f".{suffix}"), dpi=400,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    rows = read_rows(args.summary)
    output_dir = args.output_dir or args.summary.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    plot(rows, output_dir / "parameter_sensitivity_success",
         "success_percent_mean", "Δ exact reconstruction (pp)", True)
    plot(rows, output_dir / "parameter_sensitivity_routing",
         "decoded_fraction_percent_mean", "Structured decoding (%)", False)
    print(f"Figures written to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
