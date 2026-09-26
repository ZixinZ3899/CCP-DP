#!/usr/bin/env python3
"""Plot the v8.2 component ablation as supplementary heat maps."""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ORDER = [
    "guide_only", "without_cross_loo", "without_blind_ids",
    "decode_all", "local_trellis_only", "full",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    with Path(args.summary).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    datasets = list(dict.fromkeys(row["dataset"] for row in rows))
    labels = {
        row["variant_id"]: row["variant_label"] for row in rows
    }
    lookup = {(row["dataset"], row["variant_id"]): row for row in rows}
    success = np.array([
        [float(lookup[(dataset, variant)]["success_percent"])
         for dataset in datasets]
        for variant in ORDER
    ])
    runtime_ratio = np.array([
        [float(lookup[(dataset, variant)]["runtime_ratio_vs_full"])
         for dataset in datasets]
        for variant in ORDER
    ])

    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.1), constrained_layout=True)
    panels = [
        (success, "Exact reconstruction (%)", "YlGnBu", 0, 100, ".2f"),
        (runtime_ratio, "Runtime relative to Full", "YlOrBr", None, None, ".2f"),
    ]
    for panel, (ax, (matrix, cbar_label, cmap, vmin, vmax, fmt)) in enumerate(
        zip(axes, panels)
    ):
        image = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(datasets)), datasets, rotation=25, ha="right")
        ax.set_yticks(range(len(ORDER)), [labels[item] for item in ORDER])
        ax.set_title(chr(ord("a") + panel), loc="left", fontweight="bold")
        for row in range(matrix.shape[0]):
            for col in range(matrix.shape[1]):
                value = matrix[row, col]
                text = format(value, fmt)
                if panel == 1:
                    text += "×"
                ax.text(col, row, text, ha="center", va="center", fontsize=8)
        cbar = fig.colorbar(image, ax=ax, shrink=0.83)
        cbar.set_label(cbar_label)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(outdir / f"ablation_summary.{suffix}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
