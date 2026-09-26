#!/usr/bin/env python3
"""Combine the two mechanism-validation figures into a coherent 2 x 3 layout."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=HERE / "results" / "counterfactual_summary.csv")
    parser.add_argument("--pairwise", type=Path, default=HERE / "results" / "local_vs_correct_pairwise.csv")
    parser.add_argument("--allocation", type=Path, default=HERE / "results" / "adaptive_evidence_allocation_data.csv")
    parser.add_argument("--output-dir", type=Path, default=HERE / "figures")
    return parser.parse_args()


def replace_svg_font(path: Path) -> None:
    svg = path.read_text(encoding="utf-8")
    path.write_text(svg.replace("Nimbus Sans", "Arial"), encoding="utf-8")


def main() -> None:
    args = parse_args()
    counter = load_module("counterfactual_original", HERE / "plot_counterfactual.py")
    allocation = load_module("allocation_original", HERE / "plot_adaptive.py")
    counter.setup_style()

    summary = pd.read_csv(args.summary)
    pairwise = pd.read_csv(args.pairwise)
    allocation_df = pd.read_csv(args.allocation)
    allocation.validate_data(allocation_df)

    datasets = ["Srinivas-110", "Chandak-108"]
    variants = ["Local only", "Table-permuted", "Wrong-dataset", "Correct Cross-LOO"]

    fig = plt.figure(figsize=(16, 13.0))
    grid = fig.add_gridspec(
        2, 3, height_ratios=[0.78, 1.22], left=0.055, right=0.975,
        top=0.965, bottom=0.065, wspace=0.26, hspace=0.31
    )
    top_axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
    bottom_axes = [fig.add_subplot(grid[1, i]) for i in range(3)]

    counter.draw_success_panel(top_axes[0], summary, datasets, variants, panel_letter="a")
    counter.draw_mean_ed_panel(top_axes[1], summary, datasets, variants, panel_letter="b")
    counter.draw_rescue_panel(top_axes[2], pairwise, datasets, panel_letter="c")

    allocation.draw_difficulty(bottom_axes[0], allocation_df, panel_letter="d")
    allocation.draw_gain(bottom_axes[1], allocation_df, panel_letter="e")
    allocation.draw_efficiency(bottom_axes[2], allocation_df, panel_letter="f")
    for ax in bottom_axes:
        ax.set_box_aspect(1.12)

    variant_handles = [Patch(facecolor=counter.PALETTE[v], edgecolor="none", label=v) for v in variants]
    rescue_handles = [
        Patch(facecolor=counter.PALETTE["Newly rescued"], edgecolor="none", label="Newly rescued"),
        Patch(facecolor=counter.PALETTE["Newly harmed"], edgecolor="none", label="Newly harmed"),
        Patch(facecolor=counter.PALETTE["Net gain"], edgecolor="none", label="Net gain"),
    ]
    fig.legend(handles=variant_handles, loc="center", bbox_to_anchor=(0.35, 0.605), ncol=4,
               frameon=False, handlelength=1.8, columnspacing=1.2)
    fig.legend(handles=rescue_handles, loc="center", bbox_to_anchor=(0.825, 0.605), ncol=3,
               frameon=False, handlelength=1.8, columnspacing=1.15)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_dir / "combined_counterfactual_adaptive"
    fig.savefig(prefix.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight")
    replace_svg_font(prefix.with_suffix(".svg"))
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
