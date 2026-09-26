#!/usr/bin/env python3
"""Draw adaptive evidence-allocation analyses from a compact dataset summary."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


DATASET_ORDER = ["Chandak-108", "Chandak-150", "Srinivas-110", "Erlich-152"]
COLORS = {
    # Preserve the original point-figure colours.
    "Chandak-108": "#0B60A7",
    "fasta_2-110": "#E67E22",
    "Chandak-150": "#009E73",
    "Srinivas-110": "#7A3E9D",
    "Erlich-152": "#177CC5",
}
INK = "#2E2E2E"
MUTED = "#6F7780"
GRID = "#E6E8EB"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Nimbus Sans",
            "font.size": 10.5,
            "axes.titlesize": 14,
            "axes.titleweight": "bold",
            "axes.labelsize": 11.5,
            "xtick.labelsize": 10.5,
            "ytick.labelsize": 10.5,
            "axes.linewidth": 0.9,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def validate_data(df: pd.DataFrame) -> None:
    expected = set(DATASET_ORDER)
    if set(df.dataset) != expected or len(df) != len(expected):
        raise ValueError("Expected exactly one row for each of the five datasets")
    for col in ["guide_failure_pct", "structured_pct", "delta_success_pp", "time_reduction_pct"]:
        if df[col].isna().any():
            raise ValueError(f"Missing values in {col}")


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color=GRID, linewidth=0.75, zorder=0)
    ax.set_axisbelow(True)


def add_panel_letter(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.14, 1.08, letter, transform=ax.transAxes, fontsize=16, fontweight="bold", ha="left")


def scatter_points(ax: plt.Axes, df: pd.DataFrame, x: str, y: str) -> None:
    for _, row in df.iterrows():
        ax.scatter(
            row[x], row[y], s=125, color=COLORS[row.dataset], edgecolor="white",
            linewidth=1.1, zorder=4
        )


def add_fit(ax: plt.Axes, df: pd.DataFrame, x: str, y: str, xlim: tuple[float, float]) -> None:
    slope, intercept = np.polyfit(df[x].to_numpy(float), df[y].to_numpy(float), 1)
    xx = np.linspace(xlim[0], xlim[1], 100)
    ax.plot(xx, slope * xx + intercept, color="#AEB5BC", lw=1.25, linestyle=(0, (4, 4)), zorder=1)


def draw_difficulty(ax: plt.Axes, df: pd.DataFrame, panel_letter: str = "a") -> None:
    xlim = (-0.6, 14.7)
    scatter_points(ax, df, "guide_failure_pct", "structured_pct")
    add_fit(ax, df, "guide_failure_pct", "structured_pct", xlim)
    offsets = {
        "Erlich-152": (7, 8), "Chandak-150": (7, 7), "Srinivas-110": (7, -12),
        "fasta_2-110": (-82, 12), "Chandak-108": (-82, -19),
    }
    for _, row in df.iterrows():
        dx, dy = offsets[row.dataset]
        ax.annotate(row.dataset, (row.guide_failure_pct, row.structured_pct), xytext=(dx, dy),
                    textcoords="offset points", fontsize=9.3)
    ax.set_xlim(*xlim)
    ax.set_ylim(-3, 106)
    ax.set_xlabel("Guide failure rate (%)")
    ax.set_ylabel("Structured refinement (%)")
    ax.set_title("Difficulty vs activation", loc="left", pad=9)
    ax.text(0.10, 0.62, "Harder datasets\ntend to activate\nmore refinement", transform=ax.transAxes,
            color=MUTED, fontsize=9.2, fontstyle="italic")
    clean_axes(ax)
    add_panel_letter(ax, panel_letter)


def draw_gain(ax: plt.Axes, df: pd.DataFrame, panel_letter: str = "b") -> None:
    xlim = (-5, 105)
    scatter_points(ax, df, "structured_pct", "delta_success_pp")
    add_fit(ax, df, "structured_pct", "delta_success_pp", xlim)
    offsets = {
        "Erlich-152": (7, 8), "Chandak-150": (7, 7), "Srinivas-110": (7, 8),
        "fasta_2-110": (-56, 9), "Chandak-108": (-62, 10),
    }
    for _, row in df.iterrows():
        dx, dy = offsets[row.dataset]
        ax.annotate(row.dataset, (row.structured_pct, row.delta_success_pp), xytext=(dx, dy),
                    textcoords="offset points", fontsize=9.3)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.3, 10.1)
    ax.set_xlabel("Structured refinement (%)")
    ax.set_ylabel("Success gain over guide (pp)")
    ax.set_title("Activation vs accuracy gain", loc="left", pad=9)
    ax.text(0.34, 0.64, "Greater refinement\nis associated with\nlarger gain", transform=ax.transAxes,
            color=MUTED, fontsize=9.2, fontstyle="italic")
    clean_axes(ax)
    add_panel_letter(ax, panel_letter)


def draw_efficiency(ax: plt.Axes, df: pd.DataFrame, panel_letter: str = "c") -> None:
    rows = df.set_index("dataset").loc[DATASET_ORDER].reset_index()
    ys = np.arange(len(rows))
    for y, row in zip(ys, rows.itertuples()):
        value = float(row.time_reduction_pct)
        display_value = int(np.rint(value))
        if display_value == 0 and value != 0:
            display_value = 1 if value > 0 else -1
        ax.hlines(y, 0, value, color=COLORS[row.dataset], lw=2.5, zorder=2)
        ax.scatter(value, y, s=110, color=COLORS[row.dataset], edgecolor="white", linewidth=1.0, zorder=4)
        ax.text(value + (3.2 if value >= 0 else -3.2), y, f"{display_value}",
                ha="left" if value >= 0 else "right", va="center", fontsize=9.5,
                fontweight="semibold", color=COLORS[row.dataset])
    ax.axvline(0, color=INK, lw=1.25, zorder=3)
    ax.set_yticks(ys, rows.dataset)
    ax.invert_yaxis()
    ax.set_xlim(-10, 82)
    ax.set_xlabel("Time reduction vs decode-all (%)")
    # Shift the title slightly right for better alignment within the panel.
    ax.set_title("Efficiency benefit", x=0.06, ha="left", pad=9)
    ax.grid(axis="x", color=GRID, linewidth=0.75, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=6)
    add_panel_letter(ax, panel_letter)


def replace_svg_font(path: Path) -> None:
    svg = path.read_text(encoding="utf-8")
    path.write_text(svg.replace("Nimbus Sans", "Arial"), encoding="utf-8")


def main() -> None:
    args = parse_args()
    setup_style()
    df = pd.read_csv(args.data)
    validate_data(df)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_dir / "adaptive_evidence_allocation_diverging"

    fig, axes = plt.subplots(1, 3, figsize=(16, 7.8), gridspec_kw={"width_ratios": [1.05, 1.05, 0.90]})
    draw_difficulty(axes[0], df)
    draw_gain(axes[1], df)
    draw_efficiency(axes[2], df)
    fig.suptitle("Adaptive evidence allocation across representative datasets", fontsize=17,
                 fontweight="bold", y=0.975)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.80, bottom=0.15, wspace=0.28)
    fig.savefig(prefix.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight")
    replace_svg_font(prefix.with_suffix(".svg"))
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()