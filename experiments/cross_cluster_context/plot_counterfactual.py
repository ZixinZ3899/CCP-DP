#!/usr/bin/env python3
"""Draw the three-panel counterfactual validation figure from the supplied CSVs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


PALETTE = {
    # Diverging blue-red palette sampled from the supplied reference figure.
    "Local only": "#478DB5",
    "Table-permuted": "#B9D4DF",
    "Wrong-dataset": "#A33B3E",
    "Correct Cross-LOO": "#E49B80",
    "Newly rescued": "#478DB5",
    "Newly harmed": "#A33B3E",
    "Net gain": "#3065A1",
    "Accent": "#3065A1",
    "Grid": "#E6E8EB",
    "Axis": "#2E2E2E",
    "Muted": "#666666",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--pairwise", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Nimbus Sans",
            "font.size": 10.5,
            "axes.titlesize": 15,
            "axes.titleweight": "bold",
            "axes.labelsize": 11.5,
            "axes.linewidth": 0.9,
            "axes.edgecolor": PALETTE["Axis"],
            "xtick.labelsize": 10.5,
            "ytick.labelsize": 10.5,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.fontsize": 9.5,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def clean_axes(ax: plt.Axes, grid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(PALETTE["Axis"])
    ax.spines["bottom"].set_color(PALETTE["Axis"])
    if grid:
        ax.grid(axis="y", color=PALETTE["Grid"], linewidth=0.75, alpha=0.9, zorder=0)
    ax.tick_params(length=4, color=PALETTE["Axis"])


def add_group_separator(ax: plt.Axes) -> None:
    ax.axvline(0.5, color="#C9D3DE", linestyle=(0, (3, 3)), linewidth=0.9, zorder=0)


def add_panel_letter(ax: plt.Axes, letter: str) -> None:
    ax.text(
        -0.145,
        1.13,
        letter,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=19,
        fontweight="bold",
        color="#111111",
        clip_on=False,
    )


def add_shadow_bars(ax: plt.Axes, xs: np.ndarray, heights: np.ndarray, width: float, colors: list[str]):
    # A small vector shadow adds depth without raster gradients or 3-D distortion.
    shadow_dx = 0.014
    ax.bar(
        xs + shadow_dx,
        heights,
        width=width,
        color="#CDD7E2",
        alpha=0.42,
        linewidth=0,
        zorder=1,
    )
    bars = ax.bar(
        xs,
        heights,
        width=width,
        color=colors,
        edgecolor="white",
        linewidth=0.7,
        zorder=2,
    )
    return bars


def add_compare_bracket(
    ax: plt.Axes,
    x1: float,
    x2: float,
    y: float,
    y1: float,
    y2: float,
    label: str,
    log_axis: bool = False,
) -> None:
    color = PALETTE["Accent"]
    lw = 1.35
    ax.plot([x1, x1, x2], [y1, y, y], color=color, lw=lw, clip_on=False, zorder=5)
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x2, y),
        arrowprops=dict(arrowstyle="->", color=color, lw=lw, shrinkA=0, shrinkB=0),
        annotation_clip=False,
        zorder=5,
    )
    text_y = y * 1.13 if log_axis else y + 0.7
    ax.text(
        (x1 + x2) / 2,
        text_y,
        label,
        color=color,
        ha="center",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
        clip_on=False,
    )


def value_and_sd(row: pd.Series, value_col: str, sd_col: str, digits: int) -> str:
    value = float(row[value_col])
    sd = float(row[sd_col])
    base = f"{value:.{digits}f}"
    if sd > 0:
        return f"{base}\n± {sd:.{digits}f}"
    return base


def draw_success_panel(
    ax: plt.Axes,
    summary: pd.DataFrame,
    datasets: list[str],
    variants: list[str],
    panel_letter: str = "a",
) -> None:
    centers = np.arange(len(datasets), dtype=float)
    width = 0.19
    offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * width

    rows_by_key = summary.set_index(["dataset", "variant"])
    for j, variant in enumerate(variants):
        rows = [rows_by_key.loc[(dataset, variant)] for dataset in datasets]
        vals = np.array([r["success_mean_percent"] for r in rows], dtype=float)
        sds = np.array([r["success_sd_pp"] for r in rows], dtype=float)
        xs = centers + offsets[j]
        bars = add_shadow_bars(ax, xs, vals, width * 0.92, [PALETTE[variant]] * len(xs))
        for x, val in zip(xs, vals):
            if val == 0:
                ax.plot(x, 0.45, marker="_", markersize=14, markeredgewidth=2.1,
                        color=PALETTE[variant], zorder=5, clip_on=False)
        ax.errorbar(
            xs,
            vals,
            yerr=sds,
            fmt="none",
            ecolor="#171717",
            elinewidth=1.45,
            capsize=5.0,
            capthick=1.45,
            zorder=4,
        )
        for bar, row, val in zip(bars, rows, vals):
            text = value_and_sd(row, "success_mean_percent", "success_sd_pp", 2)
            y = max(val + 1.4, 1.5) if val > 0 else 1.3
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                y,
                text,
                ha="center",
                va="bottom",
                fontsize=9.3 if row["success_sd_pp"] > 0 else 8.8,
                fontweight="semibold" if row["success_sd_pp"] == 0 else "normal",
                linespacing=1.05,
                zorder=6,
            )

    local_idx = variants.index("Local only")
    correct_idx = variants.index("Correct Cross-LOO")
    for i, dataset in enumerate(datasets):
        local = float(rows_by_key.loc[(dataset, "Local only"), "success_mean_percent"])
        correct = float(rows_by_key.loc[(dataset, "Correct Cross-LOO"), "success_mean_percent"])
        # Keep the comparison bracket clearly above the bar-value labels.
        y = 104.2 if i == 0 else 103.0
        add_compare_bracket(
            ax,
            centers[i] + offsets[local_idx],
            centers[i] + offsets[correct_idx],
            y,
            local + 7.0,
            correct + 0.7,
            f"+{correct - local:.2f} pp",
        )

    ax.set_title("Exact reconstruction rate", pad=25)
    ax.set_ylabel("Success (%)")
    ax.set_xticks(centers, datasets)
    # Extra headroom is reserved for the raised brackets and their labels.
    ax.set_ylim(0, 108.5)
    ax.set_yticks(np.arange(0, 101, 20))
    add_group_separator(ax)
    clean_axes(ax)
    add_panel_letter(ax, panel_letter)


def draw_mean_ed_panel(
    ax: plt.Axes,
    summary: pd.DataFrame,
    datasets: list[str],
    variants: list[str],
    panel_letter: str = "b",
) -> None:
    centers = np.arange(len(datasets), dtype=float)
    width = 0.19
    offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * width
    rows_by_key = summary.set_index(["dataset", "variant"])

    for j, variant in enumerate(variants):
        rows = [rows_by_key.loc[(dataset, variant)] for dataset in datasets]
        vals = np.array([r["mean_ed_mean"] for r in rows], dtype=float)
        sds = np.array([r["mean_ed_sd"] for r in rows], dtype=float)
        xs = centers + offsets[j]
        bars = add_shadow_bars(ax, xs, vals, width * 0.92, [PALETTE[variant]] * len(xs))
        ax.errorbar(
            xs,
            vals,
            yerr=sds,
            fmt="none",
            ecolor="#171717",
            elinewidth=1.45,
            capsize=5.0,
            capthick=1.45,
            zorder=4,
        )
        for bar, row, val in zip(bars, rows, vals):
            text = value_and_sd(row, "mean_ed_mean", "mean_ed_sd", 3)
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                val * (1.24 if val < 1 else 1.16),
                text,
                ha="center",
                va="bottom",
                fontsize=9.3 if row["mean_ed_sd"] > 0 else 8.8,
                linespacing=1.05,
                zorder=6,
            )

    local_idx = variants.index("Local only")
    correct_idx = variants.index("Correct Cross-LOO")
    # Dynamically place each bracket well above the tallest bar and its
    # numerical annotation. Longer end caps retain a clear bracket shape.
    for i, dataset in enumerate(datasets):
        local = float(
            rows_by_key.loc[
                (dataset, "Local only"),
                "mean_ed_mean",
            ]
        )
        correct = float(
            rows_by_key.loc[
                (dataset, "Correct Cross-LOO"),
                "mean_ed_mean",
            ]
        )
        reduction = (correct / local - 1.0) * 100.0

        group_tops = []
        for variant in variants:
            row = rows_by_key.loc[(dataset, variant)]
            mean_ed = float(row["mean_ed_mean"])
            mean_sd = float(row["mean_ed_sd"])
            if not np.isfinite(mean_sd):
                mean_sd = 0.0
            group_tops.append(mean_ed + mean_sd)

        # On the logarithmic axis, multiplicative spacing is appropriate.
        # 2.20 leaves room for the two-line value/SD annotation.
        y = max(group_tops) * 2.20

        # Extend both bracket arms down to just above the corresponding
        # bar tops. Place them at the outer bar edges so that they do not
        # cross the numerical annotations.
        local_sd = float(
            rows_by_key.loc[
                (dataset, "Local only"),
                "mean_ed_sd",
            ]
        )
        correct_sd = float(
            rows_by_key.loc[
                (dataset, "Correct Cross-LOO"),
                "mean_ed_sd",
            ]
        )
        if not np.isfinite(local_sd):
            local_sd = 0.0
        if not np.isfinite(correct_sd):
            correct_sd = 0.0

        left_x = (
            centers[i]
            + offsets[local_idx]
            - width * 0.46
        )
        right_x = (
            centers[i]
            + offsets[correct_idx]
            + width * 0.46
        )

        left_bottom = (local + local_sd) * 1.08
        right_bottom = (correct + correct_sd) * 1.08

        add_compare_bracket(
            ax,
            left_x,
            right_x,
            y,
            left_bottom,
            right_bottom,
            f"{reduction:.1f}%",
            log_axis=True,
        )

    ax.set_title("Mean edit distance", pad=25)
    ax.set_ylabel("Mean ED")
    ax.set_xticks(centers, datasets)
    ax.set_yscale("log")
    ax.set_ylim(1e-2, 14)
    ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=5))
    ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=12))
    ax.grid(axis="y", which="major", color=PALETTE["Grid"], linewidth=0.75, alpha=0.9, zorder=0)
    add_group_separator(ax)
    clean_axes(ax, grid=False)
    add_panel_letter(ax, panel_letter)


def draw_rescue_panel(
    ax: plt.Axes,
    pairwise: pd.DataFrame,
    datasets: list[str],
    panel_letter: str = "c",
) -> None:
    centers = np.arange(len(datasets), dtype=float)
    height = 0.22
    offsets = np.array([-1.15, 0.0, 1.15]) * height
    measures = [
        ("correct_only_exact_rescued", "Newly rescued", 1),
        ("local_only_exact_correct_harmed", "Newly harmed", -1),
        ("net_exact_gain", "Net gain", 1),
    ]
    rows_by_dataset = pairwise.set_index("dataset")

    for j, (column, label, sign) in enumerate(measures):
        vals = np.array([float(rows_by_dataset.loc[d, column]) * sign for d in datasets])
        ys = centers + offsets[j]
        ax.barh(
            ys + 0.018,
            vals,
            height=height * 0.92,
            color="#CBD5DE",
            alpha=0.40,
            linewidth=0,
            zorder=1,
        )
        bars = ax.barh(
            ys,
            vals,
            height=height * 0.92,
            color=PALETTE[label],
            edgecolor="#FFFFFF",
            linewidth=0.75,
            zorder=2,
        )
        for bar, val in zip(bars, vals):
            if label == "Newly harmed":
                text = f"−{abs(int(val))}"
                x = val - (17 if abs(val) >= 20 else 14)
                ha = "right"
            elif label == "Net gain":
                text = f"+{int(val)}"
                x = val + 12
                ha = "left"
            else:
                text = f"{int(val)}"
                x = val + 12
                ha = "left"
            ax.text(
                x,
                bar.get_y() + bar.get_height() / 2,
                text,
                ha=ha,
                va="center",
                fontsize=9.3,
                fontweight="bold" if label == "Net gain" else "semibold",
                color=PALETTE["Accent"] if label == "Net gain" else "#171717",
                zorder=6,
            )

    ax.axvline(0, color=PALETTE["Axis"], linewidth=1.15, zorder=3)
    ax.set_title("Cluster-level exact changes", pad=25)
    ax.text(
        0.5,
        0.985,
        "Correct Cross-LOO vs Local only",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=10.5,
        color="#303640",
    )
    ax.set_xlabel("Clusters  (harmed − / rescued +)")
    ax.set_yticks(centers, datasets)
    ax.set_xlim(-75, 470)
    ax.set_xticks([-50, 0, 100, 200, 300, 400])
    ax.set_ylim(-0.55, 1.55)
    ax.invert_yaxis()
    ax.grid(axis="x", color=PALETTE["Grid"], linewidth=0.75, alpha=0.9, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(PALETTE["Axis"])
    ax.tick_params(axis="y", length=0, pad=7)
    ax.tick_params(axis="x", length=4, color=PALETTE["Axis"])
    add_panel_letter(ax, panel_letter)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(args.summary)
    pairwise = pd.read_csv(args.pairwise)

    datasets = ["Srinivas-110", "Chandak-108"]
    variants = ["Local only", "Table-permuted", "Wrong-dataset", "Correct Cross-LOO"]

    missing_summary = {
        (d, v) for d in datasets for v in variants
    } - set(zip(summary["dataset"], summary["variant"]))
    missing_pairwise = set(datasets) - set(pairwise["dataset"])
    if missing_summary or missing_pairwise:
        raise ValueError(f"Missing required rows: summary={missing_summary}, pairwise={missing_pairwise}")

    setup_style()
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(16, 7.35),
        gridspec_kw={"width_ratios": [1.02, 1.02, 1.16], "wspace": 0.24},
    )
    fig.subplots_adjust(left=0.055, right=0.985, top=0.83, bottom=0.19)

    draw_success_panel(axes[0], summary, datasets, variants)
    draw_mean_ed_panel(axes[1], summary, datasets, variants)
    draw_rescue_panel(axes[2], pairwise, datasets)

    legend_order = ["Local only", "Table-permuted", "Wrong-dataset", "Correct Cross-LOO"]
    legend_ab = [Patch(facecolor=PALETTE[v], edgecolor="none", label=v) for v in legend_order]
    fig.legend(
        handles=legend_ab,
        loc="lower center",
        bbox_to_anchor=(0.345, 0.047),
        ncol=4,
        frameon=False,
        columnspacing=1.55,
        handlelength=2.0,
        handleheight=1.05,
    )

    legend_c = [
        Patch(facecolor=PALETTE["Newly rescued"], edgecolor="none", label="Newly rescued"),
        Patch(facecolor=PALETTE["Newly harmed"], edgecolor="none", label="Newly harmed"),
        Patch(facecolor=PALETTE["Net gain"], edgecolor="none", label="Net gain"),
    ]
    fig.legend(
        handles=legend_c,
        loc="lower center",
        bbox_to_anchor=(0.825, 0.047),
        ncol=3,
        frameon=False,
        columnspacing=1.0,
        handlelength=2.0,
        handleheight=1.1,
    )

    stem = args.output_dir / "counterfactual_validation_macaron"
    fig.savefig(stem.with_suffix(".png"), dpi=400, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.08)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.08)
    # PowerPoint handles Arial text predictably. Keep SVG text editable and
    # advertise Arial while using its metric-compatible Nimbus Sans for rendering.
    svg_path = stem.with_suffix(".svg")
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(svg_text.replace("Nimbus Sans", "Arial"), encoding="utf-8")
    plt.close(fig)


if __name__ == "__main__":
    main()