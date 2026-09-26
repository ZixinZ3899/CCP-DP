from pathlib import Path
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator


ROOT = Path(__file__).resolve().parent
DATA_PATH = ROOT / "results" / "coverage_robustness_data.csv"
OUTPUT_DIR = Path(os.environ.get("COVERAGE_FIG_OUTPUT_DIR", ROOT))
OUT_PREFIX = OUTPUT_DIR / "coverage_robustness_three_panel"

# Palette aligned with the companion mechanism figure, with vivid red/green accents.
BLUE = "#4D95B8"
RED = "#E31A1C"
GREEN = "#00A651"
GOLD = "#ECAF4F"
INK = "#2E2E2E"
GRID = "#E8E8E8"

METHODS = ["CCP-DP", "BBS", "CPL", "ITR"]
COLORS = {"CCP-DP": BLUE, "BBS": RED, "CPL": GREEN, "ITR": GOLD}
MARKERS = {"CCP-DP": "o", "BBS": "s", "CPL": "^", "ITR": "D"}
DATASETS = ["Chandak-108", "Chandak-150", "Srinivas-110"]
TITLES = {
    "Chandak-108": "Chandak et al. (108 bp)",
    "Chandak-150": "Chandak et al. (150 bp)",
    "Srinivas-110": "Srinivasavaradhan et al. (110 bp)",
}


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "sans-serif"],
        "font.size": 7,
        "axes.titlesize": 7.6,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3,
        "legend.fontsize": 6.1,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.8,
        "ytick.major.size": 2.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
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


def validate_data(df):
    expected_rows = len(DATASETS) * len(METHODS) * 4
    if len(df) != expected_rows:
        raise ValueError(f"Expected {expected_rows} rows, found {len(df)}")
    if set(df.dataset) != set(DATASETS) or set(df.method) != set(METHODS):
        raise ValueError("Unexpected dataset or method labels")
    if set(df.coverage) != {5, 10, 15, 20}:
        raise ValueError("Coverage must contain 5, 10, 15 and 20")
    if not df.exact_reconstruction_percent.between(0, 100).all():
        raise ValueError("Exact reconstruction values must lie in [0, 100]")
    if (df.exact_sd < 0).any():
        raise ValueError("SD must be non-negative")
    if (df.mean_edit_distance <= 0).any():
        raise ValueError("Mean edit distance must be strictly positive before log scaling")


def style_left_axis(ax, show_ylabel):
    ax.set_xlim(4.5, 20.5)
    ax.set_xticks([5, 10, 15, 20])
    ax.set_xlabel("Read coverage (×)")
    ax.set_ylim(0, 102)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Exact reconstruction (%)" if show_ylabel else "")
    if not show_ylabel:
        ax.tick_params(axis="y", labelleft=False)
        ax.spines["left"].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.45, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(direction="out", pad=2)


def style_right_axis(ax, show_ylabel):
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 10)
    ticks = [1e-4, 1e-3, 1e-2, 1e-1, 1, 10]
    labels = ["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "1", "10"]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FixedFormatter(labels))
    ax.minorticks_off()
    ax.spines["top"].set_visible(False)
    if show_ylabel:
        ax.spines["right"].set_visible(True)
        ax.spines["right"].set_color(INK)
        ax.set_ylabel("Mean edit distance (log scale)")
        ax.tick_params(axis="y", right=True, labelright=True, colors=INK, direction="out", pad=2)
    else:
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", right=False, labelright=False)


def plot_panel(ax, frame, dataset, show_left, show_right):
    right = ax.twinx()
    for method in METHODS:
        d = frame[(frame.dataset == dataset) & (frame.method == method)].sort_values("coverage")
        color = COLORS[method]
        marker = MARKERS[method]
        emphasis = 1.45 if method == "CCP-DP" else 1.05
        ax.errorbar(
            d.coverage,
            d.exact_reconstruction_percent,
            yerr=d.exact_sd,
            color=color,
            marker=marker,
            ms=3.2,
            mfc=color,
            mec=color,
            mew=0.55,
            lw=emphasis,
            elinewidth=0.65,
            capsize=1.7,
            capthick=0.65,
            zorder=4,
        )
        right.plot(
            d.coverage,
            d.mean_edit_distance,
            color=color,
            marker=marker,
            ms=3.1,
            mfc="white",
            mec=color,
            mew=0.75,
            lw=emphasis,
            linestyle=(0, (3.0, 1.8)),
            alpha=0.9,
            zorder=3,
        )
    style_left_axis(ax, show_left)
    style_right_axis(right, show_right)
    ax.set_title(TITLES[dataset], loc="left", pad=4)


def add_panel_label(fig, ax, label):
    pos = ax.get_position()
    fig.text(pos.x0 - 0.018, pos.y1 + 0.014, label, fontsize=9, fontweight="bold", ha="right", va="bottom")


def main():
    df = pd.read_csv(DATA_PATH)
    validate_data(df)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(7.25, 2.55), sharey=True)
    for index, (ax, dataset) in enumerate(zip(axes, DATASETS)):
        plot_panel(ax, df, dataset, show_left=index == 0, show_right=index == 2)

    method_handles = [
        Line2D(
            [0],
            [0],
            color=COLORS[method],
            marker=MARKERS[method],
            lw=1.35 if method == "CCP-DP" else 1.05,
            ms=3.5,
            label=method,
        )
        for method in METHODS
    ]
    metric_handles = [
        Line2D([0], [0], color=INK, marker="o", mfc=INK, lw=1.15, ms=3.2, label="Exact reconstruction"),
        Line2D(
            [0],
            [0],
            color=INK,
            marker="o",
            mfc="white",
            lw=1.15,
            ms=3.2,
            linestyle=(0, (3.0, 1.8)),
            label="Mean edit distance",
        ),
    ]
    fig.legend(
        handles=method_handles + metric_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=6,
        handlelength=1.8,
        handletextpad=0.45,
        columnspacing=1.0,
    )

    fig.subplots_adjust(left=0.08, right=0.92, top=0.77, bottom=0.20, wspace=0.20)
    for ax, label in zip(axes, "abc"):
        add_panel_label(fig, ax, label)

    fig.savefig(OUT_PREFIX.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(OUT_PREFIX.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT_PREFIX.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(OUT_PREFIX.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


if __name__ == "__main__":
    main()
