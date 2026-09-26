from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parent
DATA_PATH = ROOT / "results" / "ccpdp_analysis_data.csv"
OUT_PREFIX = ROOT / "figures" / "ccpdp_mechanism_analysis"

# Palette adapted from the supplied reference: cool blue, peach/coral,
# teal, warm gold and restrained red.
BLUE = "#4D95B8"
BLUE_LIGHT = "#9BC5D3"
CORAL = "#E8916F"
PEACH = "#F2B092"
TEAL = "#5AA69F"
GOLD = "#ECAF4F"
GREEN = "#6AAE72"
RED = "#C83D35"
REDUCTION_FILL = "#BFDCCF"
INK = "#2E2E2E"
LIGHT = "#D9E7EC"
GRID = "#E8E8E8"

DATASETS = ["Srinivas", "Chandak-108"]
DATASET_LABELS = {"Srinivas": "Srinivas (110 bp)", "Chandak-108": "Chandak (108 bp)"}
DATASET_COLORS = {"Srinivas": BLUE, "Chandak-108": CORAL}
ERROR_ORDER = ["Substitution-only", "Paired indel", "Mixed"]
ERROR_LABELS = ["Sub.", "Paired\nindel", "Mixed"]
ERROR_COLORS = [BLUE, CORAL, TEAL]


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.titlesize": 7.5,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3,
        "legend.fontsize": 6.1,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": INK,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def panel_label(fig, ax, label, dx=0.026, dy=0.010):
    """Place panel letters in the outer gutter, clear of titles and axes."""
    pos = ax.get_position()
    fig.text(
        pos.x0 - dx,
        pos.y1 + dy,
        label,
        fontsize=9,
        fontweight="bold",
        ha="right",
        va="bottom",
    )


def clean_axis(ax, ygrid=False):
    ax.tick_params(direction="out", pad=2)
    if ygrid:
        ax.grid(axis="y", color=GRID, lw=0.45, zorder=0)
    ax.set_axisbelow(True)


def pastel(color, amount=0.48):
    rgb = np.asarray(mpl.colors.to_rgb(color))
    return tuple(rgb + (1.0 - rgb) * amount)


def label_bars(ax, bars, span, fmt="{:.1f}"):
    """Add compact value labels with enough clearance above every bar."""
    for index, bar in enumerate(bars):
        value = bar.get_height()
        # Alternate label heights within each dark/light pair so nearly equal
        # values remain separately readable.
        offset = (0.065 if index % 2 == 0 else 0.015) * span
        x_nudge = -0.025 if index % 2 == 0 else 0.025
        ax.text(
            bar.get_x() + bar.get_width() / 2 + x_nudge,
            value + offset,
            fmt.format(value),
            ha="center",
            va="bottom",
            fontsize=5.4,
            color=INK,
            clip_on=False,
            zorder=5,
        )


def grouped_bar_legend_handles():
    handles = []
    # Matplotlib fills multi-row legends down each column. Interleaving the
    # datasets therefore yields a dark Srinivas row and a light Chandak row.
    for error_label, color in zip(["Sub.", "Paired indel", "Mixed"], ERROR_COLORS):
        for dataset in DATASETS:
            is_dark = dataset == DATASETS[0]
            handles.append(
                Patch(
                    facecolor=color if is_dark else pastel(color),
                    edgecolor=INK if is_dark else color,
                    linewidth=0.4 if is_dark else 0.65,
                    label=f"{DATASET_LABELS[dataset]} · {error_label}",
                )
            )
    return handles


def plot_rank(ax, frame, dataset, letter):
    d = frame[(frame.record_type == "gate") & (frame.dataset_short == dataset)].sort_values("decile")
    total = int(d.clusters.sum())
    x = d.decile.to_numpy()
    guide = d.guide_failure_percent.to_numpy()
    final = d.final_failure_percent.to_numpy()
    ax.fill_between(
        x,
        final,
        guide,
        color=REDUCTION_FILL,
        alpha=0.62,
        linewidth=0,
        label="Failure reduction",
        zorder=1,
    )
    ax.plot(
        x,
        guide,
        color=BLUE,
        marker="o",
        ms=2.8,
        lw=1.05,
        label="Guide failure",
        zorder=3,
    )
    ax.plot(
        x,
        final,
        color=RED,
        marker="s",
        ms=2.6,
        lw=1.05,
        label="Final failure",
        zorder=3,
    )
    ax.set_xlim(1, 10)
    ax.set_xticks(np.arange(1, 11))
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel(r"Anomaly-score $a_i$ decile")
    ax.set_ylabel("Failure rate (%)")
    auc = d.auc.iloc[0]
    selected = d.activated_percent.iloc[0]
    ax.set_title(f"{DATASET_LABELS[dataset]}  |  N = {total:,}", loc="left", pad=3)
    ax.text(
        0.80,
        0.76,
        f"AUC {auc:.3f}\nDecoded {selected:.0f}%",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.3,
        linespacing=1.35,
    )
    ax.legend(loc="upper left", ncol=3, handlelength=1.2, columnspacing=0.75, borderaxespad=0)
    clean_axis(ax, ygrid=True)


def grouped_metric(ax, phase, metric, title, ylim=(0, 100), yticks=None):
    x = np.arange(len(ERROR_ORDER))
    values = {
        dataset: phase[phase.dataset_short == dataset]
        .set_index("error_type")
        .loc[ERROR_ORDER, metric]
        .to_numpy()
        for dataset in DATASETS
    }
    width = 0.32
    bars = []
    for j, color in enumerate(ERROR_COLORS):
        bars.extend(ax.bar(
            x[j] - width / 2,
            values[DATASETS[0]][j],
            width=width,
            color=color,
            edgecolor=INK,
            linewidth=0.4,
            zorder=3,
        ))
        bars.extend(ax.bar(
            x[j] + width / 2,
            values[DATASETS[1]][j],
            width=width,
            color=pastel(color),
            edgecolor=color,
            linewidth=0.65,
            zorder=3,
        ))
    ax.set_xticks(x, ERROR_LABELS)
    span = ylim[1] - ylim[0]
    ax.set_ylim(ylim[0], ylim[1] + 0.08 * span)
    ax.set_yticks(np.linspace(ylim[0], ylim[1], 3) if yticks is None else yticks)
    ax.set_title(title, loc="left", pad=3)
    for tick, color in zip(ax.get_xticklabels(), ERROR_COLORS):
        tick.set_color(color)
    label_bars(ax, bars, span)
    clean_axis(ax, ygrid=True)


def plot_routing(ax, gate):
    for dataset in DATASETS:
        d = gate[gate.dataset_short == dataset].sort_values("decile")
        ax.plot(
            d.decile,
            d.structured_decoding_percent,
            color=BLUE_LIGHT if dataset == "Srinivas" else PEACH,
            marker="o" if dataset == "Srinivas" else "s",
            ms=2.5,
            lw=1.0,
            label=DATASET_LABELS[dataset],
        )
    ax.set_xlim(1, 10)
    ax.set_xticks([1, 4, 7, 10])
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 50, 100])
    ax.set_xlabel(r"Anomaly-score $a_i$ decile")
    ax.set_ylabel("Decoded (%)")
    ax.set_title("Routing activation", loc="left", pad=3)
    ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.45), handlelength=1.3)
    clean_axis(ax, ygrid=True)


def plot_residual_share(ax, phase):
    x = np.arange(len(ERROR_ORDER))
    width = 0.32
    bars = []
    for j, (error, color) in enumerate(zip(ERROR_ORDER, ERROR_COLORS)):
        for i, dataset in enumerate(DATASETS):
            row = phase[(phase.dataset_short == dataset) & (phase.error_type == error)].iloc[0]
            bar = ax.bar(
                x[j] + (i - 0.5) * width,
                row.share_percent,
                width=width,
                color=color if i == 0 else pastel(color),
                edgecolor=INK if i == 0 else color,
                linewidth=0.4 if i == 0 else 0.65,
                zorder=3,
            )[0]
            bars.append(bar)
    ax.set_xticks(x, ERROR_LABELS)
    ax.set_ylim(0, 90)
    ax.set_yticks([0, 40, 80])
    ax.set_ylabel("Share (%)")
    ax.set_title("Residual-error composition", loc="left", pad=3)
    label_bars(ax, bars, 90)
    clean_axis(ax, ygrid=True)


def plot_hardest_decile(ax, gate):
    x = np.arange(len(DATASETS))
    d10 = gate[gate.decile == 10].set_index("dataset_short").loc[DATASETS]
    guide = d10.guide_failure_percent.to_numpy()
    final = d10.final_failure_percent.to_numpy()
    for i in range(len(DATASETS)):
        ax.plot([i, i], [final[i], guide[i]], color=LIGHT, lw=5.5, solid_capstyle="round", zorder=1)
        ax.scatter(i, guide[i], s=17, color=BLUE, marker="o", zorder=3)
        ax.scatter(i, final[i], s=16, color=RED, marker="s", zorder=3)
        reduction = 100 * (guide[i] - final[i]) / guide[i]
        ax.text(i, (guide[i] + final[i]) / 2, f"−{reduction:.0f}%", ha="center", va="center", fontsize=5.7, fontweight="bold")
    ax.set_xticks(x, ["Srinivas", "Chandak"])
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_ylabel("Failure rate (%)")
    ax.set_title("Hardest-decile rescue", loc="left", pad=3)
    ax.scatter([], [], s=17, color=BLUE, marker="o", label="Guide")
    ax.scatter([], [], s=16, color=RED, marker="s", label="Final")
    ax.legend(loc="lower left", ncol=2, handletextpad=0.35, columnspacing=0.7)
    clean_axis(ax, ygrid=True)


def plot_outcomes(ax, phase):
    labels = []
    y = []
    exact = []
    partial = []
    unchanged = []
    worsened = []
    for dataset in DATASETS:
        d = phase[phase.dataset_short == dataset].set_index("error_type").loc[ERROR_ORDER]
        prefix = "Sri" if dataset == "Srinivas" else "Cha"
        for error, short in zip(ERROR_ORDER, ["Sub", "Indel", "Mix"]):
            labels.append(f"{prefix} {short}")
            y.append(len(y))
            exact.append(d.loc[error, "exact_rescue_percent"])
            partial.append(d.loc[error, "partial_improvement_percent"])
            unchanged.append(d.loc[error, "unchanged_percent"])
            worsened.append(d.loc[error, "worsened_percent"])
    y = np.asarray(y)
    exact = np.asarray(exact)
    partial = np.asarray(partial)
    unchanged = np.asarray(unchanged)
    worsened = np.asarray(worsened)
    for values, left, color, label in [
        (exact, np.zeros_like(exact), BLUE, "Exact"),
        (partial, exact, GOLD, "Partial"),
        (unchanged, exact + partial, TEAL, "Unchanged"),
        (worsened, exact + partial + unchanged, RED, "Worsened"),
    ]:
        ax.barh(y, values, left=left, height=0.64, color=color, edgecolor="white", linewidth=0.25, label=label)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("Outcome (%)")
    ax.set_title("Phase-DP outcome composition", loc="left", pad=3)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.13), ncol=4, handlelength=1.0, columnspacing=0.75)
    clean_axis(ax, ygrid=False)


def main():
    df = pd.read_csv(DATA_PATH)
    gate = df[df.record_type == "gate"].copy()
    phase = df[df.record_type == "phase"].copy()

    fig = plt.figure(figsize=(7.25, 6.45), constrained_layout=False)
    grid = GridSpec(
        3,
        4,
        figure=fig,
        height_ratios=[1.06, 0.92, 1.06],
        wspace=0.68,
        hspace=0.62,
    )
    ax_a = fig.add_subplot(grid[0, 0:2])
    ax_b = fig.add_subplot(grid[0, 2:4])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])
    ax_e = fig.add_subplot(grid[1, 2])
    ax_f = fig.add_subplot(grid[1, 3])
    ax_g = fig.add_subplot(grid[2, 0])
    ax_h = fig.add_subplot(grid[2, 1])
    ax_i = fig.add_subplot(grid[2, 2:4])

    plot_rank(ax_a, df, "Srinivas", "a")
    plot_rank(ax_b, df, "Chandak-108", "b")
    plot_routing(ax_c, gate)
    plot_residual_share(ax_d, phase)
    grouped_metric(ax_e, phase, "exact_rescue_percent", "Exact rescue", ylim=(50, 100), yticks=[50, 75, 100])
    grouped_metric(ax_f, phase, "any_improvement_percent", "Any improvement", ylim=(65, 100), yticks=[70, 85, 100])
    grouped_metric(ax_g, phase, "ed_reduction_percent", "Edit-distance reduction", ylim=(60, 90), yticks=[60, 75, 90])
    ax_g.set_title(
        "Edit-distance reduction",
        loc="left",
        y=1.08,
        pad=0,
    )
    ax_g.set_ylabel("Reduction (%)")
    grouped_metric(ax_h, phase, "partial_improvement_percent", "Partial improvement", ylim=(0, 40), yticks=[0, 20, 40])
    ax_h.set_title(
        "Partial improvement",
        loc="left",
        y=1.08,
        pad=0,
    )
    plot_outcomes(ax_i, phase)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.965, bottom=0.085)
    for ax, letter in zip([ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g, ax_h, ax_i], "abcdefghi"):
        panel_label(fig, ax, letter)
    fig.legend(
        handles=grouped_bar_legend_handles(),
        loc="upper center",
        bbox_to_anchor=(0.64, 0.695),
        ncol=3,
        handlelength=1.0,
        handletextpad=0.45,
        columnspacing=0.85,
        labelspacing=0.45,
    )
    fig.savefig(OUT_PREFIX.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(OUT_PREFIX.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT_PREFIX.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(OUT_PREFIX.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


if __name__ == "__main__":
    main()
