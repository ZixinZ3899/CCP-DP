#!/usr/bin/env python3
"""Draw the dense Coverage x IDS main figure after merging both runs.

The main-text figure deliberately uses heatmaps only.  Absolute method curves
are retained in the exported CSV, but are not repeated in the figure.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np


EXPECTED_COVERAGES = [5, 6, 7, 8, 9, 10, 12, 15, 20]
METHODS = ["CCP-DP", "BBS", "CPL"]


def mean_sd(values: list[float]) -> tuple[float, float]:
    a = np.asarray(values, dtype=float)
    return float(a.mean()), float(a.std(ddof=1)) if len(a) > 1 else 0.0


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def aggregate_methods(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    groups: dict[tuple[int, float, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(int(row["coverage"]), float(row["total_indel_pct"]), row["method"])].append(row)
    output = []
    for (coverage, indel, method), group in sorted(groups.items()):
        success_mean, success_sd = mean_sd([float(x["success_pct"]) for x in group])
        ed_mean, ed_sd = mean_sd([float(x["mean_ed"]) for x in group])
        runtime_mean, runtime_sd = mean_sd([float(x["runtime_s"]) for x in group])
        output.append({
            "coverage": coverage,
            "total_indel_pct": indel,
            "method": method,
            "success_pct_mean": success_mean,
            "success_pct_sd": success_sd,
            "mean_ed_mean": ed_mean,
            "mean_ed_sd": ed_sd,
            "runtime_s_mean": runtime_mean,
            "runtime_s_sd": runtime_sd,
            "seeds": len(group),
        })
    return output


def paired_advantages(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    conditions: dict[tuple[int, float, int], dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        key = (int(row["coverage"]), float(row["total_indel_pct"]), int(row["seed"]))
        conditions[key][row["method"]] = row

    paired: dict[tuple[int, float], list[tuple[float, float]]] = defaultdict(list)
    for (coverage, indel, _seed), methods in conditions.items():
        if not set(METHODS).issubset(methods):
            continue
        ours, bbs, cpl = (methods[m] for m in METHODS)
        delta_success = float(ours["success_pct"]) - max(
            float(bbs["success_pct"]), float(cpl["success_pct"])
        )
        delta_ed = min(float(bbs["mean_ed"]), float(cpl["mean_ed"])) - float(ours["mean_ed"])
        paired[(coverage, indel)].append((delta_success, delta_ed))

    output = []
    for (coverage, indel), values in sorted(paired.items()):
        ds_mean, ds_sd = mean_sd([x[0] for x in values])
        de_mean, de_sd = mean_sd([x[1] for x in values])
        output.append({
            "coverage": coverage,
            "total_indel_pct": indel,
            "delta_success_mean": ds_mean,
            "delta_success_sd": ds_sd,
            "delta_ed_mean": de_mean,
            "delta_ed_sd": de_sd,
            "seeds": len(values),
        })
    return output


def lookup_matrix(
    rows: list[dict[str, object]], coverages: list[int], indels: list[float], key: str
) -> np.ndarray:
    lookup = {
        (int(r["coverage"]), float(r["total_indel_pct"])): float(r[key])
        for r in rows
    }
    return np.array(
        [[lookup.get((coverage, indel), np.nan) for coverage in coverages] for indel in indels],
        dtype=float,
    )


def annotate_heatmap(
    ax: plt.Axes, mean: np.ndarray, sd: np.ndarray,
    image: mpl.image.AxesImage, decimals: int, signed: bool,
) -> None:
    for row in range(mean.shape[0]):
        for col in range(mean.shape[1]):
            value = mean[row, col]
            if not np.isfinite(value):
                ax.text(col, row, "NA", ha="center", va="center", fontsize=6.8, color="#7A7A7A")
                continue
            rgba = image.cmap(image.norm(value))
            luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
            color = "white" if luminance < 0.54 else "#1F2933"
            sign = "+" if signed and value > 0 else ""
            ax.text(
                col, row - 0.12, f"{sign}{value:.{decimals}f}",
                ha="center", va="center", fontsize=6.9,
                fontweight="bold", color=color,
            )
            ax.text(
                col, row + 0.21, f"±{sd[row, col]:.{decimals}f}",
                ha="center", va="center", fontsize=5.6, color=color,
            )


def style_heatmap_axes(
    ax: plt.Axes, coverages: list[int], indels: list[float], title: str, ylabel: bool
) -> None:
    #ax.set_title(title, loc="left", fontsize=10.8, fontweight="bold", pad=8)
    ax.set_xticks(range(len(coverages)), [f"{x}×" for x in coverages])
    ax.set_yticks(range(len(indels)), [f"{x:g}%" for x in indels])
    ax.set_xlabel("Read coverage")
    if ylabel:
        ax.set_ylabel("Total indel rate")
    else:
        ax.tick_params(axis="y", left=False, labelleft=False)
    ax.set_xticks(np.arange(-0.5, len(coverages), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(indels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.6)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def mark_low_coverage_regime(ax: plt.Axes, coverages: list[int]) -> None:
    """Outline the prespecified low-coverage transition region (5--9x)."""
    selected = [i for i, coverage in enumerate(coverages) if 5 <= coverage <= 9]
    if not selected:
        return
    left = min(selected) - 0.48
    width = max(selected) - min(selected) + 0.96
    box = mpl.patches.Rectangle(
        (left, -0.48), width, 4.96,
        fill=False, edgecolor="#30343B", linewidth=1.05,
        linestyle=(0, (3, 2)), clip_on=False, zorder=5,
    )
    ax.add_patch(box)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--routing", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    raw = read_csv(args.metrics)
    method_summary = aggregate_methods(raw)
    advantages = paired_advantages(raw)
    routing = read_csv(args.routing)
    write_csv(args.outdir / "method_summary_dense.csv", method_summary)
    write_csv(args.outdir / "advantage_summary_dense.csv", advantages)

    found_coverages = sorted({int(x["coverage"]) for x in method_summary})
    # Keep the final nine-point layout even during a partial/preview run.
    # Missing conditions are shown as NA rather than silently removing columns.
    coverages = list(EXPECTED_COVERAGES)
    extras = [x for x in found_coverages if x not in EXPECTED_COVERAGES]
    coverages.extend(extras)
    indels = sorted({float(x["total_indel_pct"]) for x in advantages}, reverse=True)

    ds_mean = lookup_matrix(advantages, coverages, indels, "delta_success_mean")
    ds_sd = lookup_matrix(advantages, coverages, indels, "delta_success_sd")
    rt_mean = lookup_matrix(routing, coverages, indels, "decoded_fraction_pct_mean")
    rt_sd = lookup_matrix(routing, coverages, indels, "decoded_fraction_pct_sd")

    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })

    fig, axes = plt.subplots(
        1, 2, figsize=(12.2, 4.0), constrained_layout=True,
        gridspec_kw={"wspace": 0.08},
    )

    diverging = LinearSegmentedColormap.from_list(
        "gain", ["#3E6D9C", "#DCE7EE", "#F7F7F5", "#F4C1A1", "#D95843"]
    )
    diverging.set_bad("#ECEFF1")
    finite_gain = ds_mean[np.isfinite(ds_mean)]
    vmax = max(18.0, float(finite_gain.max()) if finite_gain.size else 18.0)
    vmin = min(-2.0, float(finite_gain.min()) if finite_gain.size else -2.0)
    im_gain = axes[0].imshow(
        np.ma.masked_invalid(ds_mean), aspect="auto", interpolation="nearest",
        cmap=diverging, norm=TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax),
    )
    style_heatmap_axes(
        axes[0], coverages, indels,
        "Exact-reconstruction advantage over the best baseline", True,
    )
    annotate_heatmap(axes[0], ds_mean, ds_sd, im_gain, decimals=1, signed=True)
    mark_low_coverage_regime(axes[0], coverages)
    axes[0].text(-0.105, 1.08, "a", transform=axes[0].transAxes,
                 fontsize=13, fontweight="bold", va="top")
    cb1 = fig.colorbar(im_gain, ax=axes[0], orientation="horizontal", pad=0.10, fraction=0.08)
    cb1.set_ticks([vmin, 0, 5, 10, 15, vmax])
    cb1.set_label("ΔSuccess (percentage points)")
    cb1.outline.set_visible(False)

    route_cmap = LinearSegmentedColormap.from_list(
        "routing", ["#F2F6F5", "#CBE3DD", "#6EB19F", "#137967"]
    )
    route_cmap.set_bad("#ECEFF1")
    finite_route = rt_mean[np.isfinite(rt_mean)]
    route_max = max(50.0, float(finite_route.max()) if finite_route.size else 50.0)
    im_route = axes[1].imshow(
        np.ma.masked_invalid(rt_mean), aspect="auto", interpolation="nearest",
        cmap=route_cmap, vmin=0.0, vmax=route_max,
    )
    style_heatmap_axes(axes[1], coverages, indels, "Adaptive structured-decoding allocation", False)
    annotate_heatmap(axes[1], rt_mean, rt_sd, im_route, decimals=1, signed=False)
    mark_low_coverage_regime(axes[1], coverages)
    axes[1].text(-0.105, 1.08, "b", transform=axes[1].transAxes,
                 fontsize=13, fontweight="bold", va="top")
    cb2 = fig.colorbar(im_route, ax=axes[1], orientation="horizontal", pad=0.10, fraction=0.08)
    cb2.set_ticks([0, 10, 20, 30, 40, route_max])
    cb2.set_label("Decoded clusters (%)")
    cb2.outline.set_visible(False)

    fig.text(
        0.5, -0.025,
        "Mean ± s.d. over five paired seeds. ΔSuccess = CCP-DP − max(BBS, CPL).",
        ha="center", va="top", fontsize=8, color="#4B5563",
    )

    stem = args.outdir / "coverage_ids_dense_main"
    fig.savefig(stem.with_suffix(".png"), dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Written {stem}.png/.svg/.pdf")


if __name__ == "__main__":
    main()
