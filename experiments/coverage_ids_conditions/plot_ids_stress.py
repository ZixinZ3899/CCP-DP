#!/usr/bin/env python3
"""Aggregate paired seeds and draw the Coverage x IDS stress-test heatmaps."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm


def mean_sd(values: list[float]) -> tuple[float, float]:
    a = np.asarray(values, dtype=float)
    return float(a.mean()), float(a.std(ddof=1)) if len(a) > 1 else 0.0


def read_metrics(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def build_paired(rows: list[dict[str, str]]) -> list[dict[str, float]]:
    by_condition: dict[tuple[int, float, int], dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        key = (int(row["coverage"]), float(row["total_indel_pct"]), int(row["seed"]))
        by_condition[key][row["method"]] = row

    paired = []
    required = {"CCP-DP", "BBS", "CPL"}
    for (coverage, indel, seed), methods in sorted(by_condition.items()):
        if not required.issubset(methods):
            continue
        ours = methods["CCP-DP"]
        bbs = methods["BBS"]
        cpl = methods["CPL"]
        paired.append({
            "coverage": coverage,
            "indel": indel,
            "seed": seed,
            "delta_success": float(ours["success_pct"]) - max(
                float(bbs["success_pct"]), float(cpl["success_pct"])
            ),
            "delta_ed": min(float(bbs["mean_ed"]), float(cpl["mean_ed"])) - float(ours["mean_ed"]),
        })
    if not paired:
        raise RuntimeError("No complete paired CCP-DP/BBS/CPL conditions were found")
    return paired


def aggregate(paired: list[dict[str, float]]) -> list[dict[str, float]]:
    groups: dict[tuple[int, float], list[dict[str, float]]] = defaultdict(list)
    for row in paired:
        groups[(int(row["coverage"]), row["indel"])].append(row)
    output = []
    for (coverage, indel), values in sorted(groups.items()):
        ds_mean, ds_sd = mean_sd([x["delta_success"] for x in values])
        de_mean, de_sd = mean_sd([x["delta_ed"] for x in values])
        output.append({
            "coverage": coverage,
            "indel": indel,
            "delta_success_mean": ds_mean,
            "delta_success_sd": ds_sd,
            "delta_ed_mean": de_mean,
            "delta_ed_sd": de_sd,
            "seeds": len(values),
        })
    return output


def matrix(summary: list[dict[str, float]], key: str,
           coverages: list[int], indels: list[float]) -> np.ndarray:
    lookup = {(r["coverage"], r["indel"]): r[key] for r in summary}
    return np.array([[lookup.get((c, i), np.nan) for c in coverages] for i in indels])


def symmetric_limit(values: np.ndarray) -> float:
    finite = np.abs(values[np.isfinite(values)])
    return max(float(finite.max()) if finite.size else 1.0, 1e-6)


def annotate(ax, means: np.ndarray, sds: np.ndarray, show_sd: bool) -> None:
    max_abs = symmetric_limit(means)
    for row in range(means.shape[0]):
        for col in range(means.shape[1]):
            value = means[row, col]
            if not np.isfinite(value):
                ax.text(col, row, "NA", ha="center", va="center", color="#777777")
                continue
            color = "white" if abs(value) > 0.55 * max_abs else "#20242A"
            label = f"{value:+.2f}"
            if show_sd:
                label += f"\n±{sds[row, col]:.2f}"
            ax.text(col, row, label, ha="center", va="center", color=color,
                    fontsize=8.3, fontweight="semibold")


def write_summary(path: Path, summary: list[dict[str, float]]) -> None:
    fields = list(summary[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)


def write_method_summary(path: Path, rows: list[dict[str, str]]) -> None:
    groups: dict[tuple[int, float, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(int(row["coverage"]), float(row["total_indel_pct"]), row["method"])].append(row)
    output = []
    for (coverage, indel, method), values in sorted(groups.items()):
        success_mean, success_sd = mean_sd([float(x["success_pct"]) for x in values])
        ed_mean, ed_sd = mean_sd([float(x["mean_ed"]) for x in values])
        runtime_mean, runtime_sd = mean_sd([float(x["runtime_s"]) for x in values])
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
            "seeds": len(values),
        })
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)


def write_routing_summary(path: Path, rows: list[dict[str, str]]) -> None:
    ccpdp = [r for r in rows if r["method"] == "CCP-DP"]
    if not ccpdp or "decoded_fraction_pct" not in ccpdp[0]:
        return
    groups: dict[tuple[int, float], list[dict[str, str]]] = defaultdict(list)
    for row in ccpdp:
        groups[(int(row["coverage"]), float(row["total_indel_pct"]))].append(row)
    output = []
    for (coverage, indel), values in sorted(groups.items()):
        decoded_values = [
            float(x["decoded_fraction_pct"])
            for x in values if x.get("decoded_fraction_pct", "") != ""
        ]
        changed_values = [
            float(x["changed_clusters"])
            for x in values if x.get("changed_clusters", "") != ""
        ]
        if not decoded_values:
            continue
        decoded_mean, decoded_sd = mean_sd(decoded_values)
        changed_mean, changed_sd = mean_sd(changed_values) if changed_values else (float("nan"), float("nan"))
        output.append({
            "coverage": coverage,
            "total_indel_pct": indel,
            "decoded_fraction_pct_mean": decoded_mean,
            "decoded_fraction_pct_sd": decoded_sd,
            "changed_clusters_mean": changed_mean,
            "changed_clusters_sd": changed_sd,
            "seeds": len(decoded_values),
        })
    if output:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(output[0]))
            writer.writeheader()
            writer.writerows(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--show-sd", action="store_true")
    parser.add_argument("--illustrative", action="store_true")
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = read_metrics(args.metrics)
    summary = aggregate(build_paired(rows))
    if not args.allow_partial:
        incomplete = [
            (int(r["coverage"]), float(r["indel"]), int(r["seeds"]))
            for r in summary if int(r["seeds"]) != 5
        ]
        if len(summary) != 20 or incomplete:
            raise RuntimeError(
                "Full plotting requires 20 Coverage x IDS cells with five paired "
                f"seeds each; found {len(summary)} cells, incomplete={incomplete}. "
                "Use --allow-partial only for diagnostics."
            )
    write_summary(args.outdir / "stress_summary.csv", summary)
    write_method_summary(args.outdir / "method_summary.csv", rows)
    write_routing_summary(args.outdir / "ccpdp_routing_summary.csv", rows)
    coverages = sorted({int(x["coverage"]) for x in summary})
    # High indel at the top makes the hypothesized difficult regime top-left.
    indels = sorted({float(x["indel"]) for x in summary}, reverse=True)

    success = matrix(summary, "delta_success_mean", coverages, indels)
    success_sd = matrix(summary, "delta_success_sd", coverages, indels)
    ed = matrix(summary, "delta_ed_mean", coverages, indels)
    ed_sd = matrix(summary, "delta_ed_sd", coverages, indels)

    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })
    cmap = LinearSegmentedColormap.from_list(
        "nature_diverging", ["#3C5488", "#B4D4E7", "#F7F7F5", "#F6B48F", "#E64B35"]
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.15), constrained_layout=True)

    for ax, values, sds, panel, title, cbar_label in [
        (axes[0], success, success_sd, "a", "Exact reconstruction advantage",
         r"$\Delta$Success (percentage points)"),
        (axes[1], ed, ed_sd, "b", "Residual-error advantage",
         r"$\Delta$Mean ED"),
    ]:
        limit = symmetric_limit(values)
        x_edges = np.arange(len(coverages) + 1) - 0.5
        y_edges = np.arange(len(indels) + 1) - 0.5
        image = ax.pcolormesh(
            x_edges, y_edges, values,
            cmap=cmap,
            norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
            shading="flat",
            edgecolors="white",
            linewidth=1.5,
        )
        ax.set_xlim(-0.5, len(coverages) - 0.5)
        ax.set_ylim(len(indels) - 0.5, -0.5)
        ax.set_xticks(range(len(coverages)), [f"{x}×" for x in coverages])
        ax.set_yticks(range(len(indels)), [f"{x:g}%" for x in indels])
        ax.set_xlabel("Read coverage")
        ax.set_ylabel("Total indel rate")
        ax.set_title(title, loc="left", fontsize=10.5, fontweight="semibold", pad=9)
        ax.text(-0.15, 1.08, panel, transform=ax.transAxes, fontsize=12,
                fontweight="bold", va="top")
        annotate(ax, values, sds, args.show_sd)
        cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.035)
        if getattr(cbar, "solids", None) is not None:
            cbar.solids.set_rasterized(False)
        cbar.set_label(cbar_label, fontsize=8.5)
        cbar.outline.set_linewidth(0.6)

    if args.illustrative:
        fig.text(
            0.5, 1.015, "Illustrative layout — not experimental results",
            ha="center", va="bottom", fontsize=8.5, color="#B23A2B",
            fontweight="semibold"
        )

    stem = args.outdir / "coverage_ids_stress"
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Written {stem}.png/.svg/.pdf")


if __name__ == "__main__":
    main()
