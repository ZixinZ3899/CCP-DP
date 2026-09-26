#!/usr/bin/env python3
"""Create manuscript and supplement outputs from CCP-DP sensitivity results.

The input must be sensitivity_results/summary_metrics.csv, not the smoke-test
summary.  The script produces publication figures, source-data CSV files,
LaTeX tables, and an automatically calculated textual summary.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PARAMETERS = [
    ("coverage_threshold", "Coverage switch", "Reads", 40.0, "40"),
    ("reuse_threshold", "Reuse threshold", "Threshold", 0.50, "0.50"),
    ("phase_thresholds", "Phase-load thresholds", "Low / high", 1.0, "3/6"),
    ("cross_weight_scale", "Cross-weight scale", "Scale", 1.0, "1.00"),
    ("route_min_fraction", "Minimum route fraction", "Fraction", 0.10, "0.10"),
    ("phase_band", "Phase band", "Band width", 6.0, "6"),
]
PARAM_INFO = {p[0]: p for p in PARAMETERS}
MAIN_PARAMETERS = ["coverage_threshold", "cross_weight_scale", "route_min_fraction"]
COLORS = {
    "Srinivas-110": "#0072B2",
    "Chandak-108": "#E69F00",
    "Erlich-152": "#009E73",
}
MARKERS = {"Srinivas-110": "o", "Chandak-108": "s", "Erlich-152": "D"}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in (None, "") else math.nan


def validate(rows: list[dict[str, str]], allow_single_repeat: bool) -> list[str]:
    required = {
        "dataset", "parameter", "parameter_value", "plot_x", "repeats",
        "success_percent_mean", "mean_ed_mean", "elapsed_seconds_mean",
        "elapsed_seconds_sd", "decoded_fraction_percent_mean",
    }
    if not rows:
        raise SystemExit("Input CSV is empty.")
    missing = required - set(rows[0])
    if missing:
        raise SystemExit("Input CSV is missing columns: " + ", ".join(sorted(missing)))

    datasets = list(dict.fromkeys(row["dataset"] for row in rows))
    counts = {dataset: sum(row["dataset"] == dataset for row in rows)
              for dataset in datasets}
    bad_counts = {d: n for d, n in counts.items() if n < 18}
    if bad_counts:
        details = ", ".join(f"{d}={n} rows" for d, n in bad_counts.items())
        raise SystemExit(
            "This is not the complete sensitivity summary. "
            f"Expected at least 18 rows per dataset; found {details}. "
            "You probably selected smoke_results/summary_metrics.csv. Use "
            "sensitivity_results/summary_metrics.csv instead.")

    if not allow_single_repeat:
        bad_repeats = sorted({int(float(row["repeats"])) for row in rows
                              if int(float(row["repeats"])) < 3})
        if bad_repeats:
            raise SystemExit(
                "Final manuscript output requires repeats >= 3. "
                f"Found repeat counts {bad_repeats}. Do not use smoke results.")

    for dataset in datasets:
        baseline = [row for row in rows
                    if row["dataset"] == dataset and row["parameter"] == "baseline"]
        if len(baseline) != 1:
            raise SystemExit(f"Expected exactly one baseline row for {dataset}.")
        present = {row["parameter"] for row in rows if row["dataset"] == dataset}
        absent = set(PARAM_INFO) - present
        if absent:
            raise SystemExit(f"Missing parameter scans for {dataset}: {sorted(absent)}")
    return datasets


def baseline_map(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["dataset"]: row for row in rows if row["parameter"] == "baseline"}


def expanded_points(rows: list[dict[str, str]], dataset: str,
                    parameter: str) -> list[dict[str, float | str]]:
    base = next(row for row in rows
                if row["dataset"] == dataset and row["parameter"] == "baseline")
    _, _, _, default_x, default_label = PARAM_INFO[parameter]

    def point(row: dict[str, str], x: float, label: str) -> dict[str, float | str]:
        return {
            "dataset": dataset,
            "parameter": parameter,
            "parameter_value": label,
            "plot_x": x,
            "success": number(row, "success_percent_mean"),
            "success_sd": number(row, "success_percent_sd"),
            "mean_ed": number(row, "mean_ed_mean"),
            "mean_ed_sd": number(row, "mean_ed_sd"),
            "reconstruction": number(row, "reconstruction_percent_mean"),
            "runtime": number(row, "elapsed_seconds_mean"),
            "runtime_sd": number(row, "elapsed_seconds_sd"),
            "decoded": number(row, "decoded_fraction_percent_mean"),
            "changed": number(row, "changed_clusters_mean"),
            "routing_policy": row.get("routing_policy", ""),
        }

    points = [point(base, default_x, default_label)]
    points.extend(
        point(row, number(row, "plot_x"), row["parameter_value"])
        for row in rows
        if row["dataset"] == dataset and row["parameter"] == parameter
    )
    points.sort(key=lambda item: float(item["plot_x"]))
    b_success = number(base, "success_percent_mean")
    b_ed = number(base, "mean_ed_mean")
    b_runtime = number(base, "elapsed_seconds_mean")
    for item in points:
        item["success_delta"] = float(item["success"]) - b_success
        item["mean_ed_delta"] = float(item["mean_ed"]) - b_ed
        item["runtime_delta_percent"] = 100.0 * (
            float(item["runtime"]) / b_runtime - 1.0)
    return points


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D8DEE4", linewidth=0.65, alpha=0.8)
    ax.tick_params(labelsize=8.5, width=0.8, length=3)


def dataset_style(dataset: str, index: int) -> tuple[str, str]:
    fallback_colors = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
    fallback_markers = ["o", "s", "D", "^"]
    return (COLORS.get(dataset, fallback_colors[index % len(fallback_colors)]),
            MARKERS.get(dataset, fallback_markers[index % len(fallback_markers)]))


def save_figure(fig: plt.Figure, stem: Path) -> None:
    for suffix in ("svg", "pdf", "png"):
        fig.savefig(stem.with_suffix(f".{suffix}"), dpi=600,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_one(ax: plt.Axes, rows: list[dict[str, str]], datasets: list[str],
             parameter: str, metric: str, ylabel: str | None = None,
             runtime_error: bool = False) -> None:
    _, title, xlabel, default_x, _ = PARAM_INFO[parameter]
    for index, dataset in enumerate(datasets):
        points = expanded_points(rows, dataset, parameter)
        xs = [float(p["plot_x"]) for p in points]
        ys = [float(p[metric]) for p in points]
        color, marker = dataset_style(dataset, index)
        if runtime_error:
            errors = [float(p["runtime_sd"]) for p in points]
            ax.errorbar(xs, ys, yerr=errors, color=color, marker=marker,
                        linewidth=1.65, markersize=4.8, capsize=2.2,
                        label=dataset)
        else:
            ax.plot(xs, ys, color=color, marker=marker, linewidth=1.65,
                    markersize=4.8, label=dataset)
        if parameter == "phase_thresholds":
            ax.set_xticks(xs, [str(p["parameter_value"]) for p in points])
    ax.axvline(default_x, color="#6C757D", linewidth=0.9, linestyle="--")
    if metric in {"success_delta", "mean_ed_delta", "runtime_delta_percent"}:
        ax.axhline(0.0, color="#343A40", linewidth=0.75)
    ax.set_title(title, fontsize=9.8, fontweight="semibold", pad=6)
    ax.set_xlabel(xlabel, fontsize=8.8)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9.0)
    style_axis(ax)


def add_panel_labels(axes: list[plt.Axes]) -> None:
    for label, ax in zip("abcdefghijklmnopqrstuvwxyz", axes):
        ax.text(-0.13, 1.08, label, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="top", ha="left")


def figure_main(rows: list[dict[str, str]], datasets: list[str], out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.105, top=0.86,
                        wspace=0.28, hspace=0.46)
    plot_one(axes[0, 0], rows, datasets, "coverage_threshold", "success_delta",
             "Δ exact reconstruction (pp)")
    plot_one(axes[0, 1], rows, datasets, "cross_weight_scale", "success_delta")
    plot_one(axes[1, 0], rows, datasets, "route_min_fraction", "success_delta",
             "Δ exact reconstruction (pp)")
    plot_one(axes[1, 1], rows, datasets, "route_min_fraction", "runtime",
             "Runtime (s)", runtime_error=True)
    axes[1, 1].set_title("Runtime under selective routing", fontsize=9.8,
                         fontweight="semibold", pad=6)
    add_panel_labels(list(axes.flat))
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", frameon=False,
               ncol=max(1, len(labels)), fontsize=9.0,
               bbox_to_anchor=(0.5, 0.975))
    save_figure(fig, out / "Fig_parameter_sensitivity_main")


def figure_six(rows: list[dict[str, str]], datasets: list[str], out: Path,
               metric: str, ylabel: str, stem: str) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(10.8, 6.0))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.10, top=0.86,
                        wspace=0.28, hspace=0.48)
    for index, (ax, parameter_info) in enumerate(zip(axes.flat, PARAMETERS)):
        plot_one(ax, rows, datasets, parameter_info[0], metric,
                 ylabel if index in (0, 3) else None,
                 runtime_error=(metric == "runtime"))
    add_panel_labels(list(axes.flat))
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", frameon=False,
               ncol=max(1, len(labels)), fontsize=9.0,
               bbox_to_anchor=(0.5, 0.975))
    save_figure(fig, out / stem)


def latex_escape(value: str) -> str:
    return value.replace("_", r"\_").replace("%", r"\%")


def write_default_table(path: Path, rows: list[dict[str, str]],
                        datasets: list[str]) -> None:
    bases = baseline_map(rows)
    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\caption{Default-parameter results used in the sensitivity analysis.}",
        r"\label{tab:sensitivity-default}",
        r"\small",
        r"\begin{tabular}{lrrrr}",
        r"\hline",
        r"Dataset & Success (\%) & Mean ED & Time (s) & Decoded (\%) \\",
        r"\hline",
    ]
    for dataset in datasets:
        row = bases[dataset]
        lines.append(
            f"{latex_escape(dataset)} & {number(row, 'success_percent_mean'):.2f} & "
            f"{number(row, 'mean_ed_mean'):.4f} & "
            f"{number(row, 'elapsed_seconds_mean'):.3f} $\\pm$ "
            f"{number(row, 'elapsed_seconds_sd'):.3f} & "
            f"{number(row, 'decoded_fraction_percent_mean'):.1f} \\\\")
    lines.extend([r"\hline", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines))


def robustness_rows(rows: list[dict[str, str]], datasets: list[str]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for dataset in datasets:
        for parameter, title, _, _, _ in PARAMETERS:
            points = expanded_points(rows, dataset, parameter)
            output.append({
                "dataset": dataset,
                "parameter": parameter,
                "parameter_label": title,
                "max_abs_success_delta_pp": max(abs(float(p["success_delta"])) for p in points),
                "max_abs_mean_ed_delta": max(abs(float(p["mean_ed_delta"])) for p in points),
                "runtime_min_s": min(float(p["runtime"]) for p in points),
                "runtime_max_s": max(float(p["runtime"]) for p in points),
                "decoded_min_percent": min(float(p["decoded"]) for p in points),
                "decoded_max_percent": max(float(p["decoded"]) for p in points),
            })
    return output


def write_robustness_table(path: Path, summary: list[dict[str, object]]) -> None:
    lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Parameter sensitivity over the tested ranges.}",
        r"\label{tab:sensitivity-full}",
        r"\small",
        r"\begin{tabular}{llrrrr}",
        r"\hline",
        r"Dataset & Parameter & Max $|\Delta|$ success (pp) & Max $|\Delta|$ ED & Time range (s) & Decoded range (\%) \\",
        r"\hline",
    ]
    for row in summary:
        lines.append(
            f"{latex_escape(str(row['dataset']))} & "
            f"{latex_escape(str(row['parameter_label']))} & "
            f"{float(row['max_abs_success_delta_pp']):.3f} & "
            f"{float(row['max_abs_mean_ed_delta']):.4f} & "
            f"{float(row['runtime_min_s']):.3f}--{float(row['runtime_max_s']):.3f} & "
            f"{float(row['decoded_min_percent']):.1f}--{float(row['decoded_max_percent']):.1f} \\\\")
    lines.extend([r"\hline", r"\end{tabular}", r"\end{table*}", ""])
    path.write_text("\n".join(lines))


def write_text_summary(path: Path, summary: list[dict[str, object]]) -> None:
    global_max = max(summary, key=lambda row: float(row["max_abs_success_delta_pp"]))
    by_parameter: dict[str, float] = defaultdict(float)
    for row in summary:
        by_parameter[str(row["parameter_label"])] = max(
            by_parameter[str(row["parameter_label"])],
            float(row["max_abs_success_delta_pp"]))
    ordered = sorted(by_parameter.items(), key=lambda item: item[1], reverse=True)
    lines = [
        "AUTOMATIC NUMERICAL SUMMARY",
        "",
        ("Largest observed absolute success change: "
         f"{float(global_max['max_abs_success_delta_pp']):.3f} pp "
         f"({global_max['dataset']}, {global_max['parameter_label']})."),
        "",
        "Maximum absolute success change by parameter:",
    ]
    lines.extend(f"- {name}: {value:.3f} pp" for name, value in ordered)
    lines.extend([
        "",
        "Suggested manuscript statement:",
        ("Across the tested ranges, the largest absolute change in exact "
         "reconstruction was "
         f"{float(global_max['max_abs_success_delta_pp']):.2f} percentage points. "
         "The coverage switch and cross-cluster weight produced the most visible "
         "accuracy variation, whereas parameters whose tested values did not alter "
         "the active routing regime yielded identical reconstructions. Runtime and "
         "the fraction of structured decoding should be interpreted jointly with "
         "accuracy when selecting the routing floor."),
        "",
    ])
    path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary_csv", type=Path)
    parser.add_argument("--output", type=Path,
                        default=Path("parameter_sensitivity_paper_outputs"))
    parser.add_argument("--allow-single-repeat", action="store_true",
                        help="Development only; never use for manuscript outputs")
    args = parser.parse_args()

    rows = read_rows(args.summary_csv)
    datasets = validate(rows, args.allow_single_repeat)
    preferred = list(COLORS)
    datasets = ([name for name in preferred if name in datasets] +
                [name for name in datasets if name not in preferred])
    args.output.mkdir(parents=True, exist_ok=True)

    expanded: list[dict[str, object]] = []
    main_data: list[dict[str, object]] = []
    for dataset in datasets:
        for parameter, *_ in PARAMETERS:
            points = expanded_points(rows, dataset, parameter)
            expanded.extend(points)
            if parameter in MAIN_PARAMETERS:
                main_data.extend(points)

    write_csv(args.output / "source_data_main_figure.csv", main_data)
    write_csv(args.output / "source_data_all_parameters.csv", expanded)

    bases = baseline_map(rows)
    default_rows = [{
        "dataset": dataset,
        "success_percent": number(bases[dataset], "success_percent_mean"),
        "mean_ed": number(bases[dataset], "mean_ed_mean"),
        "reconstruction_percent": number(bases[dataset], "reconstruction_percent_mean"),
        "runtime_mean_s": number(bases[dataset], "elapsed_seconds_mean"),
        "runtime_sd_s": number(bases[dataset], "elapsed_seconds_sd"),
        "decoded_percent": number(bases[dataset], "decoded_fraction_percent_mean"),
        "routing_policy": bases[dataset].get("routing_policy", ""),
    } for dataset in datasets]
    write_csv(args.output / "default_results.csv", default_rows)

    robust = robustness_rows(rows, datasets)
    write_csv(args.output / "parameter_robustness_summary.csv", robust)
    write_default_table(args.output / "table_default_results.tex", rows, datasets)
    write_robustness_table(args.output / "table_full_sensitivity.tex", robust)
    write_text_summary(args.output / "results_summary.txt", robust)

    figure_main(rows, datasets, args.output)
    figure_six(rows, datasets, args.output, "success_delta",
               "Δ exact reconstruction (pp)", "Fig_S_success_all_parameters")
    figure_six(rows, datasets, args.output, "mean_ed_delta",
               "Δ mean ED", "Fig_S_mean_ed_all_parameters")
    figure_six(rows, datasets, args.output, "decoded",
               "Structured decoding (%)", "Fig_S_routing_all_parameters")
    figure_six(rows, datasets, args.output, "runtime_delta_percent",
               "Δ runtime (%)", "Fig_S_runtime_all_parameters")

    manifest = {
        "input": str(args.summary_csv.resolve()),
        "datasets": datasets,
        "input_rows": len(rows),
        "main_parameters": MAIN_PARAMETERS,
        "supplement_parameters": [item[0] for item in PARAMETERS],
    }
    (args.output / "output_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False))
    print(f"Created manuscript outputs in {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
