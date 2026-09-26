#!/usr/bin/env python3
"""Run one-at-a-time parameter sensitivity for CCP-DP.

Every non-baseline configuration changes exactly one parameter family.  The
script evaluates exact reconstruction and Levenshtein edit distance internally,
parses routing diagnostics from the CCP-DP log, and writes resumable JSON/CSV
results.  It never uses centers during reconstruction.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[1]
SOURCE = ROOT / "main_sensitivity.cpp"
HEADER = ROOT / "ccpdp_sensitivity.hpp"
DEFAULT_BINARY = REPO_ROOT / "build" / "ccpdp_sensitivity"


@dataclass(frozen=True)
class Setting:
    key: str
    value_label: str
    plot_x: float
    tag: str
    flags: tuple[str, ...] = ()


DEFAULTS = {
    "coverage_threshold": ("40", 40.0),
    "reuse_threshold": ("0.50", 0.50),
    "phase_thresholds": ("3/6", 1.0),
    "cross_weight_scale": ("1.00", 1.0),
    "route_min_fraction": ("0.10", 0.10),
    "phase_band": ("auto (6 here)", 6.0),
}


def settings() -> list[Setting]:
    out = [Setting("baseline", "default", 0.0, "baseline")]
    out.extend(
        Setting("coverage_threshold", str(v), float(v), f"coverage_{v}",
                ("--coverage-threshold", str(v)))
        for v in (20, 30, 50, 60)
    )
    out.extend(
        Setting("reuse_threshold", f"{v:.2f}", v, f"reuse_{v:.2f}",
                ("--reuse-threshold", str(v)))
        for v in (0.30, 0.40, 0.60, 0.70)
    )
    phase_pairs = ((2.5, 5.0, 0.0), (3.5, 7.0, 2.0))
    out.extend(
        Setting("phase_thresholds", f"{low:g}/{high:g}", x,
                f"phase_{low:g}_{high:g}",
                ("--phase-low-threshold", str(low),
                 "--phase-high-threshold", str(high)))
        for low, high, x in phase_pairs
    )
    out.extend(
        Setting("cross_weight_scale", f"{v:.2f}", v, f"cross_scale_{v:.2f}",
                ("--cross-weight-scale", str(v)))
        for v in (0.75, 1.25)
    )
    out.extend(
        Setting("route_min_fraction", f"{v:.2f}", v, f"route_floor_{v:.2f}",
                ("--route-min-fraction", str(v)))
        for v in (0.05, 0.15, 0.20)
    )
    out.extend(
        Setting("phase_band", str(v), float(v), f"phase_band_{v}",
                ("--graph-phase-band", str(v)))
        for v in (4, 8)
    )
    return out


def compile_binary(binary: Path) -> None:
    if (binary.exists() and binary.stat().st_mtime >= SOURCE.stat().st_mtime
            and binary.stat().st_mtime >= HEADER.stat().st_mtime):
        return
    binary.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "g++", "-O3", "-march=native", "-std=c++17", "-fopenmp",
        "-Wall", "-Wextra", "-Wpedantic", str(SOURCE), "-o", str(binary),
    ]
    print("Compiling sensitivity binary...", flush=True)
    subprocess.run(command, check=True)


def load_sequences(path: Path) -> list[str]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        return []
    if any(line.startswith(">") for line in lines):
        sequences: list[str] = []
        current: list[str] = []
        for line in lines:
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                    current = []
            else:
                current.append(line)
        if current:
            sequences.append("".join(current).upper())
        return sequences
    return [line.upper() for line in lines]


def levenshtein_myers(pattern: str, text: str) -> int:
    """Exact Levenshtein distance using Python's arbitrary-width bit vectors."""
    if not pattern:
        return len(text)
    if not text:
        return len(pattern)
    if len(pattern) > len(text):
        pattern, text = text, pattern
    masks: dict[str, int] = {}
    for index, char in enumerate(pattern):
        masks[char] = masks.get(char, 0) | (1 << index)
    score = len(pattern)
    pv = ~0
    mv = 0
    last = 1 << (len(pattern) - 1)
    for char in text:
        eq = masks.get(char, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | ~(xh | pv)
        mh = pv & xh
        if ph & last:
            score += 1
        elif mh & last:
            score -= 1
        ph = (ph << 1) | 1
        mh <<= 1
        pv = mh | ~(xv | ph)
        mv = ph & xv
    return score


def evaluate(output_path: Path, center_path: Path, target_len: int,
             limit: int = -1) -> dict[str, Any]:
    outputs = load_sequences(output_path)
    centers = load_sequences(center_path)
    if limit > 0:
        centers = centers[:limit]
    if len(outputs) != len(centers):
        raise RuntimeError(
            f"Sequence count mismatch: output={len(outputs)}, centers={len(centers)}")
    distances = [levenshtein_myers(o, c) for o, c in zip(outputs, centers)]
    exact = sum(distance == 0 for distance in distances)
    total_ed = sum(distances)
    n = len(centers)
    return {
        "clusters": n,
        "exact": exact,
        "success_percent": 100.0 * exact / max(1, n),
        "mean_ed": total_ed / max(1, n),
        "reconstruction_percent":
            100.0 * (1.0 - total_ed / max(1, n * target_len)),
    }


def match_float(text: str, pattern: str, default: float = math.nan) -> float:
    found = re.search(pattern, text, flags=re.MULTILINE)
    return float(found.group(1)) if found else default


def match_int(text: str, pattern: str, default: int = -1) -> int:
    found = re.search(pattern, text, flags=re.MULTILINE)
    return int(found.group(1)) if found else default


def parse_log(text: str) -> dict[str, Any]:
    policy = re.search(r"^routing policy:\s*(.+)$", text, flags=re.MULTILINE)
    accepted = re.search(
        r"^cross calibration \(trust,accepted\):\s*[^,]+,(yes|no)$",
        text, flags=re.MULTILINE)
    weights = re.search(
        r"^trellis weights \(cross,local\):\s*([^,]+),([^\n]+)$",
        text, flags=re.MULTILINE)
    return {
        "elapsed_seconds": match_float(text, r"^elapsed seconds:\s*([0-9.eE+-]+)$"),
        "phase_load": match_float(text, r"^preliminary phase load:\s*([0-9.eE+-]+)$"),
        "total_error": match_float(text, r"^preliminary total error:\s*([0-9.eE+-]+)$"),
        "reuse": match_float(text, r"^cross-cluster reuse \(distinct guides\):\s*([0-9.eE+-]+)$"),
        "cross_trust": match_float(text, r"^cross calibration \(trust,accepted\):\s*([0-9.eE+-]+),"),
        "cross_accepted": accepted.group(1) if accepted else "unknown",
        "cross_weight": float(weights.group(1)) if weights else math.nan,
        "local_weight": float(weights.group(2)) if weights else math.nan,
        "decoded_clusters": match_int(text, r"^graph-decoded clusters:\s*(\d+)/"),
        "changed_clusters": match_int(text, r"^unified-graph changes:\s*(\d+)/"),
        "planned_clusters": match_int(text, r"^planned structured decoding:\s*(\d+)/"),
        "routing_policy": policy.group(1).strip() if policy else "unknown",
    }


def load_datasets(config_path: Path) -> list[dict[str, Any]]:
    data = json.loads(config_path.read_text())
    if not isinstance(data, list) or not data:
        raise RuntimeError("Config must contain a non-empty JSON list")
    for dataset in data:
        for key in ("name", "length", "clusters", "centers"):
            if key not in dataset:
                raise RuntimeError(f"Dataset entry missing {key}: {dataset}")
        for key in ("clusters", "centers"):
            path = Path(
                str(dataset[key]).format(repo_root=REPO_ROOT)
            ).expanduser()
            if "EDIT_ME" in str(path) or not path.is_file():
                raise RuntimeError(f"Please correct {key} path for {dataset['name']}: {path}")
            dataset[key] = str(path.resolve())
        dataset.setdefault("separator", "====")
    return data


def finite_or_blank(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    keys = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: finite_or_blank(row.get(key, "")) for key in keys})


def aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str, float], list[dict[str, Any]]] = {}
    for row in rows:
        key = (row["dataset"], row["parameter"], row["parameter_value"],
               float(row["plot_x"]))
        grouped.setdefault(key, []).append(row)
    metrics = (
        "success_percent", "mean_ed", "reconstruction_percent",
        "elapsed_seconds", "decoded_fraction_percent", "changed_clusters",
    )
    summary: list[dict[str, Any]] = []
    for key, group in sorted(grouped.items()):
        item: dict[str, Any] = {
            "dataset": key[0], "parameter": key[1],
            "parameter_value": key[2], "plot_x": key[3],
            "repeats": len(group),
            "routing_policy": group[0]["routing_policy"],
            "reuse": group[0]["reuse"],
            "phase_load": group[0]["phase_load"],
            "cross_trust": group[0]["cross_trust"],
            "cross_weight": group[0]["cross_weight"],
            "local_weight": group[0]["local_weight"],
        }
        for metric in metrics:
            values = [float(row[metric]) for row in group]
            item[f"{metric}_mean"] = statistics.fmean(values)
            item[f"{metric}_sd"] = statistics.stdev(values) if len(values) > 1 else 0.0
        summary.append(item)
    return summary


def run_one(binary: Path, dataset: dict[str, Any], setting: Setting,
            repeat: int, output_root: Path, jobs: int, limit: int,
            force: bool) -> dict[str, Any]:
    dataset_dir = output_root / dataset["name"] / setting.tag / f"repeat_{repeat:02d}"
    dataset_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = dataset_dir / "metrics.json"
    if metrics_path.exists() and not force:
        return json.loads(metrics_path.read_text())
    output_path = dataset_dir / "result.txt"
    diag_path = dataset_dir / "diag.csv"
    log_path = dataset_dir / "run.log"
    command = [
        str(binary), "-i", dataset["clusters"], "-l", str(dataset["length"]),
        "-s", dataset["separator"], "-o", str(output_path),
        "--diag", str(diag_path), "--mode", "auto", "--jobs", str(jobs),
        *setting.flags,
    ]
    if limit > 0:
        command.extend(("--limit", str(limit)))
    start = time.perf_counter()
    completed = subprocess.run(command, text=True, capture_output=True)
    wall = time.perf_counter() - start
    log_text = completed.stdout + completed.stderr
    log_path.write_text(log_text)
    if completed.returncode != 0:
        raise RuntimeError(f"CCP-DP failed; inspect {log_path}")
    metrics = evaluate(output_path, Path(dataset["centers"]),
                       int(dataset["length"]), limit)
    parsed = parse_log(log_text)
    metrics.update(parsed)
    metrics.update({
        "dataset": dataset["name"],
        "parameter": setting.key,
        "parameter_value": setting.value_label,
        "plot_x": setting.plot_x,
        "setting_tag": setting.tag,
        "repeat": repeat,
        "external_wall_seconds": wall,
        "decoded_fraction_percent":
            100.0 * parsed["decoded_clusters"] / max(1, metrics["clusters"]),
        "command": command,
    })
    metrics_path.write_text(json.dumps(metrics, indent=2, ensure_ascii=False))
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "datasets_server.json")
    parser.add_argument("--output", type=Path, default=ROOT / "sensitivity_results")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument(
        "--datasets", nargs="+",
        help="Run only the named datasets, e.g. Srinivas-110 Chandak-108")
    parser.add_argument("--jobs", type=int, default=32)
    parser.add_argument("--repeats", type=int, default=1,
                        help="Use 3 for final runtime error bars")
    parser.add_argument("--limit", type=int, default=-1,
                        help="Debug only; centers are truncated consistently")
    parser.add_argument("--smoke", action="store_true",
                        help="Run baseline only on the first 100 clusters")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    compile_binary(args.binary.resolve())
    datasets = load_datasets(args.config.resolve())
    if args.datasets:
        requested = set(args.datasets)
        available = {dataset["name"] for dataset in datasets}
        unknown = requested - available
        if unknown:
            raise RuntimeError(
                "Unknown dataset name(s): " + ", ".join(sorted(unknown)))
        datasets = [dataset for dataset in datasets
                    if dataset["name"] in requested]
    experiment_settings = settings()
    limit = args.limit
    if args.smoke:
        experiment_settings = experiment_settings[:1]
        limit = 100
    total = len(datasets) * len(experiment_settings) * args.repeats
    rows: list[dict[str, Any]] = []
    counter = 0
    for dataset in datasets:
        for setting in experiment_settings:
            for repeat in range(1, args.repeats + 1):
                counter += 1
                print(f"[{counter:03d}/{total:03d}] {dataset['name']} / "
                      f"{setting.tag} / repeat {repeat}", flush=True)
                row = run_one(args.binary.resolve(), dataset, setting, repeat,
                              args.output.resolve(), args.jobs, limit, args.force)
                rows.append(row)
    raw_rows = [{key: value for key, value in row.items() if key != "command"}
                for row in rows]
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "raw_metrics.csv", raw_rows)
    summary = aggregate(raw_rows)
    write_csv(args.output / "summary_metrics.csv", summary)
    (args.output / "run_manifest.json").write_text(json.dumps({
        "binary": str(args.binary.resolve()),
        "source": str(SOURCE),
        "datasets": datasets,
        "jobs": args.jobs,
        "repeats": args.repeats,
        "limit": limit,
        "defaults": DEFAULTS,
    }, indent=2, ensure_ascii=False))
    print(f"Raw metrics: {args.output / 'raw_metrics.csv'}")
    print(f"Summary:     {args.output / 'summary_metrics.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
