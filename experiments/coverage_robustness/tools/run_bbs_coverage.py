#!/usr/bin/env python3
"""Run official BBS on the same nested coverage inputs used by CCP-DP.

Each BBS process is wrapped by ``/usr/bin/time -v``. The script computes the
same exact-success, reconstruction-rate and edit-distance metrics as the
CCP-DP runner, preserves the external evaluator output, and produces long,
summary and paper-ready CSV tables.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

from run_ccpdp_coverage import (
    PASS,
    evaluate_internal,
    load_config,
    printable_command,
    read_centers,
    read_clusters,
    run_logged,
    select_datasets,
    write_csv,
    write_json,
)


def parse_duration(text: str) -> float:
    parts = text.strip().split(":")
    try:
        if len(parts) == 2:  # m:ss.xx
            return int(parts[0]) * 60.0 + float(parts[1])
        if len(parts) == 3:  # h:mm:ss.xx
            return int(parts[0]) * 3600.0 + int(parts[1]) * 60.0 + float(parts[2])
        return float(parts[0])
    except ValueError as exc:
        raise ValueError(f"Cannot parse GNU time duration: {text!r}") from exc


def parse_gnu_time(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"GNU time log was not created: {path}")
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if ": " in line:
            key, value = line.rsplit(": ", 1)
            values[key] = value.strip()

    def number(key: str) -> float | None:
        value = values.get(key)
        return float(value) if value is not None else None

    elapsed_key = "Elapsed (wall clock) time (h:mm:ss or m:ss)"
    elapsed = parse_duration(values[elapsed_key]) if elapsed_key in values else None
    rss = number("Maximum resident set size (kbytes)")
    cpu_text = values.get("Percent of CPU this job got", "").rstrip("%")
    return {
        "algorithm_time_seconds": elapsed,
        "gnu_time_elapsed_seconds": elapsed,
        "user_time_seconds": number("User time (seconds)"),
        "system_time_seconds": number("System time (seconds)"),
        "cpu_percent": float(cpu_text) if cpu_text else None,
        "max_rss_kb": int(rss) if rss is not None else None,
        "gnu_time_exit_status": (
            int(number("Exit status")) if number("Exit status") is not None else None
        ),
    }


def one_run_paths(runs_root: Path, dataset: str, seed: int, coverage: int) -> dict[str, Path]:
    run_dir = runs_root / dataset / f"seed_{seed}" / f"coverage_{coverage}x"
    return {
        "dir": run_dir,
        "result": run_dir / "result.txt",
        "stderr": run_dir / "bbs.stderr.log",
        "time": run_dir / "time_bbs.log",
        "eval_stdout": run_dir / "evaluation.stdout.log",
        "eval_stderr": run_dir / "evaluation.stderr.log",
        "metrics": run_dir / "metrics.json",
        "command": run_dir / "command.txt",
    }


def run_all(
    config: dict[str, Any], datasets: list[dict[str, Any]], coverages: list[int],
    seeds: list[int], threads: int, rerun: bool, dry_run: bool,
    keep_going: bool,
) -> None:
    executable: Path = config["_executable"]
    evaluator: Path = config["_evaluator"]
    time_value = Path(config.get("time_executable", "/usr/bin/time")).expanduser()
    time_executable = (
        time_value.resolve()
        if time_value.is_absolute()
        else (config["_root"] / time_value).resolve()
    )
    total = len(datasets) * len(seeds) * len(coverages)
    ordinal = 0
    if not dry_run:
        for required in (executable, time_executable):
            if not required.exists():
                raise FileNotFoundError(required)

    for dataset in datasets:
        for seed in seeds:
            for coverage in coverages:
                ordinal += 1
                inputs = (
                    config["_samples_root"] / dataset["name"] /
                    f"seed_{seed}" / f"coverage_{coverage}x"
                )
                clusters_path = inputs / "Clusters.txt"
                centers_path = inputs / "Centers.txt"
                paths = one_run_paths(
                    config["_runs_root"], dataset["name"], seed, coverage
                )
                if paths["metrics"].exists() and not rerun:
                    old = json.loads(paths["metrics"].read_text(encoding="utf-8"))
                    if old.get("status") == PASS:
                        print(
                            f"[{ordinal}/{total}] SKIP {dataset['name']} "
                            f"seed={seed} {coverage}x (already PASS)"
                        )
                        continue

                command = [
                    str(time_executable), "-v", "-o", str(paths["time"]),
                    str(executable), str(clusters_path),
                    "-l", str(dataset["length"]),
                    "-s", config["separator"],
                    "-t", str(threads),
                ]
                print(
                    f"\n[{ordinal}/{total}] RUN BBS | {dataset['label']} | "
                    f"seed={seed}, coverage={coverage}x\n"
                    f"{printable_command(command)} > {paths['result']}"
                )
                if dry_run:
                    continue

                paths["dir"].mkdir(parents=True, exist_ok=True)
                paths["command"].write_text(
                    f"{printable_command(command)} > {paths['result']}\n",
                    encoding="utf-8",
                )
                payload: dict[str, Any] = {
                    "status": "RUNNING",
                    "method": "BBS",
                    "dataset": dataset["name"],
                    "dataset_label": dataset["label"],
                    "length": int(dataset["length"]),
                    "seed": seed,
                    "coverage": coverage,
                    "threads": threads,
                    "clusters_path": str(clusters_path),
                    "centers_path": str(centers_path),
                    "result_path": str(paths["result"]),
                    "time_log_path": str(paths["time"]),
                }
                try:
                    if not clusters_path.exists() or not centers_path.exists():
                        raise FileNotFoundError(
                            f"Missing paired input: {clusters_path} or {centers_path}"
                        )
                    cluster_count = len(read_clusters(clusters_path))
                    center_count = len(read_centers(centers_path))
                    if cluster_count != center_count:
                        raise ValueError(
                            f"Sampled cluster/center mismatch: {cluster_count} vs "
                            f"{center_count}"
                        )
                    returncode, wrapper_wall = run_logged(
                        command, paths["result"], paths["stderr"]
                    )
                    payload.update({
                        "bbs_returncode": returncode,
                        "wrapper_wall_time_seconds": wrapper_wall,
                        "n_clusters": cluster_count,
                    })
                    if returncode != 0:
                        raise RuntimeError(
                            f"BBS exited with code {returncode}; see {paths['stderr']}"
                        )
                    payload.update(parse_gnu_time(paths["time"]))
                    payload["wall_time_seconds"] = payload["algorithm_time_seconds"]
                    payload.update(evaluate_internal(paths["result"], centers_path))

                    if evaluator.exists():
                        eval_command = [
                            sys.executable, str(evaluator),
                            "-o", str(paths["result"]),
                            "-a", str(centers_path),
                        ]
                        eval_code, eval_wall = run_logged(
                            eval_command, paths["eval_stdout"], paths["eval_stderr"]
                        )
                        payload.update({
                            "external_evaluator_returncode": eval_code,
                            "external_evaluator_wall_seconds": eval_wall,
                            "external_evaluator_command": printable_command(eval_command),
                        })
                    else:
                        payload["external_evaluator_returncode"] = None
                        payload["external_evaluator_note"] = (
                            f"Evaluator not found: {evaluator}; internal metrics were computed"
                        )
                    payload["status"] = PASS
                    rss_mib = (
                        payload["max_rss_kb"] / 1024.0
                        if payload.get("max_rss_kb") is not None else float("nan")
                    )
                    print(
                        "  success={:.4f}%  reconstruction={:.4f}%  "
                        "mean_ED={:.6f}  time={:.3f}s  maxRSS={:.1f}MiB".format(
                            payload["success_rate_percent"],
                            payload["reconstruction_rate_percent"],
                            payload["mean_edit_distance"],
                            payload["algorithm_time_seconds"],
                            rss_mib,
                        )
                    )
                except Exception as exc:
                    payload["status"] = "FAIL"
                    payload["error"] = str(exc)
                    write_json(paths["metrics"], payload)
                    print(f"  FAIL: {exc}", file=sys.stderr)
                    if not keep_going:
                        raise
                    continue
                write_json(paths["metrics"], payload)


def mean_sd_sem(values: list[float]) -> tuple[float, float, float]:
    mean = statistics.fmean(values)
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    return mean, sd, sd / math.sqrt(len(values))


def collect_results(config: dict[str, Any]) -> None:
    metrics_paths = sorted(config["_runs_root"].glob("**/metrics.json"))
    if not metrics_paths:
        raise FileNotFoundError(f"No metrics.json under {config['_runs_root']}")
    records = [json.loads(path.read_text(encoding="utf-8")) for path in metrics_paths]
    long_fields = [
        "status", "method", "dataset", "dataset_label", "length", "seed",
        "coverage", "threads", "n_clusters", "n_sequences",
        "exact_reconstructions", "success_rate_percent",
        "reconstruction_rate_percent", "mean_edit_distance",
        "total_edit_distance", "max_edit_distance", "algorithm_time_seconds",
        "user_time_seconds", "system_time_seconds", "cpu_percent", "max_rss_kb",
        "bbs_returncode", "gnu_time_exit_status", "external_evaluator_returncode",
        "result_path", "time_log_path", "error",
    ]
    output_root: Path = config["_runs_root"]
    write_csv(output_root / "results_long.csv", long_fields, records)

    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        if record.get("status") == PASS:
            grouped[(record["dataset"], int(record["coverage"]))].append(record)

    metric_names = [
        "success_rate_percent", "reconstruction_rate_percent",
        "mean_edit_distance", "algorithm_time_seconds", "max_rss_kb",
    ]
    summary_rows: list[dict[str, Any]] = []
    for (dataset, coverage), group in sorted(grouped.items()):
        row: dict[str, Any] = {
            "method": "BBS",
            "dataset": dataset,
            "dataset_label": group[0]["dataset_label"],
            "length": group[0]["length"],
            "coverage": coverage,
            "n_seeds": len(group),
            "n_clusters": group[0].get("n_clusters"),
        }
        for metric in metric_names:
            values = [
                float(record[metric]) for record in group
                if record.get(metric) is not None
            ]
            if values:
                mean, sd, sem = mean_sd_sem(values)
                row[f"{metric}_mean"] = mean
                row[f"{metric}_sd"] = sd
                row[f"{metric}_sem"] = sem
        summary_rows.append(row)

    fixed = [
        "method", "dataset", "dataset_label", "length", "coverage",
        "n_seeds", "n_clusters",
    ]
    stat_fields = [
        f"{metric}_{suffix}"
        for metric in metric_names for suffix in ("mean", "sd", "sem")
    ]
    write_csv(output_root / "results_summary.csv", fixed + stat_fields, summary_rows)

    paper_rows: list[dict[str, Any]] = []
    for row in summary_rows:
        paper_rows.append({
            "Method": "BBS",
            "Dataset": row["dataset_label"],
            "Coverage": f"{row['coverage']}x",
            "Seeds": row["n_seeds"],
            "Clusters": row["n_clusters"],
            "Success (%)": "{:.4f} ± {:.4f}".format(
                row["success_rate_percent_mean"], row["success_rate_percent_sd"]
            ),
            "Reconstruction (%)": "{:.4f} ± {:.4f}".format(
                row["reconstruction_rate_percent_mean"],
                row["reconstruction_rate_percent_sd"],
            ),
            "Mean ED": "{:.6f} ± {:.6f}".format(
                row["mean_edit_distance_mean"], row["mean_edit_distance_sd"]
            ),
            "Time-32 (s)": "{:.3f} ± {:.3f}".format(
                row["algorithm_time_seconds_mean"],
                row["algorithm_time_seconds_sd"],
            ),
            "Max RSS (MiB)": "{:.2f} ± {:.2f}".format(
                row["max_rss_kb_mean"] / 1024.0,
                row["max_rss_kb_sd"] / 1024.0,
            ),
        })
    paper_fields = [
        "Method", "Dataset", "Coverage", "Seeds", "Clusters", "Success (%)",
        "Reconstruction (%)", "Mean ED", "Time-32 (s)", "Max RSS (MiB)",
    ]
    write_csv(output_root / "paper_table.csv", paper_fields, paper_rows)
    print(f"\nCollected {len(records)} BBS runs")
    print(f"  {output_root / 'results_long.csv'}")
    print(f"  {output_root / 'results_summary.csv'}")
    print(f"  {output_root / 'paper_table.csv'}")


def add_filters(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--datasets", nargs="+")
    parser.add_argument("--coverages", type=int, nargs="+")
    parser.add_argument("--seeds", type=int, nargs="+")
    parser.add_argument("--dry-run", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=Path(__file__).with_name("datasets_bbs.json"),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run", help="run all requested BBS experiments")
    add_filters(run)
    run.add_argument("--threads", type=int)
    run.add_argument("--rerun", action="store_true")
    run.add_argument("--stop-on-error", action="store_true")
    subparsers.add_parser("collect", help="rebuild CSV tables from metrics.json")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    try:
        config = load_config(args.config.resolve())
        if args.command == "collect":
            collect_results(config)
            return
        datasets = select_datasets(config, args.datasets)
        coverages = sorted(set(args.coverages or config["coverages"]))
        seeds = sorted(set(args.seeds or config["seeds"]))
        threads = args.threads or int(config["threads"])
        run_all(
            config, datasets, coverages, seeds, threads,
            rerun=args.rerun, dry_run=args.dry_run,
            keep_going=not args.stop_on_error,
        )
        if not args.dry_run:
            collect_results(config)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
