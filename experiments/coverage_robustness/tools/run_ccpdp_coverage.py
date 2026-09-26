#!/usr/bin/env python3
"""Prepare nested coverage inputs, run CCP-DP, and collect paper-ready tables.

The runner is deliberately sequential: each reconstruction process already
uses ``--jobs 32``, so launching several processes together would oversubscribe
the machine and make the runtime comparison unreliable.

Only Python's standard library is required.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shlex
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from make_nested_coverage import read_centers, read_clusters


PASS = "PASS"
DNA = frozenset("ACGTN")
ELAPSED_RE = re.compile(
    r"elapsed\s+seconds\s*:\s*([0-9]+(?:\.[0-9]+)?)", re.IGNORECASE
)


def resolve_from_root(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def load_config(path: Path) -> dict[str, Any]:
    config = json.loads(path.read_text(encoding="utf-8"))
    root_value = Path(config["project_root"]).expanduser()
    root = (
        root_value.resolve()
        if root_value.is_absolute()
        else (path.parent / root_value).resolve()
    )
    config["_root"] = root
    for key in ("executable", "evaluator", "samples_root", "runs_root"):
        config[f"_{key}"] = resolve_from_root(root, config[key])
    for dataset in config["datasets"]:
        dataset["_clusters"] = resolve_from_root(root, dataset["clusters"])
        dataset["_centers"] = resolve_from_root(root, dataset["centers"])
    return config


def select_datasets(config: dict[str, Any], names: list[str] | None) -> list[dict[str, Any]]:
    datasets = config["datasets"]
    if not names:
        return datasets
    wanted = set(names)
    selected = [dataset for dataset in datasets if dataset["name"] in wanted]
    missing = sorted(wanted - {dataset["name"] for dataset in selected})
    if missing:
        raise ValueError(f"Unknown dataset name(s): {', '.join(missing)}")
    return selected


def printable_command(command: Iterable[str]) -> str:
    return shlex.join([str(part) for part in command])


def run_logged(command: list[str], stdout_path: Path, stderr_path: Path) -> tuple[int, float]:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with stdout_path.open("w", encoding="utf-8", newline="\n") as stdout, \
            stderr_path.open("w", encoding="utf-8", newline="\n") as stderr:
        completed = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    return completed.returncode, time.perf_counter() - started


def read_output_sequences(path: Path) -> list[str]:
    sequences: list[str] = []
    with path.open("r", encoding="utf-8-sig", newline=None) as handle:
        for line_number, raw in enumerate(handle, start=1):
            sequence = raw.strip().upper()
            if not sequence:
                continue
            illegal = sorted(set(sequence) - DNA)
            if illegal:
                raise ValueError(
                    f"{path}:{line_number}: invalid output characters {illegal}"
                )
            sequences.append(sequence)
    return sequences


def edit_distance(a: str, b: str) -> int:
    """Levenshtein distance using O(min(len(a), len(b))) memory."""
    if len(a) < len(b):
        a, b = b, a
    previous = list(range(len(b) + 1))
    for i, char_a in enumerate(a, start=1):
        current = [i]
        for j, char_b in enumerate(b, start=1):
            current.append(min(
                current[-1] + 1,
                previous[j] + 1,
                previous[j - 1] + (char_a != char_b),
            ))
        previous = current
    return previous[-1]


def evaluate_internal(result_path: Path, centers_path: Path) -> dict[str, Any]:
    outputs = read_output_sequences(result_path)
    centers = read_centers(centers_path)
    if len(outputs) != len(centers):
        raise ValueError(
            f"Output/center count mismatch: {len(outputs)} outputs versus "
            f"{len(centers)} centers"
        )
    distances = [edit_distance(output, center) for output, center in zip(outputs, centers)]
    total_bases = sum(len(center) for center in centers)
    exact = sum(distance == 0 for distance in distances)
    count = len(centers)
    return {
        "n_sequences": count,
        "exact_reconstructions": exact,
        "success_rate_percent": 100.0 * exact / count,
        "reconstruction_rate_percent": 100.0 * (1.0 - sum(distances) / total_bases),
        "mean_edit_distance": statistics.fmean(distances),
        "total_edit_distance": sum(distances),
        "max_edit_distance": max(distances, default=0),
    }


def parse_algorithm_seconds(stderr_path: Path) -> float | None:
    text = stderr_path.read_text(encoding="utf-8", errors="replace")
    matches = ELAPSED_RE.findall(text)
    return float(matches[-1]) if matches else None


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def sample_all(
    config: dict[str, Any], datasets: list[dict[str, Any]], sampler: Path,
    coverages: list[int], seeds: list[int], force: bool, dry_run: bool,
) -> None:
    for dataset in datasets:
        output_dir = config["_samples_root"] / dataset["name"]
        command = [
            sys.executable, str(sampler), "sample",
            "-i", str(dataset["_clusters"]),
            "-c", str(dataset["_centers"]),
            "-o", str(output_dir),
            "--coverages", *map(str, coverages),
            "--seeds", *map(str, seeds),
            "--min-original-coverage", str(max(coverages)),
        ]
        if force:
            command.append("--force")
        print(f"\n[SAMPLE] {dataset['label']}\n{printable_command(command)}")
        if dry_run:
            continue
        for source in (dataset["_clusters"], dataset["_centers"]):
            if not source.exists():
                raise FileNotFoundError(source)
        subprocess.run(command, check=True)


def one_run_paths(runs_root: Path, dataset: str, seed: int, coverage: int) -> dict[str, Path]:
    run_dir = runs_root / dataset / f"seed_{seed}" / f"coverage_{coverage}x"
    return {
        "dir": run_dir,
        "result": run_dir / "result.txt",
        "diag": run_dir / "diag.csv",
        "stdout": run_dir / "ccpdp.stdout.log",
        "stderr": run_dir / "ccpdp.stderr.log",
        "eval_stdout": run_dir / "evaluation.stdout.log",
        "eval_stderr": run_dir / "evaluation.stderr.log",
        "metrics": run_dir / "metrics.json",
        "command": run_dir / "command.txt",
    }


def run_all(
    config: dict[str, Any], datasets: list[dict[str, Any]], coverages: list[int],
    seeds: list[int], jobs: int, rerun: bool, dry_run: bool,
    keep_going: bool,
) -> None:
    executable: Path = config["_executable"]
    evaluator: Path = config["_evaluator"]
    total = len(datasets) * len(seeds) * len(coverages)
    ordinal = 0
    if not dry_run and not executable.exists():
        raise FileNotFoundError(executable)

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
                    str(executable),
                    "-i", str(clusters_path),
                    "-l", str(dataset["length"]),
                    "-s", config["separator"],
                    "-o", str(paths["result"]),
                    "--diag", str(paths["diag"]),
                    "--mode", "auto",
                    "--jobs", str(jobs),
                ]
                print(
                    f"\n[{ordinal}/{total}] RUN {dataset['label']} "
                    f"seed={seed}, coverage={coverage}x\n"
                    f"{printable_command(command)}"
                )
                if dry_run:
                    continue

                paths["dir"].mkdir(parents=True, exist_ok=True)
                paths["command"].write_text(
                    printable_command(command) + "\n", encoding="utf-8"
                )
                payload: dict[str, Any] = {
                    "status": "RUNNING",
                    "dataset": dataset["name"],
                    "dataset_label": dataset["label"],
                    "length": int(dataset["length"]),
                    "seed": seed,
                    "coverage": coverage,
                    "jobs": jobs,
                    "clusters_path": str(clusters_path),
                    "centers_path": str(centers_path),
                    "result_path": str(paths["result"]),
                    "diag_path": str(paths["diag"]),
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
                    returncode, wall_seconds = run_logged(
                        command, paths["stdout"], paths["stderr"]
                    )
                    payload.update({
                        "ccpdp_returncode": returncode,
                        "wall_time_seconds": wall_seconds,
                        "algorithm_time_seconds": parse_algorithm_seconds(paths["stderr"]),
                        "n_clusters": cluster_count,
                    })
                    if returncode != 0:
                        raise RuntimeError(
                            f"CCP-DP exited with code {returncode}; see {paths['stderr']}"
                        )
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
                    print(
                        "  success={:.4f}%  reconstruction={:.4f}%  "
                        "mean_ED={:.6f}  time={}".format(
                            payload["success_rate_percent"],
                            payload["reconstruction_rate_percent"],
                            payload["mean_edit_distance"],
                            (
                                f"{payload['algorithm_time_seconds']:.5f}s"
                                if payload["algorithm_time_seconds"] is not None
                                else "NA"
                            ),
                        )
                    )
                except Exception as exc:  # preserve a machine-readable failure row
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


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def collect_results(config: dict[str, Any]) -> None:
    metrics_paths = sorted(config["_runs_root"].glob("**/metrics.json"))
    if not metrics_paths:
        raise FileNotFoundError(f"No metrics.json under {config['_runs_root']}")
    records = [json.loads(path.read_text(encoding="utf-8")) for path in metrics_paths]
    long_fields = [
        "status", "dataset", "dataset_label", "length", "seed", "coverage",
        "jobs", "n_clusters", "n_sequences", "exact_reconstructions",
        "success_rate_percent", "reconstruction_rate_percent",
        "mean_edit_distance", "total_edit_distance", "max_edit_distance",
        "algorithm_time_seconds", "wall_time_seconds", "ccpdp_returncode",
        "external_evaluator_returncode", "result_path", "diag_path", "error",
    ]
    output_root: Path = config["_runs_root"]
    write_csv(output_root / "results_long.csv", long_fields, records)

    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        if record.get("status") == PASS:
            grouped[(record["dataset"], int(record["coverage"]))].append(record)

    summary_rows: list[dict[str, Any]] = []
    metrics = [
        "success_rate_percent", "reconstruction_rate_percent",
        "mean_edit_distance", "algorithm_time_seconds", "wall_time_seconds",
    ]
    for (dataset, coverage), group in sorted(grouped.items()):
        row: dict[str, Any] = {
            "dataset": dataset,
            "dataset_label": group[0]["dataset_label"],
            "length": group[0]["length"],
            "coverage": coverage,
            "n_seeds": len(group),
            "n_clusters": group[0].get("n_clusters"),
        }
        for metric in metrics:
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

    fixed = ["dataset", "dataset_label", "length", "coverage", "n_seeds", "n_clusters"]
    stat_fields = [
        f"{metric}_{suffix}" for metric in metrics for suffix in ("mean", "sd", "sem")
    ]
    write_csv(output_root / "results_summary.csv", fixed + stat_fields, summary_rows)

    paper_rows: list[dict[str, Any]] = []
    for row in summary_rows:
        paper_rows.append({
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
            "Time-32 (s)": (
                "{:.5f} ± {:.5f}".format(
                    row["algorithm_time_seconds_mean"],
                    row["algorithm_time_seconds_sd"],
                )
                if "algorithm_time_seconds_mean" in row else "NA"
            ),
        })
    paper_fields = [
        "Dataset", "Coverage", "Seeds", "Clusters", "Success (%)",
        "Reconstruction (%)", "Mean ED", "Time-32 (s)",
    ]
    write_csv(output_root / "paper_table.csv", paper_fields, paper_rows)
    print(f"\nCollected {len(records)} runs")
    print(f"  {output_root / 'results_long.csv'}")
    print(f"  {output_root / 'results_summary.csv'}")
    print(f"  {output_root / 'paper_table.csv'}")


def add_common_filters(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--datasets", nargs="+",
        help="dataset names from the config; default: all datasets",
    )
    parser.add_argument("--coverages", type=int, nargs="+")
    parser.add_argument("--seeds", type=int, nargs="+")
    parser.add_argument("--dry-run", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    default_config = Path(__file__).with_name("datasets_ccpdp.json")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=default_config)
    subparsers = parser.add_subparsers(dest="command", required=True)

    sample = subparsers.add_parser("sample", help="make paired nested inputs")
    add_common_filters(sample)
    sample.add_argument("--force", action="store_true")

    run = subparsers.add_parser("run", help="run CCP-DP and evaluate outputs")
    add_common_filters(run)
    run.add_argument("--jobs", type=int)
    run.add_argument("--rerun", action="store_true")
    run.add_argument("--stop-on-error", action="store_true")

    all_cmd = subparsers.add_parser("all", help="sample, run, then collect")
    add_common_filters(all_cmd)
    all_cmd.add_argument("--jobs", type=int)
    all_cmd.add_argument("--force", action="store_true")
    all_cmd.add_argument("--rerun", action="store_true")
    all_cmd.add_argument("--stop-on-error", action="store_true")

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
        jobs = getattr(args, "jobs", None) or int(config["jobs"])
        sampler = Path(__file__).with_name("make_nested_coverage.py")
        if args.command in ("sample", "all"):
            sample_all(
                config, datasets, sampler, coverages, seeds,
                force=args.force, dry_run=args.dry_run,
            )
        if args.command in ("run", "all"):
            run_all(
                config, datasets, coverages, seeds, jobs,
                rerun=args.rerun, dry_run=args.dry_run,
                keep_going=not args.stop_on_error,
            )
        if args.command == "all" and not args.dry_run:
            collect_results(config)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
