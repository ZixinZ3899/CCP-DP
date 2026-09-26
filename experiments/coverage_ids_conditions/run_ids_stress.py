#!/usr/bin/env python3
"""Run CCP-DP/BBS/CPL on generated stress-test conditions and evaluate them."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import shutil
import subprocess
import time
from pathlib import Path

try:
    from rapidfuzz.distance import Levenshtein

    def edit_distance(a: str, b: str) -> int:
        return Levenshtein.distance(a, b)
except ImportError:
    def edit_distance(a: str, b: str) -> int:
        """Dependency-free fallback; rapidfuzz is faster for the full run."""
        if len(a) < len(b):
            a, b = b, a
        previous = list(range(len(b) + 1))
        for i, ca in enumerate(a, start=1):
            current = [i]
            for j, cb in enumerate(b, start=1):
                current.append(min(
                    previous[j] + 1,
                    current[j - 1] + 1,
                    previous[j - 1] + (ca != cb),
                ))
            previous = current
        return previous[-1]


FIELDS = [
    "condition_id", "coverage", "total_indel_pct", "p_sub_pct", "p_ins_pct",
    "p_del_pct", "seed", "method", "clusters", "exact_count", "success_pct",
    "reconstruction_pct", "mean_ed", "runtime_s", "max_rss_kb", "output_path",
]


def read_sequences(path: Path) -> list[str]:
    return [x.strip().upper() for x in path.read_text().splitlines() if x.strip()]


def evaluate(centers_path: Path, output_path: Path) -> dict[str, float | int]:
    centers = read_sequences(centers_path)
    outputs = read_sequences(output_path)
    if len(centers) != len(outputs):
        raise RuntimeError(
            f"Output count mismatch for {output_path}: "
            f"centers={len(centers)}, outputs={len(outputs)}"
        )
    distances = [edit_distance(a, b) for a, b in zip(centers, outputs)]
    exact = sum(a == b for a, b in zip(centers, outputs))
    total_bases = sum(map(len, centers))
    return {
        "clusters": len(centers),
        "exact_count": exact,
        "success_pct": 100.0 * exact / len(centers),
        "reconstruction_pct": 100.0 * (1.0 - sum(distances) / total_bases),
        "mean_ed": sum(distances) / len(centers),
    }


def format_tokens(tokens: list[str], values: dict[str, str]) -> list[str]:
    return [token.format(**values) for token in tokens]


def parse_max_rss(path: Path) -> int:
    text = path.read_text(errors="replace")
    rss = re.search(
        r"^\s*Maximum resident set size \(kbytes\):\s*(\d+)\s*$",
        text, re.MULTILINE,
    )
    if not rss:
        raise RuntimeError(f"GNU time maximum RSS field missing from {path}")
    return int(rss.group(1))


def existing_keys(csv_path: Path) -> set[tuple[str, str]]:
    if not csv_path.exists():
        return set()
    with csv_path.open(newline="") as handle:
        return {(r["condition_id"], r["method"]) for r in csv.DictReader(handle)}


def upsert_row(csv_path: Path, row: dict[str, object]) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    prior: list[dict[str, str]] = []
    if csv_path.exists():
        with csv_path.open(newline="") as handle:
            prior = [
                old for old in csv.DictReader(handle)
                if not (
                    old["condition_id"] == str(row["condition_id"])
                    and old["method"] == str(row["method"])
                )
            ]
    temporary = csv_path.with_suffix(csv_path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(prior)
        writer.writerow(row)
    temporary.replace(csv_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--methods", default="CCP-DP,BBS,CPL")
    parser.add_argument("--max-conditions", type=int, default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    config = json.loads(args.config.read_text())
    jobs = int(config.get("jobs", 32))
    selected = {x.strip() for x in args.methods.split(",") if x.strip()}
    algorithms = [
        a for a in config["algorithms"] if a["name"] in selected and a.get("enabled", True)
    ]
    unknown = selected - {a["name"] for a in algorithms}
    if unknown:
        raise ValueError(f"Methods absent or disabled in config: {sorted(unknown)}")

    metadata_paths = sorted(args.data_root.glob("*/metadata.json"))
    if args.max_conditions is not None:
        metadata_paths = metadata_paths[: args.max_conditions]
    if not metadata_paths:
        raise RuntimeError(f"No conditions found below {args.data_root}")

    metrics_csv = args.results_root / "raw_metrics.csv"
    completed = set() if args.force else existing_keys(metrics_csv)
    time_binary = config.get("gnu_time", "/usr/bin/time")
    resolved_time = (
        str(Path(time_binary).resolve()) if Path(time_binary).exists()
        else shutil.which(time_binary)
    )
    if not resolved_time:
        print("[note] GNU time not found; runtime will use Python wall clock and RSS will be blank")

    for condition_index, metadata_path in enumerate(metadata_paths, start=1):
        meta = json.loads(metadata_path.read_text())
        condition_dir = metadata_path.parent
        result_dir = args.results_root / meta["condition_id"]
        result_dir.mkdir(parents=True, exist_ok=True)

        for algorithm in algorithms:
            key = (meta["condition_id"], algorithm["name"])
            if key in completed:
                print(f"[skip] {key[0]} / {key[1]}")
                continue

            output_path = result_dir / f"{algorithm['slug']}.txt"
            diag_path = result_dir / f"{algorithm['slug']}.diag.csv"
            log_path = result_dir / f"{algorithm['slug']}.log"
            values = {
                "clusters": str((condition_dir / "Clusters.txt").resolve()),
                "centers": str((condition_dir / "Centers.txt").resolve()),
                "output": str(output_path.resolve()),
                "diag": str(diag_path.resolve()),
                "length": str(meta["target_length"]),
                "separator": str(meta["separator"]),
                "jobs": str(jobs),
                "repo_root": str(
                    Path(__file__).resolve().parents[2]
                ),
            }
            command = format_tokens(algorithm["command"], values)
            timed_command = [resolved_time, "-v", *command] if resolved_time else command
            print(
                f"[{condition_index:03d}/{len(metadata_paths):03d}] "
                f"{meta['condition_id']} / {algorithm['name']}"
            )
            if args.dry_run:
                rendered = shlex.join(timed_command)
                if algorithm.get("cwd"):
                    rendered = (
                        f"cd {shlex.quote(algorithm['cwd'])} && " + rendered
                    )
                if algorithm.get("capture_stdout", False):
                    rendered += " > " + shlex.quote(str(output_path.resolve()))
                print("  " + rendered)
                continue

            started = time.perf_counter()
            with log_path.open("w") as log_handle:
                stdout_target = (
                    output_path.open("w") if algorithm.get("capture_stdout", False)
                    else log_handle
                )
                try:
                    proc = subprocess.run(
                        timed_command,
                        stdout=stdout_target,
                        stderr=log_handle,
                        cwd=algorithm.get("cwd") or None,
                        check=False,
                    )
                finally:
                    if hasattr(stdout_target, "close"):
                        stdout_target.close()
            if proc.returncode != 0:
                raise RuntimeError(
                    f"{algorithm['name']} failed for {meta['condition_id']}; "
                    f"inspect {log_path}"
                )
            if not output_path.exists():
                raise RuntimeError(f"Expected output was not created: {output_path}")

            python_wall = time.perf_counter() - started
            metrics = evaluate(condition_dir / "Centers.txt", output_path)
            if resolved_time:
                runtime_s, max_rss_kb = python_wall, parse_max_rss(log_path)
            else:
                runtime_s, max_rss_kb = python_wall, ""
            row = {
                "condition_id": meta["condition_id"],
                "coverage": meta["coverage"],
                "total_indel_pct": 100 * meta["total_indel"],
                "p_sub_pct": 100 * meta["p_sub"],
                "p_ins_pct": 100 * meta["p_ins"],
                "p_del_pct": 100 * meta["p_del"],
                "seed": meta["seed"],
                "method": algorithm["name"],
                **metrics,
                "runtime_s": runtime_s,
                "max_rss_kb": max_rss_kb,
                "output_path": str(output_path.resolve()),
            }
            upsert_row(metrics_csv, row)
            completed.add(key)

    if args.dry_run:
        print("Dry run finished; no algorithms were executed.")
    else:
        print(f"Raw metrics written to {metrics_csv.resolve()}")


if __name__ == "__main__":
    main()
