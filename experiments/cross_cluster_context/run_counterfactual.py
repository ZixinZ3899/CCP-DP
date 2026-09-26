#!/usr/bin/env python3
"""Run controlled Cross-LOO counterfactuals with frozen CCP-DP logic."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

try:
    from rapidfuzz.distance import Levenshtein

    def edit_distance(a: str, b: str) -> int:
        return Levenshtein.distance(a, b)
except ImportError:
    def edit_distance(a: str, b: str) -> int:
        if len(a) < len(b):
            a, b = b, a
        previous = list(range(len(b) + 1))
        for i, ca in enumerate(a, 1):
            current = [i]
            for j, cb in enumerate(b, 1):
                current.append(min(
                    previous[j] + 1,
                    current[j - 1] + 1,
                    previous[j - 1] + (ca != cb),
                ))
            previous = current
        return previous[-1]


@dataclass(frozen=True)
class Dataset:
    name: str
    clusters: Path
    centers: Path
    length: int
    separator: str = "==============================="


VARIANTS = {
    "local": "Local only",
    "table-permuted": "Table-permuted",
    "wrong-dataset": "Wrong-dataset",
    "correct": "Correct Cross-LOO",
}


def read_sequences(path: Path) -> list[str]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        raise RuntimeError(f"No sequences found in {path}")
    if any(line.startswith(">") for line in lines):
        result: list[str] = []
        current: list[str] = []
        for line in lines:
            if line.startswith(">"):
                if current:
                    result.append("".join(current).upper())
                    current = []
            else:
                current.append(line)
        if current:
            result.append("".join(current).upper())
        return result
    return [line.upper() for line in lines]


def evaluate(centers: list[str], predictions: list[str]) -> dict[str, float | int]:
    if len(centers) != len(predictions):
        raise RuntimeError(
            f"Sequence-count mismatch: centers={len(centers)}, "
            f"predictions={len(predictions)}"
        )
    distances = [edit_distance(a, b) for a, b in zip(centers, predictions)]
    exact = sum(a == b for a, b in zip(centers, predictions))
    total_bases = sum(map(len, centers))
    return {
        "n": len(centers),
        "exact": exact,
        "success": 100.0 * exact / len(centers),
        "reconstruction": 100.0 * (1.0 - sum(distances) / total_bases),
        "mean_ed": statistics.mean(distances),
    }


def sd(values: list[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else 0.0


def mcnemar_exact_two_sided(b: int, c: int) -> float:
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    probability = math.ldexp(1.0, -n)
    cumulative = probability
    for i in range(k):
        probability *= (n - i) / (i + 1)
        cumulative += probability
    return min(1.0, 2.0 * cumulative)


def run_once(
    binary: Path,
    dataset: Dataset,
    variant: str,
    repetition: int,
    outdir: Path,
    jobs: int,
    limit: int | None,
    wrong_guides: Path | None,
) -> dict[str, Path | float]:
    stem = f"{variant}_rep{repetition}"
    result = outdir / f"{stem}.txt"
    diag = outdir / f"{stem}.diag.csv"
    log = outdir / f"{stem}.log"
    guides = outdir / "frozen_guides.txt"
    command = [
        str(binary),
        "-i", str(dataset.clusters),
        "-l", str(dataset.length),
        "-s", dataset.separator,
        "-o", str(result),
        "--diag", str(diag),
        "--mode", "auto",
        "--jobs", str(jobs),
        "--decode-all",
        "--context-mode", variant,
        "--context-seed", str(2025 + repetition),
    ]
    if variant == "correct" and repetition == 1:
        command.extend(["--guide-output", str(guides)])
    if variant == "wrong-dataset":
        if wrong_guides is None:
            raise RuntimeError("wrong-dataset mode requires external guides")
        command.extend(["--context-guides", str(wrong_guides)])
    if limit is not None:
        command.extend(["--limit", str(limit)])
    start = time.perf_counter()
    with log.open("w") as handle:
        completed = subprocess.run(
            command, stdout=subprocess.DEVNULL, stderr=handle, check=False
        )
    elapsed = time.perf_counter() - start
    if completed.returncode != 0:
        raise RuntimeError(f"Run failed; inspect {log}")
    return {"result": result, "diag": diag, "log": log, "elapsed": elapsed}


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["smoke", "full"])
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=32)
    parser.add_argument("--permuted-repeats", type=int, default=5)
    args = parser.parse_args()

    root = args.project_root.resolve()
    binary = args.binary.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if not binary.exists():
        raise RuntimeError(f"Binary not found: {binary}")

    datasets = [
        Dataset(
            "Srinivas-110",
            root / "data/processed/srinivas_110/Clusters.txt",
            root / "data/processed/srinivas_110/Centers.txt",
            110,
        ),
        Dataset(
            "Chandak-108",
            root / "data/processed/chandak_108/Clusters.txt",
            root / "data/processed/chandak_108/Centers.txt",
            108,
        ),
    ]
    for dataset in datasets:
        if not dataset.clusters.exists() or not dataset.centers.exists():
            raise RuntimeError(f"Missing inputs for {dataset.name}")

    limit = 200 if args.mode == "smoke" else None
    records: dict[tuple[str, str], list[dict[str, Path | float]]] = {}

    # Correct runs first create the frozen guides used by the opposite-data
    # negative control.  All conditions use decode-all to isolate context.
    for dataset in datasets:
        dataset_dir = outdir / dataset.name
        dataset_dir.mkdir(parents=True, exist_ok=True)
        print(f"[{dataset.name}] correct")
        records[(dataset.name, "correct")] = [run_once(
            binary, dataset, "correct", 1, dataset_dir,
            args.jobs, limit, None,
        )]

    guide_paths = {
        dataset.name: outdir / dataset.name / "frozen_guides.txt"
        for dataset in datasets
    }
    opposite = {
        "Srinivas-110": guide_paths["Chandak-108"],
        "Chandak-108": guide_paths["Srinivas-110"],
    }

    for dataset in datasets:
        dataset_dir = outdir / dataset.name
        print(f"[{dataset.name}] local")
        records[(dataset.name, "local")] = [run_once(
            binary, dataset, "local", 1, dataset_dir,
            args.jobs, limit, None,
        )]
        permuted: list[dict[str, Path | float]] = []
        for repetition in range(1, args.permuted_repeats + 1):
            print(f"[{dataset.name}] table-permuted {repetition}/{args.permuted_repeats}")
            permuted.append(run_once(
                binary, dataset, "table-permuted", repetition,
                dataset_dir, args.jobs, limit, None,
            ))
        records[(dataset.name, "table-permuted")] = permuted
        print(f"[{dataset.name}] wrong-dataset")
        records[(dataset.name, "wrong-dataset")] = [run_once(
            binary, dataset, "wrong-dataset", 1, dataset_dir,
            args.jobs, limit, opposite[dataset.name],
        )]

    summary_rows: list[dict[str, object]] = []
    metrics_by_key: dict[tuple[str, str], list[dict[str, float | int]]] = {}
    for dataset in datasets:
        centers = read_sequences(dataset.centers)
        if limit is not None:
            centers = centers[:limit]
        for variant, label in VARIANTS.items():
            metrics = [
                evaluate(centers, read_sequences(Path(record["result"])))
                for record in records[(dataset.name, variant)]
            ]
            metrics_by_key[(dataset.name, variant)] = metrics
            success = [float(item["success"]) for item in metrics]
            reconstruction = [float(item["reconstruction"]) for item in metrics]
            mean_ed = [float(item["mean_ed"]) for item in metrics]
            summary_rows.append({
                "dataset": dataset.name,
                "variant": label,
                "replicates": len(metrics),
                "success_mean_percent": statistics.mean(success),
                "success_sd_pp": sd(success),
                "reconstruction_mean_percent": statistics.mean(reconstruction),
                "reconstruction_sd_pp": sd(reconstruction),
                "mean_ed_mean": statistics.mean(mean_ed),
                "mean_ed_sd": sd(mean_ed),
            })

    pairwise_rows: list[dict[str, object]] = []
    for dataset in datasets:
        centers = read_sequences(dataset.centers)
        if limit is not None:
            centers = centers[:limit]
        local = read_sequences(Path(records[(dataset.name, "local")][0]["result"]))
        correct = read_sequences(Path(records[(dataset.name, "correct")][0]["result"]))
        local_exact = [a == b for a, b in zip(centers, local)]
        correct_exact = [a == b for a, b in zip(centers, correct)]
        both_exact = sum(a and b for a, b in zip(local_exact, correct_exact))
        harmed = sum(a and not b for a, b in zip(local_exact, correct_exact))
        rescued = sum(not a and b for a, b in zip(local_exact, correct_exact))
        both_inexact = sum(not a and not b for a, b in zip(local_exact, correct_exact))
        pairwise_rows.append({
            "dataset": dataset.name,
            "both_exact": both_exact,
            "local_only_exact_correct_harmed": harmed,
            "correct_only_exact_rescued": rescued,
            "both_inexact": both_inexact,
            "net_exact_gain": rescued - harmed,
            "mcnemar_exact_two_sided_p": mcnemar_exact_two_sided(harmed, rescued),
        })

    if args.mode == "smoke":
        summary_path = outdir / "counterfactual_summary_smoke.csv"
        pairwise_path = outdir / "local_vs_correct_pairwise_smoke.csv"
    else:
        summary_path = outdir / "counterfactual_summary.csv"
        pairwise_path = outdir / "local_vs_correct_pairwise.csv"
    write_csv(summary_path, summary_rows)
    write_csv(pairwise_path, pairwise_rows)
    print(f"Summary:  {summary_path}")
    print(f"Pairwise: {pairwise_path}")


if __name__ == "__main__":
    main()
