#!/usr/bin/env python3
"""Run the six prespecified CCP-DP component ablations."""

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import subprocess
import time
from pathlib import Path

try:
    from rapidfuzz.distance import Levenshtein

    def edit_distance(left, right):
        return Levenshtein.distance(left, right)
except ImportError:
    def edit_distance(left, right):
        if len(left) < len(right):
            left, right = right, left
        previous = list(range(len(right) + 1))
        for i, ca in enumerate(left, start=1):
            current = [i]
            for j, cb in enumerate(right, start=1):
                current.append(min(
                    previous[j] + 1,
                    current[j - 1] + 1,
                    previous[j - 1] + (ca != cb),
                ))
            previous = current
        return previous[-1]


VARIANTS = [
    {
        "id": "guide_only",
        "label": "Guide only",
        "flags": ["--disable-unified-graph"],
        "joint_trellis": 0,
        "cross_loo_score": 0,
        "blind_ids_likelihood": 0,
        "sparse_gate": 0,
    },
    {
        "id": "without_cross_loo",
        "label": "Full w/o Cross-LOO score",
        "flags": ["--disable-cross-loo-score"],
        "joint_trellis": 1,
        "cross_loo_score": 0,
        "blind_ids_likelihood": 1,
        "sparse_gate": 1,
    },
    {
        "id": "without_blind_ids",
        "label": "Full w/o blind IDS likelihood",
        "flags": ["--disable-blind-ids-likelihood"],
        "joint_trellis": 1,
        "cross_loo_score": 1,
        "blind_ids_likelihood": 0,
        "sparse_gate": 1,
    },
    {
        "id": "decode_all",
        "label": "Full with decoding of all clusters",
        "flags": ["--decode-all"],
        "joint_trellis": 1,
        "cross_loo_score": 1,
        "blind_ids_likelihood": 1,
        "sparse_gate": 0,
    },
    {
        "id": "local_trellis_only",
        "label": "Local trellis only",
        "flags": [
            "--decode-all",
            "--disable-cross-loo-score",
            "--disable-blind-ids-likelihood",
        ],
        "joint_trellis": 1,
        "cross_loo_score": 0,
        "blind_ids_likelihood": 0,
        "sparse_gate": 0,
    },
    {
        "id": "full",
        "label": "Full",
        "flags": [],
        "joint_trellis": 1,
        "cross_loo_score": 1,
        "blind_ids_likelihood": 1,
        "sparse_gate": 1,
    },
]


def resolve_dataset_paths(item):
    clusters = item["clusters"] if isinstance(item["clusters"], list) else [item["clusters"]]
    centers = item["centers"] if isinstance(item["centers"], list) else [item["centers"]]
    if len(clusters) != len(centers):
        raise ValueError(f"Candidate-pair count differs for {item['name']}")
    repo_root = Path(__file__).resolve().parents[2]
    for cluster_candidate, center_candidate in zip(clusters, centers):
        cluster_path = Path(
            str(cluster_candidate).format(repo_root=repo_root)
        ).expanduser()
        center_path = Path(
            str(center_candidate).format(repo_root=repo_root)
        ).expanduser()
        if cluster_path.exists() and center_path.exists():
            return cluster_path.resolve(), center_path.resolve()
    listed = "\n".join(
        f"  - clusters={cluster}; centers={center}"
        for cluster, center in zip(clusters, centers)
    )
    raise FileNotFoundError(
        f"No complete cluster/center pair found for {item['name']}. Tried:\n{listed}"
    )


def read_sequences(path, limit=None):
    lines = [
        line.strip()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if any(line.startswith(">") for line in lines):
        sequences, current = [], []
        for line in lines:
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                    current = []
            else:
                current.append(line)
        if current:
            sequences.append("".join(current).upper())
    else:
        sequences = [line.upper() for line in lines]
    return sequences if limit is None else sequences[:limit]


def evaluate(centers, predictions):
    if len(centers) != len(predictions):
        raise RuntimeError(
            f"Sequence-count mismatch: centers={len(centers)}, "
            f"predictions={len(predictions)}"
        )
    distances = [edit_distance(a, b) for a, b in zip(centers, predictions)]
    exact = sum(a == b for a, b in zip(centers, predictions))
    total_bases = sum(len(seq) for seq in centers)
    total_ed = sum(distances)
    return {
        "sequences": len(centers),
        "exact_count": exact,
        "success_percent": 100.0 * exact / len(centers),
        "mean_ed": total_ed / len(centers),
        "reconstruction_percent": 100.0 * (1.0 - total_ed / total_bases),
    }


def evaluate_external(evaluator, centers_path, predictions_path, limit):
    command = [str(evaluator), str(centers_path), str(predictions_path)]
    if limit:
        command.append(str(limit))
    completed = subprocess.run(
        command, text=True, capture_output=True, check=True
    )
    parsed = {}
    for line in completed.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            parsed[key] = value
    return {
        "sequences": int(parsed["sequences"]),
        "exact_count": int(parsed["exact_count"]),
        "success_percent": float(parsed["success_percent"]),
        "mean_ed": float(parsed["mean_ed"]),
        "reconstruction_percent": float(parsed["reconstruction_percent"]),
    }


def extract_float(text, pattern):
    match = re.search(pattern, text, re.MULTILINE)
    return float(match.group(1)) if match else math.nan


def parse_log(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    decoded = re.search(r"^graph-decoded clusters:\s*(\d+)/(\d+)", text, re.M)
    changed = re.search(r"^unified-graph changes:\s*(\d+)/(\d+)", text, re.M)
    planned = re.search(r"^planned structured decoding:\s*(\d+)/(\d+)", text, re.M)
    weights = re.search(
        r"^trellis weights \(cross,local\):\s*([0-9.eE+-]+),([0-9.eE+-]+)",
        text,
        re.M,
    )
    policy = re.search(r"^routing policy:\s*(.+)$", text, re.M)
    return {
        "decoded_count": int(decoded.group(1)) if decoded else 0,
        "decoded_total": int(decoded.group(2)) if decoded else 0,
        "changed_count": int(changed.group(1)) if changed else 0,
        "planned_count": int(planned.group(1)) if planned else 0,
        "routing_policy": policy.group(1).strip() if policy else "unknown",
        "cross_reuse": extract_float(
            text, r"^cross-cluster reuse \(distinct guides\):\s*([0-9.eE+-]+)"
        ),
        "phase_load": extract_float(
            text, r"^preliminary phase load:\s*([0-9.eE+-]+)"
        ),
        "cross_weight": float(weights.group(1)) if weights else math.nan,
        "local_weight": float(weights.group(2)) if weights else math.nan,
        "internal_seconds": extract_float(
            text, r"^elapsed seconds:\s*([0-9.eE+-]+)"
        ),
    }


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path, rows, fieldnames=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    names = fieldnames or list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def run_one(binary, dataset, variant, repetition, outdir, jobs, limit):
    stem = f"{variant['id']}_rep{repetition}"
    result = outdir / f"{stem}.txt"
    diag = outdir / f"{stem}.diag.csv"
    log = outdir / f"{stem}.log"
    command = [
        str(binary),
        "-i", str(dataset["clusters_path"]),
        "-l", str(dataset["length"]),
        "-s", dataset["separator"],
        "-o", str(result),
        "--diag", str(diag),
        "--mode", "auto",
        "--jobs", str(jobs),
    ]
    if limit:
        command.extend(["--limit", str(limit)])
    command.extend(variant["flags"])
    start = time.perf_counter()
    with log.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(
            command,
            stdout=subprocess.DEVNULL,
            stderr=handle,
            check=False,
        )
    wall = time.perf_counter() - start
    if completed.returncode != 0:
        raise RuntimeError(
            f"Run failed ({dataset['name']} / {variant['label']}). See {log}"
        )
    parsed = parse_log(log)
    return {
        "result_path": result,
        "diag_path": diag,
        "log_path": log,
        "wall_seconds": wall,
        "sha256": file_sha256(result),
        **parsed,
    }


def render_markdown(summary_rows, path):
    headers = [
        "Dataset", "Variant", "Joint", "Cross-LOO", "Blind IDS",
        "Sparse gate", "Success (%)", "Mean ED", "Time-32 (s)",
        "Decoded (%)",
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    mark = {1: "✓", 0: "–"}
    for row in summary_rows:
        values = [
            row["dataset"],
            row["variant_label"],
            mark[int(row["joint_trellis"])],
            mark[int(row["cross_loo_score"])],
            mark[int(row["blind_ids_likelihood"])],
            mark[int(row["sparse_gate"])],
            f"{float(row['success_percent']):.2f}",
            f"{float(row['mean_ed']):.4f}",
            f"{float(row['wall_seconds_mean']):.3f} ± {float(row['wall_seconds_sd']):.3f}",
            f"{float(row['decoded_fraction_percent']):.2f}",
        ]
        lines.append("| " + " | ".join(values) + " |")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--binary", required=True)
    parser.add_argument("--evaluator")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--jobs", type=int, default=32)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    binary = Path(args.binary).resolve()
    if not binary.exists():
        raise FileNotFoundError(f"Binary not found: {binary}")
    evaluator = Path(args.evaluator).resolve() if args.evaluator else None
    if evaluator is not None and not evaluator.exists():
        raise FileNotFoundError(f"Evaluator not found: {evaluator}")
    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    datasets = []
    for item in config["datasets"]:
        current = dict(item)
        current["clusters_path"], current["centers_path"] = (
            resolve_dataset_paths(item)
        )
        datasets.append(current)

    raw_rows, summary_rows = [], []
    total = len(datasets) * len(VARIANTS) * args.repeats
    completed_count = 0
    for dataset in datasets:
        dataset_dir = outdir / dataset["name"]
        dataset_dir.mkdir(parents=True, exist_ok=True)
        centers = read_sequences(dataset["centers_path"], args.limit)
        records = {variant["id"]: [] for variant in VARIANTS}

        for repetition in range(1, args.repeats + 1):
            offset = (repetition - 1) % len(VARIANTS)
            ordered = VARIANTS[offset:] + VARIANTS[:offset]
            for variant in ordered:
                completed_count += 1
                print(
                    f"[{completed_count:03d}/{total:03d}] {dataset['name']} / "
                    f"{variant['label']} / repeat {repetition}",
                    flush=True,
                )
                record = run_one(
                    binary, dataset, variant, repetition, dataset_dir,
                    args.jobs, args.limit,
                )
                records[variant["id"]].append(record)

        for variant in VARIANTS:
            variant_records = records[variant["id"]]
            hashes = {record["sha256"] for record in variant_records}
            if len(hashes) != 1:
                raise RuntimeError(
                    f"Non-deterministic output for {dataset['name']} / "
                    f"{variant['label']}; inspect repetition outputs."
                )
            if evaluator is not None:
                metrics = evaluate_external(
                    evaluator, dataset["centers_path"],
                    variant_records[0]["result_path"], args.limit,
                )
            else:
                predictions = read_sequences(variant_records[0]["result_path"])
                metrics = evaluate(centers, predictions)
            walls = [record["wall_seconds"] for record in variant_records]
            first = variant_records[0]
            decoded_total = first["decoded_total"] or metrics["sequences"]
            decoded_fraction = 100.0 * first["decoded_count"] / decoded_total
            summary = {
                "dataset": dataset["name"],
                "length_bp": dataset["length"],
                "variant_id": variant["id"],
                "variant_label": variant["label"],
                "joint_trellis": variant["joint_trellis"],
                "cross_loo_score": variant["cross_loo_score"],
                "blind_ids_likelihood": variant["blind_ids_likelihood"],
                "sparse_gate": variant["sparse_gate"],
                "repeats": args.repeats,
                "sequences": metrics["sequences"],
                "exact_count": metrics["exact_count"],
                "success_percent": metrics["success_percent"],
                "mean_ed": metrics["mean_ed"],
                "reconstruction_percent": metrics["reconstruction_percent"],
                "wall_seconds_mean": statistics.mean(walls),
                "wall_seconds_sd": statistics.stdev(walls) if len(walls) > 1 else 0.0,
                "wall_seconds_median": statistics.median(walls),
                "decoded_fraction_percent": decoded_fraction,
                "changed_clusters": first["changed_count"],
                "routing_policy": first["routing_policy"],
                "cross_reuse": first["cross_reuse"],
                "phase_load": first["phase_load"],
                "effective_cross_weight": first["cross_weight"],
                "effective_local_weight": first["local_weight"],
            }
            summary_rows.append(summary)
            for repetition, record in enumerate(variant_records, start=1):
                raw_rows.append({
                    "dataset": dataset["name"],
                    "variant_id": variant["id"],
                    "repetition": repetition,
                    "wall_seconds": record["wall_seconds"],
                    "internal_seconds": record["internal_seconds"],
                    "decoded_count": record["decoded_count"],
                    "changed_count": record["changed_count"],
                    "sha256": record["sha256"],
                    "result_path": record["result_path"],
                    "log_path": record["log_path"],
                })

    full_by_dataset = {
        row["dataset"]: row for row in summary_rows if row["variant_id"] == "full"
    }
    for row in summary_rows:
        full = full_by_dataset[row["dataset"]]
        row["delta_success_vs_full_pp"] = (
            row["success_percent"] - full["success_percent"]
        )
        row["delta_mean_ed_vs_full"] = row["mean_ed"] - full["mean_ed"]
        row["runtime_ratio_vs_full"] = (
            row["wall_seconds_mean"] / full["wall_seconds_mean"]
        )

    write_csv(outdir / "raw_runs.csv", raw_rows)
    write_csv(outdir / "ablation_summary.csv", summary_rows)
    write_csv(
        outdir / "ablation_components.csv",
        [{key: variant[key] for key in (
            "id", "label", "joint_trellis", "cross_loo_score",
            "blind_ids_likelihood", "sparse_gate"
        )} for variant in VARIANTS],
    )
    render_markdown(summary_rows, outdir / "ablation_table.md")
    print(f"Summary: {outdir / 'ablation_summary.csv'}")
    print(f"Table:   {outdir / 'ablation_table.md'}")


if __name__ == "__main__":
    main()
