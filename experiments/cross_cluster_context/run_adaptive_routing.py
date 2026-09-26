#!/usr/bin/env python3

import argparse
import csv
import re
import statistics
import subprocess
import time
from pathlib import Path

try:
    from rapidfuzz.distance import Levenshtein

    def edit_distance(a, b):
        return Levenshtein.distance(a, b)

except ImportError:
    def edit_distance(a, b):
        """Dependency-free fallback; rapidfuzz is strongly recommended."""
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


def read_sequences(path):
    path = Path(path)
    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    if not lines:
        raise RuntimeError(f"No sequences found in {path}")

    # Supports FASTA and one-sequence-per-line files.
    if any(line.startswith(">") for line in lines):
        sequences = []
        current = []
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


def evaluate(centers, predictions):
    if len(centers) != len(predictions):
        raise RuntimeError(
            f"Sequence-count mismatch: centers={len(centers)}, "
            f"predictions={len(predictions)}"
        )

    exact = 0
    total_ed = 0
    total_bases = 0

    for center, prediction in zip(centers, predictions):
        if center == prediction:
            exact += 1
        total_ed += edit_distance(center, prediction)
        total_bases += len(center)

    n = len(centers)

    return {
        "exact_count": exact,
        "success_pct": 100.0 * exact / n,
        "mean_ed": total_ed / n,
        "reconstruction_pct": (
            100.0 * (1.0 - total_ed / total_bases)
            if total_bases else float("nan")
        ),
    }


def extract_number(text, pattern, default=float("nan")):
    match = re.search(pattern, text, flags=re.MULTILINE)
    return float(match.group(1)) if match else default


def parse_log(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")

    cluster_match = re.search(r"^clusters:\s*(\d+)", text, re.MULTILINE)
    decoded_match = re.search(
        r"^graph-decoded clusters:\s*(\d+)/(\d+)",
        text,
        re.MULTILINE,
    )
    changed_match = re.search(
        r"^unified-graph changes:\s*(\d+)/(\d+)",
        text,
        re.MULTILINE,
    )

    return {
        "clusters": int(cluster_match.group(1)) if cluster_match else None,
        "structured_count": (
            int(decoded_match.group(1)) if decoded_match else 0
        ),
        "changed_count": (
            int(changed_match.group(1)) if changed_match else 0
        ),
        "phase_load": extract_number(
            text,
            r"^preliminary phase load:\s*([0-9.eE+-]+)",
        ),
        "total_error": extract_number(
            text,
            r"^preliminary total error:\s*([0-9.eE+-]+)",
        ),
        "cross_reuse": extract_number(
            text,
            r"^cross-cluster reuse \(distinct guides\):\s*"
            r"([0-9.eE+-]+)",
        ),
        "internal_elapsed": extract_number(
            text,
            r"^elapsed seconds:\s*([0-9.eE+-]+)",
        ),
    }


def timing_stats(values):
    return {
        "median": statistics.median(values),
        "mean": statistics.mean(values),
        "sd": statistics.stdev(values) if len(values) > 1 else 0.0,
    }


def run_once(
    binary,
    clusters,
    length,
    separator,
    jobs,
    mode,
    repetition,
    outdir,
    guide_path=None,
):
    output_path = outdir / f"{mode}_rep{repetition}.txt"
    diag_path = outdir / f"{mode}_rep{repetition}.diag.csv"
    log_path = outdir / f"{mode}_rep{repetition}.log"

    command = [
        str(binary),
        "-i", str(clusters),
        "-l", str(length),
        "-s", separator,
        "-o", str(output_path),
        "--diag", str(diag_path),
        "--mode", "auto",
        "--jobs", str(jobs),
    ]

    if mode == "decode_all":
        command.append("--decode-all")

    if guide_path is not None:
        command.extend(["--guide-output", str(guide_path)])

    start = time.perf_counter()

    with log_path.open("w", encoding="utf-8") as log_file:
        completed = subprocess.run(
            command,
            stdout=subprocess.DEVNULL,
            stderr=log_file,
            check=False,
        )

    wall_seconds = time.perf_counter() - start

    if completed.returncode != 0:
        raise RuntimeError(
            f"CCP-DP failed with return code {completed.returncode}. "
            f"See {log_path}"
        )

    return {
        "result": output_path,
        "diag": diag_path,
        "log": log_path,
        "wall_seconds": wall_seconds,
    }


def update_summary(path, row):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(row.keys())
    previous = []

    if path.exists():
        with path.open(newline="", encoding="utf-8") as handle:
            previous = list(csv.DictReader(handle))

        # Replace an older result for the same dataset.
        previous = [
            old for old in previous
            if old.get("dataset") != row["dataset"]
        ]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(previous)
        writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare adaptive selective routing with decode-all "
            "using CCP-DP."
        )
    )
    parser.add_argument("--binary", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--centers", required=True)
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--separator", default="====")
    parser.add_argument("--jobs", type=int, default=32)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--summary-csv",
        default="adaptive_routing_summary.csv",
    )
    args = parser.parse_args()

    binary = Path(args.binary).resolve()
    clusters = Path(args.clusters).resolve()
    centers_path = Path(args.centers).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not binary.exists():
        raise RuntimeError(f"Binary not found: {binary}")
    if not clusters.exists():
        raise RuntimeError(f"Cluster file not found: {clusters}")
    if not centers_path.exists():
        raise RuntimeError(f"Center file not found: {centers_path}")

    guide_path = outdir / "frozen_guides.txt"
    records = {"selective": [], "decode_all": []}

    # Alternate order to reduce systematic cache/order effects.
    for repetition in range(1, args.repeats + 1):
        order = (
            ["selective", "decode_all"]
            if repetition % 2 == 1
            else ["decode_all", "selective"]
        )

        for mode in order:
            current_guide_path = (
                guide_path
                if mode == "selective" and repetition == 1
                else None
            )

            print(
                f"[{args.dataset}] {mode}, "
                f"repetition {repetition}/{args.repeats}"
            )

            record = run_once(
                binary=binary,
                clusters=clusters,
                length=args.length,
                separator=args.separator,
                jobs=args.jobs,
                mode=mode,
                repetition=repetition,
                outdir=outdir,
                guide_path=current_guide_path,
            )
            records[mode].append(record)

    centers = read_sequences(centers_path)
    guides = read_sequences(guide_path)

    selective_result = read_sequences(
        records["selective"][0]["result"]
    )
    decode_all_result = read_sequences(
        records["decode_all"][0]["result"]
    )

    guide_metrics = evaluate(centers, guides)
    selective_metrics = evaluate(centers, selective_result)
    decode_all_metrics = evaluate(centers, decode_all_result)

    selective_log = parse_log(records["selective"][0]["log"])
    decode_all_log = parse_log(records["decode_all"][0]["log"])

    n = len(centers)
    structured_count = selective_log["structured_count"]
    keep_count = n - structured_count

    selective_times = [
        record["wall_seconds"] for record in records["selective"]
    ]
    decode_all_times = [
        record["wall_seconds"] for record in records["decode_all"]
    ]

    selective_time = timing_stats(selective_times)
    decode_all_time = timing_stats(decode_all_times)

    speedup = (
        decode_all_time["median"] / selective_time["median"]
        if selective_time["median"] > 0 else float("nan")
    )
    time_reduction = (
        100.0 * (
            1.0
            - selective_time["median"] / decode_all_time["median"]
        )
        if decode_all_time["median"] > 0 else float("nan")
    )

    row = {
        "dataset": args.dataset,
        "length_bp": args.length,
        "clusters": n,
        "guide_success_pct": f"{guide_metrics['success_pct']:.4f}",
        "guide_failure_pct": f"{100-guide_metrics['success_pct']:.4f}",
        "full_success_pct": f"{selective_metrics['success_pct']:.4f}",
        "decode_all_success_pct": (
            f"{decode_all_metrics['success_pct']:.4f}"
        ),
        "delta_success_pp": (
            f"{selective_metrics['success_pct'] - guide_metrics['success_pct']:.4f}"
        ),
        "selective_minus_decode_all_pp": (
            f"{selective_metrics['success_pct'] - decode_all_metrics['success_pct']:.4f}"
        ),
        "guide_mean_ed": f"{guide_metrics['mean_ed']:.6f}",
        "full_mean_ed": f"{selective_metrics['mean_ed']:.6f}",
        "decode_all_mean_ed": f"{decode_all_metrics['mean_ed']:.6f}",
        "keep_count": keep_count,
        "keep_pct": f"{100.0 * keep_count / n:.4f}",
        "structured_count": structured_count,
        "structured_pct": f"{100.0 * structured_count / n:.4f}",
        "changed_count": selective_log["changed_count"],
        "cross_reuse": f"{selective_log['cross_reuse']:.6f}",
        "phase_load": f"{selective_log['phase_load']:.6f}",
        "selective_time_median_s": (
            f"{selective_time['median']:.6f}"
        ),
        "selective_time_sd_s": f"{selective_time['sd']:.6f}",
        "decode_all_time_median_s": (
            f"{decode_all_time['median']:.6f}"
        ),
        "decode_all_time_sd_s": f"{decode_all_time['sd']:.6f}",
        "speedup_x": f"{speedup:.4f}",
        "time_reduction_pct": f"{time_reduction:.4f}",
    }

    update_summary(args.summary_csv, row)

    print("\nCompleted")
    print(f"Dataset:                 {args.dataset}")
    print(f"Clusters:                {n}")
    print(
        f"Guide success:           "
        f"{guide_metrics['success_pct']:.2f}%"
    )
    print(
        f"Full success:            "
        f"{selective_metrics['success_pct']:.2f}%"
    )
    print(
        f"Decode-all success:      "
        f"{decode_all_metrics['success_pct']:.2f}%"
    )
    print(
        f"Structured decoding:     "
        f"{structured_count}/{n} "
        f"({100.0 * structured_count / n:.2f}%)"
    )
    print(
        f"Selective runtime:       "
        f"{selective_time['median']:.3f} s "
        f"(median)"
    )
    print(
        f"Decode-all runtime:      "
        f"{decode_all_time['median']:.3f} s "
        f"(median)"
    )
    print(f"Speedup:                 {speedup:.2f}x")
    print(f"Runtime reduction:       {time_reduction:.1f}%")
    print(f"Summary:                 {args.summary_csv}")


if __name__ == "__main__":
    main()