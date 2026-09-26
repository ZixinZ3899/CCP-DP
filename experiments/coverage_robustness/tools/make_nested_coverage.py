#!/usr/bin/env python3
"""Generate and verify deterministic nested read-coverage inputs.

For every retained cluster and seed, this script creates one deterministic
permutation of the original reads. Coverage-c input is the first c reads of
that permutation, so 5x is a prefix/subset of 10x, 10x of 15x, and 15x of
20x. The same generated Clusters.txt files should be supplied unchanged to
every reconstruction method.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import random
import sys
from pathlib import Path
from typing import Iterable


SEPARATOR = "==============================="
DNA = frozenset("ACGTN")
FORMAT_VERSION = 2


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_separator(line: str) -> bool:
    text = line.strip()
    return len(text) >= 3 and set(text) == {"="}


def dna_tokens(line: str, path: Path, line_number: int) -> list[str]:
    tokens = [token.upper() for token in line.split() if token]
    for token in tokens:
        illegal = sorted(set(token) - DNA)
        if illegal:
            raise ValueError(
                f"{path}:{line_number}: invalid DNA character(s) {illegal}"
            )
    return tokens


def read_clusters(path: Path) -> list[list[str]]:
    """Read separator-delimited clusters; a separator starts a new cluster."""
    clusters: list[list[str]] = []
    current: list[str] | None = None
    with path.open("r", encoding="utf-8-sig", newline=None) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if is_separator(line):
                if current is not None:
                    clusters.append(current)
                current = []
                continue
            if current is None:
                current = []
            current.extend(dna_tokens(line, path, line_number))
    # A final separator is treated as an end marker, not as an extra empty
    # cluster. Consecutive separators inside the file still preserve an empty
    # cluster, so cluster/center alignment problems remain visible.
    if current:
        clusters.append(current)
    if not clusters:
        raise ValueError(f"No clusters found in {path}")
    return clusters


def read_centers(path: Path) -> list[str]:
    """Read one-center-per-line text or FASTA."""
    with path.open("r", encoding="utf-8-sig", newline=None) as handle:
        lines = [line.strip() for line in handle if line.strip()]
    if not lines:
        raise ValueError(f"No centers found in {path}")
    if lines[0].startswith(">"):
        centers: list[str] = []
        sequence: list[str] = []
        for line in lines:
            if line.startswith(">"):
                if sequence:
                    centers.append("".join(sequence).upper())
                    sequence = []
            else:
                sequence.append(line)
        if sequence:
            centers.append("".join(sequence).upper())
    else:
        centers = [line.upper() for line in lines]
    for index, center in enumerate(centers):
        illegal = sorted(set(center) - DNA)
        if illegal:
            raise ValueError(
                f"{path}: center {index} has invalid DNA character(s) {illegal}"
            )
    return centers


def write_clusters(path: Path, clusters: Iterable[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for cluster in clusters:
            handle.write(SEPARATOR + "\n")
            for read in cluster:
                handle.write(read + "\n")


def write_lines(path: Path, values: Iterable[str | int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for value in values:
            handle.write(f"{value}\n")


def stable_cluster_seed(seed: int, original_cluster_index: int) -> int:
    """Derive a platform-stable per-cluster RNG seed."""
    payload = f"nested-coverage-v{FORMAT_VERSION}:{seed}:{original_cluster_index}"
    digest = hashlib.sha256(payload.encode("ascii")).digest()
    return int.from_bytes(digest[:8], byteorder="big", signed=False)


def parse_positive_unique(values: list[int], name: str) -> list[int]:
    if not values or any(value <= 0 for value in values):
        raise ValueError(f"{name} must contain positive integers")
    if len(values) != len(set(values)):
        raise ValueError(f"{name} contains duplicate values")
    return sorted(values)


def sample_command(args: argparse.Namespace) -> None:
    input_path = args.input.resolve()
    output_dir = args.output_dir.resolve()
    centers_path = args.centers.resolve() if args.centers else None
    coverages = parse_positive_unique(args.coverages, "--coverages")
    seeds = parse_positive_unique(args.seeds, "--seeds")
    required_coverage = args.min_original_coverage or max(coverages)
    if required_coverage < max(coverages):
        raise ValueError(
            "--min-original-coverage cannot be smaller than the largest coverage"
        )
    if output_dir.exists() and any(output_dir.iterdir()) and not args.force:
        raise FileExistsError(
            f"Output directory is not empty: {output_dir}\n"
            "Use a new directory or add --force to replace generated files."
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    clusters = read_clusters(input_path)
    centers = read_centers(centers_path) if centers_path else None
    if centers is not None and len(centers) != len(clusters):
        raise ValueError(
            f"Cluster/center count mismatch: {len(clusters)} clusters versus "
            f"{len(centers)} centers"
        )

    eligible_indices = [
        index for index, reads in enumerate(clusters)
        if len(reads) >= required_coverage
    ]
    if not eligible_indices:
        raise ValueError(
            f"No cluster has at least {required_coverage} original reads"
        )

    write_lines(output_dir / "eligible_original_indices.txt", eligible_indices)
    eligible_centers = (
        [centers[index] for index in eligible_indices]
        if centers is not None else None
    )
    if eligible_centers is not None:
        write_lines(output_dir / "Centers_eligible.txt", eligible_centers)

    selection_path = output_dir / "selection_manifest.csv.gz"
    run_rows: list[dict[str, object]] = []
    with gzip.open(selection_path, "wt", encoding="utf-8", newline="") as zipped:
        selection_writer = csv.DictWriter(
            zipped,
            fieldnames=[
                "seed", "eligible_cluster_index", "original_cluster_index",
                "rank", "original_read_index",
            ],
        )
        selection_writer.writeheader()
        for seed in seeds:
            sampled_at_max: list[list[str]] = []
            for eligible_index, original_index in enumerate(eligible_indices):
                reads = clusters[original_index]
                order = list(range(len(reads)))
                random.Random(stable_cluster_seed(seed, original_index)).shuffle(order)
                chosen = order[: max(coverages)]
                sampled_at_max.append([reads[index] for index in chosen])
                for rank, original_read_index in enumerate(chosen, start=1):
                    selection_writer.writerow({
                        "seed": seed,
                        "eligible_cluster_index": eligible_index,
                        "original_cluster_index": original_index,
                        "rank": rank,
                        "original_read_index": original_read_index,
                    })

            for coverage in coverages:
                output_path = (
                    output_dir / f"seed_{seed}" / f"coverage_{coverage}x" /
                    "Clusters.txt"
                )
                write_clusters(
                    output_path,
                    (reads[:coverage] for reads in sampled_at_max),
                )
                row: dict[str, object] = {
                    "seed": seed,
                    "coverage": coverage,
                    "clusters": len(eligible_indices),
                    "reads_per_cluster": coverage,
                    "relative_path": str(output_path.relative_to(output_dir)),
                    "sha256": sha256_file(output_path),
                    "centers_relative_path": "",
                    "centers_sha256": "",
                }
                if eligible_centers is not None:
                    coverage_centers_path = output_path.parent / "Centers.txt"
                    write_lines(coverage_centers_path, eligible_centers)
                    row["centers_relative_path"] = str(
                        coverage_centers_path.relative_to(output_dir)
                    )
                    row["centers_sha256"] = sha256_file(
                        coverage_centers_path
                    )
                run_rows.append(row)

    manifest_path = output_dir / "sampling_manifest.csv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(run_rows[0]))
        writer.writeheader()
        writer.writerows(run_rows)

    metadata = {
        "format_version": FORMAT_VERSION,
        "algorithm": (
            "one deterministic per-cluster permutation per seed; "
            "coverage inputs are prefixes of that permutation"
        ),
        "separator": SEPARATOR,
        "input_clusters": str(input_path),
        "input_clusters_sha256": sha256_file(input_path),
        "input_centers": str(centers_path) if centers_path else None,
        "input_centers_sha256": sha256_file(centers_path) if centers_path else None,
        "coverages": coverages,
        "seeds": seeds,
        "min_original_coverage": required_coverage,
        "original_cluster_count": len(clusters),
        "eligible_cluster_count": len(eligible_indices),
        "excluded_cluster_count": len(clusters) - len(eligible_indices),
        "output_centers": "Centers_eligible.txt" if centers is not None else None,
        "selection_manifest": selection_path.name,
        "sampling_manifest": manifest_path.name,
    }
    metadata_path = output_dir / "sampling_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    print(f"original clusters : {len(clusters)}")
    print(f"eligible clusters : {len(eligible_indices)}")
    print(f"excluded clusters : {len(clusters) - len(eligible_indices)}")
    print(f"coverages         : {', '.join(f'{x}x' for x in coverages)}")
    print(f"seeds             : {', '.join(map(str, seeds))}")
    print(f"output            : {output_dir}")
    verify_output(output_dir, verbose=True)


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def verify_output(output_dir: Path, verbose: bool) -> dict[str, object]:
    output_dir = output_dir.resolve()
    metadata_path = output_dir / "sampling_metadata.json"
    manifest_path = output_dir / "sampling_manifest.csv"
    if not metadata_path.exists() or not manifest_path.exists():
        raise FileNotFoundError(
            "sampling_metadata.json or sampling_manifest.csv is missing"
        )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    coverages = [int(value) for value in metadata["coverages"]]
    seeds = [int(value) for value in metadata["seeds"]]
    expected_clusters = int(metadata["eligible_cluster_count"])
    manifest_rows = load_manifest(manifest_path)
    manifest_by_key = {
        (int(row["seed"]), int(row["coverage"])): row
        for row in manifest_rows
    }

    failures: list[str] = []
    checked_files = 0
    root_centers: list[str] | None = None
    centers_name = metadata.get("output_centers")
    if centers_name:
        root_centers = read_centers(output_dir / centers_name)
        if len(root_centers) != expected_clusters:
            failures.append(
                f"{centers_name}: {len(root_centers)} centers, "
                f"expected {expected_clusters}"
            )
    for seed in seeds:
        sampled: dict[int, list[list[str]]] = {}
        for coverage in coverages:
            key = (seed, coverage)
            row = manifest_by_key.get(key)
            if row is None:
                failures.append(f"missing manifest row for seed={seed}, {coverage}x")
                continue
            path = output_dir / row["relative_path"]
            if not path.exists():
                failures.append(f"missing file: {path}")
                continue
            actual_hash = sha256_file(path)
            if actual_hash != row["sha256"]:
                failures.append(f"SHA-256 mismatch: {path}")
            clusters = read_clusters(path)
            sampled[coverage] = clusters
            checked_files += 1
            if len(clusters) != expected_clusters:
                failures.append(
                    f"{path}: {len(clusters)} clusters, expected {expected_clusters}"
                )
            bad_sizes = [i for i, reads in enumerate(clusters) if len(reads) != coverage]
            if bad_sizes:
                failures.append(
                    f"{path}: {len(bad_sizes)} clusters do not contain exactly "
                    f"{coverage} reads"
                )
            centers_relative_path = row.get("centers_relative_path", "")
            centers_hash = row.get("centers_sha256", "")
            if root_centers is not None:
                if not centers_relative_path:
                    failures.append(
                        f"missing Centers.txt manifest path for {path}"
                    )
                else:
                    coverage_centers_path = output_dir / centers_relative_path
                    if not coverage_centers_path.exists():
                        failures.append(f"missing file: {coverage_centers_path}")
                    else:
                        if sha256_file(coverage_centers_path) != centers_hash:
                            failures.append(
                                f"SHA-256 mismatch: {coverage_centers_path}"
                            )
                        coverage_centers = read_centers(coverage_centers_path)
                        if coverage_centers != root_centers:
                            failures.append(
                                f"{coverage_centers_path}: center rows/order do "
                                "not match Centers_eligible.txt"
                            )

        for low, high in zip(coverages, coverages[1:]):
            if low not in sampled or high not in sampled:
                continue
            if len(sampled[low]) != len(sampled[high]):
                continue
            non_nested = [
                index
                for index, (low_reads, high_reads) in enumerate(
                    zip(sampled[low], sampled[high])
                )
                if low_reads != high_reads[:low]
            ]
            if non_nested:
                failures.append(
                    f"seed={seed}: {len(non_nested)} clusters violate "
                    f"{low}x prefix-of-{high}x nesting"
                )

    report = {
        "status": "PASS" if not failures else "FAIL",
        "checked_files": checked_files,
        "expected_files": len(seeds) * len(coverages),
        "eligible_clusters": expected_clusters,
        "coverages": coverages,
        "seeds": seeds,
        "checks": [
            "file SHA-256 matches manifest",
            "identical cluster cohort and order",
            "exact reads per cluster",
            "lower coverage is an ordered prefix of higher coverage",
            "every coverage directory contains a Centers.txt file",
            "every Centers.txt matches the eligible center rows and order",
        ],
        "failures": failures,
    }
    report_path = output_dir / "verification_report.json"
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if verbose:
        print(f"verification       : {report['status']}")
        print(f"checked files      : {checked_files}/{report['expected_files']}")
        if failures:
            for failure in failures:
                print(f"  FAIL: {failure}")
    if failures:
        raise RuntimeError(f"Nested-sampling verification failed ({len(failures)} errors)")
    return report


def verify_command(args: argparse.Namespace) -> None:
    verify_output(args.output_dir, verbose=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate and verify nested 5x/10x/15x/20x cluster inputs."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    sample = subparsers.add_parser("sample", help="generate nested samples")
    sample.add_argument("-i", "--input", type=Path, required=True,
                        help="original separator-delimited Clusters file")
    sample.add_argument("-c", "--centers", type=Path,
                        help="matching centers, one per line or FASTA")
    sample.add_argument("-o", "--output-dir", type=Path, required=True)
    sample.add_argument(
        "--coverages", type=int, nargs="+", default=[5, 10, 15, 20]
    )
    sample.add_argument("--seeds", type=int, nargs="+",
                        default=[2026, 2027, 2028, 2029, 2030])
    sample.add_argument(
        "--min-original-coverage", type=int,
        help="fixed cohort threshold; default is max(--coverages)",
    )
    sample.add_argument(
        "--force", action="store_true",
        help="allow replacing files inside an existing output directory",
    )
    sample.set_defaults(func=sample_command)

    verify = subparsers.add_parser("verify", help="verify generated outputs")
    verify.add_argument("-o", "--output-dir", type=Path, required=True)
    verify.set_defaults(func=verify_command)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        args.func(args)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
