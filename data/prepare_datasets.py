#!/usr/bin/env python3
"""Prepare and validate CCP-DP benchmark datasets."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


SEPARATOR = "=" * 31


def read_nonempty_lines(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def read_clusters(path: Path) -> list[list[str]]:
    text = path.read_text(encoding="utf-8")
    raw_parts = text.split(SEPARATOR)
    clusters = [
        [line.strip() for line in part.splitlines() if line.strip()]
        for part in raw_parts
    ]
    while clusters and not clusters[0]:
        clusters.pop(0)
    while clusters and not clusters[-1]:
        clusters.pop()
    return clusters


def write_clusters(path: Path, clusters: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for cluster in clusters:
            handle.write(SEPARATOR + "\n")
            for read in cluster:
                handle.write(read + "\n")


def filter_empty(args: argparse.Namespace) -> None:
    clusters = read_clusters(args.clusters)
    centers = read_nonempty_lines(args.centers)
    if len(clusters) != len(centers):
        raise ValueError(
            f"Cluster/center count mismatch: {len(clusters)} != {len(centers)}"
        )

    kept_clusters: list[list[str]] = []
    kept_centers: list[str] = []
    removed: list[int] = []
    for index, (cluster, center) in enumerate(zip(clusters, centers), start=1):
        if cluster:
            kept_clusters.append(cluster)
            kept_centers.append(center)
        else:
            removed.append(index)

    write_clusters(args.output_clusters, kept_clusters)
    args.output_centers.parent.mkdir(parents=True, exist_ok=True)
    args.output_centers.write_text(
        "".join(center + "\n" for center in kept_centers), encoding="utf-8"
    )
    args.removed_indices.parent.mkdir(parents=True, exist_ok=True)
    args.removed_indices.write_text(
        "".join(f"{index}\n" for index in removed), encoding="utf-8"
    )
    print(f"Original clusters: {len(clusters)}")
    print(f"Removed empty clusters: {len(removed)}")
    print(f"Remaining clusters: {len(kept_clusters)}")


def is_separator(line: str) -> bool:
    return len(line) >= 3 and set(line) == {"="}


def normalize_separators(args: argparse.Namespace) -> None:
    if args.input.resolve() == args.output.resolve():
        raise ValueError("Input and output paths must differ")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cluster_open = False
    boundary_written = False
    clusters = 0
    reads = 0
    with args.input.open("r", encoding="utf-8") as source, args.output.open(
        "w", encoding="utf-8", newline="\n"
    ) as target:
        for raw_line in source:
            line = raw_line.strip()
            if not line:
                continue
            if is_separator(line):
                if cluster_open:
                    target.write(SEPARATOR + "\n")
                    clusters += 1
                    cluster_open = False
                continue
            if not boundary_written:
                target.write(SEPARATOR + "\n")
                boundary_written = True
            target.write(line + "\n")
            reads += 1
            cluster_open = True
        if cluster_open:
            target.write(SEPARATOR + "\n")
            clusters += 1
    print(f"Clusters: {clusters}")
    print(f"Reads: {reads}")
    print(f"Written: {args.output}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(args: argparse.Namespace) -> None:
    config_path = args.config.resolve()
    repo_root = config_path.parent.parent
    config = json.loads(config_path.read_text(encoding="utf-8"))
    selected = args.dataset or list(config["datasets"])
    failed = False

    for dataset_id in selected:
        if dataset_id not in config["datasets"]:
            print(f"[FAIL] unknown dataset: {dataset_id}", file=sys.stderr)
            failed = True
            continue
        item = config["datasets"][dataset_id]
        dataset_missing = False
        for kind in ("clusters", "centers"):
            path = repo_root / item[kind]
            expected = item[f"{kind}_sha256"]
            if not path.is_file():
                print(f"[SKIP] {dataset_id}/{kind}: missing {path}")
                dataset_missing = True
                continue
            actual = sha256(path)
            if actual != expected:
                print(
                    f"[FAIL] {dataset_id}/{kind}: sha256 {actual} != {expected}",
                    file=sys.stderr,
                )
                failed = True
            else:
                print(f"[PASS] {dataset_id}/{kind}")

        centers_path = repo_root / item["centers"]
        if not dataset_missing and centers_path.is_file():
            centers = read_nonempty_lines(centers_path)
            lengths = {len(center) for center in centers}
            if len(centers) != item["num_clusters"]:
                print(
                    f"[FAIL] {dataset_id}: centers={len(centers)}, "
                    f"expected={item['num_clusters']}",
                    file=sys.stderr,
                )
                failed = True
            elif lengths != {item["target_length"]}:
                print(
                    f"[FAIL] {dataset_id}: center lengths={sorted(lengths)}, "
                    f"expected={item['target_length']}",
                    file=sys.stderr,
                )
                failed = True
            else:
                print(
                    f"[PASS] {dataset_id}: {len(centers)} centers, "
                    f"{item['target_length']} bp"
                )

    if failed:
        raise SystemExit(1)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    empty_parser = subparsers.add_parser(
        "filter-empty", help="Remove empty clusters and aligned centers"
    )
    empty_parser.add_argument("--clusters", type=Path, required=True)
    empty_parser.add_argument("--centers", type=Path, required=True)
    empty_parser.add_argument("--output-clusters", type=Path, required=True)
    empty_parser.add_argument("--output-centers", type=Path, required=True)
    empty_parser.add_argument("--removed-indices", type=Path, required=True)
    empty_parser.set_defaults(func=filter_empty)

    separator_parser = subparsers.add_parser(
        "normalize-separators", help="Normalize boundaries to 31 equals signs"
    )
    separator_parser.add_argument("--input", type=Path, required=True)
    separator_parser.add_argument("--output", type=Path, required=True)
    separator_parser.set_defaults(func=normalize_separators)

    validate_parser = subparsers.add_parser(
        "validate", help="Validate present datasets against the manifest"
    )
    validate_parser.add_argument(
        "--config", type=Path, default=Path("data/datasets.json")
    )
    validate_parser.add_argument(
        "--dataset", action="append", help="Dataset ID; repeat to select several"
    )
    validate_parser.set_defaults(func=validate)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
