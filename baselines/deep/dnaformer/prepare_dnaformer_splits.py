#!/usr/bin/env python3

import argparse
import hashlib
import json
import random
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare deterministic DNAFormer train/val/test splits."
    )
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--centers", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--target-len", type=int, required=True)
    parser.add_argument("--separator", default="====")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--min-reads", type=int, default=2)
    parser.add_argument("--train-fraction", type=float, default=0.8)
    parser.add_argument("--val-fraction", type=float, default=0.1)
    parser.add_argument("--dataset-name", default="")
    return parser.parse_args()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_clusters(path, separator):
    clusters = []
    current = []

    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(separator):
                if current:
                    clusters.append(current)
                    current = []
            else:
                current.append(line.upper())

    if current:
        clusters.append(current)

    return clusters


def read_centers(path):
    with Path(path).open(encoding="utf-8") as handle:
        return [
            line.strip().upper()
            for line in handle
            if line.strip()
        ]


def write_split(output_dir, name, indices, clusters, centers):
    index_path = output_dir / f"{name}_indices.txt"
    data_path = output_dir / f"{name}.txt"

    with index_path.open("w", encoding="utf-8") as handle:
        for index in indices:
            handle.write(f"{index}\n")

    with data_path.open("w", encoding="utf-8") as handle:
        for index in indices:
            reads = "|".join(clusters[index])
            handle.write(f"{reads}:{centers[index]}\n")

    return index_path, data_path


def main():
    args = parse_args()

    if args.train_fraction <= 0 or args.val_fraction < 0:
        raise ValueError("Split fractions must be non-negative.")

    if args.train_fraction + args.val_fraction >= 1:
        raise ValueError(
            "train_fraction + val_fraction must be less than 1."
        )

    cluster_path = Path(args.clusters).resolve()
    center_path = Path(args.centers).resolve()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    clusters = read_clusters(cluster_path, args.separator)
    centers = read_centers(center_path)

    if len(clusters) != len(centers):
        raise ValueError(
            f"Cluster/center count mismatch: "
            f"{len(clusters)} != {len(centers)}"
        )

    invalid_centers = [
        index
        for index, center in enumerate(centers)
        if len(center) != args.target_len
    ]
    if invalid_centers:
        raise ValueError(
            f"{len(invalid_centers)} centers do not have target length "
            f"{args.target_len}; first index={invalid_centers[0]}"
        )

    eligible = [
        index
        for index, reads in enumerate(clusters)
        if len(reads) >= args.min_reads
    ]

    random.Random(args.seed).shuffle(eligible)

    train_end = int(len(eligible) * args.train_fraction)
    val_end = train_end + int(len(eligible) * args.val_fraction)

    splits = {
        "train": eligible[:train_end],
        "val": eligible[train_end:val_end],
        "test": eligible[val_end:],
    }

    output_files = {}
    for name, indices in splits.items():
        index_path, data_path = write_split(
            output_dir,
            name,
            indices,
            clusters,
            centers,
        )
        output_files[name] = {
            "count": len(indices),
            "indices_file": index_path.name,
            "data_file": data_path.name,
            "indices_sha256": sha256(index_path),
            "data_sha256": sha256(data_path),
        }

    center_sets = {
        name: {centers[index] for index in indices}
        for name, indices in splits.items()
    }
    duplicate_centers_across_splits = sum(
        len(center_sets[left] & center_sets[right])
        for left, right in (
            ("train", "val"),
            ("train", "test"),
            ("val", "test"),
        )
    )

    manifest = {
        "dataset": args.dataset_name,
        "seed": args.seed,
        "separator": args.separator,
        "target_length": args.target_len,
        "minimum_reads": args.min_reads,
        "train_fraction": args.train_fraction,
        "validation_fraction": args.val_fraction,
        "total_clusters": len(clusters),
        "eligible_clusters": len(eligible),
        "excluded_coverage_lt_2": len(clusters) - len(eligible),
        "duplicate_centers_across_splits":
            duplicate_centers_across_splits,
        "source": {
            "clusters": str(cluster_path),
            "centers": str(center_path),
            "clusters_sha256": sha256(cluster_path),
            "centers_sha256": sha256(center_path),
        },
        "splits": output_files,
    }

    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    print(f"Total clusters: {len(clusters)}")
    print(f"Eligible clusters: {len(eligible)}")
    print(
        "Excluded coverage < "
        f"{args.min_reads}: {len(clusters) - len(eligible)}"
    )
    for name in ("train", "val", "test"):
        print(f"{name}: {len(splits[name])}")
    print(f"Written: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
