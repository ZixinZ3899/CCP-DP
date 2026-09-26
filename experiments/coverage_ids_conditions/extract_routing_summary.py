#!/usr/bin/env python3
"""Extract CCP-DP structured-decoding counts from existing result logs."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


CONDITION_RE = re.compile(r"cov(\d+)_indel(\d+)_seed(\d+)$")
DECODE_RE = re.compile(r"graph-decoded clusters:\s*(\d+)\s*/\s*(\d+)")
CHANGE_RE = re.compile(r"unified-graph changes:\s*(\d+)\s*/\s*(\d+)")


def mean_sd(values: list[float]) -> tuple[float, float]:
    a = np.asarray(values, dtype=float)
    return float(a.mean()), float(a.std(ddof=1)) if len(a) > 1 else 0.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", action="append", required=True, type=Path)
    parser.add_argument("--raw-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    args = parser.parse_args()

    rows: dict[tuple[int, float, int], dict[str, float | int | str]] = {}
    for root in args.results_root:
        if not root.exists():
            raise FileNotFoundError(root)
        for log_path in sorted(root.glob("cov*_indel*_seed*/ccpdp.log")):
            match = CONDITION_RE.fullmatch(log_path.parent.name)
            if not match:
                continue
            coverage, indel, seed = map(int, match.groups())
            text = log_path.read_text(errors="replace")
            decoded = DECODE_RE.search(text)
            changed = CHANGE_RE.search(text)
            if not decoded:
                raise RuntimeError(f"Missing graph-decoded line in {log_path}")
            decoded_count, clusters = map(int, decoded.groups())
            changed_count = int(changed.group(1)) if changed else -1
            rows[(coverage, float(indel), seed)] = {
                "condition_id": log_path.parent.name,
                "coverage": coverage,
                "total_indel_pct": float(indel),
                "seed": seed,
                "clusters": clusters,
                "decoded_clusters": decoded_count,
                "decoded_fraction_pct": 100.0 * decoded_count / clusters,
                "changed_clusters": changed_count,
                "log_path": str(log_path.resolve()),
            }

    if not rows:
        raise RuntimeError("No CCP-DP logs containing routing statistics were found")

    raw_rows = [rows[key] for key in sorted(rows)]
    args.raw_output.parent.mkdir(parents=True, exist_ok=True)
    with args.raw_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(raw_rows[0]))
        writer.writeheader()
        writer.writerows(raw_rows)

    groups: dict[tuple[int, float], list[dict[str, float | int | str]]] = defaultdict(list)
    for row in raw_rows:
        groups[(int(row["coverage"]), float(row["total_indel_pct"]))].append(row)

    summary = []
    for (coverage, indel), group in sorted(groups.items()):
        route_mean, route_sd = mean_sd([float(x["decoded_fraction_pct"]) for x in group])
        valid_changes = [float(x["changed_clusters"]) for x in group if int(x["changed_clusters"]) >= 0]
        change_mean, change_sd = mean_sd(valid_changes) if valid_changes else (float("nan"), float("nan"))
        summary.append({
            "coverage": coverage,
            "total_indel_pct": indel,
            "decoded_fraction_pct_mean": route_mean,
            "decoded_fraction_pct_sd": route_sd,
            "changed_clusters_mean": change_mean,
            "changed_clusters_sd": change_sd,
            "seeds": len(group),
        })

    with args.summary_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print(f"Routing raw rows: {len(raw_rows)} -> {args.raw_output.resolve()}")
    print(f"Routing conditions: {len(summary)} -> {args.summary_output.resolve()}")


if __name__ == "__main__":
    main()
