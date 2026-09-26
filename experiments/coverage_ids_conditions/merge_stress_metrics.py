#!/usr/bin/env python3
"""Merge the completed 5/10/15/20x run with the new dense-coverage run."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--refined", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    base_fields, base_rows = read_rows(args.base)
    refined_fields, refined_rows = read_rows(args.refined)
    if base_fields != refined_fields:
        raise RuntimeError(
            "CSV schemas differ. Base columns:\n"
            f"{base_fields}\nRefined columns:\n{refined_fields}"
        )

    # Refined rows intentionally win if a condition was accidentally repeated.
    merged: dict[tuple[str, str], dict[str, str]] = {}
    for row in [*base_rows, *refined_rows]:
        merged[(row["condition_id"], row["method"])] = row

    method_order = {"CCP-DP": 0, "BBS": 1, "CPL": 2}
    rows = sorted(
        merged.values(),
        key=lambda r: (
            int(r["coverage"]), float(r["total_indel_pct"]), int(r["seed"]),
            method_order.get(r["method"], 99),
        ),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=base_fields)
        writer.writeheader()
        writer.writerows(rows)

    expected_coverages = {5, 6, 7, 8, 9, 10, 12, 15, 20}
    found_coverages = {int(row["coverage"]) for row in rows}
    missing = sorted(expected_coverages - found_coverages)
    print(f"Merged {len(rows)} method rows into {args.output.resolve()}")
    print(f"Coverage values: {sorted(found_coverages)}")
    if missing:
        print(f"WARNING: missing coverage values: {missing}")


if __name__ == "__main__":
    main()
