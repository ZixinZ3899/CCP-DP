#!/usr/bin/env python3
"""Build the CCP-DP mechanism-analysis table from frozen run outputs.

Ground-truth centers are used only for post-hoc evaluation. They are never
passed to the CCP-DP decoder.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

try:
    import Levenshtein as _levenshtein

    def edit_distance(left: str, right: str) -> int:
        return _levenshtein.distance(left, right)

except ImportError:
    _levenshtein = None

    def edit_distance(left: str, right: str) -> int:
        if len(left) < len(right):
            left, right = right, left
        previous = list(range(len(right) + 1))
        for row, left_base in enumerate(left, start=1):
            current = [row]
            for column, right_base in enumerate(right, start=1):
                current.append(
                    min(
                        current[-1] + 1,
                        previous[column] + 1,
                        previous[column - 1] + (left_base != right_base),
                    )
                )
            previous = current
        return previous[-1]


FIELDS = [
    "record_type",
    "dataset_short",
    "order",
    "decile",
    "clusters",
    "guide_failure_percent",
    "final_failure_percent",
    "structured_decoding_percent",
    "auc",
    "activated_percent",
    "error_type",
    "n",
    "share_percent",
    "exact_rescue_percent",
    "partial_improvement_percent",
    "unchanged_percent",
    "worsened_percent",
    "any_improvement_percent",
    "ed_reduction_percent",
]

ERROR_ORDER = ["Substitution-only", "Paired indel", "Mixed"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate the CCP-DP difficult-cluster mechanism table."
    )
    parser.add_argument(
        "--runs-dir",
        required=True,
        type=Path,
        help="Directory containing the diagnostic, guide, final, and log files",
    )
    parser.add_argument("--srinivas-centers", required=True, type=Path)
    parser.add_argument("--chandak-centers", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bins", type=int, default=10)
    parser.add_argument(
        "--compare",
        type=Path,
        help="Optional existing CSV to compare with the generated table",
    )
    parser.add_argument("--tolerance", type=float, default=1e-9)
    return parser.parse_args()


def read_sequences(path: Path):
    sequences = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            sequence = "".join(raw.strip().split()).upper()
            if not sequence or set(sequence) == {"="} or sequence.startswith("CLUSTER"):
                continue
            if not re.fullmatch(r"[ACGTN]+", sequence):
                raise ValueError(
                    f"Unexpected non-sequence line in {path}:{line_number}: {raw.rstrip()}"
                )
            sequences.append(sequence)
    return sequences


def read_diag(path: Path):
    required = {"cluster_id", "joint_anomaly", "joint_active"}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
        rows = list(reader)
    rows.sort(key=lambda row: int(row["cluster_id"]))
    observed = [int(row["cluster_id"]) for row in rows]
    if observed != list(range(1, len(rows) + 1)):
        raise ValueError(f"cluster_id must be contiguous and 1-based in {path}")
    return rows


def safe_percent(numerator, denominator):
    return 100.0 * numerator / denominator if denominator else math.nan


def binary_auc(scores, labels):
    """Tie-aware Mann-Whitney AUC; label 1 denotes a failed guide."""
    n = len(scores)
    positives = sum(labels)
    negatives = n - positives
    if positives == 0 or negatives == 0:
        return math.nan
    ordered = sorted(range(n), key=lambda index: scores[index])
    rank_sum_positive = 0.0
    start = 0
    while start < n:
        end = start + 1
        while end < n and scores[ordered[end]] == scores[ordered[start]]:
            end += 1
        average_rank = ((start + 1) + end) / 2.0
        rank_sum_positive += average_rank * sum(
            labels[ordered[index]] for index in range(start, end)
        )
        start = end
    return (
        rank_sum_positive - positives * (positives + 1) / 2.0
    ) / (positives * negatives)


def canonical_minimum_alignment(source: str, target: str):
    """Return ED, substitutions, insertions, and deletions deterministically."""
    target_length = len(target)
    previous = [(column, 0, column, 0) for column in range(target_length + 1)]
    for row, source_base in enumerate(source, 1):
        current = [(row, 0, 0, row)]
        for column, target_base in enumerate(target, 1):
            diagonal = previous[column - 1]
            if source_base == target_base:
                candidates = [diagonal]
            else:
                candidates = [
                    (diagonal[0] + 1, diagonal[1] + 1, diagonal[2], diagonal[3])
                ]
            left = current[column - 1]
            candidates.append((left[0] + 1, left[1], left[2] + 1, left[3]))
            up = previous[column]
            candidates.append((up[0] + 1, up[1], up[2], up[3] + 1))
            current.append(
                min(candidates, key=lambda item: (item[0], item[1], item[2] + item[3]))
            )
        previous = current
    return previous[-1]


def classify_guide_error(guide: str, center: str, expected_ed: int):
    if len(guide) != len(center):
        raise ValueError("Guide and center must both have the fixed target length")
    if _levenshtein is None:
        raise RuntimeError(
            "python-Levenshtein is required for the final error-type "
            "classification (pip install python-Levenshtein)"
        )
    edit_operations = _levenshtein.editops(guide, center)
    ed = len(edit_operations)
    if ed != expected_ed:
        raise ValueError(f"ED mismatch: expected={expected_ed}, recomputed={ed}")
    operation_names = [operation[0] for operation in edit_operations]
    has_substitution = "replace" in operation_names
    has_indel = "insert" in operation_names or "delete" in operation_names
    if not has_indel:
        return "Substitution-only"
    if not has_substitution:
        return "Paired indel"
    return "Mixed"


def outcome(guide_ed: int, final_ed: int):
    if final_ed == 0:
        return "exact"
    if final_ed < guide_ed:
        return "partial"
    if final_ed == guide_ed:
        return "unchanged"
    return "worsened"


def analyze_dataset(spec, runs_dir: Path, bins: int):
    diag = read_diag(runs_dir / f"{spec['prefix']}_diag.csv")
    guides = read_sequences(runs_dir / f"{spec['prefix']}_guides.txt")
    finals = read_sequences(runs_dir / f"{spec['prefix']}_final.txt")
    centers = read_sequences(spec["centers"])
    counts = {
        "diag": len(diag),
        "guides": len(guides),
        "finals": len(finals),
        "centers": len(centers),
    }
    if len(set(counts.values())) != 1:
        raise ValueError(f"Sequence/diagnostic count mismatch for {spec['key']}: {counts}")

    records = []
    for cluster_id, (row, guide, final, center) in enumerate(
        zip(diag, guides, finals, centers), 1
    ):
        if len(guide) != spec["length"] or len(final) != spec["length"]:
            raise ValueError(
                f"Unexpected reconstructed length for {spec['key']} cluster {cluster_id}"
            )
        guide_ed = edit_distance(guide, center)
        final_ed = edit_distance(final, center)
        records.append(
            {
                "cluster_id": cluster_id,
                "anomaly": float(row["joint_anomaly"]),
                "active": int(row["joint_active"]),
                "changed": int(guide != final),
                "guide": guide,
                "final": final,
                "center": center,
                "guide_ed": guide_ed,
                "final_ed": final_ed,
                "guide_failed": int(guide_ed != 0),
                "final_failed": int(final_ed != 0),
            }
        )

    ranked = sorted(
        range(len(records)), key=lambda index: (records[index]["anomaly"], index)
    )
    for rank, index in enumerate(ranked):
        records[index]["decile"] = min(bins, rank * bins // len(records) + 1)

    auc = binary_auc(
        [record["anomaly"] for record in records],
        [record["guide_failed"] for record in records],
    )
    activated_percent = safe_percent(
        sum(record["active"] for record in records), len(records)
    )

    gate_rows = []
    for bin_id in range(1, bins + 1):
        subset = [record for record in records if record["decile"] == bin_id]
        gate_rows.append(
            {
                "record_type": "gate",
                "dataset_short": spec["short"],
                "order": bin_id,
                "decile": bin_id,
                "clusters": len(subset),
                "guide_failure_percent": safe_percent(
                    sum(record["guide_failed"] for record in subset), len(subset)
                ),
                "final_failure_percent": safe_percent(
                    sum(record["final_failed"] for record in subset), len(subset)
                ),
                "structured_decoding_percent": safe_percent(
                    sum(record["active"] for record in subset), len(subset)
                ),
                "auc": auc,
                "activated_percent": activated_percent,
                "error_type": "",
                "n": "",
                "share_percent": "",
                "exact_rescue_percent": "",
                "partial_improvement_percent": "",
                "unchanged_percent": "",
                "worsened_percent": "",
                "any_improvement_percent": "",
                "ed_reduction_percent": "",
            }
        )

    selected_failures = [
        record
        for record in records
        if record["active"] == 1 and record["guide_failed"] == 1
    ]
    for record in selected_failures:
        record["error_type"] = classify_guide_error(
            record["guide"], record["center"], record["guide_ed"]
        )
        record["outcome"] = outcome(record["guide_ed"], record["final_ed"])

    phase_rows = []
    total_selected_failures = len(selected_failures)
    for order, error_type in enumerate(ERROR_ORDER, 1):
        subset = [
            record
            for record in selected_failures
            if record["error_type"] == error_type
        ]
        n = len(subset)
        outcome_counts = {
            name: sum(record["outcome"] == name for record in subset)
            for name in ["exact", "partial", "unchanged", "worsened"]
        }
        guide_ed_sum = sum(record["guide_ed"] for record in subset)
        final_ed_sum = sum(record["final_ed"] for record in subset)
        phase_rows.append(
            {
                "record_type": "phase",
                "dataset_short": spec["short"],
                "order": order,
                "decile": "",
                "clusters": "",
                "guide_failure_percent": "",
                "final_failure_percent": "",
                "structured_decoding_percent": "",
                "auc": "",
                "activated_percent": "",
                "error_type": error_type,
                "n": n,
                "share_percent": safe_percent(n, total_selected_failures),
                "exact_rescue_percent": safe_percent(outcome_counts["exact"], n),
                "partial_improvement_percent": safe_percent(
                    outcome_counts["partial"], n
                ),
                "unchanged_percent": safe_percent(outcome_counts["unchanged"], n),
                "worsened_percent": safe_percent(outcome_counts["worsened"], n),
                "any_improvement_percent": safe_percent(
                    outcome_counts["exact"] + outcome_counts["partial"], n
                ),
                "ed_reduction_percent": safe_percent(
                    guide_ed_sum - final_ed_sum, guide_ed_sum
                ),
            }
        )

    print(
        f"{spec['key']}: clusters={len(records)}, active={sum(r['active'] for r in records)}, "
        f"selected guide failures={len(selected_failures)}, AUC={auc:.12f}"
    )
    return gate_rows + phase_rows


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def compare_csv(generated: Path, expected: Path, tolerance: float):
    with generated.open(newline="", encoding="utf-8") as handle:
        actual_rows = list(csv.DictReader(handle))
    with expected.open(newline="", encoding="utf-8") as handle:
        expected_rows = list(csv.DictReader(handle))
    if len(actual_rows) != len(expected_rows):
        raise ValueError(
            f"CSV row-count mismatch: generated={len(actual_rows)}, expected={len(expected_rows)}"
        )
    mismatches = []
    for row_number, (actual, expected_row) in enumerate(
        zip(actual_rows, expected_rows), start=2
    ):
        for field in FIELDS:
            left = actual.get(field, "")
            right = expected_row.get(field, "")
            if left == right:
                continue
            try:
                equal = math.isclose(
                    float(left), float(right), rel_tol=tolerance, abs_tol=tolerance
                )
            except (TypeError, ValueError):
                equal = False
            if not equal:
                mismatches.append((row_number, field, left, right))
                if len(mismatches) >= 20:
                    break
        if len(mismatches) >= 20:
            break
    if mismatches:
        lines = ["Generated CSV differs from the expected CSV:"]
        lines.extend(
            f"  row {row}, {field}: generated={left!r}, expected={right!r}"
            for row, field, left, right in mismatches
        )
        raise ValueError("\n".join(lines))
    print(f"[PASS] Generated CSV matches {expected}")


def main():
    args = parse_args()
    if args.bins < 3:
        raise ValueError("--bins must be at least 3")
    datasets = [
        {
            "key": "srinivas_110",
            "short": "Srinivas",
            "prefix": "srinivas",
            "length": 110,
            "centers": args.srinivas_centers,
        },
        {
            "key": "chandak_108",
            "short": "Chandak-108",
            "prefix": "chandak108",
            "length": 108,
            "centers": args.chandak_centers,
        },
    ]
    rows = []
    for spec in datasets:
        rows.extend(analyze_dataset(spec, args.runs_dir, args.bins))
    write_csv(args.output, rows)
    print(f"Wrote {len(rows)} rows to {args.output}")
    if args.compare:
        compare_csv(args.output, args.compare, args.tolerance)


if __name__ == "__main__":
    main()