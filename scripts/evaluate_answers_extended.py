#!/usr/bin/env python3
"""Evaluate DNA strand reconstruction results.

The prediction file and ground-truth file must contain one sequence per line
in the same order. Empty lines are ignored.
"""

import argparse
from pathlib import Path

try:
    import Levenshtein as _levenshtein

    def levenshtein_distance(left: str, right: str) -> int:
        return _levenshtein.distance(left, right)

except ImportError:
    def levenshtein_distance(left: str, right: str) -> int:
        """Dependency-free fallback using two dynamic-programming rows."""
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


def read_sequences(path: str):
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip().upper() for line in handle if line.strip()]


def benchmark(output_file: str, answer_file: str):
    reconstructed = read_sequences(output_file)
    answers = read_sequences(answer_file)

    if not answers:
        raise ValueError("The ground-truth file is empty.")
    if len(reconstructed) != len(answers):
        raise ValueError(
            "Sequence count mismatch: "
            f"predictions={len(reconstructed)}, ground_truth={len(answers)}"
        )

    edit_distances = [
        levenshtein_distance(prediction, truth)
        for prediction, truth in zip(reconstructed, answers)
    ]

    num_strands = len(answers)
    total_bases = sum(len(sequence) for sequence in answers)
    total_edit_distance = sum(edit_distances)
    exact_correct = sum(
        prediction == truth
        for prediction, truth in zip(reconstructed, answers)
    )

    average_edit_distance = total_edit_distance / num_strands
    normalized_edit_distance = total_edit_distance / total_bases
    reconstruction_rate = 1.0 - normalized_edit_distance
    success_rate = exact_correct / num_strands

    # Kept only for compatibility with the old evaluator. Under insertions or
    # deletions this is a positional mismatch count, not alignment-based ED.
    total_positional_errors = 0
    for prediction, truth in zip(reconstructed, answers):
        overlap = min(len(prediction), len(truth))
        total_positional_errors += sum(
            prediction[index] != truth[index] for index in range(overlap)
        )
        total_positional_errors += abs(len(prediction) - len(truth))
    average_positional_distance = total_positional_errors / num_strands

    return {
        "num_strands": num_strands,
        "total_bases": total_bases,
        "exact_correct": exact_correct,
        "total_edit_distance": total_edit_distance,
        "hamming_distance": average_positional_distance,
        "edit_distance": average_edit_distance,
        "normalized_edit_distance": normalized_edit_distance,
        "reconstruction_rate": reconstruction_rate,
        "success_rate": success_rate,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark DNA sequence reconstruction results"
    )
    parser.add_argument(
        "-o", "--output_file", required=True, help="Reconstructed sequences"
    )
    parser.add_argument(
        "-a", "--answer_file", required=True, help="Ground-truth sequences"
    )
    args = parser.parse_args()

    metrics = benchmark(args.output_file, args.answer_file)

    print(Path(args.output_file).name)
    print(
        f"Sequences: {metrics['num_strands']} "
        f"(exact: {metrics['exact_correct']})"
    )
    print(f"Hamming Distance:  {metrics['hamming_distance']:.12f}")
    print(f"Edit Distance:  {metrics['edit_distance']:.12f}")
    print(f"Reconstruction Rate:  {metrics['reconstruction_rate']:.12f}")
    print(f"Success Rate:  {metrics['success_rate']:.12f}")


if __name__ == "__main__":
    main()
