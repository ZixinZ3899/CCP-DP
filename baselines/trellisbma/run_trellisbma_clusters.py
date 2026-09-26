#!/usr/bin/env python3
"""Batch runner for Microsoft's TrellisBMA implementation.

The upstream repository exposes research functions and notebooks rather than a
command-line program.  This wrapper reads separator-delimited DNA clusters,
uses the upstream uncoded (identity-code) Trellis BMA decoder, and parallelizes
independent clusters across CPU processes.
"""

import os

# One numerical-library thread per worker prevents 32 processes from creating
# hundreds of nested BLAS/OpenMP threads.
for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "NUMBA_NUM_THREADS",
):
    os.environ.setdefault(_name, "1")

import argparse
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
from tqdm import tqdm

from conv_code import conv_code
from coded_ids_multiD import coded_ids_multiD
from trellis_bma import trellis_bma


DNA_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3}
INT_TO_DNA = np.asarray(["A", "C", "G", "T"])

# These objects are constructed once in the parent and inherited by forked
# Linux workers.  Each process then owns a copy-on-write decoder instance.
_TRELLIS = None
_CODE = None
_TARGET_LEN = None
_MAX_DRIFT = None
_MAX_READS = None
_LOOKAHEAD = None
_LENGTH_POLICY = None
_SEED = None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run uncoded Trellis BMA on separator-delimited DNA clusters."
    )
    parser.add_argument("-i", "--input", required=True, help="Cluster input file")
    parser.add_argument("-o", "--output", required=True, help="One consensus per line")
    parser.add_argument("-l", "--target-len", type=int, required=True)
    parser.add_argument("-s", "--separator", default="====")
    parser.add_argument("--jobs", type=int, default=32)
    parser.add_argument(
        "--max-clusters",
        type=int,
        default=0,
        help="Process only the first N clusters; 0 means all clusters",
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=0,
        help="Use first K reads per cluster; 0 means all reads",
    )
    parser.add_argument("--max-drift", type=int, default=15)
    parser.add_argument("--p-del", type=float, default=0.01875390320814797)
    parser.add_argument("--p-ins", type=float, default=0.01674582308540412)
    parser.add_argument("--p-sub", type=float, default=0.02140018195696255)
    parser.add_argument(
        "--lookahead",
        type=int,
        choices=(0, 1),
        default=1,
        help="1: Trellis BMA with lookahead; 0: without lookahead",
    )
    parser.add_argument(
        "--length-policy",
        choices=("official-adjust", "filter"),
        default="official-adjust",
        help=(
            "How to handle reads outside target_len +/- max_drift. "
            "official-adjust reproduces upstream one_iter.py; filter discards them."
        ),
    )
    parser.add_argument("--seed", type=int, default=20260803)
    return parser.parse_args()


def load_clusters(path, separator):
    clusters = []
    current = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip().upper()
            if not line:
                continue
            if line.startswith(separator):
                if current:
                    clusters.append(current)
                    current = []
                continue
            bad = set(line) - DNA_TO_INT.keys()
            if bad:
                raise ValueError(
                    f"Invalid DNA character(s) {sorted(bad)} at input line {line_number}"
                )
            current.append(line)
    if current:
        clusters.append(current)
    if not clusters:
        raise ValueError("No non-empty clusters were found")
    return clusters


def dna_to_int(read):
    return np.fromiter((DNA_TO_INT[base] for base in read), dtype=np.int64)


def official_adjust_length(trace, rng):
    """Reproduce the out-of-drift handling used in upstream one_iter.py."""
    lower = _TARGET_LEN - _MAX_DRIFT
    upper = _TARGET_LEN + _MAX_DRIFT
    if lower <= len(trace) <= upper:
        return trace
    if len(trace) > upper:
        delete_count = len(trace) - upper
        delete_at = rng.choice(len(trace), delete_count, replace=False)
        return np.delete(trace, delete_at)

    insert_count = lower - len(trace)
    if len(trace) == 0:
        return np.zeros(lower, dtype=np.int64)
    # The official code samples without replacement.  Replacement is needed
    # only for an extreme read where more positions than its length are added.
    insert_at = rng.choice(
        len(trace), insert_count, replace=(insert_count > len(trace))
    )
    return np.insert(trace, insert_at, 0)


def prepare_traces(reads, cluster_index):
    if _MAX_READS > 0:
        reads = reads[:_MAX_READS]
    rng = np.random.default_rng(_SEED + cluster_index)
    traces = []
    for read in reads:
        trace = dna_to_int(read)
        if abs(len(trace) - _TARGET_LEN) <= _MAX_DRIFT:
            traces.append(trace)
        elif _LENGTH_POLICY == "official-adjust":
            traces.append(official_adjust_length(trace, rng))
    if not traces:
        raise ValueError("cluster has no usable reads")
    return traces


def decode_cluster(task):
    cluster_index, reads = task
    try:
        traces = prepare_traces(reads, cluster_index)
        estimate, _ = trellis_bma(
            _TRELLIS,
            traces,
            _CODE.trellis_states[0][0],
            _CODE.trellis_states[-1],
            lookahead=_LOOKAHEAD,
        )
        estimate = np.asarray(estimate, dtype=np.int64)
        if len(estimate) != _TARGET_LEN:
            raise ValueError(
                f"decoder returned length {len(estimate)}, expected {_TARGET_LEN}"
            )
        if np.any((estimate < 0) | (estimate > 3)):
            raise ValueError("decoder returned a symbol outside 0..3")
        consensus = "".join(INT_TO_DNA[estimate].tolist())
        return cluster_index, consensus, len(traces), None
    except Exception as exc:  # return context instead of losing a worker silently
        return cluster_index, None, 0, f"{type(exc).__name__}: {exc}"


def build_uncoded_decoder(args):
    global _TRELLIS, _CODE, _TARGET_LEN, _MAX_DRIFT
    global _MAX_READS, _LOOKAHEAD, _LENGTH_POLICY, _SEED

    _TARGET_LEN = args.target_len
    _MAX_DRIFT = args.max_drift
    _MAX_READS = args.max_reads
    _LOOKAHEAD = bool(args.lookahead)
    _LENGTH_POLICY = args.length_policy
    _SEED = args.seed

    code = conv_code()
    code.quar_cc(np.asarray([[1]], dtype=np.int64))
    code.make_trellis(args.target_len)
    code.make_encoder()

    trellis = coded_ids_multiD(
        4,
        4,
        code.trellis_states,
        code.trellis_edges,
        code.time_type,
        1,
        args.p_del,
        args.p_sub,
        args.p_ins,
        args.max_drift,
        input_prior=None,
    )
    _CODE = code
    _TRELLIS = trellis


def validate_args(args):
    if args.target_len <= 0 or args.target_len % 2 != 0:
        raise ValueError("target length must be a positive even integer")
    if args.jobs <= 0:
        raise ValueError("jobs must be positive")
    if args.max_reads < 0:
        raise ValueError("max-reads cannot be negative")
    if args.max_clusters < 0:
        raise ValueError("max-clusters cannot be negative")
    if args.max_drift < 0:
        raise ValueError("max-drift cannot be negative")
    probabilities = (args.p_del, args.p_ins, args.p_sub)
    if any(value < 0.0 or value >= 1.0 for value in probabilities):
        raise ValueError("IDS probabilities must be in [0, 1)")
    if sum(probabilities) >= 1.0:
        raise ValueError("p-del + p-ins + p-sub must be less than 1")
    if not args.separator:
        raise ValueError("separator cannot be empty")


def main():
    args = parse_args()
    validate_args(args)
    input_path = Path(args.input).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    clusters = load_clusters(input_path, args.separator)
    if args.max_clusters > 0:
        clusters = clusters[: args.max_clusters]
    build_uncoded_decoder(args)

    print(f"clusters: {len(clusters)}", flush=True)
    print(f"target length: {args.target_len}", flush=True)
    print(f"jobs: {args.jobs}", flush=True)
    print(f"max reads: {'all' if args.max_reads == 0 else args.max_reads}", flush=True)
    print(f"lookahead: {bool(args.lookahead)}", flush=True)
    print(
        f"IDS: del={args.p_del:.10f} ins={args.p_ins:.10f} sub={args.p_sub:.10f}",
        flush=True,
    )

    started = time.perf_counter()
    tasks = enumerate(clusters)
    if args.jobs == 1:
        results = [decode_cluster(task) for task in tqdm(tasks, total=len(clusters))]
    else:
        if "fork" not in mp.get_all_start_methods():
            raise RuntimeError("This runner requires Linux multiprocessing start method 'fork'")
        context = mp.get_context("fork")
        with context.Pool(processes=args.jobs) as pool:
            results = list(
                tqdm(
                    pool.imap(decode_cluster, tasks, chunksize=1),
                    total=len(clusters),
                    desc="Trellis BMA",
                )
            )

    errors = [(idx, message) for idx, _, _, message in results if message]
    if errors:
        print(f"ERROR: {len(errors)} cluster(s) failed", file=sys.stderr)
        for idx, message in errors[:20]:
            print(f"  cluster {idx}: {message}", file=sys.stderr)
        raise RuntimeError("No output written because one or more clusters failed")

    results.sort(key=lambda item: item[0])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    with open(temporary_path, "w", encoding="utf-8") as handle:
        for _, consensus, _, _ in results:
            handle.write(consensus + "\n")
    os.replace(temporary_path, output_path)

    elapsed = time.perf_counter() - started
    reads_used = sum(item[2] for item in results)
    print(f"reads used: {reads_used}", flush=True)
    print(f"written output: {output_path}", flush=True)
    print(f"elapsed seconds: {elapsed:.6f}", flush=True)


if __name__ == "__main__":
    main()
