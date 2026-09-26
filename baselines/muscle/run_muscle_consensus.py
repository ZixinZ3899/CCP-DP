#!/usr/bin/env python3

import os
import time
import argparse
import tempfile
import subprocess
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed


def load_clusters(path, separator):
    clusters = []
    current = []

    with open(path, "r") as f:
        for raw_line in f:
            line = raw_line.strip()

            if not line:
                continue

            # 同时兼容 ==== 和更长的等号分隔符
            if line.startswith(separator):
                if current:
                    clusters.append(current)
                    current = []
            else:
                sequence = line.upper()

                if not all(base in "ACGTN" for base in sequence):
                    raise ValueError(
                        f"发现非法DNA字符：{sequence[:80]}"
                    )

                current.append(sequence)

    if current:
        clusters.append(current)

    return clusters


def write_fasta(reads, output_path):
    with open(output_path, "w") as f:
        for index, sequence in enumerate(reads):
            f.write(f">read_{index}\n")
            f.write(sequence + "\n")


def read_aligned_fasta(path):
    sequences = []
    current = []

    with open(path, "r") as f:
        for raw_line in f:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                    current = []
            else:
                current.append(line)

    if current:
        sequences.append("".join(current).upper())

    return sequences


def majority_consensus(aligned_sequences):
    """
    将 gap 当作第五种符号参与投票。

    如果某一列 gap 的数量不少于最高碱基数量，
    删除该列；否则输出支持度最高的 A/C/G/T。

    碱基票数相同时按照 A、C、G、T 的顺序确定，
    以保证结果可复现。
    """
    if not aligned_sequences:
        return ""

    alignment_length = len(aligned_sequences[0])

    for sequence in aligned_sequences:
        if len(sequence) != alignment_length:
            raise RuntimeError("MUSCLE输出的序列长度不一致")

    consensus = []
    base_order = "ACGT"

    for column_index in range(alignment_length):
        column = [
            sequence[column_index]
            for sequence in aligned_sequences
        ]

        counts = Counter(column)

        gap_count = counts.get("-", 0)

        best_base = max(
            base_order,
            key=lambda base: (counts.get(base, 0), -base_order.index(base))
        )
        best_base_count = counts.get(best_base, 0)

        # gap获胜或与碱基打平时，删除该比对列
        if gap_count >= best_base_count:
            continue

        consensus.append(best_base)

    return "".join(consensus)


def process_cluster(task):
    (
        cluster_index,
        reads,
        muscle_path,
        muscle_threads,
        temporary_root,
    ) = task

    if len(reads) == 0:
        return cluster_index, ""

    if len(reads) == 1:
        return cluster_index, reads[0]

    with tempfile.TemporaryDirectory(
        prefix=f"cluster_{cluster_index}_",
        dir=temporary_root,
    ) as cluster_dir:

        input_fasta = os.path.join(cluster_dir, "reads.fa")
        output_fasta = os.path.join(cluster_dir, "aligned.afa")

        write_fasta(reads, input_fasta)

        command = [
            muscle_path,
            "-align", input_fasta,
            "-output", output_fasta,
            "-nt",
            "-threads", str(muscle_threads),
        ]

        environment = os.environ.copy()
        environment["OMP_NUM_THREADS"] = str(muscle_threads)

        result = subprocess.run(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            env=environment,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"cluster {cluster_index} MUSCLE运行失败：\n"
                f"{result.stderr}"
            )

        aligned_sequences = read_aligned_fasta(output_fasta)
        reconstructed = majority_consensus(aligned_sequences)

        return cluster_index, reconstructed


def main():
    parser = argparse.ArgumentParser(
        description="MUSCLE v5 cluster consensus reconstruction"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="簇文件路径",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="重建结果文件",
    )
    parser.add_argument(
        "-s", "--separator",
        default="====",
        help="簇分隔符，默认：====",
    )
    parser.add_argument(
        "--muscle",
        default="muscle",
        help="MUSCLE可执行文件路径",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=32,
        help="同时处理的簇数量，默认32",
    )
    parser.add_argument(
        "--muscle-threads",
        type=int,
        default=1,
        help="每个MUSCLE进程的线程数，默认1",
    )

    args = parser.parse_args()

    if args.jobs <= 0 or args.muscle_threads <= 0:
        raise ValueError("jobs和muscle-threads必须大于0")

    total_threads = args.jobs * args.muscle_threads

    print(f"input: {args.input}")
    print(f"muscle: {args.muscle}")
    print(f"parallel jobs: {args.jobs}")
    print(f"threads per MUSCLE: {args.muscle_threads}")
    print(f"maximum total threads: {total_threads}")

    start_time = time.perf_counter()

    clusters = load_clusters(args.input, args.separator)
    print(f"clusters: {len(clusters)}")

    reconstructed = [None] * len(clusters)

    with tempfile.TemporaryDirectory(
        prefix="muscle_consensus_"
    ) as temporary_root:

        tasks = [
            (
                index,
                reads,
                args.muscle,
                args.muscle_threads,
                temporary_root,
            )
            for index, reads in enumerate(clusters)
        ]

        with ProcessPoolExecutor(
            max_workers=args.jobs
        ) as executor:

            futures = [
                executor.submit(process_cluster, task)
                for task in tasks
            ]

            completed = 0

            for future in as_completed(futures):
                cluster_index, sequence = future.result()
                reconstructed[cluster_index] = sequence

                completed += 1

                if completed % 500 == 0:
                    print(
                        f"processed: {completed}/{len(clusters)}"
                    )

    with open(args.output, "w") as f:
        for sequence in reconstructed:
            f.write(sequence + "\n")

    elapsed = time.perf_counter() - start_time

    empty_count = sum(
        1 for sequence in reconstructed if not sequence
    )

    lengths = [
        len(sequence)
        for sequence in reconstructed
        if sequence
    ]

    print(f"written output: {args.output}")
    print(f"empty reconstructions: {empty_count}")

    if lengths:
        print(f"minimum output length: {min(lengths)}")
        print(f"maximum output length: {max(lengths)}")
        print(
            f"average output length: "
            f"{sum(lengths) / len(lengths):.4f}"
        )

    print(f"elapsed seconds: {elapsed:.6f}")


if __name__ == "__main__":
    main()
