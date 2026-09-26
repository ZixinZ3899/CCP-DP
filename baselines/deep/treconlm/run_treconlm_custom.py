#!/usr/bin/env python3
import argparse
import csv
import gc
import json
import pickle
import random
import resource
import sys
import time
from pathlib import Path

import numpy as np
import torch
from Levenshtein import distance as levenshtein_distance
from tqdm import tqdm


ADAPTER_DIR = Path(__file__).resolve().parent
REPO_ROOT = ADAPTER_DIR.parents[2]
PROJECT_ROOT = REPO_ROOT / "build" / "third_party" / "TReconLM"
sys.path.insert(0, str(PROJECT_ROOT))

from src.eval_pkg.GPT_Inference import GPT_Inference
from src.gpt_pkg.model import GPT, GPTConfig
from src.utils.helper_functions import filter_string


DNA = set("ACGT")
OUT_SEPARATOR = "==============================="


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run a released TReconLM checkpoint on a clustered DNA dataset."
    )
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--centers", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--target-len", type=int, required=True)
    parser.add_argument("--block-size", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--separator-prefix", default="====")
    parser.add_argument("--max-reads", type=int, default=10)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--device", default="cuda:0")
    return parser.parse_args()


def load_clusters(path, separator_prefix):
    clusters = []
    current = []
    with open(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip().upper()
            if not line:
                continue
            if line.startswith(separator_prefix):
                if current:
                    clusters.append(current)
                    current = []
                continue
            invalid = set(line) - DNA
            if invalid:
                raise ValueError(
                    f"Invalid characters {sorted(invalid)} in {path}, line {line_number}"
                )
            current.append(line)
    if current:
        clusters.append(current)
    return clusters


def load_centers(path, target_len):
    centers = []
    with open(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            center = raw.strip().upper()
            if not center:
                continue
            invalid = set(center) - DNA
            if invalid:
                raise ValueError(
                    f"Invalid characters {sorted(invalid)} in {path}, line {line_number}"
                )
            if len(center) != target_len:
                raise ValueError(
                    f"Center length {len(center)} != {target_len} in {path}, line {line_number}"
                )
            centers.append(center)
    return centers


def select_supported_examples(clusters, centers, max_reads, block_size, seed):
    if len(clusters) != len(centers):
        raise ValueError(
            f"Cluster/center count mismatch: {len(clusters)} vs {len(centers)}"
        )

    examples = []
    unsupported = []
    context_reduced = 0

    for index, (reads, center) in enumerate(zip(clusters, centers)):
        if len(reads) < 2:
            unsupported.append(index)
            continue

        rng = random.Random(seed + index)
        chosen = list(reads)
        if len(chosen) > max_reads:
            chosen = rng.sample(chosen, max_reads)
        else:
            rng.shuffle(chosen)

        original_n = len(chosen)
        while len(chosen) > 2 and len("|".join(chosen)) + 1 > block_size:
            chosen.pop()
        if len(chosen) != original_n:
            context_reduced += 1

        prefix_length = len("|".join(chosen)) + 1
        if prefix_length > block_size:
            raise ValueError(
                f"Cluster {index} still exceeds block size {block_size} with two reads "
                f"(prefix length {prefix_length})"
            )

        examples.append(
            {
                "original_index": index,
                "reads": chosen,
                "center": center,
                "prefix_length": prefix_length,
            }
        )

    return examples, unsupported, context_reduced


def save_selected_data(output_dir, examples, unsupported):
    with open(output_dir / "selected_clusters.txt", "w") as reads_handle, open(
        output_dir / "selected_centers.txt", "w"
    ) as centers_handle, open(output_dir / "selected_indices.txt", "w") as indices_handle:
        for item in examples:
            for read in item["reads"]:
                reads_handle.write(read + "\n")
            reads_handle.write(OUT_SEPARATOR + "\n")
            centers_handle.write(item["center"] + "\n")
            indices_handle.write(str(item["original_index"]) + "\n")

    with open(output_dir / "unsupported_single_read_indices.txt", "w") as handle:
        for index in unsupported:
            handle.write(str(index) + "\n")


def load_vocabulary():
    with open(PROJECT_ROOT / "src" / "data_pkg" / "meta_nuc.pkl", "rb") as handle:
        meta = pickle.load(handle)
    return meta["stoi"], meta["itos"]


def load_model(checkpoint_path, device):
    checkpoint = torch.load(checkpoint_path, map_location="cpu")
    model_args = checkpoint["model_args"]
    config_args = {key: value for key, value in model_args.items() if key != "model_type"}

    state_dict = {}
    for key, value in checkpoint["model"].items():
        clean_key = key.replace("_orig_mod.", "") if key.startswith("_orig_mod.") else key
        state_dict[clean_key] = value

    model = GPT(GPTConfig(**config_args))
    model.load_state_dict(state_dict, strict=True)
    model = model.half().to(device).eval()
    parameter_count = sum(parameter.numel() for parameter in model.parameters())

    del checkpoint, state_dict
    gc.collect()
    return model, parameter_count


def run_inference(args, examples, model, device):
    stoi, itos = load_vocabulary()
    decode = lambda tokens: "".join(itos[token] for token in tokens)
    encode = lambda text: [stoi.get(character, stoi.get("<unk>", 0)) for character in text]

    base_params = {
        "model": model,
        "device": device,
        "stoi": stoi,
        "itos": itos,
        "encode": encode,
        "decode": decode,
        "temperature": 1.0,
        "greedy": True,
        "ground_truth_length": args.target_len,
        "block_size": args.block_size,
        "target_type": "CPRED",
        "constrained_generation": True,
    }

    sorted_examples = sorted(examples, key=lambda item: item["prefix_length"])
    results = []

    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats(device)
    torch.cuda.synchronize(device)
    start_time = time.perf_counter()

    with torch.inference_mode():
        for start in tqdm(
            range(0, len(sorted_examples), args.batch_size), desc="TReconLM inference"
        ):
            batch = sorted_examples[start : start + args.batch_size]
            prompts = ["|".join(item["reads"]) + ":" for item in batch]
            cluster_sizes = [len(item["reads"]) for item in batch]

            params = dict(base_params)
            params["ctx"] = torch.amp.autocast(
                "cuda",
                dtype=torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float16,
            )
            output = GPT_Inference(params).inference(
                prompts, alignment_size=cluster_sizes
            )

            for item, prediction in zip(batch, output["candidate_sequences"]):
                prediction = filter_string(prediction)[: args.target_len]
                edit_distance = levenshtein_distance(item["center"], prediction)
                positional_distance = sum(
                    left != right for left, right in zip(item["center"], prediction)
                ) + abs(len(item["center"]) - len(prediction))
                results.append(
                    {
                        "original_index": item["original_index"],
                        "cluster_size": len(item["reads"]),
                        "prediction": prediction,
                        "center": item["center"],
                        "edit_distance": edit_distance,
                        "hamming_distance": positional_distance,
                        "success": int(prediction == item["center"]),
                    }
                )

    torch.cuda.synchronize(device)
    inference_seconds = time.perf_counter() - start_time
    peak_allocated_mib = torch.cuda.max_memory_allocated(device) / (1024 ** 2)
    peak_reserved_mib = torch.cuda.max_memory_reserved(device) / (1024 ** 2)

    results.sort(key=lambda item: item["original_index"])
    return results, inference_seconds, peak_allocated_mib, peak_reserved_mib


def save_results(
    args,
    output_dir,
    results,
    total_clusters,
    unsupported,
    context_reduced,
    inference_seconds,
    peak_allocated_mib,
    peak_reserved_mib,
    parameter_count,
):
    with open(output_dir / "predictions.txt", "w") as handle:
        for item in results:
            handle.write(item["prediction"] + "\n")

    with open(output_dir / "results.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    supported_count = len(results)
    success_count = sum(item["success"] for item in results)
    total_edit_distance = sum(item["edit_distance"] for item in results)
    total_hamming_distance = sum(item["hamming_distance"] for item in results)

    summary = {
        "checkpoint": str(Path(args.checkpoint).resolve()),
        "target_length": args.target_len,
        "seed": args.seed,
        "max_reads": args.max_reads,
        "batch_size": args.batch_size,
        "total_original_clusters": total_clusters,
        "supported_clusters": supported_count,
        "unsupported_single_read_clusters": len(unsupported),
        "context_reduced_clusters": context_reduced,
        "success_count": success_count,
        "success_rate_supported": success_count / supported_count,
        "success_rate_full_conservative": success_count / total_clusters,
        "total_edit_distance_supported": total_edit_distance,
        "average_edit_distance_supported": total_edit_distance / supported_count,
        "average_normalized_edit_distance_supported": total_edit_distance
        / (supported_count * args.target_len),
        "average_hamming_distance_supported": total_hamming_distance / supported_count,
        "inference_seconds": inference_seconds,
        "clusters_per_second": supported_count / inference_seconds,
        "peak_gpu_allocated_mib": peak_allocated_mib,
        "peak_gpu_reserved_mib": peak_reserved_mib,
        "model_parameters": parameter_count,
        "process_max_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
    }

    with open(output_dir / "summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)

    print("\n" + "=" * 72)
    print(json.dumps(summary, indent=2))
    print("=" * 72)
    print(f"Predictions: {output_dir / 'predictions.txt'}")
    print(f"Selected centers: {output_dir / 'selected_centers.txt'}")
    print(f"Selected clusters: {output_dir / 'selected_clusters.txt'}")


def main():
    args = parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is not available; this benchmark is intended for GPU inference")
    if not (2 <= args.max_reads <= 10):
        raise ValueError("--max-reads must be between 2 and 10")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    device = torch.device(args.device)

    clusters = load_clusters(args.clusters, args.separator_prefix)
    centers = load_centers(args.centers, args.target_len)
    examples, unsupported, context_reduced = select_supported_examples(
        clusters, centers, args.max_reads, args.block_size, args.seed
    )
    if not examples:
        raise ValueError("No supported clusters remain")

    print(f"Original clusters: {len(clusters)}")
    print(f"Supported clusters: {len(examples)}")
    print(f"Unsupported single-read clusters: {len(unsupported)}")
    print(f"Clusters reduced further to fit context: {context_reduced}")
    save_selected_data(output_dir, examples, unsupported)

    model, parameter_count = load_model(args.checkpoint, device)
    print(f"Model parameters: {parameter_count:,}")
    print(f"Device: {device} ({torch.cuda.get_device_name(device)})")

    results, seconds, peak_allocated, peak_reserved = run_inference(
        args, examples, model, device
    )
    save_results(
        args,
        output_dir,
        results,
        len(clusters),
        unsupported,
        context_reduced,
        seconds,
        peak_allocated,
        peak_reserved,
        parameter_count,
    )


if __name__ == "__main__":
    main()