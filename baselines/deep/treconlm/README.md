# TReconLM baseline

This directory contains the CCP-DP evaluation adapter and reproducibility
records for TReconLM.

## Upstream source

- Repository: MLI-lab/TReconLM
- Pinned commit: ea8a838e03e03c5eecaa6d11ceb2b9b8c5262b42
- Local source path: build/third_party/TReconLM

Fetch the pinned source with:

    ./baselines/deep/fetch_sources.sh

The adapter imports the upstream implementation from the local source path:

    baselines/deep/treconlm/run_treconlm_custom.py

## Environment

Create the recorded environment with:

    conda env create \
      -f baselines/deep/treconlm/environment.yml

Activate it with:

    conda activate treconlm

The platform-specific package lock is preserved in:

    baselines/deep/treconlm/environment-linux-64.lock.txt

The verified migration environment used an NVIDIA GeForce RTX 3090.

## Model checkpoints

Place the released checkpoints at:

    build/third_party/TReconLM/models/model_seq_len_110.pt
    build/third_party/TReconLM/models/model_var_len_50_120.pt

Verify them from the repository root:

    sha256sum -c baselines/deep/treconlm/MODEL_SHA256SUMS

Expected SHA-256 values:

- model_seq_len_110.pt:
  140ae26c126a74e936d6e247156ba7f3b639c710ddb5655a6c6c0b47e76f213a
- model_var_len_50_120.pt:
  e41876248325048edfa0258ad5683d4863b86883902480bec3e3675b775e1039

Model weights are not committed to this repository.

## Dataset path

Set the dataset root before running:

    DATA_ROOT=/path/to/clustered-nanopore-reads-dataset

## Srinivas 110 bp

Run inference with:

    mkdir -p results/reproduced/treconlm

    /usr/bin/time -v \
      -o results/reproduced/treconlm/srinivas_time.txt \
      env CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=8 \
      python baselines/deep/treconlm/run_treconlm_custom.py \
        --clusters "$DATA_ROOT/Clusters_no_empty.txt" \
        --centers "$DATA_ROOT/Centers_no_empty.txt" \
        --checkpoint \
          build/third_party/TReconLM/models/model_seq_len_110.pt \
        --target-len 110 \
        --block-size 1500 \
        --max-reads 10 \
        --seed 2026 \
        --batch-size 32 \
        --output-dir \
          results/reproduced/treconlm/srinivas_110 \
        > results/reproduced/treconlm/srinivas.log 2>&1

Evaluate supported clusters with:

    python scripts/evaluate_answers_extended.py \
      -o results/reproduced/treconlm/srinivas_110/predictions.txt \
      -a results/reproduced/treconlm/srinivas_110/selected_centers.txt

Verified original result:

- Original clusters: 9984
- Supported clusters: 9972
- Unsupported single-read clusters: 12
- Exact reconstructions: 9181
- Supported-cluster success rate: 92.0678%
- Conservative full-dataset success rate: 91.9571%
- Mean edit distance on supported clusters: 0.1679703
- Inference time: 322.02 seconds

## Chandak 108 bp

Run inference with:

    /usr/bin/time -v \
      -o results/reproduced/treconlm/chandak_time.txt \
      env CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=8 \
      python baselines/deep/treconlm/run_treconlm_custom.py \
        --clusters \
          "$DATA_ROOT/bar_lev_chandak_bbs_benchmark/Chandak_et_al/Clusters.txt" \
        --centers \
          "$DATA_ROOT/bar_lev_chandak_bbs_benchmark/Chandak_et_al/Centers.txt" \
        --checkpoint \
          build/third_party/TReconLM/models/model_var_len_50_120.pt \
        --target-len 108 \
        --block-size 2400 \
        --max-reads 10 \
        --seed 2026 \
        --batch-size 32 \
        --output-dir \
          results/reproduced/treconlm/chandak_108 \
        > results/reproduced/treconlm/chandak.log 2>&1

Evaluate supported clusters with:

    python scripts/evaluate_answers_extended.py \
      -o results/reproduced/treconlm/chandak_108/predictions.txt \
      -a results/reproduced/treconlm/chandak_108/selected_centers.txt

Verified original result:

- Original clusters: 1466
- Supported clusters: 1465
- Unsupported single-read clusters: 1
- Exact reconstructions: 706
- Supported-cluster success rate: 48.1911%
- Conservative full-dataset success rate: 48.1583%
- Mean edit distance on supported clusters: 1.3829352
- Inference time: 46.96 seconds

## Reporting rule

For comparison with methods evaluated on every cluster, report
success_rate_full_conservative from summary.json. This treats unsupported
single-read clusters as reconstruction failures.

The standalone evaluator compares predictions only against
selected_centers.txt and therefore reports the supported-cluster rate.
