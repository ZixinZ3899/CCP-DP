# DNAFormer baseline

This directory contains the CCP-DP data-preparation and fine-tuning adapters
used to evaluate DNAFormer.

## Upstream source

DNAFormer is provided inside the pinned TReconLM source tree:

- TReconLM commit: ea8a838e03e03c5eecaa6d11ceb2b9b8c5262b42
- Local source:
  build/third_party/TReconLM/DeepLearningBaselines/DNAFormer

Fetch the source with:

    ./baselines/deep/fetch_sources.sh

The local adapters are:

- prepare_dnaformer_splits.py
- dnaformer_finetune_custom.py

## Environment

DNAFormer was verified in the same CUDA environment used for TReconLM:

    conda env create \
      -f baselines/deep/treconlm/environment.yml

    conda activate treconlm

The verified hardware was an NVIDIA GeForce RTX 3090.

## Official initialization checkpoint

Place the official checkpoint at:

    build/third_party/TReconLM/DeepLearningBaselines/DNAFormer/checkpoints_official/DNAFormer_siamese.pth

Verify it with:

    sha256sum -c baselines/deep/dnaformer/MODEL_SHA256SUMS

Expected SHA-256:

    b63abe3ad041c6bfb7b5f4595809dc38e8378110e0ff2235a00e546cce8d5396

## Deterministic data preparation

Set the raw dataset root:

    DATA_ROOT=/path/to/clustered-nanopore-reads-dataset

Generate the Srinivas split:

    python baselines/deep/dnaformer/prepare_dnaformer_splits.py \
      --clusters "$DATA_ROOT/Clusters_no_empty.txt" \
      --centers "$DATA_ROOT/Centers_no_empty.txt" \
      --output-dir data/generated/dnaformer/srinivas_seed2026 \
      --target-len 110 \
      --seed 2026 \
      --dataset-name srinivas_110

Generate the Chandak split:

    python baselines/deep/dnaformer/prepare_dnaformer_splits.py \
      --clusters \
        "$DATA_ROOT/bar_lev_chandak_bbs_benchmark/Chandak_et_al/Clusters.txt" \
      --centers \
        "$DATA_ROOT/bar_lev_chandak_bbs_benchmark/Chandak_et_al/Centers.txt" \
      --output-dir data/generated/dnaformer/chandak_seed2026 \
      --target-len 108 \
      --seed 2026 \
      --dataset-name chandak_108

The preparation procedure:

1. retains clusters containing at least two reads;
2. shuffles their original zero-based indices with seed 2026;
3. assigns 80% to training and 10% to validation;
4. assigns the remaining clusters to testing;
5. writes each example as reads joined by "|" followed by ":" and the center.

Verified split sizes:

| Dataset | Train | Validation | Test | Excluded single-read |
|---|---:|---:|---:|---:|
| Srinivas | 7977 | 997 | 998 | 12 |
| Chandak | 1172 | 146 | 147 | 1 |

The generated data and index files were byte-identical to the original
experimental splits.

## Srinivas 110 bp training

Run the final max-reads-16 configuration:

    mkdir -p results/reproduced/dnaformer

    /usr/bin/time -v \
      -o results/reproduced/dnaformer/srinivas_time.txt \
      env CUDA_VISIBLE_DEVICES=0 \
          OMP_NUM_THREADS=8 \
          PYTHONUNBUFFERED=1 \
          PYTHONWARNINGS=ignore::UserWarning \
      python -u baselines/deep/dnaformer/dnaformer_finetune_custom.py \
        --train data/generated/dnaformer/srinivas_seed2026/train.txt \
        --val data/generated/dnaformer/srinivas_seed2026/val.txt \
        --test data/generated/dnaformer/srinivas_seed2026/test.txt \
        --checkpoint \
          build/third_party/TReconLM/DeepLearningBaselines/DNAFormer/checkpoints_official/DNAFormer_siamese.pth \
        --output_dir \
          results/reproduced/dnaformer/srinivas_max16_seed2026 \
        --target_len 110 \
        --input_len 132 \
        --max_reads 16 \
        --epochs 50 \
        --batch_size 8 \
        --accumulation_steps 8 \
        --backbone_lr 1e-5 \
        --new_layer_lr 1e-3 \
        --weight_decay 0 \
        --patience 6 \
        --seed 2026 \
        --num_workers 0 \
        > results/reproduced/dnaformer/srinivas.log 2>&1

Verified held-out result:

- Test clusters: 998
- Exact reconstructions: 905
- Success rate: 90.6814%
- Mean edit distance: 0.3276553
- Best epoch: 30
- Training time: 2422.12 seconds
- Test inference time: 2.94 seconds

## Chandak 108 bp training

Run the final max-reads-16 configuration:

    /usr/bin/time -v \
      -o results/reproduced/dnaformer/chandak_time.txt \
      env CUDA_VISIBLE_DEVICES=0 \
          OMP_NUM_THREADS=8 \
          PYTHONUNBUFFERED=1 \
          PYTHONWARNINGS=ignore::UserWarning \
      python -u baselines/deep/dnaformer/dnaformer_finetune_custom.py \
        --train data/generated/dnaformer/chandak_seed2026/train.txt \
        --val data/generated/dnaformer/chandak_seed2026/val.txt \
        --test data/generated/dnaformer/chandak_seed2026/test.txt \
        --checkpoint \
          build/third_party/TReconLM/DeepLearningBaselines/DNAFormer/checkpoints_official/DNAFormer_siamese.pth \
        --output_dir \
          results/reproduced/dnaformer/chandak_max16_seed2026 \
        --target_len 108 \
        --input_len 132 \
        --max_reads 16 \
        --epochs 50 \
        --batch_size 8 \
        --accumulation_steps 8 \
        --backbone_lr 1e-5 \
        --new_layer_lr 1e-3 \
        --weight_decay 0 \
        --patience 6 \
        --seed 2026 \
        --num_workers 0 \
        > results/reproduced/dnaformer/chandak.log 2>&1

Verified held-out result:

- Test clusters: 147
- Exact reconstructions: 17
- Success rate: 11.5646%
- Mean edit distance: 5.7891156
- Best epoch: 32
- Test inference time: 0.85 seconds

## Historical Chandak test-path correction

The original max-reads-16 Chandak training command accidentally supplied the
validation file through the test argument. This did not affect optimization,
early stopping, or best-checkpoint selection because the test set was used
only after training.

The saved best checkpoint was subsequently evaluated on the true held-out
test split. The preserved actual-test files are the reported Chandak result.

The reproduction command above directly supplies the correct test.txt file.

## Fine-tuned checkpoint records

Locally preserved checkpoints are expected at:

    build/models/dnaformer/srinivas_max16_seed2026.pth
    build/models/dnaformer/chandak_max16_seed2026.pth

Verify them with:

    sha256sum -c \
      baselines/deep/dnaformer/FINETUNED_MODEL_SHA256SUMS

These model files are not committed to the repository.

## Evaluation

Extract held-out centers from a generated test file:

    awk -F: '{print $NF}' \
      data/generated/dnaformer/srinivas_seed2026/test.txt \
      > results/reproduced/dnaformer/srinivas_test_centers.txt

Evaluate predictions with:

    python scripts/evaluate_answers_extended.py \
      -o results/reproduced/dnaformer/srinivas_max16_seed2026/test_predictions.txt \
      -a results/reproduced/dnaformer/srinivas_test_centers.txt

Use the corresponding Chandak test file for the Chandak evaluation.

## Reporting rule

DNAFormer values are held-out-subset results:

- Srinivas: 998 of 9984 clusters
- Chandak: 147 of 1466 clusters

They are not directly comparable with full-dataset CCP-DP results unless
CCP-DP and all other methods are evaluated using the same saved test indices.
