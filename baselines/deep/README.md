# Deep-learning baselines

This directory contains the adapters, environment records, checksum
manifests, and reproduction instructions for the deep-learning baselines
evaluated alongside CCP-DP.

Upstream source trees and model weights are stored under build/ and are not
committed to this repository.

## Status

| Method | Source | Model weights | Validation status |
|---|---|---|---|
| TReconLM | Pinned upstream commit | Two released checkpoints | Inference smoke-tested |
| DNAFormer | Pinned upstream source | Official initialization and two fine-tuned checkpoints | Data preparation and training smoke-tested |
| RobuSeqNet | Pinned upstream commit | Not included | Source snapshot only |

## Layout

- baselines/deep/treconlm: TReconLM adapter and environment records
- baselines/deep/dnaformer: DNAFormer training and data-preparation adapters
- baselines/deep/robuseqnet: source-only status documentation
- build/third_party: downloaded upstream source trees
- build/models: local model weights
- experiments/deep_baselines: preserved evaluation results and split indices

## Source setup

Run the following command after creating the source-fetching script:

    ./baselines/deep/fetch_sources.sh

## Evaluation

Evaluate reconstruction outputs with:

    python scripts/evaluate_answers_extended.py \
      -o PREDICTIONS.txt \
      -a CENTERS.txt

TReconLM does not support single-read clusters. Such clusters are counted as
failures when reporting the conservative full-dataset success rate.

DNAFormer results use deterministic held-out test subsets. They must not be
presented as full-dataset results unless every compared method is evaluated
on the same test indices.

RobuSeqNet is retained as a pinned source snapshot only. A validated
pretrained checkpoint was unavailable, so it is excluded from the reported
benchmark results.
