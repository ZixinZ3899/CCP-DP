# Deep-learning baseline results

This directory preserves the validated outputs and deterministic split
indices used for the deep-learning baseline evaluation.

## TReconLM

| Dataset | Supported clusters | Exact | Mean ED | Supported success | Conservative full success |
|---|---:|---:|---:|---:|---:|
| Srinivas 110 bp | 9972 | 9181 | 0.1679703 | 92.0678% | 91.9571% |
| Chandak 108 bp | 1465 | 706 | 1.3829352 | 48.1911% | 48.1583% |

TReconLM does not support single-read clusters. The conservative rate counts
the 12 unsupported Srinivas clusters and one unsupported Chandak cluster as
failures.

The evaluation_supported.txt files contain independent evaluations on the
supported clusters. The full conservative rates are stored in summary.json.

## DNAFormer

| Dataset | Held-out clusters | Exact | Mean ED | Success |
|---|---:|---:|---:|---:|
| Srinivas 110 bp | 998 | 905 | 0.3276553 | 90.6814% |
| Chandak 108 bp | 147 | 17 | 5.7891156 | 11.5646% |

DNAFormer results are measured on deterministic held-out subsets generated
with seed 2026. The saved indices define the exact train, validation, and test
partitions.

The original metric JSON and console logs retain historical absolute paths
for provenance. Portable evaluation outputs are stored in evaluation.txt.

The original Chandak max-reads-16 run evaluated the validation split at the
end of training. The saved best checkpoint was subsequently evaluated on the
true held-out test split; that actual-test result is the one preserved here.

## RobuSeqNet

RobuSeqNet has no result directory because only its pinned upstream source
snapshot is retained. It is not included in the validated benchmark results.

## Reproduction instructions

See:

- baselines/deep/README.md
- baselines/deep/treconlm/README.md
- baselines/deep/dnaformer/README.md
- baselines/deep/robuseqnet/README.md
