# CCP-DP experiments

This directory contains the scripts, configurations, saved source data, and
selected outputs used to evaluate CCP-DP.

Run commands from the repository root unless an experiment-specific README
states otherwise.

## Experiment index

| Directory | Purpose | Main entry point |
|---|---|---|
| `ablation/` | Component ablation on real datasets | `run_ablation.sh` |
| `coverage_robustness/` | Nested coverage-downsampling evaluation | `tools/run_ccpdp_coverage.py` |
| `coverage_ids_conditions/` | Controlled coverage-by-IDS error grid | `run_srinivas_stress.sh` |
| `cross_cluster_context/` | Counterfactual context and adaptive routing analysis | `run_adaptive_all.sh` |
| `difficult_cluster_analysis/` | Mechanism analysis on difficult clusters | `run_experiment.sh` |
| `parameter_sensitivity/` | Sensitivity to CCP-DP parameters | `run_sensitivity.sh` |
| `deep_baselines/` | Preserved deep-learning evaluation records | `README.md` |

## Shared prerequisites

Prepare the public datasets:

    ./data/download_datasets.sh --core

Build the classical baselines:

    ./baselines/build_baselines.sh

Build the main CCP-DP executable:

    mkdir -p build

    g++ \
      -std=c++17 \
      -O3 \
      -fopenmp \
      -Iinclude \
      src/main.cpp \
      -o build/ccpdp

Verify it with:

    ./build/ccpdp --help

Some experiments compile an instrumented CCP-DP executable. These builds are
written under `build/` and do not replace the main source files.

## Reproducibility levels

Experiment-specific documentation distinguishes between:

- `smoke`: a small validation run;
- `full`: the complete experiment;
- `plot`: figure generation from preserved source data.

Smoke tests verify integration and output structure. Their accuracy values
must not be reported as full benchmark results.

## Output policy

Experiment outputs fall into three categories:

1. Small tables, source data, and figures needed to document reported results
   may be committed.
2. Large generated datasets, binaries, model weights, and third-party source
   trees remain under ignored directories.
3. Temporary validation outputs should be written under `/tmp` or `build/`.

## Evaluation

A prediction file containing one reconstructed sequence per line can be
evaluated with:

    python scripts/evaluate_answers_extended.py \
      -o PREDICTIONS.txt \
      -a CENTERS.txt

The evaluator reports:

- exact reconstruction count;
- Hamming distance;
- edit distance;
- reconstruction rate;
- exact success rate.

## Experiment-specific documentation

Detailed commands, inputs, outputs, and resource notes are provided in the
README inside each experiment directory.

Full experiments can require substantially more CPU or GPU time than smoke
tests. Check the experiment-specific README before starting a full run.
