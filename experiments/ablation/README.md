# CCP-DP ablation study

This experiment measures the contribution and computational cost of the main
CCP-DP components on real DNA storage datasets.

## Evaluated components

The ablation study considers four conceptual components:

- joint phase-aware structured decoding;
- cross-cluster leave-one-out contextual score;
- blind IDS-channel likelihood;
- sparse routing gate.

The configured variants are:

| Variant | Description |
|---|---|
| Guide only | Use the frozen local guide without structured decoding |
| Full w/o Cross-LOO score | Remove the cross-cluster contextual score |
| Full w/o blind IDS likelihood | Remove the blind IDS likelihood contribution |
| Full with decoding of all clusters | Disable sparse routing and decode every cluster |
| Local trellis only | Use structured local decoding without the full evidence combination |
| Full | Complete CCP-DP method |

The exact command-line arguments for each variant are defined in
`ablation_config.json`.

## Datasets

The default configuration evaluates:

| Dataset | Length | Clusters | Centers |
|---|---:|---|---|
| Srinivas-110 | 110 bp | `data/processed/srinivas_110/Clusters.txt` | `data/processed/srinivas_110/Centers.txt` |
| Chandak-108 | 108 bp | `data/processed/chandak_108/Clusters.txt` | `data/processed/chandak_108/Centers.txt` |
| Bar-Lev-140 | 140 bp | `data/processed/bar_lev_140/Clusters.txt` | `data/processed/bar_lev_140/Centers.txt` |

Run commands from the repository root.

## Files

- `ccpdp_ablation.hpp`: instrumented CCP-DP implementation
- `main_ablation.cpp`: ablation command-line executable
- `evaluate_ablation.cpp`: standalone evaluation program
- `ablation_config.json`: datasets and ablation variants
- `run_ablation.py`: experiment runner and result collector
- `run_ablation.sh`: build-and-run entry point
- `plot_results.py`: summary figure generator
- `run_adaptive_routing.py`: shared adaptive-routing experiment helper
- `requirements.txt`: Python plotting requirements

## Prerequisites

Prepare the datasets:

    ./data/download_datasets.sh --core

Install the plotting dependencies when needed:

    python -m pip install -r experiments/ablation/requirements.txt

The shell entry point compiles the two required C++ programs automatically.

## Smoke test

The smoke test processes the first 100 clusters, uses one repeat, and runs all
six variants on all three configured datasets.

Use a small thread count for validation:

    OUTDIR=/tmp/ccpdp_ablation_smoke \
    JOBS=2 \
    experiments/ablation/run_ablation.sh smoke

A successful smoke test performs 18 runs:

    3 datasets x 6 variants x 1 repeat = 18 runs

Expected top-level outputs include:

    /tmp/ccpdp_ablation_smoke/ablation_summary.csv
    /tmp/ccpdp_ablation_smoke/ablation_table.md
    /tmp/ccpdp_ablation_smoke/raw_runs.csv

Smoke-test accuracy must not be reported as a complete benchmark result.

## Full experiment

Run the complete experiment with:

    JOBS=32 \
    REPEATS=3 \
    experiments/ablation/run_ablation.sh full

The default full output directory is:

    experiments/ablation/results/reproduced

Override it when necessary:

    OUTDIR=/path/to/output \
    JOBS=32 \
    REPEATS=3 \
    experiments/ablation/run_ablation.sh full

The complete experiment can take substantially longer than the smoke test.

## Direct runner usage

The Python runner can also be invoked directly after compiling the two
instrumented executables:

    python experiments/ablation/run_ablation.py \
      --config experiments/ablation/ablation_config.json \
      --binary build/ccpdp_ablation \
      --evaluator build/evaluate_ablation \
      --outdir /path/to/output \
      --jobs 32 \
      --repeats 3

Use `--limit N` for a small debugging run.

## Outputs

For each dataset and variant, the runner records:

- reconstructed sequences;
- diagnostic CSV output;
- execution log;
- repeat identifier;
- exact reconstruction statistics;
- reconstruction and edit-distance metrics;
- execution time.

The aggregated files are:

- `raw_runs.csv`: one row per dataset, variant, and repeat;
- `ablation_summary.csv`: metrics aggregated across repeats;
- `ablation_table.md`: compact Markdown table;
- generated figure files from `plot_results.py`.

## Rebuild figures only

If `ablation_summary.csv` already exists:

    python experiments/ablation/plot_results.py \
      --summary /path/to/ablation_summary.csv \
      --outdir /path/to/figure_directory

## Reproducibility notes

- Dataset paths are repository-relative.
- The guide-only condition does not perform structured decoding.
- Decode-all is a computational ablation and normally increases runtime.
- Timing values depend on the processor, compiler, and thread placement.
- Compare accuracy across methods using the same dataset order and centers.
- Do not compare a limited smoke run directly with a complete historical run,
  because dataset-level calibration and routing depend on the full collection
  of clusters.
