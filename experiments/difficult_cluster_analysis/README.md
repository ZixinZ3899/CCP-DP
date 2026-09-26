# Difficult-cluster mechanism analysis

This experiment studies how CCP-DP identifies and corrects difficult read
clusters.

The analysis combines:

- frozen local guides;
- final CCP-DP reconstructions;
- per-cluster diagnostic measurements;
- true centers used only for evaluation.

It examines whether the routing score concentrates structured decoding on
clusters with high local-guide failure risk and how structured decoding changes
different error types.

## Research questions

The analysis addresses three questions:

1. Does the routing score rank difficult clusters effectively?
2. Is structured decoding concentrated in the highest-risk regions?
3. Which local-guide errors are rescued, partially improved, unchanged, or
   worsened by structured decoding?

## Datasets

The experiment uses:

| Dataset | Length | Prepared clusters | Prepared centers |
|---|---:|---|---|
| Srinivasavaradhan et al. | 110 bp | `data/processed/srinivas_110/Clusters.txt` | `data/processed/srinivas_110/Centers.txt` |
| Chandak et al. | 108 bp | `data/processed/chandak_108/Clusters.txt` | `data/processed/chandak_108/Centers.txt` |

Dataset provenance and preparation instructions are provided in
[`../../data/README.md`](../../data/README.md).

## Directory contents

The principal files are:

- `run_experiment.sh`: run CCP-DP and preserve guides and diagnostics;
- `analyze_mechanism.py`: rebuild the mechanism-analysis table;
- `plot_ccpdp_mechanism.py`: generate the mechanism figure;
- `results/ccpdp_analysis_data.csv`: preserved figure source data;
- `results/ccpdp_analysis_data_rebuilt.csv`: independently rebuilt table;
- `results/raw/`: original frozen run records;
- `results/reproduced/`: reproduced run records;
- `figures/`: preserved mechanism figures.

Run all commands from the repository root.

## Prerequisites

Prepare the datasets:

    ./data/download_datasets.sh --core

Build CCP-DP:

    mkdir -p build

    g++ \
      -std=c++17 \
      -O3 \
      -fopenmp \
      -Iinclude \
      src/main.cpp \
      -o build/ccpdp

Verify the executable:

    ./build/ccpdp --help

The analysis and plotting programs require Python 3 with pandas, NumPy, and
matplotlib.

## Generated run files

For each selected dataset, `run_experiment.sh` writes:

- `PREFIX_final.txt`: final reconstructed sequences;
- `PREFIX_diag.csv`: per-cluster CCP-DP diagnostics;
- `PREFIX_guides.txt`: frozen local guides;
- `PREFIX.log`: execution log;
- `PREFIX.command.txt`: fully resolved reconstruction command.

The prefixes are:

| Dataset | Prefix |
|---|---|
| Srinivas-110 | `srinivas` |
| Chandak-108 | `chandak108` |

Existing outputs are protected by default. Use a new output directory or set
`OVERWRITE=1` when replacement is intentional.

## Smoke test

Run the first 100 clusters of both datasets with two threads:

    rm -rf /tmp/ccpdp_difficult_smoke

    OUT_DIR=/tmp/ccpdp_difficult_smoke \
    DATASETS=both \
    LIMIT=100 \
    JOBS=2 \
      experiments/difficult_cluster_analysis/run_experiment.sh

Verify the number of reconstructed sequences:

    wc -l \
      /tmp/ccpdp_difficult_smoke/srinivas_final.txt \
      /tmp/ccpdp_difficult_smoke/chandak108_final.txt

Each file should contain 100 reconstructed sequences.

Verify their lengths:

    awk '
    length($0) != 110 {
        print "Invalid Srinivas length at line", NR
        bad = 1
    }
    END {
        if (!bad) print "[PASS] Srinivas lengths are 110"
    }' /tmp/ccpdp_difficult_smoke/srinivas_final.txt

    awk '
    length($0) != 108 {
        print "Invalid Chandak length at line", NR
        bad = 1
    }
    END {
        if (!bad) print "[PASS] Chandak lengths are 108"
    }' /tmp/ccpdp_difficult_smoke/chandak108_final.txt

A limited smoke run is useful for integration testing, but it is not expected
to reproduce the complete-dataset outputs. CCP-DP estimates dataset-level
statistics from the available cluster collection, so changing the collection
can change calibration and routing.

## Full reconstruction run

Run both complete datasets with 32 threads:

    OUT_DIR=experiments/difficult_cluster_analysis/results/reproduced \
    DATASETS=both \
    JOBS=32 \
      experiments/difficult_cluster_analysis/run_experiment.sh

If the output directory already contains previous results, either select
another directory or intentionally replace them with:

    OUT_DIR=experiments/difficult_cluster_analysis/results/reproduced \
    DATASETS=both \
    JOBS=32 \
    OVERWRITE=1 \
      experiments/difficult_cluster_analysis/run_experiment.sh

The `OVERWRITE=1` option should be used only when replacing the existing
reproduced records is intended.

## Run one dataset

Run only Srinivas-110:

    OUT_DIR=/tmp/ccpdp_difficult_srinivas \
    DATASETS=srinivas_110 \
    JOBS=32 \
      experiments/difficult_cluster_analysis/run_experiment.sh

Run only Chandak-108:

    OUT_DIR=/tmp/ccpdp_difficult_chandak \
    DATASETS=chandak_108 \
    JOBS=32 \
      experiments/difficult_cluster_analysis/run_experiment.sh

Accepted `DATASETS` values are:

- `both`;
- `srinivas_110`;
- `chandak_108`.

## Override input or executable locations

The entry point accepts environment-variable overrides:

- `CCPDP_BIN`: CCP-DP executable;
- `DATA_ROOT`: prepared dataset root;
- `OUT_DIR`: output directory;
- `JOBS`: OpenMP thread count;
- `DATASETS`: selected dataset;
- `LIMIT`: positive cluster limit;
- `OVERWRITE`: whether existing results may be replaced.

For example:

    CCPDP_BIN=/path/to/ccpdp \
    DATA_ROOT=/path/to/processed_data \
    OUT_DIR=/path/to/output \
    DATASETS=both \
    JOBS=32 \
      experiments/difficult_cluster_analysis/run_experiment.sh

The custom data root must retain the expected subdirectories
`srinivas_110/` and `chandak_108/`.

## Rebuild the mechanism table

Rebuild the analysis table from complete run records:

    python \
      experiments/difficult_cluster_analysis/analyze_mechanism.py \
      --runs-dir \
      experiments/difficult_cluster_analysis/results/reproduced \
      --srinivas-centers \
      data/processed/srinivas_110/Centers.txt \
      --chandak-centers \
      data/processed/chandak_108/Centers.txt \
      --output \
      experiments/difficult_cluster_analysis/results/ccpdp_analysis_data_rebuilt.csv \
      --bins 10

The default analysis divides clusters into ten routing-score bins for each
dataset.

## Compare rebuilt and preserved tables

The analysis program can compare a newly generated table against the preserved
source data:

    python \
      experiments/difficult_cluster_analysis/analyze_mechanism.py \
      --runs-dir \
      experiments/difficult_cluster_analysis/results/reproduced \
      --srinivas-centers \
      data/processed/srinivas_110/Centers.txt \
      --chandak-centers \
      data/processed/chandak_108/Centers.txt \
      --output /tmp/ccpdp_analysis_data_check.csv \
      --bins 10 \
      --compare \
      experiments/difficult_cluster_analysis/results/ccpdp_analysis_data.csv \
      --tolerance 1e-9

The comparison checks the rebuilt numerical table within the requested
floating-point tolerance.

## Analysis table structure

The preserved table contains two record types.

### Gate records

Gate rows summarize routing-score deciles and contain:

- dataset;
- decile order;
- number of clusters;
- guide failure percentage;
- final failure percentage;
- structured-decoding percentage;
- routing AUC;
- overall activated percentage.

These rows test whether difficult guide failures and structured decoding are
concentrated in high-risk score regions.

### Error-type records

Error-type rows group guide failures into:

- substitution-only errors;
- paired-indel errors;
- mixed errors.

They report:

- number and share of errors;
- exact rescue percentage;
- partial improvement percentage;
- unchanged percentage;
- worsened percentage;
- any-improvement percentage;
- edit-distance reduction percentage.

The reference centers are used only to label outcomes and compute evaluation
metrics. They are not supplied to CCP-DP reconstruction.

## Rebuild the mechanism figure

The plotting script uses fixed repository-relative locations:

Input:

    experiments/difficult_cluster_analysis/results/
    ccpdp_analysis_data.csv

Output prefix:

    experiments/difficult_cluster_analysis/figures/
    ccpdp_mechanism_analysis

Run:

    python \
      experiments/difficult_cluster_analysis/plot_ccpdp_mechanism.py

The generated files are:

- `ccpdp_mechanism_analysis.svg`;
- `ccpdp_mechanism_analysis.pdf`;
- `ccpdp_mechanism_analysis.png`;
- `ccpdp_mechanism_analysis.tiff`.

Running the plotting script rewrites the files under `figures/`.

## Non-destructive plot validation

To validate figure generation without modifying the preserved figures, copy
the plotting inputs to a temporary directory:

    rm -rf /tmp/ccpdp_mechanism_plot_check

    mkdir -p \
      /tmp/ccpdp_mechanism_plot_check/results \
      /tmp/ccpdp_mechanism_plot_check/figures

    cp \
      experiments/difficult_cluster_analysis/plot_ccpdp_mechanism.py \
      /tmp/ccpdp_mechanism_plot_check/

    cp \
      experiments/difficult_cluster_analysis/results/ccpdp_analysis_data.csv \
      /tmp/ccpdp_mechanism_plot_check/results/

    python \
      /tmp/ccpdp_mechanism_plot_check/plot_ccpdp_mechanism.py

Expected files:

    /tmp/ccpdp_mechanism_plot_check/figures/
    ccpdp_mechanism_analysis.svg

    /tmp/ccpdp_mechanism_plot_check/figures/
    ccpdp_mechanism_analysis.pdf

    /tmp/ccpdp_mechanism_plot_check/figures/
    ccpdp_mechanism_analysis.png

    /tmp/ccpdp_mechanism_plot_check/figures/
    ccpdp_mechanism_analysis.tiff

## Reproducibility notes

- Run commands from the repository root.
- Full-dataset and limited-run outputs are not expected to be identical.
- Dataset-level calibration depends on the cluster collection being processed.
- Frozen guides are produced before dataset-level contextual estimation.
- True centers are used only by the analysis program.
- Existing experiment outputs are protected unless `OVERWRITE=1` is set.
- Historical logs may retain their original machine paths.
- Temporary validation outputs should be written under `/tmp` or `build/`.
