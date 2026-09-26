# Cross-cluster context and adaptive routing analysis

This experiment studies two questions about CCP-DP:

1. Does cross-cluster context provide meaningful reconstruction evidence?
2. Does adaptive selective routing retain decode-all accuracy while reducing
   computation?

The directory contains a counterfactual context experiment and an adaptive
evidence-allocation experiment. Their preserved source data are combined into
one multi-panel figure.

## Analysis overview

### Counterfactual context validation

The counterfactual experiment isolates the effect of cross-cluster context by
using decode-all in every condition. This prevents routing decisions from
confounding the context comparison.

It evaluates four variants:

| Variant | Description |
|---|---|
| Local only | Structured decoding without cross-cluster contextual evidence |
| Correct context | Context obtained from the correct frozen guide collection |
| Table-permuted | Context-table assignments are permuted as a negative control |
| Wrong-dataset | Context is taken from the other dataset |

Correct-context runs are executed first because their frozen guides are needed
by the wrong-dataset negative control.

The table-permuted condition is repeated five times by default. The other
conditions use one run because they are deterministic for fixed inputs.

### Adaptive evidence allocation

The adaptive-routing experiment compares:

- selective CCP-DP routing;
- decode-all structured reconstruction.

It records:

- local-guide accuracy;
- final selective accuracy;
- decode-all accuracy;
- structured-decoding fraction;
- changed-cluster count;
- cross-cluster reuse;
- phase load;
- runtime and speedup.

This analysis tests whether CCP-DP concentrates structured decoding on
difficult clusters instead of applying the expensive procedure uniformly.

## Datasets

The counterfactual experiment uses:

| Dataset | Length | Prepared directory |
|---|---:|---|
| Srinivasavaradhan et al. | 110 bp | `data/processed/srinivas_110` |
| Chandak et al. | 108 bp | `data/processed/chandak_108` |

The full adaptive-routing experiment additionally uses:

| Dataset | Length | Prepared directory |
|---|---:|---|
| Chandak et al. | 150 bp | `data/processed/chandak_150` |
| Erlich and Zielinski | 152 bp | `data/processed/erlich_152` |

An optional 110-bp `fasta_2` dataset can be included by setting both
`FASTA2_CLUSTERS` and `FASTA2_CENTERS`.

Dataset provenance and preparation are documented in
[`../../data/README.md`](../../data/README.md).

## Directory contents

The principal files are:

- `run_counterfactual.py`: counterfactual context runner;
- `run_adaptive_routing.py`: selective-versus-decode-all runner;
- `run_adaptive_all.sh`: adaptive experiment entry point;
- `plot_counterfactual.py`: counterfactual figure generator;
- `plot_adaptive.py`: adaptive-allocation figure generator;
- `plot_combined.py`: combined paper figure generator;
- `requirements.txt`: Python requirements;
- `results/counterfactual_summary.csv`: preserved counterfactual metrics;
- `results/local_vs_correct_pairwise.csv`: paired local-versus-correct counts;
- `results/adaptive_evidence_allocation_data.csv`: preserved routing summary;
- `figures/`: preserved combined figures.

Run all commands from the repository root.

## Prerequisites

Prepare the required datasets:

    ./data/download_datasets.sh --core

Install the plotting dependencies:

    python -m pip install \
      -r experiments/cross_cluster_context/requirements.txt

Build the main CCP-DP executable:

    mkdir -p build

    g++ \
      -std=c++17 \
      -O3 \
      -fopenmp \
      -Iinclude \
      src/main.cpp \
      -o build/ccpdp

Build the counterfactual executable:

    g++ \
      -std=c++17 \
      -O3 \
      -fopenmp \
      -Iinclude \
      src/counterfactual.cpp \
      -o build/ccpdp_counterfactual

Verify both executables:

    ./build/ccpdp --help

    ./build/ccpdp_counterfactual --help

## Counterfactual smoke test

The counterfactual smoke mode processes the first 200 clusters of the
Srinivas-110 and Chandak-108 datasets.

Run it with:

    rm -rf /tmp/ccpdp_counterfactual_smoke

    python \
      experiments/cross_cluster_context/run_counterfactual.py \
      smoke \
      --project-root . \
      --binary build/ccpdp_counterfactual \
      --outdir /tmp/ccpdp_counterfactual_smoke \
      --jobs 2 \
      --permuted-repeats 2

Expected summary files are:

    /tmp/ccpdp_counterfactual_smoke/counterfactual_summary_smoke.csv
    /tmp/ccpdp_counterfactual_smoke/local_vs_correct_pairwise_smoke.csv

The smoke test validates execution and output structure. Its accuracy values
must not be reported as complete results.

## Full counterfactual experiment

Run the two complete datasets with five table-permutation repetitions:

    python \
      experiments/cross_cluster_context/run_counterfactual.py \
      full \
      --project-root . \
      --binary build/ccpdp_counterfactual \
      --outdir \
      experiments/cross_cluster_context/results/reproduced_counterfactual \
      --jobs 32 \
      --permuted-repeats 5

The output directory contains reconstructed sequences, frozen guides,
diagnostic CSV files, logs, and the two aggregate tables:

- `counterfactual_summary.csv`;
- `local_vs_correct_pairwise.csv`.

The summary table reports exact success, reconstruction percentage, mean edit
distance, and standard deviations across repetitions.

The pairwise table compares local-only and correct-context outputs using:

- clusters reconstructed exactly by both;
- clusters harmed by correct context;
- clusters rescued by correct context;
- clusters reconstructed incorrectly by both;
- net exact gain;
- exact two-sided McNemar p-value.

## Adaptive-routing smoke test

The adaptive smoke mode uses:

- Srinivas-110;
- Chandak-108;
- one repeat;
- the complete selected datasets.

Run:

    experiments/cross_cluster_context/run_adaptive_all.sh \
      smoke

Outputs are written under:

    experiments/cross_cluster_context/results/adaptive_smoke

The smoke summary is:

    experiments/cross_cluster_context/results/
    adaptive_evidence_allocation_data_smoke.csv

This mode is smaller than the full adaptive experiment because it uses two
datasets and one repeat, but it still processes the complete selected
datasets.

## Full adaptive-routing experiment

Run:

    experiments/cross_cluster_context/run_adaptive_all.sh \
      full

The default full experiment uses:

- Srinivas-110;
- Chandak-108;
- Chandak-150;
- Erlich-152;
- three repeats per dataset.

Per-run outputs are written under:

    experiments/cross_cluster_context/results/adaptive_runs

The aggregated summary is:

    experiments/cross_cluster_context/results/
    adaptive_evidence_allocation_data.csv

Each reconstruction process uses 32 OpenMP threads.

## Optional additional dataset

To include an additional 110-bp dataset:

    FASTA2_CLUSTERS=/path/to/Clusters.txt \
    FASTA2_CENTERS=/path/to/Centers.txt \
      experiments/cross_cluster_context/run_adaptive_all.sh \
      full

Both environment variables must be supplied together. When neither is set,
the optional dataset is skipped.

Results obtained from the optional dataset should only be reported when its
provenance and preparation procedure are documented separately.

## Direct adaptive-runner usage

A single dataset can be evaluated directly:

    python \
      experiments/cross_cluster_context/run_adaptive_routing.py \
      --binary build/ccpdp \
      --dataset Srinivas-110 \
      --clusters data/processed/srinivas_110/Clusters.txt \
      --centers data/processed/srinivas_110/Centers.txt \
      --length 110 \
      --separator "===============================" \
      --jobs 32 \
      --repeats 3 \
      --outdir /tmp/ccpdp_adaptive_srinivas \
      --summary-csv /tmp/ccpdp_adaptive_summary.csv

The runner compares adaptive selective routing with decode-all using the same
dataset and reference centers.

## Plot the counterfactual analysis

Generate the standalone counterfactual figure from preserved results:

    rm -rf /tmp/ccpdp_counterfactual_plot
    mkdir -p /tmp/ccpdp_counterfactual_plot

    python \
      experiments/cross_cluster_context/plot_counterfactual.py \
      --summary \
      experiments/cross_cluster_context/results/counterfactual_summary.csv \
      --pairwise \
      experiments/cross_cluster_context/results/local_vs_correct_pairwise.csv \
      --output-dir /tmp/ccpdp_counterfactual_plot

The output stem is:

    counterfactual_validation_macaron

PNG, SVG, and PDF files are produced. The filename is retained for
compatibility with the existing plotting workflow and does not affect the
experiment or numerical results.

## Plot the adaptive-allocation analysis

Generate the standalone adaptive-allocation figure:

    rm -rf /tmp/ccpdp_adaptive_plot
    mkdir -p /tmp/ccpdp_adaptive_plot

    python \
      experiments/cross_cluster_context/plot_adaptive.py \
      --data \
      experiments/cross_cluster_context/results/adaptive_evidence_allocation_data.csv \
      --output-dir /tmp/ccpdp_adaptive_plot

The plotting program produces:

- `adaptive_evidence_allocation_diverging.png`;
- `adaptive_evidence_allocation_diverging.svg`;
- `adaptive_evidence_allocation_diverging.pdf`.

## Rebuild the combined figure

The combined plotting program uses the three preserved CSV files by default:

    python \
      experiments/cross_cluster_context/plot_combined.py

Default outputs are written under:

    experiments/cross_cluster_context/figures

The generated files are:

- `combined_counterfactual_adaptive.png`;
- `combined_counterfactual_adaptive.svg`;
- `combined_counterfactual_adaptive.pdf`.

For a validation build that does not alter preserved figures:

    rm -rf /tmp/ccpdp_cross_plot_check
    mkdir -p /tmp/ccpdp_cross_plot_check

    python \
      experiments/cross_cluster_context/plot_combined.py \
      --output-dir /tmp/ccpdp_cross_plot_check

## Interpretation

The counterfactual experiment should be interpreted through both aggregate
accuracy and paired cluster outcomes.

A useful contextual signal should satisfy the following expectations:

- correct context improves over local-only decoding;
- permuted context performs worse than correctly aligned context;
- wrong-dataset context does not reproduce the correct-context gain;
- rescued clusters outnumber harmed clusters.

The adaptive-routing analysis asks a different question. Selective decoding
should remain close to decode-all accuracy while using structured decoding on
a smaller fraction of clusters and reducing runtime.

Together, the two experiments show whether contextual evidence is informative
and whether CCP-DP allocates expensive decoding efficiently.

## Reproducibility notes

- Frozen guides are generated without using true centers.
- Centers are used only for evaluation.
- Every counterfactual condition uses decode-all to isolate contextual effects.
- Table permutations use repeated negative-control runs.
- Correct-context runs precede wrong-dataset runs because the latter require
  frozen guides from the opposite dataset.
- Active scripts use repository-relative dataset and executable paths.
- Historical logs may contain their original machine paths.
- Runtime comparisons should use the same thread count and an otherwise
  lightly loaded machine.
- Temporary validation outputs should be written under `/tmp` or `build/`.
