# Coverage robustness evaluation

This experiment evaluates reconstruction accuracy as read coverage decreases.
CCP-DP and the classical baselines are evaluated on exactly the same nested
coverage samples.

## Experimental design

The default evaluation uses:

- read coverages: 5x, 10x, 15x, and 20x;
- random seeds: 2026, 2027, 2028, 2029, and 2030;
- three real DNA storage datasets;
- four reconstruction methods.

The evaluated methods are:

- CCP-DP;
- BBS;
- CPL;
- ITR.

For each dataset and seed, the lower-coverage inputs are nested subsets of the
higher-coverage inputs. This reduces variation caused by independently sampled
read sets and ensures that all reconstruction methods receive identical input
clusters.

The complete design contains:

    3 datasets x 4 coverages x 5 seeds = 60 sampled inputs

Each of the four methods is evaluated on those same 60 inputs.

## Datasets

| Dataset | Target length | Prepared input |
|---|---:|---|
| Srinivasavaradhan et al. | 110 bp | `data/processed/srinivas_110` |
| Chandak et al. | 108 bp | `data/processed/chandak_108` |
| Chandak et al. | 150 bp | `data/processed/chandak_150` |

Each prepared directory contains:

- `Clusters.txt`: noisy read clusters;
- `Centers.txt`: corresponding reference sequences.

Dataset provenance and preparation instructions are provided in
[`../../data/README.md`](../../data/README.md).

## Directory contents

The principal files are:

- `tools/make_nested_coverage.py`: nested coverage sampler;
- `tools/run_ccpdp_coverage.py`: CCP-DP runner and collector;
- `tools/run_bbs_coverage.py`: BBS runner and collector;
- `tools/run_cpl_coverage.py`: CPL runner and collector;
- `tools/run_itr_coverage.py`: ITR runner and collector;
- `tools/datasets_ccpdp.json`: CCP-DP configuration;
- `tools/datasets_bbs.json`: BBS configuration;
- `tools/datasets_cpl.json`: CPL configuration;
- `tools/datasets_itr.json`: ITR configuration;
- `plot_coverage_robustness_three_panel.py`: three-panel figure generator;
- `results/coverage_robustness_data.csv`: combined figure source data.

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

Build the classical baselines:

    ./baselines/build_baselines.sh

The expected executables are:

| Method | Executable |
|---|---|
| CCP-DP | `build/ccpdp` |
| BBS | `build/baselines/bbs-target/release/bbs` |
| CPL | `build/baselines/cpl32` |
| ITR | `build/baselines/itr_omp` |

Figure generation requires pandas and matplotlib.

## Validate commands without running experiments

The following dry run resolves all CCP-DP paths but does not generate samples
or execute reconstruction:

    python \
      experiments/coverage_robustness/tools/run_ccpdp_coverage.py \
      --config \
      experiments/coverage_robustness/tools/datasets_ccpdp.json \
      all \
      --datasets srinivas_110 chandak_108 chandak_150 \
      --coverages 5 \
      --seeds 2026 \
      --jobs 2 \
      --dry-run

The classical baseline commands can be checked similarly:

    for METHOD in bbs cpl itr
    do
        python \
          "experiments/coverage_robustness/tools/run_${METHOD}_coverage.py" \
          --config \
          "experiments/coverage_robustness/tools/datasets_${METHOD}.json" \
          run \
          --datasets srinivas_110 \
          --coverages 5 \
          --seeds 2026 \
          --threads 2 \
          --dry-run
    done

## Generate the shared nested samples

Generate all nested coverage inputs once with the CCP-DP runner:

    python \
      experiments/coverage_robustness/tools/run_ccpdp_coverage.py \
      --config \
      experiments/coverage_robustness/tools/datasets_ccpdp.json \
      sample \
      --datasets srinivas_110 chandak_108 chandak_150 \
      --coverages 5 10 15 20 \
      --seeds 2026 2027 2028 2029 2030

Generated samples are written under:

    data/generated/coverage_robustness

The resulting layout follows this pattern:

    data/generated/coverage_robustness/
    └── DATASET/
        └── seed_SEED/
            └── coverage_COVERAGEx/
                ├── Clusters.txt
                └── Centers.txt

Use `--force` only when the sampled inputs must be regenerated intentionally.

## Run CCP-DP

Run CCP-DP on all existing shared samples:

    python \
      experiments/coverage_robustness/tools/run_ccpdp_coverage.py \
      --config \
      experiments/coverage_robustness/tools/datasets_ccpdp.json \
      run \
      --datasets srinivas_110 chandak_108 chandak_150 \
      --coverages 5 10 15 20 \
      --seeds 2026 2027 2028 2029 2030 \
      --jobs 32 \
      --stop-on-error

Alternatively, generate the samples, run CCP-DP, and collect its result tables
in one command:

    python \
      experiments/coverage_robustness/tools/run_ccpdp_coverage.py \
      --config \
      experiments/coverage_robustness/tools/datasets_ccpdp.json \
      all \
      --datasets srinivas_110 chandak_108 chandak_150 \
      --coverages 5 10 15 20 \
      --seeds 2026 2027 2028 2029 2030 \
      --jobs 32 \
      --stop-on-error

The runner is intentionally sequential. Each CCP-DP process already uses
multiple OpenMP threads, so running several reconstructions concurrently would
oversubscribe the machine and make runtime comparisons unreliable.

Use `--rerun` only when completed runs should be replaced.

## Run the classical baselines

The following commands reuse the nested samples produced above.

Run BBS:

    python \
      experiments/coverage_robustness/tools/run_bbs_coverage.py \
      --config \
      experiments/coverage_robustness/tools/datasets_bbs.json \
      run \
      --datasets srinivas_110 chandak_108 chandak_150 \
      --coverages 5 10 15 20 \
      --seeds 2026 2027 2028 2029 2030 \
      --threads 32 \
      --stop-on-error

Run CPL:

    python \
      experiments/coverage_robustness/tools/run_cpl_coverage.py \
      --config \
      experiments/coverage_robustness/tools/datasets_cpl.json \
      run \
      --datasets srinivas_110 chandak_108 chandak_150 \
      --coverages 5 10 15 20 \
      --seeds 2026 2027 2028 2029 2030 \
      --threads 32 \
      --stop-on-error

Run ITR:

    python \
      experiments/coverage_robustness/tools/run_itr_coverage.py \
      --config \
      experiments/coverage_robustness/tools/datasets_itr.json \
      run \
      --datasets srinivas_110 chandak_108 chandak_150 \
      --coverages 5 10 15 20 \
      --seeds 2026 2027 2028 2029 2030 \
      --threads 32 \
      --stop-on-error

ITR additionally accepts `--max-reads INT` when its per-cluster read limit
needs to be overridden.

## Rebuild result tables

Existing per-run `metrics.json` files can be collected again without rerunning
the reconstruction programs.

For CCP-DP:

    python \
      experiments/coverage_robustness/tools/run_ccpdp_coverage.py \
      --config \
      experiments/coverage_robustness/tools/datasets_ccpdp.json \
      collect

For the classical baselines:

    for METHOD in bbs cpl itr
    do
        python \
          "experiments/coverage_robustness/tools/run_${METHOD}_coverage.py" \
          --config \
          "experiments/coverage_robustness/tools/datasets_${METHOD}.json" \
          collect
    done

The collectors aggregate exact reconstruction, reconstruction rate, edit
distance, runtime, and resource-use information where available.

## Generate the coverage figure

The plotting script reads:

    experiments/coverage_robustness/results/coverage_robustness_data.csv

This combined table must contain all four methods, all three datasets, and all
four coverage levels.

Generate figures in the default experiment directory with:

    python \
      experiments/coverage_robustness/plot_coverage_robustness_three_panel.py

To write validation figures somewhere else:

    rm -rf /tmp/ccpdp_coverage_plot_check
    mkdir -p /tmp/ccpdp_coverage_plot_check

    COVERAGE_FIG_OUTPUT_DIR=/tmp/ccpdp_coverage_plot_check \
      python \
      experiments/coverage_robustness/plot_coverage_robustness_three_panel.py

The plotting script produces:

- `coverage_robustness_three_panel.svg`;
- `coverage_robustness_three_panel.pdf`;
- `coverage_robustness_three_panel.png`;
- `coverage_robustness_three_panel.tiff`.

The solid curves show exact reconstruction percentages. The dashed curves show
mean edit distance on a logarithmic scale.

## Output interpretation

The principal reported quantities are:

- exact reconstruction percentage;
- standard deviation across five seeds;
- mean edit distance;
- reconstruction rate;
- runtime;
- peak memory use where recorded.

Accuracy comparisons are valid only when methods use the same sampled
`Clusters.txt` and corresponding `Centers.txt`.

## Reproducibility notes

- Dataset, executable, sample, and result paths are repository-relative.
- Every method uses the same saved nested samples.
- Lower coverage inputs are nested within higher coverage inputs for each seed.
- Sampling seeds must be retained when regenerating results.
- The full benchmark should be run sequentially on an otherwise lightly loaded
  machine when runtime is being compared.
- `--dry-run` validates paths and commands but performs no reconstruction.
- `--force` regenerates sampled data and should be used cautiously.
- `--rerun` replaces completed reconstruction runs.
- Temporary validation figures should be written under `/tmp` or `build/`.
