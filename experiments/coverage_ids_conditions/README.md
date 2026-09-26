# Controlled coverage-by-IDS evaluation

This experiment evaluates CCP-DP under controlled combinations of read
coverage and insertion, deletion, and substitution error rates.

Unlike the coverage-robustness experiment, which downsamples real read
clusters, this experiment generates synthetic noisy read clusters from the
Srinivas 110-bp reference sequences. This makes it possible to vary coverage
and IDS-channel difficulty independently.

## Experimental question

The experiment tests whether CCP-DP provides the largest benefit when local
cluster evidence is limited or phase errors become difficult.

It compares:

- CCP-DP;
- BBS;
- CPL.

The experiment also records the fraction of clusters routed to structured
CCP-DP decoding. This connects reconstruction performance with the amount of
additional computation selected by the adaptive routing mechanism.

## Source sequences

Synthetic read clusters are generated from:

    data/processed/srinivas_110/Centers.txt

The source contains 9,984 reference sequences of length 110 bp.

Dataset provenance is documented in
[`../../data/README.md`](../../data/README.md).

## Error model

The generator controls:

- read coverage;
- total insertion-plus-deletion percentage;
- substitution percentage;
- random seed.

The substitution rate is fixed at 2% in the reported experiment.

The requested total indel percentage is divided equally between insertion and
deletion probabilities. For example, a total indel level of 6% corresponds to
3% insertion and 3% deletion.

Generated metadata record the exact coverage, error probabilities, seed,
sequence length, and generator version for every condition.

## Experimental grids

### Base grid

The base experiment uses:

- coverage: 5x, 10x, 15x, and 20x;
- total indel percentage: 2%, 4%, 6%, 8%, and 10%;
- substitution percentage: 2%;
- seeds: 2026, 2027, 2028, 2029, and 2030.

This produces:

    4 coverages x 5 indel levels x 5 seeds = 100 conditions

Three reconstruction methods are evaluated for every condition:

    100 conditions x 3 methods = 300 method runs

### Dense refinement

The dense refinement adds:

- 6x;
- 7x;
- 8x;
- 9x;
- 12x.

This produces an additional:

    5 coverages x 5 indel levels x 5 seeds = 125 conditions

and:

    125 conditions x 3 methods = 375 method runs

After merging the base and refined grids, the complete experiment contains:

    9 coverages x 5 indel levels x 5 seeds = 225 conditions

and:

    225 conditions x 3 methods = 675 method runs

## Directory contents

The main files are:

- `generate_ids_stress.py`: paired synthetic IDS-data generator;
- `run_ids_stress.py`: multi-method experiment runner;
- `experiment_config.srinivas.json`: method commands and thread settings;
- `run_srinivas_stress.sh`: base-grid entry point;
- `run_srinivas_stress_dense.sh`: dense-grid entry point;
- `merge_stress_metrics.py`: base and refined result merger;
- `extract_routing_summary.py`: CCP-DP routing-data collector;
- `plot_ids_stress.py`: base-grid plotting program;
- `plot_ids_stress_dense.py`: dense-grid plotting program;
- `requirements.txt`: Python dependencies;
- `results/`: preserved source data for the reported dense experiment.

Run all commands from the repository root.

## Prerequisites

Prepare the Srinivas dataset:

    ./data/download_datasets.sh --core

Build the classical baselines:

    ./baselines/build_baselines.sh

Install the Python dependencies:

    python -m pip install \
      -r experiments/coverage_ids_conditions/requirements.txt

The expected executables are:

| Method | Executable |
|---|---|
| CCP-DP | `build/ccpdp` |
| BBS | `build/baselines/bbs-target/release/bbs` |
| CPL | `build/baselines/cpl32` |

The base entry point compiles CCP-DP automatically when the executable is
missing or older than its source files.

## Validate resolved commands

Generate one 50-cluster condition and print the three resolved method commands
without running reconstruction:

    experiments/coverage_ids_conditions/run_srinivas_stress.sh \
      commands

This creates a small synthetic input under `smoke_data/`, then invokes the
runner with `--dry-run`.

For the 6x dense-grid condition, use:

    experiments/coverage_ids_conditions/run_srinivas_stress_dense.sh \
      commands

Dry runs validate executable, dataset, separator, output, and diagnostic
paths without executing the reconstruction methods.

## Smoke test

Run one 50-cluster condition with:

- coverage: 5x;
- total indel rate: 6%;
- substitution rate: 2%;
- seed: 2026.

Command:

    experiments/coverage_ids_conditions/run_srinivas_stress.sh \
      smoke

The three configured methods are executed. Smoke outputs are written under:

    experiments/coverage_ids_conditions/smoke_results

The dense entry point provides a corresponding 6x smoke test:

    experiments/coverage_ids_conditions/run_srinivas_stress_dense.sh \
      smoke

Its temporary outputs are written under:

    experiments/coverage_ids_conditions/smoke_results_dense

Smoke-test values validate integration only and must not be reported as
complete experimental results.

## Run the base grid

Generate the 100 paired conditions without running reconstruction:

    experiments/coverage_ids_conditions/run_srinivas_stress.sh \
      generate

The generated conditions are written under:

    experiments/coverage_ids_conditions/stress_data

Run the complete 300-evaluation base experiment:

    experiments/coverage_ids_conditions/run_srinivas_stress.sh \
      full

This command:

1. compiles CCP-DP when required;
2. generates or verifies the paired synthetic data;
3. runs CCP-DP, BBS, and CPL;
4. collects the per-run metrics.

The base results are written under:

    experiments/coverage_ids_conditions/stress_results

To continue an interrupted run without regenerating the synthetic data:

    experiments/coverage_ids_conditions/run_srinivas_stress.sh \
      resume

Check progress with:

    experiments/coverage_ids_conditions/run_srinivas_stress.sh \
      status

A complete base experiment should report:

    Conditions: 100/100
    Metric rows: 300/300

## Plot the base grid

After the base experiment has completed:

    experiments/coverage_ids_conditions/run_srinivas_stress.sh \
      plot

The default figure directory is:

    experiments/coverage_ids_conditions/stress_figure

The plotting program produces:

- `coverage_ids_stress.png`;
- `coverage_ids_stress.svg`;
- `coverage_ids_stress.pdf`;
- `stress_summary.csv`;
- `method_summary.csv`;
- `ccpdp_routing_summary.csv`.

The `--show-sd` option is enabled by the shell entry point.

## Run the dense refinement

Run only the additional 6x, 7x, 8x, 9x, and 12x conditions:

    experiments/coverage_ids_conditions/run_srinivas_stress_dense.sh \
      refine

This creates:

- 125 additional conditions;
- 375 additional method runs.

The generated inputs and results are written under:

    experiments/coverage_ids_conditions/stress_data_refined
    experiments/coverage_ids_conditions/stress_results_refined

## Merge a newly reproduced base grid

When the base grid was reproduced under `stress_results/`, merge it with the
new refined results by explicitly selecting its metric table:

    BASE_METRICS=experiments/coverage_ids_conditions/stress_results/raw_metrics.csv \
      experiments/coverage_ids_conditions/run_srinivas_stress_dense.sh \
      merge

The merged files are written under:

    experiments/coverage_ids_conditions/stress_results_dense

The merge program uses `(condition_id, method)` as its unique key. Refined rows
replace matching base rows if a condition is repeated accidentally.

## Reproduce the complete dense grid

After completing the base experiment, generate the refinement, merge the
results, and create the dense figure with:

    BASE_METRICS=experiments/coverage_ids_conditions/stress_results/raw_metrics.csv \
      experiments/coverage_ids_conditions/run_srinivas_stress_dense.sh \
      all

This is the recommended from-scratch workflow.

The dense shell script also has a default preserved base table. Explicitly
setting `BASE_METRICS` as shown above makes the provenance of a newly
reproduced run unambiguous.

## Rebuild the reported dense figure

The repository contains preserved dense-grid source data under:

    experiments/coverage_ids_conditions/results

The preserved files include:

- `raw_metrics.csv`: 675 method rows;
- `routing_raw.csv`: collected per-run routing measurements;
- `routing_summary.csv`: aggregated routing fractions;
- `method_summary_dense.csv`: aggregated method metrics;
- `advantage_summary_dense.csv`: CCP-DP advantages over the baselines.

Rebuild the reported figure in a temporary directory with:

    rm -rf /tmp/ccpdp_ids_plot_check
    mkdir -p /tmp/ccpdp_ids_plot_check

    python \
      experiments/coverage_ids_conditions/plot_ids_stress_dense.py \
      --metrics \
      experiments/coverage_ids_conditions/results/raw_metrics.csv \
      --routing \
      experiments/coverage_ids_conditions/results/routing_summary.csv \
      --outdir \
      /tmp/ccpdp_ids_plot_check

The dense plotting program produces:

- `coverage_ids_dense_main.png`;
- `coverage_ids_dense_main.svg`;
- `coverage_ids_dense_main.pdf`;
- `method_summary_dense.csv`;
- `advantage_summary_dense.csv`.

## Direct generator usage

The data generator may be called directly:

    python \
      experiments/coverage_ids_conditions/generate_ids_stress.py \
      --centers data/processed/srinivas_110/Centers.txt \
      --outdir /tmp/ccpdp_ids_data \
      --coverage 5,10,15,20 \
      --indel-pct 2,4,6,8,10 \
      --sub-pct 2 \
      --seeds 2026,2027,2028,2029,2030

Use `--limit-centers N` only for debugging or smoke testing.

## Direct runner usage

Run selected generated conditions directly with:

    python \
      experiments/coverage_ids_conditions/run_ids_stress.py \
      --config \
      experiments/coverage_ids_conditions/experiment_config.srinivas.json \
      --data-root /path/to/generated_conditions \
      --results-root /path/to/results

Useful optional arguments include:

- `--methods`: restrict the reconstruction methods;
- `--max-conditions`: restrict the number of generated conditions;
- `--dry-run`: print commands without execution;
- `--force`: replace existing completed runs.

## Recorded metrics

The combined result table records:

- condition identifier;
- coverage;
- total indel percentage;
- substitution, insertion, and deletion percentages;
- seed;
- reconstruction method;
- number of clusters;
- exact reconstruction count;
- exact success percentage;
- reconstruction percentage;
- mean edit distance;
- runtime;
- maximum resident memory;
- output path.

The primary comparison uses exact reconstruction percentage. Runtime and
routing fraction explain the computational behavior of the adaptive method.

## Reproducibility notes

- All three methods receive the same synthetic clusters for each condition.
- The same five seeds are used at every coverage and error level.
- The paired generator version is checked before a full base run.
- Generated synthetic datasets and smoke outputs are not reported results.
- Use the same thread count for runtime comparisons.
- Run the methods sequentially on an otherwise lightly loaded machine.
- Do not combine results generated with different error models.
- The historical `output_path` column may record its original run location;
  active scripts do not depend on that stored path.
