# CCP-DP parameter sensitivity

This experiment evaluates whether the main conclusions of CCP-DP remain
stable when its routing and evidence-allocation parameters are varied one
family at a time.

The sensitivity executable is separate from the production executable.
Its default settings reproduce the production CCP-DP behavior, while
additional command-line options expose selected parameters for controlled
experiments.

## Experimental principle

The experiment follows a one-at-a-time design:

- the baseline uses all production defaults;
- every non-baseline configuration changes exactly one parameter family;
- all other parameters remain fixed;
- every configuration is evaluated on the same datasets;
- full experiments use three repeats by default.

This design isolates the effect of each parameter without changing the
underlying reconstruction algorithm.

## Evaluated datasets

The default configuration is stored in `datasets_server.json`.

| Dataset | Length | Clusters | Centers |
|---|---:|---|---|
| Srinivas-110 | 110 bp | `data/processed/srinivas_110/Clusters.txt` | `data/processed/srinivas_110/Centers.txt` |
| Chandak-108 | 108 bp | `data/processed/chandak_108/Clusters.txt` | `data/processed/chandak_108/Centers.txt` |
| Erlich-152 | 152 bp | `data/processed/erlich_152/Clusters.txt` | `data/processed/erlich_152/Centers.txt` |

The standard `full` mode evaluates Srinivas-110 and Chandak-108.
Use `full-all` to include Erlich-152.

## Parameter families

The experiment varies the following parameter families.

| Parameter family | Production default | Evaluated values |
|---|---:|---|
| High/low coverage switch | 40 | 20, 30, 40, 50, 60 |
| Cross-cluster reuse threshold | 0.50 | 0.30, 0.40, 0.50, 0.60, 0.70 |
| Phase-load thresholds | 3.0 / 6.0 | 2.5 / 5.0, 3.0 / 6.0, 3.5 / 7.0 |
| Cross-evidence weight scale | 1.00 | 0.75, 1.00, 1.25 |
| Minimum routed fraction | 0.10 | 0.05, 0.10, 0.15, 0.20 |
| Phase band | automatic/default | 4, 6, 8 |

The configured grid contains 18 configurations including the baseline.

The maximum routed fraction remains fixed at 0.50 unless explicitly changed
for a separate diagnostic experiment.

## Files

- `ccpdp_sensitivity.hpp`: parameter-exposed CCP-DP implementation
- `main_sensitivity.cpp`: sensitivity executable entry point
- `datasets_server.json`: portable dataset configuration
- `run_parameter_sensitivity.py`: experiment runner and metric collector
- `run_sensitivity.sh`: build-independent experiment entry point
- `plot_parameter_sensitivity.py`: exploratory sensitivity plots
- `make_parameter_sensitivity_paper_outputs.py`: final figures, tables, and source data
- `requirements.txt`: Python analysis and plotting dependencies
- `sensitivity_results/`: preserved full-run result tables
- `smoke_results/`: preserved integration-test records
- `final_outputs/`: preserved paper figures, tables, and source data
- `VALIDATION_REPORT.md`: validation record for the sensitivity implementation

## Prerequisites

Run commands from the repository root.

Prepare the datasets:

    ./data/download_datasets.sh --core

Install the Python dependencies:

    python -m pip install \
      -r experiments/parameter_sensitivity/requirements.txt

## Build the sensitivity executable

Compile the parameter-sensitivity executable with:

    mkdir -p build

    g++ \
      -std=c++17 \
      -O3 \
      -fopenmp \
      experiments/parameter_sensitivity/main_sensitivity.cpp \
      -o build/ccpdp_sensitivity

Verify the executable:

    ./build/ccpdp_sensitivity --help

The resulting executable is written to:

    build/ccpdp_sensitivity

It does not replace the production executable at `build/ccpdp`.

## Available experiment modes

The shell entry point accepts four modes.

| Mode | Datasets | Repeats | Purpose |
|---|---|---:|---|
| `smoke` | Srinivas-110, Chandak-108 | 1 | Small integration test |
| `smoke-all` | All configured datasets | 1 | Integration test including Erlich-152 |
| `full` | Srinivas-110, Chandak-108 | 3 | Main sensitivity experiment |
| `full-all` | All configured datasets | 3 | Extended sensitivity experiment |

## Smoke test

Run a small validation with:

    JOBS=2 \
    experiments/parameter_sensitivity/run_sensitivity.sh smoke

The default smoke output directory is:

    experiments/parameter_sensitivity/results/smoke

To keep validation output outside the repository:

    OUTPUT_DIR=/tmp/ccpdp_parameter_sensitivity_smoke \
    JOBS=2 \
    experiments/parameter_sensitivity/run_sensitivity.sh smoke

The smoke mode verifies:

- executable invocation;
- dataset path resolution;
- parameter forwarding;
- prediction generation;
- diagnostic output generation;
- metric collection;
- aggregate CSV generation.

Smoke-test accuracy is not a reported benchmark result.

## Full experiment

Run the main two-dataset experiment with:

    JOBS=32 \
    experiments/parameter_sensitivity/run_sensitivity.sh full

The default output directory is:

    experiments/parameter_sensitivity/results/reproduced

The experiment contains:

    2 datasets x 18 configurations x 3 repeats = 108 runs

Run the extended three-dataset experiment with:

    JOBS=32 \
    experiments/parameter_sensitivity/run_sensitivity.sh full-all

The extended experiment contains:

    3 datasets x 18 configurations x 3 repeats = 162 runs

A full run may require substantial CPU time. Runs are executed sequentially
to avoid oversubscribing the machine and distorting runtime measurements.

## Custom output directory

Override the output location with:

    OUTPUT_DIR=/path/to/output \
    JOBS=32 \
    experiments/parameter_sensitivity/run_sensitivity.sh full

Existing completed runs are reused by default. Use the runner's `--force`
option only when intentionally replacing prior outputs.

## Direct runner usage

The Python runner can be invoked directly:

    python \
      experiments/parameter_sensitivity/run_parameter_sensitivity.py \
      --config \
      experiments/parameter_sensitivity/datasets_server.json \
      --output /path/to/output \
      --binary build/ccpdp_sensitivity \
      --datasets Srinivas-110 Chandak-108 \
      --jobs 32 \
      --repeats 3

Useful optional arguments include:

- `--limit N`: process only the first N clusters;
- `--smoke`: use the smoke-test configuration;
- `--force`: replace existing completed runs;
- `--datasets`: restrict the evaluated datasets.

Use `--limit` only for debugging. Dataset-wide calibration and routing depend
on the complete collection of clusters, so limited runs are not directly
comparable with full experiments.

## Per-run outputs

Each dataset, configuration, and repeat records:

- `result.txt`: reconstructed sequences;
- `diag.csv`: per-cluster routing and diagnostic information;
- `run.log`: program output and runtime information;
- `metrics.json`: reconstruction metrics and resolved command.

## Aggregated outputs

The runner produces:

- `raw_metrics.csv`: one row per dataset, configuration, and repeat;
- `summary_metrics.csv`: metrics aggregated across repeats;
- `run_manifest.json`: run configuration and provenance information.

The primary accuracy measurements are:

- exact reconstruction count;
- exact success rate;
- reconstruction rate;
- mean edit distance.

The experiment also records runtime and routing behavior.

## Exploratory plots

Generate the standard sensitivity plots from a completed summary table:

    python \
      experiments/parameter_sensitivity/plot_parameter_sensitivity.py \
      /path/to/summary_metrics.csv \
      --output-dir /path/to/figures

The plotting program produces success-rate and routing figures in PNG, SVG,
and PDF formats.

## Rebuild the paper outputs

The repository includes preserved full-run summary data and final source data.

Rebuild the paper-ready outputs with:

    rm -rf /tmp/ccpdp_parameter_sensitivity_outputs

    python \
      experiments/parameter_sensitivity/make_parameter_sensitivity_paper_outputs.py \
      experiments/parameter_sensitivity/sensitivity_results/summary_metrics.csv \
      --output /tmp/ccpdp_parameter_sensitivity_outputs

The paper-output program requires at least three repeats per configuration.
The option `--allow-single-repeat` is intended only for development and must
not be used for reported results.

Generated outputs include:

- `Fig_parameter_sensitivity_main.{png,svg,pdf}`;
- `Fig_S_success_all_parameters.{png,svg,pdf}`;
- `Fig_S_mean_ed_all_parameters.{png,svg,pdf}`;
- `Fig_S_runtime_all_parameters.{png,svg,pdf}`;
- `Fig_S_routing_all_parameters.{png,svg,pdf}`;
- `source_data_main_figure.csv`;
- `source_data_all_parameters.csv`;
- `default_results.csv`;
- `parameter_robustness_summary.csv`;
- `table_default_results.tex`;
- `table_full_sensitivity.tex`;
- `results_summary.txt`;
- `output_manifest.json`.

The preserved versions are stored under:

    experiments/parameter_sensitivity/final_outputs

## Interpretation

The parameter-sensitivity analysis is intended to determine whether the
reported CCP-DP behavior depends on a narrowly tuned setting.

A robust parameter family should show:

- stable exact reconstruction performance near its default;
- no abrupt degradation under moderate perturbations;
- interpretable changes in the structured-decoding fraction;
- runtime changes consistent with the amount of routed work.

Accuracy, routing, and runtime should be interpreted together.

## Reproducibility notes

- Centers are used only for evaluation, not reconstruction.
- Every non-baseline configuration changes one parameter family.
- Full reported results use three repeats.
- Use the same compiler, thread count, and machine-load conditions for runtime comparisons.
- Active scripts and dataset configurations use repository-relative paths.
- Preserved manifests may contain absolute paths from the original machine.
  These paths are historical provenance only and are not used by active scripts.
- Do not rewrite preserved historical manifests solely to replace their recorded paths.
- Generated smoke outputs should be written under `/tmp` or another ignored directory.
