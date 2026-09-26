# CCP-DP

Cross-Cluster Context-Guided Phase Dynamic Programming for DNA storage
sequence reconstruction.

CCP-DP reconstructs a fixed-length DNA sequence from each cluster of noisy
reads affected by insertion, deletion, and substitution errors. The method
combines local fixed-length guides, dataset-level cross-cluster evidence,
blind IDS-channel calibration, and selective structured decoding.

## Method overview

CCP-DP follows five stages:

1. Construct one fixed-length local guide for each read cluster.
2. Freeze all guides before estimating dataset-level statistics.
3. Estimate cross-cluster contextual evidence and the blind IDS channel.
4. Route only difficult clusters to structured phase-aware decoding.
5. Apply a final acceptance gate before replacing the local guide.

Every cluster is decoded at most once.

## Repository layout

    CCP-DP/
    ├── include/          CCP-DP implementation
    ├── src/              Main command-line programs
    ├── scripts/          Shared evaluation utilities
    ├── data/             Dataset download and preparation tools
    ├── baselines/        Classical and deep-learning baselines
    ├── experiments/      Reproduction scripts and saved results
    ├── build/            Local binaries, third-party code and weights
    └── INSTALL.md        Installation instructions

The `build/` directory, downloaded datasets, model weights, and local
smoke-test outputs are excluded from Git.

## Requirements

The main CCP-DP program requires:

- Linux x86-64
- GCC/G++ with C++17 support
- OpenMP
- Python 3 for evaluation and experiment scripts

Plotting scripts additionally use packages such as NumPy, pandas, matplotlib,
and SciPy. Deep-learning baselines use separate environments.

See [INSTALL.md](INSTALL.md) for installation details.

## Build CCP-DP

Run from the repository root:

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

Build the counterfactual executable with:

    g++ \
      -std=c++17 \
      -O3 \
      -fopenmp \
      -Iinclude \
      src/counterfactual.cpp \
      -o build/ccpdp_counterfactual

## Obtain the datasets

Dataset files are not committed to this repository.

Download and prepare the primary public benchmarks with:

    ./data/download_datasets.sh --core

Prepared data are placed under `data/processed/`.

See [data/README.md](data/README.md) for dataset provenance, citations,
download instructions, preprocessing, and expected checksums.

## Quick start

Reconstruct the Srinivas 110-bp dataset:

    mkdir -p build/quickstart

    ./build/ccpdp \
      -i data/processed/srinivas_110/Clusters.txt \
      -l 110 \
      -s "====" \
      -o build/quickstart/result.txt \
      --diag build/quickstart/diag.csv \
      --mode auto \
      --jobs 32

Evaluate the reconstruction:

    python scripts/evaluate_answers_extended.py \
      -o build/quickstart/result.txt \
      -a data/processed/srinivas_110/Centers.txt

The output contains one reconstructed sequence per line in the same cluster
order as the input.

## Main command-line options

- `-i PATH`: clustered reads
- `-l INT`: target sequence length
- `-s STRING`: cluster separator
- `-o PATH`: reconstructed sequences
- `--diag PATH`: per-cluster diagnostic table
- `--mode auto`: automatic calibration and routing
- `--jobs INT`: OpenMP thread count
- `--limit INT`: process only the first N clusters
- `--guide-output PATH`: write frozen local guides

Run `./build/ccpdp --help` for the complete option list.

## Datasets

| Dataset | Target length | Prepared directory |
|---|---:|---|
| Srinivasavaradhan et al. | 110 bp | `data/processed/srinivas_110` |
| Chandak et al. | 108 bp | `data/processed/chandak_108` |
| Bar-Lev et al. | 140 bp | `data/processed/bar_lev_140` |
| Chandak et al. | 150 bp | `data/processed/chandak_150` |
| Erlich and Zielinski | 152 bp | `data/processed/erlich_152` |

Not every experiment uses every dataset. Each experiment README specifies its
exact input files.

## Baselines

Classical baselines include BBS, MUSCLE, CPL, ITR, and TrellisBMA.

Deep-learning baselines include TReconLM and DNAFormer. RobuSeqNet is retained
as a pinned source snapshot only because a validated pretrained checkpoint
was unavailable.

See [baselines/README.md](baselines/README.md) for build and execution
instructions. Exact upstream versions are recorded in
[baselines/VERSIONS.md](baselines/VERSIONS.md).

## Experiments

Reproduction scripts are organized under `experiments/`:

- component ablation;
- coverage robustness;
- controlled coverage-by-IDS conditions;
- cross-cluster context analysis;
- difficult-cluster mechanism analysis;
- parameter sensitivity;
- deep-learning baseline records.

## Reproducibility notes

- Run commands from the repository root unless stated otherwise.
- Active dataset and executable paths are repository-relative.
- Randomized experiments record their seeds explicitly.
- Full experiments may require substantial CPU or GPU time.
- Smoke tests validate integration and are not reported benchmark results.
- Historical logs may contain paths from the original machine, but active
  scripts do not depend on those paths.

## Citation

Citation information will be added after publication.

## License

Original CCP-DP source code, repository-authored adapters, experiment
orchestration, and documentation are released under the
[MIT License](LICENSE), except where otherwise noted.

Third-party baselines retain their original licenses and attribution. See
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for the precise scope and
provenance of bundled or externally fetched components.

## Third-party software

Third-party baseline implementations retain their original attribution and
license files. Consult the corresponding directory under `baselines/` before
redistributing a baseline implementation or binary.
