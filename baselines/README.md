# Classical baselines

This directory contains the classical reconstruction baselines used in the
CCP-DP experiments. Each method has its own README with provenance, build
instructions, and an example command.

## Included methods

| Method | Implementation in this repository | Upstream version |
|---|---|---|
| BBS | Upstream Rust source | v0.2.0 |
| MUSCLE | External MUSCLE 5.3 installation plus cluster-consensus wrapper | v5.3 |
| CPL | Upstream C++ source plus OpenMP cluster driver | commit `167a682` |
| ITR | Local OpenMP adapter; upstream fetched on demand | commit `c50dec7` |
| TrellisBMA | Required upstream Python modules plus cluster runner | commit `77cb3d3` |

Exact commit hashes, licenses, binary checksums, and the verified toolchain are
recorded in [VERSIONS.md](VERSIONS.md).

## Input and output format

The examples assume a text file containing one DNA read per line. Clusters are
separated by a line beginning with a repeated equals-sign separator such as
`====`. Empty lines are ignored.

Every runner writes one reconstructed sequence per non-empty input cluster, in
the same cluster order. The matching ground-truth file must therefore contain
one center sequence per line in the same order.

## Requirements

- MUSCLE 5.3 installed separately using the supplied Bioconda environment
- GCC/G++ with C++17 and OpenMP support
- Rust and Cargo for BBS
- Python 3
- A Python 3.10 environment for the pinned TrellisBMA dependencies
- `/usr/bin/time` only when resource measurements are required

## Build

ITR does not have an explicit upstream redistribution license. Fetch its exact
upstream revision before building:

```bash
./baselines/itr/fetch_upstream.sh
```

Build BBS, CPL, and ITR and validate the external MUSCLE installation:

```bash
./baselines/build_baselines.sh
```

The generated programs are:

```text
build/baselines/bbs-target/release/bbs
build/baselines/cpl32
build/baselines/itr_omp
```

Install the TrellisBMA environment separately:

```bash
conda env create -f baselines/trellisbma/environment.yml
conda activate ccpdp-trellisbma
```

## Common variables

Set paths and dataset-specific parameters before using the example commands:

```bash
DATA_ROOT=/path/to/dataset
CLUSTERS="$DATA_ROOT/Clusters_no_empty.txt"
CENTERS="$DATA_ROOT/Centers_no_empty.txt"
LENGTH=110
SEPARATOR="===="
JOBS=32
OUT_DIR=results/baselines/srinivas_110
mkdir -p "$OUT_DIR"
```

The sequence length, separator, and filenames must be changed for other
datasets.

## Run BBS

```bash
/usr/bin/time -v -o "$OUT_DIR/bbs.time.txt" \
  ./build/baselines/bbs-target/release/bbs \
  "$CLUSTERS" -l "$LENGTH" -s "$SEPARATOR" -t "$JOBS" \
  > "$OUT_DIR/bbs.txt" \
  2> "$OUT_DIR/bbs.log"
```

## Run MUSCLE

```bash
/usr/bin/time -v \
  python baselines/muscle/run_muscle_consensus.py \
  -i "$CLUSTERS" \
  -s "$SEPARATOR" \
  -o "$OUT_DIR/muscle.txt" \
  --muscle "$(command -v muscle)" \
  --jobs "$JOBS" \
  --muscle-threads 1 \
  2> "$OUT_DIR/muscle.time.txt"
```

## Run CPL

```bash
/usr/bin/time -v -o "$OUT_DIR/cpl.time.txt" \
  ./build/baselines/cpl32 \
  -i "$CLUSTERS" -l "$LENGTH" -s "$SEPARATOR" \
  -o "$OUT_DIR/cpl.txt" -t "$JOBS" \
  > "$OUT_DIR/cpl.log" 2>&1
```

## Run ITR

```bash
OMP_NUM_THREADS="$JOBS" \
OMP_PROC_BIND=spread \
OMP_PLACES=cores \
/usr/bin/time -v -o "$OUT_DIR/itr.time.txt" \
  ./build/baselines/itr_omp \
  -i "$CLUSTERS" -l "$LENGTH" -s "$SEPARATOR" \
  -o "$OUT_DIR/itr.txt" \
  --threads "$JOBS" --max_reads 25 \
  > "$OUT_DIR/itr.log" 2>&1
```

## Run TrellisBMA

```bash
conda activate ccpdp-trellisbma

/usr/bin/time -v -o "$OUT_DIR/trellisbma.time.txt" \
  python baselines/trellisbma/run_trellisbma_clusters.py \
  -i "$CLUSTERS" -o "$OUT_DIR/trellisbma.txt" \
  -l "$LENGTH" -s "$SEPARATOR" \
  --max-clusters 0 --max-reads 0 \
  --max-drift 15 --lookahead 1 --jobs "$JOBS" \
  > "$OUT_DIR/trellisbma.log" 2>&1
```

The full TrellisBMA run can be substantially slower than the other methods.

## Evaluate a result

```bash
python scripts/evaluate_answers_extended.py \
  -o "$OUT_DIR/bbs.txt" \
  -a "$CENTERS"
```

Repeat the same command for each output file. Evaluation requires an exact
one-to-one correspondence between reconstructed sequences and ground-truth
centers.

## Reproducibility notes

- Do not commit generated binaries under `build/`.
- Do not commit private or locally downloaded datasets.
- Preserve each upstream license and attribution file.
- Use the same sampled cluster files for every method in coverage experiments.
- Record thread counts, wall time, peak RSS, sequence length, and all nondefault
  parameters with each formal benchmark.

