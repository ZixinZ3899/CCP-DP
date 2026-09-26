# ITR baseline

This directory contains `DNA_omp.cpp`, the OpenMP and command-line adapter used
for the Iterative Reconstruction (ITR) comparison.

The upstream repository is https://github.com/omersabary/Reconstruction at
commit `c50dec739bd2c7f18ac7678d921eecee02e86c6a`.

The upstream README states `License TBA`, so the upstream source is not copied
into this repository. `fetch_upstream.sh` obtains the exact revision directly
from the authors' repository under their terms.

## Fetch upstream source

```bash
./baselines/itr/fetch_upstream.sh
```

The source is placed under `build/third_party/Reconstruction`, which should
remain excluded from Git.

## Build

```bash
./baselines/build_baselines.sh
```

To build against an existing checkout instead:

```bash
ITR_UPSTREAM_ROOT=/path/to/Reconstruction \
  ./baselines/build_baselines.sh
```

## Run

```bash
CLUSTERS=/path/to/Clusters.txt
LENGTH=140
SEPARATOR="===="
JOBS=32
OUTPUT=results/itr.txt

OMP_NUM_THREADS="$JOBS" \
OMP_PROC_BIND=spread \
OMP_PLACES=cores \
/usr/bin/time -v -o results/itr.time.txt \
  ./build/baselines/itr_omp \
  -i "$CLUSTERS" -l "$LENGTH" -s "$SEPARATOR" \
  -o "$OUTPUT" \
  --threads "$JOBS" \
  --max_reads 25 \
  > results/itr.log 2>&1
```

The adapter also accepts `--seed` and `--limit`. `--limit` is intended only
for smoke tests.

