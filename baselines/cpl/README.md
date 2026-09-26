# CPL baseline

This directory contains the CPL implementation from
https://github.com/itaiorr/Deep-DNA-based-storage at commit
`167a68261719ea458e79c40fcc485ae768c1a64a`.

`main_cpl32.cpp` is the CCP-DP benchmark adapter. It reads separator-delimited
clusters, performs deterministic per-cluster reconstruction in parallel, and
writes one sequence per cluster.

See `LICENSE` and `README_UPSTREAM.md` for upstream attribution.

## Build

```bash
./baselines/build_baselines.sh
```

Manual build:

```bash
g++ -O3 -march=native -std=c++17 -fopenmp \
  -Wall -Wextra -Wpedantic \
  -include climits \
  baselines/cpl/main_cpl32.cpp \
  baselines/cpl/Cluster.cpp \
  baselines/cpl/EditDistance.cpp \
  baselines/cpl/FreqFunctions.cpp \
  baselines/cpl/Graph.cpp \
  baselines/cpl/GuessFunctions.cpp \
  baselines/cpl/Utils.cpp \
  -o build/baselines/cpl32
```

`-include climits` supplies integer-limit declarations used by upstream
`Graph.cpp` without altering that upstream source file.

## Run

```bash
CLUSTERS=/path/to/Clusters.txt
LENGTH=110
SEPARATOR="===="
JOBS=32
OUTPUT=results/cpl.txt

/usr/bin/time -v -o results/cpl.time.txt \
  ./build/baselines/cpl32 \
  -i "$CLUSTERS" -l "$LENGTH" -s "$SEPARATOR" \
  -o "$OUTPUT" -t "$JOBS" \
  > results/cpl.log 2>&1
```

The adapter defaults to at most 32 reads per cluster. Use `--max-copies N` to
change this explicitly.

