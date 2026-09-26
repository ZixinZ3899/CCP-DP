# MUSCLE baseline

This baseline aligns the reads in each cluster with MUSCLE 5.3 and constructs
a deterministic column-wise consensus. Gap characters participate in voting;
a column is removed when the gap count is at least the count of the best base.

## Licensing and installation

MUSCLE is third-party software distributed under the GNU General Public
License v3. Its executable is not committed to this repository.

Install the pinned Bioconda build with:

    conda env create \
      -f baselines/muscle/environment.yml

    conda activate ccpdp-muscle

Verify the installation:

    muscle -version

The expected version is MUSCLE 5.3.

The original formal experiment used Bioconda `muscle=5.3`, build
`h9948957_3`. Exact version and checksum records are preserved in
`../VERSIONS.md`.

The upstream license is retained in `LICENSE`. Upstream source code and
release information are available from:

    https://github.com/rcedgar/muscle

## Included files

- `run_muscle_consensus.py`: cluster parser, parallel runner, and consensus code
- `environment.yml`: pinned Bioconda MUSCLE environment
- `LICENSE`: upstream GNU General Public License v3

## Run

Run commands from the repository root:

    conda activate ccpdp-muscle

    CLUSTERS=/path/to/Clusters.txt
    SEPARATOR="===="
    JOBS=32
    OUTPUT=results/muscle.txt
    MUSCLE_BIN="$(command -v muscle)"

    /usr/bin/time -v \
      python baselines/muscle/run_muscle_consensus.py \
      -i "$CLUSTERS" \
      -s "$SEPARATOR" \
      -o "$OUTPUT" \
      --muscle "$MUSCLE_BIN" \
      --jobs "$JOBS" \
      --muscle-threads 1 \
      2> results/muscle.time.txt

`--jobs 32 --muscle-threads 1` runs independent clusters in parallel without
creating nested 32-by-32 threading.

## Alternative installation

A separately installed official MUSCLE 5.3 executable may also be used:

    python baselines/muscle/run_muscle_consensus.py \
      -i /path/to/Clusters.txt \
      -s "====" \
      -o /path/to/output.txt \
      --muscle /path/to/muscle-5.3 \
      --jobs 32 \
      --muscle-threads 1

The supplied executable must report version 5.3.

## Reproducibility note

An official Linux x86-64 MUSCLE 5.3 binary was compared with the formal
Bioconda build during migration. Both produced byte-identical output in the
10-cluster integration test. The binary itself is intentionally not
redistributed by this repository.
