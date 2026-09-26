# TrellisBMA baseline

This directory contains the modules required by the final cluster runner from
Microsoft's TrellisBMA repository:

- Upstream: https://github.com/microsoft/TrellisBMA
- Commit: `77cb3d3965548e0ff63a37ff5de13f96fa15f5a8`
- License: MIT
- Local adapter: `run_trellisbma_clusters.py`

See `LICENSE` and `README_UPSTREAM.md` for upstream attribution.

## Environment

The verified environment uses Python 3.10:

```bash
conda env create -f baselines/trellisbma/environment.yml
conda activate ccpdp-trellisbma
```

Alternatively, inside an existing Python 3.10 environment:

```bash
python -m pip install -r baselines/trellisbma/requirements.txt
```

## Run

```bash
CLUSTERS=/path/to/Clusters.txt
LENGTH=110
SEPARATOR="===="
JOBS=32
OUTPUT=results/trellisbma.txt

/usr/bin/time -v -o results/trellisbma.time.txt \
  python baselines/trellisbma/run_trellisbma_clusters.py \
  -i "$CLUSTERS" \
  -o "$OUTPUT" \
  -l "$LENGTH" \
  -s "$SEPARATOR" \
  --max-clusters 0 \
  --max-reads 0 \
  --max-drift 15 \
  --lookahead 1 \
  --jobs "$JOBS" \
  > results/trellisbma.log 2>&1
```

`--max-clusters 0` and `--max-reads 0` mean that all clusters and all reads are
used. The defaults for `--p-del`, `--p-ins`, and `--p-sub` are those used by
the migrated Srinivas benchmark runner; set them explicitly when reproducing a
different channel configuration.

Numba may emit `NumbaPendingDeprecationWarning` messages for reflected lists in
the upstream code. These warnings do not indicate a failed reconstruction.

Full-data TrellisBMA runs can take hours. Use `--max-clusters 10 --jobs 2` for a
quick integration test.

