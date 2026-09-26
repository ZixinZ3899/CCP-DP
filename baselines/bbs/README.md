# BBS baseline

This directory contains BBS v0.2.0 from
https://github.com/GZHoffie/bbs at commit
`27c39ec5b631498b3c55073038fa6a59821f2974`.

The upstream source is unmodified. See `LICENSE` and `README_UPSTREAM.md`.

## Build

From the CCP-DP repository root:

```bash
./baselines/build_baselines.sh
```

To build only BBS:

```bash
mkdir -p build/baselines
CARGO_TARGET_DIR="$PWD/build/baselines/bbs-target" \
  cargo build --release --locked \
  --manifest-path baselines/bbs/Cargo.toml
```

## Run

```bash
CLUSTERS=/path/to/Clusters.txt
LENGTH=110
SEPARATOR="===="
JOBS=32
OUTPUT=results/bbs.txt

/usr/bin/time -v -o results/bbs.time.txt \
  ./build/baselines/bbs-target/release/bbs \
  "$CLUSTERS" -l "$LENGTH" -s "$SEPARATOR" -t "$JOBS" \
  > "$OUTPUT" \
  2> results/bbs.log
```

BBS writes reconstructed sequences to standard output when `-o` is omitted.
The command above preserves the invocation used in the CCP-DP comparison.

## Check

```bash
./build/baselines/bbs-target/release/bbs --version
wc -l results/bbs.txt
```

