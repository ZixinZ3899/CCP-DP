#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
    pwd
)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
BUILD_DIR="${BUILD_DIR:-$REPO_ROOT/build/baselines}"
CXX="${CXX:-g++}"
PYTHON="${PYTHON:-python3}"
ITR_UPSTREAM_ROOT="${ITR_UPSTREAM_ROOT:-$REPO_ROOT/build/third_party/Reconstruction}"
ITR_SOURCE="$ITR_UPSTREAM_ROOT/Iterative"
MUSCLE_BIN="${MUSCLE_BIN:-}"

mkdir -p "$BUILD_DIR"

echo "[1/5] Checking MUSCLE 5.3"

if [[ -z "$MUSCLE_BIN" ]]; then
    MUSCLE_BIN=$(command -v muscle || true)
fi

if [[ -z "$MUSCLE_BIN" || ! -x "$MUSCLE_BIN" ]]; then
    echo "[ERROR] MUSCLE 5.3 executable was not found." >&2
    echo "Create the supplied environment with:" >&2
    echo "  conda env create -f baselines/muscle/environment.yml" >&2
    echo "  conda activate ccpdp-muscle" >&2
    echo "Alternatively set MUSCLE_BIN=/path/to/muscle-5.3." >&2
    exit 1
fi

MUSCLE_VERSION=$(
    "$MUSCLE_BIN" -version 2>&1 |
    head -n 1
)

echo "$MUSCLE_VERSION"

if [[ "$MUSCLE_VERSION" != *"5.3"* ]]; then
    echo "[ERROR] Expected MUSCLE 5.3, but found:" >&2
    echo "        $MUSCLE_VERSION" >&2
    exit 1
fi

echo "[2/5] Building BBS"

CARGO_TARGET_DIR="$BUILD_DIR/bbs-target" \
    cargo build \
    --release \
    --locked \
    --manifest-path "$REPO_ROOT/baselines/bbs/Cargo.toml"

echo "[3/5] Building CPL"

"$CXX" -O3 -march=native -std=c++17 -fopenmp \
    -Wall -Wextra -Wpedantic \
    -include climits \
    "$REPO_ROOT/baselines/cpl/main_cpl32.cpp" \
    "$REPO_ROOT/baselines/cpl/Cluster.cpp" \
    "$REPO_ROOT/baselines/cpl/EditDistance.cpp" \
    "$REPO_ROOT/baselines/cpl/FreqFunctions.cpp" \
    "$REPO_ROOT/baselines/cpl/Graph.cpp" \
    "$REPO_ROOT/baselines/cpl/GuessFunctions.cpp" \
    "$REPO_ROOT/baselines/cpl/Utils.cpp" \
    -o "$BUILD_DIR/cpl32"

echo "[4/5] Building ITR"

if [[ ! -f "$ITR_SOURCE/Cluster2.cpp" ]]; then
    echo "[ERROR] ITR upstream source was not found at:" >&2
    echo "        $ITR_SOURCE" >&2
    echo "Run ./baselines/itr/fetch_upstream.sh first, or set" >&2
    echo "ITR_UPSTREAM_ROOT=/path/to/Reconstruction." >&2
    exit 1
fi

"$CXX" -O3 -march=native -std=c++17 -fopenmp \
    -Wall -Wextra -Wpedantic \
    -I"$ITR_SOURCE" \
    "$REPO_ROOT/baselines/itr/DNA_omp.cpp" \
    "$ITR_SOURCE/Clone.cpp" \
    "$ITR_SOURCE/Cluster2.cpp" \
    "$ITR_SOURCE/CommonSubstring2.cpp" \
    "$ITR_SOURCE/DividerBMA.cpp" \
    "$ITR_SOURCE/EditDistance.cpp" \
    "$ITR_SOURCE/LCS2.cpp" \
    "$ITR_SOURCE/LongestPath.cpp" \
    -o "$BUILD_DIR/itr_omp"

echo "[5/5] Checking TrellisBMA Python syntax"

"$PYTHON" -m py_compile \
    "$REPO_ROOT/baselines/trellisbma/coded_ids_multiD.py" \
    "$REPO_ROOT/baselines/trellisbma/coded_ids_multiD_jit.py" \
    "$REPO_ROOT/baselines/trellisbma/conv_code.py" \
    "$REPO_ROOT/baselines/trellisbma/helper_functions.py" \
    "$REPO_ROOT/baselines/trellisbma/trellis_bma.py" \
    "$REPO_ROOT/baselines/trellisbma/run_trellisbma_clusters.py"

echo
echo "[PASS] Baseline build completed"
echo "BBS: $BUILD_DIR/bbs-target/release/bbs"
echo "CPL: $BUILD_DIR/cpl32"
echo "ITR: $BUILD_DIR/itr_omp"
echo "MUSCLE: $MUSCLE_BIN"
echo "TrellisBMA: Python runner; install its environment separately"
