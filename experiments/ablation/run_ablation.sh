#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
MODE="${1:-full}"
CXX="${CXX:-g++}"

mkdir -p "${BUILD_DIR}"

"${CXX}" -O3 -march=native -std=c++17 -fopenmp \
  -Wall -Wextra -Wpedantic \
  "${SCRIPT_DIR}/main_ablation.cpp" \
  -o "${BUILD_DIR}/ccpdp_ablation"

"${CXX}" -O3 -march=native -std=c++17 \
  -Wall -Wextra -Wpedantic \
  "${SCRIPT_DIR}/evaluate_ablation.cpp" \
  -o "${BUILD_DIR}/evaluate_ablation"

case "${MODE}" in
  smoke)
    OUTDIR="${OUTDIR:-${SCRIPT_DIR}/results/smoke}"
    REPEATS=1
    LIMIT_ARGS=(--limit 100)
    ;;
  full)
    OUTDIR="${OUTDIR:-${SCRIPT_DIR}/results/reproduced}"
    REPEATS="${REPEATS:-3}"
    LIMIT_ARGS=()
    ;;
  *)
    echo "Usage: $0 [smoke|full]" >&2
    exit 2
    ;;
esac

python3 "${SCRIPT_DIR}/run_ablation.py" \
  --config "${SCRIPT_DIR}/ablation_config.json" \
  --binary "${BUILD_DIR}/ccpdp_ablation" \
  --evaluator "${BUILD_DIR}/evaluate_ablation" \
  --outdir "${OUTDIR}" \
  --jobs "${JOBS:-32}" \
  --repeats "${REPEATS}" \
  "${LIMIT_ARGS[@]}"

python3 "${SCRIPT_DIR}/plot_results.py" \
  --summary "${OUTDIR}/ablation_summary.csv" \
  --outdir "${OUTDIR}"
