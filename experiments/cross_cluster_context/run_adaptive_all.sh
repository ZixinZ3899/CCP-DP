#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-smoke}"
if [[ "$MODE" != "smoke" && "$MODE" != "full" ]]; then
  echo "Usage: $0 [smoke|full]" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BINARY="${REPO_ROOT}/build/ccpdp"
RUNNER="${SCRIPT_DIR}/run_adaptive_routing.py"
BASE="${SCRIPT_DIR}/results"

if [[ "$MODE" == "smoke" ]]; then
  RESULTS="$BASE/adaptive_smoke"
  SUMMARY="$BASE/adaptive_evidence_allocation_data_smoke.csv"
  REPEATS=1
else
  RESULTS="$BASE/adaptive_runs"
  SUMMARY="$BASE/adaptive_evidence_allocation_data.csv"
  REPEATS=3
fi

mkdir -p "$RESULTS"

run_dataset() {
  local name="$1"
  local length="$2"
  local clusters="$3"
  local centers="$4"

  [[ -f "$clusters" ]] || { echo "Missing clusters: $clusters" >&2; exit 1; }
  [[ -f "$centers" ]] || { echo "Missing centers: $centers" >&2; exit 1; }

  python "$RUNNER" \
    --binary "$BINARY" \
    --dataset "$name" \
    --clusters "$clusters" \
    --centers "$centers" \
    --length "$length" \
    --separator "===============================" \
    --jobs 32 \
    --repeats "$REPEATS" \
    --outdir "$RESULTS/$name" \
    --summary-csv "$SUMMARY"
}

run_dataset \
  "Srinivas-110" 110 \
  "${REPO_ROOT}/data/processed/srinivas_110/Clusters.txt" \
  "${REPO_ROOT}/data/processed/srinivas_110/Centers.txt"

run_dataset \
  "Chandak-108" 108 \
  "${REPO_ROOT}/data/processed/chandak_108/Clusters.txt" \
  "${REPO_ROOT}/data/processed/chandak_108/Centers.txt"

if [[ "$MODE" == "full" ]]; then
  run_dataset \
    "Chandak-150" 150 \
    "${REPO_ROOT}/data/processed/chandak_150/Clusters.txt" \
    "${REPO_ROOT}/data/processed/chandak_150/Centers.txt"

  run_dataset \
    "Erlich-152" 152 \
    "${REPO_ROOT}/data/processed/erlich_152/Clusters.txt" \
    "${REPO_ROOT}/data/processed/erlich_152/Centers.txt"

  if [[ -n "${FASTA2_CLUSTERS:-}" || -n "${FASTA2_CENTERS:-}" ]]; then
    if [[ -z "${FASTA2_CLUSTERS:-}" || -z "${FASTA2_CENTERS:-}" ]]; then
      echo "Set both FASTA2_CLUSTERS and FASTA2_CENTERS, or neither." >&2
      exit 1
    fi
    run_dataset \
      "fasta_2-110" 110 \
      "$FASTA2_CLUSTERS" \
      "$FASTA2_CENTERS"
  else
    echo "WARNING: fasta_2-110 skipped because its original inputs were not found." >&2
  fi
fi

echo "Adaptive summary: $SUMMARY"
