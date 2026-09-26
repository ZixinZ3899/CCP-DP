#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CENTERS="${REPO_ROOT}/data/processed/srinivas_110/Centers.txt"
SOURCE="${REPO_ROOT}/src/main.cpp"
HEADER="${REPO_ROOT}/include/ccpdp.hpp"
BINARY="${REPO_ROOT}/build/ccpdp"
CONFIG="${SCRIPT_DIR}/experiment_config.srinivas.json"

# Use dedicated directories to keep smoke and full results separate.
SMOKE_DATA="${SCRIPT_DIR}/smoke_data"
SMOKE_RESULTS="${SCRIPT_DIR}/smoke_results"
STRESS_DATA="${SCRIPT_DIR}/stress_data"
STRESS_RESULTS="${SCRIPT_DIR}/stress_results"
STRESS_FIGURE="${SCRIPT_DIR}/stress_figure"

usage() {
  echo "Usage: $0 {compile|commands|smoke|generate|full|resume|plot|status}"
  echo "  compile   compile the CCP-DP adaptive-routing binary"
  echo "  commands  generate a 50-cluster smoke input and print resolved commands"
  echo "  smoke     compile and run one 50-cluster condition with all methods"
  echo "  generate  generate the 100 paired Coverage x IDS conditions only"
  echo "  full      compile, generate/verify data, and run all 300 evaluations"
  echo "  resume    resume the 300 evaluations without regenerating data"
  echo "  plot      draw PNG/SVG/PDF from the completed results"
  echo "  status    report generated conditions and completed metric rows"
}

compile_ccpdp() {
  [[ -f "${SOURCE}" ]] || { echo "Missing source: ${SOURCE}" >&2; exit 1; }
  [[ -f "${HEADER}" ]] || { echo "Missing header: ${HEADER}" >&2; exit 1; }
  if [[ ! -x "${BINARY}" || "${SOURCE}" -nt "${BINARY}" || "${HEADER}" -nt "${BINARY}" ]]; then
    echo "Compiling ${BINARY}"
    g++ -O3 -march=native -std=c++17 -fopenmp \
      -Wall -Wextra -Wpedantic \
      -I"${REPO_ROOT}/include" "${SOURCE}" -o "${BINARY}"
  else
    echo "Binary is current: ${BINARY}"
  fi
}

generate_smoke() {
  python3 "${SCRIPT_DIR}/generate_ids_stress.py" \
    --centers "${CENTERS}" \
    --outdir "${SMOKE_DATA}" \
    --coverage 5 --indel-pct 6 --sub-pct 2 --seeds 2026 \
    --limit-centers 50
}

generate_full() {
  python3 "${SCRIPT_DIR}/generate_ids_stress.py" \
    --centers "${CENTERS}" \
    --outdir "${STRESS_DATA}" \
    --coverage 5,10,15,20 \
    --indel-pct 2,4,6,8,10 \
    --sub-pct 2 \
    --seeds 2026,2027,2028,2029,2030
}

verify_full_data() {
  local count
  if [[ -d "${STRESS_DATA}" ]]; then
    count="$(find "${STRESS_DATA}" -mindepth 2 -maxdepth 2 -name metadata.json | wc -l)"
  else
    count=0
  fi
  if [[ "${count}" -ne 100 ]]; then
    echo "Expected 100 conditions in ${STRESS_DATA}, found ${count}." >&2
    echo "Run: $0 generate" >&2
    exit 1
  fi
  if ! grep -q '"generator_version": "coverage-ids-v2-paired"' \
      "${STRESS_DATA}/cov05_indel02_seed2026/metadata.json"; then
    echo "The data directory was not produced by the corrected paired generator." >&2
    echo "Use a fresh ${STRESS_DATA} directory and run: $0 generate" >&2
    exit 1
  fi
}

run_full() {
  verify_full_data
  python3 "${SCRIPT_DIR}/run_ids_stress.py" \
    --config "${CONFIG}" \
    --data-root "${STRESS_DATA}" \
    --results-root "${STRESS_RESULTS}"
}

case "${1:-}" in
  compile)
    compile_ccpdp
    ;;
  commands)
    compile_ccpdp
    generate_smoke
    python3 "${SCRIPT_DIR}/run_ids_stress.py" \
      --config "${CONFIG}" \
      --data-root "${SMOKE_DATA}" \
      --results-root "${SMOKE_RESULTS}" \
      --dry-run
    ;;
  smoke)
    compile_ccpdp
    generate_smoke
    python3 "${SCRIPT_DIR}/run_ids_stress.py" \
      --config "${CONFIG}" \
      --data-root "${SMOKE_DATA}" \
      --results-root "${SMOKE_RESULTS}" \
      --force
    ;;
  generate)
    generate_full
    verify_full_data
    ;;
  full)
    compile_ccpdp
    generate_full
    run_full
    ;;
  resume)
    compile_ccpdp
    run_full
    ;;
  plot)
    python3 "${SCRIPT_DIR}/plot_ids_stress.py" \
      --metrics "${STRESS_RESULTS}/raw_metrics.csv" \
      --outdir "${STRESS_FIGURE}" \
      --show-sd
    ;;
  status)
    if [[ -d "${STRESS_DATA}" ]]; then
      conditions="$(find "${STRESS_DATA}" -mindepth 2 -maxdepth 2 -name metadata.json | wc -l)"
    else
      conditions=0
    fi
    rows=0
    if [[ -f "${STRESS_RESULTS}/raw_metrics.csv" ]]; then
      rows="$(( $(wc -l < "${STRESS_RESULTS}/raw_metrics.csv") - 1 ))"
    fi
    echo "Conditions: ${conditions}/100"
    echo "Metric rows: ${rows}/300"
    ;;
  *)
    usage
    exit 2
    ;;
esac
