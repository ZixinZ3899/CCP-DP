#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CENTERS="${REPO_ROOT}/data/processed/srinivas_110/Centers.txt"
CONFIG="${SCRIPT_DIR}/experiment_config.srinivas.json"

REFINED_DATA="${SCRIPT_DIR}/stress_data_refined"
REFINED_RESULTS="${SCRIPT_DIR}/stress_results_refined"
DENSE_RESULTS="${SCRIPT_DIR}/stress_results_dense"
FIGURE_DIR="${SCRIPT_DIR}/stress_figure_dense"

# Override this on the command line if your completed 5/10/15/20x file lives elsewhere:
# BASE_METRICS=/path/to/raw_metrics.csv ./run_srinivas_stress_dense.sh merge
if [[ -z "${BASE_METRICS:-}" ]]; then
  BASE_METRICS="${SCRIPT_DIR}/results/raw_metrics.csv"
fi
BASE_RESULTS_ROOT="$(dirname "${BASE_METRICS}")"

usage() {
  echo "Usage: $0 {commands|smoke|refine|merge|plot|all}"
  echo "  commands  print one resolved 6x smoke-test command"
  echo "  smoke     run 6x/6% indel on 50 clusters"
  echo "  refine    run only 6x,7x,8x,9x,12x (125 conditions; 375 method runs)"
  echo "  merge     merge new rows with the completed 5x,10x,15x,20x experiment"
  echo "  plot      merge, extract routing, and draw the dense main figure"
  echo "  all       refine + merge + plot"
}

generate_refined() {
  python3 "${SCRIPT_DIR}/generate_ids_stress.py" \
    --centers "${CENTERS}" \
    --outdir "${REFINED_DATA}" \
    --coverage 6,7,8,9,12 \
    --indel-pct 2,4,6,8,10 \
    --sub-pct 2 \
    --seeds 2026,2027,2028,2029,2030
}

run_refined() {
  python3 "${SCRIPT_DIR}/run_ids_stress.py" \
    --config "${CONFIG}" \
    --data-root "${REFINED_DATA}" \
    --results-root "${REFINED_RESULTS}"
}

merge_results() {
  if [[ ! -f "${BASE_METRICS}" ]]; then
    echo "Base metrics not found: ${BASE_METRICS}" >&2
    echo "Set BASE_METRICS=/absolute/path/to/base/raw_metrics.csv" >&2
    exit 1
  fi
  if [[ ! -f "${REFINED_RESULTS}/raw_metrics.csv" ]]; then
    echo "Refined metrics not found; run '$0 refine' first." >&2
    exit 1
  fi
  mkdir -p "${DENSE_RESULTS}"
  python3 "${SCRIPT_DIR}/merge_stress_metrics.py" \
    --base "${BASE_METRICS}" \
    --refined "${REFINED_RESULTS}/raw_metrics.csv" \
    --output "${DENSE_RESULTS}/raw_metrics.csv"

  python3 "${SCRIPT_DIR}/extract_routing_summary.py" \
    --results-root "${BASE_RESULTS_ROOT}" \
    --results-root "${REFINED_RESULTS}" \
    --raw-output "${DENSE_RESULTS}/routing_raw.csv" \
    --summary-output "${DENSE_RESULTS}/routing_summary.csv"
}

plot_results() {
  merge_results
  python3 "${SCRIPT_DIR}/plot_ids_stress_dense.py" \
    --metrics "${DENSE_RESULTS}/raw_metrics.csv" \
    --routing "${DENSE_RESULTS}/routing_summary.csv" \
    --outdir "${FIGURE_DIR}"
}

case "${1:-}" in
  commands)
    python3 "${SCRIPT_DIR}/generate_ids_stress.py" \
      --centers "${CENTERS}" --outdir "${SCRIPT_DIR}/smoke_data_dense" \
      --coverage 6 --indel-pct 6 --sub-pct 2 --seeds 2026 --limit-centers 50
    python3 "${SCRIPT_DIR}/run_ids_stress.py" \
      --config "${CONFIG}" --data-root "${SCRIPT_DIR}/smoke_data_dense" \
      --results-root "${SCRIPT_DIR}/smoke_results_dense" --dry-run
    ;;
  smoke)
    python3 "${SCRIPT_DIR}/generate_ids_stress.py" \
      --centers "${CENTERS}" --outdir "${SCRIPT_DIR}/smoke_data_dense" \
      --coverage 6 --indel-pct 6 --sub-pct 2 --seeds 2026 --limit-centers 50
    python3 "${SCRIPT_DIR}/run_ids_stress.py" \
      --config "${CONFIG}" --data-root "${SCRIPT_DIR}/smoke_data_dense" \
      --results-root "${SCRIPT_DIR}/smoke_results_dense"
    ;;
  refine)
    generate_refined
    run_refined
    ;;
  merge)
    merge_results
    ;;
  plot)
    plot_results
    ;;
  all)
    generate_refined
    run_refined
    plot_results
    ;;
  *)
    usage
    exit 2
    ;;
esac
