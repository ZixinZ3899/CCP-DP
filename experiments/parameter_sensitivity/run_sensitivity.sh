#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${root_dir}/../.." && pwd)"
mode="${1:-full}"
config="${2:-${root_dir}/datasets_server.json}"
python_cmd="${PYTHON:-python}"
binary="${CCPDP_SENSITIVITY_BINARY:-${repo_root}/build/ccpdp_sensitivity}"
jobs="${JOBS:-32}"

if [[ "$mode" == "smoke" ]]; then
  output="${OUTPUT_DIR:-${root_dir}/results/smoke}"
  "$python_cmd" "${root_dir}/run_parameter_sensitivity.py" \
    --config "$config" --output "$output" \
    --binary "$binary" --jobs "$jobs" --smoke \
    --datasets Srinivas-110 Chandak-108
elif [[ "$mode" == "smoke-all" ]]; then
  output="${OUTPUT_DIR:-${root_dir}/results/smoke_all}"
  "$python_cmd" "${root_dir}/run_parameter_sensitivity.py" \
    --config "$config" --output "$output" \
    --binary "$binary" --jobs "$jobs" --smoke
elif [[ "$mode" == "full" ]]; then
  output="${OUTPUT_DIR:-${root_dir}/results/reproduced}"
  "$python_cmd" "${root_dir}/run_parameter_sensitivity.py" \
    --config "$config" --output "$output" \
    --binary "$binary" --jobs "$jobs" --repeats 3 \
    --datasets Srinivas-110 Chandak-108
  "$python_cmd" "${root_dir}/plot_parameter_sensitivity.py" \
    "$output/summary_metrics.csv"
elif [[ "$mode" == "full-all" ]]; then
  output="${OUTPUT_DIR:-${root_dir}/results/reproduced_all}"
  "$python_cmd" "${root_dir}/run_parameter_sensitivity.py" \
    --config "$config" --output "$output" \
    --binary "$binary" --jobs "$jobs" --repeats 3
  "$python_cmd" "${root_dir}/plot_parameter_sensitivity.py" \
    "$output/summary_metrics.csv"
else
  echo "Usage: $0 {smoke|smoke-all|full|full-all} [config.json]" >&2
  exit 2
fi
