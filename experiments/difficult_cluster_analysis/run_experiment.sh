#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

BIN="${CCPDP_BIN:-${REPO_ROOT}/build/ccpdp}"
DATA_ROOT="${DATA_ROOT:-${REPO_ROOT}/data/processed}"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/results/reproduced}"

JOBS="${JOBS:-32}"
DATASETS="${DATASETS:-both}"
OVERWRITE="${OVERWRITE:-0}"

if [[ ! -x "${BIN}" ]]; then
    echo "CCP-DP executable not found or not executable: ${BIN}" >&2
    exit 1
fi

if [[ "${DATASETS}" != "both" &&
      "${DATASETS}" != "srinivas_110" &&
      "${DATASETS}" != "chandak_108" ]]; then
    echo "DATASETS must be: both, srinivas_110, or chandak_108" >&2
    exit 2
fi

LIMIT_ARGS=()

if [[ -n "${LIMIT:-}" ]]; then
    if [[ ! "${LIMIT}" =~ ^[1-9][0-9]*$ ]]; then
        echo "LIMIT must be a positive integer" >&2
        exit 2
    fi
    LIMIT_ARGS=(--limit "${LIMIT}")
fi

mkdir -p "${OUT_DIR}"

run_one() {
    local key="$1"
    local prefix="$2"
    local input="$3"
    local length="$4"
    local separator="$5"

    local final_file="${OUT_DIR}/${prefix}_final.txt"
    local diag_file="${OUT_DIR}/${prefix}_diag.csv"
    local guide_file="${OUT_DIR}/${prefix}_guides.txt"
    local log_file="${OUT_DIR}/${prefix}.log"
    local command_file="${OUT_DIR}/${prefix}.command.txt"

    if [[ ! -f "${input}" ]]; then
        echo "Input dataset not found: ${input}" >&2
        exit 1
    fi

    if [[ "${OVERWRITE}" != "1" ]] &&
       [[ -e "${final_file}" ||
          -e "${diag_file}" ||
          -e "${guide_file}" ||
          -e "${log_file}" ]]; then
        echo "Outputs already exist for ${key}: ${OUT_DIR}" >&2
        echo "Use another OUT_DIR or set OVERWRITE=1." >&2
        exit 1
    fi

    local cmd=(
        "${BIN}"
        -i "${input}"
        -l "${length}"
        -s "${separator}"
        -o "${final_file}"
        --diag "${diag_file}"
        --guide-output "${guide_file}"
        --mode auto
        --jobs "${JOBS}"
    )

    cmd+=("${LIMIT_ARGS[@]}")

    printf '%q ' "${cmd[@]}" > "${command_file}"
    printf '\n' >> "${command_file}"

    echo "[RUN] ${key}"
    echo "      input:  ${input}"
    echo "      output: ${final_file}"

    "${cmd[@]}" 2> "${log_file}"

    echo "[DONE] ${key}"
    grep -E \
        "build:|clusters:|graph-decoded clusters:|planned structured decoding:|written output:|written diag:|elapsed seconds:" \
        "${log_file}" || true
}

if [[ "${DATASETS}" == "both" ||
      "${DATASETS}" == "srinivas_110" ]]; then
    run_one \
        "srinivas_110" \
        "srinivas" \
        "${DATA_ROOT}/srinivas_110/Clusters.txt" \
        110 \
        "==============================="
fi

if [[ "${DATASETS}" == "both" ||
      "${DATASETS}" == "chandak_108" ]]; then
    run_one \
        "chandak_108" \
        "chandak108" \
        "${DATA_ROOT}/chandak_108/Clusters.txt" \
        108 \
        "===="
fi

echo "All requested CCP-DP runs completed."
echo "Results: ${OUT_DIR}"