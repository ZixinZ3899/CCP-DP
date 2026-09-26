#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
RAW_DIR="${SCRIPT_DIR}/raw"
PROCESSED_DIR="${SCRIPT_DIR}/processed"
MODE="${1:---core}"

CNR_COMMIT="6938f44796185902a08381943c2895782886c5c3"
ZENODO_BASE="https://zenodo.org/records/16959565/files"

require_command() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "Missing required command: $1" >&2
        exit 1
    }
}

download_file() {
    local url=$1
    local output=$2
    if [[ -s "$output" ]]; then
        echo "[KEEP] $output"
        return
    fi
    echo "[DOWNLOAD] $url"
    curl --fail --location --retry 5 --retry-delay 2 \
        --output "${output}.tmp" "$url"
    mv "${output}.tmp" "$output"
}

clone_or_update() {
    local url=$1
    local directory=$2
    if [[ ! -d "${directory}/.git" ]]; then
        git clone "$url" "$directory"
    else
        git -C "$directory" fetch origin
    fi
}

prepare_core() {
    require_command git
    require_command curl
    require_command python3
    mkdir -p "$RAW_DIR" "$PROCESSED_DIR"

    local cnr_dir="${RAW_DIR}/srinivas_cnr"
    clone_or_update \
        "https://github.com/microsoft/clustered-nanopore-reads-dataset.git" \
        "$cnr_dir"
    git -C "$cnr_dir" fetch origin "$CNR_COMMIT"
    git -C "$cnr_dir" checkout --detach "$CNR_COMMIT"

    mkdir -p "${PROCESSED_DIR}/srinivas_110"
    python3 "${SCRIPT_DIR}/prepare_datasets.py" filter-empty \
        --clusters "${cnr_dir}/Clusters.txt" \
        --centers "${cnr_dir}/Centers.txt" \
        --output-clusters "${PROCESSED_DIR}/srinivas_110/Clusters.txt" \
        --output-centers "${PROCESSED_DIR}/srinivas_110/Centers.txt" \
        --removed-indices "${PROCESSED_DIR}/srinivas_110/removed_empty_cluster_ids.txt"

    local zenodo_dir="${RAW_DIR}/zenodo_16959565"
    mkdir -p "$zenodo_dir"
    download_file \
        "${ZENODO_BASE}/oligo0_UnderlyingClusters.txt?download=1" \
        "${zenodo_dir}/oligo0_UnderlyingClusters.txt"
    download_file \
        "${ZENODO_BASE}/oligo0refs.txt?download=1" \
        "${zenodo_dir}/oligo0refs.txt"
    download_file \
        "${ZENODO_BASE}/BinnedNanoporeTwoFlowcells_clusters_subsampled.txt?download=1" \
        "${zenodo_dir}/BinnedNanoporeTwoFlowcells_clusters_subsampled.txt"
    download_file \
        "${ZENODO_BASE}/BinnedNanoporeTwoFlowcells_centers_subsampled.txt?download=1" \
        "${zenodo_dir}/BinnedNanoporeTwoFlowcells_centers_subsampled.txt"

    (
        cd "$zenodo_dir"
        printf '%s  %s\n' \
          'fb790f5b42262cbf76d49f9db99845d7' 'oligo0_UnderlyingClusters.txt' \
          '8f506d9cbb59650aaf3833fcf0e0c660' 'oligo0refs.txt' \
          '991af24a2e2b5dc9e2a0a85a0c0afac1' 'BinnedNanoporeTwoFlowcells_clusters_subsampled.txt' \
          'cb9f451b28c6425fe4575b637538c7cb' 'BinnedNanoporeTwoFlowcells_centers_subsampled.txt' \
          | md5sum --check --strict -
    )

    mkdir -p \
        "${PROCESSED_DIR}/chandak_108" \
        "${PROCESSED_DIR}/bar_lev_140"
    cp "${zenodo_dir}/oligo0_UnderlyingClusters.txt" \
       "${PROCESSED_DIR}/chandak_108/Clusters.txt"
    cp "${zenodo_dir}/oligo0refs.txt" \
       "${PROCESSED_DIR}/chandak_108/Centers.txt"
    cp "${zenodo_dir}/BinnedNanoporeTwoFlowcells_clusters_subsampled.txt" \
       "${PROCESSED_DIR}/bar_lev_140/Clusters.txt"
    cp "${zenodo_dir}/BinnedNanoporeTwoFlowcells_centers_subsampled.txt" \
       "${PROCESSED_DIR}/bar_lev_140/Centers.txt"

    python3 "${SCRIPT_DIR}/prepare_datasets.py" validate \
        --config "${SCRIPT_DIR}/datasets.json" \
        --dataset srinivas_110 \
        --dataset chandak_108 \
        --dataset bar_lev_140
}

fetch_public_sources() {
    require_command git
    mkdir -p "$RAW_DIR"
    clone_or_update \
        "https://github.com/shubhamchandak94/LDPC_DNA_storage.git" \
        "${RAW_DIR}/LDPC_DNA_storage"
    clone_or_update \
        "https://github.com/shubhamchandak94/LDPC_DNA_storage_data.git" \
        "${RAW_DIR}/LDPC_DNA_storage_data"
    clone_or_update \
        "https://github.com/TeamErlich/dna-fountain.git" \
        "${RAW_DIR}/dna-fountain"
    echo "Public source repositories downloaded."
    echo "Large ENA reads and the Grass supplementary data are not downloaded."
}

case "$MODE" in
    --core)
        prepare_core
        ;;
    --sources)
        fetch_public_sources
        ;;
    --all-public)
        prepare_core
        fetch_public_sources
        ;;
    --help|-h)
        cat <<'EOF'
Usage: ./data/download_datasets.sh [--core|--sources|--all-public]

  --core        Download and prepare the three primary CCP-DP benchmarks.
  --sources     Clone public source repositories for additional datasets.
  --all-public  Perform both operations.

Raw and processed data remain under data/ and are ignored by Git.
EOF
        ;;
    *)
        echo "Unknown option: $MODE" >&2
        exit 2
        ;;
esac

