#!/usr/bin/env bash
set -euo pipefail

ROOT=$(
    cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
    pwd
)

THIRD_PARTY="$ROOT/build/third_party"
mkdir -p "$THIRD_PARTY"

TRECONLM_COMMIT="ea8a838e03e03c5eecaa6d11ceb2b9b8c5262b42"
ROBUSEQNET_COMMIT="ff7bb0039e1667884ab9b6ba05af732084c4f052"

SCHEME="https"
SEPARATOR="://"

TRECONLM_URL="${SCHEME}${SEPARATOR}github.com/MLI-lab/TReconLM.git"
ROBUSEQNET_URL="${SCHEME}${SEPARATOR}github.com/qinyunnn/RobuSeqNet.git"

fetch_repository() {
    local name=$1
    local url=$2
    local commit=$3
    local destination="$THIRD_PARTY/$name"

    if [[ -e "$destination" && ! -d "$destination/.git" ]]; then
        echo "[ERROR] $destination exists but is not a Git repository" >&2
        exit 1
    fi

    if [[ ! -d "$destination/.git" ]]; then
        echo "[CLONE] $name"
        git clone "$url" "$destination"
    fi

    if ! git -C "$destination" cat-file -e "${commit}^{commit}" 2>/dev/null; then
        echo "[FETCH] $name $commit"
        git -C "$destination" fetch origin "$commit"
    fi

    git -C "$destination" checkout --detach "$commit"

    actual=$(
        git -C "$destination" rev-parse HEAD
    )

    if [[ "$actual" != "$commit" ]]; then
        echo "[ERROR] $name commit mismatch" >&2
        exit 1
    fi

    echo "[PASS] $name is pinned to $actual"
}

fetch_repository \
    TReconLM \
    "$TRECONLM_URL" \
    "$TRECONLM_COMMIT"

fetch_repository \
    RobuSeqNet \
    "$ROBUSEQNET_URL" \
    "$ROBUSEQNET_COMMIT"

echo
echo "[PASS] Third-party source setup completed"
echo "Model weights must be obtained separately."
