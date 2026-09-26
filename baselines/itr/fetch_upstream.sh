#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
    pwd
)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)

UPSTREAM_URL="https://github.com/omersabary/Reconstruction.git"
UPSTREAM_COMMIT="c50dec739bd2c7f18ac7678d921eecee02e86c6a"
DESTINATION="${1:-$REPO_ROOT/build/third_party/Reconstruction}"

if [[ -e "$DESTINATION" && ! -d "$DESTINATION/.git" ]]; then
    echo "[ERROR] Destination exists but is not a Git repository:"
    echo "        $DESTINATION"
    exit 1
fi

mkdir -p "$(dirname -- "$DESTINATION")"

if [[ ! -d "$DESTINATION/.git" ]]; then
    echo "[CLONE] $UPSTREAM_URL"
    git clone "$UPSTREAM_URL" "$DESTINATION"
fi

echo "[FETCH] $UPSTREAM_COMMIT"
git -C "$DESTINATION" fetch origin "$UPSTREAM_COMMIT"

echo "[CHECKOUT] $UPSTREAM_COMMIT"
git -C "$DESTINATION" checkout --detach "$UPSTREAM_COMMIT"

ACTUAL_COMMIT=$(git -C "$DESTINATION" rev-parse HEAD)
if [[ "$ACTUAL_COMMIT" != "$UPSTREAM_COMMIT" ]]; then
    echo "[ERROR] Unexpected ITR commit: $ACTUAL_COMMIT"
    exit 1
fi

echo "[PASS] ITR upstream is pinned to $ACTUAL_COMMIT"
echo "[PATH] $DESTINATION"
echo
echo "Note: the upstream README currently states 'License TBA'."
echo "The downloaded source remains under the upstream authors' terms."

