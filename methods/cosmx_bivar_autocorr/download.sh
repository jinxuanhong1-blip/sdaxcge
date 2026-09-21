#!/usr/bin/env bash
# Figshare 25976224: CellCharter-clustered He et al. 2022 CosMx NSCLC h5ad.
# https://doi.org/10.6084/m9.figshare.25976224
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEST="${1:-$ROOT/data/cosmx_human_nsclc_clustered.h5ad}"
mkdir -p "$(dirname "$DEST")"
if [[ -f "$DEST" ]]; then
  echo "[skip] $DEST"
  exit 0
fi
echo "[get] figshare 46841842 -> $DEST"
curl -fL --retry 5 --retry-all-errors --retry-delay 5 \
  -o "$DEST" "https://ndownloader.figshare.com/files/46841842"
ls -lh "$DEST"
