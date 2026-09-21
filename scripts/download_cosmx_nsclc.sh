#!/usr/bin/env bash
# Download official CosMx NSCLC tarballs and extract exprMat + metadata only.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${1:-$ROOT/data/cosmx}"
mkdir -p "$DEST"
TMP="${TMPDIR:-/tmp}/cosmx_tar"
mkdir -p "$TMP"

SAMPLES=(
  Lung5_Rep1
  Lung5_Rep2
  Lung5_Rep3
  Lung6
  Lung9_Rep1
  Lung9_Rep2
  Lung12
  Lung13
)

BASE="https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed"

have_sample() {
  local s="$1"
  find "$DEST" -type f -name "${s}_exprMat_file.csv" | grep -q .
}

download_one() {
  local s="$1"
  if have_sample "$s"; then
    echo "[skip] $s already extracted"
    return 0
  fi
  local url="${BASE}/${s}/${s}+SMI+Flat+data.tar.gz"
  local tarpath="${TMP}/${s}.tar.gz"
  echo "[get] $s -> $tarpath"
  curl -fL --globoff --retry 4 --retry-all-errors --retry-delay 8 --connect-timeout 30 \
    -o "$tarpath" "$url"
  echo "[list] $s"
  mapfile -t members < <(tar -tzf "$tarpath" | grep -E "${s}_(exprMat|metadata)_file\.csv$")
  if [[ ${#members[@]} -lt 2 ]]; then
    echo "[fail] $s: expected exprMat+metadata, got: ${members[*]:-none}" >&2
    return 1
  fi
  echo "[extract] ${members[*]}"
  tar -xzf "$tarpath" -C "$DEST" "${members[@]}"
  rm -f "$tarpath"
  if ! have_sample "$s"; then
    echo "[fail] $s: exprMat missing after extract" >&2
    return 1
  fi
  echo "[ok] $s"
}

if [[ $# -ge 2 ]]; then
  shift
  for s in "$@"; do
    download_one "$s"
  done
else
  for s in "${SAMPLES[@]}"; do
    download_one "$s"
  done
fi

echo "Extracted files:"
find "$DEST" -type f \( -name '*exprMat_file.csv' -o -name '*metadata_file.csv' \) | sort
