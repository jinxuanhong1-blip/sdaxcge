#!/usr/bin/env bash
# Public processed cell-bin H5ADs for GSE328481 (Stereo-XCR-seq LUAD).
set -euo pipefail
DEST="${1:-data/GSE328481}"
mkdir -p "$DEST"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9684nnn"
declare -A FILES=(
  [GSM9684206]=GSM9684206_D06053D2.h5ad
  [GSM9684207]=GSM9684207_D06047C3.h5ad
  [GSM9684208]=GSM9684208_D06047F6.h5ad
  [GSM9684209]=GSM9684209_D06047E1.h5ad
  [GSM9684210]=GSM9684210_D06050A2.h5ad
  [GSM9684211]=GSM9684211_D06047A2.h5ad
  [GSM9684212]=GSM9684212_D06050C2.h5ad
  [GSM9684213]=GSM9684213_D06047D4.h5ad
  [GSM9684214]=GSM9684214_D06047E2.h5ad
  [GSM9684215]=GSM9684215_D06050D4.h5ad
  [GSM9684216]=GSM9684216_D06050E4.h5ad
)
for acc in "${!FILES[@]}"; do
  fn="${FILES[$acc]}"
  out="${DEST}/${fn}"
  if [[ -s "$out" ]]; then
    echo "exists $fn"
    continue
  fi
  curl -L --fail --retry 4 --retry-delay 8 -C - -o "$out" "${BASE}/${acc}/suppl/${fn}"
done
ls -lh "$DEST"
