#!/usr/bin/env bash
# Public processed GSE189357 10x MTX only. No FASTQ. No GSE189487. No GSE148071.
set -euo pipefail
OUT="${1:-/tmp/gse189357}"
mkdir -p "$OUT"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar"
DEST="$OUT/GSE189357_RAW.tar"
if [[ -f "$DEST" ]] && [[ "$(stat -c%s "$DEST")" -gt 1000000 ]]; then
  echo "exists $DEST ($(stat -c%s "$DEST") bytes)"
else
  curl -L --retry 5 --retry-delay 4 -o "$DEST.partial" "$URL"
  mv "$DEST.partial" "$DEST"
  echo "wrote $DEST ($(stat -c%s "$DEST") bytes)"
fi
ls -lh "$DEST"
