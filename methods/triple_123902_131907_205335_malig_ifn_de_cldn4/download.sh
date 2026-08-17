#!/usr/bin/env bash
# Public processed GEO only. GSE123902 RAW.tar is 90.4 MB (<2 GB).
# GSE131907 / GSE205335 malignant UMI-sums are already in data/ (PR #456).
set -euo pipefail
DEST="${1:-/tmp/geo_dl}"
mkdir -p "$DEST"
curl -fL --max-time 300 -o "$DEST/GSE123902_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar"
ls -lh "$DEST/GSE123902_RAW.tar"
