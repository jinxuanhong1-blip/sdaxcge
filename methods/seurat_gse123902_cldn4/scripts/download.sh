#!/usr/bin/env bash
# Public processed GSE123902 (GEO dense UMI CSVs). Author 36.5 GB H5 is not required.
set -euo pipefail
OUT="${1:-/tmp/gse123902_seurat}"
mkdir -p "$OUT"
cd "$OUT"
if [[ ! -s GSE123902_RAW.tar ]]; then
  curl -fL --retry 4 --retry-delay 4 -o GSE123902_RAW.tar.partial \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar
  mv GSE123902_RAW.tar.partial GSE123902_RAW.tar
fi
curl -fsSL --retry 4 -o GSE123902_GEO_README.rtf \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_GEO_README.rtf || true
mkdir -p csv
if ! ls csv/GSM*_dense.csv.gz >/dev/null 2>&1; then
  tar -xf GSE123902_RAW.tar -C csv
fi
ls -lh GSE123902_RAW.tar csv | head
