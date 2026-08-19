#!/usr/bin/env bash
# Download public GSE221733 GeoMx CTA RNA matrices + series matrix metadata.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DEST="${GSE221733_DATA:-$ROOT/data/gse221733}"
mkdir -p "$DEST"

GEO_CGI="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE221733&format=file&file="
FTP_SUPPL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221733/suppl"
FTP_MATRIX="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221733/matrix"

fetch() {
  local url="$1" out="$2"
  if [[ -s "$out" ]]; then
    echo "exists: $out"
    return 0
  fi
  echo "GET $url"
  curl -fL --retry 5 --retry-delay 4 -o "$out" "$url"
}

fetch "${GEO_CGI}GSE221733_4301_CTA_norm.xlsx" "$DEST/GSE221733_4301_CTA_norm.xlsx"
fetch "${GEO_CGI}GSE221733_4301_CTA_QC.xlsx" "$DEST/GSE221733_4301_CTA_QC.xlsx"
fetch "${GEO_CGI}GSE221733_4301_CTA_initial.csv.gz" "$DEST/GSE221733_4301_CTA_initial.csv.gz"
fetch "$FTP_MATRIX/GSE221733_series_matrix.txt.gz" "$DEST/GSE221733_series_matrix.txt.gz"

# FTP fallback if the GEO CGI wrote an empty file
for f in GSE221733_4301_CTA_norm.xlsx GSE221733_4301_CTA_QC.xlsx GSE221733_4301_CTA_initial.csv.gz; do
  if [[ ! -s "$DEST/$f" ]]; then
    fetch "$FTP_SUPPL/$f" "$DEST/$f"
  fi
done

if [[ ! -s "$DEST/GSE221733_series_matrix.txt" ]]; then
  python3 - << PY
import gzip
from pathlib import Path
src = Path("$DEST") / "GSE221733_series_matrix.txt.gz"
dst = Path("$DEST") / "GSE221733_series_matrix.txt"
with gzip.open(src, "rb") as fh, open(dst, "wb") as out:
    out.write(fh.read())
print("unzipped", dst, "bytes", dst.stat().st_size)
PY
fi

echo "GSE221733 files in $DEST"
ls -lh "$DEST"
