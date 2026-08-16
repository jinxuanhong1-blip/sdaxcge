#!/usr/bin/env bash
# Download public GSE221322 GeoMx protein DSP tables (Monkman et al.).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DEST="${GSE221322_DATA:-$ROOT/data/gse221322}"
mkdir -p "$DEST"

BASE="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE221322&format=file&file="
FTP_SUPPL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221322/suppl"
FTP_MATRIX="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221322/matrix"

fetch() {
  local url="$1" out="$2"
  if [[ -s "$out" ]]; then
    echo "exists: $out"
    return 0
  fi
  echo "GET $url"
  curl -fsSL -o "$out" "$url"
}

fetch "${BASE}GSE221322_4301_protein_QC.csv.gz" "$DEST/GSE221322_4301_protein_QC.csv.gz"
fetch "${BASE}GSE221322_4301_protein_norm.csv.gz" "$DEST/GSE221322_4301_protein_norm.csv.gz"
fetch "$FTP_MATRIX/GSE221322_series_matrix.txt.gz" "$DEST/GSE221322_series_matrix.txt.gz"

# Fallback if the GEO download CGI is empty
for f in GSE221322_4301_protein_QC.csv.gz GSE221322_4301_protein_norm.csv.gz; do
  if [[ ! -s "$DEST/$f" ]]; then
    fetch "$FTP_SUPPL/$f" "$DEST/$f"
  fi
done

python3 - << PY
from pathlib import Path
import gzip
dest = Path("$DEST")
for name in [
    "GSE221322_4301_protein_QC.csv.gz",
    "GSE221322_4301_protein_norm.csv.gz",
    "GSE221322_series_matrix.txt.gz",
]:
    p = dest / name
    out = dest / name.replace(".gz", "")
    if out.exists() and out.stat().st_size > 0:
        print("unzipped exists:", out)
        continue
    with gzip.open(p, "rb") as src, open(out, "wb") as dst:
        dst.write(src.read())
    print("unzipped", out, "bytes", out.stat().st_size)
PY

echo "GSE221322 files in $DEST"
ls -lh "$DEST"
