#!/usr/bin/env bash
# Public GEO processed 10x MTX only. No SRA / FASTQ. No private 8-KL object.
set -euo pipefail
OUT="${1:-/tmp/gse180963}"
mkdir -p "$OUT"
TAR="$OUT/GSE180963_RAW.tar"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180963/suppl/GSE180963_RAW.tar"
SMX="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180963/matrix/GSE180963_series_matrix.txt.gz"

if [[ ! -s "$TAR" ]]; then
  echo "downloading $URL"
  curl -L --fail --retry 5 --retry-delay 4 -o "$TAR" "$URL"
fi
if [[ ! -s "$OUT/GSE180963_series_matrix.txt.gz" ]]; then
  curl -L --fail --retry 5 --retry-delay 4 -o "$OUT/GSE180963_series_matrix.txt.gz" "$SMX" || true
fi

python3 - <<PY
import tarfile
from pathlib import Path
out = Path("$OUT")
tar = out / "GSE180963_RAW.tar"
with tarfile.open(tar) as tf:
    tf.extractall(out)
for name in ("GSM5481386_K.tar.gz", "GSM5481387_KL.tar.gz"):
    inner = out / name
    with tarfile.open(inner) as tf:
        tf.extractall(out)
need = [
    out / "K" / "matrix.mtx",
    out / "K" / "genes.tsv",
    out / "K" / "barcodes.tsv",
    out / "KL" / "matrix.mtx",
    out / "KL" / "genes.tsv",
    out / "KL" / "barcodes.tsv",
]
missing = [str(p) for p in need if not p.exists()]
if missing:
    raise SystemExit("missing processed MTX files: " + ", ".join(missing))
print("ok processed MTX under", out)
PY
