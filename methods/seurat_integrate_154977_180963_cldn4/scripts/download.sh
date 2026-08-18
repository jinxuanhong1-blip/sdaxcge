#!/usr/bin/env bash
# Public GEO processed matrices only. Pair = GSE154977 + GSE180963.
# No SRA / FASTQ. No GSE267321. No private 8-KL. No human.
set -euo pipefail
OUT="${1:-/tmp/pair_154977_180963}"
mkdir -p "$OUT/GSE154977" "$OUT/GSE180963"

# ---- GSE180963 10x MTX ----
TAR="$OUT/GSE180963/GSE180963_RAW.tar"
URL180="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180963/suppl/GSE180963_RAW.tar"
if [[ ! -s "$TAR" ]]; then
  echo "downloading $URL180"
  curl -L --fail --retry 5 --retry-delay 4 -o "$TAR" "$URL180"
fi
python3 - <<PY
import tarfile
from pathlib import Path
out = Path("$OUT") / "GSE180963"
tar = out / "GSE180963_RAW.tar"
need = [
    out / "K" / "matrix.mtx",
    out / "K" / "genes.tsv",
    out / "K" / "barcodes.tsv",
    out / "KL" / "matrix.mtx",
    out / "KL" / "genes.tsv",
    out / "KL" / "barcodes.tsv",
]
if all(p.exists() and p.stat().st_size > 0 for p in need):
    print("ok GSE180963 MTX already extracted")
else:
    with tarfile.open(tar) as tf:
        tf.extractall(out)
    for name in ("GSM5481386_K.tar.gz", "GSM5481387_KL.tar.gz"):
        inner = out / name
        with tarfile.open(inner) as tf:
            tf.extractall(out)
    missing = [str(p) for p in need if not p.exists()]
    if missing:
        raise SystemExit("missing GSE180963 MTX: " + ", ".join(missing))
    print("ok GSE180963 MTX under", out)
PY

# ---- GSE154977 custom COO h5 + tables ----
GEO154="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl"
for f in \
  GSE154977_mmLung10x_cis_dSp_rawCount.h5 \
  GSE154977_mmLung10x_cis_smpTable.csv.gz \
  GSE154977_mmLung10x_cis_geneTable.csv.gz \
  GSE154977_mmLung10x_cis_dZ_annot_annot_smpTable.csv.gz \
  GSE154977_mmLung10x_cis_dZ_QCstat_QCstat_smpTable.csv.gz
do
  dest="$OUT/GSE154977/$f"
  if [[ ! -s "$dest" ]]; then
    echo "downloading $GEO154/$f"
    curl -L --fail --retry 5 --retry-delay 4 -o "$dest" "$GEO154/$f"
  fi
done
python3 - <<PY
from pathlib import Path
d = Path("$OUT") / "GSE154977"
need = [
    d / "GSE154977_mmLung10x_cis_dSp_rawCount.h5",
    d / "GSE154977_mmLung10x_cis_smpTable.csv.gz",
    d / "GSE154977_mmLung10x_cis_geneTable.csv.gz",
]
missing = [str(p) for p in need if not p.exists() or p.stat().st_size < 1000]
if missing:
    raise SystemExit("missing GSE154977 files: " + ", ".join(missing))
print("ok GSE154977 processed files under", d)
PY

echo "pair download complete under $OUT"
