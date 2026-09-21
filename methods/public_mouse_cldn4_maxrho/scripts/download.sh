#!/usr/bin/env bash
# Public GEO processed matrices for the integrate cohorts only.
# GSE154977, GSE180963, GSE165641. No private 8 KL. No SRA/FASTQ.
set -euo pipefail

OUT="${1:-/tmp/kpkl_10x}"
mkdir -p "$OUT"

retry() {
  local dest="$1"
  local url="$2"
  if [[ -s "$dest" ]]; then
    echo "have $dest"
    return 0
  fi
  echo "downloading $url"
  curl -L --fail --retry 5 --retry-delay 4 -o "$dest" "$url"
}

G154="$OUT/GSE154977"
mkdir -p "$G154"
BASE154="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl"
for f in \
  GSE154977_mmLung10x_cis_dSp_rawCount.h5 \
  GSE154977_mmLung10x_cis_smpTable.csv.gz \
  GSE154977_mmLung10x_cis_geneTable.csv.gz
do
  retry "$G154/$f" "$BASE154/$f"
done

G180="$OUT/GSE180963"
mkdir -p "$G180"
retry "$G180/GSE180963_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180963/suppl/GSE180963_RAW.tar"
python3 - <<PY
import tarfile
from pathlib import Path
out = Path("$G180")
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
    raise SystemExit("GSE180963 missing: " + ", ".join(missing))
print("GSE180963 OK")
PY

G165="$OUT/GSE165641"
mkdir -p "$G165"
retry "$G165/GSM5047302_KL1_count.tar.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5047nnn/GSM5047302/suppl/GSM5047302_KL1_count.tar.gz"
retry "$G165/GSM5047303_KL2_count.tar.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5047nnn/GSM5047303/suppl/GSM5047303_KL2_count.tar.gz"
python3 - <<PY
import tarfile
from pathlib import Path
out = Path("$G165")
for name in ("GSM5047302_KL1_count.tar.gz", "GSM5047303_KL2_count.tar.gz"):
    with tarfile.open(out / name) as tf:
        tf.extractall(out)
hits = list(out.rglob("filtered_feature_bc_matrix"))
if len(hits) < 2:
    raise SystemExit("GSE165641 filtered matrices not found")
print("GSE165641 OK", len(hits))
PY

echo "matrices ready under $OUT"
