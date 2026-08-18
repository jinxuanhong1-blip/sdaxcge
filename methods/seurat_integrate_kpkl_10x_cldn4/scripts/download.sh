#!/usr/bin/env bash
# Public GEO processed matrices only. No SRA / FASTQ. No private 8-KL.
# Drop any series whose processed matrix is missing after download.
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

# ---- GSE154977 KP 30w (processed COO h5 + tables) ----
G154="$OUT/GSE154977"
mkdir -p "$G154"
BASE154="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl"
for f in \
  GSE154977_mmLung10x_cis_dSp_rawCount.h5 \
  GSE154977_mmLung10x_cis_smpTable.csv.gz \
  GSE154977_mmLung10x_cis_geneTable.csv.gz \
  GSE154977_mmLung10x_cis_dZ_annot_annot_smpTable.csv.gz \
  GSE154977_mmLung10x_cis_dZ_QCstat_QCstat_smpTable.csv.gz
do
  retry "$G154/$f" "$BASE154/$f"
done
if [[ -s "$G154/GSE154977_mmLung10x_cis_dSp_rawCount.h5" ]]; then
  echo "OK GSE154977 processed h5"
  echo "yes" > "$G154/HAS_MATRIX"
else
  echo "DROP GSE154977 — no processed matrix"
  echo "no" > "$G154/HAS_MATRIX"
fi

# ---- GSE180963 K / KL (processed 10x MTX in RAW.tar) ----
G180="$OUT/GSE180963"
mkdir -p "$G180"
retry "$G180/GSE180963_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180963/suppl/GSE180963_RAW.tar"
python3 - <<PY
import tarfile
from pathlib import Path
out = Path("$G180")
tar = out / "GSE180963_RAW.tar"
ok = False
if tar.exists() and tar.stat().st_size > 1000:
    with tarfile.open(tar) as tf:
        tf.extractall(out)
    for name in ("GSM5481386_K.tar.gz", "GSM5481387_KL.tar.gz"):
        inner = out / name
        if inner.exists():
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
    ok = all(p.exists() for p in need)
(out / "HAS_MATRIX").write_text("yes\n" if ok else "no\n")
print("GSE180963 processed MTX", "OK" if ok else "DROP")
PY

# ---- GSE165641 KL1 / KL2 (Cell Ranger filtered MTX) ----
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
ok = True
for name, dest in (
    ("GSM5047302_KL1_count.tar.gz", "KL1_count"),
    ("GSM5047303_KL2_count.tar.gz", "KL2_count"),
):
    inner = out / name
    if not inner.exists() or inner.stat().st_size < 1000:
        ok = False
        continue
    with tarfile.open(inner) as tf:
        tf.extractall(out)
    mtx_dir = out / dest / "filtered_feature_bc_matrix"
    if not mtx_dir.exists():
        # some tars unpack a nested folder
        hits = list(out.rglob("filtered_feature_bc_matrix"))
        if hits:
            print(name, "unpacked to", hits[0])
        else:
            ok = False
(out / "HAS_MATRIX").write_text("yes\n" if ok else "no\n")
print("GSE165641 Cell Ranger MTX", "OK" if ok else "DROP")
PY

echo "==== matrix presence ===="
for s in GSE154977 GSE180963 GSE165641; do
  echo -n "$s "; cat "$OUT/$s/HAS_MATRIX"
done
