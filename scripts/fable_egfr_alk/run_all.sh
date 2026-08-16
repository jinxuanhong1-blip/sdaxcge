#!/usr/bin/env bash
# Reproduce the full TACSTD2/CLDN4 vs immune/response analysis slice.
# Raw inputs are downloaded into data/fable_egfr_alk/ (git-ignored).
set -euo pipefail
cd "$(dirname "$0")"

DATA="../../data/fable_egfr_alk"
mkdir -p "$DATA"

echo "==> [1/6] Download public data (idempotent)"
dl () { # url  outfile
  if [ ! -s "$DATA/$2" ]; then
    echo "    fetch $2"
    curl -fsSL --retry 4 --retry-delay 4 -o "$DATA/$2" "$1"
  else
    echo "    have  $2"
  fi
}
GDC="https://gdc-hub.s3.us-east-1.amazonaws.com/download"
for C in LUAD LUSC; do
  dl "$GDC/TCGA-$C.star_tpm.tsv.gz"            "TCGA-$C.star_tpm.tsv.gz"
  dl "$GDC/TCGA-$C.somaticmutation_wxs.tsv.gz" "TCGA-$C.somaticmutation_wxs.tsv.gz"
  dl "$GDC/TCGA-$C.survival.tsv.gz"            "TCGA-$C.survival.tsv.gz"
done
GEO="https://ftp.ncbi.nlm.nih.gov/geo/series"
dl "$GEO/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz"                 "GSE126044_counts.txt.gz"
dl "$GEO/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz"        "GSE126044_series_matrix.txt.gz"
dl "$GEO/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz" "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
dl "$GEO/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz"        "GSE135222_series_matrix.txt.gz"

echo "==> [2/6] Prepare TCGA-LUAD";  python3 00_prepare_tcga.py LUAD
echo "==> [3/6] Prepare TCGA-LUSC";  python3 00_prepare_tcga.py LUSC
echo "==> [4/6] TCGA immune axis";   python3 01_tcga_immune.py LUAD; python3 01_tcga_immune.py LUSC
echo "==> [5/6] ICI response axis";  python3 02_ici_response.py
echo "==> [6/6] Verification";       python3 03_verify.py
echo "==> DONE. Outputs in results/fable_egfr_alk/"
