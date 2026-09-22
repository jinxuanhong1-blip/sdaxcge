#!/usr/bin/env bash
# Public GEO only. No private 8KL. No SRA/FASTQ.
set -euo pipefail
OUT="${1:-/workspace/data/paper_funnel_public_mouse}"
mkdir -p "$OUT" "$OUT/extracted/GSE165641" "$OUT/extracted/GSE180963" \
  "$OUT/extracted/GSE264739" "$OUT/extracted/GSE179501" \
  "$OUT/extracted/GSE165641/x10" "$OUT/extracted/GSE180963/x10"
cd "$OUT"

retry() {
  local dest="$1" url="$2"
  if [[ -s "$dest" ]]; then echo "have $dest"; return 0; fi
  echo "GET $dest"
  curl -fL --retry 5 --retry-delay 3 -A "Mozilla/5.0" -o "$dest" "$url"
}

retry GSE137244_counts.fpkm.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/GSE137244_counts.fpkm.csv.gz"
retry GSE165641_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165641/suppl/GSE165641_RAW.tar"
retry GSE180963_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180963/suppl/GSE180963_RAW.tar"
retry GSE179501_XTR_scRNAseq_matrix.mtx.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179501/suppl/GSE179501_XTR_scRNAseq_matrix.mtx.gz"
retry GSE179501_XTR_scRNAseq_features.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179501/suppl/GSE179501_XTR_scRNAseq_features.tsv.gz"
retry GSE179501_XTR_scRNAseq_barcodes.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179501/suppl/GSE179501_XTR_scRNAseq_barcodes.tsv.gz"
retry GSE264739_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264739/suppl/GSE264739_RAW.tar"
retry GSE154977_mmLung10x_cis_dSp_rawCount.h5 \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl/GSE154977_mmLung10x_cis_dSp_rawCount.h5"
retry GSE154977_mmLung10x_cis_geneTable.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl/GSE154977_mmLung10x_cis_geneTable.csv.gz"
retry GSE154977_mmLung10x_cis_smpTable.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl/GSE154977_mmLung10x_cis_smpTable.csv.gz"

# Unpack
tar -xf GSE165641_RAW.tar -C extracted/GSE165641
tar -xf GSE180963_RAW.tar -C extracted/GSE180963
tar -xf GSE264739_RAW.tar -C extracted/GSE264739
tar -xzf extracted/GSE165641/GSM5047302_KL1_count.tar.gz -C extracted/GSE165641/x10
tar -xzf extracted/GSE165641/GSM5047303_KL2_count.tar.gz -C extracted/GSE165641/x10
tar -xzf extracted/GSE180963/GSM5481386_K.tar.gz -C extracted/GSE180963/x10
tar -xzf extracted/GSE180963/GSM5481387_KL.tar.gz -C extracted/GSE180963/x10
cp -f GSE179501_XTR_scRNAseq_matrix.mtx.gz extracted/GSE179501/matrix.mtx.gz
cp -f GSE179501_XTR_scRNAseq_features.tsv.gz extracted/GSE179501/features.tsv.gz
# barcodes deposit is plain text mislabeled .gz
cp -f GSE179501_XTR_scRNAseq_barcodes.tsv.gz extracted/GSE179501/barcodes.tsv
echo "download done -> $OUT"
