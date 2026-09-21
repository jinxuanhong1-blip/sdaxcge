#!/usr/bin/env bash
# Public GEO processed matrices only. No SRA, no FASTQ, no private 8-KL.
set -euo pipefail

OUT="${1:-/tmp/public_mouse_scrna}"
mkdir -p "$OUT"
cd "$OUT"

retry() {
  local dest="$1"
  local url="$2"
  if [[ -s "$dest" ]]; then
    echo "have $dest"
    return 0
  fi
  echo "GET $url"
  curl -fL --retry 5 --retry-delay 4 -A "Mozilla/5.0" -o "$dest" "$url"
}

retry GSE154989_mmLungPlate_fQC_dSp_rawCount.h5 \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/GSE154989_mmLungPlate_fQC_dSp_rawCount.h5"
retry GSE154989_mmLungPlate_fQC_geneTable.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/GSE154989_mmLungPlate_fQC_geneTable.csv.gz"
retry GSE154989_mmLungPlate_fQC_smpTable.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/GSE154989_mmLungPlate_fQC_smpTable.csv.gz"
retry GSE154989_mmLungPlate_fQC_dZ_annot_smpTable.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/GSE154989_mmLungPlate_fQC_dZ_annot_smpTable.csv.gz"

retry GSE154977_mmLung10x_cis_dSp_rawCount.h5 \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl/GSE154977_mmLung10x_cis_dSp_rawCount.h5"
retry GSE154977_mmLung10x_cis_geneTable.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl/GSE154977_mmLung10x_cis_geneTable.csv.gz"
retry GSE154977_mmLung10x_cis_smpTable.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl/GSE154977_mmLung10x_cis_smpTable.csv.gz"
retry GSE154977_mmLung10x_cis_dZ_annot_annot_smpTable.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl/GSE154977_mmLung10x_cis_dZ_annot_annot_smpTable.csv.gz"

retry GSE179501_XTR_scRNAseq_matrix.mtx.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179501/suppl/GSE179501_XTR_scRNAseq_matrix.mtx.gz"
retry GSE179501_XTR_scRNAseq_features.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179501/suppl/GSE179501_XTR_scRNAseq_features.tsv.gz"
retry GSE179501_XTR_scRNAseq_barcodes.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179501/suppl/GSE179501_XTR_scRNAseq_barcodes.tsv.gz"

retry GSE179502_XTR_sorted_scRNAseq_matrix.mtx.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179502/suppl/GSE179502_XTR_sorted_scRNAseq_matrix.mtx.gz"
retry GSE179502_XTR_sorted_scRNAseq_features.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179502/suppl/GSE179502_XTR_sorted_scRNAseq_features.tsv.gz"
retry GSE179502_XTR_sorted_scRNAseq_barcodes.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179502/suppl/GSE179502_XTR_sorted_scRNAseq_barcodes.tsv.gz"

retry GSE165641_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165641/suppl/GSE165641_RAW.tar"
retry GSE180963_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180963/suppl/GSE180963_RAW.tar"
retry GSE149813_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149813/suppl/GSE149813_RAW.tar"

retry GSE267321_Normalized_expression_matrix_02122026.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267321/suppl/GSE267321_Normalized_expression_matrix_02122026.csv.gz"

retry GSE127465_mouse_counts_normalized_15939x28205.mtx.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/GSE127465_mouse_counts_normalized_15939x28205.mtx.gz"
retry GSE127465_gene_names_mouse_28205.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/GSE127465_gene_names_mouse_28205.tsv.gz"
retry GSE127465_mouse_cell_metadata_15939x12.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/GSE127465_mouse_cell_metadata_15939x12.tsv.gz"

retry GSE136246_mouse_normalized_gene_counts.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136246/suppl/GSE136246_mouse_normalized_gene_counts.txt.gz"
retry GSE136246_mouse_cell_groupings.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136246/suppl/GSE136246_mouse_cell_groupings.csv.gz"

echo "download done -> $OUT"
