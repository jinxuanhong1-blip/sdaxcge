#!/usr/bin/env bash
# Download processed matrices (+ series_matrix metadata) for the ICI bulk lung analysis.
# Processed data only, no FASTQ. All files < 2 GB. Staged in /tmp/ici_bulk_data.
set -euo pipefail
DEST="${1:-/tmp/ici_bulk_data}"
mkdir -p "$DEST"; cd "$DEST"

dl(){ url="$1"; out="$2"; for i in 1 2 3 4; do
  curl -sSL -o "$out" "$url" && [ -s "$out" ] && { echo "OK  $out"; return 0; }
  echo "retry $i $out"; sleep $((2**i)); done; echo "FAIL $out"; return 1; }

# series_matrix metadata (response labels)
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz" GSE126044_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz" GSE135222_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136961/matrix/GSE136961_series_matrix.txt.gz" GSE136961_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/matrix/GSE166449_series_matrix.txt.gz" GSE166449_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93157/matrix/GSE93157_series_matrix.txt.gz"   GSE93157_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/matrix/GSE207422_series_matrix.txt.gz" GSE207422_series_matrix.txt.gz

# processed expression matrices
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz" GSE126044_counts.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz" GSE135222_exp.tsv.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136961/suppl/GSE136961_TPM.tsv.gz" GSE136961_TPM.tsv.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz" GSE166449_TPM.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93157/suppl/GSE93157_raw_data_values.txt.gz" GSE93157_raw.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz" GSE207422_log2TPM.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx" GSE207422_metadata.xlsx

# round-2 matching cohorts: neoadjuvant durvalumab NSCLC (correlation-only)
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/matrix/GSE253564_series_matrix.txt.gz" GSE253564_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/matrix/GSE248378_series_matrix.txt.gz" GSE248378_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz" GSE253564_pre_FPKM.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/suppl/GSE248378_Durva_Post_FPKMs.txt.gz" GSE248378_post_FPKM.txt.gz

echo "Done. Files in $DEST"
echo "Note: cbioportal_analysis.py fetches TCGA-LUAD/LUSC + OncoSG live from cbioportal.org."
