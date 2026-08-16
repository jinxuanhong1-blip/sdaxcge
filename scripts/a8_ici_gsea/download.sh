#!/usr/bin/env bash
# Download processed public ICI lung-tumor bulk matrices (no FASTQ, <2 GB).
# Leftover of the TCGA-only A8 slide: pre-treatment open GEO series with TACSTD2.
set -euo pipefail
DEST="${1:-/tmp/a8_ici_gsea_data}"
mkdir -p "$DEST"
cd "$DEST"

dl() {
  url="$1"
  out="$2"
  if [ -s "$out" ]; then
    echo "cached $out ($(wc -c < "$out") bytes)"
    return 0
  fi
  for i in 1 2 3 4; do
    if curl -fsSL --retry 2 -A "a8-ici-gsea/1.0" -o "$out.tmp" "$url" && [ -s "$out.tmp" ]; then
      mv "$out.tmp" "$out"
      echo "OK  $out ($(wc -c < "$out") bytes)"
      return 0
    fi
    echo "retry $i $out"
    sleep $((2 ** i))
  done
  echo "FAIL $out from $url"
  return 1
}

# --- GSE126044 Cho 2019 anti-PD-1 NSCLC counts ---
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz" GSE126044_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz" GSE126044_counts.txt.gz

# --- GSE135222 Jung/Kim 2019–2020 anti-PD-(L)1 NSCLC TPM ---
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz" GSE135222_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz" GSE135222_exp.tsv.gz

# --- GSE166449 Hwang 2021 pembrolizumab LUAD TPM ---
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/matrix/GSE166449_series_matrix.txt.gz" GSE166449_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz" GSE166449_TPM.txt.gz

# --- GSE253564 Altorki leftover neoadjuvant durvalumab pre-treatment FPKM ---
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/matrix/GSE253564_series_matrix.txt.gz" GSE253564_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz" GSE253564_pre_FPKM.txt.gz

# --- GSE190265 France3 anti-PD-1 NSCLC TPM (France4 lacks TACSTD2) ---
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/matrix/GSE190265_series_matrix.txt.gz" GSE190265_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/suppl/GSE190265_TPM_France3.csv.gz" GSE190265_TPM_France3.csv.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190265/suppl/GSE190265_samples_info_France3.csv.gz" GSE190265_samples_info_France3.csv.gz

# --- GSE283829 Lindberg 2025 leftover ICI NSCLC raw counts ---
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283829/matrix/GSE283829_series_matrix.txt.gz" GSE283829_series_matrix.txt.gz
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283829/suppl/GSE283829_raw_express_matrix_all_samples.txt.gz" GSE283829_raw_counts.txt.gz

# Ensembl -> HGNC for GSE135222 / GSE283829
dl "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.tsv" hgnc_complete_set.tsv

echo "Done. Files in $DEST"
ls -lh "$DEST"
