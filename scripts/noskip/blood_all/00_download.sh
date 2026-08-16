#!/usr/bin/env bash
# Public GEO files used by 03_analyze.py / 04_analyze_gse285888.py.
# Not committed. Default dest: /tmp/geo_blood
set -euo pipefail
DEST="${GEO_DIR:-/tmp/geo_blood}"
mkdir -p "$DEST" "$DEST/gse285888"
cd "$DEST"

get() { local out="$1" url="$2"; if [[ -s "$out" ]]; then echo "have $out"; else wget -q -c -O "$out" "$url"; echo "got $out"; fi; }

get GSE111414_gene_counts.csv.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111414/suppl/GSE111414_gene_counts.csv.gz
get GSE111414_series_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111414/matrix/GSE111414_series_matrix.txt.gz
get GSE202417_series_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202417/matrix/GSE202417_series_matrix.txt.gz
get GSE249262_series_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249262/matrix/GSE249262_series_matrix.txt.gz
get GSE235048_24CFlow_TPM.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235048/suppl/GSE235048_24CFlow_TPM.txt.gz
get GSE235048_series_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235048/matrix/GSE235048_series_matrix.txt.gz
get GSE225620_pre_vs_post_featureCounts.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225620/suppl/GSE225620_pre_vs_post_featureCounts.txt.gz
get GSE225620_series_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225620/matrix/GSE225620_series_matrix.txt.gz
get GSE305086_Expression_matrix_final.csv.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305086/suppl/GSE305086_Expression_matrix_final.csv.gz
get GSE305086_series_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305086/matrix/GSE305086_series_matrix.txt.gz
get GSE152590_TPM_matrix.xlsx https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152590/suppl/GSE152590_TPM_matrix.xlsx
get GSE216297_TEP_Count_Matrix.RData.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE216nnn/GSE216297/suppl/GSE216297_TEP_Count_Matrix.RData.gz
get GSE216297_series_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE216nnn/GSE216297/matrix/GSE216297_series_matrix.txt.gz
get GSE285888_RAW.tar https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285888/suppl/GSE285888_RAW.tar
if [[ ! -s gse285888/GSM8712033_matrix.txt.gz ]]; then
  tar -xf GSE285888_RAW.tar -C gse285888
fi
echo done
