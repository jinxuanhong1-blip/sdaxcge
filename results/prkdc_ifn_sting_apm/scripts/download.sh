#!/usr/bin/env bash
# Fetch the public GEO matrices scored by score_meta.py.
# Usage: PRKDC_DATA=/tmp/geo/data bash scripts/download.sh
set -euo pipefail
DEST="${PRKDC_DATA:-/tmp/geo/data}"
mkdir -p "$DEST"
cd "$DEST"

fetch() {
  local url="$1"
  local name
  name="$(basename "$url")"
  if [[ -s "$name" ]]; then
    echo "have $name"
    return
  fi
  echo "get $name"
  curl -fL --retry 3 -o "$name" "$url"
}

fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273409/suppl/GSE273409_salmon_gene_counts_normalized.csv.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242255/suppl/GSE242255_AREK_DNAPK_DESeq2_normalised_counts.csv.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116765/suppl/GSE116765_DNAPK_AllTreats_countTable.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116765/suppl/GSE116765_22Rv1_DNAPKi_CountTable.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167956/suppl/GSE167956_RAW_COUNTS-gene_expressions_in_24_samples.xlsx"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180581/suppl/GSE180581_all_samples_counts.xlsx"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268415/suppl/GSE268415_gene_count.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285698/suppl/GSE285698_raw_counts.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129436/suppl/GSE129436_norm_matrix.csv.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE315nnn/GSE315862/suppl/GSE315862_processed_counts.xlsx"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319513/suppl/GSE319513_KU_24HVSCON_24H_Gene_differential_expression.xlsx"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319513/suppl/GSE319513_KU_48HVSCON_48H_Gene_differential_expression.xlsx"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93074/suppl/GSE93074_non-normalized-nu7441.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93075/suppl/GSE93075_non-normalized-siPRKDCTC32.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85202/suppl/GSE85202_non-normalized.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63480/matrix/GSE63480_series_matrix.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz"
fetch "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt"
fetch "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
fetch "https://bioconductor.org/packages/release/data/annotation/src/contrib/hugene20sttranscriptcluster.db_8.8.0.tar.gz"

if [[ ! -d hugene20sttranscriptcluster.db ]]; then
  tar -xzf hugene20sttranscriptcluster.db_8.8.0.tar.gz
fi
echo "ready in $DEST"
