#!/usr/bin/env bash
# Download every public file used by this slice. Raw data stay under
# $OPUS_TLS_DATA (default /tmp/opus_tls_data) and are NOT committed.
# Patil 2019 / OAK (EGAD00001002517, EGAD00001008548) are EGA-restricted
# and are listed only — see notes/opus_tls/DATA_ACCESS.md.
set -euo pipefail
ROOT="${OPUS_TLS_DATA:-/tmp/opus_tls_data}"
mkdir -p "$ROOT"/{tcga,geo,meta,processed}

XENA=https://gdc-hub.s3.us-east-1.amazonaws.com/download
GDC=https://api.gdc.cancer.gov/data
GEO=https://ftp.ncbi.nlm.nih.gov/geo

dl() {  # dest url
  local dest="$1" url="$2"
  if [[ -s "$dest" ]]; then
    echo "[skip] $dest"
    return
  fi
  echo "[get]  $dest"
  curl -fsSL --retry 4 --retry-delay 4 -o "$dest" "$url"
}

echo "== TCGA Xena GDC hub =="
for f in \
  TCGA-LUAD.star_tpm.tsv.gz TCGA-LUSC.star_tpm.tsv.gz \
  TCGA-LUAD.clinical.tsv.gz TCGA-LUSC.clinical.tsv.gz \
  TCGA-LUAD.survival.tsv.gz TCGA-LUSC.survival.tsv.gz \
  gencode.v36.annotation.gtf.gene.probemap
do
  dl "$ROOT/tcga/$f" "$XENA/$f"
done

echo "== GDC PanCanAtlas immune / purity =="
dl "$ROOT/meta/absolute_purity.txt"          "$GDC/4f277128-f793-4354-a13d-30cc7fe9f6b5"
dl "$ROOT/meta/leukocyte_fraction.tsv"       "$GDC/6f75c9d7-5134-4ed1-b8f3-72856c98a4e8"
dl "$ROOT/meta/cibersort_relative.tsv"       "$GDC/b3df502e-3594-46ef-9f94-d041a20a0b9a"
dl "$ROOT/meta/mutation_load.txt"            "$GDC/ff3f962c-3573-44ae-a8f4-e5ac0aea64b6"

echo "== GEO expression + series matrices =="
# ICI
dl "$ROOT/geo/GSE135222_exp.tsv.gz"          "$GEO/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
dl "$ROOT/geo/GSE135222_matrix.txt.gz"       "$GEO/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz"
dl "$ROOT/geo/GSE126044_counts.txt.gz"       "$GEO/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz"
dl "$ROOT/geo/GSE126044_matrix.txt.gz"       "$GEO/series/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz"
dl "$ROOT/geo/GSE207422_bulk_log2TPM.txt.gz" "$GEO/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz"
dl "$ROOT/geo/GSE207422_bulk_meta.xlsx"      "$GEO/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx"
dl "$ROOT/geo/GSE207422_matrix.txt.gz"       "$GEO/series/GSE207nnn/GSE207422/matrix/GSE207422_series_matrix.txt.gz"
# Atlases
dl "$ROOT/geo/GSE81089_fpkm.tsv.gz"          "$GEO/series/GSE81nnn/GSE81089/suppl/GSE81089_FPKM_cufflinks.tsv.gz"
dl "$ROOT/geo/GSE81089_matrix.txt.gz"        "$GEO/series/GSE81nnn/GSE81089/matrix/GSE81089_series_matrix.txt.gz"
dl "$ROOT/geo/GSE72094_matrix.txt.gz"        "$GEO/series/GSE72nnn/GSE72094/matrix/GSE72094_series_matrix.txt.gz"
dl "$ROOT/geo/GPL15048.txt"                  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15048&targ=self&form=text&view=data"
dl "$ROOT/geo/GSE190265_tpm.csv.gz"          "$GEO/series/GSE190nnn/GSE190265/suppl/GSE190265_TPM_France3.csv.gz"
dl "$ROOT/geo/GSE190265_matrix.txt.gz"       "$GEO/series/GSE190nnn/GSE190265/matrix/GSE190265_series_matrix.txt.gz"

echo "[done] public downloads in $ROOT"
echo "[note] Patil/OAK (EGA) were NOT downloaded. See notes/opus_tls/DATA_ACCESS.md"
