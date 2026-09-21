#!/bin/bash
# Public GEO files for the mouse-unit mediation. Nothing here is the private 8 KL set.
set -euo pipefail
mkdir -p /tmp/geo_dl
cd /tmp/geo_dl
curl -fL --retry 4 -o GSE137244_counts.fpkm.csv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/GSE137244_counts.fpkm.csv.gz"
curl -fL --retry 4 -o GSE264739_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264739/suppl/GSE264739_RAW.tar"
curl -fL --retry 4 -o GSE295824_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295824/suppl/GSE295824_RAW.tar"
curl -fL --retry 4 -o GSE201247_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201247/suppl/GSE201247_RAW.tar"
curl -fL --retry 4 -o GSE165641_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165641/suppl/GSE165641_RAW.tar"
curl -fL --retry 4 -o GSE180963_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180963/suppl/GSE180963_RAW.tar"
curl -fL --retry 4 -o GSE179501_matrix.mtx.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179501/suppl/GSE179501_XTR_scRNAseq_matrix.mtx.gz"
curl -fL --retry 4 -o GSE179501_features.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179501/suppl/GSE179501_XTR_scRNAseq_features.tsv.gz"
curl -fL --retry 4 -o GSE179501_barcodes.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179501/suppl/GSE179501_XTR_scRNAseq_barcodes.tsv.gz"
for s in d10_1 d10_2 dTom_1 dTom_2; do
  for k in matrix.mtx features.tsv barcodes.tsv; do
    curl -fL --retry 4 -o "GSE266323_${s}_${k}.gz" \
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266323/suppl/GSE266323_${s}_MMT_${k}.gz"
  done
done
