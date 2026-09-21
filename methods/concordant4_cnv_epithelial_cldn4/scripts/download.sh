#!/usr/bin/env bash
# Public concordant-4 matrices only. No GSE148071. No FASTQ / EGA / SRA.
# GSE205335 is the double-gzipped GEO RDS (the only UMI matrix on GEO).
set -euo pipefail
OUT="${1:-/tmp/concordant4_raw}"
mkdir -p "$OUT"
dl() {
  local dest="$1" url="$2"
  if [[ -s "$dest" ]]; then
    echo "HAVE $dest"
    return 0
  fi
  echo "GET $url"
  curl -fL --retry 6 --retry-delay 5 -o "${dest}.partial" "$url"
  mv "${dest}.partial" "$dest"
  echo "OK $dest"
}
dl "$OUT/GSE123902_RAW.tar" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar
dl "$OUT/GSE131907_ann.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz
dl "$OUT/GSE131907_umi.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz
dl "$OUT/GSE205335_ident.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz
dl "$OUT/GSE205335_family.soft.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz
dl "$OUT/GSE205335_umi.rds.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz
dl "$OUT/GSE189357_RAW.tar" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar
dl "$OUT/refGene.txt.gz" \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
echo "refused: GSE148071 GSE127465 GSE207422 GSE154826 GSE200563 E-MTAB-13526"
