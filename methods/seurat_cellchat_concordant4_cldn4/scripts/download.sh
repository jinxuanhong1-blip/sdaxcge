#!/usr/bin/env bash
# Public processed concordant-four matrices. No GSE148071. No FASTQ/EGA/SRA.
# Prefer MTX / dense CSV / UMI TXT. GSE205335 has only an RDS on GEO.
set -euo pipefail
OUT="${1:-/tmp/concordant4_raw}"
mkdir -p "$OUT"/GSE123902 "$OUT"/GSE131907 "$OUT"/GSE205335 "$OUT"/GSE189357

dl() {
  local dest="$1" url="$2"
  if [[ -s "$dest" ]] && [[ "$(stat -c%s "$dest")" -gt 1000 ]]; then
    echo "HAVE $dest ($(stat -c%s "$dest") bytes)"
    return 0
  fi
  echo "GET $url"
  curl -fL --retry 6 --retry-delay 8 --retry-all-errors -o "${dest}.partial" "$url"
  mv "${dest}.partial" "$dest"
  echo "OK $dest ($(stat -c%s "$dest") bytes)"
}

# GSE123902: Laughney dense UMI CSVs (90 MB). Skip 36.5 GB H5.
dl "$OUT/GSE123902/GSE123902_RAW.tar" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar

# GSE131907: raw UMI TXT + author annotation. Skip 2.9 GB log2TPM text and both RDS.
dl "$OUT/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz
dl "$OUT/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz

# GSE205335: GEO ships only an RDS UMI matrix. Identity + SOFT are TXT.
dl "$OUT/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz
dl "$OUT/GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz
dl "$OUT/GSE205335/GSE205335_family.soft.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz

# GSE189357: 10x MTX/TSV TAR.
dl "$OUT/GSE189357/GSE189357_RAW.tar" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar

echo "refused: GSE148071 GSE127465 7-pool dual-high FASTQ EGA SRA 36.5GB H5 2.9GB log2TPM"
ls -lh "$OUT"/*/*
