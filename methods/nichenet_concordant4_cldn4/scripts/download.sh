#!/usr/bin/env bash
# Public processed concordant-four matrices + NicheNet v2 prior.
# Does not download GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, or GSE207422.
set -euo pipefail
OUT="${1:-/tmp/concordant4_raw}"
PRI="${2:-/tmp/nichenet_prior}"
mkdir -p "$OUT"/GSE123902 "$OUT"/GSE131907 "$OUT"/GSE205335 "$OUT"/GSE189357 "$PRI"

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

dl "$PRI/ligand_target_matrix_nsga2r_final.rds" \
  "https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final.rds?download=1"
dl "$PRI/lr_network_human_21122021.rds" \
  "https://zenodo.org/records/7074291/files/lr_network_human_21122021.rds?download=1"

dl "$OUT/GSE123902/GSE123902_RAW.tar" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar
dl "$OUT/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz
dl "$OUT/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz
dl "$OUT/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz
dl "$OUT/GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz
dl "$OUT/GSE189357/GSE189357_RAW.tar" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar

# NicheNet source used by activity.R, pinned.
SHA=66f90d5eeafef280b2b2f339b3fd70ffec1781dd
SRC="${3:-/tmp/nichenet_src}"
mkdir -p "$SRC"
for f in supporting_functions.R evaluate_model_target_prediction.R evaluate_model_ligand_prediction.R; do
  dl "$SRC/$f" "https://raw.githubusercontent.com/saeyslab/nichenetr/${SHA}/R/${f}"
done

echo "refused: GSE148071 GSE127465 GSE154826 GSE200563 E-MTAB-13526 GSE207422"
