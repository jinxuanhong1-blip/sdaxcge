#!/usr/bin/env bash
set -u
BASE="https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files"
OUT="data/emtab13530"
mkdir -p "$OUT"
# grab SDRF for sample annotation
curl -sS -m 120 -o "$OUT/E-MTAB-13530.sdrf.txt" "$BASE/E-MTAB-13530.sdrf.txt"
mapfile -t SAMPLES < <(grep 'filtered_feature_bc_matrix.h5' /tmp/emtab_files.txt | sed 's/-filtered_feature_bc_matrix.h5//')
echo "found ${#SAMPLES[@]} samples"
fetch() {
  local s="$1"
  local h5="$OUT/${s}-filtered_feature_bc_matrix.h5"
  local sp="$OUT/${s}-spatial.tar"
  [ -s "$h5" ] || curl -sS -m 300 -o "$h5" "$BASE/${s}-filtered_feature_bc_matrix.h5"
  [ -s "$sp" ] || curl -sS -m 300 -o "$sp" "$BASE/${s}-spatial.tar"
  echo "done $s"
}
export -f fetch; export BASE OUT
printf "%s\n" "${SAMPLES[@]}" | xargs -P 6 -I {} bash -c 'fetch "$@"' _ {}
echo EMTAB_ALL_DONE
