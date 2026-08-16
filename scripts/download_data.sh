#!/usr/bin/env bash
# Download the tumor ("T") 10x Visium sections of E-MTAB-13530
# (human NSCLC lesions and non-involved tissue) from BioStudies/ArrayExpress.
#
# For each section we fetch:
#   <sample>-filtered_feature_bc_matrix.h5   gene x spot count matrix (in-tissue)
#   <sample>-spatial.tar                     tissue_positions_list.csv + scalefactors
#
# Data are written to data/E-MTAB-13530/ and are NOT committed (see .gitignore).
set -euo pipefail

DEST="data/E-MTAB-13530"
BASE="https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files"
mkdir -p "$DEST"

SAMPLES="P10_T1 P10_T2 P10_T3 P10_T4 \
P11_T1 P11_T2 P11_T3 P11_T4 \
P15_T1 P15_T2 P16_T1 P16_T2 \
P17_T1 P17_T2 P19_T1 P19_T2 \
P24_T1 P24_T2 P25_T1 P25_T2"

for s in $SAMPLES; do
  for suf in filtered_feature_bc_matrix.h5 spatial.tar; do
    f="${s}-${suf}"
    if [ ! -s "$DEST/$f" ]; then
      echo "downloading $f"
      curl -sS -m 300 --retry 4 -o "$DEST/$f" "$BASE/$f"
    fi
  done
done
echo "done -> $DEST"
