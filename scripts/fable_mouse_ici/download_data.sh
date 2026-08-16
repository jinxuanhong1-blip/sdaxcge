#!/usr/bin/env bash
# Download processed supplementary files for the mouse-lung ICI slice into a
# temp cache (NOT committed to the repo; only derived tables/plots are).
# All target files are < 2 GB (verified in notes/file_inventory.json), so none
# are skipped for size.
set -euo pipefail
D=/tmp/fable_ici_data
mkdir -p "$D"
cd "$D"
UA="Mozilla/5.0 (fable-mouse-ici)"

dl() { # url outname
  if [ -f "$2" ]; then echo "cached: $2"; return; fi
  echo "GET $1"
  curl -fSL --retry 4 --retry-delay 3 -A "$UA" -o "$2" "$1"
}

B=https://ftp.ncbi.nlm.nih.gov/geo/series
dl $B/GSE239nnn/GSE239485/suppl/GSE239485_Processed_data.xlsx GSE239485_Processed_data.xlsx
dl $B/GSE241nnn/GSE241978/suppl/GSE241978_2020-07-21_Sherr_analysis_CMT_KO_vs_Cas9Ctrl.xlsx GSE241978.xlsx
dl $B/GSE197nnn/GSE197260/suppl/GSE197260_RNAseqTPM_MM_EGFR-TKI-CD8.txt.gz GSE197260_TPM.txt.gz
dl $B/GSE297nnn/GSE297630/suppl/GSE297630_processed_data.xlsx GSE297630_processed_data.xlsx
dl $B/GSE297nnn/GSE297630/suppl/GSE297630_RAW.tar GSE297630_RAW.tar
dl $B/GSE133nnn/GSE133604/suppl/GSE133604_RAW.tar GSE133604_RAW.tar
dl $B/GSE129nnn/GSE129297/suppl/GSE129297_RAW.tar GSE129297_RAW.tar
dl $B/GSE222nnn/GSE222158/suppl/GSE222158_RAW.tar GSE222158_RAW.tar
dl $B/GSE330nnn/GSE330658/suppl/GSE330658_RAW.tar GSE330658_RAW.tar
dl $B/GSE297nnn/GSE297632/suppl/GSE297632_RAW.tar GSE297632_RAW.tar

# ArrayExpress / BioStudies
dl https://www.ebi.ac.uk/biostudies/files/E-MTAB-13704/GEMMS_raw_counts.csv E-MTAB-13704_GEMMS_raw_counts.csv
dl https://www.ebi.ac.uk/biostudies/files/E-MTAB-13704/E-MTAB-13704.sdrf.txt E-MTAB-13704.sdrf.txt

echo "=== downloaded ==="
ls -lh "$D"
