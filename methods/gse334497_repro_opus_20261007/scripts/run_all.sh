#!/usr/bin/env bash
# Full GSE334497 reanalysis, in dependency order. Outputs go to tables/ and figures/.
set -euo pipefail
cd "$(dirname "$0")"
bash 00_fetch.sh
python3 01_prepare_qc.py
Rscript 02_de.R
python3 03_scores.py
python3 04_gsea.py
python3 05_related_datasets.py
python3 07_tissue_contamination.py
python3 06_figures.py
