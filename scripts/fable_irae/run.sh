#!/usr/bin/env bash
# Reproduce the fable_irae slice end-to-end.
# Requires: python3 with pandas, numpy, scipy, statsmodels, matplotlib, anndata, h5py.
#   pip install --user pandas numpy scipy statsmodels matplotlib anndata h5py GEOparse
# Large matrices (~2 GB) are downloaded to results/fable_irae/data/raw/ (gitignored).
set -euo pipefail
cd "$(dirname "$0")"

python3 01_discover_geo.py          # GEO search -> candidate tables
python3 02_download.py              # download processed matrices (<2 GB)
gunzip -kf ../../results/fable_irae/data/raw/GSE206300_ircolitis-tissue-epithelial.h5ad.gz
python3 03_inspect_h5ad.py          # structural sanity check of the two h5ad files
python3 04_analyze_colon.py         # GSE206300 colon epithelium (colitis)
python3 05_analyze_blood.py         # GSE319496 whole-blood bulk (irAE Y/N)
python3 06_analyze_balf.py          # GSE277136 BALF (pneumonitis, lung)
python3 07_verify.py                # controls, robustness, master summary + JSON verdict
python3 08_summary_figure.py        # cross-dataset measurability figure
echo "Done. See results/fable_irae/tables/master_summary.tsv and notes/fable_irae/WRITEUP.md"
