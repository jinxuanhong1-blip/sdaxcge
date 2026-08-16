#!/usr/bin/env bash
# Reproduce the LUSC / ICI TACSTD2-CLDN4 slice from public GEO + Xena files.
set -euo pipefail
cd "$(dirname "$0")"
python3 01_download.py
python3 02_build_cohorts.py
python3 03_tcga_reference.py
python3 04_histology.py
python3 05_analyze.py
python3 06_figures.py
python3 07_size_check.py
