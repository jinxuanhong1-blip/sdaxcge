#!/usr/bin/env bash
# Reproduce the full fable_stk11 slice end-to-end.
# All outputs land under results/fable_stk11/ (data, tables, figures).
set -euo pipefail
cd "$(dirname "$0")"

echo "[1/5] fetch public data (cBioPortal + GEO)"
python3 01_fetch_data.py

echo "[2/5] TCGA: genotype <-> TACSTD2/CLDN4 expression"
python3 02_analyze_tcga.py

echo "[3/5] ICI cohorts: genotype <-> response (DCB/PFS)"
python3 03_analyze_ici_genotype.py

echo "[4/5] GEO GSE135222: expression <-> response"
python3 04_analyze_geo_expression.py

echo "[5/5] verification"
python3 05_verify.py

echo "DONE. See results/fable_stk11/ and notes/fable_stk11/WRITEUP.md"
