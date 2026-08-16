#!/usr/bin/env bash
# Reproduce the full fable_paired ICI-lung TACSTD2/CLDN4 analysis end to end.
#
# Raw GEO downloads land in $FABLE_RAW_DIR (default /tmp/fable_raw), OUTSIDE the
# repo, so committed processed artifacts stay well under the 2 GB budget.
set -euo pipefail
cd "$(dirname "$0")"

export FABLE_RAW_DIR="${FABLE_RAW_DIR:-/tmp/fable_raw}"
echo "Raw data dir: $FABLE_RAW_DIR"

python3 -m pip install --quiet -r requirements.txt

python3 download_data.py
python3 01_scrna_gse207422.py       # core: pre vs post, epithelial, response
python3 02_baseline_bulk.py         # baseline verification (3 NSCLC cohorts)
python3 03_paired_pre_on_gse91061.py  # orthogonal within-patient pre->on dynamics
python3 06_paired_gse248249.py      # A7 analog: same-patient NSCLC pre vs acquired resistance
python3 07_mouse_gse246922.py       # TISMO stand-in: mouse lung ICB lines
python3 04_summary.py               # master table + effect-size forest + FDR
python3 05_verify.py                # positive-control sanity checks

echo "Done. See results/fable_paired/ for tables and figures."
