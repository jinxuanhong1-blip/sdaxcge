#!/usr/bin/env bash
# Reproduce the full slice. Total download < 2 GB (processed data only).
# SCLC_DATA_DIR (default /tmp/sclc_data) holds downloads; the repository
# only receives small tables and figures under results/fable_sclc/.
set -euo pipefail
cd "$(dirname "$0")"

python3 01_download_data.py "$@"
python3 02_bulk_subtypes.py
python3 03_geomx_ici.py
python3 04_scrna_atlas.py
