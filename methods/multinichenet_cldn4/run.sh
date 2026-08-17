#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 scripts/00_download.py
python3 scripts/01_convert_prior.py
python3 scripts/02_extract_panels.py
python3 scripts/03_analyze.py
