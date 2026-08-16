#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
python3 scripts/w200/C4_GSE334497/analyze.py
