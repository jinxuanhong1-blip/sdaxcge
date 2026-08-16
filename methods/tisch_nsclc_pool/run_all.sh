#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
python3 methods/tisch_nsclc_pool/download.py
python3 methods/tisch_nsclc_pool/analyze.py
python3 methods/tisch_nsclc_pool/pool.py
python3 methods/tisch_nsclc_pool/write_report.py
