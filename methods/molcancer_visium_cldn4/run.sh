#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
python3 methods/molcancer_visium_cldn4/count_fastq.py
python3 methods/molcancer_visium_cldn4/analyze.py
