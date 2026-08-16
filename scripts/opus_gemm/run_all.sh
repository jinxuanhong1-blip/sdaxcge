#!/usr/bin/env bash
# Reproduce the KRAS/LKB1 + NSCLC GEMM ICI Tacstd2/Cldn4 slice.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

python3 scripts/opus_gemm/genes.py --out results/opus_gemm/gene_ids.tsv
python3 scripts/opus_gemm/fetch_data.py
python3 scripts/opus_gemm/build_matrices.py
python3 scripts/opus_gemm/analyze.py
python3 scripts/opus_gemm/plot.py
python3 scripts/opus_gemm/write_verification.py
echo "done"
