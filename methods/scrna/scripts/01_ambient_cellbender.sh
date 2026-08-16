#!/usr/bin/env bash
# Ambient-RNA removal with CellBender, per capture lane / 每个测序lane单独运行。
# Input MUST be the RAW (unfiltered) 10x .h5 (all barcodes), not filtered.
# 输入必须是未过滤的 raw .h5。
#
# Install: pip install cellbender    (GPU strongly recommended / 建议GPU)
# Docs:    https://cellbender.readthedocs.io
set -euo pipefail

RAW_H5=${1:?"usage: 01_ambient_cellbender.sh raw_feature_bc_matrix.h5 out_prefix [expected_cells] [total_droplets]"}
OUT=${2:?"output prefix"}
EXPECTED_CELLS=${3:-5000}
TOTAL_DROPLETS=${4:-20000}

cellbender remove-background \
    --input "${RAW_H5}" \
    --output "${OUT}_cellbender.h5" \
    --expected-cells "${EXPECTED_CELLS}" \
    --total-droplets-included "${TOTAL_DROPLETS}" \
    --fpr 0.01 \
    --epochs 150 \
    --cuda

# QC: inspect ${OUT}_cellbender.pdf (learning curve, cell-probability, ambient profile).
# Downstream: load ${OUT}_cellbender_filtered.h5 with scanpy.read_10x_h5(..., gex_only=True).
# 校验：查看输出PDF；下游读取 *_filtered.h5。
echo "CellBender done -> ${OUT}_cellbender_filtered.h5"
