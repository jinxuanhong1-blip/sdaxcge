#!/usr/bin/env bash
# Download the GSE271689 GeoMx DSP-WTA demo data (NSCLC, PD-1 immunotherapy).
# 下载 GSE271689 GeoMx DSP-WTA 演示数据（NSCLC，PD-1 免疫治疗）。
#
# The GEO supplementary archive GSE271689_RAW.tar is ~36 MB and contains 586
# per-AOI .dcc.gz count files -- well under the 2 GB processed-data budget.
# GEO 补充档 GSE271689_RAW.tar 约 36 MB，含 586 个 AOI 的 .dcc.gz 计数文件，
# 远小于 2 GB 的处理数据预算。
#
# NOT included in the GEO archive, and required to build a GeoMxSet:
# GEO 档案中不含、但构建 GeoMxSet 所必需：
#   1) The WTA probe kit config  Hs_R_NGS_WTA_v1.0.pkc  (from NanoString).
#      WTA 探针配置文件（来自 NanoString）。
#   2) An AOI annotation sheet mapping each DCC id -> slide/ROI/compartment/
#      patient (+ clinical OS columns for template 09). Reconstruct from the
#      GEO series matrix / sample characteristics.
#      将每个 DCC id 映射到 切片/ROI/区室/患者（+模板09所需临床OS列）的注释表；
#      可由 GEO series matrix / 样本特征重建。
#
# Usage: bash methods/spatial/demo/download_geomx_gse271689.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${HERE}/data/GSE271689"
mkdir -p "${DEST}"

TAR="${HERE}/data/GSE271689_RAW.tar"
URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE271689&format=file"

if [[ ! -f "${TAR}" ]]; then
  echo "[download] fetching GSE271689_RAW.tar (~36 MB) ..."
  curl -fSL -o "${TAR}" "${URL}"
fi

echo "[download] extracting DCC files ..."
tar -xf "${TAR}" -C "${DEST}"

N=$(find "${DEST}" -name '*.dcc.gz' | wc -l | tr -d ' ')
echo "[download] extracted ${N} .dcc.gz files into ${DEST}"
echo "[next] obtain Hs_R_NGS_WTA_v1.0.pkc and an AOI annotation sheet, then run:"
echo "       Rscript methods/spatial/templates/06_geomx_qc_normalization.R \\"
echo "           --dcc ${DEST} --pkc <path/Hs_R_NGS_WTA_v1.0.pkc> \\"
echo "           --annotation <path/geomx_annotation.csv> --out methods/spatial/demo/out/geomx"
