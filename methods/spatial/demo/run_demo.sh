#!/usr/bin/env bash
# End-to-end GeoMx demo runner for GSE271689 (guarded, no fabricated outputs).
# GSE271689 的端到端 GeoMx 演示运行器（带守卫，不伪造任何输出）。
#
# Behaviour: each prerequisite is checked; if anything required is missing the
# script explains exactly what to supply and exits WITHOUT producing fake
# results. It only runs a stage when its real inputs and tools are present.
# 行为：逐项检查前置条件；若缺任何必需项，脚本明确说明需补充什么并退出，绝不产生
# 伪结果。仅当某阶段的真实输入与工具齐备时才运行该阶段。
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
DCC="${HERE}/data/GSE271689"
PKC="${PKC:-${HERE}/data/pkc/Hs_R_NGS_WTA_v1.0.pkc}"
ANNOT="${ANNOT:-${HERE}/data/geomx_annotation.csv}"
OUT="${HERE}/out/geomx"
T="${ROOT}/methods/spatial/templates"

echo "== GeoMx demo (GSE271689) =="

# 1) Data present? 数据是否就绪？
if ! ls "${DCC}"/*.dcc.gz >/dev/null 2>&1 && ! ls "${DCC}"/*.dcc >/dev/null 2>&1; then
  echo "[step] downloading demo DCC files ..."
  bash "${HERE}/download_geomx_gse271689.sh"
fi
N=$(find "${DCC}" -name '*.dcc*' 2>/dev/null | wc -l | tr -d ' ')
echo "[ok] found ${N} DCC files."

MISSING=0
# 2) Tooling present? 工具是否安装？
if ! command -v Rscript >/dev/null 2>&1; then
  echo "[missing] Rscript not found. Install R >= 4.3, then:"
  echo "          Rscript methods/spatial/environment/install-geomx.R"
  MISSING=1
else
  if ! Rscript -e 'quit(status = as.integer(!requireNamespace("GeomxTools", quietly=TRUE)))' 2>/dev/null; then
    echo "[missing] R package GeomxTools. Install with:"
    echo "          Rscript methods/spatial/environment/install-geomx.R"
    MISSING=1
  fi
fi

# 3) PKC + annotation present? These are NOT in the GEO archive. PKC 与注释是否就绪？
if [[ ! -f "${PKC}" ]]; then
  echo "[missing] WTA PKC: ${PKC}"
  echo "          Obtain Hs_R_NGS_WTA_v1.0.pkc from NanoString and place it there (or set \$PKC)."
  MISSING=1
fi
if [[ ! -f "${ANNOT}" ]]; then
  echo "[missing] AOI annotation: ${ANNOT}"
  echo "          Build it from the GSE271689 series matrix / sample characteristics (or set \$ANNOT)."
  echo "          See methods/spatial/demo/README.md for the required columns."
  MISSING=1
fi

if [[ "${MISSING}" -ne 0 ]]; then
  echo
  echo "[stop] Prerequisites missing (see above). No analysis run; nothing fabricated."
  echo "       Supply the items above and re-run this script."
  exit 0
fi

mkdir -p "${OUT}"
echo "[run] 06 QC + Q3 normalization ..."
Rscript "${T}/06_geomx_qc_normalization.R" --dcc "${DCC}" --pkc "${PKC}" \
  --annotation "${ANNOT}" --out "${OUT}" || { echo "[fail] step 06"; exit 1; }

echo "[run] 07 compartment model (edit --compartment/--contrast for your sheet) ..."
Rscript "${T}/07_geomx_compartment_models.R" --rds "${OUT}/geomx_target_qnorm.rds" \
  --compartment segment --contrast Tumor,Immune --subject patient \
  --engine mixed --out "${OUT}" || echo "[warn] step 07 needs compartment/contrast matching your annotation"

echo "[done] outputs in ${OUT}"
echo "[next] deconvolution (08) needs a cell-profile matrix; survival (09) runs only if OS labels exist."
