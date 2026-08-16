#!/usr/bin/env bash
# Pointer/downloader for the E-MTAB-13530 Visium demo (human NSCLC lesions and
# non-involved lung; 36 sections; deconvolved with cell2location in the source
# study). ArrayExpress/BioStudies hosts the files.
# E-MTAB-13530 Visium 演示数据的下载指引（人 NSCLC 病灶与非受累肺；36 张切片；
# 原研究用 cell2location 去卷积）。文件托管于 ArrayExpress/BioStudies。
#
# SIZE WARNING: full raw Visium (all sections + images) can exceed the 2 GB
# processed-data budget. Prefer a single processed section (or a study-provided
# .h5ad) for a quick demo; the QC template accepts either a spaceranger 'outs'
# directory or an .h5ad. Inspect sizes on the study page before bulk download.
# 体积提醒：完整原始 Visium（所有切片+图像）可能超过 2 GB。快速演示建议只取单张
# 处理后的切片（或研究提供的 .h5ad）；QC 模板可接受 spaceranger 'outs' 目录或
# .h5ad。批量下载前请先在研究页面核对体积。
#
# Study landing pages (browse the File list to pick specific sections):
#   https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530
#   FTP:  https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files/
#
# Usage: bash methods/spatial/demo/download_visium_emtab13530.sh [FILE_RELATIVE_PATH]
# Pass an exact file path (relative to the study Files/ root) to fetch just that
# file; with no argument the script only lists guidance and exits.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${HERE}/data/EMTAB13530"
BASE="https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files"
mkdir -p "${DEST}"

if [[ $# -lt 1 ]]; then
  cat <<'EOF'
[info] No file specified. This dataset is large; choose specific files.
[info] 1) Open the study page and review the Files list + sizes:
          https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530
[info] 2) Re-run with a single file path relative to the Files/ root, e.g.:
          bash download_visium_emtab13530.sh <relative/path/to/section.h5ad>
[info] Keep the total under the 2 GB processed-data budget for the demo.
EOF
  exit 0
fi

REL="$1"
echo "[download] fetching ${REL} ..."
curl -fSL --create-dirs -o "${DEST}/$(basename "${REL}")" "${BASE}/${REL}"
echo "[download] saved to ${DEST}/$(basename "${REL}")"
echo "[next] run QC:  python methods/spatial/templates/01_visium_qc.py \\"
echo "           --input ${DEST}/<outs_dir_or.h5ad> --sample <id>"
