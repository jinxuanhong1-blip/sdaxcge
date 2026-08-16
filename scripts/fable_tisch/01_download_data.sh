#!/usr/bin/env bash
# Download open lung cancer scRNA-seq objects from TISCH2 (http://tisch.comp-genomics.org).
# Both datasets are well under 2 GB and ship with uniform MAESTRO cell-type annotations,
# including an explicit "Malignant" major lineage, which is what we need for the
# tumor-vs-immune TACSTD2/CLDN4 comparison.
#
# NSCLC_GSE131907  Kim et al. 2020, Nat Commun  (LUAD, ~208k cells, tumor+normal+metastasis)
# NSCLC_EMTAB6149  Lambrechts et al. 2018, Nat Med  (NSCLC, ~52k cells)
#
# NOTE: we deliberately do NOT download the 12 GB extended LuCA atlas.

set -euo pipefail

DATA_DIR="${1:-data_fable_tisch}"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

BASE="https://tisch.compbio.cn/static/data"   # canonical redirect target of biostorage.aws.tisch.comp-genomics.org

for ds in NSCLC_GSE131907 NSCLC_EMTAB6149; do
  for f in expression.h5 CellMetainfo_table.tsv; do
    out="${ds}_${f}"
    if [[ ! -s "$out" ]]; then
      echo ">> downloading $out"
      curl -sL --retry 3 -o "$out" "${BASE}/${ds}/${ds}_${f}"
    fi
    ls -lh "$out"
  done
done
