#!/usr/bin/env bash
# Download open lung cancer scRNA-seq objects from TISCH2 (http://tisch.comp-genomics.org).
# Both datasets are well under 2 GB and ship with uniform MAESTRO cell-type annotations,
# including an explicit "Malignant" major lineage, which is what we need for the
# tumor-vs-immune TACSTD2/CLDN4 comparison.
#
# NSCLC_GSE131907  Kim et al. 2020, Nat Commun  (LUAD, ~208k cells; subset to tumor tissue)
# NSCLC_EMTAB6149  Lambrechts et al. 2018, Nat Med  (NSCLC, ~40k cells, explicit Malignant)
# NSCLC_GSE127465  Zilionis et al. 2019, Immunity  (~31k cells, explicit Malignant)
# NSCLC_GSE148071  Wu et al. 2021, Nat Commun  (~82k cells, explicit Malignant + leftover epithelium)
#
# NOTE: we deliberately do NOT download the 12 GB extended LuCA atlas.
# GSE131907 expression.h5 is ~2.5 GB on disk but is subsetted to tumor-tissue
# cells (~65k) before any statistics are computed.

set -euo pipefail

DATA_DIR="${1:-data_fable_tisch}"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

BASE="https://tisch.compbio.cn/static/data"   # canonical redirect target of biostorage.aws.tisch.comp-genomics.org

for ds in NSCLC_GSE131907 NSCLC_EMTAB6149 NSCLC_GSE127465 NSCLC_GSE148071; do
  for f in expression.h5 CellMetainfo_table.tsv; do
    out="${ds}_${f}"
    if [[ ! -s "$out" ]]; then
      echo ">> downloading $out"
      curl -sL --retry 3 -o "$out" "${BASE}/${ds}/${ds}_${f}"
    fi
    ls -lh "$out"
  done
done
