#!/usr/bin/env bash
# End-to-end demo: resources -> GEO ICI + TCGA NSCLC -> scores -> TACSTD2/CLDN4 correlations.
# Run from anywhere; paths are resolved relative to methods/bulk_immune/.
set -euo pipefail

HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}"

python3 scripts/00_fetch_resources.py
python3 scripts/01_prepare_geo.py
python3 scripts/02_prepare_tcga.py

python3 scripts/03_score_cohort.py \
  --name GSE126044 \
  --log data/GSE126044.log2cpm.tsv.gz \
  --linear data/GSE126044.counts.tsv.gz \
  --linear-is-counts \
  --pheno data/GSE126044.pheno.tsv

python3 scripts/03_score_cohort.py \
  --name GSE135222 \
  --log data/GSE135222.log2fpkm.tsv.gz \
  --linear data/GSE135222.fpkm.tsv.gz \
  --pheno data/GSE135222.pheno.tsv

python3 scripts/03_score_cohort.py \
  --name TCGA_NSCLC \
  --log data/TCGA_NSCLC.log2rsem.tsv.gz \
  --linear data/TCGA_NSCLC.linear.tsv.gz \
  --pheno data/TCGA_NSCLC.pheno.tsv

python3 scripts/04_correlate.py \
  --indir results/GSE126044 \
  --group-col response --positive responder \
  --batch-col tissue_preservation \
  --purity-col estimate:ESTIMATEScore

python3 scripts/04_correlate.py \
  --indir results/GSE135222 \
  --time-col pfs_time --event-col pfs_event \
  --purity-col estimate:ESTIMATEScore

python3 scripts/04_correlate.py \
  --indir results/TCGA_NSCLC \
  --batch-col cohort \
  --strata-col cohort \
  --purity-col estimate:ESTIMATEScore

python3 scripts/06_make_figures.py --root "$HERE"
python3 scripts/07_concordance.py --root "$HERE"

echo "demo complete. see $HERE/results/"
