#!/usr/bin/env bash
# Fetch all inputs for the hunt_pacific analyses into data/hunt_pacific/.
# All URLs verified 2026-08-16.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/hunt_pacific
cd data/hunt_pacific

# Open durvalumab NSCLC trial FPKM matrices (GEO supplementary files)
curl -sO "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"
curl -sO "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/suppl/GSE248378_Durva_Post_FPKMs.txt.gz"

# ESTIMATE stromal/immune gene sets (official R package, R-Forge)
curl -sL -o estimate.tar.gz "http://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz"
tar xzf estimate.tar.gz estimate/inst/extdata/SI_geneset.gmt -O > SI_geneset.gmt

# TCGA PanCanAtlas ABSOLUTE purity (official GDC publication-page file)
curl -s -o TCGA_mastercalls.abs_tables_JSedit.fixed.txt \
  "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"

ls -la
