# scripts

```bash
# GSE243013 (author TCR table + residual scores)
python3 run_gse243013.py \
  --tcr data/public/gse243013/GSE243013_T_with_TCR_annotation.csv.gz \
  --meta data/public/gse243013/GSE243013_NSCLC_immune_scRNA_metadata.csv.gz \
  --residual methods/scrna_tcr/data/gse243013_residual_tacstd2_cldn4.tsv \
  --out methods/scrna_tcr/results/GSE243013

# Caushi public processed (VDJ + RDS + GSE176022). No EGA.
python3 run_caushi.py \
  --vdj-dir data/public/caushi/vdj \
  --series data/public/caushi/GSE176021_series_matrix.txt.gz \
  --cd3-rds data/public/caushi/GSE176021_CD3_annotations.rds \
  --gse176022-tar data/public/caushi/GSE176022_RAW.tar \
  --out methods/scrna_tcr/results/Caushi

# Leftover GSE179994
python3 run_leftover.py \
  --tcr data/public/leftover/GSE179994_all.scTCR.tsv.gz \
  --meta data/public/leftover/GSE179994_Tcell.metadata.tsv.gz \
  --series data/public/leftover/GSE179994_series_matrix.txt.gz \
  --out methods/scrna_tcr/results/leftover
```

Unit is always the sample/patient. Cells are never the MPR replicate.
