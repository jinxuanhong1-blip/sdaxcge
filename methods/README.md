# Methods · GSE253013 extra LUAD scRNA figure

Public-only analysis of GEO [GSE253013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253013) (Sze/Xiang *Cancer Research* 2024, PMID 38335304).

This is **additive** paper content. User A3 on GSE207422 is taken as given and is **not** re-analyzed here.

## What this series actually is

- 89 GSM 10x lanes from **9 treatment-naïve** LUAD patients (tumor + adjacent non-tumor lung).
- Processed object: `GSE253013_all_luad_garnett_temp.rds.gz` (9.3 GB; double-gzipped XDR RDS).
- Series matrix is metadata only (6.7 KB; no expression).
- **Not** a neoadjuvant / ICI response cohort. GEO has no MPR / pCR / R / NR fields.

## Reproduce

```bash
python3 methods/gse253013_download.py --outdir data/gse253013
python3 methods/gse253013_extract.py \
  --rds data/gse253013/GSE253013_all_luad_garnett_temp.rds.gz \
  --outdir data/gse253013/extracted
python3 methods/gse253013_assemble.py --extracted data/gse253013/extracted
python3 methods/gse253013_analyze.py \
  --extracted data/gse253013/extracted \
  --outdir results/gse253013
```

The RDS does not fit in 16 GB RAM as a Seurat object. `gse253013_extract.py` walks the XDR stream and writes large sparse slots to disk, then keeps only the gene panel.

## Primary numbers (after the public RDS run)

Tumor, marker malignant-like, patient unit, n=9: TACSTD2 vs T/NK ρ=−0.72, p=0.030; CLDN4 ρ=−0.33, p=0.38. No MPR/R labels. See `results/gse253013/` and `paper/extra_gse253013_luad_scrna.md`.

## Scoring

- Lineage = argmax of mean log1p marker scores (epithelial / T / NK / B / myeloid / fibroblast / endothelial).
- Malignant-like = epithelial AND near-zero normal-lung markers (`SFTPA2`, `AGER`, `SCGB1A1`, `SCGB3A1`, `TPPP3`). This is a marker proxy, not CopyKAT.
- T/NK fraction and TACSTD2 / CLDN4 scores are aggregated **per patient** on tumor lanes (eligible: ≥10 malignant-like and ≥20 T/NK cells).
- Tests are Spearman ρ / p at the patient unit. Cell-level correlations are not reported.
- Response / MPR / R contrasts are omitted when labels are absent (this series).
