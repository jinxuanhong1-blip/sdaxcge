# B1 analog — TCGA-KIRP TACSTD2 × CLDN4 surface-gene rank

Claim B1 (user PPT) says that in a **TCGA pan-cancer** surfaceome ranking,
**CLDN4** is the **top** co-expression partner of **TACSTD2 (TROP2)**.
This slice is the **KIRP-only analog**: same question, one cohort, honest rank.

It is the papillary renal-cell counterpart of `results/w200/B1_BRCA`.

## Methods (fixed before looking at the rank)

- **Expression:** UCSC Xena TCGA hub `TCGA.KIRP.sampleMap/HiSeqV2`,
  log2(norm_count+1), HGNC symbols. Same matrix family as the BRCA analog.
- **Samples:** primary tumours only (barcode sample-type `01`).
- **Surfaceome:** Bausch-Fluck et al., PNAS 2018, “in silico surfaceome only”
  sheet of table S3 (2,886 proteins). Universe = those genes present in the
  matrix, with TACSTD2 itself removed.
- **Primary metric:** Spearman correlation of each surface gene vs TACSTD2.
  Pearson is reported as a sensitivity check. FDR is Benjamini–Hochberg
  across the surfaceome universe.
- **w200:** reporting window for `top200.csv`. The full ranking is always
  written to `coexpression_TACSTD2_surfaceome.csv`.

## How to reproduce

```bash
python3 scripts/w200/B1_KIRP/download_data.py
python3 scripts/w200/B1_KIRP/run_analysis.py
```

Outputs land in `results/w200/B1_KIRP/`.

## Result

Filled after the ranking is computed. See `results/w200/B1_KIRP/summary.md`.
Do not invent a rank.

## What this does **not** show

- Not pan-cancer. A KIRP rank cannot confirm or refute the pan-cancer B1 claim.
- Bulk RNA-seq. Shared epithelial / stromal content can inflate two surface
  epithelial genes. No purity adjustment in this slice (that is a B1-PAAD
  extra, not part of the BRCA analog).
- Not protein, not IHC, not ICI outcome.
