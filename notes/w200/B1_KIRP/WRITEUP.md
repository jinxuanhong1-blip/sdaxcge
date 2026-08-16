# B1 analog — TCGA-KIRP TACSTD2 × CLDN4 surface-gene rank

Claim B1 (user PPT) says that in a **TCGA pan-cancer** surfaceome ranking,
**CLDN4** is the **top** co-expression partner of **TACSTD2 (TROP2)**.
This slice is the **KIRP-only analog**: same question, one cohort, honest rank.

It is the papillary renal-cell counterpart of `results/w200/B1_BRCA`.

## Methods (fixed before looking at the rank)

- **Expression:** UCSC Xena TCGA hub `TCGA.KIRP.sampleMap/HiSeqV2`,
  log2(norm_count+1), HGNC symbols. Same matrix family as the BRCA analog.
- **Samples:** 290 primary tumours (barcode sample-type `01`). Matrix has
  323 columns / 20,530 genes before sample-type filter.
- **Surfaceome:** Bausch-Fluck et al., PNAS 2018, “in silico surfaceome only”
  sheet of table S3. 2,618 surface genes present after removing TACSTD2;
  10 zero-variance genes dropped (undefined correlation); **2,608 ranked**.
- **Primary metric:** Spearman correlation of each surface gene vs TACSTD2.
  Pearson is a sensitivity check. FDR is Benjamini–Hochberg across the
  ranked surfaceome.
- **w200:** reporting window for `top200.csv`. Full ranking:
  `coexpression_TACSTD2_surfaceome.csv`.

## How to reproduce

```bash
python3 scripts/w200/B1_KIRP/download_data.py
python3 scripts/w200/B1_KIRP/run_analysis.py
```

Outputs land in `results/w200/B1_KIRP/`.

## Result (computed, not invented)

**No. CLDN4 is not the top TACSTD2 surface partner in TCGA-KIRP.**

| quantity | value |
| --- | --- |
| CLDN4 Spearman rank | **#57 of 2,608** (97.9th percentile) |
| CLDN4 Spearman rho | 0.435 (FDR q = 2.68e-13) |
| CLDN4 Pearson rank | #80 of 2,608 (rho = 0.386) |
| Actual #1 | **MUC1** (rho = 0.651) |
| Actual #2 | **CLDN7** (rho = 0.621) |

CLDN4 is co-expressed with TACSTD2 (rho is real and FDR-significant) but
it is **far from rank 1**. Another claudin, CLDN7, outranks it by a wide
margin. Compare the BRCA analog, where CLDN4 was #4 of 2,618 (rho = 0.348).
KIRP does **not** reproduce a “CLDN4 = #1 surface partner” readout.

## What this does **not** show

- Not pan-cancer. A KIRP rank cannot confirm or refute the pan-cancer B1 claim.
- Bulk RNA-seq. Shared epithelial content can inflate two surface epithelial
  genes. No purity adjustment here (that extra is in the PAAD analog, not
  in the BRCA B1 protocol this slice follows).
- Not protein, not IHC, not ICI outcome.
