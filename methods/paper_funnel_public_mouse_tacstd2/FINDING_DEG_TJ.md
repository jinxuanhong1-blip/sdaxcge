# Tacstd2-high epithelial DEG → tight junction (public integrate mice)

Public processed counts only: GSE154977, GSE180963, GSE165641. Private 8 KL matrices were not read and were not merged. Epithelium is the locked gate from the integrate (Epcam∧structural∧Ptprc−, or AT2-FACS Ptprc− with Epcam∨structural). **Tacstd2 is not an epithelial caller.**

Within each mouse, epithelial cells with Tacstd2 count > 0 are compared to epithelial cells with Tacstd2 count = 0. Mice below the high/low cell floor are inventoried and dropped from DEG.

## Module scores (primary)

Frozen TJ lists from the GSE137244 transfer signatures (`gene_sets/signatures.gmt`). Cell-level mean log1p module score; one-sided Mann–Whitney that Tacstd2-high > Tacstd2-low; binomial test on the sign of the mouse-level delta.

| set | n mice | Δ>0 | Δ>0 and p<0.05 | mean Δ | binomial p (direction) |
|---|---:|---:|---:|---:|---:|
| CLDN4_TJ_EDGE | 7 | 7 | 7 | +0.1652 | 0.007812 |
| TJ_EPITHELIAL | 7 | 7 | 7 | +0.1439 | 0.007812 |
| TJ_TISMO | 7 | 7 | 7 | +0.2608 | 0.007812 |

**Headline (TJ_TISMO 7-gene):** 7/7 mice have Tacstd2-high epithelium above Tacstd2-low (mean Δ=+0.2608; one-sided binomial p=0.007812). 7/7 also have within-mouse MW p<0.05.

## Ranked gene enrichment (secondary)

Mean across-mouse pseudobulk delta (Tacstd2-high − low) for genes present in ≥2 mice. TJ_TISMO rank Mann–Whitney (set deltas greater than background) p=3.601e-09. Strongest hypergeometric among top-N cutoffs: top 500, k=5/7, enrichment=28.58, p=1.926e-07.

Cell-level p-values inside a mouse are descriptive for calling the module direction; the mouse is the inferential unit for the binomial and for the across-mouse gene ranking.

Mice scored: GSE154977_KP_30w_Cis72_m5, GSE154977_KP_30w_Cis72_m6, GSE154977_KP_30w_ND_m3, GSE154977_KP_30w_ND_m4, GSE180963_KL, GSE165641_KL1, GSE165641_KL2.
Private 8 KL mice used: 0.
