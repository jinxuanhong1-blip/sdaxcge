# Methods — winning-pair real Slingshot, CLDN4 only

ADDITIVE. Does **not** rewrite `methods/scrna_paga_cldn4/` (PR #325) or
`methods/winpair_131907_205335_slingshot_cldn4/` (PR #449 DPT fallback).
GSE207422 and GSE148071 are not added. TACSTD2 does not define groups.
PR #320 T/NK Q4 vs Q1 r=−0.705 is given and is not re-audited.

**Question.** On the winning public pair GSE131907+GSE205335 **malignant/epithelial
cells**, where does **CLDN4** sit on a Slingshot lineage relative to a
CLDN4-excluded barrier/keratin program and an IFN program? Root is AT2-like or
lowest-CLDN4, never CLDN4-high.

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix | Kept | Dropped |
| --- | --- | --- | --- | --- |
| GSE131907 | Kim et al., *Nat Commun* 2020, PMID 32385277 | raw UMI text + author annotation | `Cell_type == Epithelial cells` and `Sample_Origin` ∈ {tLung, nLung, tL/B, mLN, mBrain} | PE unlabeled epithelium; 2.86 GB log2TPM; EGA FASTQ |
| GSE205335 | Ahn / Lee, *eLife* 2024 | UMI dgCMatrix + author identity + SOFT | `lineage.total == Epithelial cells` on non-normal tissues | Normal Lung / LN / Brain; EGA FASTQ |

## Trajectory clock (real Slingshot, not DPT)

Slingshot (Street et al. 2018) is the primary clock.

1. Try Bioconductor `slingshot` via `Rscript` if present.
2. If R is missing, use **pyslingshot-bio** (`pyslingshot.core.slingshot`), a
   documented Python port of Street et al. 2018
   (https://github.com/omicverse/py-Slingshot; reported Spearman ≥0.99 vs R on a
   Y-shaped fixture). Package `__init__` plotting (`ggplot2_py`) is not required;
   the core MST + principal-curve fit is imported directly.
3. If that import fails, a local Street-2018 port in `scripts/slingshot_py.py`
   (MST over centroids → root-to-leaf lineages → LOWESS principal curve) is used.

Diffusion pseudotime is **not** the primary clock. Palantir is an optional companion.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA 30 → Harmony (`batch = dataset`) → k-NN 30 |
| Clusters | Leiden 0.4; PAGA on Leiden (geometry) |
| Root | Leiden with most GSE131907 nLung author AT2, else highest AT2 among below-median-CLDN4 clusters. **Max-CLDN4 cluster is excluded.** |
| Score | CLDN4 continuous |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| IFN | STAT1, IRF1, IRF7, ISG15, IFIT1/2/3, OAS1/2, MX1, IFI27, IFI44, IFI44L, RSAD2, IFI6, GBP1, IDO1, CXCL9/10/11 (**no CLDN4**) |
| Inferential n | GSE131907 `Sample` + GSE205335 `patient`. Cell-level ρ is descriptive. |
| Cap | ≤280 cells / unit after protecting nLung AT2 |
| Unused | GSE207422; GSE148071; dual-high TACSTD2∩CLDN4; ICI / MPR / RECIST as a trajectory label |

## Tests

Primary: sample-level Spearman of mean CLDN4 vs mean Slingshot PT / AT2 / barrier / IFN.
BH inside that list only. Along-curve binned means use the primary lineage
(terminal = highest mean CLDN4 among root-started lineages).
