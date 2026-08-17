# Methods — GSE207422 PAGA / DPT scored by CLDN4 vs MPR

Additive. Does not rewrite `methods/scrna_paga/`, `methods/scrna_paga_cldn4/`,
`methods/gse148071_paga_cldn4/`, or the given GSE207422 T/NK / dual-high slices.

**Question.** Among malignant-like epithelial cells on public GSE207422
(neoadjuvant PD-1 + platinum), does **CLDN4** sit on a leftover/AT2-rooted
diffusion trajectory that also separates MPR vs NMPR?

**Not a TACSTD2 redo.** TACSTD2 is a companion Spearman only. No dual-high gate.
GSE207422-only T/NK was flat in the given slice and is not re-audited.

| Item | Choice |
| --- | --- |
| Cohort | GSE207422 (Hu et al. 2023); 15 scRNA samples, tests on 12 post-treatment patients |
| Counts | Public GEO processed UMI only. Not HRA001033. |
| Labels | Marker lineage (Hu canonical argmax). Author CopyKAT barcodes are not public. |
| A3-malignant-like | Epithelial AND zero UMI for SFTPA2 / AGER / SCGB1A1 / SCGB3A1 / TPPP3 |
| Graph | log1p(CP10k) → Seurat HVG 3000 → PCA → k-NN (30/30) on epithelium |
| Clusters | Leiden 0.6; PAGA on Leiden |
| Cap | max 500 cells/sample if n>8000 (seed 0) |
| Root | Leftover epithelium, not CLDN4-high, median AT2 score |
| Score | CLDN4 continuous + tertiles |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| Inferential n | Patient. Cell-level ρ is descriptive. |
| Primary tests | Patient-mean CLDN4 vs DPT (A3-malignant, ≥10 cells); DPT vs MPR (exact MWU) |
| Extra figure | Within-sample CLDN4-high vs low program scores (min 8 cells/arm) |
| Unused | Slingshot (R absent). Dual-high. T/NK re-audit. |

**Honest n.** This set is small (12 post patients; MPR n=4 including pCR). Several
MPR residuals have 0 A3-malignant-like cells. Do not write the graph cell count
as the test n. `SFTPC` is absent from the public UMI; the AT2 score uses the
remaining locked AT2 genes.
