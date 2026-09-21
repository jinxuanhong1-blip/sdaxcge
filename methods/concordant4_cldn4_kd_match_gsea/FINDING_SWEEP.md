# SWEEP — tumor-intrinsic IFN on the concordant-4 malignant pseudobulk

**NOT a true KD.** Same malignant UMI-sums and the same unpurified IFN-γ NES as the quartile GSEA. No new cohort. No gene was added to the IFN sets after looking at NES. The STING core score reads CGAS from MB21D1 where that is the symbol, and STING1 from TMEM173. Those aliases are not forced into the prerank, so the IFN rank stays the published one. Reactome STING GSEA therefore misses CGAS and STING1 when those official symbols are absent from the shared matrix.

The unpurified IFN-γ leading edge contained T/B/NK markers (GZMA, CD69, CD86, LCP2). Those genes, and the rest of the locked lineage panel in `data/tnk_lineage_markers.tsv`, are removed in the purified arms. Myeloid genes that are not on that panel stay in the rank.

Primary purified arm: **Q1 vs Q4, full T/B/NK panel dropped** from the rank and from every gene set. Honest n = **18 lowest vs 16 highest**. Markers present in the matrix: 66. Spearman of that T/B/NK score vs CLDN4 %pos inside this contrast: -0.673.

## Call

| readout | unpurified Q4 | primary purified (panel dropped) |
|---|---|---|
| Hallmark IFN-γ NES | +3.789 (p 0.001) | +3.776 (p 0.001, thesis FDR 0.002) |
| direction vs private KD (IFN up after loss) | same direction | **same direction** |

Highest Hallmark IFN-γ NES among the purified arms: **+3.776** on `Q4 drop T/B/NK panel` (p 0.001). Lowest purified IFN-γ NES: **+1.245** on `continuous OLS + T/B/NK covariate`. Panel removal does not beat the unpurified NES. The sweep maximum is the panel-drop quartile, which is the least stringent purification. The stricter test is the T/B/NK score covariate, because that score tracks CLDN4 (Spearman -0.673 on the quartile tails). Quote the covariate row next to the maximum. Do not quote the maximum alone.

The purified leading edge still opens with CCL5 and CSF2RB. Those are not on the T/B/NK panel, so they were not removed.

Purified IFN-γ leading edge: `CCL5,CSF2RB,GBP4,TAP1,IL10RA,EPSTI1,IL15RA,STAT1,SAMHD1,FGL2,JAK2,SAMD9L,PSME2,PSMB9,PFKP,STAT4,OAS2,IFI35,LAP3,RTP4,TNFAIP3,VCAM1,STAT2,CASP1,XAF1`

## Sweep

Positive NES = higher when CLDN4 is lower. Thesis-aligned IFN and STING are positive. Thesis-aligned NHEJ is negative. Thesis FDR is BH inside IFN-γ, APM, Reactome STING, and Reactome NHEJ with histones removed.

| arm | n low/high | IFN-γ NES (p) | intrinsic ISG NES (p) | STING NES (p) | NHEJ no-histone NES (p) |
|---|---:|---:|---:|---:|---:|
| Q4 unpurified | 18/16 | +3.789 (0.001) | +3.291 (0.001) | +1.392 (0.037) | -0.865 (0.523) |
| Q4 drop leading-edge markers | 18/16 | +3.639 (0.001) | +3.335 (0.001) | +1.440 (0.024) | -0.860 (0.498) |
| Q4 drop T/B/NK panel | 18/16 | +3.776 (0.001) | +3.528 (0.001) | +1.468 (0.022) | -0.857 (0.540) |
| Q4 + T/B/NK covariate | 18/16 | +2.656 (0.001) | +3.113 (0.001) | -0.893 (0.432) | -0.947 (0.426) |
| Q5 drop T/B/NK panel | 15/14 | +3.447 (0.001) | +3.387 (0.001) | +1.097 (0.107) | -0.810 (0.576) |
| Q5 + T/B/NK covariate | 15/14 | +2.380 (0.001) | +3.071 (0.001) | -0.973 (0.334) | -0.801 (0.580) |
| continuous OLS, panel dropped | n=64 | +2.865 (0.001) | +2.878 (0.001) | +0.980 (0.166) | +0.732 (0.293) |
| continuous Spearman, panel dropped | n=64 | +2.986 (0.001) | +2.551 (0.001) | +1.043 (0.140) | -0.888 (0.424) |
| continuous OLS + T/B/NK covariate | n=64 | +1.245 (0.003) | +1.883 (0.001) | -0.990 (0.317) | +1.138 (0.078) |

### Primary purified sets

| set | ES | NES | nom p | thesis FDR | n |
|---|---:|---:|---:|---:|---:|
| IFN-γ | +0.669 | +3.776 | 0.001 | 0.002 | 191 |
| IFN-α | +0.681 | +3.485 | 0.001 | — | 95 |
| MHC-I / APM | +0.764 | +2.639 | 0.001 | 0.002 | 21 |
| intrinsic ISG | +0.753 | +3.528 | 0.001 | — | 59 |
| Reactome STING | +0.496 | +1.468 | 0.022 | 0.029 | 12 |
| Reactome NHEJ (no histone) | -0.264 | -0.857 | 0.540 | 0.540 | 30 |
| Reactome NHEJ (full) | -0.281 | -0.891 | 0.496 | — | 31 |
| KEGG TJ (CLDN4 out) | -0.252 | -1.040 | 0.370 | — | 147 |

Reactome NHEJ v2023.2.Hs is 68 genes and is mostly histone symbols. The no-histone set is the pre-specified GSEA module (histone symbols H2BC*, H3-*, H4C*, H2AX removed). The full set is the companion row. The 7-gene c-NHEJ core and the 5-gene STING core are module scores, not NES, because both are below the size floor of 8.

### Core module scores (low vs high)

Positive t = higher when CLDN4 is lower.

| arm | module | genes | low-vs-high t | p |
|---|---|---:|---:|---:|
| Q4 unpurified | cNHEJ_7 | 7 | -0.964 | 0.343 |
| Q4 unpurified | STING_5 | 5 | +2.205 | 0.036 |
| Q4 drop leading-edge markers | cNHEJ_7 | 7 | -0.964 | 0.343 |
| Q4 drop leading-edge markers | STING_5 | 5 | +2.205 | 0.036 |
| Q4 drop T/B/NK panel | cNHEJ_7 | 7 | -0.964 | 0.343 |
| Q4 drop T/B/NK panel | STING_5 | 5 | +2.205 | 0.036 |
| Q4 + T/B/NK covariate | cNHEJ_7 | 7 | -1.252 | 0.221 |
| Q4 + T/B/NK covariate | STING_5 | 5 | +0.764 | 0.451 |
| Q5 drop T/B/NK panel | cNHEJ_7 | 7 | -0.693 | 0.495 |
| Q5 drop T/B/NK panel | STING_5 | 5 | +1.931 | 0.065 |
| Q5 + T/B/NK covariate | cNHEJ_7 | 7 | -1.220 | 0.235 |
| Q5 + T/B/NK covariate | STING_5 | 5 | +0.553 | 0.586 |
| continuous OLS, panel dropped | cNHEJ_7 | 7 | -0.117 | 0.908 |
| continuous OLS, panel dropped | STING_5 | 5 | +1.111 | 0.271 |
| continuous OLS + T/B/NK covariate | cNHEJ_7 | 7 | -0.014 | 0.989 |
| continuous OLS + T/B/NK covariate | STING_5 | 5 | -0.633 | 0.529 |

## AUCell IFN, malignant cells only

This is not a NES. AUCell was already run on author-malignant cells. GSE148071 is not in this table. GSE123902 and GSE189357 have no cell-level AUCell in that file, so they are absent here rather than imputed. The AUCell IFN gene list is the intrinsic ISG list (no TCR, BCR, or NK markers).

| cohort | patients | cells | Spearman CLDN4 vs IFN AUCell | p | Q1−Q4 median IFN AUC | MWU p |
|---|---:|---:|---:|---:|---:|---:|
| GSE131907 (Q1 8 / Q4 8) | 31 | 31131 | -0.127 | 0.495 | +0.0013 | 0.574 |
| GSE205335 (Q1 6 / Q4 6) | 22 | 28512 | +0.098 | 0.665 | -0.0004 | 0.818 |
| GSE131907+GSE205335 within-dataset ranks | 53 | 59643 | +0.012 | 0.931 | — | — |

A negative Spearman would be the thesis direction (IFN activity higher when malignant CLDN4 is lower). The pooled rank correlation is the AUCell result for this sweep. It is not a NES and it is not a hit.

## What this sweep does not do

- It does not relabel the proxy as a knockdown.
- It does not drop GSE205335, or any other cohort, to raise the NES.
- It does not remove myeloid genes that were not on the T/B/NK panel.
- It does not bring GSE148071 into the AUCell pool.
- It does not treat the sweep maximum as the only number.

Reproduce: `python3 methods/concordant4_cldn4_kd_match_gsea/sweep_intrinsic.py`
