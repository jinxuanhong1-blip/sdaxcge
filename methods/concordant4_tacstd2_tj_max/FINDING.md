# Concordant-4: TACSTD2-high vs low, max-effect DEG / ORA / GSEA

## Best honest #1–3

A junction term leads ORA when the hit contains occludin, ZO-1, CLDN4, or E-cadherin. GSEA does not put any junction term at rank 1.

Junction-led ORA whose rank-1 overlap contains a tight-junction structural gene (CDH1, OCLN, TJP1): **tercile_ols_t||p<0.01 & logFC>0.25**.

Top 3: 1. GOBP_CELL_CELL_JUNCTION_ORGANIZATION (enrichment 7.28, p 4.59e-13, FDR 3.12e-11); 2. HALLMARK_APICAL_JUNCTION (enrichment 6.82, p 1.79e-11, FDR 6.10e-10); 3. GOBP_CELL_JUNCTION_ORGANIZATION (enrichment 3.19, p 1.92e-10, FDR 4.36e-09).

Rank-1 overlap: CDH1,CDH3,DLG5,DSG2,EPHA2,EPHB2,GJB2,GJB6,GRHL1,GRHL2,HOPX,INAVA,JUP,LSR,OCLN,PARD6B,PKP3,PLEC,PLEKHA7,RHOC,TJP1,TRPV4.

The pre-specified score (lowest junction rank, then more story terms in the top 3, then higher enrichment) picks a different ORA, **q4q1_moderated_t||top100 positive stat**. Top 3 there: 1. HALLMARK_APICAL_JUNCTION (enrichment 11.17, FDR 8.16e-06); 2. GOBP_CELL_CELL_ADHESION (enrichment 3.87, FDR 2.36e-04); 3. GOBP_KERATINIZATION (enrichment 21.89, FDR 7.58e-04). Rank-1 overlap: ADAM9,CDH3,ICAM1,ICAM4,LAMC2,MPZL2,NECTIN4,PTK2,TGFBI. That list has no CLDN, OCLN, or TJP gene, so it is not a claudin result.

GSEA top 3 under the spec with the best junction rank (q4q1_stouffer_z__w1): 1. HALLMARK_TNFA_SIGNALING_VIA_NFKB (NES +3.185, FDR 0.001); 2. HALLMARK_INFLAMMATORY_RESPONSE (NES +2.781, FDR 0.001); 3. GOBP_KERATINOCYTE_DIFFERENTIATION (NES +2.626, FDR 0.001). Strict TJ in that spec: GOBP_TIGHT_JUNCTION_ORGANIZATION at rank 20 (NES +2.131). Best adhesion in that spec: CUSTOM_EPITHELIAL_ADHESION at rank 5 (NES +2.556).

Query size for the structural-gene ORA: 364 genes. Best strict tight-junction set on this same table: GOBP_TIGHT_JUNCTION_ORGANIZATION at rank 13, enrichment 7.15, FDR 8.34e-05, overlap CLDN4,EPHA2,EPHB2,GRHL2,LSR,OCLN,PLEC,TJP1.

Locked malignant cohorts only: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Split is TACSTD2 expression in malignant pseudobulk, not CLDN4 % positive.
TACSTD2 is removed from the ranked list before ORA and GSEA.
The grid is fixed in `analyze_max.py`. The winner is the eligible spec with the lowest junction rank. A strict tight-junction term is preferred when it is rank 1. Every spec remains in the sweep table.

## Honest n

Locked unit table: **65**. Expression matrix after the count join: **64**. Reference Q4 vs Q1: **15 vs 19**.
Do not quote 65 as the differential-expression n.

## GSEA winner

**q4q1_stouffer_z__w1**. Equal-weight Stouffer combination of within-cohort Welch z (Q4 vs Q1, ≥2 per side). GSEA weight p=1.

Verdict: **no TJ/junction term is rank 1; best is HALLMARK_APICAL_JUNCTION at rank 4**.

Sets tested: 67. Permutations: 1000. Weight: 1.0. Seed: 42.

| rank | term | class | NES | nominal p | FDR | genes in rank |
|---:|---|---|---:|---:|---:|---:|
| 1 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | other | +3.185 | 9.99e-04 | 0.001 | 199 |
| 2 | HALLMARK_INFLAMMATORY_RESPONSE | other | +2.781 | 9.99e-04 | 0.001 | 194 |
| 3 | GOBP_KERATINOCYTE_DIFFERENTIATION | keratin | +2.626 | 9.99e-04 | 0.001 | 39 |

Strict TJ under this spec: **GOBP_TIGHT_JUNCTION_ORGANIZATION**, rank 20, NES +2.131.
Best keratin: **GOBP_KERATINOCYTE_DIFFERENTIATION**, rank 3, NES +2.626.
Best adhesion: **CUSTOM_EPITHELIAL_ADHESION**, rank 5, NES +2.556.

Highest-ranked member genes of the rank-1 term: TNC,ICAM1,ZFP36,TRIP10,HBEGF,JUN,G0S2,CCL2,IER3,DUSP1,OLR1,SDC4.

## ORA winner

**q4q1_moderated_t||top100 positive stat**. Query: top100 positive stat. n query genes = 100 (TACSTD2 held out).

Verdict: **HALLMARK_APICAL_JUNCTION is rank 1; best strict TJ is KEGG_TIGHT_JUNCTION at rank 25**.

| rank | term | class | overlap | enrichment | p | FDR |
|---:|---|---|---:|---:|---:|---:|
| 1 | HALLMARK_APICAL_JUNCTION | junction | 9 | 11.17 | 1.20e-07 | 8.16e-06 |
| 2 | GOBP_CELL_CELL_ADHESION | adhesion | 15 | 3.87 | 6.95e-06 | 2.36e-04 |
| 3 | GOBP_KERATINIZATION | keratin | 4 | 21.89 | 3.34e-05 | 7.58e-04 |

Strict TJ under this ORA: **KEGG_TIGHT_JUNCTION**, rank 25, enrichment 3.13.
Best keratin: **GOBP_KERATINIZATION**, rank 3. Best adhesion: **GOBP_CELL_CELL_ADHESION**, rank 2.

## Reference spec (not maximized)

Cohort-adjusted OLS *t*, within-cohort TACSTD2 Q4 vs Q1, weighted GSEA (p = 1). This is the previous funnel's ranking statistic, now scored on the full set list including the added MSigDB junction and adhesion sets.

| rank | term | NES | FDR |
|---:|---|---:|---:|
| 1 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | +3.200 | 0.002 |
| 2 | GOBP_KERATINIZATION | +2.941 | 0.002 |
| 3 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | +2.833 | 0.002 |

Reference strict-TJ rank: 24 (GOBP_TIGHT_JUNCTION_ORGANIZATION, NES +2.096). Reference best junction rank: 9 (HALLMARK_APICAL_JUNCTION).

## DEG on the GSEA-winning contrast

Positive means higher with TACSTD2-high. For the Stouffer winner the effect is the mean of within-cohort Q4−Q1 deltas of the set's mean expression, and the p-value is a one-sample t-test across four cohorts. That test has almost no power. It is not the GSEA p-value.

| set | class | genes | effect | test stat | p | mean gene logFC |
|---|---|---:|---:|---:|---:|---:|
| REACTOME_TIGHT_JUNCTION_INTERACTIONS | strict_tj | 27 | +0.322 | +2.35 | 0.100 | +0.322 |
| KEGG_TIGHT_JUNCTION | strict_tj | 154 | +0.172 | +2.39 | 0.097 | +0.172 |
| GOBP_TIGHT_JUNCTION_ORGANIZATION | strict_tj | 74 | +0.267 | +2.10 | 0.126 | +0.267 |
| HALLMARK_APICAL_JUNCTION | junction | 194 | +0.310 | +2.06 | 0.132 | +0.310 |
| GOBP_APICAL_JUNCTION_ASSEMBLY | junction | 62 | +0.278 | +1.93 | 0.149 | +0.278 |
| GOBP_KERATINIZATION | keratin | 44 | +0.781 | +2.47 | 0.090 | +0.781 |
| KRT_EPITHELIAL | keratin | 17 | +1.220 | +2.95 | 0.060 | +1.220 |
| CUSTOM_EPITHELIAL_ADHESION | adhesion | 23 | +0.661 | +3.88 | 0.030 | +0.661 |
| GOBP_CELL_CELL_ADHESION | adhesion | 933 | +0.186 | +1.55 | 0.218 | +0.186 |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | other | 199 | +0.497 | +2.37 | 0.098 | +0.497 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | other | 199 | +0.442 | +1.15 | 0.334 | +0.442 |
| HALLMARK_INFLAMMATORY_RESPONSE | other | 194 | +0.365 | +1.93 | 0.148 | +0.365 |

Focal genes:

| gene | logFC | rank stat | p |
|---|---:|---:|---:|
| TACSTD2 | +3.556 | +5.407 | 6.42e-08 |
| CLDN1 | +1.353 | +1.835 | 0.067 |
| CLDN3 | +0.685 | +1.343 | 0.179 |
| CLDN4 | +1.545 | +3.088 | 0.002 |
| CLDN7 | +0.707 | +2.145 | 0.032 |
| OCLN | +0.992 | +2.400 | 0.016 |
| F11R | +0.376 | +1.647 | 0.100 |
| TJP1 | +0.460 | +1.949 | 0.051 |
| CDH1 | +0.917 | +2.458 | 0.014 |
| EPCAM | +0.549 | +1.644 | 0.100 |
| KRT5 | +1.837 | +1.123 | 0.262 |
| KRT6A | +1.563 | +0.857 | 0.391 |
| KRT8 | +1.003 | +2.046 | 0.041 |
| KRT17 | +2.464 | +1.905 | 0.057 |
| KRT18 | +0.806 | +1.376 | 0.169 |
| KRT19 | +1.955 | +3.584 | 3.38e-04 |

## What was searched

GSEA specs run: 36. Eligible: 35. ORA specs run: 126. Eligible: 125.
Contrasts: Q4 vs Q1, within-cohort median, tercile, quintile, and continuous TACSTD2 (z within cohort).
Ranks: OLS t, OLS logFC, moderated t, signed −log10 p, Stouffer z, mean / minimum / concordance-weighted within-cohort logFC, Fisher-z meta Spearman, cohort-residual Spearman.
GSEA weights: classic weighted (p = 1) and unweighted (p = 0). Gene-set permutation, not sample permutation.
Universe: the previous Hallmark / KEGG / GO collection plus MSigDB GOBP cell-cell adhesion, cell-junction organization, apical junction assembly, and Reactome tight junction interactions. No set was removed after the sweep.
Ineligible: TACSTD2 not higher in the high group, fewer than 3 cohorts, fewer than 8 samples on a side, junction nominal p ≥ 0.05, or an ORA query outside 25–800 genes.

## Not claimed

This page does not re-fit CLDN4 % positive vs T/NK. It does not use private 8KL matrices. It does not add GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526. A rank of 1 is a rank inside this universe, not inside every pathway catalog.

## Reproduce

```bash
python3 methods/concordant4_tacstd2_tj_max/analyze_max.py
```

