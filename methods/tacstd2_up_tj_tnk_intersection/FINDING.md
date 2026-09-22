# TACSTD2-high UP ∩ TJ ∩ low T/NK

Story gap after PR #741. Question: among genes up in TACSTD2-high
malignant cells on the #741 list, which are tight-junction members
**and** individually predict low T/NK?

Locked cohorts only: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Split and the UP list are TACSTD2, not CLDN4. This page does not
replace the locked CLDN4 %pos vs T/NK result (ρ = −0.531, n = 65).

## Call

| piece | rule |
|---|---|
| UP list | #741 ORA query: p<0.01 and logFC>0.25 on the TACSTD2 Q4 vs Q1 OLS. TACSTD2 held out. FDR arm of that DEG is empty, which is why #741 used this threshold. |
| TJ member | #741 TJ family: KEGG_TIGHT_JUNCTION ∪ GOBP_TIGHT_JUNCTION_ORGANIZATION ∪ CUSTOM_EPITHELIAL_ADHESION. TACSTD2 held out. |
| Predicts low T/NK | Within-cohort Spearman of malignant log2(TMM-CPM+1) vs frac_tnk, DerSimonian–Laird meta, **ρ < 0 and p < 0.05**. BH-FDR is reported inside the screened set and is not required for the nominal call. |
| Unit | Patient / donor / sample. Expression n = 64 (P4001 has T/NK and CLDN4 %pos but no count column). |

A sensitivity UP list (p<0.05 and logFC>0) is reported so CLDN4, which
misses the #741 threshold at p = 0.013, is not dropped in silence.
There is no per-cell matrix here, so this screen is pseudobulk expression.
The only %pos column in the locked tables is CLDN4.

## Set sizes

| set | n |
|---|---:|
| #741 UP list (p<0.01, logFC>0.25, TACSTD2 out) | 274 |
| nominal UP (p<0.05, logFC>0, TACSTD2 out) | 966 |
| TJ family defined (TACSTD2 out) | 224 |
| TJ family genes in the expression matrix | 206 |
| #741 UP ∩ TJ | 11 |
| #741 UP ∩ classical TJ (KEGG ∪ GO) | 8 |
| #741 UP ∩ TJ ∩ predicts low T/NK | 0 |
| nominal UP ∩ TJ | 27 |
| nominal UP ∩ TJ ∩ predicts low T/NK | 0 |
| TJ family with nominal low T/NK (ρ<0, p<0.05) | 9 |
| of those, also on the #741 UP list | 0 |

## Primary table: #741 UP ∩ TJ, each gene vs T/NK

Positive logFC = higher in TACSTD2 Q4. Positive ρ = higher malignant expression with higher T/NK.

| gene | logFC | DE p | TJ source | meta ρ | meta p | 95% CI | cohorts ρ<0 | I² | predicts low T/NK |
|---|---:|---:|---|---:|---:|---|---:|---:|---|
| GRHL2 | +1.084 | 0.004 | GOBP_TJ_ORG | -0.051 | 0.713 | -0.312 to +0.217 | 2/4 | 0.0% | no |
| LSR | +0.924 | 0.006 | GOBP_TJ_ORG | -0.015 | 0.913 | -0.279 to +0.251 | 3/4 | 0.0% | no |
| DSG2 | +1.828 | 0.003 | CUSTOM_ADHESION | -0.012 | 0.962 | -0.455 to +0.436 | 2/4 | 65.1% | no |
| ERBB2 | +1.053 | 0.003 | KEGG_TJ | +0.007 | 0.980 | -0.474 to +0.485 | 1/4 | 70.7% | no |
| JUP | +0.894 | 0.008 | CUSTOM_ADHESION | +0.015 | 0.914 | -0.251 to +0.279 | 1/4 | 0.0% | no |
| EPHA2 | +1.712 | 0.001 | GOBP_TJ_ORG | +0.038 | 0.784 | -0.230 to +0.300 | 2/4 | 0.0% | no |
| CLDN23 | +1.038 | 0.003 | KEGG_TJ|GOBP_TJ_ORG | +0.039 | 0.882 | -0.439 to +0.500 | 2/4 | 69.3% | no |
| AMOTL2 | +1.198 | 0.007 | KEGG_TJ | +0.159 | 0.320 | -0.154 to +0.442 | 1/4 | 22.4% | no |
| SYNPO | +1.703 | 9.98e-04 | KEGG_TJ | +0.170 | 0.215 | -0.099 to +0.417 | 1/4 | 0.0% | no |
| PLEC | +1.421 | 0.009 | GOBP_TJ_ORG | +0.172 | 0.499 | -0.319 to +0.591 | 2/4 | 68.6% | no |
| NECTIN4 | +1.513 | 0.002 | CUSTOM_ADHESION | +0.184 | 0.279 | -0.149 to +0.479 | 1/3 | 0.0% | no |

**Three-way intersection (#741 UP ∩ TJ ∩ predicts low T/NK) = 0 genes.** FDR < 0.05 inside this screen: 0.

No gene that is on the #741 UP list and in the TJ family individually predicts low T/NK. The 11 meta ρ values sit between -0.051 and +0.184; the smallest meta p is 0.215.

Restricting TJ to KEGG ∪ GO organization (dropping adhesion-only members) leaves 8 genes and 0 low-T/NK calls. The empty intersection does not depend on the custom adhesion panel.

The same call on the entire #741 UP list, TJ or not: 0 / 274 genes have meta ρ < 0 and p < 0.05. The most negative meta ρ inside the UP list is KRT13 (ρ = -0.359, p = 0.169). Significant UP-list associations that do exist are in the other direction (14 genes with ρ > 0 and p < 0.05; best FDR among those is 0.008).

NECTIN4 (and NECTIN2 on the nominal list) have no rank correlation in GSE131907 because expression does not vary there, so those metas use 3 cohorts.

## Sensitivity: nominal UP (p<0.05, logFC>0) ∩ TJ

n = 27. Genes with meta ρ < 0 and p < 0.05: **0**. FDR < 0.05 and ρ < 0: **0**.
CLDN4 is in this nominal list and is absent from the #741 UP list.

| gene | on #741 UP list | logFC | DE p | meta ρ | meta p | cohorts ρ<0 | predicts low T/NK |
|---|---|---:|---:|---:|---:|---:|---|
| CCND1 | no | +1.558 | 0.024 | -0.313 | 0.300 | 3/4 | no |
| OCLN | no | +0.989 | 0.013 | -0.240 | 0.188 | 3/4 | no |
| JUN | no | +0.742 | 0.019 | -0.209 | 0.369 | 3/4 | no |
| CLDN4 | no | +1.632 | 0.013 | -0.207 | 0.131 | 4/4 | no |
| CDH1 | no | +1.032 | 0.025 | -0.178 | 0.295 | 3/4 | no |
| SCRIB | no | +0.467 | 0.042 | -0.168 | 0.443 | 2/4 | no |
| PRKCZ | no | +0.731 | 0.010 | -0.089 | 0.591 | 3/4 | no |
| PRKCE | no | +0.744 | 0.015 | -0.076 | 0.735 | 3/4 | no |
| GRHL2 | yes | +1.084 | 0.004 | -0.051 | 0.713 | 2/4 | no |
| TJP1 | no | +0.593 | 0.037 | -0.018 | 0.922 | 2/4 | no |
| LSR | yes | +0.924 | 0.006 | -0.015 | 0.913 | 3/4 | no |
| DSG2 | yes | +1.828 | 0.003 | -0.012 | 0.962 | 2/4 | no |
| MSN | no | +1.147 | 0.022 | +0.003 | 0.990 | 1/4 | no |
| ERBB2 | yes | +1.053 | 0.003 | +0.007 | 0.980 | 1/4 | no |
| JUP | yes | +0.894 | 0.008 | +0.015 | 0.914 | 1/4 | no |
| EPHA2 | yes | +1.712 | 0.001 | +0.038 | 0.784 | 2/4 | no |
| CLDN23 | yes | +1.038 | 0.003 | +0.039 | 0.882 | 2/4 | no |
| MICALL2 | no | +0.892 | 0.019 | +0.086 | 0.533 | 3/4 | no |
| NECTIN2 | no | +0.620 | 0.028 | +0.119 | 0.484 | 2/3 | no |
| PARD6B | no | +0.592 | 0.042 | +0.141 | 0.499 | 2/4 | no |
| AMOTL2 | yes | +1.198 | 0.007 | +0.159 | 0.320 | 1/4 | no |
| SYNPO | yes | +1.703 | 9.98e-04 | +0.170 | 0.215 | 1/4 | no |
| PLEC | yes | +1.421 | 0.009 | +0.172 | 0.499 | 2/4 | no |
| NECTIN4 | yes | +1.513 | 0.002 | +0.184 | 0.279 | 1/3 | no |
| EZR | no | +0.554 | 0.048 | +0.221 | 0.226 | 1/4 | no |
| CLDN1 | no | +1.603 | 0.039 | +0.243 | 0.073 | 0/4 | no |
| CLDN16 | no | +0.863 | 0.047 | +0.308 | 0.029 | 1/4 | no |

CLDN16 is the one gene in this nominal table with meta p < 0.05 (ρ = +0.308, p = 0.029). The direction is higher T/NK, so it stays out of the low-T/NK intersection.

## Claudin rows, and the CLDN4 %pos reference

Every claudin in the #741 claudin panel that is present in the expression matrix is listed. CLDN4 %pos is a different measurement (cell-level percent positive), recomputed on these 64 units with the same meta. It is not a member of the UP list.

| gene | on #741 UP list | logFC | DE p | meta ρ | meta p | cohorts ρ<0 | readout |
|---|---|---:|---:|---:|---:|---:|---|
| CLDN4 | no | — | — | -0.523 | 2.86e-05 | 4/4 | %pos |
| CLDN19 | no | -0.018 | 0.703 | -0.328 | 0.014 | 4/4 | pseudobulk |
| CLDN5 | no | +0.357 | 0.592 | -0.288 | 0.033 | 3/4 | pseudobulk |
| CLDN3 | no | +0.652 | 0.336 | -0.283 | 0.122 | 2/4 | pseudobulk |
| CLDN15 | no | -0.090 | 0.787 | -0.227 | 0.095 | 3/4 | pseudobulk |
| CLDN7 | no | +0.713 | 0.089 | -0.219 | 0.387 | 2/4 | pseudobulk |
| CLDN4 | no | +1.632 | 0.013 | -0.207 | 0.131 | 4/4 | pseudobulk |
| CLDN22 | no | +0.102 | 0.532 | -0.189 | 0.303 | 3/4 | pseudobulk |
| CLDN11 | no | -0.513 | 0.267 | -0.187 | 0.173 | 4/4 | pseudobulk |
| CLDN14 | no | -0.152 | 0.626 | -0.166 | 0.227 | 3/4 | pseudobulk |
| CLDN12 | no | +0.166 | 0.573 | -0.139 | 0.408 | 3/4 | pseudobulk |
| CLDN6 | no | +1.055 | 0.151 | -0.129 | 0.478 | 2/4 | pseudobulk |
| CLDN9 | no | +0.830 | 0.189 | -0.030 | 0.827 | 3/4 | pseudobulk |
| CLDN2 | no | -0.308 | 0.392 | -0.009 | 0.951 | 2/4 | pseudobulk |
| CLDN23 | yes | +1.038 | 0.003 | +0.039 | 0.882 | 2/4 | pseudobulk |
| CLDN18 | no | -0.363 | 0.382 | +0.064 | 0.644 | 2/4 | pseudobulk |
| CLDN10 | no | +0.043 | 0.964 | +0.111 | 0.420 | 2/4 | pseudobulk |
| CLDN8 | no | +0.481 | 0.337 | +0.119 | 0.550 | 2/4 | pseudobulk |
| CLDN1 | no | +1.603 | 0.039 | +0.243 | 0.073 | 0/4 | pseudobulk |
| CLDN16 | no | +0.863 | 0.047 | +0.308 | 0.029 | 1/4 | pseudobulk |

## TJ genes with a nominal low-T/NK call are not the UP list

Across 206 TJ-family genes in the expression matrix, 9 have meta ρ < 0 and p < 0.05. Of those, 0 are on the #741 UP list and 0 are on the nominal UP list. None survives BH-FDR across the TJ family (smallest FDR = 0.090). They are recorded so the empty intersection is visible as disjoint sets, not as an untested corner. They are not an exclusion signature.

| gene | logFC | DE p | on #741 UP | meta ρ | meta p | FDR in TJ family | cohorts ρ<0 |
|---|---:|---:|---|---:|---:|---:|---:|
| PRKAA2 | -0.732 | 0.143 | no | -0.418 | 0.001 | 0.090 | 4/4 |
| DLG3 | +0.029 | 0.920 | no | -0.339 | 0.011 | 0.293 | 4/4 |
| PPP2R2C | -1.012 | 0.027 | no | -0.339 | 0.011 | 0.293 | 3/4 |
| CLDN19 | -0.018 | 0.703 | no | -0.328 | 0.014 | 0.293 | 4/4 |
| SLC39A9 | +0.120 | 0.591 | no | -0.323 | 0.016 | 0.293 | 4/4 |
| RAMP2 | -0.256 | 0.717 | no | -0.306 | 0.034 | 0.440 | 3/4 |
| CDC42 | +0.115 | 0.461 | no | -0.297 | 0.027 | 0.440 | 4/4 |
| CLDN5 | +0.357 | 0.592 | no | -0.288 | 0.033 | 0.440 | 3/4 |
| TGFBR1 | -0.263 | 0.414 | no | -0.270 | 0.046 | 0.552 | 3/4 |

## Lead-in

The TACSTD2-high malignant program and the low-T/NK direction do not meet on a TJ gene. The claudin that clears the #741 UP list is CLDN23 (logFC +1.038, DE p = 0.003), and its meta ρ vs T/NK is +0.039 (p = 0.882). CLDN4 misses that UP list (logFC +1.632, DE p = 0.013). On the same pseudobulk scale CLDN4 is negative in 4/4 cohorts (ρ = -0.207, p = 0.131, I² = 0.0%) and does not meet p < 0.05. The measurement that predicts low T/NK is CLDN4 %pos on these units (ρ = -0.523, p = 2.86e-05, negative in 4/4, I² = 0.0%). The screen this intersection points at is cell-level CLDN4 %pos.

CLDN1 is nominally up (DE p = 0.039) and its meta ρ vs T/NK is +0.243 (0/4 cohorts negative). A claudin that is up with TACSTD2 is not automatically a low-T/NK gene.

## Not claimed

- A TJ gene from the #741 UP list excludes T/NK. None does, at this rule.
- Pseudobulk CLDN4 expression is a significant T/NK predictor. It is directionally negative in 4/4 and not significant.
- The nine nominal TJ hits above are an exclusion list. They miss the UP list and miss FDR.
- Cell-level %pos for any gene other than CLDN4. Those columns are not in the locked tables.
- Private 8KL, KD coculture, or a Visium spatial claim.
- Pooling non-concordant accessions.

## Reproduce

```bash
python3 methods/tacstd2_up_tj_tnk_intersection/analyze.py
```

