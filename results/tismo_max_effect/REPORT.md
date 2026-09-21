# TISMO Tacstd2/Cldn4 vs ICI response and immune infiltrate

## Paper sentence

Across TISMO baseline tumors, the largest pan-cancer study-level association with at least five distinct values on each axis is **Cldn4 vs Granulocytes_mMCPcounter**, Spearman ρ = **0.869** (95% bootstrap CI 0.666 to 0.961; n = **23 studies**; leave-one-study-out ρ 0.850 to 0.906). The largest study-level odds ratio is **Tacstd2 vs Granulocytes_mMCPcounter**, OR = **100** (Woolf 95% CI 5.5 to 1830; Fisher p = 0.000346; high/high, high/low, low/high, low/low = 10/1/1/10; n = 22 studies). Cliff's δ for high vs low Tacstd2 on Granulocytes_mMCPcounter is **0.901** (n_high = 11, n_low = 11). All three point the same way: higher bulk Tacstd2/Cldn4 tracks higher granulocyte infiltrate. Restricting to studies with at least 8 baseline mice leaves ρ = 0.765 (n = 14) and the Tacstd2 granulocyte OR at 36 (6/1/1/6, n = 14).

小鼠单位上的最大 |ρ| 是 Cldn4 vs Neutrophils_mMCPcounter，ρ = 0.637（n = 299 只小鼠）。这个小鼠相关嵌在细胞系差异里。扣掉细胞系中位数之后，最大 |ρ| 是 0.336（Tacstd2 vs T CD4 Th2_xCell，ρ = -0.336，n = 294 只小鼠，17 个细胞系）。ICI 应答在全癌种 arm 单位上接近空：Tacstd2 OR = 0.866，Cldn4 OR = 0.682。

## Units

Two units were eligible for the maximum.

- **Mouse.** One unique SRX run. Baseline mice that TISMO copies onto two arms sharing an isotype control are kept once.
- **Study.** One source series (GEO `GSE` or ArrayExpress `ERP`/`RU` accession). The gene and the infiltrate score are the median of that series' baseline mice. A 3-mouse series counts the same as a 91-mouse series. The sensitivity below repeats the winners after dropping series with fewer than 8 baseline mice.

Responder labels are an arm property. In the Tacstd2 export, 62 of 64 treated arms are uniformly responders or uniformly non-responders. Only CT26 anti-PD-1 (GSE139475, 5 vs 4 mice) and YUMM1.7 BRAF-inhibitor + anti-PD-L2 (GSE103725, 2 vs 2 mice) carry mouse-level labels. A mouse-level response test would repeat the arm label. It is not used as a maximized effect.

Cell-line medians and within-line residuals are reported so the study-level correlation is not mistaken for a within-mouse effect. TISMO has no KL or KP lung line. LLC is the only lung carcinoma in this ICB table.

## Immune infiltrate, study unit

Search space: 2 genes × 84 immune scores at the study unit. 166 rows had enough non-missing studies to score, and 154 of those had at least five distinct values on each axis. Stromal, endothelial, fibroblast, and progenitor scores were excluded before ranking.

| Metric | Winner | Estimate | n | Direction |
|---|---|---|---|---|
| max \|ρ\| | Cldn4 vs Granulocytes_mMCPcounter | ρ = 0.869 (p = 7.54e-08) | 23 studies | higher gene, higher infiltrate |
| max \|Cliff δ\| | Tacstd2 vs Granulocytes_mMCPcounter | δ = 0.901 (CI 0.669 to 1.000) | 11 vs 11 studies | high-gene studies have higher infiltrate |
| max OR | Tacstd2 vs Granulocytes_mMCPcounter | OR = 100 (10/1/1/10) | 22 studies | high gene co-occurs with high infiltrate |

Odds ratios use one pre-specified cut: the median of each variable inside that spec. A split is ineligible when either side has fewer than 3 observations, when either side is under 20% of n, or when any cell of the 2×2 is 0. Zero-inflated scores such as mMCPcounter neutrophils often have a median of 0, so their OR is left undefined even when Spearman ρ is large. No outcome-guided cutpoint was scanned.

## Immune infiltrate, mouse unit

| Metric | Winner | Estimate | n |
|---|---|---|---|---|
| max \|ρ\| | Cldn4 vs Neutrophils_mMCPcounter | ρ = 0.637 (p = 1.74e-35) | 299 mice |
| max \|Cliff δ\| | Cldn4 vs Neutrophils_mMCPcounter | δ = 0.627 | 150 vs 149 mice |
| max OR | Cldn4 vs Dendritic resting_CIBERSORT_abs | OR = 9.24 (113/37/37/112; Fisher p = 8.24e-19) | 299 mice |

Mouse-level p-values treat SRX runs as independent. They are descriptive. After subtracting each cell line's median, the largest |ρ| is 0.336 (Tacstd2 vs T CD4 Th2_xCell, ρ = -0.336, n = 294 mice across 17 lines). The large infiltrate association is a between-study fact.

## CD8 reference

CD8 scores stayed in the search. At the study unit the strongest is Tacstd2 vs T CD8_CIBERSORT_abs, ρ = +0.738 (n = 22). It is positive and smaller than the granulocyte ρ. This table does not support CD8 exclusion.

| Gene | Score | Study ρ | n studies | Mouse ρ | n mice |
|---|---|---|---|---|---|
| Tacstd2 | T CD8_TIMER | 0.390 | 22 | -0.005 | 294 |
| Tacstd2 | T CD8_CIBERSORT_abs | 0.738 | 22 | 0.253 | 294 |
| Tacstd2 | T CD8_EPIC | 0.336 | 22 | 0.006 | 294 |
| Tacstd2 | T CD8_quanTIseq | 0.639 | 22 | 0.246 | 294 |
| Tacstd2 | T CD8_xCell | 0.120 | 22 | -0.006 | 294 |
| Tacstd2 | CD8 T_mMCPcounter | 0.666 | 22 | 0.315 | 287 |
| Cldn4 | T CD8_TIMER | 0.433 | 23 | -0.048 | 299 |
| Cldn4 | T CD8_CIBERSORT_abs | 0.651 | 23 | 0.196 | 299 |
| Cldn4 | T CD8_EPIC | 0.407 | 23 | 0.053 | 299 |
| Cldn4 | T CD8_quanTIseq | 0.653 | 23 | 0.152 | 299 |
| Cldn4 | T CD8_xCell | 0.173 | 23 | -0.092 | 299 |
| Cldn4 | CD8 T_mMCPcounter | 0.598 | 23 | 0.228 | 292 |

## ICI response

Pure arms only. An arm is one TISMO comparison group. Mixed arms are held out of this table and listed below.

| Gene | Unit | Stratum | n | ρ (gene vs response) | Cliff δ (responders vs non-responders) | OR for response given high gene | 2×2 R/NR among high, R/NR among low |
|---|---|---|---|---|---|---|---|
| Tacstd2 | arm | pan-cancer | 62 | 0.052 | 0.064 | 0.866 | 20/11/21/10 |
| Tacstd2 | arm | Mammary | 29 | 0.008 | 0.010 | 1.125 | 9/6/8/6 |
| Tacstd2 | arm | Melanoma | 13 | 0.275 | 0.333 | 1.250 | 5/2/4/2 |
| Tacstd2 | study | pan-cancer | 22 | -0.059 | NA | 3.062 | 7/4/4/7 |
| Cldn4 | arm | pan-cancer | 63 | -0.078 | -0.095 | 0.682 | 20/12/22/9 |
| Cldn4 | arm | Mammary | 29 | -0.192 | -0.225 | 0.635 | 8/7/9/5 |
| Cldn4 | arm | Melanoma | 13 | 0.581 | 0.667 | NA | 9/4/0/0 |
| Cldn4 | study | pan-cancer | 23 | -0.052 | NA | 3.500 | 8/4/4/7 |

The largest eligible response |Cliff's δ| is Cldn4 in Melanoma arms: δ = 0.667 (responders minus non-responders), 9 responder arms vs 4 non-responder arms, Spearman ρ = 0.581. Positive δ means responder arms have the higher gene median. That is the opposite of a high-gene resistance marker. Pan-cancer arm ORs sit near 1.

Mouse-level labels inside the only two mixed arms:

| Gene | Arm | Mice R vs NR | Median R | Median NR | Cliff δ |
|---|---|---|---|---|---|
| Tacstd2 | CT26_GSE139475_antiPD1 | 5 vs 4 | 1.143 | 1.395 | -0.200 |
| Tacstd2 | YUMM1.7_GSE103725_BRAF.inhibitor_antiPDL2 | 2 vs 2 | 0.194 | 0.047 | 0.500 |
| Cldn4 | CT26_GSE139475_antiPD1 | 5 vs 4 | 0.766 | 0.890 | -0.400 |
| Cldn4 | YUMM1.7_GSE103725_BRAF.inhibitor_antiPDL2 | 2 vs 2 | 0.034 | 0.073 | 0.000 |

## What was not done

- The locked Tacstd2 49/64 ICB-versus-control count was not recomputed and is not this result.
- LLC is not labeled KL. No KP lung line is in this table.
- Cutpoints were not tuned to maximize OR.
- Mouse-level response tests that repeat an arm label were not treated as independent mice.

## Methods

Expression: TISMO `POST /rtismo/gene/downVivoExprn` for Tacstd2 and Cldn4, all six ICB treatment classes and all vivo models. Infiltrate: `POST /rtismo/gene/downICBTreated` with `type=3`. That export returns every deconvolution score (TIMER, CIBERSORT absolute, EPIC, quanTIseq, xCell, mMCPcounter) rather than the one score named in the request. Scores were joined on SRX. Spearman ρ is the rank correlation. A ρ is eligible for the maximum only when both variables have at least five distinct values, so a tie pattern such as several melanoma series at zero does not rank as ρ = 1. Cliff's δ is `2U/(n1 n2) − 1` from the Mann–Whitney U of the high-gene group versus the low-gene group. OR is the median-split odds ratio with a Woolf interval and a two-sided Fisher exact test. Bootstrap intervals are percentile intervals from 2000 resamples (seed 0).

Code: `scripts/tismo_max_effect.py`. Tables: `results/tismo_max_effect/tables/`.

