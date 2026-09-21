# TISMO: Tacstd2 versus infiltrate and ICI, partialling Cldn4

## Paper sentence

Tacstd2 and Cldn4 are one bulk axis in the 22 paired TISMO series (study Spearman ρ = 0.863; mouse ρ = 0.830; 294 baseline mice). On that axis, the study-level Tacstd2–granulocyte link shrinks and the largest Tacstd2–CD8 link does not. ICI response has nothing to attenuate.

Study unit, mMCPcounter granulocytes: Tacstd2 ρ = 0.860 falls to a partial ρ of 0.473 given Cldn4 (resample CI -0.219 to 0.818). Shrinkage Δ = 0.387 (CI 0.095 to 1.001). The shrinkage interval stays above zero. The interval for the partial correlation includes zero. Cldn4 partialled on Tacstd2 moves from 0.852 to 0.426 (CI 0.033 to 0.779). The two leftover partials are 0.473 and 0.426. Tacstd2’s interval includes zero. Cldn4’s interval stays off zero. That is not a clean assignment of the granulocyte link to either gene.

Study unit, CIBERSORT CD8, the largest pre-specified CD8 score: Tacstd2 ρ = 0.738, partial given Cldn4 ρ = 0.556 (CI 0.044 to 0.800), shrinkage Δ = 0.183 (CI -0.061 to 0.587). The point estimate is closer to zero. The interval for the shrinkage includes zero. The partial stays positive. This is not CD8 exclusion.

Mouse unit, the two myeloid scores diverge. Tacstd2 versus Neutrophils_mMCPcounter falls from 0.527 to 0.014 (study-cluster CI -0.149 to 0.285; shrinkage CI 0.156 to 0.668). Cldn4 versus that neutrophil score, partialled on Tacstd2, stays at 0.401 (CI 0.045 to 0.545). The same mice versus Granulocytes_mMCPcounter do the opposite: Tacstd2 moves from 0.510 to a partial of 0.286 (CI 0.201 to 0.432). Within line, the neutrophil residual correlation is already -0.056. ICI, study unit: Tacstd2 versus responder fraction ρ = -0.059, partial given Cldn4 ρ = 0.055 (CI -0.391 to 0.477, n = 22 series).

小鼠和 series 分开报，不把两个分数说成同一个答案。22 个 series 上 Tacstd2 与 Cldn4 的 ρ = 0.863。series 单位 mMCP 粒细胞：Tacstd2 从 0.860 到偏相关 0.473（区间 -0.219 to 0.818），缩小量 0.387（区间 0.095 to 1.001）。Cldn4 偏掉 Tacstd2 之后是 0.426，两个剩余相关点估计差不多，不能把粒细胞链单独记到 Cldn4 上。CD8 最大的 CIBERSORT 从 0.738 到 0.556，缩小量的区间含 0，偏相关仍为正，不是排斥。小鼠单位 mMCP 中性粒从 0.527 到 0.014，但 mMCP 粒细胞的小鼠偏相关仍是 0.286；扣掉细胞系中位数后中性粒边际相关是 -0.056。ICI 应答在 series 单位上是 -0.059 → 0.055。没有重算 49/64。LLC 不是 KL。

## Units

A mouse is one SRX run. A study is one source series. 20 of 22 baseline series contain a single cell line. GSE124821 pools four mammary lines (91 baseline mice) and GSE159344 pools two melanoma lines (6 mice). The study correlation gives every series one vote. The mouse correlation gives GSE124821 91 of 294 votes. Within-line residuals subtract the line median before the partial correlation, so that row is the within-model association.

The paired table keeps mice measured for both genes. Tacstd2 baseline mice are a subset of the Cldn4 export: 5 Cldn4 baseline mice have no Tacstd2 value (GSE146027 n=2, RU31562203 n=3). RU31562203 (MOC22) is the series that drops out. Cldn4 versus Granulocytes_mMCPcounter on the unpaired Cldn4 export is ρ = 0.869 (n = 23 series). On the paired series it is ρ = 0.852.

Responder labels are constant inside 62 of 64 treated arms. A mouse-level response test would repeat the arm label. It is not used. The mixed arms are listed below and are too small for a partial correlation.

TISMO has no KL or KP lung line. LLC is the only lung carcinoma in this table, and it is not called KL. The locked Tacstd2 49/64 ICB-versus-control count is not recomputed.

## Collinearity

Study-level Spearman of the two genes is 0.863 (n = 22; resample CI 0.599 to 0.950; leave-one-study-out 0.842 to 0.900). Mouse-level Spearman is 0.830 (n = 294; study-cluster CI 0.684 to 0.904). With that shared axis, a partial correlation is a noisy split of one program. The variance inflation from ρ = 0.863 is 3.91.

## Granulocytes

Study unit, Tacstd2 partialled on Cldn4. Granulocytes mMCP: marginal ρ = 0.860 (resample CI 0.672 to 0.953), partial ρ = 0.473 (CI -0.219 to 0.818), shrinkage Δ = 0.387 (CI 0.095 to 1.001). The shrinkage interval stays above zero. The interval for the partial correlation includes zero.

Cldn4 partialled on Tacstd2 is ρ = 0.426 (CI 0.033 to 0.779), from a marginal ρ of 0.852. Leaving out GSE148856 (402230 sarcoma, 8 mice) moves the Tacstd2 partial from 0.473 to 0.239. Leave-one-study-out partials span 0.239 to 0.540, so no single series creates the positive point estimate. The bootstrap interval is wider than that leave-one-out range and includes zero. On z-scored ranks the Tacstd2 coefficient in granulocytes ~ Tacstd2 + Cldn4 is 0.490 (CI -0.221 to 0.893), down from the simple coefficient 0.860. The semipartial correlation, which keeps the granulocyte rank intact, is 0.248.

Studies with at least 8 baseline mice (n = 14): marginal ρ = 0.849, partial ρ = 0.625 (CI -0.351 to 0.967), Δ = 0.224 (CI 0.001 to 1.170).

Neutrophil scores are the rest of the pre-specified granulocyte family. They are not a second discovery set. At the study unit, Neutrophils_mMCPcounter moves from 0.606 to 0.180 (CI -0.246 to 0.463; shrinkage CI 0.005 to 0.863). The shrinkage interval stays above zero. The interval for the partial correlation includes zero.

| Score | n | Tacstd2 ρ | partial given Cldn4 (CI) | shrinkage Δ (CI) | Cldn4 ρ | Cldn4 partial given Tacstd2 |
|---|---|---|---|---|---|---|
| Granulocytes mMCP | 22 | 0.860 | 0.473 (-0.219 to 0.818) | 0.387 (0.095 to 1.001) | 0.852 | 0.426 |
| Neutrophils mMCP | 22 | 0.606 | 0.180 (-0.246 to 0.463) | 0.426 (0.005 to 0.863) | 0.620 | 0.242 |
| Neutrophil TIMER | 22 | 0.328 | 0.134 (-0.354 to 0.608) | 0.194 (-0.275 to 0.623) | 0.306 | 0.047 |
| Neutrophils CIBERSORT | 22 | 0.582 | 0.266 (-0.285 to 0.629) | 0.316 (-0.053 to 0.827) | 0.544 | 0.102 |
| Neutrophil quanTIseq | 22 | 0.509 | 0.014 (-0.436 to 0.384) | 0.495 (0.077 to 1.004) | 0.584 | 0.332 |
| Neutrophil xCell | 22 | 0.531 | 0.300 (-0.164 to 0.642) | 0.231 (-0.108 to 0.671) | 0.459 | 0.003 |

Mouse unit, neutrophils. Neutrophils mMCP: marginal ρ = 0.527 (resample CI 0.275 to 0.616), partial ρ = 0.014 (CI -0.149 to 0.285), shrinkage Δ = 0.512 (CI 0.156 to 0.668). The shrinkage interval stays above zero. The interval for the partial correlation includes zero. Cldn4 versus the same score, partialled on Tacstd2, stays at 0.401 (marginal 0.628; cluster CI 0.045 to 0.545). Dropping GSE124821 leaves the Tacstd2 neutrophil partial at 0.195 (cluster CI 0.041 to 0.312; marginal 0.508, n = 203). Within line, the residual Tacstd2–neutrophil correlation is -0.056 and the partial is -0.076 (n = 294 mice, 17 lines; line-cluster CI for the partial -0.234 to 0.072). The mouse-level attenuation is a between-line fact. There is no within-line Tacstd2–neutrophil link to explain.

Mouse-level granulocytes, the same score as the study result: marginal ρ = 0.510, partial ρ = 0.286 (cluster CI 0.201 to 0.432), shrinkage Δ = 0.224 (CI 0.049 to 0.423). The shrinkage interval stays above zero, and the partial correlation stays on the same side of zero. Cldn4 versus mouse granulocytes, partialled on Tacstd2, is ρ = 0.040 from a marginal ρ of 0.442 (cluster CI -0.132 to 0.235). Within line, the residual granulocyte correlation is 0.210 and the partial is 0.170 (line-cluster CI 0.049 to 0.293; shrinkage CI -0.046 to 0.131). mMCPcounter granulocytes and mMCPcounter neutrophils are not interchangeable in this table.

Mouse unit, Tacstd2 partialled on Cldn4. Intervals are study-cluster bootstrap intervals.

| Score | n | Tacstd2 ρ | partial given Cldn4 (CI) | shrinkage Δ (CI) | Cldn4 ρ | Cldn4 partial given Tacstd2 |
|---|---|---|---|---|---|---|
| Granulocytes mMCP | 294 | 0.510 | 0.286 (0.201 to 0.432) | 0.224 (0.049 to 0.423) | 0.442 | 0.040 |
| Neutrophils mMCP | 294 | 0.527 | 0.014 (-0.149 to 0.285) | 0.512 (0.156 to 0.668) | 0.628 | 0.401 |
| Neutrophil TIMER | 294 | -0.005 | 0.113 (-0.136 to 0.267) | -0.117 (-0.397 to 0.385) | -0.081 | -0.139 |
| Neutrophils CIBERSORT | 294 | 0.358 | 0.304 (0.156 to 0.392) | 0.054 (-0.062 to 0.314) | 0.233 | -0.123 |
| Neutrophil quanTIseq | 294 | 0.372 | 0.083 (-0.041 to 0.397) | 0.289 (-0.049 to 0.410) | 0.397 | 0.171 |
| Neutrophil xCell | 294 | 0.251 | 0.099 (-0.025 to 0.330) | 0.152 (-0.016 to 0.426) | 0.238 | 0.055 |

## CD8

CD8 was pre-specified as the six scores in the earlier TISMO table. Study-level CIBERSORT CD8: marginal ρ = 0.738 (resample CI 0.446 to 0.882), partial ρ = 0.556 (CI 0.044 to 0.800), shrinkage Δ = 0.183 (CI -0.061 to 0.587). The point estimate is closer to zero. The interval for the shrinkage includes zero. Cldn4 versus that score moves from 0.594 to -0.126 after Tacstd2 (CI -0.522 to 0.376; shrinkage CI 0.243 to 1.045). On CIBERSORT, Cldn4’s CD8 link shrinks and Tacstd2’s shrinkage interval includes zero. quanTIseq is the other way around: the Tacstd2 shrinkage interval stays above zero and the partial interval includes zero. TIMER, EPIC, xCell, and mMCPcounter CD8 have Tacstd2 partial intervals that include zero. No pre-specified study-level CD8 score changes from a positive Tacstd2 correlation to a negative partial correlation. This is not CD8 exclusion.

| Score | n | Tacstd2 ρ | partial given Cldn4 (CI) | shrinkage Δ (CI) | Cldn4 ρ | Cldn4 partial given Tacstd2 |
|---|---|---|---|---|---|---|
| CD8 TIMER | 22 | 0.390 | 0.144 (-0.321 to 0.568) | 0.246 (-0.155 to 0.692) | 0.374 | 0.080 |
| CD8 CIBERSORT | 22 | 0.738 | 0.556 (0.044 to 0.800) | 0.183 (-0.061 to 0.587) | 0.594 | -0.126 |
| CD8 EPIC | 22 | 0.336 | 0.115 (-0.476 to 0.687) | 0.220 (-0.262 to 0.698) | 0.325 | 0.075 |
| CD8 quanTIseq | 22 | 0.639 | 0.284 (-0.326 to 0.648) | 0.355 (0.010 to 0.919) | 0.608 | 0.147 |
| CD8 xCell | 22 | 0.120 | 0.144 (-0.396 to 0.506) | -0.025 (-0.414 to 0.503) | 0.054 | -0.098 |
| CD8 mMCP | 22 | 0.666 | 0.433 (-0.194 to 0.722) | 0.234 (-0.072 to 0.811) | 0.563 | -0.032 |

Mouse-level CIBERSORT CD8 is ρ = 0.253, partial ρ = 0.200 (cluster CI 0.049 to 0.399). The within-line residual marginal correlation for that score is -0.067. Mouse xCell CD8 is a null marginal correlation (ρ = -0.006) whose partial is 0.177 (CI 0.003 to 0.283). That is not a marginal CD8 link that failed to attenuate. Mouse CD8 correlations are small once the line is held fixed.

| Score | n | Tacstd2 ρ | partial given Cldn4 (CI) | shrinkage Δ (CI) | Cldn4 ρ | Cldn4 partial given Tacstd2 |
|---|---|---|---|---|---|---|
| CD8 TIMER | 294 | -0.005 | 0.069 (-0.064 to 0.214) | -0.074 (-0.223 to 0.233) | -0.053 | -0.087 |
| CD8 CIBERSORT | 294 | 0.253 | 0.200 (0.049 to 0.399) | 0.053 (-0.109 to 0.349) | 0.172 | -0.070 |
| CD8 EPIC | 294 | 0.006 | -0.057 (-0.270 to 0.360) | 0.063 (-0.130 to 0.349) | 0.045 | 0.072 |
| CD8 quanTIseq | 294 | 0.246 | 0.260 (0.116 to 0.391) | -0.014 (-0.239 to 0.393) | 0.123 | -0.150 |
| CD8 xCell | 294 | -0.006 | 0.177 (0.003 to 0.283) | -0.183 (-0.314 to 0.183) | -0.126 | -0.216 |
| CD8 mMCP | 287 | 0.315 | 0.265 (0.112 to 0.435) | 0.050 (-0.199 to 0.364) | 0.205 | -0.107 |

## Within-line residuals

After subtracting each line's median, the largest |marginal ρ| in the pre-specified list is 0.225 (Tacstd2 vs Neutrophil_quanTIseq, ρ = 0.225, partial 0.217, n = 294 mice across 17 lines). Line medians themselves (n = 17 lines) are a smaller between-model table. Tacstd2 versus granulocytes at the line unit is ρ = 0.664, partial ρ = 0.455 (line-resample CI -0.052 to 0.807). That interval is wide because seventeen lines are collinear on Tacstd2 and Cldn4 (line-level ρ = 0.697).

## ICI response

Treated mice with both genes: 311 mice, 62 pure arms, 2 mixed arms, 22 series. GSE124821 contributes 19 of the 62 pure arms, so an arm-level correlation gives one mammary series nineteen votes. The study-level partial resamples series and is the unit used here.

Study-level Tacstd2 versus responder fraction: ρ = -0.059 (CI -0.484 to 0.436), partial given Cldn4 ρ = 0.055 (CI -0.391 to 0.477), Δ = -0.114 (CI -0.500 to 0.346). The marginal correlation is already near zero, so there is no link to attenuate. Cldn4 versus the same responder fraction is ρ = -0.117, partial given Tacstd2 ρ = -0.115.

Arm-level Tacstd2, descriptive because of the GSE124821 pile-up: ρ = 0.052, partial ρ = 0.217 (series-cluster CI -0.096 to 0.396). Median-split OR for response given high Tacstd2 is 0.866 (Woolf 0.302 to 2.481; table high-R/high-NR/low-R/low-NR = 20/11/21/10). Stratifying that split on the Cldn4 median gives a Mantel–Haenszel OR of 1.143 (Robins–Breslow–Greenland 0.323 to 4.044; stratum tables 15/9/4/3;5/2/17/7). The crude OR and the stratified OR both sit near 1.

Mouse-level labels inside the only mixed arms:

| Gene | Arm | Mice R vs NR | Median R | Median NR | Cliff δ |
|---|---|---|---|---|---|
| Tacstd2 | CT26_GSE139475_antiPD1 | 5 vs 4 | 1.143 | 1.395 | -0.200 |
| Cldn4 | CT26_GSE139475_antiPD1 | 5 vs 4 | 0.766 | 0.890 | -0.400 |
| Tacstd2 | YUMM1.7_GSE103725_BRAF.inhibitor_antiPDL2 | 2 vs 2 | 0.194 | 0.047 | 0.500 |
| Cldn4 | YUMM1.7_GSE103725_BRAF.inhibitor_antiPDL2 | 2 vs 2 | 0.034 | 0.073 | 0.000 |

## What this does not say

- Partialling Cldn4 is an observational split of two correlated bulk measurements. It is not evidence that Cldn4 mediates a Tacstd2 effect, and it is not a knockdown.
- The other immune scores in the TISMO export were not scanned for the largest attenuation.
- Mouse-level analytic p-values are not reported as evidence. Mice are nested in lines and series. The mouse intervals above are study-cluster intervals.
- The locked Tacstd2 49/64 count was not recomputed. LLC is not a KL line. No KP lung line is in this table.
- Cutpoints for the response odds ratio are the medians. They were not tuned.

## Methods

Expression is the TISMO `value` column from `POST /rtismo/gene/downVivoExprn` for Tacstd2 and Cldn4. Infiltrate scores are `POST /rtismo/gene/downICBTreated` with `type=3`, joined on the sample id. The files are the same export used for the study-level granulocyte table. A partial Spearman correlation is the Pearson correlation of ranks after each rank is residualized on the control rank, with an intercept. On z-scored ranks the simple coefficient equals the Spearman correlation, and the coefficient in the two-gene model is `(ρ_xy − ρ_xz ρ_yz) / (1 − ρ_xz²)`. The semipartial correlation correlates the control-residualized exposure rank with the raw outcome rank. Bootstrap intervals are percentile intervals from 2000 resamples. Study rows are resampled as rows. Mouse rows are resampled by series (every mouse from a drawn series), which stops one 91-mouse series from being treated as 91 independent tumors. Within-line intervals resample lines. Leave-one-series-out is the influence check. The response odds ratio is the median split. The stratified odds ratio is the Mantel–Haenszel estimate across the Cldn4 median split, with the Robins–Breslow–Greenland interval. An analytic partial-correlation p-value is stored in the table. It is not the interval interpreted above. The bootstrap interval is. Seed sequence salt is `0` plus the analysis name.

Code: `scripts/tismo_partial_cldn4.py`. Tables: `results/tismo_partial_cldn4/tables/`.

