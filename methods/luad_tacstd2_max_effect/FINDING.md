# OncoSG + GEO LUAD: maximum story-fit TACSTD2 |ρ| and TJ/junction NES

Public bulk LUAD only. Every ρ and NES below is written by `analyze.py` from the downloaded matrices. The maximum is the extreme of a pre-declared grid, so its p-value is smaller than a single pre-specified test. Sweep q is Benjamini–Hochberg across that grid.

Story-fit rule, fixed before the maximum was read: keep TACSTD2–CD8/ImmuneScore correlations only when ρ is negative, and keep TJ/junction enrichments only when NES is positive in TACSTD2-high. GSE19804 (tissue labeled lung cancer) and GSE233774 (n = 30) are computed and excluded from the headline. GSE10072, GSE11969, and GSE248378 stay closed.

## Maximum |ρ|

Eligible cells: 192 of 236 primary-family tests (612 rows stored, including companions and ineligible cohorts). Eligible |ρ| median 0.191, 95th percentile 0.419.

**Search maximum:** GSE273377_validation, CD8A vs TACSTD2, subset `all`, covariate `ESTIMATE_cosine_TumorPurity` (cosine_purity). n = 60. Spearman ρ = **-0.443** (p = 0.000443, primary-family q = 0.00227). Genes in the endpoint: 1 (CD8A).

In that same cohort the unadjusted CD8A ρ is -0.423 and the stromal9 partial is -0.442. The search maximum is 0.020 above the unadjusted CD8A ρ. Cosine TumorPurity and ESTIMATEScore give the same partial ρ in this cohort, so the cosine transform is monotone with ESTIMATEScore here and is the same adjustment. OncoSG high-tumor cytotoxic6 (CD8A, GZMA, GZMB, PRF1, NKG7, GNLY z-mean), no further covariate, is ρ = -0.438 (n = 87, q = 0.000169).

**Cross-cohort maximum:** cytotoxic6_zmean, subset `all`, covariate class `none`. Pooled ρ = **-0.308** (DL random, 95% CI -0.361 to -0.253, p = 1.15e-25, I² = 15%, k = 6, n sum = 1398). Cohorts: GSE273377_validation,GSE282774,GSE31210,GSE68465,GSE72094,OncoSG.

Locked single-gene reference, recomputed in this run. The partial column is the composition covariate: published PURITY on OncoSG, ESTIMATE StromalScore on arrays, and the 9-gene stromal z-mean on RNA-seq GEO.

| cohort | n | CD8A unadjusted ρ | CD8A composition partial ρ | immune8 unadjusted ρ | immune8 composition partial ρ |
|---|---:|---:|---:|---:|---:|
| OncoSG | 169 | -0.380 | -0.309 | -0.387 | -0.318 |
| GSE31210 | 226 | -0.289 | -0.194 | -0.302 | -0.208 |
| GSE72094 | 442 | -0.227 | -0.153 | -0.228 | -0.156 |
| GSE68465 | 443 | -0.152 | -0.166 | -0.226 | -0.250 |
| GSE273377_discovery | 103 | 0.039 | -0.028 | 0.096 | 0.030 |
| GSE273377_validation | 60 | -0.423 | -0.442 | -0.330 | -0.355 |
| GSE282774 | 58 | -0.206 | -0.205 | -0.179 | -0.178 |
| GSE19804 | 60 | -0.024 | 0.104 | -0.133 | -0.043 |
| GSE233774_tumor | 30 | -0.437 | -0.266 | -0.439 | -0.268 |

Array partials on ESTIMATEScore (stromal ssGSEA + immune ssGSEA), the GSE31210 covariate in PR 736: GSE31210 -0.107; GSE72094 -0.084; GSE68465 -0.170. Cosine TumorPurity matches those three partials exactly in this run.

CD8A, all samples, unadjusted, pooled across cohorts where that correlation is negative: ρ = -0.263 (k = 6, I² = 54%). immune8 z-mean on the same slice: ρ = -0.265 (k = 6, I² = 11%).

GSE273377 discovery (n = 103) is strict LUAD and is left out of the negative pools where its sign is positive. All-sample unadjusted ρ: CD8A +0.039; cytotoxic6_zmean +0.166; immune8_zmean +0.096; estimate_immune141 -0.039. The pooled ρ values are therefore pooled over the cohorts that stay negative.

Ayers GEP18 is a companion, outside the CD8/ImmuneScore maximum. Its largest negative ρ on strict LUAD with n ≥ 50 is -0.432 (OncoSG, high_tumor, none, n = 87).

## Maximum TJ/junction NES in TACSTD2-high

Library: 19 name-filtered sets (`data/junction_sets.gmt`). Tests in the strict-LUAD n≥50 family: 467. Positive NES (story-fit): 419.

**Positive-NES maximum:** GSE273377_validation, `A8::HALLMARK_APICAL_JUNCTION`, contrast `tercile_welch`. NES = **2.893** (ES = 0.557, nominal p = 0.000999, family q = 0.00161, genes in rank = 200, arms 20/20). Leading edge: ACTB,CDH3,CLDN4,CDH1,JUP,MAPK13,MYH9,VCL,ADAM9,TGFBI,ACTN4,SDC3,MYL12B,NECTIN1,MMP2,VAV2,CAP1,ITGB1,CD276,CNN2.

That Hallmark apical-junction specification across strict LUAD: GSE273377_discovery +2.25; GSE273377_validation +2.89; GSE282774 +1.64; GSE31210 -1.23; GSE68465 +2.09; GSE72094 +1.71; OncoSG -1.59. OncoSG and GSE31210 are negative, so this cell is the positive-NES extreme and does not carry a cross-cohort apical-junction claim. Locked KEGG tight junction, Q4 vs Q1 Welch, recomputed here: OncoSG -1.01; GSE31210 +1.94 (PR 736 reported OncoSG −1.00 and GSE31210 +2.00 under the same ranking statistic and an independent permutation stream).

**Locked A8 set with the largest positive NES** (KEGG tight junction, GO tight-junction organization, GO bicellular tight-junction assembly, Hallmark apical junction): GSE273377_validation, `A8::HALLMARK_APICAL_JUNCTION`, `tercile_welch`, NES = **2.893** (nominal p = 0.000999, q within the four-set family = 0.00145).

**Positive in all 7 strict LUAD cohorts:** `ENRICHR_GOCC2023::Apical Junction Complex (GO:0043296)`, `q4_vs_q1_welch`, largest NES = **2.535** in GSE273377_validation (nominal p = 0.000999, family q = 0.00161). The smallest NES in that specification is 1.210 in OncoSG. The specification with the highest floor is `ENRICHR_REACTOME2022::Tight Junction Interactions R-HSA-420029`, `spearman_vs_tacstd2`, minimum NES **1.823** (OncoSG) and maximum 2.271 (GSE282774).

A8 sets at the Q4 vs Q1 Welch contrast (positive NES is the story direction):

| cohort | KEGG tight junction | GO TJ organization | GO bicellular TJ assembly | Hallmark apical junction |
|---|---:|---:|---:|---:|
| OncoSG | -1.006 | 1.069 | 1.477 | -1.501 |
| GSE31210 | 1.941 | 1.787 | 1.717 | -1.158 |
| GSE72094 | 1.920 | 1.943 | 1.880 | 1.833 |
| GSE68465 | 2.012 | 1.804 | 1.827 | 2.148 |
| GSE273377_discovery | 2.574 | 2.148 | 1.985 | 2.398 |
| GSE273377_validation | 2.370 | 2.436 | 2.243 | 2.761 |
| GSE282774 | 2.221 | 2.097 | 2.029 | 1.643 |

OncoSG Q4 vs Q1 for the four locked sets: GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY NES 1.48; GOBP_TIGHT_JUNCTION_ORGANIZATION NES 1.07; HALLMARK_APICAL_JUNCTION NES -1.50; KEGG_TIGHT_JUNCTION NES -1.01. 2 of 4 are positive in TACSTD2-high. Sets at or below zero: HALLMARK_APICAL_JUNCTION, KEGG_TIGHT_JUNCTION.

## What this is

These are bulk RNA associations in surgical or resected LUAD. They are not a spatial exclusion measurement and not an ICI-response result. OncoSG remains the East-Asian z-score matrix with n = 169 (portal RNA list 181 is not the matrix n). ESTIMATE is not applied to OncoSG z-scores. Selecting the grid maximum inflates |effect| relative to the locked CD8A and KEGG rows above.

## Reproduce

```bash
pip install -r methods/luad_tacstd2_max_effect/requirements.txt
python3 methods/luad_tacstd2_max_effect/analyze.py
```

