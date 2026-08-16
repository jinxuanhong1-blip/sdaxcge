# Public lung scRNA: TACSTD2/CLDN4 vs T/NK, split by histology then meta

**Verdict.** In public tumor scRNA, **LUAD** malignant/epithelial TACSTD2 vs T/NK is a weak negative that does **not** reach 0.05 after honest pooling: random-effects ρ = **−0.21** (95% CI −0.44 to +0.05), **p = 0.11**, k = 5, n = 72, I² = 0%. LUAD CLDN4 is null (ρ = +0.02, p = 0.89). **LUSC cannot be meta-analysed**: no series has LUSC n ≥ 6 with this estimand. GSE207422 LUSC is n = 5 (ρ = 0.00 / −0.10). GSE205335 squamous is n = 3. GSE241934 is LUAD/ASC, not LUSC. GSE291670 has no public LUAD/LUSC field.

This is the larger-n scRNA counterpart of the bulk TLS histology interaction. It does **not** recover a LUSC-specific anti-T/NK package. The only adequately powered histology arm is LUAD, and LUAD is not ρ ≈ −0.4 to −0.5.

## 中文摘要

公开肺癌肿瘤 scRNA 按作者组织学拆开后：LUAD 恶性/上皮 TACSTD2 与 T/NK 比例的随机效应 ρ = −0.21（n = 72，p = 0.11），CLDN4 为 0。LUSC **没有** n≥6 的可合并队列。不能把 bulk LUSC TLS 结果写成单细胞已验证。

## Estimand

Sample/patient Spearman of malignant (else epithelial) TACSTD2 or CLDN4 mean vs T/NK fraction. Fisher-z DL meta **within** LUAD and **within** LUSC. Primary n ≥ 6. ASC / SCLC / NUT / unlabeled NSCLC held out.

## Primary LUAD (n ≥ 6)

| Cohort | n | TACSTD2 ρ | p | CLDN4 ρ | p |
| --- | ---: | ---: | ---: | ---: | ---: |
| GSE205335 ADC, ≥20 mal | 14 | −0.23 | 0.44 | −0.15 | 0.60 |
| GSE131907 tLung epi | 11 | +0.09 | 0.79 | −0.05 | 0.89 |
| GSE253013 LUAD malig-like | 9 | −0.72 | 0.030 | −0.33 | 0.38 |
| GSE241934 IIT LUAD | 10 | −0.08 | 0.83 | −0.05 | 0.88 |
| GSE241934 RWC LUAD | 28 | −0.17 | 0.39 | +0.22 | 0.26 |
| **RE meta** | **72** | **−0.21** | **0.11** | **+0.02** | **0.89** |

GSE253013 is the only LUAD cohort with p < 0.05. Removing it (not pre-specified) would make TACSTD2 even closer to zero. I² = 0% for both genes; FE and RE coincide.

## Reported, not pooled

| Cohort | Histology | n | TACSTD2 ρ | CLDN4 ρ |
| --- | --- | ---: | ---: | ---: |
| GSE207422 post, ≥10 mal-like | LUAD | 4 | +0.40 | +0.40 |
| GSE207422 post, ≥10 mal-like | LUSC | 5 | 0.00 | −0.10 |
| GSE205335 SQ, ≥20 mal | LUSC | 3 | — | — |
| GSE241934 IIT/RWC ASC | ASC | 1+1 | — | — |
| GSE131907 tumor-site epi (sensitivity) | LUAD | 36 | +0.17 | −0.21 |
| GSE131907 tumor-site mal (mostly mets) | LUAD | 21 | +0.08 | −0.40 (p=0.075) |
| GSE241934 LUAD IIT+RWC combined | LUAD | 38 | −0.17 | +0.23 |

## Catalogued, not in the ρ table

| Series | Why |
| --- | --- |
| GSE291670 | TACSTD2/CLDN4 vs T/NK exist (n=6) but GEO histology is **NSCLC**, not LUAD/LUSC |
| GSE148071 | Mixed LUAD/LUSC in the paper; GEO SOFT has age/sex only |
| GSE154826 | CD45-enriched; compact objects lack usable epithelium |
| LuCA / HLCA integrated | Public but >2 GB objects; skipped |

## Read against the bulk TLS interaction

Bulk: TACSTD2–TLS/B/CD8 after purity **holds in TCGA-LUSC**, not as TLS in TCGA-LUAD.  
scRNA: the large public n is **LUAD** (n=72) and is **not** a significant anti-T/NK correlation. LUSC scRNA n is 5+3. The histology interaction seen in bulk **cannot be tested at single-cell with the open matrices listed here**.

## Figures

- `figures/forest_TACSTD2_by_histology.png`
- `figures/forest_CLDN4_by_histology.png`

Open circles = n < 6, shown, not pooled. Numbers from `cohort_effects.tsv` and `meta_by_histology.tsv`.
