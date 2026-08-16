# TACSTD2 and CLDN4 vs DCB in GSE190265 / GSE190266

**Verdict: no support for TACSTD2. CLDN4 is a weak, unreplicated trend (unadjusted p = 0.052 in France4 only). Do not treat either gene as a DCB marker in these cohorts.**

Anti-PD1 monotherapy, advanced NSCLC, CGFL Dijon (Limagne / Ghiringhelli; GSE190270 sub-series; PMID 35051357). Tumor RNA-seq TPM from GEO. DCB is derived from deposited PFS, not from RECIST labels (none were uploaded).

## What was tested

| item | rule |
|---|---|
| DCB | PFS time ≥ 6 months |
| NDB | progression/death event at PFS < 6 months |
| excluded | censored before 6 months (benefit unknown) |
| expression | provider TPM, `log2(TPM+1)` |
| primary test | two-sided Wilcoxon rank-sum, DCB vs NDB |
| pooled test | cohort-stratified van Elteren (not a raw TPM merge) |
| multiplicity | Benjamini–Hochberg over the four gene × cohort Wilcoxon tests |

The 6-month cutoff is the conventional ICI DCB rule. It was not chosen to fit these genes.

## Data limitations (read these first)

1. **GSE190266 TACSTD2 is missing.** The deposited `GSE190266_TPM_France4.csv.gz` has 16,383 gene columns and stops at `MTMR14`. That is exactly an Excel-style 2^14 column cap. TACSTD2 sorts after T and is not in the file. Raw SRA is available but was not re-quantified here. **Any TACSTD2 claim from France4 would be fabricated.**
2. **GSE190266 PFS is already capped at 6 months** (`pfs_time (6 months)` / `pfs_evt (6 months)` in the series matrix). DCB in that cohort means “still progression-free at the 6-month cap.” Continuous PFS tests there are biased and are reported only as sensitivity.
3. **GSE190265 has uncapped PFS** in `GSE190265_samples_info_France3.csv.gz` (max 46.4 months). That is the only cohort where Spearman / log-rank vs time is interpretable.
4. **Platforms differ** (France3: HiSeq/NextSeq, 33,976 genes; France4: NovaSeq, truncated matrix). TPM values are not on a common scale.
5. **No RECIST response field** is deposited. DCB ≠ ORR.
6. **Histology is incomplete** on GEO (France3: 18/43 labeled unknown; France4: 2 unknown, plus one `disease state: 5` treated as missing).

## Cohort composition

| cohort | platform | TPM n | DCB | NDB | excluded |
|---|---|---|---|---|---|
| GSE190265 France3 | GPL18573 | 43 | 14 | 29 | 0 |
| GSE190266 France4 | GPL24676 | 70 | 17 | 52 | 1 (`A16-1489-1`, PFS 0.03 mo, event 0) |

## Primary results

| cohort | gene | n DCB / NDB | median TPM DCB vs NDB | AUC | Wilcoxon p | BH p |
|---|---|---|---|---|---|---|
| GSE190265 | TACSTD2 | 14 / 29 | 72.5 vs 61.1 | 0.46 | 0.67 | 0.93 |
| GSE190265 | CLDN4 | 14 / 29 | 22.4 vs 17.1 | 0.49 | 0.93 | 0.93 |
| GSE190266 | TACSTD2 | — | **not in GEO matrix** | — | — | — |
| GSE190266 | CLDN4 | 17 / 52 | 55.9 vs 21.7 | 0.66 | **0.052** | 0.15 |

AUC is P(a random DCB sample has higher expression than a random NDB sample). 0.5 = no difference.

**TACSTD2 (France3 only):** medians almost overlap, AUC on the wrong side of 0.5, p = 0.67. Spearman vs uncapped PFS ρ = −0.03 (p = 0.84). Median-split log-rank p = 0.88. **Null.**

**CLDN4:** null in France3 (AUC 0.49, p = 0.93). In France4 the DCB median is higher (AUC 0.66) but the two-sided Wilcoxon is 0.052 and BH-adjusted 0.15. Spearman vs *capped* PFS ρ = 0.20 (p = 0.10). Median-split log-rank p = 0.31. **A one-cohort trend, not a finding.**

## Pooled (CLDN4 only is real pooling)

| test | TACSTD2 | CLDN4 |
|---|---|---|
| van Elteren Z (p) | −0.44 (0.66) — France3 only | +1.67 (0.095) |
| within-cohort z-score Wilcoxon p | 0.67 | 0.27 |

Pooling does not rescue CLDN4. The France3 null dilutes the France4 trend.

## Histology (confounder check, exploratory)

DCB is almost entirely non-squamous in France4 (15/17). CLDN4 is slightly higher in non-squamous tumors, not significantly (France3 p = 0.12; France4 p = 0.35).

Restricted to non-squamous (tiny n in France3):

| cohort | gene | n DCB / NDB | AUC | p |
|---|---|---|---|---|
| GSE190265 | TACSTD2 | 4 / 10 | 0.55 | 0.84 |
| GSE190265 | CLDN4 | 4 / 10 | 0.85 | 0.054 |
| GSE190266 | CLDN4 | 15 / 39 | 0.66 | 0.066 |

Same pattern: CLDN4 point estimates favor DCB, every p > 0.05, France3 non-squamous DCB n = 4. This is not confirmation.

## Honest conclusion

- **TACSTD2 / TROP2:** no association with 6-month DCB in the only cohort where the gene is actually present.
- **CLDN4:** one borderline unadjusted p-value in France4 that fails BH correction, fails log-rank, and is absent in France3. Exploratory non-squamous splits stay above 0.05.
- These datasets are small (14 and 17 DCB events). They can rule out a *large* effect for TACSTD2 in France3; they cannot prove a small CLDN4 effect.
- Do not cite this folder as evidence that TROP2 or claudin-4 expression predicts PD-1 benefit.

## How to reproduce

```bash
python3 scripts/fetch_geo.py
python3 scripts/analyze_tacstd2_cldn4_dcb.py
```

Outputs in this directory: `per_cohort_stats.csv`, `stratified_van_elteren.csv`, `pooled_zscored_stats.csv`, `nonsquamous_only_stats.csv`, `expression_long.csv`, `clinical_dcb_labels.csv`, `TACSTD2_CLDN4_DCB_boxplots.png`.
