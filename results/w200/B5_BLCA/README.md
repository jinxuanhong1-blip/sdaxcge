# B5_BLCA — CLDN4-high vs ICI objective response, urothelial carcinoma

**Verdict: claim NOT supported (`inconsistent`).**
Claimed OR = 0.42. Observed pooled OR = **1.31 (95 % CI 0.82–2.10)** across
n = 406 response-evaluable patients in 3 public cohorts. Opposite direction;
0.42 is rejected at **p = 2.1 × 10⁻⁶**. `did_we_tune_to_0.42: false`.

Read `WRITEUP.md` (bilingual EN/中文) for the full report and `analysis_plan.md`
for the prespecification, which was committed before any result was computed.

## Headline numbers

| | CLDN4-high ORR | CLDN4-low ORR | OR (95 % CI) | p |
|---|---|---|---|---|
| IMvigor210 (n = 298) | 26.2 % | 19.5 % | 1.47 (0.82–2.64) | 0.214 |
| BACI / GSE176307 (n = 87) | 18.6 % | 18.2 % | 1.03 (0.30–3.53) | 1.000 |
| Snyder 2017 (n = 21) | 30.0 % | 36.4 % | 0.76 (0.08–6.52) | 1.000 |
| **Pooled, random effects** | — | — | **1.31 (0.82–2.10)** | 0.258 |

I² = 0 %, τ² = 0, Q = 0.71 (df 2, p = 0.700).

Why the null is informative rather than just underpowered:

- **Positive control** (CD8 T-effector, published to track response in IMvigor210):
  pooled OR 1.76 (1.09–2.83), p = 0.021; per SD in IMvigor210 OR 1.46, p = 0.0069.
- **Negative control** (housekeeping genes): pooled OR 0.88 (0.55–1.41), p = 0.592.
- **Power** to detect OR = 0.42 had it been real: 0.70 (20 000 simulations).

Every prespecified sensitivity analysis agrees — five cutoff rules, continuous
exposure, Mann–Whitney, ITT endpoint, bladder-site-only, DESeq normalisation,
an independent processing pipeline, and multivariable adjustment. None lands near
0.42.

## Files

| File | Contents |
|---|---|
| `analysis_plan.md` | Prespecification, committed before results |
| `WRITEUP.md` | Full bilingual report (EN + 中文) |
| `summary.json` | Machine-readable results, incl. `honest_verdict` block |
| `cohort_effects.csv` | Per-cohort primary 2×2 tables and ORs |
| `meta_summary.csv` | Pooled estimates for all four models + heterogeneity |
| `cohort_flow.csv` | n at every stage: available → RNA → evaluable → analysed |
| `validation.json` | Checks of the mirrored IMvigor210 against published values |
| `controls.csv`, `controls_continuous.csv` | Positive/negative pipeline controls |
| `sensitivity_cutoffs.csv` | S1 cutoff robustness |
| `continuous_logistic.csv`, `mannwhitney.csv` | S2, S3 |
| `exploratory_genes.csv` | Exploratory gene family with BH q-values (m = 11) |
| `per_sample_expression.csv` | Auditable per-patient table with `expression_scale` |
| `download_manifest.json` | URLs, byte counts, SHA-256, declared scales |
| `fig_forest_meta.png` | Forest plot with the claimed 0.42 marked |
| `fig_orr_by_cldn4.png` | ORR by CLDN4 median split, with counts |
| `fig_cldn4_by_response.png` | CLDN4 distribution in responders vs non-responders |
| `fig_cutoff_robustness.png` | Pooled OR across five cutoff rules |
| `fig_controls.png` | Controls vs the CLDN4 primary |

## Reproduce

```bash
python3 scripts/w200/B5_BLCA/download.py         # ~125 MB, writes manifest + catalog.tsv
python3 scripts/w200/B5_BLCA/prepare_cohorts.py  # runs extract_imvigor210.R, harmonises cohorts
python3 scripts/w200/B5_BLCA/analyze.py          # statistics, figures, summary.json
```

Requires Python 3 with pandas/scipy/statsmodels/matplotlib, and R (used only to
read the archived `cds.RData` S4 object).

## Data sources

- **IMvigor210** — Mariathasan et al., *Nature* 2018 (doi:10.1038/nature25501),
  atezolizumab in metastatic urothelial carcinoma. The official
  `IMvigor210CoreBiologies` URL now returns HTTP 404; `cds.RData` was taken from a
  GitHub mirror and accepted only after reproducing the published RECIST
  (CR 25 / PR 43 / SD 63 / PD 167 / NE 50) and immune-phenotype
  (76/134/74/64) distributions exactly, with CLDN4 agreeing at ρ = 0.991 with an
  independently harmonised copy of the same trial.
- **BACI** — GEO `GSE176307`, Robertson et al. 2021, real-world metastatic
  urothelial carcinoma on atezolizumab/pembrolizumab/other ICI.
- **Snyder 2017** — *PLoS Med* 14:e1002309, via the bhklab PredictIO compendium
  (Zenodo 7058399; harmonised by Bareche et al., *Ann Oncol* 2022, PMID 36055464).

## Scope

Claim B5 refers to 11 cohorts across GBM/NSCLC/melanoma/RCC/urothelial. This slice
tests **only the urothelial arm**, and it cannot by itself confirm or refute a
pooled 11-cohort figure. Separately, no published CLDN4 × ICI meta-analysis could
be located, so OR = 0.42 is treated as an unsourced assertion under test rather
than a result being reproduced. We do not name cohorts we cannot enumerate.
