# A4 TISMO colon/CRC — Tacstd2 after ICB

**Claim (user, all-cancer):** Tacstd2 up in 49/64 ICI-treated mouse models (p=5.8e-5).

**This slice:** colon/CRC models only. Filters were not tuned to recover 49/64.

**Honest verdict: `8/10` up. The 49/64 statistic is not a colon result.**

NOT the 49/64 result. Colon-only is **8/10** cohorts up. Paired t / Wilcoxon are only barely <0.05 if the 10 cohorts are treated as independent; the sign test, study collapse, cell-line collapse, CT26-only, MC38-only, and strict-ICB subsets are all null.

## Primary contrast

TISMO in-vivo Gene module, ICB regimens vs the study's own `Baseline=1` control arm.
One pair per TISMO cohort (`cellline_study_condition_regimen`).
Values are TISMO quantile-normalised, ComBat-corrected log-scale TPM as served.

| Subset | n cohorts | n lines | n studies | up/n | mean Δ | paired t p | Wilcoxon p | sign-test p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| colon ICB cohorts (TISMO definition) | 10 | 2 | 5 | 8/10 | 0.160 | 0.0479 | 0.0488 | 0.109 |
| CT26 only | 7 | 1 | 3 | 6/7 | 0.213 | 0.0627 | 0.109 | 0.125 |
| MC38 only | 3 | 1 | 2 | 2/3 | 0.036 | 0.397 | 0.5 | 1 |
| strict ICB regimens only | 6 | 2 | 5 | 5/6 | 0.115 | 0.322 | 0.219 | 0.219 |
| ICB plus other partner | 4 | 2 | 2 | 3/4 | 0.227 | 0.0697 | 0.25 | 0.625 |
| collapsed to one pair per study | 5 | 2 | 5 | 4/5 | 0.084 | 0.517 | 0.438 | 0.375 |
| collapsed to one pair per cell line | 2 | 2 | 5 | 2/2 | 0.124 | 0.394 | 0.5 | 0.5 |

Primary n = **10** paired cohorts, **46** baseline + **50** ICB samples, **2** cell lines, **5** studies.

- mean baseline → ICB: **0.406 → 0.565** (mean Δ = 0.160)
- paired t = 2.29, p = 0.0479
- Wilcoxon p = 0.0488
- sign test (8 up / 2 down vs 0.5) p = 0.109
- Cohen's d_z = 0.72

User-reported 49/64 and p=5.8e-05 are **all-cancer**. This download recovers that all-cancer count (49/64 up). It is not a colon result.

## Per-cohort table (every colon pair TISMO serves)

| Model | Study | Regimen | n B/ICB | mean B | mean ICB | Δ | dir | strict ICB |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| CT26 | ERP114266 | antiPDL1; antiCTL4; | 5/5 | 0.243 | 0.656 | +0.413 | up | yes |
| CT26 | ERP114266 | antiPDL1; antiCTL4; | 5/5 | 1.215 | 1.533 | +0.318 | up | yes |
| CT26 | GSE139475 | antiPD1 | 5/9 | 1.629 | 1.313 | -0.317 | down | yes |
| CT26 | GSE153239 | antiGARP:TGFb FcD; antiPD1 FcS; | 5/5 | 0.262 | 0.620 | +0.358 | up | no |
| CT26 | GSE153239 | antiGARP:TGFb WT; antiPD1 FcS; | 5/5 | 0.189 | 0.466 | +0.276 | up | no |
| CT26 | GSE153239 | antiPD1_FcS | 5/5 | 0.220 | 0.375 | +0.155 | up | yes |
| CT26 | GSE153239 | antiTGFb; antiPD1 FcS; | 5/5 | 0.270 | 0.555 | +0.285 | up | no |
| MC38 | GSE172162 | antiPD1 | 4/4 | 0.016 | 0.036 | +0.020 | up | yes |
| MC38 | GSE172162 | antiPD1; exercise | 4/4 | 0.013 | 0.000 | -0.013 | down | no |
| MC38 | GSE93017 | antiPDL1 | 3/3 | 0.000 | 0.100 | +0.100 | up | yes |

`strict ICB` = treated arm is only anti-PD-1 / PD-L1 / CTLA-4 (no TGF-β, GARP, exercise, or other partners).

## Independence (why p≈0.05 is not a colon confirmation)

Cohorts cluster inside two cell lines and five studies. Treating n=10 as independent overstates precision.

| Cluster | n clusters | mean Δ | 95% bootstrap CI | two-sided p |
| --- | ---: | ---: | --- | ---: |
| cell_line | 2 | 0.160 | NA | NA | k=2 clusters: bootstrap p/CI not reported (resamples collapse to the cluster means)
| study | 5 | 0.160 | [-0.074, 0.285] | 0.0846 |

Cell-line k=2 is too small for a bootstrap p. Study-level k=5 is the usable cluster check.

- Collapsed to one pair per **study** (n=5): 4/5 up; paired t p and Wilcoxon p are null (see table).
- Collapsed to one pair per **cell line** (n=2): both lines go up (CT26 Δ larger; MC38 near zero). n=2 cannot support a p-value.
- **MC38** Tacstd2 is essentially off (cohort means 0.00–0.10). A 2/3 “up” count there is noise around zero.
- **CT26 only** (n=7, 6/7 up): paired t p ≈ 0.06, Wilcoxon / sign test null.

## Coverage (honest missingness)

- `cellLineMeta` colorectal carcinoma models: 1638N-T1, CMT93, CT26, MC38.
- Of those, the ICB Gene-module model list contains: CT26, MC38.
- Colorectal models with **no** ICB Gene-module expression: 1638N-T1, CMT93.
- Samples appearing in more than one colon cohort: 0 (should be 0).
- All-cancer Tacstd2 ICB cohorts TISMO returned in this download: 64 (49 up / 15 down). Colon is 10 of those 64, not a separate 64.

## What this is not

- A reproduction of 49/64 or p=5.8e-5 in colon. That numerator/denominator is all-cancer.
- Evidence that Tacstd2 induction after ICB is a CRC-specific finding.
- Human CRC, MSI/MSS, or TROP2-ADC outcome.
- A Cldn4 analysis.

## Caveats

1. TISMO units are site-normalised, not raw TPM. Direction within a study is usable; absolute values are not portable.
2. Several CT26 arms are ICB + TGF-β / GARP. The strict-ICB subset is null (see table).
3. ERP114266 day-7 and day-14 are the same study at two time points.
4. n=10 is small. A p-value that sits on 0.05 under the most liberal pairing is not a confirmation.
5. 1638N-T1 and CMT93 are colorectal in TISMO metadata but are absent from the ICB Gene module, so they cannot enter the denominator.

## Rerun

```bash
python3 scripts/w200/A4_colon/download.py
python3 scripts/w200/A4_colon/analyze.py
```

See `summary.json`, `paired_summaries.csv`, and `cohort_pairs.csv`.
