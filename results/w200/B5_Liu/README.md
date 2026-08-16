# B5 Liu 2019: CLDN4 expression versus anti-PD-1 response

## Bottom line

Open Liu 2019 melanoma anti-PD-1 RNA cohort. **No evidence that higher
pretreatment `CLDN4` predicts response.**

Primary OR/n/p (CR/PR vs PD; `CLDN4` above vs at-or-below the primary-cohort
median TPM):

- **OR = 0.70**
- **n = 103** (47 responders, 56 progressors; 2x2 = 21/30/26/26)
- **p = 0.431**
- 95% CI 0.33–1.52

The point estimate is in the opposite direction from a positive biomarker
(high `CLDN4` had *lower* response odds). The interval includes 1. Detection
and continuous logistic models are the same story.

## OR / n / p

| analysis | n | OR | 95% CI | p |
|---|---:|---:|---|---:|
| Primary median-split high vs low | 103 | 0.700 | 0.326–1.524 | 0.431 |
| Primary detected vs undetected | 103 | 0.688 | 0.316–1.518 | 0.421 |
| Primary top-quartile vs rest | 103 | 0.676 | 0.282–1.677 | 0.496 |
| Primary logistic per +1 log2(TPM+1) | 103 | 0.931 | 0.302–2.864 | 0.900 |
| All RNA median-split (CR/PR vs PD/SD/MR) | 121 | 0.630 | 0.306–1.318 | 0.264 |

2x2 cells are high-CLDN4 responder / high-CLDN4 nonresponder / low-CLDN4
responder / low-CLDN4 nonresponder. Odds ratios are Fisher exact except the
logistic row. Confidence intervals for 2x2 tests use the Haldane–Anscombe
correction.

Ipilimumab-stratified and melanoma-subtype median-splits are also null (all
OR 95% CIs include 1; all Holm-adjusted p = 1). The lowest point estimate is
skin/occult (OR = 0.51, n = 91, p = 0.144). That is not a validated subgroup
finding.

## Continuous test

The paper's published single-gene test is two-sided Mann–Whitney on expression
in CR/PR versus PD. That result is reproduced exactly:

- Responders: n = 47, median = 0.0261 TPM
- Progressors: n = 56, median = 0.0356 TPM
- U = 1178, two-sided p = 0.347
- Rank-biserial = -0.105 (bootstrap 95% CI -0.315 to 0.114)

`CLDN4` is low and zero-inflated (detected in 26/47 responders and 36/56
progressors). The data do not support `CLDN4` as a standalone ICI-response
biomarker in this cohort.

## Cohort

Public patient-by-gene TPM matrix and Supplementary Table 1 from Liu,
Schilling, et al. Response is best RECIST response to anti-PD-1. The
prespecified primary contrast follows the paper: CR/PR versus PD, excluding
16 SD and 2 MR cases among 121 RNA-profiled patients.

## Files

- `or_n_p.csv` / `or_n_p.md`: odds ratio, n, p, 2x2 counts, and CIs
- `or_n_p.png` / `.svg`: forest plot of binary ORs
- `cldn4_vs_response.png` / `.svg`: primary expression comparison
- `primary_result.json`: headline OR/n/p plus the Mann–Whitney result
- `sensitivity_results.csv`: alternative outcome and subgroup rank tests
- `cldn4_patient_data.csv`: public sample-level values
- `source_manifest.json`: source URLs, SHA-256 hashes, software, validation
- `analyze_cldn4.py`: complete reproducible analysis

## Reproduce

```bash
python3 -m pip install -r requirements.txt
python3 analyze_cldn4.py
```

The script downloads the two public supplements to
`/tmp/liu2019_cldn4_source`, verifies their SHA-256 hashes, and regenerates all
tables and figures. It uses a fixed seed and 20,000 stratified bootstrap
replicates.

## Source

Liu D, Schilling B, Liu D, et al. *Integrative molecular and clinical modeling
of clinical outcomes to PD1 blockade in patients with metastatic melanoma.*
Nature Medicine. 2019;25:1916–1927.
[doi:10.1038/s41591-019-0654-5](https://doi.org/10.1038/s41591-019-0654-5).

This is a single-cohort, retrospective association analysis. It is not a
clinical validation study, and a null result does not establish absence of a
small effect.
