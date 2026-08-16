# B5 Liu 2019: CLDN4 expression versus anti-PD-1 response

## Bottom line

There is **no evidence that higher pretreatment tumor `CLDN4` RNA predicts
response** in this cohort. In the paper's primary comparison (CR/PR versus PD),
responders had slightly lower, not higher, `CLDN4`, but the difference was
small and compatible with no association:

- Responders: n = 47, median = 0.0261 TPM
- Progressors: n = 56, median = 0.0356 TPM
- Mann–Whitney U = 1178, two-sided p = 0.347
- Rank-biserial correlation = -0.105 (bootstrap 95% CI -0.315 to 0.114)
- AUC if higher `CLDN4` predicts response = 0.448

`CLDN4` was low and zero-inflated: it was detected above zero in 26/47
responders and 36/56 progressors (Fisher p = 0.421). The data therefore do not
support `CLDN4` as a standalone response biomarker here.

## Cohort and analysis

The analysis uses the authors' public patient-by-gene TPM matrix and clinical
Supplementary Table 1 from Liu, Schilling, et al. Response is best RECIST
response to anti-PD-1. The prespecified primary contrast follows the paper:
CR/PR ("responders") versus PD ("progressors"), excluding 16 SD and 2 MR cases
among the 121 RNA-profiled patients.

The two-sided Mann–Whitney test is appropriate for the strongly skewed,
zero-inflated expression values. The test on `log2(TPM + 1)` has the same ranks
and p-value as a test on TPM. The computed p-value exactly reproduces the
published `CLDN4` entry in Supplementary Table 4.

Sensitivity analyses do not change the main conclusion. The conventional
CR/PR-versus-PD/SD/MR contrast gives p = 0.216. Ipilimumab-stratified and
skin/occult analyses are also null. A post hoc skin-only subset gives a nominal
p = 0.027 in the opposite direction (lower `CLDN4` in responders), but it does
not survive Holm adjustment over the five listed sensitivity checks
(adjusted p = 0.136). It should not be treated as a validated subgroup finding.

## Files

- `cldn4_vs_response.png` / `.svg`: primary comparison
- `primary_result.json`: primary test, effect size, and bootstrap intervals
- `sensitivity_results.csv`: alternative outcome and subgroup checks
- `cldn4_patient_data.csv`: public sample-level values used in the analysis
- `source_manifest.json`: source URLs, SHA-256 hashes, software, and validation
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
