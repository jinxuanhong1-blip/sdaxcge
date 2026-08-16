# B5 analog: open melanoma ICI RNA, CLDN4 versus response

## Bottom line

Public melanoma ICI transcriptomes **do not support** the B5 claim that
CLDN4-high patients have worse ICI response (user-reported pooled OR = 0.42).

The named analog, **Riaz 2017 / GSE91061 pretreatment**, is null:
responder-higher AUC 0.487 (bootstrap 95% CI 0.328–0.649), Mann–Whitney
p = 0.911, median-split OR 1.52 (Fisher p = 0.725). That is no
discrimination, not a negative-predictive signal.

The two other open melanoma RNA+response cohorts that are large enough to
test are also non-significant. A Mantel–Haenszel pool of Riaz + Hugo + Liu
gives OR 0.82 (0.47–1.44). The interval includes 1 and excludes 0.42.

## What was tested

B5 (user PPT) states that across an 11-cohort ICI meta-analysis, CLDN4-high
is associated with worse response (OR = 0.42). This folder tests the
**open melanoma** slice of that claim: public RNA plus a reconstructable
response label. It does not reconstruct the unpublished 11-cohort list.

Primary analog: Riaz et al., *Cell* 2017, GEO
[GSE91061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061),
nivolumab, pretreatment biopsy, GEO `PRCR` vs `SD`/`PD`.

Added open melanoma ICI RNA cohorts:

| Cohort | Access | Endpoint | Result |
|---|---|---|---|
| Riaz 2017 GSE91061 | GEO FPKM + response | Pre, PRCR vs SD/PD | Null |
| Hugo 2016 GSE78220 | GEO FPKM + RECIST | Pre, CR/PR vs PD | Null |
| Liu 2019 *Nat Med* | GitHub raw counts + Nature Table 1 | Pre-PD1 RNA, CR/PR vs SD/PD | NS trend toward lower CLDN4 in responders |
| Auslander 2018 GSE115821 | GEO counts + R/NR | Pre, patient-collapsed | n=7 (2 responders); not interpretable |

## Riaz primary result

| Analysis | R / NR | AUC (R higher) | 95% bootstrap CI | MW p | CLDN4-high OR | Fisher p |
|---|---:|---:|---:|---:|---:|---:|
| Pretreatment (primary) | 10 / 39 | 0.487 | 0.328–0.649 | 0.911 | 1.52 | 0.725 |
| On-treatment (exploratory) | 13 / 43 | 0.604 | 0.438–0.760 | 0.264 | 2.65 | 0.205 |
| Paired on − pre (exploratory) | 9 / 33 | 0.633 | 0.434–0.818 | 0.232 | 2.22 | 0.454 |

On-treatment and paired changes are post-baseline and cannot establish a
pretreatment predictor. They are not significant.

## Open melanoma B5-style odds ratios

CLDN4-high = at or above the cohort-internal median. OR is the odds of
response in CLDN4-high versus CLDN4-low, with Haldane–Anscombe 0.5
correction. OR < 1 is the claimed direction (CLDN4-high worse).

| Cohort | R / NR | AUC | MW p | OR | 95% CI | Fisher p |
|---|---:|---:|---:|---:|---:|---:|
| Riaz 2017 pre | 10 / 39 | 0.487 | 0.911 | 1.52 | 0.39–5.87 | 0.725 |
| Hugo 2016 pre | 15 / 13 | 0.538 | 0.747 | 1.31 | 0.31–5.51 | 1.000 |
| Liu 2019 pre-PD1 RNA | 41 / 61 | 0.420 | 0.160 | 0.57 | 0.25–1.27 | 0.213 |
| MH pool (those three) | 66 / 113 | — | — | **0.82** | **0.47–1.44** | — |
| Auslander 2018 pre | 2 / 5 | 0.700 | 0.561 | 7.00 | 0.22–219 | 0.429 |

Auslander is shown in tables only. Two pretreatment responders cannot test
B5.

Liu is the only cohort whose point estimate is in the claimed direction, and
it is not significant. Liu's own supplementary Table 4 already reported
CLDN4 Mann–Whitney p = 0.347 for responders versus PD; that published
number is also non-significant. The analysis here uses CR/PR versus SD/PD
on pre-PD1 RNA samples (n=102 after excluding mixed response).

## Not used, with reasons

| Cohort | Why not analyzed as a CLDN4 test |
|---|---|
| Prat 2017 GSE93157 | 25 melanoma samples, but the nCounter immune panel has no CLDN4 |
| Gide 2019 PRJEB23709 | Public deposit is raw ENA RNA; no processed matrix + response table |
| Van Allen 2015 | RNA is controlled-access |
| Chen 2016 GSE67501 | Anti-PD-1 RNA is RCC, not melanoma |

This is **not** the claimed 11-cohort meta-analysis. Missing melanoma RNA
(Gide processed, Van Allen) and all non-melanoma ICI cohorts are outside
this analog.

## Limitations

1. Riaz and Hugo are small. Wide intervals are expected.
2. GEO Riaz labels collapse CR and PR as `PRCR`.
3. Liu uses raw counts, not TPM. The rank test is valid within cohort; values
   are not comparable across cohorts.
4. Median split is the B5 recipe, not an optimized cutoff. Continuous
   Mann–Whitney is the better single-gene test and is also null.
5. No purity, subtype, or prior-ipilimumab adjustment except Liu's
   pre-PD1 biopsy filter. Riaz GEO metadata cannot drop uveal cases.
6. A three-cohort melanoma pool is not a substitute for the unpublished
   11-cohort list. It is enough to say the open melanoma RNA does not
   reproduce OR = 0.42.

## Reproduction

```bash
python3 -m pip install -r requirements.txt
python3 analyze.py
```

The script downloads public GEO, Nature, and GitHub inputs when absent.
`provenance.tsv` stores URLs and SHA-256 hashes. Large source files stay
untracked under `data/`.

Outputs:

- `cldn4_vs_response.png` / `.svg` — Riaz primary and exploratory plots
- `open_melanoma_forest.png` / `.svg` — B5-style OR forest
- `summary.csv` — Riaz tests
- `open_melanoma_summary.csv` / `open_melanoma_pooled.csv`
- `inventory.csv` — analyzed and skipped cohorts
- `sample_level.csv`, `hugo_sample_level.csv`, `liu_sample_level.csv`,
  `auslander_*.csv`, `paired_changes.csv`

## References

- Riaz N, et al. *Cell*. 2017;171:934–949. GEO GSE91061.
- Hugo W, et al. *Cell*. 2016;165:35–44. GEO GSE78220.
- Liu D, et al. *Nat Med*. 2019;25:1916–1927. Supplementary Table 1 and
  [vanallenlab/schadendorf-pd1](https://github.com/vanallenlab/schadendorf-pd1)
  `addData.zip`.
- Auslander N, et al. *Nat Med*. 2018;24:1545–1549. GEO GSE115821.
- Prat A, et al. *Cancer Res*. 2017;77:3540–3550. GEO GSE93157 (no CLDN4).
