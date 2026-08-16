# Open SCLC ICI RNA: TACSTD2/CLDN4 versus outcome

## Result

In the two analyzable open cohorts, patient-level pretreatment **TACSTD2 was not
associated with RECIST objective response**. Among 38
CR/PR and 18 SD/PD patients, median log2-normalized
expression was 7.46 versus
7.5 (difference
-0.0418; rank-biserial
-0.105; two-sided Mann–Whitney
P=0.533). The cohort-adjusted odds ratio per
1-SD higher TACSTD2 was 0.826 (95% CI
0.442–1.54;
P=0.548).

**CLDN4 was not measured by the deposited GeoMx CTA/custom panel and was
skipped.** No proxy gene was substituted.

The secondary long-term-benefit comparison was also null: median TACSTD2 was
7.43 in 12
patients with TTP ≥12 months and 7.47 in
43 evaluable patients without it (Mann–Whitney
P=1). As a continuous check, Spearman
ρ for TACSTD2 versus TTP was -0.0497
(P=0.711).

## Included open cohorts

- **GSE261345 / CANTABRICO:** 26 ES-SCLC patients, pretreatment spatial RNA,
  first-line durvalumab plus platinum–etoposide.
- **GSE261348 / IMfirst:** 32 ES-SCLC patients, pretreatment spatial RNA,
  first-line atezolizumab plus platinum–etoposide.
- Expression and outcomes were taken from the pinned open harmonization
  [SCAPeSCLC v1.3.5](https://doi.org/10.5281/zenodo.21614171), derived from
  those GEO deposits. Analysis is at the patient level, avoiding ROI
  pseudoreplication.

## Interpretation limits

- Both studies are single-arm chemoimmunotherapy cohorts. This analysis is an
  **outcome association**, not evidence that TACSTD2 predicts benefit specific
  to PD-L1 blockade rather than chemotherapy or prognosis.
- Objective response (CR/PR versus SD/PD) was the primary endpoint. The
  prespecified study-style secondary endpoint was time to progression at least
  12 months; patients censored before 12 months were excluded from that binary
  endpoint.
- Tests are exploratory, two-sided, and unadjusted for clinical covariates
  beyond cohort in logistic models. `association_tests.tsv` also reports a
  conservative Holm adjustment across every inferential row shown.
- The assay is a targeted ~1,800-gene panel, not whole-transcriptome RNA-seq.
  Patient values are harmonized patient-level summaries of spatial ROIs.

## Controlled/request-only cohorts skipped

See `cohort_audit.tsv`. IMpower133 was deliberately excluded because its
expression and linked clinical data are EGA-controlled. Other cohorts were
skipped when sample-level expression linked to response was request-only,
controlled, or not deposited; plots and aggregate tables were not digitized.

## Reproduce

```bash
python3 scripts/analyze_sclc_open.py --output results/w200/SCLC_open
```

The script pins the open record version and verifies source-file MD5 checksums.
Search/audit cutoff: 2026-08-16.
