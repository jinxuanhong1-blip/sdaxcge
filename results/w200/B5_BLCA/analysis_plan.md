# B5_BLCA — Prespecified analysis plan

**Status: PRESPECIFIED.** This file was written and committed *before* any
CLDN4-versus-response statistic was computed. Its git commit precedes the commit
that adds `summary.json` and the results tables. Nothing in this plan was revised
after seeing an outcome-associated result; if a deviation became necessary it is
recorded in the "Deviations" section of `WRITEUP.md`, not silently edited here.

## 1. The claim under test

Claim **B5** as stated by the user, verbatim:

> Reconstruct 11-cohort ICI meta CLDN4-high OR=0.42 across GBM/NSCLC/melanoma/RCC/urothelial.

This slice, **B5_BLCA**, tests the **urothelial/bladder arm only**:

> In patients with urothelial carcinoma treated with immune-checkpoint inhibitors,
> tumours with high CLDN4 expression have lower odds of objective response than
> tumours with low CLDN4 expression, with an odds ratio of approximately 0.42.

The user-asserted number is fixed as a module constant before any computation:

```python
USER_OR = 0.42
EXACT_DECIMALS = 2     # "matches" means equal after rounding to 2 dp
NEARBY_ABS    = 0.05   # "nearby" band, declared in advance
```

**Anti-anchoring commitment.** No cutoff, cohort definition, filter, normalisation
or covariate set will be chosen because it moves the estimate toward 0.42. The
primary specification below is fixed in advance; every alternative specification
is reported as a declared sensitivity, whatever it shows. `summary.json` will
carry `did_we_tune_to_0.42: false`.

**Provenance of the claim.** A prior agent's literature search for a CLDN4 × ICI
meta-analysis returned nothing (PubMed `"claudin 4" AND ("checkpoint inhibitor" OR
"anti-PD-1")` → 0 hits). We therefore treat OR = 0.42 as an **unsourced assertion**
to be tested, not as a published result to be reproduced. We will not name or
imply the existence of 11 cohorts we cannot enumerate.

## 2. Cohort inclusion criteria (fixed in advance)

A cohort is eligible if it is (a) public, (b) urothelial/bladder carcinoma,
(c) treated with an anti-PD-1 or anti-PD-L1 agent, (d) has pre-treatment tumour
RNA expression including CLDN4, and (e) has RECIST best overall response
sufficient to derive ORR.

| # | Cohort | Source | Treatment | RNA | Response field |
|---|--------|--------|-----------|-----|----------------|
| 1 | IMvigor210 (Mariathasan 2018) | `IMvigor210CoreBiologies` `cds.RData` | atezolizumab | counts → TPM | `Best Confirmed Overall Response` |
| 2 | BACI (Robertson 2021) | GEO GSE176307 | atezolizumab / pembrolizumab / other ICI | salmon TPM | `io.response` |
| 3 | Snyder 2017 | bhklab PredictIO, Zenodo 7058399 | atezolizumab | log2(TPM+0.001) | `recist` |

**No double counting.** IMvigor210 contributes **exactly one** estimate to the
meta-analysis, taken from the original `cds.RData`. The PredictIO-harmonised copy
of the same trial is used only as a processing-concordance check and is explicitly
excluded from pooling.

Cohorts considered and excluded are listed in `WRITEUP.md` with the reason.

## 3. Primary analysis (one prespecified test per cohort, one pooled estimate)

**Exposure.** CLDN4-high = expression *strictly greater than* the cohort-specific
median; ties at the median go to the low group. The split is computed within
cohort, on the cohort's own declared expression scale, and **without reference to
response**.

**Endpoint.** ORR from RECIST best overall response:
`responder = {CR, PR}`, `non-responder = {SD, PD}`.
Non-evaluable / not-assessed / missing records are **excluded** from the primary
analysis. Counts dropped at each stage are reported per cohort.

**Test.** Two-sided Fisher exact test on the 2×2 table (CLDN4 high/low ×
responder/non-responder). Effect measure: odds ratio of *response* for CLDN4-high
versus CLDN4-low, so **OR < 1 means CLDN4-high responds less often** — the
direction the claim asserts. Per-cohort OR is the conditional maximum-likelihood
estimate with an exact (Fisher) 95% CI.

**Pooling.** For inverse-variance pooling, log OR and its standard error use the
Woolf estimator with a Haldane–Anscombe 0.5 correction applied only when a cell is
zero. Primary pooled estimate: **DerSimonian–Laird random effects**, chosen in
advance because the cohorts differ in platform, agent and response ascertainment.
Because k is small, **Hartung–Knapp** confidence limits are reported alongside.
Fixed-effect inverse-variance and Mantel–Haenszel estimates are also reported.
Heterogeneity: Cochran Q with df and p, I², τ².

**Comparison to the claim.** Reported as an explicit `honest_verdict` block:
`MATCH_AT_2DP` / `NEARBY_NOT_EXACT` / `DOES_NOT_MATCH`, plus the separate and more
informative statement of whether 0.42 falls inside the pooled 95% CI. Being
"nearby" is not treated as confirmation.

## 4. Prespecified sensitivity analyses

All are reported whatever their result. None can promote or demote the primary
conclusion; they characterise robustness.

- **S1 Cutoff robustness** — upper vs lower tertile; Q4 vs Q1; upper quartile vs
  the rest; 60/40 split.
- **S2 Continuous exposure** — logistic regression of response on CLDN4, reporting
  OR per 1 unit log2 expression and per 1 within-cohort SD. Avoids dichotomisation
  entirely.
- **S3 Distribution-level** — Mann–Whitney U of CLDN4 in responders versus
  non-responders, with rank-biserial correlation.
- **S4 ITT-style endpoint** — non-evaluable patients counted as non-responders.
- **S5 Biopsy site** — IMvigor210 restricted to bladder-site biopsies.
- **S6 Normalisation** — IMvigor210 repeated on DESeq size-factor-normalised
  counts instead of TPM.
- **S7 Pipeline concordance** — CLDN4 from our own IMvigor210 extraction versus the
  PredictIO harmonisation of the same samples (Spearman correlation, and the
  primary test recomputed on the harmonised values).
- **S8 Covariate adjustment** — multivariable logistic regression within
  IMvigor210 adjusting for the prespecified set below.

Prespecified adjustment set (IMvigor210, chosen for prognostic relevance and
availability, not for effect on the estimate): `FMOne mutation burden per MB`
(log-transformed), `Immune phenotype`, `Baseline ECOG Score`, `Received platinum`,
and biopsy `Tissue` collapsed to bladder / non-bladder. Records missing an
adjustment variable are dropped from S8 only, with n reported.

## 5. Pipeline calibration controls

Prespecified, and reported even if they behave unexpectedly — they test whether
the pipeline can detect a signal at all, and whether it manufactures one.

- **Positive control:** CD8 T-effector / IFN-γ programme, mean of per-gene
  within-cohort z-scores of `CD8A, GZMA, GZMB, IFNG, EOMES, CXCL9, CXCL10, TBX21,
  PRF1`. This programme was reported to associate with atezolizumab response in
  IMvigor210 itself, so a **positive** association (OR > 1) is the expected
  direction. Failure here would indicate the pipeline is underpowered or broken.
- **Negative control:** mean within-cohort z of housekeeping genes `ACTB, GAPDH,
  TBP, RPL13A, PGK1`. Expected null (OR ≈ 1). A strong association here would
  indicate a normalisation or leakage artefact.

Controls use the identical median-split + Fisher machinery as the primary test.

## 6. Exploratory analyses (labelled exploratory, BH-adjusted)

Family: `CLDN1, CLDN2, CLDN3, CLDN7, CLDN18, TACSTD2, EPCAM, CDH1, OCLN, TJP1,
CD274` — m = 11 genes, tested by the same median-split Fisher procedure in
IMvigor210. Benjamini–Hochberg adjustment within this family, threshold q < 0.05,
reported as `q-value` (never as "FDR"). CLDN4 itself is **not** in this family; it
is the primary hypothesis and is not BH-adjusted against these genes.

Also exploratory: Spearman correlation of CLDN4 with the CD8 T-effector score and
its distribution across the IMvigor210 immune phenotype (desert / excluded /
inflamed), which speaks to the project's immune-exclusion thesis but is not a test
of the B5 claim.

## 7. Units, QC and reproducibility

- Expression scale is declared per cohort and carried in an `expression_scale`
  column of every per-sample output. IMvigor210 and BACI: `log2(TPM+1)`.
  PredictIO/Snyder: `log2(TPM+0.001)` as deposited — verified empirically, since
  the observed matrix minimum equals log2(0.001).
- CLDN4 is resolved by Ensembl `ENSG00000189143` where IDs are available and by
  exact symbol match otherwise; a symbol resolving to multiple rows is summed at
  TPM level, and the number of matched rows is recorded.
- Filtering never references CLDN4 expression or response. No low-expression
  filter is applied to CLDN4.
- One record per patient. `n` is reported per cohort at every stage: available →
  RNA present → response evaluable → analysed.
- Seed fixed in `repro/seed.txt`; environment captured in `repro/pip-freeze.txt`;
  every downloaded file recorded in `catalog.tsv` and
  `results/w200/B5_BLCA/download_manifest.json` with URL, byte count and SHA-256.

## 8. Interpretation rules agreed in advance

- The urothelial arm alone cannot confirm or refute a pooled 11-cohort number. The
  most this slice can say is whether the urothelial evidence is *compatible* with
  OR = 0.42.
- Conclusion vocabulary is restricted to **supported / suggestive / inconsistent /
  not testable**. "Supported" requires concordant direction in ≥2 independent
  cohorts and a pooled interval excluding 1.
- A null result is a result and will be reported as prominently as a positive one.
- Neither an unexpected direction nor a null will be handled by adding cohorts,
  moving cutoffs, or dropping a cohort after the fact.
