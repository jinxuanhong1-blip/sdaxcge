# Open ICI cohorts: TACSTD2 after CLDN4

Patient is the unit. Cohorts are not pooled. No quantile search. Odds ratios are benefit in the high or per-SD group versus the reference. A zero cell is left undefined; nothing was filled in with 0.5.

OAK and POPLAR were not analyzed. The public EGA page for `EGAD00001008391` names the OAK log2(TPM+1) matrix and DAC `EGAC00001002120` (study `EGAS00001005013`; POPLAR log2(TPM+1) is `EGAD00001008390`). No EGA login was used and no matrix was downloaded, so this file contains no OAK or POPLAR odds ratio. Bessede et al., Clin Cancer Res 2024, reported a TACSTD2 result on those trials. That result was not recomputed here.

## What was scored

| Cohort | n in the response model | Endpoint | Cancer and drug |
|---|---:|---|---|
| IMvigor210 | 298 (68 CR/PR) | Best confirmed response, CR/PR vs SD/PD | Metastatic urothelial, atezolizumab, single arm |
| GSE218989 | 355 (168 responders) | Deposited responder vs non-responder | Lung, PD-1/PD-L1 inhibitor |
| GSE190265 | 43 (14 with PFS ≥ 6) | PFS time ≥ 6 in the deposited units | NSCLC, anti-PD-1 monotherapy |
| GSE135222 | 27 (7 with PFS ≥ 180 days) | PFS time ≥ 180 days | NSCLC, anti-PD-1/PD-L1 |
| GSE207422 | 24 (9 MPR) | Pathologic MPR, including pCR, vs NMPR | NSCLC, neoadjuvant anti-PD-1 plus chemotherapy |
| GSE166449 | 22 (7 responders) | Title responder vs non-responder | Advanced lung, immunotherapy |
| GSE126044 | 16 (5 responders) | GEO responder vs non-responder | NSCLC, anti-PD-1 |

GSE190266 is the Dijon anti-PD-1 series whose deposited TPM header has 16,384 fields and stops at MTMR14. TACSTD2 is not in that file. It was not scored.

IMvigor210 `binaryResponse` is CR/PR versus SD/PD. It is not durable clinical benefit. The cohort is single-arm, so the odds ratio is an on-treatment association, not an atezolizumab-versus-chemotherapy interaction. GSE135222 is a PFS landmark, not the authors' RECIST durable-benefit flag. GSE207422 is not ICI monotherapy.

Expression is log2 within each cohort (counts as log2(CPM+1); IMvigor210 as log2(count/sizeFactor+1); GSE166449 and GSE207422 were already on a log scale and were not logged again), then z-scored inside the complete-case sample. The primary contrast is a logistic odds ratio per 1 SD. The median split is a 2×2 odds ratio with a Woolf interval and a Fisher exact p value. The median for that 2×2 is taken inside the analyzed sample.

The indirect effect is a regression imputation on the probability scale for a +1 SD shift in TACSTD2: CLDN4 is linear in TACSTD2, the outcome is logistic in both, and each patient's residual is kept. Bootstrap intervals are percentile intervals from 2,000 patient resamples (seed 814), re-fitting the z-score and both regressions each time. Logistic odds-ratio intervals in the table are Wald intervals. A proportion mediated is reported only when the total-effect interval and the indirect-effect interval both exclude 0 and have the same sign.

## Benefit

Per 1 SD of TACSTD2. "After CLDN4" is the same model with z(CLDN4) added. CD8A is the positive control in the same patients, per 1 SD, unadjusted.

| Cohort | TACSTD2 OR (95% CI) | p | After CLDN4 | p | Median OR (Fisher p) | CD8A OR | p |
|---|---|---:|---|---:|---|---|---:|
| IMvigor210 | 1.06 (0.80–1.39) | 0.70 | 0.97 (0.69–1.37) | 0.88 | 1.36 (0.33) | 1.35 (1.01–1.79) | 0.042 |
| GSE218989 | 1.04 (0.84–1.28) | 0.73 | 1.22 (0.93–1.58) | 0.15 | 1.08 (0.75) | 1.40 (1.13–1.74) | 0.0023 |
| GSE190265 | 0.86 (0.46–1.62) | 0.64 | 0.91 (0.44–1.88) | 0.79 | 1.43 (0.75) | 2.75 (1.20–6.33) | 0.017 |
| GSE135222 | 1.05 (0.43–2.56) | 0.92 | 0.96 (0.26–3.58) | 0.95 | 0.61 (0.68) | 2.37 (0.80–7.03) | 0.12 |
| GSE207422 | 0.72 (0.31–1.69) | 0.46 | 0.31 (0.041–2.29) | 0.25 | 0.33 (0.40) | 2.86 (0.94–8.74) | 0.065 |
| GSE166449 | 1.43 (0.55–3.74) | 0.46 | 1.47 (0.45–4.79) | 0.52 | 3.75 (0.36) | 2.16 (0.76–6.12) | 0.15 |
| GSE126044 | 0.45 (0.14–1.46) | 0.19 | 0.80 (0.14–4.45) | 0.80 | 0.56 (1.0) | not quoted |  |

Every TACSTD2 benefit interval, before and after CLDN4, includes 1. The natural indirect effects on the probability scale also include 0. There is no TACSTD2–benefit association in these cohorts for CLDN4 to carry. The TACSTD2×CLDN4 product term was not significant in any benefit model. The smallest p was 0.27, in GSE207422.

Five NSCLC median tables match the unadjusted median TACSTD2 rows from the earlier public sweep, with the same cell counts and the same odds ratios: GSE126044, GSE135222, GSE166449, GSE190265, and GSE207422. GSE218989 was not in that sweep. On IMvigor210, a median cut taken on all 348 tumors and then applied to the 298 response-evaluable patients gives OR 1.27 (Fisher p = 0.41). That is the split from the earlier TACSTD2-only extraction. The primary row above uses the median inside the 298.

CD8A is associated with benefit in IMvigor210, GSE218989, and GSE190265, so those endpoints can detect an immune association of that size. In the 234 IMvigor210 patients with both a response and FMOne TMB, log2(TMB+1) has OR 2.57 per SD (95% CI 1.79–3.70, p = 3.6×10⁻⁷). With TACSTD2 in the model that TMB odds ratio is 2.72 (p = 1.9×10⁻⁷). TACSTD2 given TMB is 0.79 (0.57–1.08, p = 0.14). TACSTD2 given TMB and CLDN4 is 0.77 (0.53–1.12, p = 0.17).

GSE126044's CD8A per-SD logistic fit separated and is not quoted. In that same cohort, 154 of 2,000 benefit bootstraps did not converge, so the median odds ratio is the stable summary. Overall survival on IMvigor210 (348 patients, 232 deaths) is null: HR per SD of TACSTD2 1.00 (0.88–1.12, p = 0.96), and 1.01 (0.87–1.18, p = 0.89) after CLDN4. The bootstrap interval for the change in the log hazard ratio is −0.11 to 0.10.

## TMB and immune scores

TMB is present as a per-sample column only in IMvigor210 (`FMOne mutation burden per MB`, 272 tumors with TACSTD2 and CLDN4). The other GEO series scored here have no TMB field. FMOne TMB and neoantigen burden correlate at Spearman ρ = 0.69 (n = 218, p = 7.3×10⁻³²), which is a check that the column is the mutation-burden field, not a second endpoint.

Standardized beta is the change in the z-scored outcome per 1 SD of TACSTD2. Intervals are bootstrap percentile intervals.

| Cohort | Outcome | n | Spearman ρ (p) | Partial ρ given CLDN4 (p) | Beta | After CLDN4 |
|---|---|---:|---|---|---|---|
| IMvigor210 | log2(TMB+1) | 272 | +0.234 (9.8×10⁻⁵) | +0.117 (0.053) | +0.227 (0.118–0.337) | +0.144 (0.013–0.274) |
| IMvigor210 | CD8A | 348 | −0.226 (2.2×10⁻⁵) | −0.154 (0.0041) | −0.203 (−0.302 to −0.105) | −0.141 (−0.272 to −0.023) |
| IMvigor210 | Teff | 348 | −0.198 (2.1×10⁻⁴) | −0.133 (0.013) | −0.187 (−0.284 to −0.092) | −0.127 (−0.261 to −0.008) |
| IMvigor210 | Phenotype order | 284 | −0.132 (0.026) | −0.089 (0.14) | −0.123 (−0.233 to −0.011) | −0.069 (−0.197 to 0.063) |
| GSE218989 | CD8A | 355 | −0.262 (5.3×10⁻⁷) | −0.202 (1.3×10⁻⁴) | −0.247 (−0.356 to −0.133) | −0.211 (−0.344 to −0.068) |

Teff is the mean of within-cohort z-scores of CD8A, GZMA, GZMB, IFNG, EOMES, CXCL9, CXCL10, and TBX21. Phenotype order is desert = 0, excluded = 1, inflamed = 2 (76 / 134 / 74).

The TMB association shrinks after CLDN4, and the bootstrap interval for that drop just excludes 0 (0.001 to 0.168 on the standardized scale). The adjusted coefficient is still above 0. The indirect effect is 0.37 of the unadjusted coefficient. The rank correlation after CLDN4 is ρ = 0.117, p = 0.053. A median split of TMB, which is this cohort's median and not 10 mut/Mb, gives a TACSTD2 odds ratio of 1.52 per SD (1.16–2.00, p = 0.0022) and 1.30 (0.95–1.78, p = 0.095) after CLDN4. The interval for that change in the log odds ratio includes 0. The unadjusted median-split odds ratio is 1.98 (1.22–3.21, Fisher p = 0.0075); after CLDN4 it is 1.46 (0.85–2.54, p = 0.17).

The CD8A association does not come out. In IMvigor210 and in GSE218989 the partial correlation and the adjusted coefficient stay negative, and the interval for the drop includes 0. The same pattern holds for the IMvigor210 Teff score. Binary inflamed versus desert or excluded is null before CLDN4 (OR 0.88, 0.68–1.13, p = 0.31) and after it (OR 1.04, 0.76–1.43, p = 0.82).

The smaller NSCLC series do not have a TACSTD2–CD8A or TACSTD2–Teff association whose interval excludes 0, except GSE207422 CD8A (ρ = −0.50, p = 0.013, n = 24). There TACSTD2 and CLDN4 correlate at Pearson r = 0.89. The adjusted coefficient moves further from 0, and the interval for that move includes 0. That is collinearity in a two-predictor model with 24 patients, not evidence that CLDN4 strengthens the CD8A link.

Pearson r of TACSTD2 with CLDN4 in the benefit samples is 0.61 (IMvigor210), 0.59 (GSE218989), 0.48 (GSE190265), 0.74 (GSE135222), 0.89 (GSE207422), 0.58 (GSE166449), and 0.72 (GSE126044).

## What this does and does not say

In these open cohorts the TACSTD2 association with benefit is null, so conditioning on CLDN4 has nothing to attenuate and the indirect effect includes 0. Where TACSTD2 is inversely associated with CD8A, that association is still present after CLDN4. Where TACSTD2 is positively associated with TMB, in IMvigor210 only, part of the linear association is shared with CLDN4 and a residual association remains.

This is an observational regression. CLDN4 may be a mediator, a confounder, or both. The imputation assumes a linear CLDN4 model, no unmeasured confounding of the TACSTD2–CLDN4, CLDN4–outcome, and TACSTD2–outcome paths, and no interaction. It does not identify a causal path, and it does not re-test the controlled OAK/POPLAR analysis.

Figure: `figures/fig_response_or.png`, `figures/fig_tmb_immune.png`. Rows: `tables/response_or.tsv`, `tables/continuous.tsv`, `tables/per_sample.tsv`.

```bash
python3 methods/ici_tacstd2_cldn4_mediation/analyze.py
```

The script downloads GEO supplements and the IMvigor210 `cds.RData` mirror into `data/cache/` (gitignored) and refits every model. It needs Python with numpy, scipy, pandas, statsmodels, and matplotlib, plus R with the survival package.
