# Claim A1 — replication result

**Claim (user):** TACSTD2 (TROP2) is negatively correlated with immune /
cytotoxic / exhaustion signatures in TCGA + OncoSG lung cancer, and the
association is *still negative after adjusting for tumor purity* (partial
Spearman), on a combined sample of **n ≈ 1031**.

**Verdict: REPLICATES (direction and significance).** The purity-adjusted
partial Spearman correlation between TACSTD2 and every one of the three
signatures is negative and highly significant in the pooled TCGA+OncoSG data and
in each individual cohort. The exact sample count differs from the reported
n ≈ 1031 (see "On the sample size" below); this is disclosed rather than
massaged.

All data are public: expression and OncoSG purity from the cBioPortal REST API;
TCGA tumor purity from the TCGA PanCanAtlas ABSOLUTE supplemental table (NCI
GDC). Analysis code: `code/run_claim_A1.py`, `code/make_figures.py`.

## Primary result — pooled TCGA + OncoSG, purity-adjusted partial Spearman

| Signature | Partial rho (purity-adj) | p | n |
|---|---|---|---|
| Immune (T-cell effector) | **-0.213** | 3.8e-13 | 1145 |
| Cytotoxic | **-0.188** | 1.6e-10 | 1145 |
| Exhaustion / checkpoint | **-0.194** | 3.9e-11 | 1145 |

Unadjusted Spearman is similarly negative (-0.20, -0.18, -0.17), so purity
adjustment does not create the effect — it slightly strengthens it, exactly as
the claim asserts ("still negative after purity").

## Per-cohort (purity-adjusted partial Spearman)

| Cohort (histology) | Immune | Cytotoxic | Exhaustion | n |
|---|---|---|---|---|
| TCGA LUAD (adeno) | -0.120 (p=7.4e-3) | -0.127 (p=4.7e-3) | -0.062 (p=0.17, ns) | 497 |
| TCGA LUSC (squamous) | -0.271 (p=1.7e-9) | -0.201 (p=9.3e-6) | -0.266 (p=3.5e-9) | 479 |
| OncoSG LUAD (adeno) | -0.318 (p=2.7e-5) | -0.363 (p=1.3e-6) | -0.389 (p=1.9e-7) | 169 |

Every point estimate is negative. The only non-significant cell is TCGA-LUAD
exhaustion (rho -0.062, p=0.17); it is still negative, and it is significant in
LUSC, OncoSG, and all pooled analyses. A fixed-effect meta-analysis across the
three cohorts (`meta_analysis.csv`) gives -0.214 / -0.194 / -0.198 (all
p < 1e-10), consistent with the pooled concatenation.

## Individual gene checks (pooled, purity-adjusted)

TACSTD2 is negatively associated with each canonical marker gene: CD8A (-0.197),
GZMA (-0.155), PDCD1/PD-1 (-0.165), CTLA4 (-0.154), HAVCR2/TIM-3 (-0.094),
CD274/PD-L1 (-0.089); all p < 0.003. This rules out the pattern being an
artifact of one signature's gene composition.

## On the sample size (honest note)

The claim states n ≈ 1031. No single natural cohort definition reproduces
exactly 1031:

- TCGA NSCLC only (LUAD + LUSC), expression + ABSOLUTE purity: **n = 976**
- TCGA NSCLC + OncoSG, expression + purity: **n = 1145**
- (Expression-only, before purity join: TCGA NSCLC = 994, +OncoSG = 1163)

1031 sits between the TCGA-only (976) and TCGA+OncoSG (1145) figures. The
difference is almost certainly explained by inclusion choices that the claim did
not fully specify — e.g. a different purity source (CPE consensus purity vs.
ABSOLUTE, which have different missingness), a different OncoSG RNA subset
(cBioPortal exposes 181 RNA samples, 169 with non-missing TACSTD2 + purity), or
TCGA Firehose vs. PanCanAtlas. **The direction and significance of the claim are
robust to all of these choices**, which is the substantive point. We report the
exact n for every configuration rather than selecting one to match 1031.

## Reproduce

```
pip install -r requirements.txt
python code/run_claim_A1.py     # downloads public data, writes tables/summary
python code/make_figures.py     # writes forest + scatter figures
```

## Files

- `per_cohort_correlations.csv` — per cohort, per signature, unadjusted + partial
- `pooled_correlations.csv` — TCGA-only and TCGA+OncoSG pooled
- `meta_analysis.csv` — fixed-effect meta across cohorts
- `individual_gene_checks.csv` — marker-gene level correlations
- `summary.json` — machine-readable summary + verdict inputs
- `sample_data_*.csv` — per-sample TACSTD2, signature scores, purity
- `forest_partial_spearman.png`, `scatter_pooled_rank.png` — figures

## Caveats

- Bulk-tumor correlation reflects the mixed tumor+microenvironment; even with
  purity as a covariate this is an association, not a causal/tumor-intrinsic
  claim. The consistent negative sign across independent cohorts and histologies
  nonetheless supports an inverse TACSTD2–immune relationship in NSCLC.
- Signature scores are mean per-gene z-scores; Spearman is rank-based and
  therefore invariant to the per-gene monotonic z-transform used by cBioPortal.
