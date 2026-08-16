# GSE166449: TACSTD2 and CLDN4 vs pembrolizumab response

## Bottom line

This small cohort provides **no persuasive evidence** that pretreatment `TACSTD2`,
`CLDN4`, or their equal-weight two-gene score is associated with pembrolizumab
response. Point estimates are weak and uncertain. If anything, `TACSTD2` is
numerically higher in responders, opposite to a simple “high TROP2 means
resistance” hypothesis, but the confidence interval includes effects in both
directions.

## Cohort and endpoint

- 22 advanced lung adenocarcinoma patients: 7 responders and 15 non-responders.
- All received pembrolizumab; all biopsies were pretreatment.
- The source paper assessed response using RECIST 1.1 and defined CR/PR as
  responder and SD/PD as non-responder.
- GEO titles provide only the binary responder label, not patient-level RECIST
  categories, follow-up, survival, treatment line, PD-L1, histology details, or
  clinical covariates.
- The deposited matrix values are `log2(TPM + 1)` (although the source filename
  says TPM). Tests use these deposited log-scale values.

## Results

| Feature | Mean R | Mean NR | Difference R−NR (95% CI) | Exact Wilcoxon/Mann–Whitney p | BH q (2 genes) | AUC, high predicts R (bootstrap 95% CI) |
|---|---:|---:|---:|---:|---:|---:|
| TACSTD2 | 3.295 | 2.775 | 0.520 (-1.204, 2.244) | 0.407 | 0.814 | 0.619 (0.333, 0.867) |
| CLDN4 | 1.539 | 1.402 | 0.136 (-0.733, 1.006) | 0.945 | 0.945 | 0.514 (0.238, 0.781) |

The equal-weight score (mean of within-cohort z-scored `TACSTD2` and `CLDN4`)
also does not separate groups: AUC 0.610
(bootstrap 95% CI 0.333–0.857);
exact label-permutation p=0.546. This score is
exploratory and was not externally validated.

`TACSTD2` and `CLDN4` are moderately correlated in this cohort (Spearman
rho=0.475, p=0.026), so treating them as independent
signals would overstate the information available.

## Interpretation

- `TACSTD2`: exact two-sided p=0.407;
  AUC=0.619. This is compatible with
  chance and is not evidence of predictive utility.
- `CLDN4`: exact two-sided p=0.945;
  AUC=0.514, essentially random.
- Neither gene survives even a minimal two-gene multiplicity correction.
- These are treatment-outcome associations in a single-arm cohort. They cannot
  establish a treatment-specific predictive biomarker because there is no
  untreated/control arm; prognostic effects cannot be separated from
  pembrolizumab interaction effects.
- With only 7 responders, estimates and bootstrap intervals are wide. Cutpoint
  searching, multivariable fitting, or reporting an in-sample optimized model
  would be especially prone to overfitting, so none was used.
- Bulk RNA-seq mixes tumor expression with purity and cell-composition effects.
  No available covariates permit adjustment for those factors or FFPE versus
  fresh tissue.

## Files

- `sample_level_expression.csv`: response labels and both genes per sample.
- `association_statistics.csv`: descriptive statistics, effect sizes, tests,
  AUCs, and univariable odds ratios.
- `summary.json`: machine-readable results and explicit honest verdict.
- `expression_by_response.*` and `tacstd2_cldn4_scatter.*`: plots.
- `analysis_manifest.txt`: input checksums and software versions.
- `analyze.py`: complete reproducible analysis.

## Sources

1. GEO [GSE166449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166449).
2. Lee et al. *Cell* (2021), PMID 33857424,
   DOI [10.1016/j.cell.2021.03.030](https://doi.org/10.1016/j.cell.2021.03.030).
