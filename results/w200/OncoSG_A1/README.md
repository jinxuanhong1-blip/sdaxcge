# OncoSG LUAD: TACSTD2 versus immune scores after purity adjustment

## Result

The public cBioPortal/Datahub matrix was available, so the cohort was analyzed rather
than skipped. There were 169 expression samples; 169
had purity and 169 had all eight
IMSIG immune scores. Complete-case counts for each test are reported in
`associations.tsv`.

The strongest purity-adjusted association was Neutrophils
(partial Spearman rho = -0.456,
95% CI -0.569 to -0.328,
BH q = 4.05e-09; n = 169), a negative
association. Immune associations passing 5% BH FDR: B cells, Plasma cells, Monocytes, T cells, Macrophages, NK cells, Neutrophils.

TACSTD2 expression itself correlated with purity (Spearman rho =
0.238, p = 0.00184), supporting the requested
purity adjustment.

## Method

- Source: cBioPortal Datahub study `luad_oncosg_2020` (OncoSG, Nat Genet 2020).
- Exposure: TACSTD2 (Entrez 4070) z-score from log RNA-seq V2 RSEM, standardized
  relative to all samples.
- Outcomes: the eight sample-level IMSIG immune scores supplied with the study.
  Proliferation and translation were analyzed as controls but excluded from immune-test
  multiple-testing correction.
- Primary statistic: partial Spearman correlation. TACSTD2, each score, and purity were
  rank transformed; each ranked variable was residualized on ranked purity; the
  residual correlation was tested with df = n - 3.
- Multiplicity: Benjamini-Hochberg correction across the eight immune outcomes.
- Confidence intervals: Fisher-z approximation for the partial correlation.
- Missing values: per-outcome complete cases; no imputation.

This is an association analysis of bulk tumors. Purity adjustment reduces one source of
mixture confounding but does not establish cell-intrinsic regulation or causality.

## Files

- `associations.tsv`: raw and purity-adjusted results.
- `analysis_samples.tsv`: joined public values used in the tests.
- `partial_spearman_forest.png`: primary estimates and 95% intervals.
- `purity_adjusted_residuals.png`: rank-residual diagnostic panels.
- `provenance.json`: source URLs, checksums, and counts.
- `analyze.py`: complete reproduction script.

Run `python3 analyze.py` from this directory. Downloads are cached outside the
repository by default; use `--cache-dir PATH` to select another cache.
