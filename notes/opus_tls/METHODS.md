# Methods

## Question

Bessede and colleagues (and related TROP2 / ADC literature) propose that
**TACSTD2-high lung tumours are TLS-low / B-cell-low / T-cell-low**, and
that this is not only a purity artefact. This slice tests that claim in
every open lung ICI or atlas RNA matrix that could be downloaded here
(< 2 GB, no EGA).

Primary genes: **TACSTD2** (TROP2), **CLDN4**.
Primary immune programmes: TLS (Cabrita 8-gene, 12-chemokine, Meylan imprint),
B-cell, plasma-cell, Tfh, CD8 T-cell, IFN-γ (Ayers).

## Expression

- TCGA LUAD/LUSC: Xena GDC STAR `log2(TPM+1)`, GENCODE v36 symbols.
  Primary tumours only (`-01A/B`); one sample per barcode.
- GEO RNA-seq: TPM or FPKM → `log2(x+1)`; GSE126044 counts → log2(CPM+1);
  GSE207422 already log2 TPM; GSE72094 already log2 microarray.
- Duplicate symbols: keep the highest-mean transcript.

## Signatures

Mean of gene-wise z-scores (sample-wise). A rank-percentile alternative is
computed as a sensitivity score (`scores_rank`) and is not the primary
statistic. Gene lists and coverage are in
`results/opus_tls/tables/signature_gene_coverage.csv`.

TLS gene sets are published lymphoid-aggregate programmes, not a
histology TLS call. No spatial / IHC TLS label is available in the open
matrices used here.

## Confounding

TACSTD2 and CLDN4 are epithelial. In bulk RNA, immune signatures fall as
tumour-cell fraction rises. Every association is therefore reported as:

1. crude Spearman ρ (bootstrap 95% CI, 1000 resamples);
2. **partial Spearman** after rank-transform + residualising on
   - **ABSOLUTE purity** (TCGA; Carter / PanCanAtlas mastercalls), and
   - **pan-epithelial score** (EPCAM/KRT8/KRT18/KRT19/CDH1/KRT7) in every
     cohort, and
   - ABSOLUTE + methylation **leukocyte fraction** in TCGA;
3. position of TACSTD2 in the **transcriptome-wide** distribution of the
   same partial correlation (is it special, or a typical epithelial gene?).

Partial Spearman: rank x, y and covariates; OLS-residualise; Pearson on
residuals; p from a t statistic with n − 2 − k df. Key TACSTD2 partial
correlations also have a 2 000-resample bootstrap 95% CI
(`key_partial_ci.csv`).

A call of **HOLDS** requires: n ≥ 40, ρ_adj < 0, p_adj < 0.05.
Small ICI cohorts (n = 16–27) are reported but labelled UNDERPOWERED.

## Other statistics

- Quartile (Q4 vs Q1) Mann–Whitney with rank-biserial r and Hodges–Lehmann
  shift; Stouffer combination across purity tertiles.
- Multivariable OLS of inverse-normal signature on inverse-normal TACSTD2
  ± purity ± age/sex/stage, HC3 SEs.
- Orthogonal TCGA readouts: CIBERSORT relative B-naive / B-memory / plasma
  / Tfh / CD8, methylation leukocyte fraction, nonsilent TMB.
- Meta-analysis of Fisher z with Spearman-specific variance
  (Bonett & Wright 2000) and DerSimonian–Laird random effects. High I² is
  expected if the effect is histology-specific; that is a result, not a
  bug.
- ICI endpoints: Cox PH for GSE135222 PFS; Mann–Whitney for DCB / RECIST /
  MPR. These are descriptive.
- Follow-up (script 06): five extra public matrices (GSE4573, GSE17710,
  GSE19188, GSE50081, GSE103584) scored the same way. GSE17710 is also
  residualised on pathologist `tumor_percent`. Formal TACSTD2 × LUSC
  interaction in pooled TCGA (inverse-normal signature ~ TACSTD2 + LUSC +
  TACSTD2×LUSC + ABSOLUTE, HC3). OS Cox in each TCGA histology.

No p-hacking of the gene sets after seeing the data: lists were locked in
`opus_tls_lib.py` before the first correlation table was written.
FDR (BH) is applied within cohort × gene across signatures.

## Software

Python 3.12, pandas, numpy, scipy, statsmodels, lifelines, matplotlib.
Scripts are under `scripts/opus_tls/`. Seed `20260816`.
