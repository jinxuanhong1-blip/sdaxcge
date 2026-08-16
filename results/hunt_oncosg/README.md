# OncoSG LUAD: purity-partial TACSTD2 versus immune programs

## What was found

A public OncoSG LUAD expression cohort exists in cBioPortal/Datahub as
`luad_oncosg_2020` (Chen et al., *Nature Genetics* 2020, PMID
[32015526](https://pubmed.ncbi.nlm.nih.gov/32015526/)). A targeted search found
no GEO or Synapse accession for this cohort. The former OncoSG URL
`https://src.gisapps.org/OncoSG/` returned HTTP 404 on 2026-08-16.

The public Datahub file is **not raw expression**. It contains gene-wise
z-scores calculated from log RNA Seq V2 RSEM expression relative to all
samples. TACSTD2 has no missing values in the downloadable matrix. Because a
gene-wise z-score is a monotone linear transformation, TACSTD2 ranks—and hence
its Spearman correlations—are unchanged from the underlying log-expression
ranks, absent additional rounding ties.

## Cohort audit

- cBioPortal currently reports 305 total study samples and 181 samples with
  RNA Seq V2 data.
- The downloadable all-sample expression matrix actually contains 169 sample
  columns.
- All 169 matrix samples match the clinical sample file and have nonmissing
  TACSTD2, source-provided `PURITY`, and all ten IMSIG fields.
- The count mismatch is retained as a limitation; results use the 169 observed
  complete cases and do not claim 181.
- Checksums and exact source URLs are in `source_manifest.tsv`; raw source
  files remain outside git under `/tmp/hunt_oncosg_data` by default.

## Analysis and result

For each score, TACSTD2 and the score were rank-transformed, each was linearly
residualized on ranked purity (with an intercept), and the residuals were
Pearson-correlated. This is a one-covariate partial Spearman correlation.
Two-sided t-test P values use `n - 3` degrees of freedom. BH correction was
applied only across the prespecified family of eight immune metrics. IMSIG
proliferation and translation are reported as nonimmune controls and are not
included in that FDR family.

Seven of eight immune associations remain inverse at BH q < 0.05. The largest
absolute adjusted association is neutrophils (partial rho = -0.456,
95% CI -0.568 to -0.328, P = 5.06e-10, BH q = 4.05e-9). Interferon is the
exception (partial rho = -0.063, P = 0.416, q = 0.416). See
`tacstd2_immune_partial.csv` for all estimates.

## Honest interpretation

These are observational bulk-tumor associations, not evidence that TACSTD2
causes immune exclusion. IMSIG scores and TACSTD2 come from the same expression
matrix, so compositional effects and signature overlap can contribute to the
correlations. Purity is a source-provided numeric field, but its estimation
method is not encoded in the cBioPortal clinical metadata used here; this
analysis therefore does not label it as a specific algorithm's estimate.
Residualizing on purity reduces linear rank confounding by that estimate but
does not deconvolve malignant-cell expression or remove all technical and
clinical confounding.

## Files

- `tacstd2_immune_partial.csv`: primary and control association estimates.
- `analysis_samples.csv`: exact 169-row analysis input extracted from the
  downloaded public files.
- `data_audit.json`: observed counts, profile identifier, checksums, and
  software versions.
- `source_manifest.tsv`: download URLs, byte sizes, and SHA-256 checksums.
- `summary.json`: compact machine-readable result.

## Reproduce

```bash
python3 -m pip install -r scripts/hunt_oncosg/requirements.txt
python3 scripts/hunt_oncosg/download.py
python3 scripts/hunt_oncosg/analyze.py
```

Set `ONCOSG_DATA_DIR` or `ONCOSG_RESULTS_DIR` to override the default data and
output locations.
