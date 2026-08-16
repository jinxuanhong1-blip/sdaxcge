# Riaz 2017 melanoma ICI: CLDN4 versus response

## Bottom line

Pretreatment tumor **CLDN4 does not separate nivolumab responders from
nonresponders in this public cohort**. Among 49 evaluable pretreatment samples,
the responder-higher AUC was 0.487 (bootstrap 95% CI 0.328–0.649) and the
two-sided Mann–Whitney p-value was 0.911. This is a null result, not evidence
that CLDN4 predicts response.

## Cohort and endpoint

- Source: Riaz et al., *Cell* 2017, GEO
  [GSE91061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061).
- Treatment: nivolumab (anti-PD-1) in advanced melanoma.
- Public RNA data: 109 biopsies from 65 patients (51 pretreatment, 58
  on-treatment), quantified by GEO as FPKM against hg19 knownGene identifiers.
- Gene: **CLDN4**, NCBI Entrez Gene ID **1364**.
- Primary analysis: one pretreatment biopsy per patient with a known response.
- Response definition: GEO `PRCR` = responder; `SD` or `PD` = nonresponder.
  `UNK` was excluded.
- Expression scale: `log2(FPKM + 1)`. The Mann–Whitney result is unchanged by
  this monotonic transformation.

## Results

| Analysis | R / NR | Median R | Median NR | AUC (R higher) | 95% bootstrap CI | MW p |
|---|---:|---:|---:|---:|---:|---:|
| Pretreatment (primary) | 10 / 39 | 0.246 | 0.205 | 0.487 | 0.328–0.649 | 0.911 |
| On-treatment (exploratory) | 13 / 43 | 0.278 | 0.134 | 0.604 | 0.438–0.760 | 0.264 |
| Paired on − pre change (exploratory) | 9 / 33 | 0.212 | 0.015 | 0.633 | 0.434–0.818 | 0.232 |

Medians are on the log2(FPKM + 1) scale; the paired row reports median changes.
For the primary analysis, CLDN4 was nonzero in 10/10 responders and 35/39
nonresponders (two-sided Fisher p=0.569). Neither expression magnitude nor
detection supports association with response.

The on-treatment and paired analyses are post-baseline and therefore cannot
establish pretreatment predictive value. Their confidence intervals are wide
and their p-values are not significant; they should not be presented as
positive findings.

## Important limitations

1. The effective primary sample is small and imbalanced (10 versus 39), so the
   result is imprecise despite being centered near no discrimination.
2. GEO combines complete and partial responses as `PRCR`; response depth cannot
   be reconstructed from these public annotations.
3. The public sample metadata do not expose prior ipilimumab status, melanoma
   subtype, biopsy site, tumor purity, or other patient covariates needed for a
   credible adjusted model. In particular, ocular/uveal samples cannot be
   excluded from GEO metadata alone.
4. FPKM is suitable for this within-gene, between-sample comparison but is not a
   raw-count differential-expression analysis. CLDN4 is low and right-skewed,
   which is why the prespecified rank test is emphasized.
5. This is an association analysis of one marker. It does not test whether
   CLDN4 is specifically predictive of nivolumab benefit versus prognostic in
   untreated melanoma.

## Reproduction

From this directory:

```bash
python3 -m pip install -r requirements.txt
python3 analyze.py
```

The script downloads the two public GEO inputs when absent, checks dimensions
and sample-name concordance, and regenerates all tables and figures.
`provenance.tsv` records input URLs and SHA-256 hashes.

Outputs:

- `cldn4_vs_response.png` / `.svg`: primary and exploratory plots
- `summary.csv`: tests, effect sizes, confidence intervals, and detection rates
- `sample_level.csv`: parsed public annotations and CLDN4 values
- `paired_changes.csv`: matched-patient changes
- `provenance.tsv`: source URLs and checksums

## References

- Riaz N, Havel JJ, Makarov V, et al. Tumor and Microenvironment Evolution
  during Immunotherapy with Nivolumab. *Cell*. 2017;171:934–949.e16.
  [doi:10.1016/j.cell.2017.09.028](https://doi.org/10.1016/j.cell.2017.09.028)
- NCBI GEO accession
  [GSE91061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061)
