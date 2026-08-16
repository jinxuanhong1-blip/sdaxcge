# GSE135222: TACSTD2 (TROP2) and CLDN4 vs durable clinical benefit (DCB)

**Bottom line: no significant association for either gene.** In this small
anti-PD-1/PD-L1 NSCLC cohort (n = 27; 7 DCB vs 20 NDB), neither TACSTD2 nor
CLDN4 tumor expression differed between DCB and NDB patients, and neither
discriminated DCB (all AUC 95% CIs span 0.5).

## Cohort and data

- GEO [GSE135222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222)
  (Jung et al., *Nat Commun* 2019, PMID 31537801): RNA-seq (RSEM TPM, hg19) of
  27 advanced NSCLC patients treated with anti-PD-1/PD-L1.
- Genes: TACSTD2 = ENSG00000184292, CLDN4 = ENSG00000189143; analyzed as
  log2(TPM + 1).

## DCB definition (derived, not provided by GEO)

GEO provides only the PFS event indicator and `pfs.time` (days), not an
explicit DCB label. DCB was defined as **PFS ≥ 183 days (6 months)**, the
standard convention. This is unambiguous here: every censored patient has
follow-up ≥ 205 days, so no patient's 6-month status is unknown. Result: 7 DCB
(6 censored ≥ 205 d plus 1 progressor at 257 d), 20 NDB.

## Results

| Feature | Median DCB | Median NDB | MWU p (two-sided) | BH q (2 genes) | AUC (DCB) | AUC 95% CI |
|---|---|---|---|---|---|---|
| TACSTD2 (log2 TPM+1) | 7.34 | 7.86 | 0.61 | 0.69 | 0.43 | 0.22–0.66 |
| CLDN4 (log2 TPM+1) | 7.69 | 7.13 | 0.69 | 0.69 | 0.56 | 0.31–0.80 |
| Combined z (exploratory) | 0.31 | 0.42 | 0.61 | — | 0.43 | 0.22–0.65 |

Secondary Spearman correlations with PFS time (ignoring censoring, so
interpret with caution): TACSTD2 r = −0.16 (p = 0.44), CLDN4 r = −0.14
(p = 0.48), combined z r = −0.28 (p = 0.15). Directionally the (non-significant)
trend is toward *lower* combined TACSTD2/CLDN4 in longer-PFS patients — e.g.,
the longest-PFS patient (NSCLC947, censored at 618 d) has near-floor expression
of both genes — but this is driven by a few samples and is not significant.

## Honest caveats

- **Severely underpowered**: 7 vs 20. At this size, only very large effects
  (AUC ≳ 0.85) would reach p < 0.05; a true modest effect could easily be
  missed. Equally, the null here does not "prove" no association.
- The AUC confidence intervals are wide and all include 0.5.
- DCB labels are derived from PFS, not author-assigned (though the derivation
  has no ambiguous cases).
- Bulk TPM; no adjustment for tumor purity, histology, or PD-L1 status
  (not available in GEO metadata).
- Two genes plus one exploratory composite were tested; BH correction across
  the two primary genes is reported, and nothing survives (or even approaches)
  significance before correction anyway.

## Files

- `sample_data.csv` — per-patient clinical (sex, age, PFS event/time, DCB) and
  expression values.
- `stats.csv` — group medians, Mann-Whitney U, p, BH q, AUC with bootstrap
  (2000×, stratified) 95% CI.
- `spearman_pfs.csv` — secondary Spearman correlations with PFS time.
- `boxplots_dcb.png`, `roc_dcb.png` — figures.
- Code: `scripts/gse135222_tacstd2_cldn4_dcb.py` (downloads data from GEO,
  fully reproducible; seed 20260816).
