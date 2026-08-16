# Kim 2018 gastric ICI: CLDN4 versus pembrolizumab response

## Bottom line

Public pretreatment tumor **CLDN4 does not associate with pembrolizumab
response** in the open Kim metastatic gastric cohort. Median-split high versus
low CLDN4: **OR 1.17, n=45, Fisher p=1.0** (95% CI 0.32–4.25). This is a null
result, not evidence that CLDN4 predicts gastric ICI benefit.

The same public labels and matrix recover the expected immune-gene signal
(CD274 / CD8A / CXCL9), so the CLDN4 null is not a broken-file artifact.

## Cohort and endpoint

- Source: Kim et al., *Nat Med* 2018;24:1449–1458
  ([doi:10.1038/s41591-018-0101-z](https://doi.org/10.1038/s41591-018-0101-z);
  PMID 30013197; NCT02589496).
- Treatment: pembrolizumab (anti-PD-1) salvage therapy in metastatic gastric
  cancer.
- Public RNA + sample-level RECIST: ORCESTRA / PredictIO `ICB_Kim` TSV release,
  Zenodo [10.5281/zenodo.7058399](https://doi.org/10.5281/zenodo.7058399)
  (`ICB_Kim.zip`). Raw FASTQ is ENA [PRJEB25780](https://www.ebi.ac.uk/ena/browser/view/PRJEB25780)
  and was not reprocessed.
- n=45 pretreatment tumors with gene-level log2 TPM (Gencode v19/v40-style
  Ensembl IDs in the release). Floor near log2(0.001) is already applied.
- Gene: **CLDN4**, Ensembl **ENSG00000189143.10**.
- Primary binary endpoint: PredictIO `recist` **CR versus SD+PD**. In this
  release every objective responder is coded `CR` (13 CR, 14 SD, 18 PD); there
  are no `PR` labels. PredictIO `response` is R for CR, NR for PD, and missing
  for SD.
- Primary exposure: cohort-median split of CLDN4 (high ≥ median).
- Primary effect: two-sided Fisher exact **odds ratio, n, and p**. Woolf logit
  95% CI is reported beside the Fisher OR. Continuous Mann–Whitney AUC is
  supporting only.

## Primary result (OR / n / p)

| Gene | Contrast | High R/NR | Low R/NR | OR | 95% CI | n | Fisher p |
|---|---|---:|---:|---:|---|---:|---:|
| **CLDN4 (primary)** | CR vs SD+PD | 7/16 | 6/16 | **1.17** | 0.32–4.25 | **45** | **1.0** |
| CLDN4 | CR vs PD (drop SD) | 7/9 | 6/9 | 1.17 | 0.28–4.87 | 31 | 1.0 |
| CLDN4 | CR+SD vs PD | 13/10 | 14/8 | 0.74 | 0.22–2.46 | 45 | 0.76 |
| TACSTD2 | CR vs SD+PD | 7/16 | 6/16 | 1.17 | 0.32–4.25 | 45 | 1.0 |
| EPCAM | CR vs SD+PD | 7/16 | 6/16 | 1.17 | 0.32–4.25 | 45 | 1.0 |
| CD274 (control) | CR vs SD+PD | 10/13 | 3/19 | 4.87 | 1.12–21.2 | 45 | 0.047 |

CLDN4 high is not enriched for response. The 2×2 is essentially balanced
(7 versus 6 responders on each side of the median). The confidence interval
includes both a 3-fold protective effect and a 4-fold harmful effect.

## Supporting continuous tests

| Gene | Contrast | R / NR | Median R | Median NR | AUC (R higher) | 95% bootstrap CI | MW p |
|---|---|---:|---:|---:|---:|---|---:|
| CLDN4 | CR vs SD+PD | 13 / 32 | 7.33 | 6.54 | 0.618 | 0.434–0.790 | 0.225 |
| CLDN4 | CR vs PD | 13 / 18 | 7.33 | 6.94 | 0.611 | 0.403–0.803 | 0.307 |
| TACSTD2 | CR vs SD+PD | 13 / 32 | 6.06 | 5.85 | 0.526 | 0.332–0.716 | 0.793 |
| EPCAM | CR vs SD+PD | 13 / 32 | 8.45 | 8.99 | 0.519 | 0.322–0.711 | 0.851 |
| CD274 | CR vs SD+PD | 13 / 32 | 2.86 | 0.99 | 0.781 | 0.619–0.918 | 0.0035 |

Medians are on the release log2 TPM scale. Bootstrap CIs use 20,000
patient-level resamples (seed 25780). Direction for CLDN4 is slightly
responder-higher and is not significant. TACSTD2 and EPCAM are also null.
CD274 (PD-L1 RNA) is higher in responders, matching the original paper's
PD-L1 story at the RNA level.

## Important limitations

1. **n=45 is small and the OR is uninformative.** Thirteen responders cannot
   support a biomarker claim. A null here is the expected outcome for a modest
   effect, not proof of no association.
2. **PredictIO collapsed CR and PR.** The public `recist` field has 13 `CR` and
   0 `PR`. The original trial reported mixed CR/PR among a 61-patient cohort;
   this RNA slice cannot reconstruct response depth. We treat public `CR` as
   objective response because that is what the open table contains.
3. **MSI-H and EBV drive response in this trial and are not in the public
   TSV.** Kim et al. reported ORR 85.7% in MSI-H and 100% in EBV-positive
   tumors. The Zenodo metadata columns for TMB, CIN, and histology are empty
   for `ICB_Kim`. An unadjusted CLDN4 test cannot separate a tight-junction
   effect from those molecular subtypes.
4. **CLDN4 is an epithelial / tight-junction gene.** High CLDN4 may index
   epithelial content rather than immune exclusion. EPCAM is likewise null;
   immune genes are not. That pattern argues against a hidden global
   normalization failure, not for a CLDN4 mechanism.
5. **Single open gastric ICI RNA cohort.** This is the public gastric
   pembrolizumab RNA set used by PredictIO (`predictio.ca` dataset 22). There
   is no second public gastric ICI RNA matrix in this slice to test sign
   concordance.
6. **Processed matrix only.** Results are tied to the ORCESTRA kallisto / gene
   TPM release. They are not interchangeable with TIDE quantile log10 FPKM or
   a local STAR recount of PRJEB25780.

## Reproduction

From this directory:

```bash
python3 -m pip install -r requirements.txt
python3 analyze.py
```

The script downloads `ICB_Kim.zip` when absent, checks the published MD5,
extracts the three needed TSVs, verifies 45 samples and the CR/SD/PD counts,
and regenerates tables and figures. `provenance.tsv` records checksums.

Outputs:

- `or_n_p.csv`: Fisher OR, Woolf 95% CI, 2×2 counts, n, p
- `summary.csv`: OR tests plus Mann–Whitney / AUC
- `sample_level.csv`: public RECIST and gene values
- `cldn4_vs_response.png` / `.svg`
- `provenance.tsv`

## References

- Kim ST, Cristescu R, Bass AJ, et al. Comprehensive molecular
  characterization of clinical responses to PD-1 inhibition in metastatic
  gastric cancer. *Nat Med*. 2018;24:1449–1458.
  [doi:10.1038/s41591-018-0101-z](https://doi.org/10.1038/s41591-018-0101-z)
- ORCESTRA ICB TSV files, Zenodo
  [10.5281/zenodo.7058399](https://doi.org/10.5281/zenodo.7058399)
- ENA project [PRJEB25780](https://www.ebi.ac.uk/ena/browser/view/PRJEB25780)
