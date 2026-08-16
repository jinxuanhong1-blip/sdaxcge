# Methods — B1 analog, TCGA-PRAD

## Question

Is `CLDN4` the top co-expression partner of `TACSTD2` (TROP2) among
cell-surface genes in TCGA prostate adenocarcinoma?

## Why this design

This is the prostate analog of `results/w200/B1_BRCA`. The B1 readout is
an *honest surface rank*: every surfaceome gene is scored against the
anchor; the focus gene's true position is reported. We do not restrict
the neighbourhood to claudins or to a hand-picked panel.

`w200` is only the reporting window (`top200.csv`). The ranking always
covers the full surface-gene universe.

## Data

| Input | Source | Notes |
| --- | --- | --- |
| Expression | UCSC Xena `TCGA.PRAD.sampleMap/HiSeqV2` | log2(norm_count+1), HGNC symbol rows |
| Surfaceome | Bausch-Fluck et al., PNAS 2018, table S3 | "in silico surfaceome only" sheet; UniProt gene symbols |
| Surfaceome file | [steveneschrich/surfaceome](https://github.com/steveneschrich/surfaceome) `data-raw/surfy/table_S3_surfaceome.xlsx` | Same workbook as the BRCA analog. The Wollscheid lab host now returns HTML / Git LFS pointers. |

Primary tumours only (TCGA sample-type code `01`). This matrix has 497
code-01 aliquots from 497 patients, plus 52 adjacent-normal (`11`) and
1 metastatic (`06`) aliquot that are not used.

## Statistics

- Primary metric: Spearman ρ of each surface gene vs TACSTD2.
- Secondary: Pearson r (same samples).
- Two-sided p from the correlation t approximation; Benjamini–Hochberg FDR
  across the ranked universe.
- Genes with undefined correlation (zero variance) are dropped, not ranked.
- Stability: 1,000 case-resampling bootstraps (seed 20260816) of the
  **full** surfaceome ranking. Report the 95% percentile CI on CLDN4's
  Spearman ρ and the fraction of resamples in which CLDN4 is rank #1.

## What this is not

- Not protein, IHC, or single-cell.
- Not a purity-adjusted partial correlation (ABSOLUTE / ESTIMATE were
  not applied; bulk prostate tumours are typically high-purity, but
  stroma still mixes the signal).
- Not a claim that CLDN4 is the unique partner: a 0.001 Spearman gap
  is a near-tie.

## Reproduce

```bash
pip install -r requirements.txt
python scripts/download_data.py
python scripts/coexpression.py
```
