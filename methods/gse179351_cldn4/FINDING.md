# GSE179351 — CLDN4 vs immune / CD274 (additive leftover)

**Additive only.** Public leftover slice. This is **not** a lung ICI cohort. Parikh et al., *Nat Cancer* 2021 ([PMID 35122060](https://pubmed.ncbi.nlm.nih.gov/35122060/); GEO [GSE179351](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179351)). Phase II radiation + ipilimumab + nivolumab in metastatic **MSS CRC and PDAC**. Out-of-field metastatic biopsies. No prior page in this repo scored CLDN4 here.

The series matrix is public **metadata**. Expression is **not** in that matrix (`!series_matrix_table` has 0 rows). Tests use the deposited supplementary DESeq2-normalized counts. This is not an empty finding.

## Honest n

| item | n | source |
|---|---:|---|
| GEO samples | **54** | GSM5416019–GSM5416072; Public on Jul 06 2021 |
| series-matrix **expression rows** | **0** | public file; table is header then end |
| DESeq2 count columns | **54** | `GSE179351_DESeq2_NormalizedCountsForAllSamples.txt.gz` |
| titles that join GEO ↔ counts | **53** | `20-Pre-xRT` is GEO-only; `22-Pre-xRT` is count-only |
| unique patients after join | **25** | title patient ID; patient 20 does not join |
| **primary (one baseline / patient)** | **25** | Pre-Tx if present, else Pre-xRT |
| Pre-Tx only | **23** | sensitivity; drops Pre-xRT-only patients 7 and 51 |
| MSS CRC / PDAC at baseline | **14 / 11** | GEO `cancer type` on the baseline row |
| CLDN4 / CD274 / CD8A / IFNG | 25 | protein-coding rows present |
| T/NK signature | 25 | **16/16** locked genes |
| IFN (Ayers-6) | 25 | **6/6**; **IFNG present** |
| RECIST / DCR / ORR on GEO | **0** | not characteristics; not recovered here |
| PFS / OS on GEO | **0** | not characteristics |

Do not quote n=54 or n=53 as the test n. Those are GSM / biopsy counts. Several patients have 2–3 timepoints. The unit is the **patient**. Response labels used in the paper (SD/PR/CR vs PD; pretreatment 5 vs 7) are **not on GEO** (n=0).

Patient 52 is labeled MSS CRC on Pre-Tx and PDAC on later GEO rows. Baseline uses the Pre-Tx (CRC) row. `20-Pre-xRT` / `22-Pre-xRT` are **not** remapped.

## One-row table

| dataset | n | matrix | CLDN4–CD274 ρ (p) | CLDN4–T/NK ρ (p) | CLDN4–IFN ρ (p) | CLDN4–CD8A ρ (p) | response n |
|---|---:|---|---|---|---|---|---:|
| GSE179351 baseline DESeq2 | **25** | suppl. DESeq2 (series matrix empty) | **−0.060 (0.77)** | **−0.425 (0.034)** | **−0.462 (0.020)** | −0.210 (0.31) | **0** |

Full numeric row: `tables/one_row.tsv`.

## Matrix and labels (nothing invented)

| field | public? | n | what is there |
|---|---|---:|---|
| series-matrix expression | empty | 0 | RNA-seq placeholder table |
| DESeq2 normalized counts | yes | 54 columns | HTSeq → DESeq2; HG38 / Gencode |
| cancer type | yes | 53 joined | MSS CRC / PDAC |
| time point | yes | 53 | Pre-Tx 23 / Pre-xRT 19 / Post-xRT 11 |
| RECIST / DCR / ORR | **no** | 0 | not deposited |
| PFS / OS | **no** | 0 | not deposited |

Paper mapping (not a GEO field): biopsies are out-of-field metastases around ICI ± radiation. Timepoint is not a response label.

## Scores

`log2(DESeq2_normalized + 1)`. Signatures = mean of gene-wise z-scores on the **25 baseline** columns (z computed on all 54 count columns, then subset).

| signature | locked genes | used |
|---|---|---|
| T/NK | CD8A, CD8B, CD2, CD3D, CD3E, CD3G, NKG7, GNLY, PRF1, GZMA, GZMB, GZMK, KLRD1, KLRK1, KLRB1, NCR1 | **16/16** |
| IFN | Ayers-6: IFNG, STAT1, CXCL9, CXCL10, IDO1, HLA-DRA | **6/6** |
| MHC-I | HLA-A, HLA-B, HLA-C | 3/3 |

## Primary Spearman (baseline n=25)

| pair | n | ρ | p |
|---|---:|---:|---:|
| **CLDN4 vs CD274** | 25 | **−0.060** | **0.77** |
| **CLDN4 vs T/NK** | 25 | **−0.425** | **0.034** |
| **CLDN4 vs IFN (6/6)** | 25 | **−0.462** | **0.020** |
| CLDN4 vs CD8A | 25 | −0.210 | 0.31 |
| CLDN4 vs NKG7 | 25 | **−0.645** | **0.00050** |
| CLDN4 vs MHC-I | 25 | −0.358 | 0.079 |
| CLDN4 vs PDCD1 | 25 | +0.005 | 0.98 |
| CLDN4 vs TACSTD2 | 25 | +0.559 | 0.0037 |
| CLDN4 vs EPCAM | 25 | +0.779 | 4.4×10⁻⁶ |

**CLDN4 vs CD274 does not hold.** The leftover that holds is **CLDN4-high = T/NK-low / IFN-low**, strongest on NKG7. CD8A alone is the same sign and is not significant.

## After EPCAM (epithelial content)

Partial Spearman = Pearson of rank residuals on EPCAM. n=25.

| pair | partial ρ \| EPCAM | p |
|---|---:|---:|
| CLDN4 vs CD274 | −0.120 | 0.57 |
| CLDN4 vs T/NK | **−0.575** | **0.0027** |
| CLDN4 vs IFN | **−0.530** | **0.0064** |
| CLDN4 vs CD8A | **−0.443** | **0.026** |
| CLDN4 vs NKG7 | **−0.631** | **0.00073** |
| CLDN4 vs MHC-I | **−0.616** | **0.0010** |

CD274 stays null after EPCAM. T/NK and IFN stay negative. CD8A and MHC-I become significant only after the EPCAM residual — do not upgrade those to primary claims.

## Median split (CLDN4-high vs CLDN4-low)

Median CLDN4 log2(DESeq2+1) = 7.783. High n=13, low n=12. Cliff δ = P(high>low) − P(high<low). Negative δ = CLDN4-high is immune-low.

| axis | nH / nL | δ | MWU p |
|---|---|---:|---:|
| CD274 | 13 / 12 | −0.22 | 0.36 |
| T/NK | 13 / 12 | −0.42 | 0.077 |
| IFN | 13 / 12 | **−0.53** | **0.028** |
| CD8A | 13 / 12 | −0.17 | 0.49 |
| NKG7 | 13 / 12 | **−0.69** | **0.0035** |
| MHC-I | 13 / 12 | −0.36 | 0.13 |

Median-split CD274 is null. IFN and NKG7 hold; the composite T/NK split is the same sign and p=0.077 at 13 vs 12.

## Sensitivity

| cut | n | CLDN4–CD274 ρ (p) | CLDN4–T/NK ρ (p) | CLDN4–IFN ρ (p) |
|---|---:|---|---|---|
| Pre-Tx only | 23 | +0.001 (1.00) | **−0.480 (0.020)** | **−0.522 (0.011)** |
| all joined biopsies | 53 | −0.081 (0.56) | −0.445 (0.00084) | −0.335 (0.014) |
| baseline MSS CRC | 14 | −0.200 (0.49) | −0.446 (0.11) | −0.402 (0.15) |
| baseline PDAC | 11 | +0.273 (0.42) | −0.491 (0.13) | **−0.609 (0.047)** |

n=53 is **pseudoreplication** (repeat biopsies). CRC n=14 and PDAC n=11 are honest small strata, not confirmations. CD274 is null in every cut.

## How to read this

- **Matrix is public, primary n=25, not empty.** The empty objects are the series-matrix expression table and GEO response labels (n=0).
- **CLDN4 vs CD274 does not hold** (ρ=−0.06, p=0.77). Do not claim a PD-L1 axis on this leftover.
- **CLDN4 vs T/NK and IFN hold** at n=25 (ρ=−0.43 / −0.46). They remain negative after EPCAM. NKG7 is the strongest single immune gene.
- **CD8A alone does not hold** on the continuous primary test.
- **Response is empty here** (n=0 on GEO). This page does not score DCR / ORR / PFS.
- This is MSS CRC + PDAC ICI + radiation, not lung. Extra analog weight only.

## Files

- `analyze.py` — download series matrix + DESeq2 counts to `/tmp/gse179351`, then score
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `per_patient_baseline.tsv`, `label_inventory.tsv`, `gene_coverage.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd274.png`
- `figures/fig2_cldn4_vs_tnk.png`
- `figures/fig3_cldn4_vs_ifn.png`
- `figures/fig4_spearman_forest.png`

```bash
python3 methods/gse179351_cldn4/analyze.py
```
