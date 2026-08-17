# GSE248378 — CLDN4 vs T/NK or IFN (Durva neoadjuvant leftover)

**Additive only.** Prior leftover scores of this series vs recurrence / PFS (PRs 175 / 224) and the leftover-bulk “given” tag (PR 296) are **not re-audited**. This folder only asks: is a public processed matrix here, and if so what is **CLDN4 vs T/NK or IFN**.

Altorki et al., *Nat Commun* 2023 ([PMID 38114518](https://pubmed.ncbi.nlm.nih.gov/38114518/); GEO [GSE248378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378)). NCT02904954 neoadjuvant durvalumab ± SBRT, stages I–III NSCLC. Unit is the **resected FFPE tumor** (one GSM / one FPKM column). No slide was re-scored.

## Honest n

| item | n | source |
|---|---:|---|
| GEO samples | **29** | GSM7912309–GSM7912337 |
| series-matrix **expression rows** | **0** | public file; table is header then `!series_matrix_table_end` |
| supplementary FPKM samples | **29** | `GSE248378_Durva_Post_FPKMs.txt.gz` |
| supplementary FPKM genes | **15165** | HGNC symbols; Cufflinks FPKM, Gencode v19 |
| CLDN4 / TACSTD2 / CD8A / NKG7 | 29 | present on FPKM |
| T/NK signature | 29 | 16/16 locked genes present |
| IFN (Ayers-6) | 29 | **5/6**; **IFNG absent** |
| MPR / recurrence / PFS on GEO | **0** | not characteristics; not recovered here |
| Adenocarcinoma / Squamous / other | 20 / 7 / 2 | GEO `cell type` (NOS=1, sarcomatoid=1) |
| Arm1 / Arm2 | 18 / 11 | GEO `treatment` |

The series matrix is public **metadata**. Expression is **not** in that matrix (`!Sample_data_row_count` = 0). Tests use the public supplementary FPKM. Titles join 29/29 to GSM IDs. This is not empty.

Do not treat n as a pretreatment or MPR cohort. These are post-treatment resections. MPR and recurrence are not GEO fields (n=0).

## One-row table

| dataset | n | matrix | CLDN4–T/NK ρ (p) | CLDN4–IFN ρ (p) | CLDN4–CD8A ρ (p) | T/NK genes | IFN genes | IFNG | MPR n |
|---|---:|---|---|---|---|---|---|---|---:|
| GSE248378 post-durva FPKM | 29 | suppl. FPKM (series matrix empty) | **−0.578 (0.0010)** | −0.270 (0.16) | **−0.464 (0.011)** | 16/16 | 5/6 | **absent** | 0 |

Full numeric row: `tables/one_row.tsv`.

## Matrix and labels (nothing invented)

| field | public? | n | what is there |
|---|---|---:|---|
| series-matrix expression | empty | 0 | RNA-seq placeholder table |
| FPKM | yes | 29 × 15165 | Twist-exome capture RNA-seq |
| treatment | yes | 29 | Arm1 18 / Arm2 11 |
| histology | yes | 29 | Adenocarcinoma 20 / Squamous 7 / NOS 1 / sarcomatoid 1 |
| MPR | **no** | 0 | not deposited |
| recurrence / PFS / OS | **no** | 0 | not deposited |
| IFNG | **no** | 0 | not a row on the FPKM; IFNAR1/2 and IFNGR1/2 are present |

Paper mapping (not a GEO field): Arm1 = durvalumab alone; Arm2 = durvalumab + SBRT. Arm is treatment assignment, not response.

## Scores

`log2(FPKM+1)`. Signatures = mean of gene-wise z-scores on the 29 columns.

| signature | locked genes | used |
|---|---|---|
| T/NK | CD8A, CD8B, CD2, CD3D, CD3E, CD3G, NKG7, GNLY, PRF1, GZMA, GZMB, GZMK, KLRD1, KLRK1, KLRB1, NCR1 | **16/16** |
| IFN | Ayers-6: IFNG, STAT1, CXCL9, CXCL10, IDO1, HLA-DRA | **5/6** (no IFNG) |
| MHC-I | HLA-A, HLA-B, HLA-C | 3/3 |

KLRF1 was never in the locked T/NK list (it is also absent). IFN is the Ayers list minus the missing gene, not a substitute IFNG.

## Primary Spearman (n=29)

| pair | n | ρ | p |
|---|---:|---:|---:|
| **CLDN4 vs T/NK** | 29 | **−0.578** | **0.0010** |
| CLDN4 vs IFN (5/6) | 29 | −0.270 | 0.16 |
| CLDN4 vs CD8A | 29 | **−0.464** | **0.011** |
| CLDN4 vs NKG7 | 29 | **−0.555** | **0.0018** |
| CLDN4 vs MHC-I | 29 | +0.049 | 0.80 |
| CLDN4 vs TACSTD2 | 29 | +0.766 | 1.3×10⁻⁶ |
| CLDN4 vs EPCAM | 29 | +0.665 | 8.4×10⁻⁵ |
| TACSTD2 vs T/NK (companion) | 29 | −0.837 | 1.5×10⁻⁸ |
| TACSTD2 vs IFN (companion) | 29 | −0.512 | 0.0045 |
| T/NK vs IFN (sanity) | 29 | +0.787 | 4.2×10⁻⁷ |
| CD8A vs T/NK (sanity; CD8A is in T/NK) | 29 | +0.894 | 6.4×10⁻¹¹ |

CLDN4 tracks **T/NK-low** on this public FPKM. The IFN axis is the same direction and is **not** significant (and is missing IFNG). MHC-I is near zero.

## After EPCAM (epithelial content)

Partial Spearman = Pearson of rank residuals on EPCAM. n=29, df = 26.

| pair | partial ρ \| EPCAM | p |
|---|---:|---:|
| CLDN4 vs T/NK | **−0.434** | **0.019** |
| CLDN4 vs IFN | −0.133 | 0.49 |
| CLDN4 vs CD8A | −0.313 | 0.099 |
| CLDN4 vs NKG7 | −0.359 | 0.056 |
| CLDN4 vs MHC-I | +0.328 | 0.083 |

The T/NK anti-correlation is not only EPCAM. IFN goes to near zero after EPCAM.

## Median split (CLDN4-high vs CLDN4-low)

Median CLDN4 log2(FPKM+1) = 5.744. High n=15, low n=14. Cliff δ = P(high>low) − P(high<low). Negative δ = CLDN4-high is immune-low.

| axis | nH / nL | δ | MWU p |
|---|---|---:|---:|
| T/NK | 15 / 14 | **−0.65** | **0.0032** |
| NKG7 | 15 / 14 | **−0.59** | **0.0073** |
| CD8A | 15 / 14 | **−0.50** | **0.025** |
| IFN | 15 / 14 | −0.31 | 0.16 |
| MHC-I | 15 / 14 | −0.067 | 0.78 |

## Histology and arm (GEO labels only)

| stratum | n | CLDN4–T/NK ρ (p) | CLDN4–IFN ρ (p) |
|---|---:|---|---|
| Adenocarcinoma | 20 | **−0.800 (2.3×10⁻⁵)** | −0.340 (0.14) |
| Squamous | 7 | −0.857 (0.014) | −0.571 (0.18) |
| Arm1 | 18 | **−0.606 (0.0077)** | −0.393 (0.11) |
| Arm2 | 11 | −0.473 (0.14) | −0.055 (0.87) |

Squamous n=7 and Arm2 n=11 are honest small strata, not confirmations. NOS and sarcomatoid are n=1 each and are not tested alone.

## How to read this

- **Matrix is public, n=29, not empty.** The empty object is the series-matrix expression table. The FPKM is the usable matrix.
- **CLDN4 vs T/NK holds** (ρ=−0.58, p=0.001; median δ=−0.65, p=0.003). It remains negative after EPCAM (ρ=−0.43, p=0.019).
- **CLDN4 vs IFN does not hold** on the genes that exist (ρ=−0.27, p=0.16). **IFNG is not on the matrix**, so this is not a full Ayers-6 test.
- **MPR / recurrence are empty here** (n=0 on GEO). This page does not re-score those endpoints.
- TACSTD2 is a companion (vs T/NK ρ=−0.84). It is not a gate.

This is one post-treatment Twist-exome FPKM, n=29. It is extra weight for **CLDN4-high = T/NK-low** in leftover neoadjuvant durvalumab RNA. It is not an IFN or MPR claim.

## Files

- `analyze.py` — download series matrix + FPKM to `/tmp/gse248378`, then score
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `label_inventory.tsv`, `gene_coverage.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_tnk.png`
- `figures/fig2_cldn4_vs_ifn.png`
- `figures/fig3_cldn4_vs_cd8a.png`
- `figures/fig4_spearman_forest.png`

```bash
python3 methods/gse248378_cldn4/analyze.py
```
