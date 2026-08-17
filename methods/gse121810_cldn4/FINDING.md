# GSE121810 — CLDN4 vs immune / CD274 (GBM pembrolizumab leftover)

**Additive only.** Prior B5 analog (PR 202) already scored **CLDN4 vs OS/PFS** on this series. That survival audit is **not re-run**. This folder only asks: is a public processed matrix here, and if so what is **CLDN4 vs immune / CD274**.

Cloughesy et al., *Nat Med* 2019 ([PMID 30742122](https://pubmed.ncbi.nlm.nih.gov/30742122/); GEO [GSE121810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121810)). NCT02852655 recurrent, surgically resectable glioblastoma; pembrolizumab neoadjuvant (one dose before resection) vs adjuvant-only. Unit is the **resected tumor** (one GSM / one HUGO column). No slide was re-scored.

## Honest n

| item | n | source |
|---|---:|---|
| GEO samples | **29** | GSM3447008–GSM3447036 |
| series-matrix **expression rows** | **0** | public file; table is header then `!series_matrix_table_end` |
| supplementary HUGO samples | **29** | `GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx` |
| supplementary HUGO genes | **37036** | HUGO symbols; htseq-count, GRCh38 |
| CLDN4 / CD274 / CD8A / PTPRC rows | 29 | present on HUGO matrix |
| **CLDN4 CPM ≥ 1** | **3** | median CPM **0.481** (IQR 0.303–0.640); max 8.03 |
| T/NK signature | 29 | **16/16** locked genes present |
| IFN (Ayers-6) | 29 | **6/6** rows present; **IFNG usable n = 0** (0/29 CPM ≥ 1; max 0.21) |
| RECIST / ORR on GEO | **0** | not a characteristic |
| OS / PFS in this leftover | **0** | not scored here |
| Neoadjuvant (A) / adjuvant (B) | 14 / 15 | GEO `therapy` |

The series matrix is public **metadata**. Expression is **not** in that matrix (`!Sample_data_row_count` is an empty table). Tests use the public supplementary HUGO counts. Patient numbers join 29/29 to GSM titles (xlsx `Pt3_A` = GEO “Patient 3 tumor”). This is **not** an empty accession.

It is also **not** a CLDN4-high epithelial cohort. GBM is glial. CLDN4 sits at the detection floor (28/29 nonzero raw counts, but only **3/29** reach CPM ≥ 1). GFAP median is 11,344 CPM in the same matrix. Spearman on n=29 is a rank test of mostly-floor CLDN4. The CPM ≥ 1 subset is **n=3** and is not tested (Spearman requires n≥5).

Do not treat n=29 as a pretreatment lung ICI cohort. Arm A RNA is **on-treatment** (after one pembrolizumab dose). Arm B is pretreatment. RECIST is empty here (n=0).

## One-row table

| dataset | n | usable CLDN4 (CPM≥1) | matrix | CLDN4–CD274 ρ (p) | CLDN4–T/NK ρ (p) | CLDN4–IFN ρ (p) | CLDN4–CD8A ρ (p) | T/NK genes | IFN genes | IFNG | RECIST n |
|---|---:|---:|---|---|---|---|---|---|---|---|---:|
| GSE121810 recurrent GBM HUGO | 29 | **3** | suppl. counts (series matrix empty) | +0.294 (0.12) | −0.115 (0.55) | +0.186 (0.33) | −0.060 (0.76) | 16/16 | 6/6 | **floor (0/29 CPM≥1)** | 0 |

Full numeric row: `tables/one_row.tsv`.

## Matrix and labels (nothing invented)

| field | public? | n | what is there |
|---|---|---:|---|
| series-matrix expression | empty | 0 | RNA-seq placeholder table |
| HUGO counts | yes | 29 × 37036 | htseq-count; log2(CPM+1) for tests |
| therapy | yes | 29 | neoadjuvant 14 / adjuvant 15 |
| RECIST / ORR | **no** | 0 | not deposited |
| OS / PFS | not used | 0 | prior B5 analog only |
| CLDN4 | yes, floor | 29 row / **3** CPM≥1 | Pt33_A 3.04; Pt15_B 8.03; Pt20_B 1.32 |
| CD274 | yes | 29 (27 CPM≥1) | median CPM 2.49 |
| IFNG | row only | **0** usable | 6/29 nonzero; max CPM 0.21 |
| EPCAM | yes, floor | 5 CPM≥1 | not an epithelial-content control here |

Paper mapping (GEO field is `therapy`, not “arm”): neoadjuvant = one dose before surgery (on-treatment RNA); adjuvant = post-surgical PD-1 only (pretreatment RNA). Therapy is assignment, not response.

## Scores

`log2(CPM+1)` from library-size CPM. Signatures = mean of gene-wise z-scores on the 29 columns.

| signature | locked genes | used |
|---|---|---|
| T/NK | CD8A, CD8B, CD2, CD3D, CD3E, CD3G, NKG7, GNLY, PRF1, GZMA, GZMB, GZMK, KLRD1, KLRK1, KLRB1, NCR1 | **16/16** |
| IFN | Ayers-6: IFNG, STAT1, CXCL9, CXCL10, IDO1, HLA-DRA | **6/6** (IFNG is a row, not a usable gene) |
| MHC-I | HLA-A, HLA-B, HLA-C | 3/3 |

KLRF1 is on the matrix and was never in the locked T/NK list. IFN includes the IFNG row (honest 6/6); that gene is empty and does not drive the score (STAT1 / HLA-DRA / CXCL9 / CXCL10 are the expressed members).

## Primary Spearman (n=29)

| pair | n | ρ | p |
|---|---:|---:|---:|
| **CLDN4 vs CD274** | 29 | +0.294 | 0.12 |
| **CLDN4 vs T/NK** | 29 | −0.115 | 0.55 |
| CLDN4 vs IFN (6/6; IFNG floor) | 29 | +0.186 | 0.33 |
| CLDN4 vs CD8A | 29 | −0.060 | 0.76 |
| CLDN4 vs NKG7 | 29 | −0.276 | 0.15 |
| CLDN4 vs PTPRC (CD45) | 29 | +0.019 | 0.92 |
| CLDN4 vs MHC-I | 29 | +0.381 | 0.041 |
| CLDN4 vs PDCD1 | 29 | +0.087 | 0.65 |
| CLDN4 vs TACSTD2 (companion; also floor) | 29 | +0.227 | 0.24 |
| CLDN4 vs EPCAM (floor) | 29 | +0.086 | 0.66 |
| CLDN4 vs GFAP | 29 | +0.028 | 0.89 |
| CD274 vs T/NK (sanity) | 29 | +0.499 | 0.0059 |
| T/NK vs IFN (sanity) | 29 | +0.712 | 1.5×10⁻⁵ |
| T/NK vs PTPRC (sanity) | 29 | +0.673 | 6.3×10⁻⁵ |
| CD8A vs T/NK (sanity; CD8A is in T/NK) | 29 | +0.902 | 2.4×10⁻¹¹ |

**Requested leftover is null.** CLDN4 does not track CD274, T/NK, IFN, or CD8A at n=29. Immune-axis sanity checks are intact (CD274, T/NK, IFN, and PTPRC co-vary with each other), so the null is not a broken matrix.

MHC-I is a nominal p=0.041 in the **wrong direction for an immune-low story** (positive, not T/NK-low) and is **not** the leftover question. It is reported so it is not hidden. It is not a CD274 or T/NK claim.

## After PTPRC / GFAP

Partial Spearman = Pearson of rank residuals. n=29.

| pair | partial ρ \| PTPRC | p | partial ρ \| GFAP | p |
|---|---:|---:|---:|---:|
| CLDN4 vs CD274 | +0.343 | 0.069 | +0.302 | 0.11 |
| CLDN4 vs T/NK | −0.173 | 0.37 | −0.139 | 0.47 |
| CLDN4 vs IFN | +0.236 | 0.22 | +0.186 | 0.34 |
| CLDN4 vs CD8A | −0.108 | 0.58 | −0.076 | 0.70 |
| CLDN4 vs MHC-I | +0.471 | 0.0099 | +0.389 | 0.037 |

CD274 / T/NK / IFN stay null after CD45 or GFAP. MHC-I stays nominally positive; that is still not the leftover axis.

## Median split (CLDN4-high vs CLDN4-low)

Median CLDN4 log2(CPM+1) = 0.567 (raw median 15 counts). High n=15, low n=14. Cliff δ = P(high>low) − P(high<low). This split is **floor vs slightly-less-floor**, not a biological CLDN4-high group.

| axis | nH / nL | δ | MWU p |
|---|---|---:|---:|
| CD274 | 15 / 14 | +0.38 | 0.085 |
| T/NK | 15 / 14 | −0.038 | 0.88 |
| IFN | 15 / 14 | +0.32 | 0.14 |
| CD8A | 15 / 14 | +0.11 | 0.62 |
| NKG7 | 15 / 14 | −0.11 | 0.63 |
| MHC-I | 15 / 14 | +0.41 | 0.064 |
| PTPRC | 15 / 14 | +0.019 | 0.95 |

## Floor-only and arm (honest small n)

| stratum | n | CLDN4–CD274 ρ (p) | CLDN4–T/NK ρ (p) | CLDN4–IFN ρ (p) | CLDN4–CD8A ρ (p) |
|---|---:|---|---|---|---|
| CLDN4 CPM < 1 | 26 | +0.350 (0.079) | +0.015 (0.94) | +0.285 (0.16) | +0.045 (0.83) |
| CLDN4 CPM ≥ 1 | **3** | **empty** (n<5) | empty | empty | empty |
| Neoadjuvant A (on-tx) | 14 | +0.341 (0.23) | −0.130 (0.66) | +0.064 (0.83) | −0.064 (0.83) |
| Adjuvant B (pre-tx) | 15 | +0.282 (0.31) | −0.154 (0.58) | +0.432 (0.11) | +0.039 (0.89) |

The three tumors above CPM 1 are Pt33_A (3.04), Pt15_B (8.03), and Pt20_B (1.32). They are the same leverage samples the B5 OS audit already flagged. They are not a CLDN4-high immune contrast.

## How to read this

- **Matrix is public, n=29, not an empty GEO record.** The empty object is the series-matrix expression table. The HUGO xlsx is the usable matrix.
- **CLDN4 itself is almost empty.** Usable expression is **n=3** (CPM ≥ 1). Do not write this as a CLDN4-high vs immune test with n=29.
- **CLDN4 vs CD274 does not hold** (ρ=+0.29, p=0.12). Same for T/NK, IFN, and CD8A.
- **IFNG is a row, not a signal** (usable n=0). The Ayers score is 6/6 present and is still null vs CLDN4.
- **RECIST / OS are empty here** (n=0 on this page). This leftover does not re-score survival.
- MHC-I p=0.041 is a nominal extra row, positive, and is not upgraded to a leftover claim.

This is one recurrent-GBM pembrolizumab count matrix. It adds **no** CLDN4–immune or CLDN4–CD274 weight. The honest leftover result is public n=29 / usable CLDN4 n=3 / requested correlations null.

## Files

- `analyze.py` — download series matrix + HUGO xlsx to `/tmp/gse121810`, then score
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `label_inventory.tsv`, `gene_coverage.tsv`, `detectability.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd274.png`
- `figures/fig2_cldn4_vs_tnk.png`
- `figures/fig3_cldn4_vs_cd8a.png`
- `figures/fig4_spearman_forest.png`
- `figures/fig5_detectability.png`

```bash
python3 methods/gse121810_cldn4/analyze.py
```
