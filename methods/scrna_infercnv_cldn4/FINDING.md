# inferCNV-proxy malignant CLDN4 vs same-patient T/NK

ADDITIVE. **CLDN4 only.** Public GEO UMI, GSE207422 (Hu et al., *Genome Medicine* 2023, PMID 36869384). Author CopyKAT barcodes are not on GEO. This is a window-smoothed expression CNV proxy, not the R `infercnv` or `copykat` packages.

Malignant = marker epithelium **and** CNV score > 95th percentile of stromal (fibroblast + endothelial) reference. T/NK = lineage T ∪ NK. Score = mean `log1p(CP10k)` CLDN4, with % UMI>0 as a second metric. Unit = patient. Samples below the cell floor drop out.

## Verdict

CNV-aneuploid epithelium is CLDN4-high versus the same patient’s T/NK in every post-treatment pair. That is a compartment contrast, not an ICI-response claim.

- **inferCNV-like stromal p95, post-treatment, ≥1 cell/arm: n=12** (NMPR 8 / MPR 4, including pCR P06). Malignant median 1.429 vs T/NK 0.012. 12/12 patients malignant > T/NK. Wilcoxon signed-rank (greater) p=2.4×10⁻⁴.
- **Same definition, ≥5 cells/arm: n=11.** Drops P13 (3 inferCNV-malignant cells). 11/11 still malignant > T/NK. p=4.9×10⁻⁴.
- CopyKAT-like p95 keeps all **n=12** at the ≥5 floor (median 1.706 vs 0.012, p=2.4×10⁻⁴).
- Versus T/NK **fraction** the Spearman is a null (inferCNV n=12: ρ=−0.140, p=0.665; n=11: ρ=−0.055, p=0.873).
- NMPR vs MPR malignant CLDN4 is **not** NMPR>MPR (inferCNV n=12: 1.288 vs 1.768, p=0.214). Do not reuse the TACSTD2 A3 direction for CLDN4.

GSE131907 was not required; GSE207422 has public UMI, paired T/NK, and MPR labels.

## Honest n

| Malignant def | Floor | n (NMPR/MPR) | Dropped | Why dropped |
|---|---|---|---|---|
| inferCNV-like p95 | ≥1 | 12 (8/4) | — | — |
| inferCNV-like p95 | ≥5 | 11 (7/4) | P13 | 3 CNV-malignant cells |
| CopyKAT-like p95 | ≥1 or ≥5 | 12 (8/4) | — | — |
| marker epithelium (sensitivity) | ≥1 or ≥5 | 12 (8/4) | — | not a CNV call |

Labeled post-treatment samples = 12. Pre-treatment biopsies (P01, P05, P08) are excluded from the primary pair. Residual epithelium is thin after MPR: P06 has 5 inferCNV-malignant cells; P13 has 3.

## Paired CLDN4 (same patient)

Post-treatment. Wilcoxon signed-rank, alternative = malignant > T/NK.

| Malignant def | n | median mal | median T/NK | # mal>T/NK | p | %pos mal | %pos T/NK |
|---|---|---|---|---|---|---|---|
| inferCNV-like p95 | 12 | 1.429 | 0.012 | 12/12 | 2.4×10⁻⁴ | 86.5 | 0.66 |
| inferCNV-like p95 | 11 | 1.428 | 0.014 | 11/11 | 4.9×10⁻⁴ | 84.4 | 0.72 |
| CopyKAT-like p95 | 12 | 1.706 | 0.012 | 12/12 | 2.4×10⁻⁴ | 97.7 | 0.66 |
| marker epithelium | 12 | 1.529 | 0.012 | 12/12 | 2.4×10⁻⁴ | 85.5 | 0.66 |

Cell-level (descriptive; n is huge): post inferCNV-malignant 1,178 cells, median CLDN4 1.34, 81.8% pos vs T/NK 33,760 cells, median 0, 1.2% pos.

## Per-patient post-treatment (inferCNV-like)

| Patient | Response | Histology | n mal | mal CLDN4 | mal %pos | n T/NK | T/NK CLDN4 | T/NK %pos |
|---|---|---|---|---|---|---|---|---|
| P02 | NMPR | Adeno | 15 | 1.428 | 53.3 | 2899 | 0.014 | 0.72 |
| P03 | MPR | Squamous | 194 | 0.978 | 69.1 | 5060 | 0.016 | 0.91 |
| P04 | NMPR | Adeno | 26 | 1.103 | 50.0 | 5616 | 0.005 | 0.32 |
| P06 | MPR (pCR) | Adeno | 5 | 2.486 | 100 | 2924 | 0.004 | 0.21 |
| P07 | NMPR | Squamous | 493 | 1.328 | 94.1 | 1087 | 0.217 | 9.66 |
| P09 | NMPR | Squamous | 115 | 1.429 | 88.7 | 3283 | 0.046 | 3.17 |
| P10 | NMPR | Squamous | 35 | 0.705 | 40.0 | 4080 | 0.014 | 0.86 |
| P11 | MPR | Adeno | 115 | 1.991 | 98.3 | 1690 | 0.046 | 2.72 |
| P12 | NMPR | Adeno | 95 | 0.731 | 34.7 | 1386 | 0.009 | 0.51 |
| P13 | NMPR | Squamous | 3 | 1.931 | 100 | 1347 | 0.002 | 0.15 |
| P14 | MPR | Squamous | 50 | 1.618 | 96.0 | 2863 | 0.010 | 0.59 |
| P15 | NMPR | Adeno | 32 | 1.648 | 84.4 | 1525 | 0.006 | 0.33 |

P13 is in the n=12 row and out of the n=11 row. P07 T/NK CLDN4 is the highest immune value (ambient / doublet risk) and is still below that patient’s malignant mean.

## Secondary: vs T/NK fraction, and NMPR vs MPR

These are reported so they are not silently reused.

| Malignant def | n | ρ vs T/NK frac | p(ρ) | mean NMPR | mean MPR | Δ | p |
|---|---|---|---|---|---|---|---|
| inferCNV-like p95 | 12 | −0.140 | 0.665 | 1.288 | 1.768 | −0.480 | 0.214 |
| inferCNV-like p95 | 11 | −0.055 | 0.873 | 1.196 | 1.768 | −0.572 | 0.164 |
| CopyKAT-like p95 | 12 | +0.210 | 0.513 | 1.711 | 1.699 | +0.012 | 1.000 |
| marker epithelium | 12 | +0.035 | 0.914 | 1.480 | 1.502 | −0.023 | 1.000 |

n=12 is small. A true ρ=−0.45 has two-sided Spearman p≈0.14 here. Direction and p are separate.

## Methods (short)

- Matrix: 24,292 genes × 92,330 barcodes. Lineage = argmax of `log1p(CP10k)` modules (T, NK, B, plasma, myeloid, neutrophil, mast, epithelial, fibroblast, endothelial).
- CNV keep-set = epithelium + stroma + 2,500 T/myeloid cells. Reference = 1,303 stromal cells. 11,572 genes with coordinates, mean UMI>0.05, ≥20 expressing cells. `log1p(CP10k)`, subtract reference (mean for inferCNV-like, median for CopyKAT-like), chromosome-wise window=25, score = mean|smooth| or sum|smooth|. Aneuploid = epithelium and score > stromal p95 (1,517 and 6,080 cells).
- T/NK fraction denominator = all cells in that sample.
- Figure: `fig_cldn4_infercnv_vs_tnk.png`.

## Files

| File | Role |
|---|---|
| `FINDING.md` | this note |
| `analyze.py` / `download.py` | GEO download + CNV proxy + CLDN4 tests |
| `per_sample.tsv` | patient-level CLDN4 and T/NK |
| `summary.json` | all n / ρ / p including pre-treatment |
| `cell_calls.tsv.gz` | barcode flags and CLDN4 |
| `sample_metadata.tsv` | GEO sheet (pCR grouped with MPR) |
| `assets/gene_chr.tsv` | gene coordinates for the CNV window |
