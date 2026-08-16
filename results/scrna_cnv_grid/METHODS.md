# Methods — malignant-definition grid

Public processed matrices only. No FASTQ. No private author RDS.

## GSE207422

- Hu et al., *Genome Medicine* 2023, PMID 36869384.
- GEO UMI: 24,292 genes × 92,330 barcodes.
- Sample sheet: 15 rows. Primary contrast = 12 post-treatment resections. pCR (P06) grouped with MPR (n=4) vs NMPR (n=8). Three pre-treatment biopsies are excluded from the MPR test.
- Author CopyKAT / epithelium barcode IDs are not deposited. DRMref Seurat labels (30,877 cells; 2,051 `Malignant cells`) are used as a **third-party** axis, not as author CopyKAT.

## GSE241934

- Zhang et al., *Cell Reports Medicine* 2024 (NEOTIDE/CTONG2104).
- Public MTX + `major.cell.type` + `Pathological Response`.
- IIT: 11 EGFR-mutant trial tumors (MPR 4, non-MPR 7). Real-world: 34 EGFR-WT tumors (pCR grouped with MPR).
- Author epithelium = `major.cell.type == Epi` (IIT 1,699 cells; paper reports CopyKAT malignant N=1,669).

## Marker lineage

Canonical `log1p(CP10k)` module scores. Assigned lineage = argmax. T vs NK broken by CD3E when scores are close. Normal-lung program: SFTPA1/2, SFTPB/C, AGER, SCGB1A1, SCGB3A1/2, TPPP3, FOXJ1, CAPS.

## CopyKAT / inferCNV-style score

Not the R `copykat` or `infercnv` packages. Window-smoothed expression CNV:

1. Diploid reference = marker (GSE207422) or author (GSE241934) fibroblasts + endothelia. Stromal reference is subsampled to ≤2,500 cells on the larger real-world matrix.
2. Genes with genomic coordinates, mean UMI > 0.05, ≥20 expressing cells in the keep-set.
3. `log1p(CP10k)`, subtract the median of reference cells per gene.
4. Order by chromosome, start. Sliding-window mean, window = 25 genes. No smoothing across chromosomes.
5. Per-cell score = sum of |smoothed values|.
6. Aneuploid call = epithelial (author and/or marker) AND score > 95th percentile of reference-cell scores. p90 and 1-D 2-means cuts are reported as sensitivities.

## TACSTD2 metrics (malignant-only, per sample)

- mean `log1p(CP10k)` (primary)
- % cells with UMI ≥ 1
- pseudobulk CPM (`sum TACSTD2 UMI / sum library UMI × 10⁶`)

Samples with zero cells in a definition contribute NaN and drop out of that test. A sensitivity requires ≥5 malignant cells.

## T/NK fraction (denominator = all cells in the sample)

- lineage T ∪ NK
- UMI gate: CD3E≥1 OR CD8A≥1 OR NKG7≥1
- GSE207422 DRMref: CD8+ T + CD4+ T + NK among DRMref-annotated cells
- GSE241934 author: `major.cell.type` in {T, NK}

## Tests

- **NMPR vs MPR:** exact two-sided Wilcoxon by enumerating assignments when C(n,k) ≤ 30,000; otherwise asymptotic Mann–Whitney. Unit = sample.
- **vs T/NK:** Spearman ρ of the sample-level TACSTD2 metric vs each T/NK fraction.

## Direction (pre-specified)

User direction = **NMPR mean > MPR mean** and **ρ < 0**. A definition **recovers the user direction** when both are true on the primary metric (mean log1p CP10k, no cell-count floor). p-values are reported honestly; direction is not a significance claim.
