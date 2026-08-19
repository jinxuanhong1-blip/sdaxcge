# RESULTS — additive public Visium, CLDN4-only

Public processed Visium only. **CLDN4-only** (no dual-high / TACSTD2 gate). No private 8-KL. No invented accessions or statistics.

## Sources

| ID | Accession | Design | Slide chemistry used for coordinates | Status |
|---|---|---|---|---|
| A | [GSE277206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277206) | 10x Visium CytAssist FFPE, never-smoker LUAD progression (MIA-034, MIA-039; small n) | visium-v5 / 11 mm (barcode whitelist; GEO uploaded H5 only) | ran (2 sections) |
| B | [Zenodo 13337961](https://zenodo.org/records/13337961) / [10.3389/fimmu.2024.1430163](https://doi.org/10.3389/fimmu.2024.1430163) | lepidic vs solid LUAD spatial transcriptome | visium-v4 / CytAssist 6.5 mm (record text says 11 mm; uploaded barcodes are the 4,992-spot 6.5 mm set) | ran (2 sections) |

## Locked methods

- QC: spots with ≥200 detected genes and a matched Visium array coordinate.
- Expression: log1p(CP10K).
- Same-spot test: Spearman **CLDN4 vs CD8A** on all QC spots.
- Epithelial-like: section Q3+ of the mean of available `EPCAM`, `KRT8`, `KRT18`, `KRT19`, `CDH1`, `KRT7` (KRT18 is absent from the Visium Human Transcriptome Probe Set v2.0 matrices used here).
- CLDN4-high / low: Q4 / Q1 of CLDN4 **among epithelial-like spots**.
- CD8A-high: section Q4 of CD8A among all QC spots **if Q3 > 0**. If Q3 == 0 (zero-inflated CD8A; Q4 is not a high tail), CD8A-high = CD8A > 0.
- Nearest CD8A-high distance: hex-aware Euclidean distance in µm (100 µm Visium pitch). Self is excluded when a query spot is itself CD8A-high.
- Neighbor CD8A: mean CD8A of hex ring-1 neighbors `(row, col±2)` and `(row±1, col±1)`; require ≥3 neighbors.
- KRT8 residualization: (i) partial Spearman of CLDN4 vs CD8A controlling for KRT8 ranks; (ii) linear residual of log CLDN4 on log KRT8, then the same Q4/Q1 distance and neighbor tests on residual quartiles among epithelial-like spots.
- Q4 vs Q1: two-sided Mann–Whitney U. Section is the unit; n=2 per source is too small for a signed-rank across sections.
- Maps use array coordinates recovered from the public 10x barcode inclusion list. Section-specific H&E / `tissue_positions` pixel maps were **not** deposited with these count matrices.

## Per-section results (computed)

| Source | Section | Histology | n QC | n epi-like | Spearman CLDN4 vs CD8A | Partial \| KRT8 | Nearest CD8A-high µm (Q4 vs Q1) | Neighbor CD8A (Q4 vs Q1) |
|---|---|---|---:|---:|---|---|---|---|
| GSE277206 | GSE277206_MIA-034 | MIA | 10674 | 2669 | ρ=-0.001, p=0.8845 (n=10674) | ρ=0.118, p=3.32e-34 | med 173.2 vs 100.0 (Δ=73.2, p=5.55e-15) | med 0.000 vs 0.234 (Δ=-0.234, p=6.74e-31) |
| GSE277206 | GSE277206_MIA-039 | MIA | 4758 | 1190 | ρ=0.077, p=9.71e-08 (n=4758) | ρ=0.087, p=1.72e-09 | med 100.0 vs 173.2 (Δ=-73.2, p=2.43e-08) | med 0.077 vs 0.000 (Δ=0.077, p=0.1709) |
| zenodo_13337961 | Zenodo13337961_lepidic | lepidic | 4991 | 1248 | ρ=0.097, p=7.41e-12 (n=4991) | ρ=0.029, p=0.04304 | med 100.0 vs 100.0 (Δ=-0.0, p=4.86e-10) | med 0.192 vs 0.123 (Δ=0.069, p=1.98e-12) |
| zenodo_13337961 | Zenodo13337961_solid | solid | 3701 | 926 | ρ=-0.202, p=1.59e-35 (n=3701) | ρ=-0.141, p=5.63e-18 | med 173.2 vs 100.0 (Δ=73.2, p=0.06775) | med 0.033 vs 0.046 (Δ=-0.013, p=0.0003634) |

### KRT8-residual Q4 vs Q1 (epithelial-like residual quartiles)

| Source | Section | Resid nearest µm Δ (Q4−Q1) | p | Resid neighbor CD8A Δ (Q4−Q1) | p |
|---|---|---:|---:|---:|---:|
| GSE277206 | GSE277206_MIA-034 | 73.2 | 1.08e-05 | -0.206 | 1.61e-12 |
| GSE277206 | GSE277206_MIA-039 | -73.2 | 2.45e-07 | 0.076 | 0.1562 |
| zenodo_13337961 | Zenodo13337961_lepidic | 0.0 | 6.65e-07 | 0.071 | 1.01e-12 |
| zenodo_13337961 | Zenodo13337961_solid | 73.2 | 0.04334 | -0.012 | 0.001444 |

## What this does not claim

- It does not use private 8-KL or any controlled/unpublished matrix.
- It does not treat n=2 sections as a multi-cohort meta-analysis.
- Same-spot CLDN4 vs CD8A is composition-sensitive; the KRT8 residual / partial is the control, not a proof of a CLDN4-specific barrier.
- Array maps are barcode-whitelist coordinates, not the unpublished section H&E.

## Reproduce

```bash
python3 methods/visium_public_luad/download.py
python3 methods/visium_public_luad/analyze.py
```

## Files

- `methods/visium_public_luad/analyze.py`
- `methods/visium_public_luad/download.py`
- `methods/visium_public_luad/coords/visium_v4_cytassist_6p5mm.csv`
- `methods/visium_public_luad/coords/visium_v5_cytassist_11mm.csv`
- `results/visium_public_luad/tables/section_stats.tsv`
- `results/visium_public_luad/tables/summary.json`
- `results/visium_public_luad/maps/`

