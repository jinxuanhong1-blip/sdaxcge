# GSE207422 scVI/Harmony latent — malignant CLDN4 (not TACSTD2 redo)

Additive methods only. Public GEO UMI. **CLDN4 is the score.** TACSTD2 is a companion gene and is **never a gate**. The given A3 TACSTD2 analysis and the marker-only CLDN4 dualhigh slice are not re-argued.

**Verdict (honest n):** after a real scVI latent on GSE207422 (batch = 10x sample), Leiden-cluster malignant CLDN4 vs T/NK in the 12 post-treatment patients is ρ=-0.03, p=0.93, n=10. Marker-lineage malignant CLDN4 vs T/NK is ρ=0.07, p=0.83, n=11. A3-malignant CLDN4 vs T/NK is ρ=0.27, p=0.42, n=11. None of these is a significant anti-correlation. NMPR vs MPR on cluster-malignant CLDN4 is mean 1.398 vs 1.483 (Δ=-0.085); exact p=0.67; n=7 vs 3.

## What this is (and is not)

| Existing slice | This slice |
|---|---|
| scVI combo GSE207422+GSE241934 IIT scored **TACSTD2** (n=26) | GSE207422 only; scored **CLDN4** |
| Marker-only GSE207422 malignant CLDN4 (no latent) | Same public UMI, but malignant calls also from the **latent Leiden** |
| A3 TACSTD2 NMPR/MPR / T/NK | Taken as given; not re-run as the primary |

## Data and n

- Public GEO UMI only: **92330** cells × **24292** genes. After QC (≥200 UMI, ≥200 genes): **92330** cells. Raw GSA-Human HRA001033 was not used. Author CopyKAT barcodes are not on GEO.
- **Header n = 12 post-treatment patients** (MPR n=4 including pCR P06; NMPR n=8). The three pre-treatment biopsies (P01/P05/P08) are excluded from primary tests. All 15 samples are in the tables.
- Eligible for a Spearman: ≥10 malignant-like cells and ≥20 T/NK cells and a finite CLDN4 mean. Cell-level p-values are not reported.
- Marker lineage after QC: {"B/Plasma": 11295, "Endothelial": 2510, "Epithelial": 14067, "Fibroblast": 1114, "Mast": 1085, "Myeloid": 23808, "T/NK": 38451}.
- Leiden cluster lineage after scVI: {"B/Plasma": 11081, "Endothelial": 240, "Epithelial": 12118, "Fibroblast": 868, "Mast": 842, "Myeloid": 28627, "T/NK": 38554}.
- Cluster-malignant ineligible post samples (n_malignant < 10): **P06 / BD_immune06 (pCR; 5 cells)** and **P13 / BD_immune13 (NMPR; 3 cells)**. Primary Spearman n=10/12.
- Marker-malignant ineligible post: **P11 / BD_immune11 (MPR; 2 cells)**. Spearman n=11/12.
- A3-malignant ineligible post: **P11 / BD_immune11 (2 cells)**. Spearman n=11/12. This A3 call uses the same zero-UMI normal-lung panel as the given A3 TACSTD2 slice, but the epithelial parent is the broader combo-style argmax (not the Hu-threshold lineage). It is a sensitivity, not a dualhigh redo.
- Leiden epithelial clusters: 4, 11, 16 called malignant (n=7,369 + 2,069 + 732 = 10,170). Cluster 12 is epithelial but normal-lung-high (mean 1.87) and is **not** malignant.

## Integration

- Method: **scVI**
- Batch key: `sample` (the 15 10x libraries)
- HVGs: 2000 (`seurat_v3`, batch-aware); latent dim 10
- Cells in the latent: 92330
- scVI version / note: 1.5.0.post1

CLDN4 and T/NK are scored on raw log1p(CP10k) / lineage fractions, not on the latent coordinates. The latent is used to build a joint embedding and a Leiden malignant call. That is the additive piece versus the marker-only CLDN4 slice.

## Definitions

| Item | Rule |
|---|---|
| Marker malignant-like | argmax lineage = Epithelial AND normal-lung score ≤ epithelial 75th percentile |
| Cluster malignant | Leiden cluster on `X_int` whose mean marker score is Epithelial and whose mean normal-lung score is ≤ the epithelial-cluster 75th percentile |
| A3-malignant | epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (same rule as the given A3 TACSTD2 slice; **not** CopyKAT). Sensitivity only. |
| CLDN4 mean | sample mean `log1p(CP10k)` inside the named malignant call |
| T/NK | matching lineage fraction of all QC cells |
| CXCL13+ | (T/NK) AND CXCL13 UMI ≥ 1 / all QC cells |
| cyto-high | T/NK AND mean log1p(GZMB, GZMA, PRF1, IFNG, NKG7) ≥ pooled T/NK 75th percentile / all QC cells |
| TACSTD2 | companion only; not a gate |

## Primary — post n=12, malignant CLDN4 vs immune

| Malignant call | vs T/NK | vs CXCL13+ | vs cyto-high |
|---|---|---|---|
| Cluster (latent Leiden) | ρ=-0.03, p=0.93, n=10 | ρ=-0.27, p=0.45, n=10 | ρ=-0.13, p=0.73, n=10 |
| Marker lineage | ρ=0.07, p=0.83, n=11 | ρ=0.31, p=0.35, n=11 | ρ=0.34, p=0.31, n=11 |
| A3-malignant (sensitivity) | ρ=0.27, p=0.42, n=11 | ρ=0.41, p=0.21, n=11 | ρ=0.22, p=0.52, n=11 |

## Primary — NMPR vs MPR (post only)

| Score | Result |
|---|---|
| Cluster-malignant CLDN4 mean | mean 1.398 vs 1.483 (Δ=-0.085); exact p=0.67; n=7 vs 3 |
| Marker-malignant CLDN4 mean | mean 0.633 vs 0.324 (Δ=+0.309); exact p=0.38; n=8 vs 3 |
| A3-malignant CLDN4 mean | mean 0.716 vs 0.323 (Δ=+0.393); exact p=0.38; n=8 vs 3 |
| Cluster T/NK fraction | mean 0.404 vs 0.537 (Δ=-0.132); exact p=0.68; n=8 vs 4 |
| Marker T/NK fraction | mean 0.406 vs 0.535 (Δ=-0.129); exact p=0.68; n=8 vs 4 |

## Sensitivity — all 15 samples (3 pre + 12 post)

Pre-treatment biopsies are not a second clinical unit. They are shown so the 15-library latent is not silently subset.

| Malignant call | vs T/NK | vs CXCL13+ | vs cyto-high |
|---|---|---|---|
| Cluster (latent Leiden) | ρ=-0.13, p=0.68, n=12 | ρ=-0.41, p=0.19, n=12 | ρ=-0.20, p=0.53, n=12 |
| Marker lineage | ρ=-0.11, p=0.71, n=14 | ρ=0.10, p=0.72, n=14 | ρ=0.13, p=0.67, n=14 |

## Companion TACSTD2 (not a gate; not a redo)

Reported so CLDN4 can be read next to TACSTD2 on the **same latent malignant call**. This does not replace the given A3 TACSTD2 write-up or the scVI combo TACSTD2 score.

- Cluster-malignant CLDN4 vs TACSTD2 (post eligible): ρ=0.50, p=0.14, n=10
- Cluster-malignant TACSTD2 vs T/NK (post eligible): ρ=-0.15, p=0.68, n=10
- Marker-malignant TACSTD2 vs T/NK (post eligible): ρ=0.20, p=0.56, n=11

## Honest limits

1. Header n=12 (4 MPR vs 8 NMPR) is the cohort. Spearman |ρ|≈0.45 has two-sided p≈0.14 at n=12. Eligible n can be lower when a malignant call is empty.
2. A3-malignant-like empties some MPR residual tumors (normal-lung program). That is reported, not patched.
3. Leiden-cluster malignant is a latent annotation, not Hu et al. CopyKAT. Residual unmarked epithelium can leak.
4. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is not a gate.
5. Expression is log1p(CP10k) from the public UMI. The latent does not change the CLDN4 numbers; it changes who is called malignant.
6. This is GSE207422 only. It is not the GSE207422+GSE241934 scVI combo.

## Files

- `sample_table_cluster.tsv` / `sample_table_marker.tsv` / `sample_table_a3.tsv`
- `cluster_annotation.tsv` / `spearman.tsv` / `nmpr_vs_mpr.tsv` / `summary.json`
- Figures: cluster and marker CLDN4 vs T/NK; UMAPs (sample, marker lineage, cluster lineage)

## Reproduce

```bash
python3 methods/scrna_scvi_cldn4/download.py
python3 methods/scrna_scvi_cldn4/prepare.py
python3 methods/scrna_scvi_cldn4/integrate_score.py
```

scVI is preferred. If training fails, the script falls back to Harmony and records the error.

