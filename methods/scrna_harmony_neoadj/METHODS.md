# METHODS — public Harmony integration of neoadjuvant lung scRNA

Public processed UMI / author matrices only. No FASTQ. User A3 (GSE207422 CopyKAT) is taken as given and was not re-run.

## Series considered

| Accession | Public object | Why in / out |
|---|---|---|
| GSE207422 | `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` (gene × cell UMI) + sample xlsx | In. Epithelium + immune. Neoadjuvant MPR/NMPR/pCR. LUAD (Adeno) and LUSC (Squamous). |
| GSE241934 | IIT + Real MTX + author meta (`major.cell.type`, MPR) | In. Author Epi / T / NK. Neoadjuvant IO+chemo. Histology LUAD / ASC (no LUSC). |
| GSE291670 | Six 10x MTX in `GSE291670_RAW.tar` | In. Tumor nuclei, MPR vs Non-MPR in titles. NSCLC not split on GEO. |
| GSE205335 | `GSE205335_Lung_IO_CellIdentity.txt.gz` inventoried; UMI is `*_UMI_matrix.rds.gz` | **Out of the Harmony object.** Author epithelium / malignant cells exist, but the series is palliative ICI with RECIST, not neoadjuvant MPR. Mixing RECIST into an MPR model is an endpoint swap. RDS-only UMI was not pulled. |

## Shared gene space

Gene symbols were intersected across the three included matrices. `TACSTD2` and `CLDN4` are required. A core lineage/target panel plus highly variable genes among the intersection (capped) were used for PCA / Harmony. Library size is author `nCount_RNA` (GSE241934) or the sum of deposited gene UMIs (GSE207422, GSE291670).

## Integration

CP10k → log1p → gene z-score → PCA (30 PCs) → `harmonypy` with batch = dataset (`GSE207422`, `GSE241934_IIT`, `GSE241934_REAL`, `GSE291670`). Leiden / UMAP on Harmony PCs (sample-stratified subsample if needed; labels transferred by kNN). Clusters labeled by marker means. Malignant-like = epithelial clusters with a low normal-lung score. T/NK = T + NK clusters.

## Statistics

Patient / post-treatment sample is the unit. GSE207422 MPR tests use the 12 post-treatment surgery samples (pCR counted as MPR). Minimum 20 malignant-like cells; T/NK pairing also ≥20 T/NK.

1. Malignant-like mean log1p(CP10k) `TACSTD2` and `CLDN4`: two-sided MWU (exact when the assignment count is small) and OLS `score ~ C(MPR) + C(cohort)`.
2. Spearman of that score vs joint-embedding T/NK fraction; partial Spearman after residualizing cohort. Unadjusted ρ is reported but is not interpreted without the cohort residual.
3. LUAD vs LUSC where labeled (GSE207422 Pathology; GSE241934 Histology). ASC is not LUSC.

GSE291670 is nuclear RNA; patient-level CP10k scores are not Harmony-corrected, so the cohort term is required. Sensitivities at ≥5 and ≥10 malignant-like cells are stored in `tables/sensitivity_floors.json`.

No cell-level p-values for patient-level labels.
