# Methods — Seurat / Harmony integration of the triple that differs

## Design (locked)

| Piece | Choice |
|---|---|
| Cohorts | GSE123902 + GSE131907 + GSE205335 only |
| Gate | PR #459 locked units (n=56): 13 donors + 21 samples + 22 patients |
| Split | Within-dataset Seurat malignant CLDN4 %pos Q4 vs Q1 |
| Malignant | GSE123902 marker-malignant (`EPCAM\|KRT8\|KRT18\|KRT19>0` and `PTPRC==0`); GSE131907 author `Cell_subtype` in {Malignant cells, tS1, tS2, tS3}; GSE205335 author `lineage.sub == Malignant cells` |
| T/NK | Locked full-sample `frac_tnk` (not the capped object) |
| Integration | Seurat v5 `IntegrateLayers(HarmonyIntegration)` on `dataset` layers |
| Scores | `AddModuleScore` on joined RNA: Hallmark IFNα∪IFNγ and custom MHC-I + MHC-II; CLDN4 excluded |
| Unit | **patient** (GSE123902 donor; GSE131907 GEO Sample = patient on this gate; GSE205335 patient). Cells are not n. |
| Dataset covariate | Residual Spearman (both axes residualized on dataset) and OLS `y ~ CLDN4_z + dataset` |
| Cap | ≤80 malignant + ≤80 T/NK cells per patient (memory). Means use the capped malignant cells. |

## Reproduce

```bash
export R_LIBS_USER=/tmp/r_libs
bash methods/seurat_integrate_123902_131907_205335_cldn4/scripts/00_download.sh
Rscript methods/seurat_integrate_123902_131907_205335_cldn4/scripts/01_keep_cells.R
python3 methods/seurat_integrate_123902_131907_205335_cldn4/scripts/02_extract.py
Rscript methods/seurat_integrate_123902_131907_205335_cldn4/scripts/03_seurat_integrate.R
```

Python is extract-only for the GSE131907 gene × cell TSV and the GSE123902 dense CSVs.
Primary integration and tests are R + Seurat / Harmony.
If Seurat cannot install, stop. Do not fall back to a Python-only primary.
