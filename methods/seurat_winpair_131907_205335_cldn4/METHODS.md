# Methods — Seurat integration of the winning pair

## Design (locked)

| Piece | Choice |
|---|---|
| Cohorts | GSE131907 + GSE205335 only |
| Gate | Author malignant CLDN4, n_mal ≥ 20 (PR #320 tables) |
| Split | Within-cohort Seurat malignant CLDN4 %pos Q4 vs Q1 |
| Malignant | Author `Cell_subtype == Malignant cells` / `lineage.sub == Malignant cells` |
| T/NK | Author T/NK labels; **fraction is the locked full-sample `frac_tnk`** (not the capped object) |
| Integration | Seurat v5 `IntegrateLayers(RPCAIntegration)` on dataset layers |
| Scores | `AddModuleScore` on joined RNA: Hallmark IFNα∪IFNγ and custom MHC-I + MHC-II; CLDN4 excluded from both |
| Unit | GSE131907 = GEO Sample; GSE205335 = patient (tumor GSMs collapsed). Cells are not n. |
| Cap | ≤120 malignant + ≤120 T/NK cells per unit (memory). Means use the capped malignant cells. |

## Reproduce

```bash
bash methods/seurat_winpair_131907_205335_cldn4/scripts/00_download.sh
Rscript methods/seurat_winpair_131907_205335_cldn4/scripts/01_keep_cells.R
python3 methods/seurat_winpair_131907_205335_cldn4/scripts/02_extract_gse131907.py
Rscript methods/seurat_winpair_131907_205335_cldn4/scripts/03_seurat_integrate.R
```

Python is extract-only for the GSE131907 gene × cell TSV. Primary integration
and tests are R + Seurat.
