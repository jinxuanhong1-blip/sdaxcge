# FINDING — Seurat / Harmony triple GSE123902 + GSE131907 + GSE205335, CLDN4-only

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. Not a bigger merge.
The triple that differs is given (PR #459: malignant CLDN4 %pos vs T/NK,
Spearman ρ=−0.522, N=56; Q4 vs Q1 r=−0.735, n=28). Thesis is already
correct and is **not** re-derived: CLDN4-high malignant cells sit with
lower T/NK and lower own IFN/MHC. This folder adds a **Seurat v5
`IntegrateLayers` (Harmony)** object of the same three public GEO sets.

Primary engine: **R + Seurat / Harmony**. Python is extract-only (GSE131907
UMI TSV; GSE123902 dense CSVs). No Python-only primary. No GSE148071.
No GSE189357.

Honest unit = **patient** (GSE123902 donor; GSE131907 GEO Sample, unique
patient numbers on the gated tumor set; GSE205335 patient). Cells are
counts, not replicates.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE123902 + GSE131907 + GSE205335 only | Not 148071 / 189357 / 207422 |
| Integration | Seurat v5 `IntegrateLayers` Harmony on `dataset` layers | Integration is for UMAP; scores use joined RNA log-norm |
| CLDN4 | Single gene; %pos = counts > 0 in malignant cells | Capped ≤80 malignant cells / patient |
| T/NK fraction | Locked full-sample `frac_tnk` (PR #459) | Not recomputed from the capped object |
| IFN / MHC | `AddModuleScore`; Hallmark IFNα∪IFNγ; custom MHC-I + MHC-II | CLDN4 excluded from both modules |
| Quartiles | Within-dataset Seurat CLDN4 %pos | Mid quartiles unused in Q4 vs Q1 |
| Dataset covariate | Residual Spearman + OLS `y ~ CLDN4_z + dataset` | Linear; descriptive |
| Model | Within-dataset Spearman; DerSimonian–Laird on Fisher-z; stacked MW | Three datasets; p-values are descriptive |

## Honest n

_Filled after the Seurat run. See `results/tables/n_honest.tsv`._

## Results

_Pending `03_seurat_integrate.R`. Do not quote numbers from this stub._

## What this is not

- Not a dual-high TACSTD2×CLDN4 gate.
- Not GSE148071 or GSE189357 or any fourth cohort.
- Not a cell-level Wilcoxon sold as n.
- Not evidence that CLDN4 *causes* T/NK loss or IFN/MHC drop.
- Not a Python Scanpy primary.
- Not a re-derivation of the thesis.

## Files

- `scripts/03_seurat_integrate.R` — Seurat merge, Harmony, module scores, patient tests
- `results/tables/patient_scores.tsv`
- `results/tables/patient_level_tests.tsv`
- `results/tables/dataset_covariate_ols.tsv`
- `results/figures/`

## Reproduce

```bash
export R_LIBS_USER=/tmp/r_libs
bash methods/seurat_integrate_123902_131907_205335_cldn4/scripts/00_download.sh
Rscript methods/seurat_integrate_123902_131907_205335_cldn4/scripts/01_keep_cells.R
python3 methods/seurat_integrate_123902_131907_205335_cldn4/scripts/02_extract.py
Rscript methods/seurat_integrate_123902_131907_205335_cldn4/scripts/03_seurat_integrate.R
```
