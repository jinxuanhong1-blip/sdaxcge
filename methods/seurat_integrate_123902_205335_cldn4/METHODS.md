# Methods — Seurat / Harmony pair GSE123902 + GSE205335 (CLDN4-only)

ADDITIVE. **CLDN4 only.** Thesis already correct. Honest unit = **patient**
(GSE123902 donor / LX ID; GSE205335 patient). Do **not** add GSE148071 or
GSE127465. Concordant-4 four-way merge is a different agent.

## Design (locked)

| Piece | Choice |
|---|---|
| Cohorts | GSE123902 + GSE205335 only |
| Integration | Seurat v5 `IntegrateLayers(HarmonyIntegration)` on `dataset` layers. Fallback: `harmony::RunHarmony` on PCA. |
| CLDN4 | Single gene. Mean = log1p CP10k; %pos = counts > 0. Never a lineage marker. |
| GSE123902 malignant | Four-way marker argmax on log1p CP10k (epithelial vs T/NK vs myeloid vs B). Keep if top ≥ 0.12 and top ≥ 1.15 × second. Malignant = epithelial **in tumor** (primary or metastasis). Matched normal dropped. |
| GSE205335 malignant | Author `lineage.sub == Malignant cells` in tumor tissue (normal lung/LN/brain dropped). |
| T/NK fraction | **Full-sample** tumor cells of that patient (GSE123902 marker T/NK; GSE205335 author `lineage.total == T/NK cells`). Not the Harmony-capped object. |
| IFN / MHC / TJ | Mean log1p CP10k of locked genes on **all** malignant cells of that patient (before the UMAP cap). CLDN4 held out of TJ. |
| Quartiles | Within-cohort malignant CLDN4 **%pos** Q4 vs Q1. Mid unused. |
| Unit | Patient. Cells are counts, not n. Eligible: ≥20 malignant and ≥20 T/NK. |
| Harmony cap | ≤150 malignant + ≤150 T/NK per patient (memory / UMAP only). |
| Dual-high | Not defined. No TACSTD2∩CLDN4 gate. |

## Tests (patient unit)

- Spearman of malignant CLDN4 mean or %pos vs full-sample T/NK fraction (Fisher-z 95% CI), within cohort and DerSimonian–Laird on Fisher-z.
- Spearman of malignant CLDN4 vs IFN / MHC-I/APM / TJ (CLDN4 held out).
- If eligible n ≥ 8: Q4 vs Q1 Mann–Whitney on IFN / MHC / TJ / T/NK (within-cohort quartiles; stacked).

p-values are descriptive. Cell-level tests are not the claim.

## Reproduce

```bash
bash methods/seurat_integrate_123902_205335_cldn4/scripts/00_download.sh
Rscript methods/seurat_integrate_123902_205335_cldn4/scripts/01_seurat_harmony_integrate.R
```

Primary engine: **R + Seurat + Harmony**. If Seurat cannot be installed, stop
(do not fall back to a Python-only primary). Python is used only to peel the
GEO gzip-of-gzip RDS.
