# Harvested tables — existing public scores / matrices

Additive CLDN4-only slice. No TACSTD2 gate. No dual-high.

| File | n | Source branch | Source path | What is used |
|---|---:|---|---|---|
| `GSE285029_sample_scores.csv` | 234 | `cursor/a11-gse285029-pdl1-4243` | `results/w200/A11_GSE285029/sample_scores.csv` | CLDN4, CD274, IFN_compact, IFNG, HALLMARK_IFNG |
| `GSE218989_per_patient.tsv` | 355 | `cursor/gse218989-cldn4-ici-a090` | `methods/gse218989_cldn4_ici/tables/per_patient.tsv` | CLDN4, IFNG, Ayers6, GEO response |
| `../data/GSE218989_gene_extract.tsv` | 11 genes × 355 | GEO `GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz` | slim extract only | CD274 + IFN-compact genes (not on the harvested table) |
| `GSE126044_per_sample_scores.csv` | 16 | `cursor/w200-b4-gse126044-recompute-0fc6` | `results/w200/B4_GSE126044/tables/per_sample_scores.csv` | CLDN4 log2 CPM, response, sample type |
| `GSE126044_clinical.csv` | 16 | same | `results/w200/B4_GSE126044/data/GSE126044_clinical.csv` | labels |
| `../data/GSE126044_counts.txt.gz` | 18,747 × 16 | same | `results/w200/B4_GSE126044/data/GSE126044_counts.txt.gz` | CD274 + IFN-compact on log2 CPM |
| `GSE166449_sample_level_expression.csv` | 22 | `cursor/gse166449-pembro-analysis-731b` | `results/w200/GSE166449/tables/sample_level_expression.csv` | CLDN4 (deposited scale), response |
| `../data/GSE166449_Raw_gene_TPM_matrix.txt.gz` | 20,345 × 22 | same | `results/w200/GSE166449/data/GSE166449_Raw_gene_TPM_matrix.txt.gz` | CD274 + IFN-compact; scale is the deposited matrix (used as-is) |

**Not invented:** IFN-compact gene list is the A11 a-priori set (`IFNG STAT1 IRF1 CXCL9 CXCL10 CXCL11 IDO1 GBP1`). Ayers6 on GSE218989 is the already-stored column. Hallmark IFN-γ is used only where it already exists (GSE285029).

**Not a dual-high analysis:** TACSTD2 is never a gate, never a co-filter, never a quartile partner.
