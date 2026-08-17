# Harvested tables — existing public matrices only

No new GEO/Xena/CPTAC download. Each file is a copy of a sample-level
table already computed on a public matrix in this repo.

| File | n (rows) | Source branch | Source path |
|---|---:|---|---|
| GSE218989_per_patient.tsv | 355 | `cursor/gse218989-cldn4-ici-a090` | `methods/gse218989_cldn4_ici/tables/per_patient.tsv` |
| GSE285029_sample_scores.csv | 234 | `cursor/a11-gse285029-pdl1-4243` | `results/w200/A11_GSE285029/sample_scores.csv` |
| TCGA_LUAD_xena_star_tpm.tsv | 502 | `cursor/xena-luad-tacstd2-immune-purity-ddb5` | `results/w200/Xena_LUAD/tables/analysis_table_star_tpm.tsv` |
| TCGA_LUSC_xena_sample_table.tsv | 502 | `cursor/xena-lusc-tacstd2-immune-purity-c1cc` | `results/w200/Xena_LUSC/sample_table.tsv` |
| CPTAC_LUAD_sample_table.tsv | 110 | `cursor/cptac-luad-trop2-cldn4-6a8b` | `results/grok_cptac_luad/tables/sample_table.tsv` |
| CPTAC_LSCC_sample_scores.tsv | 108 | `cursor/cptac-lusc-cldn4-protein-b9ed` | `methods/cptac_lusc_cldn4_protein/tables/sample_scores.tsv` |
| CPTAC_LSCC_sample_level_features.tsv | 108 | `cursor/cptac-lscc-protein-rna-immune-4215` | `results/grok_cptac_lscc/tables/sample_level_features.tsv` |
| GSE68465_samples.tsv | 443 | `cursor/a1-luad-wave2-6484` | `results/rework/A1_extra_luad/samples_GSE68465.tsv` |
| GSE31210_samples.tsv | 226 | `cursor/a1-luad-wave2-6484` | `results/rework/A1_extra_luad/samples_GSE31210.tsv` |
| OncoSG_sample_table.tsv | 169 | `cursor/oncosg-a1-cd8-gep-c2f6` | `results/rework/OncoSG_A1/sample_table.tsv` |

**OncoSG CLDN4:** the public cBioPortal z-score matrix used in that table
does not carry CLDN4. TACSTD2-only recut.

**Not invented:** IFN / MHC-I / GEP18 / OS are used only when the harvested
table already has that column. No signature was recomputed from raw counts.
