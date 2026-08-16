# 01 · Cohort and data QC

Release: **DepMap Public 24Q4** (figshare DOI `10.6084/m9.figshare.27993248`, retrieved 2026-08-16T19:30:56Z). All files are open release files; md5 checksums are recorded in `data/opus_depmap/MANIFEST.json`.

## Cohort definition

* Lung lineage = `Model.csv` rows with `OncotreeLineage == "Lung"`.
* Non-cancerous lung models (`OncotreePrimaryDisease == "Non-Cancerous"`) are dropped.
* `LungGroup` collapses Oncotree annotations into NSCLC / SCLC / Other lung (the latter is mostly SMARCA4-deficient undifferentiated thoracic tumours and non-SCLC neuroendocrine models).

## Counts

| metric | value |
| --- | --- |
| release | DepMap Public 24Q4 |
| figshare_doi | 10.6084/m9.figshare.27993248 |
| models_total | 2105 |
| lung_lineage_cancer_models | 254 |
| lung_with_CRISPR | 126 |
| lung_with_expression | 208 |
| lung_with_CRISPR_and_expression | 123 |
| all_models_with_CRISPR | 1178 |
| genes_in_CRISPR_matrix | 17916 |
| all_models_with_expression | 1673 |
| genes_in_expression_matrix | 19193 |
| lung_NSCLC_models | 164 |
| lung_NSCLC_with_CRISPR | 98 |
| lung_NSCLC_with_expression | 143 |
| lung_SCLC_models | 81 |
| lung_SCLC_with_CRISPR | 25 |
| lung_SCLC_with_expression | 59 |
| lung_Other lung_models | 9 |
| lung_Other lung_with_CRISPR | 3 |
| lung_Other lung_with_expression | 6 |

## Lung subtypes with CRISPR screens

| OncotreeSubtype | n with CRISPR |
| --- | --- |
| Lung Adenocarcinoma | 53 |
| Small Cell Lung Cancer | 25 |
| Lung Squamous Cell Carcinoma | 21 |
| Large Cell Lung Carcinoma | 10 |
| Non-Small Cell Lung Cancer | 6 |
| Giant Cell Carcinoma of the Lung | 3 |
| SMARCA4-deficient undifferentiated tumor | 2 |
| Lung Adenosquamous Carcinoma | 2 |
| NUT Carcinoma of the Lung | 1 |
| Lung Carcinoid | 1 |
| Mucoepidermoid Carcinoma of the Lung | 1 |
| Poorly Differentiated Non-Small Cell Lung Cancer | 1 |

## Measurement notes

* `CRISPRGeneEffect.csv` is the integrated Chronos gene effect matrix, scaled so that the median common essential is -1 and the median non-essential is 0.
* `CRISPRGeneDependency.csv` gives the posterior probability that a model is dependent on a gene; DepMap's convention is that probability > 0.5 counts as a dependent line.
* Expression is `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, i.e. log2(TPM+1) from the GTEx-style RNA-seq pipeline, model-level.
* Duplicate gene symbols (same symbol, different Entrez id) are de-duplicated by keeping the first occurrence; neither TACSTD2 nor CLDN4 is affected.
