# Residual TACSTD2 / CLDN4 scores (not a tumor score)

`gse243013_residual_tacstd2_cldn4.tsv` is the patient-level detection table
computed from the public CD45+ count matrix
`GSE243013_NSCLC_immune_scRNA_counts.mtx.gz` (see PR
`cursor/gse243013-mpr-analysis-e70e`).

These values are **immune-compartment residual / ambient detection**.
They are **not** tumor-epithelial TACSTD2/CLDN4 and must not be used as E2.
