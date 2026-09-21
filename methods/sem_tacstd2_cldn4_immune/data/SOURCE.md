# Inputs

## Concordant-4 patient table

`concordant4_patient_scores.tsv` is the patient-level table from the concordant-4 ELF3 / TACSTD2 / CLDN4 analysis (branch `cursor/concordant4-elf3-module-tnk-9b4c`, file `methods/concordant4_elf3_regulon/results/tables/patient_scores.tsv`).

- Unit: patient / donor / sample, n = 65.
- Cohorts: GSE123902 (13), GSE131907 (21), GSE205335 (22), GSE189357 (9).
- `pct_TACSTD2` and `pct_CLDN4`: percent of malignant cells with a detected UMI.
- `mean_TACSTD2` and `mean_CLDN4`: mean log1p(UMI) in malignant cells.
- `frac_tnk`: T/NK fraction of all cells in the unit.
- `n_malignant`: malignant-cell count. Every unit has at least 27.

The analysis script checks that the DerSimonian–Laird Spearman of `pct_CLDN4` versus `frac_tnk` equals the locked value −0.5311678045689989.

## TCGA gene table

`tcga_primary01_genes.tsv.gz` is a patient-level extract of UCSC Xena GDC STAR log2(TPM+1) (`gdc-hub.s3.us-east-1.amazonaws.com`, GENCODE v36). Primary solid tumor only (sample-type code 01). Replicate aliquots are averaged. Ensembl IDs and the SHA-256 of the extract are in `tcga_provenance.json`.

Cohorts: LUAD 516, LUSC 501, BRCA 1095, CESC 304, KIRC 533, STAD 412, BLCA 406, PAAD 178.

The analysis script checks the LUAD partial Spearman of TACSTD2 and of CLDN4 versus the CD8 score after KRT8/KRT18/KRT19 against the values published with the earlier keratin analysis (−0.103504 and −0.084776).
