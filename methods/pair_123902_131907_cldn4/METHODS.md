# Methods — pair GSE123902 + GSE131907, CLDN4 only

Patient is the unit. CLDN4 is the only scoring gene. TACSTD2 is loaded only as an audit column and is not used to define groups.

## Matrices

- **GSE123902:** GEO `GSE123902_RAW.tar` (90.4 MB). Per-sample SEQC dense UMI, cells × genes. Full-library nUMI is the row sum of all ~16k genes; CP10k uses that denominator. The 36.5 GB author H5 is skipped.
- **GSE131907:** GEO raw UMI `txt.gz` (389.8 MB) + `GSE131907_Lung_Cancer_cell_annotation.txt.gz`. The 2.9 GB log2TPM text is skipped.

## Compartments

- **GSE123902:** lineage = argmax of mean log1p(CP10k) among epithelial (EPCAM, KRT8, KRT18, KRT19, KRT7), T/NK (CD3D, CD3E, CD8A, NKG7, GNLY, KLRD1), myeloid (LYZ, CD14, CSF1R, AIF1), B (MS4A1, CD79A). Winner must be ≥0.12 and ≥1.15× the runner-up. Malignant = epithelial in PRIMARY_TUMOUR or METASTASIS. CLDN4 is not in the assigner.
- **GSE131907:** malignant = author `Cell_subtype` ∈ {Malignant cells, tS1, tS2, tS3}. T/NK = author `Cell_type` ∈ {T lymphocytes, NK cells}. Tumor-origin samples only (tLung, tL/B, mLN, mBrain). One row per patient.

## Test

Spearman of malignant mean CLDN4 vs T/NK fraction. Eligibility: ≥20 malignant and ≥20 T/NK. Combo: DerSimonian–Laird on Fisher-z. Difference rule (CellChat trigger): opposite sign with both |ρ|≥0.10, or |Δρ|≥0.25, or one arm p<0.05 and |Δρ|≥0.20, or I²≥50%.
