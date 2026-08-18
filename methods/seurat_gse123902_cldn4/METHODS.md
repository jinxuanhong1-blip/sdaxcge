# Methods — Seurat GSE123902 CLDN4-only

ADDITIVE. Human LUAD/NSCLC. CLDN4 only. Patient/donor is the unit.

## Input

Public processed GSE123902 SEQC dense unnormalized UMI CSVs (`GSE123902_RAW.tar`). Cells × genes, one file per sample. Gene universes differ by sample; the union is used and missing genes are 0. The 36.5 GB author annotated H5 is not used.

CSVs are converted to a 10x-style MTX (`matrix.mtx` + `barcodes.tsv` + `features.tsv`). Seurat `ReadMtx` then `CreateSeuratObject` is the object constructor. If Seurat cannot be installed, the analysis stops.

## Lineage (CLDN4 never used)

On Seurat `LogNormalize` (log1p CP10k): four-way marker mean for epithelial (`EPCAM`, `KRT8`, `KRT18`, `KRT19`, `KRT7`), T/NK (`CD3D`, `CD3E`, `CD8A`, `NKG7`, `GNLY`, `KLRD1`), myeloid (`LYZ`, `CD14`, `CSF1R`, `AIF1`), B (`MS4A1`, `CD79A`). Assign the winner if top ≥ 0.12 and top ≥ 1.15 × second. Malignant = epithelial in tumor (primary or metastasis). Matched normal is not malignant.

## Donor tests (mouse-style)

Eligible donor: ≥20 malignant and ≥20 T/NK cells in that donor’s tumor samples.

- Spearman of malignant CLDN4 mean or %pos vs T/NK fraction (Fisher-z 95% CI).
- Spearman of malignant CLDN4 mean vs IFN / MHC-I/APM / TJ. TJ holds CLDN4 out.
- If eligible n ≥ 8: Q4 vs Q1 on malignant CLDN4 mean for IFN / MHC / TJ (Mann–Whitney).

Family score = mean log1p CP10k of locked genes (not Seurat control-pool module scores).

## Not done

No TACSTD2∩CLDN4 dual-high. No GSE148071 merge. No Python-only primary. No ICI / MPR / survival claim.
