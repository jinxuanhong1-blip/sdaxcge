# NSCLC ICI scRNA / spatial — GEO survey (this slice)

Question: in **malignant / epithelial cells**, is **TACSTD2** (and **CLDN4**) associated with T/NK (or CD8 / TLS) infiltration and ICI **response**?

Only processed matrices. Files >2 GB skipped. T-sorted-only series skipped (GSE176022, GSE99254). LuCA h5ad skipped. No invented accessions.

## Used

| Accession | Design | Processed file | Size | ICI label |
|---|---|---|---|---|
| **GSE207422** | Neoadjuvant PD-1 + chemo NSCLC (Hu et al. *Genome Med* 2023, PMID 36869384) | `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` + sample `metadata.xlsx` | 175 MB + 11 KB | **MPR / pCR / NMPR** and RECIST. 15 samples (3 pre-tx biopsy, 12 post-tx surgery). No per-cell CopyKAT labels on GEO. |
| **GSE205335** | Palliative lung ICI atlas (Park/Ahn/Lee) | `GSE205335_Lung_IO_UMI_matrix.rds.gz` + `CellIdentity.txt.gz` + SOFT | 499 MB + 719 KB | **RECIST** only (PR/SD/PD/NE). Not MPR. Author `lineage.sub` includes `Malignant cells`. |
| **GSE271689** | GeoMx DSP, ICI, OS (spatial, not single-cell) | `GSE271689_RAW.tar` (35 MB) | 35 MB | OS / ICI. Secondary; not required for the “done” gate. |

## Surveyed, not used as primary ICI-response tests

| Accession | Why not primary |
|---|---|
| GSE131907 | NSCLC atlas (Kim/Lee). Treatment-naive / mixed; **not ICI-response labeled**. Raw UMI txt 390 MB OK; log2TPM txt **2.9 GB skipped**. |
| GSE154826 | LCAM (Leader). MARS-seq amp_batch tarballs; not an ICI-response trial. |

## Skipped by instruction

- GSE176022, GSE99254 — T-cell sorted only.
- Any file >2 GB (GSE131907 normalized log2TPM txt; LuCA h5ad).

## GEO URLs (verified 2026-08-16)

- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/suppl/
