# GEO verification — human ICI scRNA-seq slice

Verified 2026-08-16 directly against NCBI GEO (accession viewer + family SOFT
files downloaded from `ftp.ncbi.nlm.nih.gov`).

## GSE207422 — neoadjuvant anti-PD-1 + chemotherapy NSCLC (MPR endpoint)

- Title: "Tumor microenvironment remodeling after neoadjuvant immunotherapy in
  non-small cell lung cancer revealed by single-cell RNA sequencing"
- Status: Public on Mar 21 2023; PubMed 36869384 (and 41007733).
- Platform GPL24676 (Illumina NovaSeq 6000, Homo sapiens); 39 GSM samples
  (GSM6287372–GSM6287410; scRNA + bulk).
- Contact: Peng Zhang, Shanghai Pulmonary Hospital. Raw data in GSA HRA001033
  (not on GEO); processed data on the Series record.
- Processed supplementary files used (both well under the 2 GB limit):
  - `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` — 175.5 MB, dense genes ×
    cells UMI count matrix, 92,330 cells, barcodes prefixed by sample
    (`BD_immune01_...` … `BD_immune15_...`).
  - `GSE207422_NSCLC_scRNAseq_metadata.xlsx` — 11.1 KB, 15 scRNA samples with
    Patient, Resource (pre-treatment biopsy n=3 / post-treatment surgery
    n=12), PD-1 antibody, chemo, Pathologic Response (MPR / NMPR / pCR / NE),
    Residual Tumor fraction, RECIST.
  - Bulk RNA-seq files exist on the Series but were not needed for this
    scRNA-only slice.
- No per-cell cell-type annotation is provided on GEO, so cell lineages were
  assigned in-house from a marker panel (see WRITEUP methods and
  `results/fable_scrna_ici/gse207422_lineage_marker_qc.tsv`).

## GSE205335 — lung cancer under immune checkpoint inhibitors (RECIST endpoint)

- Title: "Single-cell transcriptome profiles of tumor tissues from lung cancer
  patients receiving immune checkpoint inhibitors."
- Status: Public on Nov 05 2024. Contact: Hae-Ock Lee, The Catholic University
  of Korea (contributors Myung-Ju Ahn, Hae-Ock Lee). Raw data in EGA
  EGAD00001008703 (not on GEO).
- Platforms GPL16791 / GPL24676; 33 GSM samples (GSM6210624–GSM6210656) from
  26 patients: metastatic LN, liver, effusion, lung/bronchus tumors plus a few
  normal tissues (Normal Lung / Normal LN / Normal Brain).
- Per-sample characteristics in family SOFT: patient, tissue, tumor stage,
  cancer subtype (ADC / SQ / SCLC / NUT), RECIST (PR / SD / PD / NE), 10x
  platform (3' or 5').
- Processed supplementary files used (both under 2 GB):
  - `GSE205335_Lung_IO_UMI_matrix.rds.gz` — 499.5 MB download; NOTE: the file
    is double-gzipped (a gzip of an internally gzip-compressed RDS); one
    `zcat` pass yields a readable RDS containing a `dgCMatrix` of
    33,714 genes × 96,505 cells.
  - `GSE205335_Lung_IO_CellIdentity.txt.gz` — 718.7 KB; per-cell annotation
    (96,505 cells): `orig.ident` (sample), `lineage.total`, `lineage.sub`
    (incl. "Malignant cells", "CD8+ T cells"), fine `celltype`
    (incl. "CD8+ TEX").
- GSM ↔ `orig.ident` mapping was reconstructed from sample titles + platform
  (e.g. "P1006 EBUS_06" + 3' → `EBUS-06-3P`) and verified to match 33/33 with
  no unmatched `orig.ident`
  (`results/fable_scrna_ici/gse205335_sample_table.tsv`).

## Size policy

All downloads were processed matrices/annotations only; the largest file was
523.7 MB on disk (GSE205335 RDS). Nothing over 2 GB was downloaded; no raw
FASTQ/BAM touched.
