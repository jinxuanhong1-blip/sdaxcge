# GEO 2015–2016 leftover: lung × PD-1/PD-L1 × TACSTD2/CLDN4

Honest slice of NCBI GEO Series (GSE) with publication date **2015-01-01 to 2016-12-31**.
Accessions are whatever `esearch` returns. None are invented.

## Query

`db=gds`, `gse[Entry Type]`, `2015/01/01:2016/12/31[PDAT]`, plus

`lung` AND (`PD-1` OR `PD1` OR `PD-L1` OR `PDL1` OR `PDCD1` OR `CD274` OR `B7-H1` OR `programmed death` OR `programmed cell death`).

Script: `scripts/geo_2015_2016_lung_pd1_pdl1.py`.

## What is leftover vs already-easy

First pass already pulled array series matrices whose GPL tables have a gene-symbol column.
This slice finishes the leftovers from that same NCBI hit list:

| Accession | Why leftover | TACSTD2 / CLDN4 |
|---|---|---|
| GSE84789 | featureCounts Geneid is Gencode v19 `ENSG*.version` | mapped via Ensembl REST `xrefs/symbol` (not guessed) |
| GSE76356 | per-sample tables use `ENSMUSG*_chr`; value is `Unique_RPKM` | Tacstd2 present; Cldn4 **absent from the file** |
| GSE84797 | NanoString GPL19965 (`#ID = Official Symbol`); PanCancer Immune panel | **neither gene is on the panel** |
| GSE81258 | SuperSeries, no series-level matrix | RNA-seq SubSeries **GSE81257** named in `Series_relation` (DESeq table, not counts) |

Array / FPKM series that already had symbols (GSE43453, GSE57133, GSE65041, GSE67501, GSE78220) are re-extracted the same way so the folder is complete.

## Do not over-claim

`lung_pd_flags.tsv` only checks title+summary text. Several NCBI hits are **not lung tumors**:

- GSE65041, GSE78220: melanoma PD-1/PD-L1
- GSE67501: RCC anti-PD-1
- GSE76356: macrophage / DICER + PD-1, not a lung cohort
- GSE57133: mouse lung, but title/summary do not mention PD-1/PD-L1
- GSE81258 / GSE81257: Nfib metastasis RNA-seq (KP lung model). Tacstd2 `baseMean` is 0; Cldn4 is present as DESeq stats only

True lung × PD-1/PD-L1 expression leftovers in this window are **GSE84789** (NSCLC B-cell scRNA-seq counts) and **GSE84797** (NSCLC NanoString immune panel, genes not measured).

Ensembl IDs used for mapping are written to `ensembl_id_lookup.tsv` as returned by Ensembl REST. If a gene is missing from a platform or supplementary file, the status line says so. No imputed values.

## Outputs

- `series_metadata.tsv` — NCBI esummary for every hit
- `lung_pd_flags.tsv` — lung / PD term flags from title+summary
- `ensembl_id_lookup.tsv` — live Ensembl IDs
- `processing_status.tsv` — what was done per GSE
- `tacstd2_cldn4_summary.tsv` — per-series numeric summary
- `GSE*_TACSTD2_CLDN4_expression.tsv` — extracted values (empty table = gene not on platform)
