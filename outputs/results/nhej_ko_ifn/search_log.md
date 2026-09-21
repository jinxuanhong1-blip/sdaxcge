# GEO search log

Date: 2026-09-21. Database: NCBI GEO DataSets (`db=gds`) through E-utilities.

## Queries

Perturbation block used with each gene:

`knockout OR knock-out OR knock out OR CRISPR OR gene deletion OR deficient OR deficiency OR null OR sgRNA OR loss-of-function OR loss of function OR deleted`

| Query | Series count |
|---|---:|
| LIG4 / Lig4 / DNA ligase IV / DNA ligase 4, plus perturbation, GSE | 44 |
| XRCC4 / Xrcc4, plus perturbation, GSE | 27 |
| PRKDC / Prkdc / DNA-PKcs / DNA PKcs, plus perturbation, GSE | 38 |
| TP53BP1 / Trp53bp1 / 53BP1, plus perturbation, GSE | 50 |
| LIG4 or ligase in the title, GSE | 37 |
| XRCC4 in the title, GSE | 17 |
| PRKDC / DNA-PKcs / DNA-PK in the title, GSE | 57 |
| TP53BP1 / 53BP1 in the title, GSE | 53 |
| Union | 217 |

Follow-up queries that recovered series the broad text match had under-called because the genotype was only in the design (GSE135274 XRCC4(-/-), GSE145148 XRCC4KO):

| Query | Count |
|---|---:|
| XRCC4 plus knockout/CRISPR/deficient, expression profiling by high-throughput sequencing | 6 |
| LIG4 knockout / Lig4-/- / LIG4 KO, GSE | 6 |
| DNA-PK(cs) knockout / PRKDC knockout / DNA-PKcs-deficient, GSE | 4 |
| 53BP1 knockout / 53BP1-/- / 53BP1 KO, GSE | 15 |

Among the 217 series, 73 were expression profiling by high-throughput sequencing (including mixed experiment types). Titles, summaries, and sample titles were screened for a genetic knockout of LIG4, XRCC4, PRKDC, or TP53BP1. Inhibitor, siRNA, and shRNA experiments were recorded and not scored as knockouts.

## Eligibility

A row enters the primary cancer-line table only when all of these hold:

- processed RNA-seq matrix on GEO
- human or mouse
- cultured cancer cell line
- engineered knockout of one of the four genes, not a knockdown or a drug
- a matched control in the same series
- for the primary row, no drug and no hypoxia

Tumor RNA-seq of a knockout cancer line is scored separately because immune infiltrate can create an interferon signature. Non-malignant lines are scored separately.

## Matrices used

| Accession | File | Quantification |
|---|---|---|
| GSE135274 | `GSE135274_all_sample_human_RNA_max.txt.gz` | counts; exact-match rows only; technical column pairs summed |
| GSE154443 | `GSE154443_gene.fpkm.matrix.xlsx` | FPKM, one column per genotype |
| GSE285698 | `GSE285698_raw_counts.txt.gz` | counts; library ids from GEO `Sample_description` |
| GSE84986 | `GSE84986_RAW.tar` | gene counts; Ensembl symbols from Ensembl REST |
| GSE145148 | `GSE145148_count_matrix.csv.gz` | counts |
| GSE280049 | `GSE280049_53BP1_KO_RNA-seq.xlsx` | integer counts |
| GSE237615 | four `Raw_Counts_*.xls.txt.gz` files | counts, joined on gene symbol |

Ensembl symbol lookup was `https://rest.ensembl.org/lookup/symbol/{species}` on 2026-09-21.
