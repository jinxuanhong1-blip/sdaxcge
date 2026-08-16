# Search log — C4_more

## Question

C4 public analog: only mouse lung GSE50927 opens IFN; cancer lines close it.
Hunt **any additional public CLDN4 KD/KO transcriptome** that is **not**
GSE50927, GSE207704, GSE22493, or GSE245459. If none, say none.

## Eligibility

Included only if all of the following are true:

1. Direct CLDN4/Cldn4 loss: knockdown, knockout, siRNA, shRNA, or CRISPR.
2. Transcriptome-scale assay (RNA-seq or expression microarray), not qPCR-only.
3. Public accession with a downloadable gene-level matrix or author DE table.
4. Accession is not one of the four excluded IDs.

Excluded:

- Observational CLDN4 expression, overexpression-only, or ceRNA/tumor-normal.
- Ligand / peptide downregulation (C-CPE) without genetic KD/KO.
- Perturbation of another claudin, or complete Cldn-null (all family members).
- TACSTD2/TROP2 KD/KO (including GSE245459 and GSE334497).
- Papers that mention RNA-seq but deposit no public accession or matrix.
- Methylation, ChIP, ATAC, 4C, or other non-expression assays.

## What was queried (2026-08-16)

NCBI GEO DataSets (`gds`), GSE entry type, via E-utilities
(`scripts/rework_C4_more/search_public.py` → `live_search.json`):

| Query | GSE-level hits |
|---|---|
| `CLDN4[All Fields] AND gse[Entry Type]` | 21 |
| `Cldn4[All Fields] AND gse[Entry Type]` | 21 (same UID set) |
| `"claudin 4"[All Fields] AND gse[Entry Type]` | 10 |
| `claudin-4[All Fields] AND gse[Entry Type]` | 10 |
| `CLDN4` + knockdown/knockout/siRNA/shRNA/CRISPR/silencing + GSE | 4 |
| `Cldn4` + knockdown/knockout/KO/siRNA/shRNA/CRISPR + GSE | 5 |
| `Cldn-null AND gse[Entry Type]` | 1 (GSE274940) |

Union of those GSE series, plus explicit lookups of GSE22421, GSE274940,
GSE245459, and GSE334497, gave **29 unique GSE records**. Each title,
summary, assay type, and sample titles were read from `esummary`.

SRA text searches for `CLDN4 knockout RNA-seq`, `Cldn4 knockout transcriptome`,
`CLDN4 knockdown RNA-seq`, and `CLDN4 shRNA RNA-seq` each returned **0** runs.

PubMed combinations of CLDN4/Cldn4/claudin-4 with knockout/CRISPR/silencing
and RNA-seq/transcriptome were used to catch papers that sequenced but did
not index the gene name in the GEO record. Two papers report CLDN4-loss
RNA-seq without a public accession (PMID 41016339, PMID 40892111).

ArrayExpress / BioStudies:

- `E-MTAB AND CLDN4` → 0 hits.
- CLDN4 + knockdown/knockout hits are E-GEOD mirrors of GSE22493 / GSE50927
  plus literature records, not a new unique study.

OmicsDI search for CLDN4 knockout/knockdown returned the same known GEO
mirrors (E-GEOD-22493, E-GEOD-50927) and literature, not a new series.

This is a systematic indexed-resource hunt. It cannot prove that an
unindexed lab website, controlled-access archive, or unreleased matrix
does not exist. No accession in this folder was invented.
