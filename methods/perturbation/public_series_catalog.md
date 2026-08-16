# Public CLDN4 / TACSTD2 perturbation series (verified)

**Only accessions returned live by NCBI E-utilities / GEO FTP are listed.
Nothing here is invented. Re-query before citing — GEO changes.**

Queried 2026-08-16 via `templates/geo_query.py` and `esearch/esummary` on `db=gds`.
Suppl file sizes from HTTP `Content-Length` on NCBI GEO FTP.

## RNA-seq series used in this folder

| Accession | Organism | n | Perturbation | Processed suppl | Size | Run |
|---|---|---:|---|---|---:|---|
| **GSE207704** | human | 8 (1/condition) | CLDN4 CRISPR KO vs WT (MCF7, T47D) | `GSE207704_CLDN4_RNAseq.txt.gz` (FPKM) | 1.0 MB | [`example_run/`](example_run/) |
| **GSE245459** | human | 12 (n=3) | shTACSTD2 vs shNC ± cisplatin (SKOV3) | `GSE245459_fpkm.anno.txt.gz` (FPKM) | 14.0 MB | [`example_gse245459/`](example_gse245459/) |
| **GSE334497** | mouse | 10 (n=5) | Trop2/Tacstd2 KO vs WT (4T1 tumors) | `GSE334497_normalized_counts.csv.gz` | 1.2 MB | [`example_gse334497/`](example_gse334497/) |

All three are **public** and **< 2 GB**. None of the processed matrices are raw
integer counts, so none of the example runs claim a DESeq2/edgeR p-value.

## Other verified public series (catalogued, not re-run here)

| Accession | Organism | n | Assay | Why it is listed | Why it was not re-run |
|---|---|---:|---|---|---|
| **GSE22493** | human | 3 | microarray | CLDN4-silencing vs CLDN4-overexpressing SKOV-3 | array (not RNA-seq); `GSE22493_RAW.tar` is CEL files (13.0 MB) needing an Affymetrix pipeline |

## Searches that did **not** yield a dedicated series

These queries were run; they returned either zero Series or only keyword-coincidence
hits (the gene name appears in a sample/description but the perturbation is of
something else). **Absence is reported, not papered over.**

| Query | Dedicated perturbation RNA-seq series? |
|---|---|
| `TROP2 knockdown` | none |
| `TROP2 shRNA` | none |
| `TACSTD2 CRISPR` as a standalone human cell-line KO RNA-seq | none found beyond GSE334497 (mouse in-vivo KO) |
| additional human **CLDN4** KD/KO RNA-seq besides GSE207704 | none found in this search (GSE22493 is array) |

## How to refresh this catalog

```bash
python templates/geo_query.py --suppl \
    "CLDN4 knockout" "CLDN4 knockdown" "CLDN4 shRNA" "CLDN4 CRISPR" \
    "TACSTD2 knockdown" "TACSTD2 shRNA" "TROP2 knockout"
```

Fielded GEO web query (RNA-seq Series only):

```
CLDN4[Description] AND "expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]
TACSTD2[Description] AND "expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]
```
