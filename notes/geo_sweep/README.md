# Remaining GEO lung series sweep (TACSTD2 / TROP2 / CLDN4 / claudin-4 / ICI)

Exhaustive search of NCBI GEO DataSets (`db=gds`) for **OPEN human/mouse lung Series (GSE)** whose record text mentions any of:

`TACSTD2`, `TROP2`, `CLDN4`, `claudin-4`, `ICI`

Terms were **quoted** so the search matches literal text (title / summary / sample fields), not platform gene-annotation presence. An unquoted `TACSTD2` query expands to `"tacstd2 protein, human"` and returns hundreds of array platforms that merely list the gene.

Filter: `lung AND (Homo sapiens[ORGN] OR Mus musculus[ORGN]) AND gse[ETYP]`.

## Classic set excluded (not in the catalog)

`GSE126044`, `GSE135222`, `GSE136961`, `GSE166449`, `GSE93157`, `GSE207422`, `GSE205335`

Of these, only `GSE205335` was returned by the quoted-term search and was dropped.

## Outputs

| Path | What it is |
| --- | --- |
| `notes/geo_sweep/catalog.tsv` | Every verified extra accession (80 GSE). No invented accessions. |
| `results/geo_sweep/marker_stats.tsv` | **The single results table.** One row per (accession, gene) for TACSTD2 and CLDN4, plus four extra rows for the GSE50927 DE contrasts. |
| `scripts/geo_sweep/` | Reproducible search / probe / download / extract / table-build scripts. |

Raw series matrices and supplementary blobs live under `results/geo_sweep/matrices/` and `results/geo_sweep/supp/` (git-ignored).

## Catalog counts

- 80 extra OPEN series, all with a processed series matrix < 2 GB (downloaded).
- 63 ICI-term hits; 16 KD/KO-classified (regex on title/summary/sample text).
- 8 series carry an embedded expression table; the rest are metadata-only matrices (typical for RNA-seq).

## Marker analysis (no fabricated statistics)

`marker_stats.tsv` numeric fields are filled only when a real processed table contained the gene.

- **Embedded matrices** (Affymetrix / Illumina / Clariom): probe → symbol via GEO platform annotation. Group means use GEO sample titles. Mann–Whitney is reported only when exactly two groups have n≥2.
- **Supplementary gene-level tables**: symbol, Ensembl gene, Ensembl transcript, or Entrez id. Mouse *Cldn4* is `ENSMUSG00000047501` (verified via Ensembl REST 2026-08-16). `ENSMUSG00000024959` is *Bad* and was not used.
- **GSE50927** (Cldn4 KO lungs ± VILI): all four author edgeR contrasts are reported (logFC / P / FDR from the deposited CSVs).
- **GSE182261** (NIK knockdown, ICI+KD): raw counts, shGFP (n=2) vs shNIK (n=4).
- **GSE169688** (ICI, LLC1 scRNA): marker rows stream-extracted from the 682 MB 3-day matrix (under 2 GB).
- **GSE235122** (ICI, melanoma methylation): 320 MB processed matrix streamed; TACSTD2/CLDN4 rows absent.
- **GSE152590** (ICI, CD8 T-cell TPM): gene list has neither TACSTD2 nor CLDN4.

Blank numeric cells mean the value was not computed. `marker_status` records why: `found`, `not_in_table`, `not_on_platform`, `no_processed_table`, `too_large_or_absent`.

Usable new ICI / KD sets with stats in the table include GSE182261, GSE238006, GSE266364, GSE274960, GSE285029, GSE288083, GSE317011, GSE330941, GSE333285, GSE33348, GSE48443, and GSE50927 (plus ICI-only series that had a gene-level table).
