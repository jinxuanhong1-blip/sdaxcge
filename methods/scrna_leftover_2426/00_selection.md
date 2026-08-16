# Selection log — leftover 2024–2026 lung ICI / neoadjuvant scRNA

Live GEO / ArrayExpress metadata pulled 2026-08-16. No accessions invented. Core set excluded: GSE207422, GSE205335, GSE241934, GSE291670, GSE243013, GSE253013, GSE131907.

## Candidate list (user)

| Accession | Keep / skip | Why |
|-----------|-------------|-----|
| GSE274584 | inventory only | EV RNA count table. Not scRNA. |
| GSE274588 | skip scoring | PBMC 10x. No epithelium. |
| GSE274595 | **score** | Tumor tissue + nuclei MTX, 8 patients, 181 MB. No ICI labels. |
| GSE302113 | skip scoring | scATAC + mtDNA (36 GB). Not RNA. |
| GSE267108 | **score** | 8 LUAD tumors, processed MTX 249 MB. PD-L1 groups inferred from published cell counts. Treatment-naïve. |
| GSE176021 | skip scoring | GEO processed = lymphocytes / CD3+CD8 RDS. Raw FASTQ is EGA. Has MPR but no epithelium. |
| GSE176022 | skip scoring | Bulk TCR / MANAFEST. Not scRNA expression. |
| GSE186446 | bulk proxy only | scRNA is T cells (3 ICB pts). GEO file `fcount_aggr.txt.gz` is bulk. |
| GSE337519 | **descriptive** | n=1 neoadjuvant chemo-IO. 34 MB MTX. |
| GSE308745 | skip scoring | PBMC cohort. Skip if no epithelium. |
| E-MTAB-13526 | **score extra n** | Cvejic/De Zuani atlas. Not ICI. Tumor h5ad 58 GB streamed remotely. |

## What was downloaded

GEO processed / RAW MTX for GSE274595, GSE267108, GSE337519. GSE176021 annotation RDS (confirm T-only). GSE186446 fcount. GSE274584 EV table. E-MTAB tumor h5ad was **not** fully downloaded; obs + two gene columns were streamed over HTTPS.

## Not used

EGA (GSE176021 FASTQ), dbGaP, GSE185206 SuperSeries RAW.tar (TCR/scRNA T cells; not required once T-only is established), GSE274588 1.1 GB PBMC, GSE302113 36 GB fragments, GSE308745 1.6 GB PBMC.
