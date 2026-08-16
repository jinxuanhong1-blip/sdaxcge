# Data access (this slice)

Public files used here are downloaded by `scripts/opus_tls/00_download.sh`
into `$OPUS_TLS_DATA` (default `/tmp/opus_tls_data`). They are **not**
committed. Every file is < 200 MB except the two TCGA STAR TPM matrices
(~160–168 MB each) and the unused GSE207422 scRNA UMI matrix (184 MB,
downloaded then ignored).

## Used (open)

| Cohort | Accession / hub | Why | n used | Size |
|---|---|---|---|---|
| TCGA-LUAD | UCSC Xena GDC hub `TCGA-LUAD.star_tpm.tsv.gz` | LUAD atlas + ABSOLUTE purity | 528 primary tumours | 168 MB |
| TCGA-LUSC | same, `TCGA-LUSC.star_tpm.tsv.gz` | LUSC atlas + ABSOLUTE purity | 501 primary tumours | 160 MB |
| PanCanAtlas purity | GDC `4f277128-…` ABSOLUTE mastercalls | tumour-cell fraction (DNA) | — | 0.9 MB |
| PanCanAtlas leukocytes | GDC `6f75c9d7-…` | methylation leukocyte fraction | — | 0.6 MB |
| PanCanAtlas CIBERSORT | GDC `b3df502e-…` | orthogonal B / plasma / Tfh / CD8 | — | 3.8 MB |
| GSE72094 | GEO series matrix + GPL15048 | LUAD microarray atlas, OS | 442 | 106 MB |
| GSE81089 | GEO FPKM | Uppsala mixed NSCLC atlas, OS | 197 tumours | 24 MB |
| GSE135222 | GEO TPM | advanced NSCLC, anti-PD-1/PD-L1, PFS | 27 | 1.6 MB |
| GSE126044 | GEO counts | NSCLC anti-PD-1, RECIST | 16 | 0.6 MB |
| GSE207422 | GEO log2 TPM + xlsx | neoadjuvant PD-1 + chemo, MPR | 24 bulk (pre+post) | 5.6 MB |
| GSE190265 | GEO TPM (France3) | NSCLC biopsies | 26 | 4.1 MB |

## Listed only — not downloaded (EGA / controlled)

These are the lung ICI RNA-seq sets that would be the most relevant
replication of Bessede *et al.* (POPLAR/OAK-style atezolizumab) and of
the Patil TLS paper. They sit behind EGA DAC and cannot be fetched in
this environment.

| Dataset | EGA / dbGaP | Why it matters | Approx. size |
|---|---|---|---|
| **Patil et al. 2022 TLS in IMpower150 / IMpower010 / OAK / POPLAR** | EGAD00001008548 (OAK RNA), related IMpower EGA studies | histologically scored TLS + atezolizumab RNA in NSCLC | multi-GB, controlled |
| **OAK (NCT02008227) atezolizumab vs docetaxel** | EGAD00001002517 / EGAS00001005013 (study-dependent) | largest 2L NSCLC ICI RNA-seq with OS/PFS | multi-GB, controlled |
| **POPLAR (NCT01903993)** | EGA (Roche/Genentech) | companion randomised ICI RNA | controlled |
| **IMpower150 / 130 / 132** | EGA (Roche) | 1L chemo-IO RNA | controlled |
| **TRACERx / PEACE** | EGAD00001004585 and related | multi-region LUAD, TLS spatial | controlled |
| **Stand Up To Cancer / SU2C-MARK** | dbGaP phs002916 / related | metastatic NSCLC ICI RNA | controlled |
| **MSK-IMPACT NSCLC ICI** | dbGaP (Hellmann / Rizvi) | WES+RNA subset, not uniformly open | controlled |

If a later slice obtains EGA access, the same signature / partial-Spearman
code in `scripts/opus_tls/` can be pointed at those matrices without
changing the statistics.

## Inspected and not used as a primary test

| Dataset | Reason |
|---|---|
| GSE207422 scRNA UMI (184 MB) | single-cell; this slice is bulk signatures |
| GSE93157 NanoString 730-gene | panel lacks a full TLS / plasma set |
| GSE166449 | not a lung-ICI or lung-atlas cohort (PRECISION oncology mix) |
| GSE72094 RAW.tar | redundant with the series matrix |
