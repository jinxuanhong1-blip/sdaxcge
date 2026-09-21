# Methods

## Catalogs

Live queries on 2026-09-21. Script: `scripts/00_inventory.py`. Log: `results/search_log.json`, `results/inventory_summary.json`.

PRIDE project search (`/pride/ws/archive/v2/search/projects`):

| Keyword | Projects |
| --- | ---: |
| CLDN4 | 0 |
| claudin-4 | 0 |
| Cldn4 | 0 |
| claudin | 30 |
| claudin knockdown | 0 |
| claudin interactome | 0 |
| BioID claudin | 0 |

Each of the 30 projects was re-read from `/projects/{accession}`. A project was marked in scope when the title or description contained CLDN4 / claudin-4 together with a knockdown or interactome term, or when a deposited filename contained `Cldn4`. That filename rule is what recovers [PXD031094](https://www.ebi.ac.uk/pride/archive/projects/PXD031094): the project text says “pan-claudin” and “claudins”, and the CLDN4 CoIP is the SEARCH file `proteinGroups_Cldn4.txt`.

MassIVE ProXI `keywords=CLDN4` does not filter. The first compact record is MSV000065795. The script therefore downloaded `https://massive.ucsd.edu/ProteoSAFe/datasets_json.jsp` (20,228 datasets) and filtered titles and descriptions for `claudin` or `CLDN` plus a digit. One dataset mentions claudins: MSV000095998 (CLDN10B). None mention CLDN4 or claudin-4.

## Files used

| File | Origin | Check |
| --- | --- | --- |
| `data/proteinGroups_Cldn4.txt` | PRIDE SEARCH for PXD031094 | SHA-1 `4dbde6c0d827105412ef1f34cbb969b8dd8acb36`, equal to the live PRIDE file checksum. 3,950,383 bytes. 2,788 protein groups. |
| `data/plos2015_s2_all_proteins.xlsx` | PLOS ONE 2015 S2 Table, DOI [10.1371/journal.pone.0117074.s006](https://doi.org/10.1371/journal.pone.0117074.s006) | Complete BioID identifications. BL-Cldn4 names are column M. |
| `data/plos2015_s3_enriched.xlsx` | PLOS ONE 2015 S3 Table, DOI [10.1371/journal.pone.0117074.s007](https://doi.org/10.1371/journal.pone.0117074.s007) | Three-fold enriched subset. Shared strings scanned. |
| `data/crc24_suppst1_bioid.docx` | Figshare file 53292048, MD5 `70bcbba8b14311102973f11eda6cc41d` | Villagómez et al. 2024 Supplementary Table 1, “BioID Data - Claudin-4 proximal proteins”. |

`ftp.pride.ebi.ac.uk` returned a TLS unexpected EOF from this environment, and FTP timed out. `proteinGroups_Cldn4.txt` was taken from the earlier public checkout of the same SEARCH file after the SHA-1 matched the live API checksum. `proteinGroups_PRISMA.txt` (66,679,807 bytes, PRIDE SHA-1 `770554ea6737df2a850b08fa2404d0f5be0f75b6`) was not retrieved and was not scanned. RAW files were not downloaded. `scripts/02_download.py` writes `results/download_log.json`.

## Partner scan

`scripts/01_scan_partners.py` searches gene symbols and protein-name phrases:

`PRKDC`, `XRCC4`, `XRCC5`, `XRCC6`, `STING1`, `TMEM173`, `CGAS`, `MB21D1`, DNA-PKcs, Ku70, Ku80, Ku86, “DNA-dependent protein kinase”, “X-ray repair cross”, “stimulator of interferon genes”.

The token `Ku` alone is not used. In the canine FASTA it matches Kunitz-type SPINT2.

For PXD031094 the fields searched are `Gene names`, `Protein names`, `Fasta headers`, and `Majority protein IDs`. LFQ columns are `LFQ intensity Cldn4_1..4` and `LFQ intensity GFPctrl_1..4`.

The 2024 BioID docx is parsed as repeating rows of gene, control sample 1, CLDN4 sample 2, and log2(CLDN4/control). That table is one control column and one CLDN4 column.

## Out of scope

CosMx, concordant-4 scRNA, GSE137244, TCGA keratin correction, and TISMO are unchanged. CPTAC abundance tables and ICI response proteomes are a different question. CLDN6-deficient PXD066158 and CLDN10B MSV000095998 are different genes and were not scanned for DNA-PKcs.
