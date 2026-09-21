# FINDING — CLDN4 knockdown proteomes and claudin-4 partners of DNA-PKcs / Ku / STING

**The open CLDN4 CoIP and the two open claudin-4 BioID tables contain zero rows for PRKDC, XRCC5, XRCC6, STING1, or TMEM173.** No PRIDE or MassIVE project is a CLDN4 knockdown proteome. The only in-scope PRIDE project is the pan-claudin CoIP [PXD031094](https://www.ebi.ac.uk/pride/archive/projects/PXD031094).

This does not revise the locked CosMx, concordant-4, GSE137244, TCGA, or TISMO results.

---

## What was searched

| Catalog | Query | Result |
| --- | --- | --- |
| PRIDE keyword `CLDN4`, `claudin-4`, `Cldn4` | project search | 0 projects |
| PRIDE keyword `claudin` | project search | 30 projects, classified in `results/inventory_pride.tsv` |
| PRIDE in scope | text plus a `Cldn4` filename | **PXD031094 only** |
| MassIVE catalog | 20,228 datasets, title and description | **0** mention CLDN4 or claudin-4 |
| MassIVE claudin mention | same catalog | MSV000095998, CLDN10B in renal cancer |

PRIDE keyword `claudin` returns 11 projects whose text does not describe a claudin protein (keyword false positives), 12 other claudin mentions, 5 other-claudin knockdown or deficiency designs, and 1 “claudin-type” AMPA auxiliary-subunit project (PXD044621, GSG1L). Those are listed in the inventory and were not treated as CLDN4 datasets.

Nearest knockdown proteomes that are a different gene:

- PXD066158, CRISPR CLDN6-deficient germ-cell tumor cells.
- MSV000095998, CLDN10B in renal cancer.

## Open tables that were scanned

Targets: DNA-PKcs `PRKDC`, Ku80 `XRCC5`, Ku70 `XRCC6`, plus `XRCC4`, STING `STING1`/`TMEM173`, cGAS `CGAS`/`MB21D1`. Counts are from `results/scan_summary.json`.

| Table | Design | Rows scanned | Target hits |
| --- | --- | ---: | ---: |
| PXD031094 `proteinGroups_Cldn4.txt` | Full-length human CLDN4 GFP-Trap CoIP in MDCK-C7, n=4 vs 4 GFP. MaxQuant LFQ. SHA-1 matches the live PRIDE checksum. | 2,788 protein groups | **0** |
| Fredriksson 2015 PLOS S2, BL-Cldn4 column | BioID, biotin ligase on the N-terminus of claudin-4, MDCK II. Complete identification list. Not deposited in PRIDE or MassIVE. | 953 proteins | **0** |
| Fredriksson 2015 PLOS S3 | Three-fold enriched subset of the same BioID. Shared strings. | 1,820 strings | **0** |
| Villagómez 2024 Figshare Supplementary Table 1 | BioID “Claudin-4 proximal proteins” (OVCAR3 experiment reported in Heliyon 2022). One control column and one CLDN4 column. MD5 matches Figshare. Raw MS was not in PRIDE; the paper says raw data are available from the authors. | 411 genes | **0** |

Of the 411 BioID genes, 309 have the printed log2(CLDN4/control) greater than 1. That column is a single paired ratio, not a statistical partner call. The Heliyon 2022 text reports 52 proteins after a stricter filter (difference greater than 10 and log2 greater than 1 across four experiments). The stricter supplement PDF was not retrieved here (PMC `mmc1.pdf` returned 404). The Figshare table still contains the proteins those papers name from this BioID (RAB7A, SLC1A5, YES1, NDRG1, TRIP11, GOLGB1) and does not contain the DNA-PKcs / Ku / STING symbols.

The same MaxQuant file recovers the bait and a previously ranked partner, which shows that the scan is reading FASTA gene symbols:

| Gene | UniProt | Unique peptides | Cldn4 LFQ (4 replicates) | GFP LFQ |
| --- | --- | ---: | --- | --- |
| CLDN4 | O14493 | 3 | 1.99e9, 1.51e9, 3.17e9, 3.69e9 | 0, 0, 0, 0 |
| OAS1 | F1PLW6 | 5 | 1.04e8, 9.99e7, 1.04e8, 9.08e7 | 0, 0, 0, 0 |

PARP1 (J9NXE3; F1Q2M3) is in the MaxQuant file with 2 unique peptides. Its only nonzero LFQ is GFP replicate 3 (`38298000`). All four CLDN4 LFQ values are 0.

## Not scanned

`proteinGroups_PRISMA.txt` is the C-terminal peptide-matrix screen in the same PRIDE project (66,679,807 bytes, API SHA-1 `770554ea6737df2a850b08fa2404d0f5be0f75b6`). The PRIDE FTP host did not complete a TLS handshake from this environment, so that file was not read. RAW files (including the four `CoIP_Cldn4` RAW files) were not downloaded.

## CLDN4 knockdown protein abundance

No open PRIDE or MassIVE dataset is a CLDN4 knockdown or knockout proteome. Papers that change CLDN4 and discuss DNA repair or STING do not deposit that experiment in those archives:

- Yamamoto et al. 2022, PMID 35373300. Claudin-4 knockdown and PARP-inhibitor sensitivity. DNA-repair readouts include 53BP1 and XRCC1. The data-availability text points to public RNA and a prior ovarian cell-line mass spectrometry resource, and says RPPA data are available from the authors. No PXD accession.
- Villagómez et al. 2024, PMID 38867360. Claudin-4 knockdown with metabolomics (their Supplementary Table 2) and the BioID table scanned above. Data availability: data are in the paper and supplement; flow-cytometry files from the authors. No PXD accession.
- Villagómez et al. 2025, PMID 41214101. Claudin-4 knockdown or overexpression with immunoblot and microscopy of pSTING. Public data named in the paper are TCGA (dbGaP PHS000178). The BioID cited there is the earlier OVCAR3 experiment, not a new STING immunoprecipitation.

## 中文

公开的 CLDN4 CoIP（PXD031094，2788 个蛋白组）和两份 claudin-4 BioID 表（MDCK 953 个蛋白；卵巢细胞 BioID 411 个基因）里，PRKDC、XRCC5、XRCC6、STING1、TMEM173 的命中数都是 0。PRIDE 用 CLDN4 / claudin-4 检索到 0 个项目；MassIVE 目录 20228 个数据集里没有 CLDN4。没有找到 CLDN4 敲低后的公开蛋白组。PRISMA 大表这次没有下到，未扫描。没有改动已锁定的 CosMx、concordant-4、GSE137244、TCGA 和 TISMO 结论。
