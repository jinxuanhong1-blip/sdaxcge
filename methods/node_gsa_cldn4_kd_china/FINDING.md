# CLDN4 knockdown in NODE / GSA: lung, ovarian, breast

**0 Chinese accessions are a CLDN4 knockdown in lung, ovarian, or breast cancer. 0 files were downloaded.**

NODE, GSA-Human (HRA), and OMIX each return 0 records for `CLDN4`. GSA’s 11 `CLDN4` hits are INSDC mirrors (SRX), not CRA projects. The only Chinese BioProject whose title contains CLDN4 is not a knockdown and is controlled.

## Counts

| Archive | Query | Records | CLDN4 knockdown, lung/ovary/breast |
|---|---|---:|---|
| NODE | `CLDN4`, `claudin-4` | 0 | 0 |
| NODE | `knockdown` / `shRNA` / `siRNA` | 63 / 43 / 16 | 0 mention CLDN4 |
| NODE | `CLDN` | 31 | all CLDN6 colorectal (OEP00006646), not CLDN4 |
| GSA-Human | `CLDN4` | 0 | 0 |
| GSA-Human | `knockdown` | 105 | 0 mention CLDN4 |
| OMIX | `CLDN4` | 0 | 0 |
| OMIX | `knockdown` | 87 | 0 mention CLDN4 |
| OMIX | `claudin` / `CLDN` | 7 / 5 | CLDN18.2 gastric or trial biomarker tables, not CLDN4 |
| GSA | `CLDN4` | 11 | INSDC experiments, 0 CRA / 0 PRJCA |
| BioProject | `CLDN4` | 16 | 1 PRJCA, and that project is not a knockdown |

## Nearest Chinese lung deposit (not downloaded)

HRA006761 / PRJCA023797. Shanghai Pulmonary Hospital (Caicun Zhou). PMID 38629624. Sixteen baseline malignant pleural effusion scRNA-seq samples from advanced NSCLC (8 recurrent, 8 non-recurrent). CLDN4 is a marker in the paper, not a silenced gene. GSA-Human lists the study as **controlled**. The public catalog has 16 runs (HRR1529727–HRR1529742), 32 FASTQ files, 100.53 GB, Illumina NextSeq 550. NODE has no copy. Metadata xlsx is restricted to authorized users. Sequence was not requested.

## Chinese papers that perturb CLDN4, with no omics deposit

- PMID 26058359 (2015, Guangdong): CLDN4 shRNA and overexpression in MCF-7. Proliferation, apoptosis, migration, xenograft. No GEO, HRA, OMIX, or OEP.
- PMID 25871476 (2015, Beijing): CLDN4 antibody blockade of vasculogenic mimicry in MDA-MB-231. Not a deposited knockdown transcriptome.

Papers that knock down a different gene and read CLDN4 out also have no NODE/GSA accession: PMID 30808546 (PAK4 knockdown, breast), PMID 33006362 (CRAD knockdown microarray in A549/H1299; GEO has no matching series; NODE query `CRAD` is 0), PMID 32779991 (ELFN1-AS1 knockdown, ovarian).

## Open CLDN4-knockdown transcriptomes that are not Chinese deposits

These appear in the CNCB BioProject/GSA index because CNCB mirrors INSDC. They were left at GEO/SRA.

- GSE22493 / PRJNA128653. SKOV-3 CLDN4 siRNA microarray, n=3. Harvard Medical School, USA.
- GSE207704 / PRJNA856719. T47D and MCF7 CLDN4-knockout RNA-seq. Fukushima Medical University, Japan. PMID 37059993. GSA lists SRX16067609 and SRX16067610.
- GSE50927 / PRJNA219385. Mouse lung, WT vs Cldn4 knockout after ventilator-induced lung injury. University of Southern California. Not cancer. PMID 25106430.

GSA also lists 4C-seq at the CLDN4 locus in gastric lines YCC21 and SNU16 (SRX8386834–SRX8386842). That is capture of the locus, not knockdown of CLDN4, and not lung, ovary, or breast.
