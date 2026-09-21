# Sci Rep 2025 CLDN4 CRISPRi / pSTING: raw RNA-seq and proteomics hunt

**0 open raw RNA-seq files. 0 open raw proteomics files.** Nothing from this paper was deposited in GEO, SRA, or PRIDE. The open supplements were downloaded and read. They are uncropped immunoblots, a flow-antibody list, and two small TCGA cBioPortal tables. They are not CRISPRi count matrices and they are not mass-spec spectra.

No differential-expression or protein-abundance analysis was run. There is no CRISPRi transcriptome or proteome to score.

Paper: Villagomez et al. (Bitler senior), *Scientific Reports* 2025;15:39257. DOI [10.1038/s41598-025-23137-1](https://doi.org/10.1038/s41598-025-23137-1). PMID 41214101. PMC12603150. Published 10 November 2025.

CRISPRi in the paper: OVCAR3 and OVCA429, dCas9-KRAB-MeCP2 (Addgene 110824) plus gRNA `GCTGGCTTGCGCATCAGGAC`. Claudin-4 overexpression is in OVCAR8. pSTING is confocal microscopy and immunoblot (Figure 2), not sequencing.

---

## Data availability, as printed

> All materials are available upon request of the corresponding author. All data is presented in the manuscript, and all raw data is available upon request of the corresponding author.

Methods list immunoblot, immunofluorescence, autophagy-flux flow, ISRE luciferase, HIS-BRGS PDX flow cytometry, and public-database mining (cBioPortal TCGA ovarian serous cystadenocarcinoma, PanCancer Atlas, accessed 14 April 2020; dbGaP PHS000178 is the TCGA study, not a new deposit; TIMER; STRING). There is no RNA-seq library kit, no FASTQ, and no mass spectrometer.

The BioID result cited for CLDN4–Rab7 is a previous paper (Neville et al., *Heliyon* 2022), not a proteome generated here.

---

## Repositories checked (21 Sep 2026)

| Source | Query | Result |
| --- | --- | --- |
| GEO | `Villagomez[Author]` | 11 series. None are this paper (EHMT/PRMT5 ovarian or MCL; GBM chromatin). |
| GEO | `CLDN4 AND (CRISPRi OR CRISPR) AND (ovarian OR OVCAR)` | 0 |
| GEO | `Bitler[Author]` plus CLDN4 / claudin / STING / CRISPRi / OVCAR / CMP | 4 hits, all unrelated after MeSH expansion: GSE186618/619 colorectal ERV RNA-seq; GSE182089/197 SETDB1 CRISPR screen. |
| GEO elink | PubMed 41214101 | no linked GEO |
| SRA | PMID, and Bitler/Villagomez organization plus claudin / CLDN4 / STING / CRISPRi | 0 |
| BioProject | `Villagomez FR[Author] AND claudin` | 0 |
| PRIDE v2 | DOI, Bitler, Villagomez, `claudin-4` | 0 projects. Keyword `ovarian` returns PXD083384, so the endpoint answered. |
| MassIVE ProXI | `filter=claudin-4` | no datasets |
| BioStudies | S-EPMC12603150 | literature record only. 0 GSE / SRP / PRJNA / PXD |
| OmicsDI | this DOI | one literature hit, omics type Unknown |
| Figshare | DOI and `Villagomez claudin-4` | 0 |
| Zenodo | DOI | API returned HTTP 403 from this environment. A web search did not find a Zenodo record for this DOI. |

Full query strings are in `search_log.tsv`.

Same CRISPRi gRNA is used in the two sister *Cancer Research Communications* papers (CRC-24-0558, CRC-24-0240). BioStudies S-EPMC11705808 and S-EPMC11218812 also contain no GSE or PXD. The autophagy preprint says generated data are in the manuscript and supplement, and that datamining and flow raw files are available on request.

Not this experiment: GSE22493 (2010-12-15, SKOV-3 CLDN4 siRNA microarray).

---

## Supplements downloaded

Springer MOESM files for this DOI. SHA-256 and byte sizes are in `supplement_inventory.tsv`. Blot images were not committed.

| File | What it is |
| --- | --- |
| MOESM1 PDF, 2.3 MB, 9 pages | Uncropped whole membranes for Figure 2, Figure 6, and Supplementary Figures 4–10. Extracted text is the figure titles only. |
| MOESM2 DOCX, 15 KB | Supplementary Table 2. TCGA ovarian cBioPortal, 617 serous cystadenocarcinomas. Five CLDN4–ISG rows. Reprinted as `moesm2_tcga_isg_cbioportal.tsv`. Spelling “Co-ocurrence” is as printed. |
| MOESM3 XLSX, 12 KB, 36 rows | Supplemental Table 1. Flow antibodies and reagents (CD4, CD8, CD247/TCRz, IFNg, TNFa, and the rest of the HIS panel). Not a proteome. |
| MOESM4 DOCX, 8.9 MB | Supplementary figure legends and 21 PNG panels. No RNA-seq, proteomics, GEO, or PRIDE text. |
| MOESM5 DOCX, 19 KB | Supplementary Table 3. TCGA PanCancer Atlas ovarian protein values stratified by CLDN4 mRNA z-score. 18 proteins. Reprinted as `moesm5_tcga_protein_vs_cldn4.tsv`. |

Table 3’s Word header prints `CLDN4: EXP>1` and `CLDN4: EXP<-1` twice. The TSV keeps both columns and labels them by column index. The caption does not name those four columns, so they were not relabeled as mean or standard deviation. STING / pSTING is not in the table. This is a public TCGA protein summary already in the supplement, not raw spectra from the CRISPRi lines.

---

## What was not done

- No FASTQ, count matrix, or PRIDE raw/mzML download. None are open.
- No new TCGA or CPTAC reanalysis. The locked TCGA keratin result is unchanged.
- No private 8-KL, HIS, or PDX proteome was used. The HIS-PDX experiment in this paper is spectral flow, and its raw FCS is author-on-request.
