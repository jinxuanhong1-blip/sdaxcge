# Bitler lab claudin-4 papers: accessions, RPPA, supplements

Public files only. This folder inventories every Bitler-lab (University of Colorado Anschutz) paper that studies claudin-4 / CLDN4, records every omics accession and the RPPA experiment, and stores the supplementary files that could be downloaded.

No Bitler claudin-4 paper deposited a GEO, SRA, PRIDE, MetaboLights, or RPPA-portal accession of its own. The lab-generated omics that exist in public form are supplementary tables (BioID, metabolomics, TCGA reanalysis, drug screen). Those tables are copied here as the original files and as TSV.

## Papers

| Folder | Paper | PMID / PMCID | Own omics deposit |
|---|---|---|---|
| `mct2022_yamamoto/` | Yamamoto et al. Loss of Claudin-4 Reduces DNA Damage Repair and Increases Sensitivity to PARP Inhibitors. *Mol Cancer Ther* 2022. DOI 10.1158/1535-7163.MCT-21-0827 | 35373300 / PMC8988515 | RPPA run, matrix not released |
| `scirep2025/` | Villagomez et al. Claudin-4 as a dual regulator of genome stability and immune evasion in high grade serous ovarian cancer. *Sci Rep* 2025. DOI 10.1038/s41598-025-23137-1 | 41214101 / PMC12603150 | none (TCGA reanalysis + supplements) |
| `aacr2025_lb387/` | Villagomez et al. AACR 2025 abstract LB387. DOI 10.1158/1538-7445.AM2025-LB387 | abstract | none |
| `crc2024_autophagy/` | Villagomez et al. Claudin-4 Modulates Autophagy via SLC1A5/LAT1 as a Mechanism to Regulate Micronuclei. *Cancer Res Commun* 2024. DOI 10.1158/2767-9764.CRC-24-0240 | 38867360 / PMC11218812 | BioID + metabolomics in supplements |
| `crc2025_genome/` | Villagomez et al. Claudin-4 Stabilizes the Genome via Nuclear and Cell-Cycle Remodeling. *Cancer Res Commun* 2025. DOI 10.1158/2767-9764.CRC-24-0558 | 39625235 / PMC11705808 | none (figure supplements only) |
| `heliyon2022_neville/` | Neville, Webb, Baumgartner, Bitler. Claudin-4 localization in epithelial ovarian cancer. *Heliyon* 2022. DOI 10.1016/j.heliyon.2022.e10862 | 36237976 / PMC9552118 | BioID table in the PDF supplement |
| `mcr2019_breed/` | Breed et al. Ovarian Tumor Cell Expression of Claudin-4 Reduces Apoptotic Response to Paclitaxel. *Mol Cancer Res* 2019. DOI 10.1158/1541-7786.MCR-18-0451 | 30606772 / PMC6916652 | none; reanalysis of GSE18521 |

Preprints of the two Cancer Research Communications papers (10.1101/2024.01.18.576263 and 10.1101/2024.09.04.611120) state the same data policy as the journals. bioRxiv returned HTTP 429 during this download, so journal supplements were used.

A 2017 mammary-gland claudin paper (Baumgartner et al., PMID 28455726) is not an ovarian claudin-4 omics study and is not included. The 2022 SGO abstract (Gynecol Oncol, DOI 10.1016/s0090-8258(22)01453-6) is a meeting version of the MCT paper and has no separate data.

## RPPA (MCT 2022) — no accession, no matrix

Methods and Figure 3A: protein from OVCAR3 shCtrl (n=3) and shCLDN4#1 (n=3) was run on a reverse-phase protein array at the MD Anderson RPPA Core (NCI CA16672; Yiling Lu R50CA221675). The heatmap shows proteins with p<0.05 and FDR<15%. The text names XRCC1 and 53BP1 as decreased after claudin-4 knockdown and confirms them by immunoblot.

What is public:

- `mct2022_yamamoto/rppa/figure3_panelA_is_rppa_heatmap.jpg` — published Figure 3. Panel A is that heatmap. Panels B–G are immunoblots, DepMap correlations, the drug screen, and the Qian single-cell plot. This JPEG is not a numeric matrix.
- No supplementary xlsx is the RPPA table. Supplements 6–9 are the TCGA transcript table, cell-line BRCA status, the drug/DepMap table, and the GTFB tumor table.
- Data-availability statement: data are available upon request from the corresponding author.

## Lab-generated omics that are in the supplements

**Heliyon 2022 BioID.** OVCAR3 and OVCAR8 expressing claudin-4–biotin ligase. Averaged intensities are Supplementary Table 1 inside `mmc1`. Parsed to `heliyon2022_neville/tables/supp_table_1_bioid_averaged.tsv` (54 symbols). No PXD.

**CRC 2024 BioID.** Supplementary Table 1, control vs CLDN4 spectral counts and log2 fold change, 412 proteins: `crc2024_autophagy/tables/supp_table_1_bioid_cldn4_proximal_proteins.tsv`. No PXD. This is a separate table from the Heliyon averages.

**CRC 2024 metabolomics.** Supplementary Table 2, four blocks, WT vs CLDN4 KD (n=3):

- `supp_table_2_metabolomics_OVCA429_cells.tsv` — sample IDs DS2-041-001+ to 006+
- `supp_table_2_metabolomics_OVCA429_supernatant.tsv` — DS2-041-007+ to 012+
- `supp_table_2_metabolomics_OVCAR3_cells.tsv` — DS2-026-001+ to 006+
- `supp_table_2_metabolomics_OVCAR3_supernatant.tsv` — DS2-026-007+ to 012+

No MetaboLights or Metabolomics Workbench accession is stated. Flow-cytometry raw data are upon request.

## Accessions they used but did not generate

Full list with URLs is `ACCESSIONS.tsv`.

| Accession | Where it is used | In this folder |
|---|---|---|
| GSE18521 (GPL570) | MCR 2019, CLDN4 in 10 HOSE vs 53 HGSOC. Probe `218182_s_at`. | series matrix + CLDN4 row |
| phs000178 | TCGA OV, via cBioPortal (Firehose Legacy on 2020-12-20 in MCT; PanCancer Atlas on 2020-04-14 in Sci Rep and CRC 2025) | reanalysis tables only, not the controlled raw files |
| PXD003668 (MSV000080814) | Coscia 2016 ovarian-line proteome, MCT Figure S3 | not copied (RAW archive) |
| E-MTAB-8107 (+ E-MTAB-6149, E-MTAB-6653) | Qian 2020 blueprint scRNA, MCT Figure 3G, via blueprint.lambrechtslab.org on 2021-04-09 | not copied (large count matrices) |
| DepMap / GDSC2 | MCT drug-AUC correlations, accessed 2021-12-09 | the paper's exported correlation table |

## Supplements downloaded

Journal supplements are under each paper's `supplements/` directory. TSV exports of the numeric tables are under `tables/`.

MCT 2022 (PMC author manuscript, via the Internet Archive because live PMC bin URLs require reCAPTCHA): all 9 supplements, including the four xlsx tables and supplementary figures 1–4 (png).

Sci Rep 2025: MOESM1–MOESM5 (PDF, docx, xlsx). MOESM3 is the spectral-flow antibody panel.

CRC 2024: supplementary tables 1–4, supplementary figures 1–4, movies 1–4.

CRC 2025: supplementary figures 1–8 only.

Heliyon 2022: `mmc1` PDF.

MCR 2019: 3 of 8 supplements (Sup Fig 1, 5, 6). The other five were not archived and could not be fetched through the reCAPTCHA wall. URLs are in `mcr2019_breed/MISSING_SUPPLEMENTS.txt`.

AACR LB387: no file to download.

## Checksums

`checksums.sha256` covers every file in this directory except itself.
