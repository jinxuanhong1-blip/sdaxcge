# AACR abstracts and bioRxiv 2023–2026: CLDN4 CRISPR/KD omics accessions

Screen date: 2026-09-21. Public sources only. No private matrices. No new differential statistics.

## Harvest

No GEO, SRA, BioProject, PRIDE, MetaboLights, Metabolomics Workbench, or ArrayExpress accession from 2023–2026 is a CLDN4 CRISPR, CRISPRi, shRNA, or siRNA contrast tied to an AACR abstract or a bioRxiv preprint.

One processed metabolomics matrix of that contrast is public, as a journal supplementary table on AACR figshare. It belongs to the same study as AACR abstract B072 and bioRxiv `10.1101/2024.01.18.576263`.

| Accession | What it is | CLDN4 KD contrast |
|---|---|---|
| `10.1158/2767-9764.28688934` | Supplementary Table 2, Villagomez 2024. Targeted metabolomics intensities. Figshare article 28688934. File `https://ndownloader.figshare.com/files/53292045` | **Yes.** OVCA429 and OVCAR3, cells and supernatant, WT n=3 vs CLDN4 KD n=3 |
| `10.1158/2767-9764.c.7312139` | Parent figshare collection for the Cancer Research Communications paper | Collection, not a separate matrix |
| `10.1158/2767-9764.28688937` | Supplementary Table 1. BioID proximal proteins (411 gene symbols; control sample vs CLDN4-bait sample, one each) | **No.** CLDN4 bait versus control, not the knockdown contrast |

The preprint data-availability line says the data are in the manuscript and its supplement. The stable public IDs above are the published figshare records (collection posted 2024-07-02; current table files dated 2025-03-28 on the figshare API). The Cancer Research Communications data-availability statement does not name a GEO or PRIDE accession.

## Design inventory of the KD matrix

Parsed from the docx. Intensities were not reanalyzed.

| Block | WT sample IDs | CLDN4 KD sample IDs | Compounds |
|---|---|---|---|
| OVCA429 cells | DS2-041-001+, 002+, 003+ | DS2-041-004+, 005+, 006+ | 155 |
| OVCA429 supernatant | DS2-041-007+, 008+, 009+ | DS2-041-010+, 011+, 012+ | 149 |
| OVCAR3 cells | DS2-026-001+, 002+, 003+ | DS2-026-004+, 005+, 006+ | 157 |
| OVCAR3 supernatant | DS2-026-007+, 008+, 009+ | DS2-026-010+, 011+, 012+ (labeled KD) | 147 |

CRISPRi in the paper and preprint: OVCA429 and OVCAR3, dCas9-KRAB-MeCP2, gRNA `GCTGGCTTGCGCATCAGGAC`. OVCAR8 in the same papers is overexpression, not knockdown.

## Same-study records with no repository accession

| Record | Perturbation | Deposit |
|---|---|---|
| AACR abstract B072, DOI `10.1158/1538-7445.ovarian23-b072` | CRISPRi of claudin-4 in OVCA429, OVCAR3, and OVCAR8, plus CMP | Abstract text describes imaging, flow, and immunoblot. The metabolomics table is on the later paper, not cited as a repository ID in the abstract |
| bioRxiv `10.1101/2024.09.04.611120` (2024-09-04) | Claudin-4 downregulation by CRISPRi in the full text; published as CRC `10.1158/2767-9764.CRC-24-0558` | Preprint: data are in the manuscript and supplement. Journal: TCGA dbGaP `PHS000178` reused; other data upon request. Figshare items for this paper are figures, not an omics table |
| AACR 2025 abstract LB387, DOI `10.1158/1538-7445.am2025-lb387` | Retrieved abstract text: CMP, Rab7, STING. CRISPR is not in that text. Companion Sci Rep `10.1038/s41598-025-23137-1` uses the same CRISPRi system | Sci Rep: raw data upon request. No new accession |

## Outside this wave

`GSE207704` is CLDN4 CRISPR RNA-seq in T47D and MCF7 (Murakami, Breast Cancer Research 2023; GEO public 2023-04-19). It is not an AACR abstract and no bioRxiv record was found for it. It was already scored in PR #193. Not reanalyzed here.

`CNP0006650` (Cell Reports Medicine 2025) is patient HCC scRNA-seq. That paper also uses CLDN4 shRNA, but the deposited single-cell data are patient tumors, not the knockdown. Not bioRxiv.

`PRJNA1424384` and Zenodo `10.5281/zenodo.18557106` (Nature 2026) are a genome-wide CRISPR screen that identified claudin-4 as a toxin receptor. The deposit is the screen, not a CLDN4-knockout transcriptome. Not an AACR abstract or bioRxiv.

## Files

- `accessions.tsv` — harvested IDs
- `sample_manifest.tsv` — KD metabolomics sample IDs
- `screened.tsv` — records opened in this wave
- `search_log.md` — queries
