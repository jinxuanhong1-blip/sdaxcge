# Source and scope notes

Search frozen on 2026-08-16 (UTC). “Direct” means that repository metadata or the
primary paper explicitly says that the samples/data came from the named trial.
Records that only cite a trial, or independent cohorts using a similar regimen,
were excluded.

## Primary identifiers

| Trial | Registry |
|---|---|
| IMpower150 | [NCT02366143](https://clinicaltrials.gov/study/NCT02366143) |
| IMpower110 | [NCT02409342](https://clinicaltrials.gov/study/NCT02409342) |
| IMpower130 | [NCT02367781](https://clinicaltrials.gov/study/NCT02367781) |
| PACIFIC | [NCT02125461](https://clinicaltrials.gov/study/NCT02125461) |

## Direct IMpower150 deposits

1. [EGA EGAS00001006703](https://www.ega-archive.org/studies/EGAS00001006703)
   links four datasets: ctDNA features/mutation calls
   (EGAD00001009725), clinical data (EGAD00001009726), analysis/model artifacts
   (EGAD00001009764), and a 311-gene list (EGAD50000000273). The
   [Nature Medicine paper](https://doi.org/10.1038/s41591-023-02226-6)
   explicitly states that IMpower150 clinical and ctDNA data and documented R
   code were deposited under this study accession.
2. [EGA EGAS50000001272](https://www.ega-archive.org/studies/EGAS50000001272)
   and dataset
   [EGAD50000001814](https://www.ega-archive.org/datasets/EGAD50000001814)
   are explicitly described as baseline IMpower150 tumor transcriptomes and
   relevant clinical metadata. The three processed CSVs total 325,295,332
   bytes.
3. The official
   [EGA public metadata API](https://ega-archive.org/discovery/metadata/public-metadata-api/)
   reports `access_type: controlled` for every dataset above. File size below
   2 GB does not make a controlled EGA object openly downloadable.

## Cross-trial open processed result

[Khan et al., Genome Medicine 2023](https://doi.org/10.1186/s13073-023-01193-4)
analyzed WGS across 14 Roche trials, explicitly including IMpower110,
IMpower130, and IMpower150. Raw WGS cannot be placed in a public/controlled
archive under the consents; the paper instead makes pooled common-variant GWAS
summary statistics public:

- [LocusZoom 74850](https://my.locuszoom.org/gwas/74850), whole pooled cohort,
  direct file 304,844,876 bytes.
- [LocusZoom 743668](https://my.locuszoom.org/gwas/743668), taxane-treated
  subcohort, direct file 279,782,143 bytes.

These are processed, public, and below 2 GiB, so the download script permits
them. They are **not trial-specific datasets** and cannot support an
IMpower110-versus-130-versus-150 comparison.

The pooled hypothyroidism genetics study also includes IMpower130 and exposes
three public PGS Catalog scoring files
([PGP000164](https://www.pgscatalog.org/publication/PGP000164/);
PGS000759, PGS000760, and PGS000761). These are reusable score weights, not
IMpower130 genotypes or trial-stratified results.

## Direct open IMpower150 derivatives

Three BioStudies/Europe PMC publication packages expose small processed
workbooks derived directly from IMpower150:

- [S-EPMC10115641](https://www.ebi.ac.uk/biostudies/studies/S-EPMC10115641):
  ctDNA cohort/sample and analysis summary tables, 182,673 bytes.
- [S-EPMC11316765](https://www.ebi.ac.uk/biostudies/studies/S-EPMC11316765):
  figure-level ctDNA source data, 250,646 bytes.
- [S-EPMC12775477](https://www.ebi.ac.uk/biostudies/studies/S-EPMC12775477):
  figure-level transcriptomic subtype/biomarker source data, 59,235 bytes.

These workbooks are truly open and processed, so they are downloaded. They
contain summaries/source values, not the controlled patient-level matrices.
Narrative PDFs and peer-review/reporting files were excluded.

## Negative findings and request-only routes

- A complete exact-name scan of the official EGA public API covered 21,321
  datasets and 10,502 studies. It found direct records only for IMpower150.
  Generic `Pacific` hits (Pacific populations and Pacific Biosciences) were
  excluded. The reproducible scan is `scripts/gpt_impower/audit_ega.py`.
- The 2024
  [MOSBY paper](https://doi.org/10.1038/s41598-024-69198-6) says IMpower110
  and, at that time, IMpower150 data were not public pending primary biomarker
  manuscripts. IMpower150 has since acquired a 2025 EGA transcriptomic deposit;
  no analogous IMpower110 record was found in the 2026 EGA scan.
- The
  [IMpower130 primary report](https://doi.org/10.1016/S1470-2045(19)30167-6)
  directs qualified researchers to sponsor-mediated participant-level data
  requests rather than an omics archive.
- Original PACIFIC participant-level data are requestable under the
  [AstraZeneca/Vivli route](https://vivli.org/ourmember/astrazeneca/), but no
  direct public omics accession for the original PACIFIC trial was found.
- PACIFIC BioStudies records such as S-EPMC8412232 and S-EPMC9015199 are
  publication packages, not omics experiments. Their aggregate IHC/survival
  supplements were not classified as processed omics.
- [PRJNA1026052](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1026052) is an
  open TCR-sequencing study of a small, independent durvalumab-after-CRT
  cohort. It cites PACIFIC as the treatment precedent but is not the PACIFIC
  trial, so it is deliberately excluded from the accession table and download.

## Important disambiguations

- EGAS00001004888 and its IMpower133 datasets are **IMpower133 (SCLC)**, not
  IMpower130. Search snippets that label EGAS00001004888 as IMpower130 are
  wrong.
- IMvigor130 is a urothelial-cancer trial and is not IMpower130.
- PACIFIC-R, PACIFIC-2, PACIFIC-5, PACIFIC-6, and regimen-matched institutional
  cohorts are not the original PACIFIC randomized trial.
