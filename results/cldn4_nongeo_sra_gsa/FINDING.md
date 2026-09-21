# CLDN4 shRNA / CRISPR RNA-seq outside GEO

Queried 2026-09-21. Question: is there a public CLDN4 or Cldn4 shRNA, siRNA, or CRISPR RNA-seq in SRA/BioProject or CNCB (GSA, GSA-Human, OMIX) that is **not** already a GEO series?

**No open non-GEO matrix was found, so nothing was downloaded.**

The only public CLDN4-loss transcriptomes remain the ones already in GEO:

| Accession | What it is | Openness |
|---|---|---|
| GSE207704 / PRJNA856719 | Human MCF7 and T47D CRISPR CLDN4 knockout RNA-seq | SRA public (~1.8 GB for the four CLDN4-/- runs) |
| GSE50927 / PRJNA219385 | Mouse lung germline Cldn4 knockout RNA-seq | SRA public (~34 GB, 5 runs) |
| GSE22493 / PRJNA128653 | SKOV-3 CLDN4 siRNA microarray | GEO matrix; 0 SRA runs |

Those three are excluded here because they are in GEO. Their FASTQ was not pulled again.

## What was searched

NCBI E-utilities on 2026-09-21:

- BioProject text `CLDN4` / `Cldn4` / `claudin-4` / `claudin 4`: 31 projects.
- SRA text for the same tokens: 45 experiments, all of them already tied to a GEO series or to the 4C-seq near-miss below.
- Phrase queries `"CLDN4 shRNA"`, `"CLDN4 knockout"`, `"CLDN4 knockdown"`, `shCLDN4`, `sgCLDN4` in SRA: **0 runs**.
- BioProject `"CLDN4 knockout"`: only PRJNA856719. `"CLDN4 KO"` / `"Cldn4 KO"`: only PRJNA219385.
- PubMed union of CLDN4/Cldn4 knockout, knockdown, CRISPR, shRNA, and siRNA phrases plus the RNA-seq abstract query: 87 papers. Europe PMC full text was opened when a PMCID existed.

CNCB on the same day:

- GSA search `CLDN4` and `Cldn4`: 11 items, all INSDC mirrors of GSE50927, GSE207704, or the 4C-seq study. **0 CRA accessions.**
- GSA-Human `finished.json`: 7,312 studies. Strict CLDN4 / claudin-4 / Cldn4 text match: **0**.
- OMIX `releaseList`: 10,741 accessions in the complete HTML. **0** CLDN4 / claudin-4 strings.

BioStudies `CLDN4 AND (knockout OR knockdown OR shRNA OR CRISPR OR siRNA)` returned literature records and the ArrayExpress mirrors E-GEOD-22493 and E-GEOD-50927, not a new unique study.

## Papers that look like the missing dataset

**PMID 41016339** (Kashiwagi et al., BBRC 2025) does report RNA-seq of CRISPR CLDN4-knockout H1688 SCLC cells, with SAA1 up after knockout. No GSE, PRJNA, PRJDB, CRA, HRA, or OMIX accession is attached in Europe PMC, and a BioProject query for Kashiwagi plus CLDN4/H1688/SAA1/SCLC returned 0. The reads are not publicly downloadable.

**PMID 40892111** (Zheng et al., Functional & Integrative Genomics 2025) is not a CLDN4-knockdown RNA-seq. The RNA-seq is cerulein acute pancreatitis versus control; CLDN4 was picked from that list and then knocked down. Knockdown readouts are qPCR, Western blot, ELISA, and histology. The data-availability statement says the corresponding author will share data on reasonable request. The only open supplement (`10142_2025_1683_MOESM1_ESM.docx`, 14.3 KB) is a primer table (GAPDH, CLDN4, GPX4, ACSL4). It was inspected and not committed.

**PMID 41697223** (Cldn4 deletion, abdominal sepsis, 2026) reports permeability and survival, not RNA-seq, and has no BioProject.

## Open records that are not this experiment

These are public and not in GEO, but they are not CLDN4 shRNA/CRISPR RNA-seq, so the FASTQ/sgRNA files were not downloaded:

- **PRJNA1424384** — six public sgRNA-seq runs (~229 MB) from a genome-wide Avana screen of HT29 cells under Bacteroides fragilis toxin. The toxin binds claudin-4. The perturbed variable is the sgRNA library, not a CLDN4 RNA-seq contrast. PMID 42020735.
- **PRJNA634686** — six public 4C-seq runs (~4.8 GB) with the viewpoint at the CLDN4 locus in MSH2-knockout gastric lines. Chromatin, not expression.
- **HRA013141** — open GSA-Human ATAC-seq of NCI-H1688 after XPO1 knockdown.

## Bottom line

There is still no second public CLDN4-loss RNA-seq that lives only in SRA or in a CNCB archive. The H1688 CRISPR RNA-seq is real in the paper and absent from the public archives checked here.
