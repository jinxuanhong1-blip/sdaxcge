# Lung-cancer CLDN4 siRNA / shRNA / CRISPR, 2018–2026

Public literature and deposit hunt only. No private matrices were used. No new accession was found, so no matrix was downloaded and no IFN score was computed.

Search date: 2026-09-21.

## Answer

**No public omics deposit of a lung-cancer CLDN4 knockdown or knockout was found.**

Two papers in 2018–2026 actually perturb CLDN4 in lung-cancer cells. Neither deposited a CLDN4-perturbation transcriptome in GEO, SRA, GSA, GSA-Human, NODE (via OmicsDI), or CNCB BioProject.

| PMID | Year | Model | Perturbation | Omics of the CLDN4 arm | Deposit |
|---|---|---|---|---|---|
| [41016339](https://pubmed.ncbi.nlm.nih.gov/41016339/) | 2025 | H1688 SCLC | CRISPR knockout | RNA-seq (SAA1 called as a downstream gene) | **None found** |
| [33006362](https://pubmed.ncbi.nlm.nih.gov/33006362/) | 2020 | A549 and H1299 NSCLC | lentiviral shRNA (rescue in CRAD-silenced cells) | none (qPCR / viability). Microarray is CRAD knockdown, not CLDN4 | **None.** “available from the corresponding author upon request.” |

A third A549 CLDN4 siRNA paper ([39502213](https://pubmed.ncbi.nlm.nih.gov/39502213/), 2024) studies LPS acute lung injury, not lung cancer, and also has no public matrix.

IFN scoring was not run. There is no new count matrix to score.

## What counted

A paper counted if all of the following were true:

1. Publication date 2018–2026.
2. The perturbed gene is CLDN4 / claudin-4 (siRNA, shRNA, or CRISPR / knockout), not another gene whose readout happens to include CLDN4.
3. The cells or tumours are lung cancer (NSCLC, SCLC, LUAD, LUSC, or a lung-cancer line used as such).

A deposit counted if a genome-wide omics matrix or read set was in GEO, SRA, GSA, GSA-Human, NODE, or CNCB, whether or not the paper cited the accession.

## Papers without a deposit

### PMID 41016339 — Kashiwagi et al., Biochem Biophys Res Commun 2025

CLDN4 knockout was made in H1688 small-cell lung cancer cells. Knockout increased proliferation. The abstract states that RNA-seq identified SAA1 as upregulated after CLDN4 knockout.

The full text is not in PMC (Europe PMC: `inPMC=N`, `hasSuppl=N`, `hasTMAccessionNumbers=N`, `hasDbCrossReferences=N`). Searches on 2026-09-21 found no matching series:

- GEO: `H1688 AND (CLDN4 OR SAA1)` = 0; `CLDN4 AND (SCLC OR "small cell") AND gse` = 0.
- SRA: `H1688 AND CLDN4` = 0; `Dokkyo AND (CLDN4 OR H1688 OR SAA1)` = 0.
- Author GEO under Yazawa is other experiments (GSE310370 plasma miRNA after nivolumab; GSE281782 STMN1 polyamide in SCLC), not this knockout.
- GSA keyword `CLDN4`: 11 records, all previously known non-matches (below).
- GSA-Human `CLDN4`: 0 records.
- OmicsDI `H1688 AND CLDN4`: 0 datasets.
- CNCB BioProject `CLDN4`: 16 projects, none an H1688 CLDN4 knockout.

### PMID 33006362 — Cui et al., Biosci Rep 2020

Lentiviral shRNA against claudin-4 was used in A549 cells as a rescue after CRAD silencing. The deposited-style assay in the paper is a microarray of CRAD knockdown, not of CLDN4 knockdown. Data availability: “The data used to support the findings of this study are available from the corresponding author upon request.” GEO `CRAD AND (A549 OR NSCLC OR lung)` returned 0 series. The full text contains no GSE, SRA, GSA, or CNCB accession.

## A549 siRNA that is not a lung-cancer study

PMID 39502213 (Zheng et al., Heliyon 2024) transfected siRNA against claudin-4 in A549 cells. The experiment is LPS acute lung injury. No transcriptome was generated. Data availability: “Data are available upon reasonable request.”

## Already scored elsewhere — not this wave

These are real CLDN4-loss transcriptomes and were not re-downloaded or re-scored:

| Accession | Why it is outside this wave |
|---|---|
| GSE50927 / PRJNA219385 | Mouse whole-lung germline Cldn4 knockout ± ventilator injury. Published 2014. Not a lung tumour. |
| GSE207704 / PRJNA856719 | MCF7 and T47D breast CLDN4 CRISPR RNA-seq. |
| GSE22493 / PRJNA128653 | SKOV-3 ovarian CLDN4 siRNA microarray. |

## Checked and rejected

H1688 RNA-seq that is not a CLDN4 perturbation:

- PRJNA1005054 / SRP455878 (SRR25694432–SRR25694437, libraries S1_1–S1_3 and S2_1–S2_3). Description: oridonin on SCLC. Hangzhou First People’s Hospital. Not CLDN4.
- GSE / SRA shKLF9 in NCI-H1688 (GSM9043561–GSM9043564). Not CLDN4.

Other lung papers that mention CLDN4 next to a knockdown of a different gene, with no CLDN4 siRNA/shRNA/CRISPR arm: CGN/FOXO1 in A549 (PMID 38338691), ELF3 in DMS53 (PMID 36840413, GSE190618), Cx26 in Calu-3 (PMID 40380298; “no datasets were generated”). Observational CLDN4+ malignant pleural effusion scRNA (PMID 38629624; CNCB PRJCA023797) does not perturb CLDN4.

Arabi et al. 2024 (PMID 39364319) wrote that, to their knowledge, no study had examined a mechanistic role of claudin-4 in LUAD. That review predates the 2025 H1688 knockout and does not cancel the 2020 A549 rescue shRNA above.

## IFN

No new accession, so no download and no IFN / MHC score. A score from the author-only RNA-seq in PMID 41016339 cannot be computed from public data.
