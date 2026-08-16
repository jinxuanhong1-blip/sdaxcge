# Data Catalog — TROP2 / CLDN ADC ± ICI in Lung Cancer

Parallel slice `fable_adc`. This catalog lists **only public resources whose accessions/DOIs were
verified during this work** (GEO record pages, journal DOIs, dbGaP study pages, press releases).
No accession is invented. Where patient-level omics are **not** public, that is stated explicitly.

Scope: TROP2-directed antibody–drug conjugates (ADCs) — sacituzumab govitecan (SG),
datopotamab deruxtecan (Dato-DXd) — with or without immune-checkpoint inhibitors (ICI) in
lung cancer; plus CLDN-targeted agents (CLDN18.2) where any public omics touch lung.

Last verified: 2026-08-16.

---

## A. Analyzable public omics used in this slice (open, patient-level)

These are fully open GEO records with processed expression **and** outcome/response annotation.
None is a TROP2-ADC-treated cohort (see §C) — they are **ICI ± chemo NSCLC cohorts** used as
public surrogates to characterize the ADC targets (TROP2 = *TACSTD2*; CLDN18.2 gene = *CLDN18*),
the SG payload target (*TOP1*), and PD-L1 (*CD274*) versus immunotherapy outcome.

| Accession | Design / treatment | n (used) | Outcome field | Files used | Platform |
|---|---|---|---|---|---|
| [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044) | NSCLC, anti–PD-1 monotherapy | 16 | responder / non-responder | `GSE126044_counts.txt.gz` (gene-symbol raw counts) + series matrix | Illumina HiSeq 2500 (GPL16791) |
| [GSE135222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222) | Advanced NSCLC, anti–PD-1/PD-L1 | 27 | PFS event + PFS time (days) | `GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz` (Ensembl TPM) + series matrix | Illumina HiSeq 2500 (GPL16791) |
| [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) | Neoadjuvant chemo + anti–PD-1 (toripalimab / camrelizumab / sintilimab), **pre-treatment** biopsy | 29 bulk | pathologic response (MPR incl. pCR / NMPR), residual-tumor fraction, RECIST | `GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz` + `..._metadata.xlsx` | Illumina NovaSeq 6000 (GPL24676) |

Target gene identifiers used:

| Gene (symbol) | Meaning for this slice | Ensembl (for GSE135222) |
|---|---|---|
| `TACSTD2` | TROP2 — SG / Dato-DXd target | ENSG00000184292 |
| `CLDN18` | Claudin-18 (isoform 18.2 is the therapeutic epitope; bulk RNA-seq cannot resolve isoform) | ENSG00000066405 |
| `TOP1` | Topoisomerase-1 — target of SG payload SN-38 / Dato-DXd payload DXd | ENSG00000198900 |
| `CD274` | PD-L1 — ICI context / SG+pembro combination rationale | ENSG00000120217 |

> Note on GSE207422: this is the closest public analog to "TROP2 target ± ICI in lung," because it is a
> chemo-**immunotherapy** neoadjuvant NSCLC cohort with pre-treatment transcriptomes and graded
> pathologic response. It also carries a matched scRNA-seq matrix (not used here).

## B. Additional public omics resources (cataloged, not analyzed here)

| Resource | Relevance | Access |
|---|---|---|
| [TCGA LUAD / LUSC](https://xenabrowser.net/datapages/) via UCSC Xena (RRID:SCR_018938) | Large-cohort prognostic context for *TACSTD2* / *CLDN18* mRNA vs OS in NSCLC | Open (bulk download) |
| [GSE302284](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302284) | TROP2 in EGFR-mutant NSCLC drug-tolerant persisters; SG vs TROP2 CAR-T in models | Open (GEO) |
| [PMC11999141](https://pmc.ncbi.nlm.nih.gov/articles/PMC11999141/) | TROP2 mRNA/protein prevalence in NSCLC using TCGA LUAD/LUSC (reference for target biology) | Open (article) |

## C. TROP2-ADC clinical trials — biomarker findings WITHOUT public patient-level omics

Cataloged from abstracts / peer-reviewed papers / company press releases. **No public omics
accession exists** for these; do not fabricate one. Patient-level data are controlled or unreleased.

| Trial (ID) | Agent(s) | Setting | Public biomarker finding | Source |
|---|---|---|---|---|
| TROPION-Lung01 ([NCT04656652](https://clinicaltrials.gov/study/NCT04656652)) | Dato-DXd vs docetaxel | 2L+ advanced/metastatic NSCLC | TROP2 QCS **normalized membrane ratio (NMR)**: QCS-NMR+ (≥75% of tumor cells with ratio ≤0.56) enriched PFS/ORR benefit for Dato-DXd | [JCO 10.1200/JCO-24-01544](https://ascopubs.org/doi/10.1200/JCO-24-01544); [IASLC WCLC24 PL02.11](https://www.iaslc.org/iaslc-news/press-release/normalized-membrane-ratio-trop2-quantitative-continuous-scoring-predictive); [Daiichi Sankyo PR](https://www.daiichisankyo.com/files/news/pressrelease/pdf/202409/20240908_E.pdf) |
| EVOKE-01 ([NCT05089734](https://clinicaltrials.gov/study/NCT05089734)) | SG vs docetaxel | Previously treated metastatic NSCLC | TROP2 IHC H-score highly expressed but **not** predictive of OS benefit; possibly negative prognostic | [AACR 2025 Abstract LB260](https://doi.org/10.1158/1538-7445.am2025-lb260) |
| EVOKE-02 ([NCT05186974](https://clinicaltrials.gov/study/NCT05186974)) | SG + pembrolizumab ± carboplatin | 1L metastatic NSCLC | Efficacy (ORR/PFS) **independent** of TROP2 H-score (median cut) | [AACR 2025 Abstract LB399](https://doi.org/10.1158/1538-7445.am2025-lb399) |
| EVOKE-03 | SG + pembrolizumab | 1L NSCLC, PD-L1 TPS ≥50% | Phase III ongoing (context) | referenced in LB399 |

## D. CLDN18.2-targeted agents — lung relevance & omics status

| Item | Relevance to lung | Public omics? | Source |
|---|---|---|---|
| CLDN18.2 IHC in lung invasive mucinous adenocarcinoma (IMA) | IMA can express CLDN18.2 → potential zolbetuximab / CLDN18.2-ADC / CAR-T candidates | IHC only (no open omics accession) | [Cancer Diagn Progn article 583](https://www.cancerdiagnosisprognosis.org/article/583/expression-of-cldn182-in-invasive-mucinous-adenocarcinomas-of-the-lung) |
| CLDN18.2 pan-cancer analysis | Uses TCGA (UCSC Xena) + TCGASpliceSeq; *CLDN18* expression/PSI, prognosis | Open via TCGA/Xena | [Front. Pharmacol. 10.3389/fphar.2024.1494131](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2024.1494131/full) |
| Agents: zolbetuximab (mAb), CMG-901 / SOT102 / CPO102 (ADC), CT041 (CAR-T) | CLDN18.2 therapeutic landscape; lung is a minor/mucinous-restricted indication | No lung-specific trial omics public | landscape reviews |

In this slice, CLDN18.2 is analyzed at the **gene level (`CLDN18`)** in the three NSCLC cohorts above
to describe its (expectedly low/rare) expression in unselected NSCLC and any outcome association.
Bulk RNA-seq cannot distinguish the 18.1 vs 18.2 isoform — this is stated as a caveat.

---

## Provenance / integrity notes
- All GEO file names and sizes were confirmed by direct download from the NCBI FTP mirror.
- dbGaP `phs002822.v1.p1` (SU2C-MARK ICI NSCLC) exists and is **controlled-access** — not used for
  analysis here; listed only for completeness.
- No accession, DOI, or NCT number in this catalog was invented; each is linked to its source.
