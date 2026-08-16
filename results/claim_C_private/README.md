# Claim C private — public analog hunt

**Scope.** C1–C10 are private. This folder records **public analogs only**. Private claim text is not copied or reconstructed.

**Date.** 2026-08-16. Nothing invented: every “FOUND” row has a DOI/PMID that was opened or API-confirmed.

## Found vs none

| # | Analog class (as requested) | Verdict | Public analog (if any) |
| --- | --- | --- | --- |
| 1 | TROP2–CLDN4 colocalization papers / supplements | **FOUND** | Nakatsukasa 2010 PLA in cornea; Van Itallie 2015 CLDN4-BioID lists TROP-2 |
| 2 | TROP2-ADC (SG / Dato / SKB264) proteomics showing tight-junction down | **NONE** | No SG / Dato-DXd / SKB264 treatment proteome reports claudin / TJ decrease |
| 3 | Co-IP **or** AlphaFold2 of TROP2–CLDN4, published | **NONE** | Published Co-IP is **negative**; no AF2 TROP2–CLDN4 complex paper |
| 4 | Cell-line IC50 vs TROP2 (SG / Dato / SKB264) | **FOUND** | Multiple preclinical papers with IC50/GI50 stratified by TROP2 |

Do not treat near-misses as hits. They fail the analog class as written.

## 1. Colocalization — FOUND

### Nakatsukasa et al., *Am J Pathol* 2010

- DOI `10.2353/ajpath.2010.100149`, PMID 20651236.
- **PLA (Fig. 2B–C):** positive in situ proximity-ligation signals between TACSTD2 and **CLDN4** (also CLDN1, CLDN7, TJP1) on **normal human cornea**; signals sit on the plasma membrane (desmoplakin overlay). OCLN PLA was negative.
- **Not Co-IP:** the same paper states CLDN4 / TJP1 / OCLN **did not** co-immunoprecipitate with TACSTD2 at any detergent concentration tested.
- Authors’ own reading: CLDN4 “closely resides” next to TACSTD2 (possibly via CLDN1/7) but **does not bind**.
- Journal supplements: Fig. S1–S3 at the AJP site (isotype controls; GDLD stains; desmosome PLA negative). Not a numeric coloc table.

### Van Itallie et al., *PLOS ONE* 2015

- DOI `10.1371/journal.pone.0117074`, PMID 25789658, PMC4366163.
- BioID of occludin and **claudin-4** in MDCK II. Text + **Table 2** name tumor-associated calcium signal transducer 2 (**TROP-2**) as a signaling protein **enriched only around Cldn4**.
- Assay is proximity biotinylation + MS, not dual immunofluorescence. Supporting information on PLOS is the MS inventory.

**Not found under this lane:** a published dual-IF / multiplex IF paper or supplement that reports a TROP2∩CLDN4 overlap fraction or paired H-scores in lung/NSCLC (see sibling hunt `results/hunt_ihc`). Europe PMC title intersection `TROP2` ∩ `CLDN4` = 0.

## 2. TROP2-ADC proteomics, tight junction down — NONE

Searched Europe PMC, PRIDE v2, and ProteomeXchange for sacituzumab govitecan, datopotamab / DS-1062, SKB264 / sacituzumab tirumotecan / MK-2870 plus proteome / mass spectrometry / claudin / tight junction.

**No hit** that is (ADC treatment) AND (global or phospho-proteomics) AND (TJ / claudin decrease).

Near-misses (not this analog):

| Item | Why it fails |
| --- | --- |
| Zhao et al. *JITC* 2026, 10.1136/jitc-2025-012265, PMID 41932810 | TROP2 **KO** and **hRS7** (SG antibody, not the ADC) lower **claudin-7** / occludin by IF/IHC in TNBC. Not SG/Dato/SKB264 proteomics; not CLDN4. |
| ASCO 2026 abstract e13006, 10.1200/jco.2026.44.16_suppl.e13006 | Paired pre/post-**SG** MS (mPROBE) of Trop2, TOPO1, HER2/3, Nectin4, SLFN11. **No TJ/claudin panel reported.** |
| PRIDE / ProteomeXchange ADC-name queries | 0 projects titled or keyworded sacituzumab / datopotamab / SKB264 / govitecan / DS-1062. |
| PXD039272 | Trop2 xenograft prostate proteome (not an ADC-treatment TJ study). |

## 3. Co-IP or AlphaFold2 TROP2–CLDN4 — NONE

| Test | Result |
| --- | --- |
| Co-IP TROP2–CLDN4 | **None positive.** Nakatsukasa 2010 is an explicit **negative** Co-IP. IntAct has Co-IP only for TACSTD2–CLDN1 and TACSTD2–CLDN7 from that paper; CLDN4 is curated as **PLA proximity** (EBI-9302237, MI score 0.27), not Co-IP. BioPlex 2021 (Huttlin, PMID 33961781) reports tagged Co-IP TACSTD2–**CLDN7**, not CLDN4. UniProt P09758 INTERACTION list does not include CLDN4. |
| AlphaFold2 / Multimer / AF3 of TROP2–CLDN4 | **None.** Europe PMC `(TROP2 OR TACSTD2) AND (CLDN4 OR claudin-4) AND (AlphaFold…)` returned no TROP2–CLDN4 complex paper. Published AF2/AF3 uses for TROP2 are nanobody / dimer models, not CLDN4. AF2 papers on **CLDN4–CPE** are the wrong partner. |

Vendor pages that list TACSTD2 as a CLDN4 “interactor” (e.g. Creative BioMart) are not publications and were not counted.

## 4. Cell-line IC50 vs TROP2 — FOUND

Public preclinical papers report IC50/GI50 for TROP2-ADCs on lines with measured TROP2 (flow / IHC), including TROP2-low/negative comparators.

| ADC | Paper | DOI | PMID | IC50 / GI50 vs TROP2 (as stated) |
| --- | --- | --- | --- | --- |
| Dato-DXd | Okajima 2021 *Mol Cancer Ther* | 10.1158/1535-7163.mct-21-0206 | 34413126 | TROP2-high GI50 **0.48–7.8 nM**; TROP2-low LK-2, Calu-6 **>100 nM** (Table 1) |
| Dato-DXd | Bellone 2025 *Cancer Res Commun* | 10.1158/2767-9764.crc-25-0057 | 40299780 | TROP2 3+ USC IC50 **0.11 µg/mL** vs control ADC 30–49 µg/mL; TROP2 0/1+ not selective |
| SKB264 | Cheng 2022 *Front Oncol* | 10.3389/fonc.2022.951589 | 36620535 | IC50 **1.281–18.83 nM** on TROP2+ lines; NCI-H23 TROP2+ > parental TROP2− |
| SG / IMMU-132 | Cardillo 2011 *Clin Cancer Res* | 10.1158/1078-0432.ccr-10-2939 | 21372224 | Cytotoxicity across Trop-2+ epithelial lines (later cited IC50 ~2–23 nM) |
| SG | Goldenberg 2015 *Oncotarget* | 10.18632/oncotarget.4318 | 26101915 | Trop-2 target characterization + efficacy in Trop-2+ models |
| SG | Lopez 2020 *Front Oncol* | 10.3389/fonc.2020.00118 | 32117765 | Trop-2+ EOC more sensitive than Trop-2− OVA14 vs control ADC |
| SG | Zeybek 2020 *Sci Rep* | 10.1038/s41598-020-58009-3 | 31969666 | Trop-2+ cervical IC50 **0.18–0.26 nM**; Trop-2− ADX-2 not selective |
| SG | Hoff 2022 *Cancers* | 10.3390/cancers14194789 | 36230712 | Eso26 (TROP2-high) IC50 **~6.5 nM** vs FLO-1 (TROP2-neg) **~100 nM** |
| SG | Cardillo 2019 *Oncotarget* | 10.18632/oncotarget.27766 | 33196706 | MDA-MB-231 Trop-2 transfection clones C13/C39 vs parental |
| SG | Liu 2023 *npj Breast Cancer* | 10.1038/s41523-023-00573-8 | 37567892 | TROP2 overexpression **lowers SG IC50** in BCX lines |

Payload resistance can dominate among TROP2+ lines (e.g. AACR 2025 abs. 2930, pancreatic Dato-DXd). That does not cancel the published TROP2-stratified IC50 tables above.

## Private claims (out of scope)

C1–C10 remain private. See `excluded_private.tsv`. This hunt does not quote, infer, or score those claims.

## Files

| File | Contents |
| --- | --- |
| `VERDICT.txt` | One-screen found vs none |
| `found_vs_none.tsv` | Same table, machine-readable |
| `found.tsv` | Qualifying public analogs |
| `none.tsv` | Lanes with no public analog |
| `near_misses.tsv` | Related papers that fail the analog as written |
| `excluded_private.tsv` | C1–C10 marked private |
| `search_log.md` | Query-by-query log |
