# Search log — public analogs of private C1–C10

Date: 2026-08-16. Private claim text was not used. Four analog classes only.

## Europe PMC (`https://www.ebi.ac.uk/europepmc/webservices/rest/search`)

| Query | hitCount | Qualifying analog |
| --- | ---: | --- |
| TITLE:(TROP2 OR Trop-2 OR TACSTD2) AND TITLE:(CLDN4 OR "claudin-4" OR "claudin 4") | 0 | none |
| ABSTRACT:(TROP2 OR Trop-2 OR TACSTD2) AND ABSTRACT:(CLDN4 OR "claudin-4" OR "claudin 4") | 13 | no dual-IF coloc paper; Chen 2024 scRNA MPE; patents; EMT cell-line paper |
| (TROP2 OR Trop-2 OR TACSTD2) AND (CLDN4 OR "claudin-4") AND (colocali* OR PLA OR "dual immunofluorescence" OR multiplex) | 70 | Nakatsukasa 2010 (PMID 20651236) is the PLA hit; Van Itallie 2015 (PMID 25789658) is BioID; rest are token-overlap |
| (TROP2 OR Trop-2 OR TACSTD2) AND (CLDN4 OR "claudin-4") AND (immunoprecipit* OR "co-IP" OR interactome OR pulldown) | 42 | Nakatsukasa 2010 **negative** Co-IP for CLDN4; no positive Co-IP |
| (TROP2 OR TACSTD2) AND (AlphaFold OR AlphaFold2 OR AlphaFold-Multimer) AND (CLDN4 OR claudin) | 24 | no TROP2–CLDN4 complex |
| (TROP2 OR TACSTD2) AND (CLDN4 OR "claudin-4") AND (AlphaFold OR "predicted structure" OR multimer) | 5 | none are TROP2–CLDN4 AF |
| ("sacituzumab govitecan" OR datopotamab OR DS-1062 OR SKB264 OR "sacituzumab tirumotecan" OR MK-2870) AND (proteom* OR "mass spectrometry" OR phosphoproteom*) | 635 | reviews / PK LC-MS; no ADC-treatment TJ proteome |
| ("sacituzumab govitecan" OR datopotamab OR SKB264 OR "sacituzumab tirumotecan" OR DS-1062) AND (claudin OR "tight junction" OR occludin OR CLDN) | 264 | Zhao JITC 2026 (CLDN7 / hRS7); reviews; no ADC MS TJ-down |
| ADC names AND (PRIDE OR ProteomeXchange OR "proteomic profiling") | 54 | no SG/Dato/SKB264 TJ dataset |
| ("sacituzumab govitecan" AND (IC50 OR "half maximal") AND (TROP2 OR Trop-2) AND (cell OR line)) | 447 | multiple SG IC50 papers (see found.tsv) |
| (datopotamab AND (IC50 OR "IC 50") AND (TROP2 OR Trop-2)) | 60 | Okajima 2021; Bellone 2025; AACR abstracts |
| ((SKB264 OR "sacituzumab tirumotecan" OR MK-2870) AND (IC50 OR "IC 50") AND (TROP2 OR Trop-2)) | 36 | Cheng 2022 is the primary preclinical IC50 source |
| (GDLD OR "gelatinous drop") AND (TACSTD2 OR TROP2) AND claudin | 27 | Nakatsukasa 2010 primary; later GDLD functional papers |

DOI lookups confirmed: 10.2353/ajpath.2010.100149, 10.1371/journal.pone.0117074, 10.1158/1535-7163.mct-21-0206, 10.3389/fonc.2022.951589, 10.1158/1078-0432.ccr-10-2939, 10.18632/oncotarget.4318, 10.3389/fonc.2020.00118, 10.1038/s41598-020-58009-3, 10.3390/cancers14194789, 10.18632/oncotarget.27766, 10.1038/s41523-023-00573-8, 10.1158/2767-9764.crc-25-0057, 10.1136/jitc-2025-012265.

ASCO abstract 10.1200/jco.2026.44.16_suppl.e13006 did not index in Europe PMC; abstract text was read from the JCO page (Trop2/TOPO/HER/Nectin4 MS only).

## IntAct (`/intact/ws/interaction/findInteractions/TACSTD2`)

Claudin-related binary records:

| Intact ac | Pair | Method | PMID | Note |
| --- | --- | --- | --- | --- |
| EBI-9302190 | TACSTD2–CLDN7 | anti-bait Co-IP | 20651236 | positive |
| EBI-9302190 | TACSTD2–CLDN1 | anti-bait Co-IP | 20651236 | positive |
| EBI-9302231 | CLDN7–TACSTD2 | PLA | 20651236 | proximity |
| EBI-9302215 | CLDN1–TACSTD2 | PLA | 20651236 | proximity |
| EBI-9302237 | **CLDN4–TACSTD2** | **PLA only** | 20651236 | MI score 0.27; **not Co-IP** |
| EBI-54945048 | CLDN7–TACSTD2 | anti-tag Co-IP | 33961781 | BioPlex; not CLDN4 |

## UniProt P09758 (TACSTD2)

INTERACTION comment lists KRTAP1-3, KRTAP10-8, KRTAP5-2, LCE2C. **No CLDN4.**

## PRIDE / ProteomeXchange

| Store | Query | Result |
| --- | --- | --- |
| PRIDE WS v2 `/search/projects?keyword=` | sacituzumab, datopotamab, SKB264, govitecan, "TROP2 ADC", DS-1062 | 0 projects each |
| ProteomeXchange `GetDataset?filterstr=` | sacituzumab, datopotamab | `[]` |
| ProteomeXchange | TROP2 | unrelated / no-publication datasets (e.g. PXD075364); not ADC TJ-down |

## Web / publisher pages inspected

- Nakatsukasa 2010 full text (Semantic/AJP extract): Fig 2A negative Co-IP; Fig 2B–C positive PLA including CLDN4; discussion paragraph on proximity vs binding.
- Van Itallie 2015 PLOS HTML: Table 2 + paragraph naming TROP-2 around Cldn4.
- Okajima 2021 PMC9398094: Table 1 GI50 vs TROP2 copies/cell.
- Cheng 2022 PMC9817100: Fig 2D IC50 range; TROP2+ vs parental NCI-H23.
- Zhao 2026 JITC PDF/HTML: Co-IP and IF are **claudin 7**; hRS7 not SG/Dato/SKB264 proteomics.
- JCO 2026 e13006 abstract: SG paired MS without TJ proteins.
- Creative BioMart CLDN4 page: vendor interactor list; discarded.

## Conclusion of the log

- Coloc lane: **FOUND** (PLA + BioID). Not a lung dual-IF table.
- ADC proteomics TJ-down: **NONE**.
- Co-IP or AF2 TROP2–CLDN4: **NONE** (Co-IP published negative; AF2 absent).
- IC50 vs TROP2: **FOUND** for SG, Dato-DXd, and SKB264.
