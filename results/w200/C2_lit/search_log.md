# Search log — C2 public literature (CLDN4 lysosomal degradation after ADC)

Date: 2026-08-16. Catalog only. Private C2 text was not used.

Analog as written in the request: **CLDN4 lysosomal degradation after ADC**.

A row is FOUND only if a public paper shows all three: (1) an ADC, (2) CLDN4 protein loss or lysosomal colocalization, (3) lysosome dependence (chloroquine / bafilomycin / LAMP1-CLDN4, or equivalent). Near-misses fail at least one of those.

## Europe PMC (`https://www.ebi.ac.uk/europepmc/webservices/rest/search`)

| Query | hitCount | Qualifying analog |
| --- | ---: | --- |
| `CLDN4 AND (lysosom* OR "lysosomal degradation") AND (ADC OR "antibody-drug" OR sacituzumab OR datopotamab OR SKB264)` | 19 | none; reviews and wrong-claudin ADCs (CLDN6, CLDN2, EpCAM/CLDN3) |
| `(CLDN4 OR "claudin-4") AND (lysosom*) AND (internali* OR endocyt* OR degrad*)` | 602 | token overlap; no ADC + CLDN4 protein-fate paper in the first pages. Starvation/DON/EGF papers pulled as near-misses. |
| `(TROP2 OR TACSTD2 OR Trop-2) AND (CLDN4 OR "claudin-4") AND (lysosom* OR degrad* OR internali*)` | 116 | no TROP2-ADC CLDN4 lysosome experiment. Sun 2026 RNA co-expression only. |
| `("sacituzumab govitecan" OR datopotamab OR SKB264 OR "sacituzumab tirumotecan") AND (CLDN4 OR "claudin-4" OR claudin) AND (lysosom* OR degrad*)` | 143 | Zhao JITC 2026 (CLDN7 / hRS7); ADC reviews. No SG/Dato/SKB264 + CLDN4 lysosome. |
| `(CLDN4 OR "claudin-4") AND (chloroquine OR bafilomycin OR "lysosomal inhibitor") AND (degrad*)` | 110 | Ikari 2019 chloroquine raises CLDN4 (constitutive). No ADC. |
| `(TROP2 OR TACSTD2) AND (claudin) AND (lysosom*)` | 28+ | Wu 2020 Cells (CLDN1/7 + chloroquine) is the closest TROP2–claudin–lysosome paper. |
| `(CLDN4 OR "claudin-4") AND (CPE OR "C-CPE" OR enterotoxin) AND (lysosom* OR endocyt* OR internali*)` | 254 | Saeki 2011 CPE-rGel (cargo trapped); Hashimoto 2013 CPE-ETA'. Ligand, not ADC. |
| `hRS7 AND (claudin OR CLDN)` | 32 | Zhao JITC 2026 is the only primary hRS7–claudin experiment (CLDN7). |
| `(CLDN4 OR "claudin-4") AND ("antibody-drug conjugate" OR ADC) AND (internali* OR lysosom*)` | 57 | no CLDN4-ADC protein-fate paper. CLDN6-23-ADC and reviews. |
| `AUTHOR:Wu AND matriptase AND (EpCAM OR TROP2) AND claudin` | 6 | Wu 2017 JCI; Wu 2020 Cells. |
| `TITLE:(CLDN4 OR "claudin-4") AND (lysosom* OR endocyt* OR internali* OR trafficking OR degrad*)` | 63 | Li 2021 FASEB/DON; Ogawa 2012; Cong 2015 (via related queries). |
| `TITLE:("endocytic regulation of intestinal tight junction")` | 1 | Li 2021 FASEB J, PMID 33484473. |
| `TITLE:(deoxynivalenol) AND TITLE:(endocytosis) AND claudin-4` | 4 | Li 2021 Arch Toxicol, PMID 33847777. |
| `DOI:10.1242/jcs.165878` | 1 | Cong 2015, PMID 25948584. |
| `TITLE:(Epidermal growth factor modulates claudins) AND AUTHOR:Ogawa` | 1 | Ogawa 2012, PMID 22544349. |
| `TITLE:(Immunotoxin-mediated targeting of claudin-4)` | 1 | Hashimoto 2013, PMID 23563899. |
| `("CPE") AND gelonin AND claudin` | 5 | Saeki/Yuan 2011 BMC Cancer, PMID 21303546. |
| `("antibody-drug conjugate") AND CLDN4 AND (MMAE OR SN-38 OR DXd OR payload)` | 17 | CLDN6 ADC and reviews; no CLDN4-ADC paper. |

DOI lookups confirmed on Europe PMC or publisher pages:  
`10.3390/cells9041027`, `10.1172/jci88428`, `10.1136/jitc-2025-012265`, `10.1158/1535-7163.mct-21-0206`, `10.1158/1078-0432.ccr-10-2939`, `10.18632/oncotarget.4318`, `10.2353/ajpath.2010.100149`, `10.1371/journal.pone.0117074`, `10.1242/jcs.165878`, `10.1096/fj.202002098r`, `10.1007/s00204-021-03044-w`, `10.1007/s00418-012-0956-x`, `10.1038/s41598-019-46250-4`, `10.1186/1471-2407-11-61`, `10.3892/ijo.2013.1881`, `10.1158/1078-0432.ccr-22-2981`.

ASCO 2026 `10.1200/jco.2026.44.16_suppl.e13006` did not index in Europe PMC; counted from the sibling public-analog hunt (no TJ/claudin panel).

## PRIDE / ProteomeXchange

Same empty ADC-name result as `results/claim_C_private`: no sacituzumab / datopotamab / SKB264 treatment proteome that reports claudin / TJ decrease. Keyword `CLDN4 ADC` on PRIDE WS v2 returned an empty list.

## Web / publisher pages inspected

- Wu 2020 Cells PMC7226414: chloroquine rescues claudin-1 and claudin-7 after HAI-1/2 knockdown; cleaved TROP2/EpCAM accumulate. CLDN4 is not in that rescue figure.
- Zhao 2026 JITC HTML: hRS7 lowers claudin-7; combination with anti-PD-1. Not SG/Dato/SKB264; not CLDN4.
- Okajima 2021 PMC9398094: Dato-DXd lysosomal DXd release; GI50 vs TROP2. No CLDN4 blot.
- Nakatsukasa 2010: PLA TACSTD2–CLDN4 positive; Co-IP CLDN4 negative.
- Cong 2015 JCS abstract/publisher extract: CLDN4 S195 → β-arrestin2 → clathrin; ubiquitin-dependent degradation after carbachol.
- Li 2021 FASEB abstract: lysosome-dependent loss of CLDN1/3/4 in starved IPEC-J2; dynamin-dependent for CLDN3/4.
- Saeki 2011 BMC Cancer PDF extract: CPE-rGel enters vesicles and is degraded/trapped in the endosomal/lysosomal compartment unless R9 is added.
- Hashimoto 2013: CPE-ETA' cartoon is LRP1/ETA translocation (Golgi/ER), not lysosomal ADC release.
- Patent WO2016175551: CLDN4 antibody–dendron–doxorubicin; claimed acidic endosome/lysosome release. Not counted as a paper.

## What was not counted as FOUND

- Reviews that restate ADC lysosomal payload release in general.
- CLDN6 / CLDN18.2 / CLDN2 ADCs.
- Vendor pages and patents.
- RNA co-expression of TACSTD2 and CLDN4.
- CLDN4 endocytosis after starvation, DON, EGF, or carbachol.

## Conclusion of the log

- Direct analog (CLDN4 lysosomal degradation after ADC): **NONE**.
- Closest public mechanism papers are Wu 2020 (TROP2/claudin lysosome, wrong claudin) and the CLDN4 endocytosis set (right protein, no ADC).
- TROP2-ADC papers document lysosomal trafficking of **TROP2/ADC**, not of **CLDN4**.
