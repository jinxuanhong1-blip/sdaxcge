# Search log — TROP2+CLDN4 lung IHC tables

Date: 2026-08-16. Sources limited to Figshare, Zenodo, and journal supplements for a **hit**. Other indexes used only to find candidate supplements.

## Zenodo API (`https://zenodo.org/api/records`)

| Query | Total (API) | Dual TROP2+CLDN4 lung H-score / coloc table |
| --- | --- | --- |
| TROP2 AND CLDN4 | 0 (strict AND, first pass) / 21 (OR-ish `TROP2+CLDN4`) | none |
| TACSTD2 AND CLDN4 | 0 | none |
| "Trop-2" AND claudin AND lung | 0 | none |
| TROP2 AND lung AND "H-score" | 0 | none |
| CLDN4 AND lung AND "H-score" | 0 | none |
| TROP2 AND CLDN4 AND immunohistochemistry | 0 | none |
| TACSTD2 AND claudin-4 AND IHC | 0 | none |
| TACSTD2 | 6 | none |
| TROP2 lung | 15560 (token-split; mostly unrelated "lung") | none meeting criteria |
| CLDN4 lung immunohistochemistry | 17217 (token-split) | none meeting criteria |

Inspected records that mentioned TROP2 and lung:

- 10.5281/zenodo.18543127 / 18494664 — TROP-2 + cMET H-scores; **restricted**; no CLDN4.
- 10.5281/zenodo.18941275 — `TROP2_paper_public_bundle.zip`; title/metadata do not claim CLDN4 lung IHC tables.
- Other TROP2 Zenodo hits: ADC reviews, EMPD, breast, ELISA methods — not lung dual IHC.

## Figshare API (`https://api.figshare.com/v2/articles/search`)

| Query | n returned (page_size 8–20) | Dual table |
| --- | --- | --- |
| TROP2 CLDN4 | 0 | none |
| TACSTD2 CLDN4 | 1 (colorectal cohesin supplement; gene-list mention) | no |
| TROP2 CLDN4 lung | 0 | none |
| TROP2 H-score lung | 0 | none |
| CLDN4 H-score lung | 0 | none |
| claudin-4 TROP2 immunohistochemistry | 4 (testicular germ-cell surface-protein tables) | not lung |
| Trop-2 claudin-4 NSCLC | 0 | none |
| CLDN4 lung immunohistochemistry | 0 | none |
| claudin-4 NSCLC H-score | 0 | none |
| claudin-4 lung H-score supplementary | 0 | none |
| TROP2 NSCLC supplementary H-score | 0 | none |
| TACSTD2 immunohistochemistry lung table | 0 | none |
| Moldvay claudin lung | 0 | none |
| Paschoud claudin lung | 0 | none |
| `:title: TROP2 OR TACSTD2 AND CLDN4` | 18 TROP2-related files | none paired with CLDN4 in lung |

Inspected Figshare objects:

- 10.6084/m9.figshare.19665135 — Dum et al. TROP2 TMA supplement: three PDFs of **figures**, not a case-level H-score table; no CLDN4.
- 10.3389/fonc.2025.1638054.s001 / .s002 — TROP2 lung review data sheets (bioinformatics / trial text), not dual IHC H-scores.
- 10.3389/fimmu.2025.1709316.s001 — GPC3 + TROP2 LUSC (n=10); **no CLDN4**.

## Europe PMC

| Query | hitCount | Dual lung IHC table |
| --- | --- | --- |
| TITLE:"TROP2" AND TITLE:"CLDN4" | 0 | none |
| ABSTRACT:(TROP2 AND CLDN4 AND (lung OR NSCLC) AND (IHC OR immunohistochemistry)) | 0 | none |
| ABSTRACT:(TROP2 OR Trop-2 OR TACSTD2) AND ABSTRACT:(CLDN4 OR "claudin-4" OR "claudin 4") | 13 | scRNA MPE (Chen 2024); ovarian expression patents; no lung dual IHC table |
| TITLE:(TROP2 OR Trop-2 OR TACSTD2) AND TITLE:(CLDN4 OR claudin OR claudins) | 8 | breast TROP2/CLDN7; biliary CLDN18.2+Trop-2; ovarian CLDN6+Trop-2; mouse/cell claudin papers — **not lung TROP2+CLDN4 H-scores** |

## Journal supplements checked by paper

| Paper | Supplement / data statement | Dual table |
| --- | --- | --- |
| Ahmed et al. PLOS ONE 10.1371/journal.pone.0321555 | Supporting figs S1–S7 (Trop-2 / PD-L1 mRNA and demographics). Patient-level IHC on request to Gilead. | no |
| Dum et al. Pathobiology 10.1159/000522206 | Figshare 19665135: Suppl Fig 1–3 PDFs | no |
| Liu et al. J Transl Med 10.1186/s12967-025-06888-3 | Trop-2 H-score methods + in-article tables; data availability is article/supplement of **TROP2 only** | no |
| Báthori et al. 10.3389/pore.2023.1611328 | CLDN H-scores in rare lung tumors; no TROP2 | no |
| Sun et al. 10.1038/s41698-026-01523-w | RNA co-expression of TACSTD2 with CLDN4 | no |
| Zhao et al. 10.1136/jitc-2025-012265 | TROP2 / claudin-7 in breast | no |

## Negative API / web notes

- No Figshare or Zenodo record title/description claimed “TROP2 and CLDN4” H-scores or co-localization in lung.
- Tokenized Zenodo searches for “TROP2 lung” and “CLDN4 lung immunohistochemistry” return thousands of unrelated records; titles inspected on the first pages do not include a dual-marker H-score table.
- Dryad / OSF / Mendeley Data web hits for this query were immune/mIF or other markers, not TROP2+CLDN4.

## Conclusion of the log

Zero qualifying public tables. Result recorded as NONE.
