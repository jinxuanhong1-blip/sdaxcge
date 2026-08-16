# B7 search log

Date: 2026-08-16

## Queries

Europe PMC (JSON API, `resultType=lite`):

1. `(TITLE:"claudin-4" OR TITLE:"claudin 4" OR TITLE:CLDN4 OR TITLE:Cldn4) AND (TEER OR "transepithelial electrical resistance" OR "transepithelial resistance" OR TER)` — 62 hits (many gut/skin/kidney).
2. `(ABSTRACT:"claudin-4" OR ABSTRACT:CLDN4) AND (ABSTRACT:TEER OR ABSTRACT:"transepithelial electrical resistance") AND (lung OR alveolar OR airway OR pulmonary OR A549 OR H441 OR Calu-3 OR NHBE OR NSCLC)` — 47 hits; most mention CLDN4 as a TJ marker without CLDN4 perturbation.
3. `(TITLE:"claudin-4" OR TITLE:CLDN4) AND (lung OR alveolar) AND (barrier OR TEER OR TER OR permeability)`
4. Author-targeted: Wray; Mitchell + alveolar; Kage + knockout; Van Itallie + conductance.
5. `CPE OR "C-CPE" AND claudin-4 AND (TEER OR TER) AND (lung OR cancer)`

Web / PMC HTML / OA PDF follow-up for:

- Wray 2009 PMC2742793
- Mitchell 2011 PMC3129905 (Europe PMC PDF)
- Kage 2014 PMC4187039
- Rokkam 2011 PMC3157178
- Lee 2018 PMC5705480 (e-aair PDF)
- Van Itallie 2001 JCI PDF
- Litkouhi 2007 PMC1854850
- Nishiguchi 2019 PMC6481338
- Kuwada 2019 PMC6825989
- Schlingmann 2015 PMC4562902 (review)
- Overgaard 2012 PMC3375852 (review)
- Arabi 2024 PMC11446878 (lung-cancer review)
- Fujiwara-Tani 2023 PMC10051602 (cancer-therapy review)

Additional targeted web searches:

- `"claudin-4" siRNA TEER (A549 OR H441 OR Calu-3 OR NSCLC)`
- Rokkam alveolar fluid clearance
- Kage knockout TEER

## Inclusion

Paper is **in** if it (a) perturbs or scores CLDN4 and (b) reports a
barrier endpoint, or is a review that quotes those numbers.

## Exclusion / near-misses (not counted as lung-cancer TEER hits)

| Item | Why not a lung-cancer CLDN4-TEER hit |
| --- | --- |
| A549 / H441 / Calu-3 TGF-β TEER papers | TEER of lung lines without CLDN4 KD/OE |
| Park 2021 COPD (PMC8502106) | Plasma/sputum CLDN4 vs lung function; no TEER |
| Zhang 2024 Clin Transl Med MPE scRNA | CLDN4+ cells in recurrent MPE; no TEER |
| Moldvay / Soini NSCLC IHC | Expression only |
| PIKfyve YM201636 NSCLC (PMC7877718) | Mentions TEER conceptually; no CLDN4 TEER |
| Sun 2024 PRRSV endothelial CLDN4/8 | Porcine microvascular endothelium, not epithelium |
| Chen 2017 C2-ceramide alveolar TEER | CLDN4 protein down with ceramide; not a CLDN4 gene perturbation; OA PDF not retrieved here |
| Shlyonsky 2005 H441 dexamethasone | CLDN4 mRNA + dome / *R*t context; not a CLDN4 KD/OE TEER curve |
| Tokumasu 2017 MDCK CLDN4 KO (PLoS ONE) | Kidney line; CLDN4 dispensable in WT MDCK |
| Watari 2017 thiostrepton (Sci Rep) | Caco-2 TER, not lung |

## Lung-cancer TEER gap

No public primary paper was located that reports TEER/TER after CLDN4
siRNA, shRNA, KO, or overexpression in a lung-cancer cell line.
Arabi et al. 2024: “To our knowledge, there have been no studies
examining the mechanistic role of claudin-4 in LUAD development.”

## Number rules

- Copy printed means, SEM/SD, *n*, and *P* strings only.
- Figure-only TEER without a printed Ω value → `value_text=figure_only`.
- Reviews that restate a primary number are tagged `quoted_from_review`
  and the primary PMID is given.
- No *P* values were computed from figures or re-tested.
