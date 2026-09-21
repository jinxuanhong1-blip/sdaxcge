# Search log

Date: 2026-09-21. Window: first publication or meeting year 2023 through 2026. Question: AACR abstracts and bioRxiv preprints that perturb CLDN4 by CRISPR, CRISPRi, shRNA, or siRNA and that deposit omics.

## Queries that returned the screened set

Europe PMC, synonym off:

- Preprints with CLDN4 / Cldn4 / claudin-4 in the title, `SRC:PPR`, `FIRST_PDATE:[2023-01-01 TO 2026-12-31]`. Eight hits. Listed in `screened.tsv`.
- Preprint abstracts containing CLDN4 and CRISPR, CRISPRi, knockdown, knockout, siRNA, or shRNA in that window. Five hits. Two are Research Square (Connexin26; not bioRxiv). One is a claudin-3 knockout mouse preprint. One is the kidney-stone CLDN4 variant. One is `10.1101/2024.01.18.576263`.
- Exact phrases `"CLDN4 knockdown"`, `"CLDN4 knockout"`, `"CLDN4 shRNA"`, `"Cldn4 knockout"`, `"Cldn4 knockdown"`, `"claudin-4 knockdown"`, `"claudin-4 knockout"`, `"claudin-4 CRISPR"`, `"CLDN4 CRISPR"` over the same dates, all sources. Used to find journal papers that might also exist as bioRxiv. `"CLDN4 CRISPR"` returned 0.

NCBI:

- GEO DataSets: `(CLDN4 OR Cldn4 OR claudin-4) AND (CRISPR OR knockdown OR knockout OR siRNA OR shRNA) AND gse[ETYP]`. Series: GSE207704, GSE84742, GSE60885, GSE22493. Only GSE207704 is a CLDN4 knockout RNA-seq series, and it is outside this wave’s venues.
- GEO title search for CLDN4 series with publication date 2022–2026: GSE207704 only.
- BioProject: CLDN4 plus CRISPR/knockdown/knockout/siRNA/shRNA, 2023–2026: 0.
- PubMed: CLDN4 in title/abstract plus CRISPR/knockdown/knockout in Cancer Research, Clinical Cancer Research, Molecular Cancer Therapeutics, or Cancer Research Communications, 2023–2026: PMID 38867360 (the autophagy paper). Meeting abstracts with DOI prefix `10.1158/1538-7445` were not in this PubMed result.

Repositories:

- PRIDE archive keyword `claudin-4`: empty project list.
- Metabolomics Workbench `study_title/claudin/summary`: `[]`.
- MetaboLights study search `claudin-4`: no study.
- Figshare API collection 7312139 and article search `Claudin-4 Stabilizes the Genome`.

Publisher text:

- bioRxiv `10.1101/2024.01.18.576263`: CRISPRi methods; data availability says manuscript plus supplement; flow cytometry and data mining upon request.
- bioRxiv HTML `10.1101/2024.09.04.611120v1.full`: “All data is provided within the manuscript and supplemental information.”
- Europe PMC full text for PMC11218812, PMC11705808, PMC12603150, PMC10105442, PMC12281411, PMC13253352. Accessions recorded in `screened.tsv`.

## Limits

AbstractsOnline (`abstractsonline.com`) returned HTTP 403, so AACR annual-meeting and special-conference abstracts were taken from publisher pages and web results that quote the abstract, not from a complete program dump. A meeting abstract that never reached a DOI landing page or a search snippet could have been missed. Europe PMC did not index `10.1158/1538-7445.ovarian23-b072` or `10.1158/1538-7445.am2025-lb387`.

Compound counts and sample IDs were read from the figshare docx (`word/document.xml` table rows), not copied from a paper’s narrative.
