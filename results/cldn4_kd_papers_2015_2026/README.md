# CLDN4 siRNA / shRNA / CRISPR in cancer, 2015–2026

Master table: `master_table.tsv` (32 papers). Non-cancer experiments that used the same perturbations are in `appendix_noncancer.tsv` and are not counted above.

`data_pointer` is the place to get the CLDN4 perturbation data:

- a repository accession, only when that accession is the CLDN4 knockdown or knockout profile
- otherwise `SUPPLEMENT:` plus the PMCID or article DOI, when the journal files hold the experiment
- otherwise `UNDEPOSITED`

`cldn4_perturbation_profile` is `DEPOSITED:GSE207704` or `UNDEPOSITED`. `also_deposited` lists other accessions from the same paper. Those are not CLDN4-knockdown matrices.

## What was deposited

One cancer CLDN4-knockout expression matrix falls in this window.

| Pointer | Paper | What it is |
|---|---|---|
| GSE207704 (BioProject PRJNA856719) | Murakami-Nishimagi 2023, PMID 37059993 | RNA-seq of T47D and MCF-7, parental vs CLDN4-/-, two replicates each (8 samples) |
| PRJNA1424384 and 10.5281/zenodo.18557106 | White 2026, PMID 42020735 | Genome-wide CRISPR screen in which claudin-4 was the hit. HT29 CLDN4-knockout phenotyping is in the paper, not a knockout transcriptome |

No NODE (OEP/OEX) accession was found. CNP0006650 (PMID 40592346) is single-cell RNA-seq of hepatocellular carcinoma, not the CLDN4 knockdown.

Lung and ovarian knockdown/knockout profiles in this window are UNDEPOSITED. The SCLC paper (Kashiwagi 2025, PMID 41016339) describes RNA-seq of CLDN4-knockout H1688 cells and is not linked to GEO or SRA. The ovarian shRNA and CRISPRi series (Hicks/Breed, Yamamoto, Neville, Villagomez, Bitler; PMIDs 27724921, 30606772, 35373300, 36237976, 38867360, 39625235, 41214101) are journal supplements or “upon request”.

## Priority counts

| Group | Papers | Perturbation profile deposited |
|---|---|---|
| lung | 2 | 0 |
| ovarian | 7 | 0 |
| breast | 5 | 1 (GSE207704) |
| other cancer | 18 | 0 |

PMID 37427351 is included under other cancer with an explicit note: the CRISPR edit is the Twist1 E-box in the CLDN4 promoter, not a coding-sequence knockout.

PMID 39502213 uses A549 and a claudin-4 siRNA, but the disease model is LPS acute lung injury. It is in the appendix, not the cancer table.

## Outside this window

These are real CLDN4-loss expression sets and they are not 2015–2026 cancer papers:

- GSE22493 (public 2010): microarray of CLDN4 siRNA versus CLDN4-high SKOV-3 ovarian cells. No PMID on the GEO series.
- GSE50927 (PMID 25106430, public 2014): mouse lung RNA-seq, wild-type versus Cldn4 knockout, after ventilator-induced lung injury.

## How the list was built

Searches were run on 21 Sep 2026.

- Europe PMC, publication date 2015-01-01 to 2026-12-31, 81 exact phrases (`CLDN4 siRNA`, `claudin-4 knockdown`, `shCLDN4`, `knockdown of CLDN4`, and the same forms for knockout, silencing, depletion, and Cldn4). Hyphen-slash forms such as `CLDN4-/-` were dropped because the index expanded them into thousands of unrelated hits. Unique records: 139.
- PubMed title/abstract query for CLDN4 plus siRNA/shRNA/CRISPR/knockdown/silencing in the same years: 122 PMIDs, unioned with the phrase set.
- Full text for 124 PMC articles was mined for a CLDN4-directed perturbation sentence, the data-availability section, repository regexes (GEO, SRA, ENA, ArrayExpress, GSA, NODE, CNGB, figshare, Zenodo, Dryad), and `supplementary-material`.
- GEO E-utilities were searched for CLDN4 plus knockdown/knockout/siRNA/shRNA/CRISPR. PubMed-to-GEO links were checked for the closed papers.

A row is a paper that itself silenced or edited CLDN4 in a cancer model. Reviews that only cite such an experiment are not rows. Knockdown of a different gene that changes CLDN4 (for example PAK4, PMID 30808546; FOXA1, PMID 37397926; cingulin/FOXO1, PMID 38338691) is not a row.

Closed full text (no PMCID) was called UNDEPOSITED only after a negative PubMed–GEO link and, where a PDF could be opened, a negative scan of that PDF. A paywalled supplement can still exist for those PMIDs; the note on the row says so.
