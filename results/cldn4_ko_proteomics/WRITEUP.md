# CLDN4 knockout proteomics: ArrayExpress/BioStudies and PRIDE

Search date: 2026-09-21. Live counts are in `search_counts.tsv`. The script is `scripts/cldn4_ko_proteomics/search_cldn4_ko_proteomics.py`.

## Verdict

No public CLDN4 or Cldn4 knockout or knockdown proteome was found for a human epithelial line. MHC-I, type-I IFN, and type-II IFN proteins were not scored. The score table has 90 set-membership rows (82 gene symbols) and every row is `not_scored`, with `log2fc_ko_over_control` left as `NA`.

This is a repository result. It does not re-score the existing CLDN4-loss RNA series, and it does not convert a different proteomic contrast into a CLDN4-loss contrast.

## What counted as a hit

All three had to be true:

1. The perturbed gene is CLDN4 or Cldn4 (knockout, knockdown, CRISPR, siRNA, or shRNA).
2. The measurement is a proteome (protein abundance table), not RNA, methylation, metabolomics, or a western blot.
3. The cells are a human epithelial line. Any epithelial lineage was eligible.

## PRIDE

Quoted EBI Search queries `"CLDN4"`, `"Cldn4"`, `"claudin-4"`, `"claudin-4 knockdown"`, and `"claudin-4 knockout"` each returned 0 projects. The boolean query `CLDN4 AND (knockout OR knockdown OR siRNA OR shRNA OR CRISPR)` also returned 0. The PRIDE v2 keyword search returned 0 for `CLDN4`, `Cldn4`, and `claudin-4`.

The unquoted query `CLDN4` returned 2 projects, PXD051838 and PXD051839. Both are mouse Foxf2 endothelial knockouts. Their descriptions name Cldn5, not Cldn4 (`pride_false_positive_tokens.tsv`). They were not scored.

The keyword `claudin` returned 30 projects (`pride_keyword_claudin_accessions.txt`). The closest records are a different gene, a subtype label, or an interactome:

| Accession | Why it is not a CLDN4-loss proteome |
|---|---|
| PXD066158 | CLDN6-deficient human germ-cell tumor proteome |
| PXD005292 | Glycoproteome of claudin-low breast lines |
| PXD031094 | Pan-claudin co-IP interactome |
| PXD003651 | CLDN3 in prostate-cancer exosomes |
| PXD010908 | Phosphorylation of claudin-11 |

OmicsDI query `CLDN4` returned 166 datasets. All 166 were fetched. Source counts: BioStudies literature 114, ENA project 22, GEO 20, ArrayExpress 7, iProX 1, BioStudies other 1, BioImages 1. PRIDE datasets fetched: 0.

## ArrayExpress / BioStudies

The ArrayExpress collection API query `CLDN4` returned 4 studies. None is a proteome (`CLDN4 AND proteomics` = 0):

| Accession | What it is |
|---|---|
| E-GEOD-22493 | Human ovarian SKOV3, CLDN4 siRNA versus overexpression, microarray RNA (GSE22493) |
| E-GEOD-50927 | Mouse lung Cldn4 knockout RNA-seq (GSE50927) |
| E-GEOD-60885 | Trophoblast DNA methylation; claudin-4 is a named gene |
| E-GEOD-84742 | Mouse colon differentiation RNA-seq; not a CLDN4 knockout |

EBI Search over `biostudies-arrayexpress` returned 7 for the same word. The extra three (E-GEOD-16876, E-GEOD-4328, E-GEOD-8337) do not contain CLDN4, Cldn4, or claudin-4 in the title or description returned with the hit. They are unrelated RNA arrays.

BioStudies query `CLDN4` returned 116 records, mostly literature mirrors. The one deposited study with a protein file is S-BSST1967.

## S-BSST1967 is not a CLDN4-loss proteome

S-BSST1967 (human, released 2025-05-28) is the BioStudies deposit for the CLDN4 palmitoylation paper (PMC12281411, doi:10.1016/j.xcrm.2025.102208). Its files include `ExpAll.xlsx` (258,294 bytes, file description "protein omics") and an RNA-seq DESeq2 table.

The paper's mass-spectrometry experiments are:

- proteomics of lenvatinib-tolerant versus mock primary HCC cells
- label-free palmitoyl-proteomics of membrane fractions
- LC-MS/MS of CLDN4 immunoprecipitates

CLDN4 shRNA in MHCC97H and PLC/PRF/5 is a functional knockdown, not the proteomic contrast. The paper's data-availability statement deposits single-cell RNA at CNGB CNP0006650 and does not list a PRIDE accession. `ftp.ebi.ac.uk` did not serve `ExpAll.xlsx` in this session (FTP 421, HTTPS handshake closed), so MHC/IFN proteins were not read from that workbook. They would still be the wrong contrast if the workbook is the lenvatinib or palmitoyl table.

## Human epithelial CLDN4 loss that stayed off the proteome archives

These are real CLDN4 perturbations. They do not add a protein matrix.

| Record | Line | Measurement |
|---|---|---|
| E-GEOD-22493 | SKOV3 ovarian | microarray RNA |
| GSE207704 | MCF7 and T47D breast | CRISPR RNA-seq; no PRIDE project |
| PMID 41214101 (PMC12603150) | high-grade serous ovarian models | knockdown, ISRE reporter, immunoblot, flow; no PXD or GSE in the full text |
| PMID 35373300 | ovarian models | knockdown; DNA-repair immunoblots; no proteome accession found |
| PMID 38867360 (PMC11218812) | ovarian models | CRISPRi; UHPLC-MS metabolomics, not a proteome |
| PMID 41016339 | SCLC H1688 | knockout RNA-seq; no public matrix |

PubMed title/abstract query for CLDN4 plus knockout/knockdown/siRNA/shRNA/CRISPR plus proteome/"mass spectrometry" returned 2 papers: PMID 30353739 (Caco-2 treated with Bifidobacterium factors) and PMID 20511395 (dog MDCK plasma-membrane proteomics during EMT). Neither is a human CLDN4 knockout.

## MHC / IFN proteins

The panel is the same manual set used for the CLDN4/TACSTD2 knockdown RNA wave: `IFN_ALPHA_TYPE1`, `IFN_GAMMA_TYPE2`, and `ANTIGEN_PRESENTATION_MHC1`. Membership is written on every row of `mhc_ifn_protein_scores.tsv`.

No qualifying protein matrix was available, so no log2 fold-change, replicate count, or p-value was calculated. `qualifying_hits.tsv` has a header and no data rows.
