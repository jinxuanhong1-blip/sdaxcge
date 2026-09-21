# Search log — lung-cancer CLDN4 KD match wave

Date: 2026-09-21.

Question: which 2018–2026 lung-cancer papers used CLDN4 siRNA, shRNA, or CRISPR, which of those deposited omics in GEO, SRA, GSA, NODE, or CNCB, and, if a new accession exists, download it and score IFN.

## Paper screen

PubMed E-utilities, publication date 2018–2026. Exact phrases included:

`knockdown of CLDN4`, `knockdown of claudin-4`, `CLDN4 knockdown`, `claudin-4 knockdown`, `CLDN4 siRNA`, `CLDN4 shRNA`, `siCLDN4`, `CLDN4 knockout`, `Cldn4 knockout`, `silencing of CLDN4`, `CLDN4 was silenced`, `CLDN4 was knocked`, `knockout of CLDN4`, `knockout of claudin-4`, `deletion of CLDN4`, `CLDN4 CRISPR`, `siRNA targeting CLDN4`, `shRNA targeting CLDN4`, `siRNA against CLDN4`, `shRNA against CLDN4`, `interfering CLDN4`, `CLDN4 gene silencing`.

Hyphenated queries that NCBI parsed as Boolean NOT (`CLDN4-specific`, `CLDN4−/−`) were discarded. After that filter, 101 PMIDs remained. Abstracts were fetched and flagged if they mentioned lung, NSCLC, SCLC, A549, H1299, H1688, Calu-3, pulmonary, bronchial, or alveolar tissue. Those full texts that are in PMC were opened on Europe PMC.

Europe PMC full-text proximity (`CLDN4` / `claudin-4` within 15 words of siRNA, shRNA, knockdown, knockout, silencing, CRISPR, or sgRNA, plus a lung term) was used as a second net. Co-occurrence hits were lung injury, other claudins, or citations of non-lung knockdowns, except the papers in `papers.tsv`.

## Deposit screen

| Source | Query | Result |
|---|---|---|
| GEO gds | `CLDN4 AND (knockdown OR knockout OR siRNA OR shRNA OR CRISPR) AND gse`, limited to 2018–2026 | 1 hit: GSE207704 (breast) |
| GEO gds | `CLDN4 AND (lung OR NSCLC OR A549) AND (knockdown OR siRNA OR shRNA OR CRISPR OR knockout) AND gse` | 0 |
| GEO gds | `H1688 AND (CLDN4 OR SAA1)` | 0 |
| GEO gds | `CRAD AND (A549 OR NSCLC OR lung)` | 0 |
| GEO gds | `Yazawa[Author] AND (CLDN4 OR lung OR SCLC)` | GSE310370, GSE281782, GSE233214 (vasculitis spatial). None is a CLDN4 knockout |
| SRA | `H1688 AND CLDN4` | 0 |
| SRA | `H1688 AND (knockout OR knockdown OR CRISPR OR SAA1 OR CLDN4)` | shKLF9 only (GSM9043563, GSM9043564) |
| SRA | `Dokkyo AND (CLDN4 OR SAA1 OR H1688)` | 0 |
| BioProject | `CLDN4 AND (knockout OR knockdown OR CRISPR) AND (lung OR SCLC OR H1688)` | 0 |
| OmicsDI | `CLDN4 AND (knockout OR knockdown OR siRNA OR CRISPR) AND (lung OR SCLC OR NSCLC)` | E-GEOD-50927 and literature mirrors, no lung-cancer KD |
| OmicsDI | `H1688 AND CLDN4` | 0 |
| OmicsDI | `CLDN4 AND source:node` | 0 |
| GSA HTML search | `CLDN4` | 11 items: GSE50927, GSE207704, gastric 4C-seq at the CLDN4 locus |
| GSA HTML search | `H1688` | 19 items. Relevant RNA-seq is PRJNA1005054 (oridonin) plus shKLF9 and a methylation sample |
| CNCB `/search/api/specific?db=gsa-human&q=CLDN4` | GSA-Human | recordsTotal 0 |
| CNCB same API `db=gsa` | GSA | 11, same set as the HTML search |
| CNCB same API `db=bioproject` | BioProject | 16. Lung-related projects are observational, Cldn18/CEBPA, or GSE50927. No H1688 CLDN4 knockout. PRJCA023797 is CLDN4-positive effusion scRNA |
| CNCB same API `db=biosample` | BioSample | 5: MCF7 CLDN4-/- and Cldn4 KO VILI samples only |
| NODE | `https://www.biosino.org/node/api/public/search?query=CLDN4` and `/api/v1/search` | HTTP 200 body `404 NOT_FOUND`. Public search HTML is a shell with no embedded accessions |

PRJNA1005054 was opened. Six H1688 RNA-seq runs (S1_1–S1_3, S2_1–S2_3). Project description is oridonin antitumor activity in SCLC. Not scored.

## IFN

Not run. No new lung-cancer CLDN4-perturbation matrix was available to download.
