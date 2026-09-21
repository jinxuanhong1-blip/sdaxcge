# CLDN4 KD match wave (2024–2026 refresh)

Search date: **2026-09-21**.

**New public CLDN4/Cldn4 KD/KO/CRISPR/shRNA/siRNA expression datasets outside the catalog: 0.**

IFN/APM overlap was **not scored**. There is no new matrix to download. No overlap count, logFC, or p-value was imputed.

Catalog held out, not re-scored: **GSE22493, GSE50927, GSE207704, GSE274940**.

## What counted as a hit

A hit is a public expression profile (RNA-seq, microarray, or similar gene-expression matrix) in which CLDN4 or Cldn4 itself was knocked down, knocked out, or targeted by CRISPR, shRNA, or siRNA, released or first public in 2024–2026, and not one of the four catalog accessions.

Repositories searched: NCBI GEO (GDS), NCBI SRA, EBI BioStudies/ArrayExpress, NGDC GSA (including GSA-Human / HRA, CRA, OMIX, and BioProject). Europe PMC 2024–2026 was used only to find papers that might have deposited into those repositories.

## Repository result

| Source | Query shape | What came back |
|---|---|---|
| GEO | CLDN4/Cldn4 plus knockdown, knockout, CRISPR, shRNA, or siRNA | Only catalog **GSE207704** (MCF7/T47D CLDN4−/−) and **GSE22493** (SKOV-3 siRNA), plus older non-KD claudin series |
| GEO | same perturbation terms, PDAT 2024–2026 | **0** series |
| GEO | CLDN4 or Cldn4, PDAT 2024–2026, any design | 11 series; none is a CLDN4/Cldn4 genetic loss expression contrast |
| GEO | exact `CLDN4-/-` / `CLDN4 knockout` / `shCLDN4` / `siCLDN4` / `sgCLDN4` | **GSE207704** only |
| SRA | CLDN4 plus knockout/knockdown/CRISPR/shRNA | 3 experiments, all **4C-seq at the CLDN4 locus in MSH2-KO cells** (not CLDN4 loss, not expression) |
| SRA | same terms, PDAT 2024–2026 | **0** |
| BioStudies | CLDN4 plus knockdown/knockout/CRISPR/shRNA/siRNA (84 records) | Expression accessions: **E-GEOD-22493, E-GEOD-50927** only. E-GEOD-60885 and E-GEOD-84742 are not CLDN4-loss transcriptomes |
| GSA / GSA-Human | `CLDN4` and `Cldn4` | GSA: 11 INSDC mirrors of **GSE50927**, **GSE207704**, and the 4C-seq locus experiments. CRA: 0. OMIX: 0. HRA: 0 for CLDN4/Cldn4 |
| NGDC BioProject | `CLDN4` | 16 projects. The only CLDN4-loss expression projects are **PRJNA128653** (GSE22493) and **PRJNA219385** (GSE50927). **PRJNA856719** is GSE207704 |

GSE274940 (EpH4 WT vs pan-claudin-null, public 2025-09-24) is inside the stated catalog. It was not treated as a new hit and was not re-scored. Prior catalog note stands: it is not a Cldn4-specific loss (Cldn4 RNA is not down).

## Literature checked because the abstract says CLDN4 was perturbed

These papers are 2024–2026 and describe CLDN4/Cldn4 loss or a CLDN4-targeted CRISPR experiment. None added a new expression accession in GEO, ArrayExpress, SRA, or GSA. PubMed records for all eight have an empty DataBank list (checked 2026-09-21).

| PMID | What was done | Expression matrix in the four repositories |
|---|---|---|
| 41016339 | CRISPR CLDN4 knockout in H1688 SCLC; authors report RNA-seq and SAA1 up | **No.** No GEO series under the perturbation, under H1688+CLDN4, or under SAA1+CLDN4 in 2024–2026. Not scored |
| 40892111 | shRNA CLDN4 in cerulein AP (266-6 cells and mice). RNA-seq in the abstract is the AP-vs-control discovery contrast used to pick CLDN4 | **No.** No ferroptosis+CLDN4 GEO series. KD readouts in the abstract are qPCR/Western/ELISA, not a KD transcriptome |
| 41214101 | Stable CLDN4 overexpression and knockdown in HGSC; ISRE reporter | **No.** OA full text: data available on request. No GSE/SRA/HRA accession. RNA-seq mention is a cited T-cell paper |
| 38867360 | CRISPRi and a claudin-mimic peptide in ovarian cancer | **No.** OA full text: data in the manuscript/supplement. No RNA-seq sentence and no accession |
| 41697223 | Germline Cldn4 deletion, abdominal sepsis (permeability, cytokines, survival) | **No.** No GEO series for claudin-4 plus sepsis |
| 42159105 | shRNA CLDN4 in ox-LDL HUVECs and an AP rat model (viability, ROS, histology) | **No** databank accession |
| 42020735 | Genome-wide CRISPR screen that identifies claudin-4 as the BFT receptor | **SRA PRJNA1424384** is the screen (Zenodo 18557106). It is not a CLDN4-knockout expression profile |
| 40592346 | CLDN4 palmitoylation in HCC. scRNA-seq of tumors and RNA-seq of CPP-S4 peptide-treated MHCC97H | **Not in GEO/ArrayExpress/SRA/GSA.** CNGB CNSA **CNP0006650**. Peptide treatment is not CLDN4 KD/KO/CRISPR/shRNA/siRNA |

Author spot checks: Sato-Yazawa plus CLDN4/SAA1/RNA-seq returns no GEO series. Kashiwagi GEO series from 2024–2026 are other genes (TDP-43, AMH, ITGB3), not H1688 CLDN4 knockout.

## Other records opened and dropped

| Accession | Why it is not a new hit |
|---|---|
| GSE22493, GSE50927, GSE207704, GSE274940 | Catalog |
| E-GEOD-22493, E-GEOD-50927, PRJNA128653, PRJNA219385, PRJNA856719, SRX352051/053/054, SRX16067609/610 | Same catalog studies, ArrayExpress or GSA/SRA mirrors |
| SRX8386834 and related | 4C-seq at the CLDN4 locus; the knockout is MSH2 |
| HRA003416 | Name contains “claudin-low” subtype of ccRCC (2022). Not a CLDN4 knockdown |
| PRJCA023797 | CLDN4-positive cancer cells in malignant pleural effusion. Observational, not a perturbation |
| GSE310539, GSE309751, GSE330007, GSE309894, GSE279689, GSE282371, GSE244820, GSE268909, GSE247271, GSE247130, GSE264098 | GEO 2024–2026 text hits for the string CLDN4/Cldn4. Designs are injury, leukoplakia, ALK, fish gonadotropins, atopic dermatitis, coronavirus, hyperglycemia, or CEBPA. Not CLDN4 loss |
| GSE303714, GSE297058 | 266-6 / pancreatitis RNA-seq. Not a Cldn4 knockdown |
| GSE174462 | H1688 shFOXM1, 2022. Not CLDN4 |

## IFN/APM

Prespecified panels from the existing catalog audit (claim C4), not re-fit here:

- CLAIM6: IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A
- IFN/ISG panel (39 genes) and MHC-I/APM panel (18 genes)

They were not applied. A set-level overlap on zero new contrasts would be a fabricated number.

## Reproduce

Queries, counts, and the exclusion table are in `queries.tsv` and `inventory.tsv` next to this file. The search script reprints the same queries:

```bash
python3 scripts/cldn4_kd_match_wave/search_public.py
```
