# CLDN4 KD match wave (2024–2026 refresh)

Search date: **2026-09-21**.

**New public CLDN4/Cldn4 KD/KO/CRISPR/shRNA/siRNA expression datasets outside the catalog: 0.**

The second pass (NODE, CNGB, Figshare, Zenodo, author GitHub, and full-text accession mining) also found **no new matrix**.

IFN/APM and NHEJ/STING/IFN were **not scored**. There is no new count matrix, FPKM table, or DEG table to download. No overlap count, logFC, or p-value was imputed.

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

## Expanded hunt (NODE, CNGB, Figshare, Zenodo, GitHub, full text)

Searched 2026-09-21 after the GEO/ArrayExpress/SRA/GSA pass. Still **0** new CLDN4/Cldn4-loss expression matrices.

| Source | Query | Result |
|---|---|---|
| NODE `POST /node/api/app/browse/search` | `queryWord` = CLDN4, Cldn4, claudin-4, "CLDN4 knockout", "H1688 CLDN4" | **0** each. The same endpoint returns 7,888 hits for `lung`, so the empty CLDN4 result is not a dead API |
| CNGB `search/ajax/project` | CLDN4 | 5 hits: catalog PRJNA128653 (GSE22493), three pre-2024 ceRNA projects, and **CNP0006650** |
| CNGB project | H1688, Kashiwagi, Sato-Yazawa | **0** |
| Figshare article search | "H1688 CLDN4", "CLDN4 CRISPR RNA-seq", "CLDN4 knockout" | **0** |
| Zenodo title `CLDN4` | metadata.title:CLDN4 | 1 dataset, not a CLDN4 knockdown (see below) |
| GitHub | repo search `CLDN4`; code search `CLDN4 knockout`, `H1688 CLDN4` | one public repo, variant tables, not expression |
| Europe PMC annotations | PMID 41016339, 40892111, 41214101, 38867360, 41697223, 42159105 and the OA PMCs | no GSE/SRA/CNP/OEP/Zenodo/Figshare expression accession on the H1688 paper |

### Records opened

| ID | What it actually is |
|---|---|
| PMID 41016339 (10.1016/j.bbrc.2025.152710) | Unpaywall: no OA PDF. Elsevier `mmc1`–`mmc5` URLs for PII S0006291X25014263 returned 404. Europe PMC text-mined accessions are only Addgene plasmids (RRID:Addgene_98293, Addgene_52961). No expression deposit to score |
| PMID 40892111 supplement `10142_2025_1683_MOESM1_ESM.docx` | Primer sequences (GAPDH, CLDN4, GPX4, ACSL4). Data-availability line: available from the corresponding author on request. Not a matrix |
| Zenodo 16886088 (10.5281/zenodo.16886088), 2025-09-24 | Source data for Science Advances adx7431 (27-claudin EpH4 paper). The paper states the RNA-seq is **GSE274940** (catalog). Figure 4 and figure S9 source workbooks are conductance / dilution-potential / strand counts for single-Cldn addbacks, not gene-expression matrices |
| Zenodo 21976731, 2026-08-17 | Title names CLDN4 as a gastric-cancer biomarker with KRT18, GPRC5A, EPCAM, and CLDN7. Zip central directories: qPCR/wound images for **si-KRT18**, plus a re-pack of public **GSE163558** and TCGA TPM. Not a CLDN4 knockdown transcriptome |
| CNP0006650 / CSE0000463 | CNGB project already tied to PMID 40592346 (HCC scRNA and CPP-S4 peptide RNA-seq). Peptide treatment is not CLDN4 KD/KO/CRISPR/shRNA/siRNA |
| github.com/poojascis/CLDN4_Data | `Supple_Table1_nsSNV.csv` is missense variants on ENST00000340958.4. No count matrix |

Figshare hits for the string "Sato-Yazawa" are other people with the same surname (amorphous alumina, spider silk), not the Dokkyo CLDN4 paper.

## IFN/APM and NHEJ/STING/IFN

Prespecified panels, not re-fit on a result:

- CLAIM6: IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A
- IFN/ISG (39 genes) and MHC-I/APM (18 genes), as in the catalog audit
- NHEJ: XRCC6, XRCC5, PRKDC, LIG4, XRCC4, NHEJ1, DCLRE1C, PAXX, POLL, POLM
- STING axis: CGAS, STING1, TBK1, IRF3, IFI16, TREX1, CCL5, CXCL10, IFNB1

They were not applied. A set-level score on zero new contrasts would be a fabricated number.

## Reproduce

Queries, counts, and the exclusion table are in `queries.tsv` and `inventory.tsv` next to this file. The search script reprints the same queries:

```bash
python3 scripts/cldn4_kd_match_wave/search_public.py
```
