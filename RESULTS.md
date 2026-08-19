# PhenoCycler / CODEX hunt: lung NSCLC CLDN4 protein + CD8 + XY

**Verdict: no public lung PhenoCycler / CODEX / CODEX-like panel includes CLDN4 protein.** CD8 + XY exist in several public lung datasets. CLDN4 protein does not. No claudin proxy was used. Nearest-CD8 µm, radius counts, and mixing were **not computed** because there is no CLDN4-protein channel to score.

Scope: additive multiplex protein imaging only (PhenoCycler / CODEX / CODEX-like). IMC / MIBI were not re-hunted. Private 8-KL material was not used.

---

## 1. Required lung hit (CLDN4 protein AND CD8 protein AND XY)

**None found.**

A lung dataset is a hit only if the **same** multiplex protein assay measures **CLDN4 protein** (not RNA, not IHC-on-a-serial-section, not another claudin) **and** **CD8 protein**, with cell centroids in XY. No public lung panel met that bar.

---

## 2. Panels checked (lung)

Each row is a panel whose antibody list was read from the paper supplement, vendor/OMAP table, or deposited channel file. CLDN4 = **absent** unless stated.

### 2.1 Takano / Suzuki Visium + PhenoCycler LUAD (PMC11621540)

**Paper:** Takano et al., *Nat Commun* 2024;15:10637. PMC11621540. DOI 10.1038/s41467-024-54671-7.

**What was done:** Visium (30 LUAD cases: 8 IA + 22 AIS/MIA) plus PhenoCycler on representative cases plus Xenium. This is the Visium+PhenoCycler LUAD development paper named in the request.

**Supplement / data (retrieved, not assumed):**

| Resource | Accession / URL | What it is |
| --- | --- | --- |
| Supplement PDF | [41467_2024_54671_MOESM1_ESM.pdf](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-54671-7/MediaObjects/41467_2024_54671_MOESM1_ESM.pdf) (102 MB; Tables S4–S5) | Full PhenoCycler antibody + cycle tables |
| Supplement Data S1 | [41467_2024_54671_MOESM4_ESM.xlsx](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-54671-7/MediaObjects/41467_2024_54671_MOESM4_ESM.xlsx) | 302-gene **Xenium RNA** panel (not protein) |
| Processed public | [DBKERO project page](https://kero.hgc.jp/Ad-SpatialAnalysis_2024.html) | Per-case Visium / PhenoCycler / Xenium downloads |
| Code | [Ad-SpatialAnalysis_2024](https://github.com/asuzuki-asuzuki/Ad-SpatialAnalysis_2024), [LUAD-Spatial](https://github.com/akinaka-dd/LUAD-Spatial) | Seurat + Visium–PhenoCycler alignment |
| Raw images / seq | JGAS000613, JGAS000677 (JGA / NBDC; controlled) | Personally identifiable; not open |

There is **no Zenodo** deposit for this paper. Raw PhenoCycler QPTIFF is controlled-access JGA. Processed PhenoCycler measurement tables **are** open on DBKERO.

**Table S4 — FFPE IA PhenoCycler (LUAD No. 2–4), 35 antibodies.**  
Inventoried: CD31, CD20, beta-actin, Pan-Cytokeratin, CD44, E-cadherin, CD45RO, β-catenin1, Podoplanin, **CD8 (C8/144B, BX026)**, Mac2/Galectin-3, Ki67, CD4, LIF, CD68, CD11c, IDO1, CD3e.  
Custom: α-SMA/ACTA2, Collagen-IV, HLA-A, HLA-DPB1, TP53, CD163, ICOS, TIGIT, DC-LAMP, PD-L1, FOXP3, CD56, LAG3, CXCL13, CTLA4, TTF-1/NKX2-1, PD-1.  
**CLDN4: not listed.**

**Table S5 — FF AIS/MIA PhenoCycler (TSU-23/24/28/30/31/35).**  
CD19, **CD8 (SK1, BX004)**, PD-1, CD68, CD45RO, Pan-Cytokeratin, CD4, Podoplanin, CD11c, CD31, Mac2/Galectin-3, Ki67, TP53, CD163, CXCL13, CTLA4, PD-L1, TTF-1/NKX2-1, FOXP3, α-SMA/ACTA2, CD56.  
**CLDN4: not listed.**

**Public cell table check (DBKERO, not a proxy claim):**  
Downloaded `TSU-23_FF.tar.gz` from DBKERO. File `TSU-23_FF_measurements_mod.csv.gz` has **56,357 cells** and columns:

`Centroid X um`, `Centroid Y um`, `CD8: Cell: Mean`, plus DAPI, PD-1, CD19, CD45RO, CD4, Pan-Cytokeratin, Podoplanin, CD11c, TP53, CD31, NKX2-1, CD163, Ki67, a-SMA, FOXP3.

**No CLDN4 / Claudin column.** CD8 + XY are real; CLDN4 protein is not.

**Xenium (RNA, listed only to close the loop):** Data S1 has **CLDN25**, CD8A. No CLDN4. Visium TME-score gene lists use CLDN18 (well-differentiated) and CLDN5 (endothelium). Those are RNA signatures, not PhenoCycler protein, and they are **not CLDN4**.

### 2.2 Haga et al. 2023 AIS/MIA PhenoCycler (prior LUAD paper)

**Paper:** Haga et al., *Nat Commun* 2023;14:8375. DOI 10.1038/s41467-023-43732-y. Cited by PMC11621540 as the earlier AIS/MIA study.

**Supplement Table S7** (41467_2023_43732_MOESM1_ESM.pdf, retrieved): CD19, **CD8 (SK1)**, LIF, PD-1, CD68, CD45RO, Pan-Cytokeratin, β-catenin1, CD4, Podoplanin, CD11c, CD31, Mac2/Galectin-3, Ki67, Caveolin, CD163, CXCL13, CTLA4, PD-L1, TTF-1, FOXP3, α-SMA.  
**CLDN4: not listed.**

### 2.3 Monkman / Kulasinghe CODEX NSCLC (public cells + XY)

**Paper:** Monkman et al., *J Transl Med* 2024. PMC10910756. DOI 10.1186/s12967-024-05035-8.

**Data:** Zenodo [10.5281/zenodo.10258578](https://doi.org/10.5281/zenodo.10258578) — `4301_cells.csv`, `anndata_4301_v4_annotated.h5ad`, `4301_channelnames.csv`, plus raw OME-TIFFs. Open.

**Panel (Supplementary Table S1 + deposited `4301_channelnames.csv`, both read):** DAPI, CD31, Siglec8, CD4, CD15, CD44, CD107a, CD20, CD38, CD68, CD34, CD45RO, CD45, CD141, CD11b, CD11c, PanCK, Podoplanin, CD197, RORgammaT, **CD8 (C8/144B)**, GranzymeB, CD25, CD21, SPP1, FoxP3, CD56, HLA-DR, Vimentin, Ki67, CD117, CD14, CD163, CD183, CD45RA, CD3e, PGP9.5.  
**CLDN4: not listed.** CD8 + XY are deposited.

### 2.4 HuBMAP lung CODEX OMAPs (healthy / reference lung, not NSCLC; still checked)

| Panel | Source | CD8 | CLDN4 | Notes |
| --- | --- | --- | --- | --- |
| OMAP-16 lung CODEX v1.2 | [lod.humanatlas.io/omap/16-lung-codex](https://lod.humanatlas.io/omap/16-lung-codex/); graph mapped via catalog/RRID to OMAP-25 symbols | Yes (CD8A, C8/144B, 4250012) | **No** | 34 Abs: CDH1, COL1A1, MKI67, ACTA2, SFTPC, HLA-DRA, CD68, PTPRC, CD1C, CD3E, MS4A1, PECAM1, PF4, LYVE1, PanCK, CD8A, AGER, TP63, TPSAB1, CD14, ATP1A1, CD163, COL4A1, KRT5, FOXP3, ITGAX, MPO, SCEL, TUBB3, SCGB1A1, SCGB3A2, MUC5AC, plus AF2727 (RRID AB_2170716; not CLDN4) |
| OMAP-25 lung CODEX v1.1 | [cdn CSV](https://cdn.humanatlas.io/digital-objects/omap/25-lung-codex/latest/assets/omap-25-lung-codex.csv); DOI 10.48539/HBM497.LJTV.953 | Yes (CD8A) | **No** | 45-plex FFPE lung; SFTPC, AGER, SCGB1A1, MUC5AC, PanCK, immune set |
| OMAP-38 lung CODEX v1.0 | [cdn CSV](https://cdn.humanatlas.io/digital-objects/omap/38-lung-codex/latest/assets/omap-38-lung-codex.csv); DOI 10.48539/HBM652.FCHM.439 | Yes (CD8A) | **No** | Senescence-expanded lung panel (GDF15, GLB1, H2AX, LMNB1, CDKN2A, CDKN1A, …) |

protocols.io **813.1** PhenoCycler-Fusion FFPE lung example cocktail (Purkerson / Pryhuber) matches the OMAP-25 style list (SMA, HLADR, FOXP3, ColIV, **CD8**, CD68, CD45, CD4, CD3e, CD31, CD20, CD163, CD14, CD11c, E-cadherin, TPSAB1, SFTPC, SCGB1A1, …). **No CLDN4.**

These OMAPs are **not lung cancer** cohorts. They are listed so the lung CODEX antibody space is covered.

### 2.5 Akoya / Quanterix inventoried and kit panels (lung-cancer–tested reagents)

| Source | CLDN4 | CD8 | Note |
| --- | --- | --- | --- |
| [PhenoCycler Antibody Database](https://www.quanterix.com/phenocycler-antibody-database/) (page retrieved) | **No inventoried CLDN4** | Inventoried (human FFPE 4250012; human FF 4450004 / 4150004) | Only claudin entry: **Claudin-5**, status **Screened**, “Contact for more info” — not CLDN4, not inventoried |
| PhenoCode Discovery IO60 (product biomarker list) | **No** | Yes | Immune typing + architecture (E-cadherin, PanCK, SMA, …). No tight-junction claudin |
| PhenoCode Discovery Tissue Architecture module | **No** | No (architecture module) | E-cadherin, SMA, Vimentin, ColIV, CD31, β-catenin, … |
| PhenoCode Signature 5-plex kits (Immuno-Contexture / Immune Profile / Activated TIL / T Cell Status) | **No** | Yes | PanCK / immune; add-in list is CD163, CD20, CD4, CD45RO, CD68, FoxP3, PD-L1, PD-1, SMA — no CLDN4 |

CLDN4 is not an Akoya inventoried PhenoCycler barcode. A custom conjugation would be required; no public lung study that did so was found.

### 2.6 Other lung CODEX / PhenoCycler reports (panel stated or not inspectable)

| Study | Platform | Public protein table? | CLDN4 | CD8 |
| --- | --- | --- | --- | --- |
| AACR 2023 abs. 4627 Monkman / Kulasinghe chemo + ICI NSCLC | PhenoCycler CODEX, 38 markers | No (abstract; later JTM 36-plex is the open one above) | Not claimed | Yes (T-cell panel) |
| SITC 2023 222-E Enable / Kulasinghe NSCLC TMA n=42 | PhenoCycler-Fusion | No public cell table found | Not claimed | Yes |
| SITC 2022 120 PhenoCycler-Fusion 57-Ab + PhenoCode Signature on NSCLC | Poster / abstract | No | Not claimed | Yes |
| Visiopharm / Akoya IO60 + Metabolic spike-in NSCLC poster | PhenoCycler-Fusion | No public deposit | Not claimed (GLUT1, LDHA, HK1, …) | Yes |
| Zenodo [20263843](https://zenodo.org/records/20263843) / concept [19477006](https://doi.org/10.5281/zenodo.19477006) “CODEX and H&E … lung cancer TMAs” | CODEX | **Files restricted** (login / request). Description does not name antibodies | **Unknown — panel not inspectable** | Unknown |

Restricted Zenodo lung CODEX is **not** treated as a CLDN4 hit.

---

## 3. Optional secondary: epithelial solid-tumor CODEX with CLDN4 + CD8 (labeled not-lung)

**No public hit.** Epithelial CODEX OMAPs and the PDAC multi-omic paper were checked so a non-lung fallback is not missed.

| Dataset / panel | Tissue | CD8 | CLDN4 protein in CODEX | Label |
| --- | --- | --- | --- | --- |
| OMAP-2 intestine CODEX | Intestine | Yes | **No** | not-lung |
| OMAP-15 intestine CODEX | Intestine | Yes | **No** | not-lung |
| OMAP-6 pancreas CODEX | Pancreas | No | **No** | not-lung |
| OMAP-13 pancreas CODEX | Pancreas | Yes | **No** | not-lung |
| OMAP-34 esophagus CODEX | Esophagus | Yes | **No** | not-lung |
| OMAP-29 mouth mucosa CODEX | Oral mucosa | Yes | **No** | not-lung |
| Trevino / Willmann et al. *Nat Commun* 2024 PMC11686138 (PDAC Visium + CODEX + IHC) | PDAC | CD45 / immune in 51-plex CODEX (“as in ref. 50”) | **CLDN4 protein is IHC + Visium RNA, not a CODEX channel.** Fig. 2 is “High-resolution CLDN4 immunohistochemistry.” Methods: “immune and fibroblast populations are profiled using … CODEX and IHC.” | not-lung; **not a dual-protein CODEX hit** |
| Knott TNBC CODEX Zenodo [10045066](https://zenodo.org/records/10045066) | Breast (TNBC) | Yes (celltype CD8T; XY present) | **No** CLDN4 column in deposited phenotype table | not-lung |

Do not treat PDAC IHC CLDN4 or Visium *CLDN4* mRNA as a CODEX protein channel.

---

## 4. Metrics (nearest CD8 µm, radius counts, mixing)

**Not run.** There is no lung (or secondary not-lung) public CODEX/PhenoCycler table with a CLDN4 protein column plus CD8 plus XY. Computing distances on PanCK, E-cadherin, or another claudin would be a proxy and was forbidden.

What *would* be computable on existing public lung CODEX/PhenoCycler tables (CD8 + XY only, **no CLDN4**):

- PMC11621540 DBKERO CSVs: `Centroid X um`, `Centroid Y um`, `CD8: Cell: Mean` (e.g. TSU-23, 56,357 cells).
- Monkman Zenodo 10258578: `4301_cells.csv` + annotated h5ad with CD8 and centroids.

Those are CD8-only spatial tables. They do not answer CLDN4–CD8 mixing.

---

## 5. What was *not* used as a substitute

- E-cadherin, Pan-Cytokeratin, β-catenin, TTF-1, or EpCAM as a “CLDN4-like” epithelial score.
- CLDN5 (Akoya screened; endothelium), CLDN18 / CLDN6 / CLDN25 (RNA in Visium/Xenium of PMC11621540).
- Single-plex CLDN4 IHC (diagnostic LUAD vs mesothelioma literature; not multiplex + CD8 + XY).
- IMC / MIBI (other hunt).
- Private 8-KL.

---

## 6. Practical implication

A public additive multiplex protein dataset that can measure **CLDN4-high tumor cells vs nearest CD8 (µm)** in lung NSCLC **does not currently exist** in the PhenoCycler/CODEX corpus checked above. Closing that gap requires a **custom-conjugated CLDN4** barcode on PhenoCycler (not in the Akoya inventoried catalog as of this hunt) plus CD8, on public lung tissue with deposited cell tables.

Closest **open** lung protein+XY resources for a *future* CLDN4 add-on, not for current CLDN4 metrics:

1. **DBKERO PhenoCycler CSVs** from PMC11621540 (LUAD IA FFPE + AIS/MIA FF; CD8 + µm centroids; QPTIFF under JGA).
2. **Monkman Zenodo 10258578** (NSCLC ICI TMA; 36-plex including CD8; open images + anndata).
