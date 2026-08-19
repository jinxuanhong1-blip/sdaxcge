# CLDN4 protein + CD8 protein spatial proteomics hunt

**0 public datasets had CLDN4 protein and CD8 protein on the same multiplex panel with cell XY coordinates. 0 of those were lung.**

No spatial-metric numbers were computed. Running nearest-CD8 distance, radius neighbor counts, or a mixing score on a panel that lacks CLDN4 protein would be fabrication. This note is the hunt result, not a claim-failed analysis.

## What was required

Same-section, cell-level protein table with:

1. CLDN4 / Claudin-4 **protein** intensity (not *CLDN4* RNA, not “claudin-low” subtype, not CLDN18.2)
2. CD8 / CD8a **protein** intensity
3. XY (or centroid) coordinates
4. Public download without EGA/dbGaP access

Lung first. If lung is empty, the best non-lung epithelial tumor set with **both** proteins. No private 8-KL material.

## Headline counts

| Filter | n | Source |
|---|---|---|
| Public series / kits / papers screened | 41 | `hunt/screened_inventory.csv` |
| Had **CD8 protein** + a public cell-XY table | 14 | same file, `cd8_protein=yes` and `public_cell_xy=yes` |
| Had **CLDN4 protein** on that same panel | **0** | same file, `qualified=yes` |
| Of those 0, lung | **0** | |
| Non-lung epithelial tumor fallback with both proteins | **0** | |
| Spatial metrics run | **0** | |

Many multiplex panels have CD8 and a generic epithelial marker (panCK, E-cadherin, EpCAM, KRT8/18). That is expected and is **not** a CLDN4-protein series.

## Lung (empty)

These public lung spatial-protein resources have CD8 (or an immune panel that includes it) and **do not** list CLDN4 protein:

| Resource | Tech | CLDN4 protein? | CD8 protein? | Why not run |
|---|---|---|---|---|
| Sorin et al. 2023 *Nature* LUAD IMC (416 pts; Zenodo 7760826) | IMC 35-plex | No (panCK / lineage, not CLDN4) | Yes | Missing CLDN4 protein |
| Sorin et al. 2025 *Nat Commun* driver-mutation LUAD IMC | same 35-plex | No | Yes | Missing CLDN4 protein |
| Cords et al. early LUAD IMC (Zenodo 14827075) | IMC | Not a public cell table (files restricted) | Unknown public | Restricted; not downloaded |
| Bruker CosMx NSCLC demo / Liu CAF reanalysis | CosMx **RNA** | *CLDN4* RNA only | CD8 RNA/protein mixed | Not CLDN4 **protein** |
| CosMx Human IO 64-plex protein (commercial kit used on lung) | CosMx protein | No (full published target list) | Yes | Missing CLDN4 |
| Xenium protein 27-plex (tumor + immune subpanels) | Xenium protein | No | Yes (CD8A) | Missing CLDN4 |
| GeoMx IO protein core/modules | GeoMx DSP | No | Yes | No cell XY; no CLDN4 |
| GSE271689 NSCLC “spatial multi-omics” | mostly sequencing | Not a CLDN4+CD8 protein cell table | — | Wrong assay class |

GEO eutils queries that combined `CLDN4`/`Claudin-4` with IMC/CODEX/MIBI/PhenoCycler/CosMx/Xenium returned 9 GDS IDs. All 9 were RNA/ATAC/CyTOF-blood false positives (GSE200036, GSE253013, GSE255031, GSE271689, GSE285215/287/288, GSE295969, GSE316782). None are multiplex imaging cell tables with both proteins.

## Non-lung epithelial fallback (also empty)

Best GI/breast/pancreas public protein tables were inspected. **None** carry CLDN4 protein + CD8 on one panel:

| Resource | Tissue | Tech | CD8 | CLDN4 protein | Evidence |
|---|---|---|---|---|---|
| Hickey et al. 2023 *Nature* / Zenodo 7311360 | healthy intestine | CODEX 57-plex | Yes (column `CD8`) | **No** | Header downloaded; markers are MUC2…CD161 + OLFM4/FAP/CD25/CollIV/CK7/MUC6 |
| Schürch et al. 2020 *Cell* CRC | CRC | CODEX 56-plex | Yes | **No** | Immune/structural panel; PMC text has no CLDN/Claudin |
| Jackson et al. 2020 *Nature* breast | breast | IMC | Yes | **No** | Standard breast IMC (ECAD/CK/ER/HER2), not CLDN4 |
| Keren et al. MIBI TNBC | breast | MIBI-TOF | Yes | **No** | Published 36-plex lacks CLDN4 |
| IMMUcan / Bodenmiller example (Zenodo 5949116 `panel.csv`) | mixed (incl. NSCLC, CRC, H&N) | IMC | Yes (`CD8a`) | **No** | Panel file downloaded |
| OMAP-2 intestine (HuBMAP) | intestine | CODEX | Yes | **No** | OMAP text: pCK, MUC, immune; no CLDN4 |
| OMAP-34 esophagus (Zenodo 17613542) | esophagus | CODEX 31-plex | Yes | **No** | KRT5/13/14/16, IVL, PPL, MUC1, CD8; no CLDN4 |
| PhenoCode Discovery IO60 | multi-tumor kit | PhenoCycler | Yes | **No** | Published biomarker list (PanCK/EPCAM/KRTs, no CLDN4) |
| CosMx Human IO 64-plex | multi-tumor kit | CosMx protein | Yes | **No** | Published target list |
| Scheuermann et al. 2024 MACSima 118-plex | HCC / solid TIME | cyclic IF | Yes | **No** | Frontiers Table 1: no CLDN/Claudin string |
| Karlsen et al. 2024 breast CAF IMC 42-plex | breast | IMC | Yes | **No** | Published table: PanCK/ECAD/CK5/CK8/18, no CLDN4 |
| Wang/Ferrara et al. 2024 *Nat Commun* PDAC | pancreas | Visium + IHC + 51-plex CODEX | CODEX has immune markers | CLDN4 is **Visium RNA + single-plex IHC**, not on the 51-plex | Same-panel rule fails |
| Bolen CosMx sigmoid protein (Zenodo 14851272) | colon | CosMx protein | Yes (kit) | **No** (kit) | Standard 64-plex |

The closest **CLDN4 protein** public imaging is the Ferrara/Wang PDAC **single-plex IHC** (and PET peptide work). That is not multiplex with CD8 and has no public CLDN4+CD8 cell table.

## Repositories queried

| Source | Query style | Usable CLDN4+CD8+XY hit |
|---|---|---|
| NCBI GEO / GDS eutils | CLDN4 or Claudin-4 ∩ IMC/CODEX/MIBI/PhenoCycler/CosMx/Xenium | 0 |
| Zenodo API | quoted Claudin-4 ∩ CODEX/IMC/MIBI/PhenoCycler | 0 (title/description quoted search total 0) |
| Europe PMC | Claudin-4 ∩ IMC / PhenoCycler / CosMx protein / mIF+CD8 | papers mention CLDN4 RNA or single-plex IHC, not dual-protein cell tables |
| Figshare API | CLDN4 CODEX/IMC/multiplex | 0 relevant datasets |
| Dryad | Hickey intestine CODEX (same as Zenodo 7311360) | CD8 only |
| HuBMAP / OMAP | intestine, esophagus, pancreas, kidney OMAPs | no CLDN4 |
| HTAN docs / published atlases | multiplex microscopy + PDAC Visium | no public CLDN4+CD8 protein cell table |
| Synapse | syn61831984 (Ferrara Stanford) | Visium **transcriptome**, not protein |
| ImmPort | CLDN4 text search | no multiplex imaging cell table with both proteins |
| CosMx / Xenium / GeoMx / Akoya IO kit lists | official marker tables | CD8 yes, CLDN4 no |
| EGA / dbGaP | skipped (no access), per instructions | — |

## Why this is empty (not a failed run)

CLDN4 is a tight-junction epithelial marker used in single-plex IHC and in spatial **transcriptomics**. Off-the-shelf IO multiplex protein panels were built for immune lineage + pan-epithelial (panCK/EPCAM/ECAD), not claudins. Custom CLDN4 metal/oligo conjugates exist as reagents; they are not on the public lung/GI/breast cell tables we could open.

## What was not done

- No nearest-CD8 distance, radius counts, or mixing score
- No substitution of panCK/EPCAM/ECAD for CLDN4
- No *CLDN4* RNA CosMx/Xenium/Visium as protein
- No private 8-KL
- No EGA/dbGaP download

## Figure

`figures/fig1_hunt_funnel.png` — hunt funnel from `hunt/screened_inventory.csv`: 41 screened → 14 CD8+public-XY → 0 CLDN4+CD8 same panel → 0 lung.

## Reproducibility

- `hunt/scripts/query_public_catalogs.py` — GEO / Zenodo / Europe PMC / Figshare queries
- `hunt/panels_checked.md` — panel-level yes/no
- `hunt/evidence/hickey_codex_header.txt` — first line of Zenodo 7311360 cell table
- `hunt/evidence/immucan_panel.csv` — Zenodo 5949116 IMC panel
