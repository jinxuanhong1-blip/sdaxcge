# Dataset selection notes — `fable_irae` slice

Goal: measure TACSTD2 (TROP2) and CLDN4 (claudin-4) in **public** immune checkpoint
inhibitor (ICI) immune-related adverse event (irAE) transcriptomic datasets, focused on
lung-relevant toxicities (pneumonitis, colitis, etc.), keeping the processed footprint
small (< 2 GB) and everything reproducible.

## How candidates were found

`scripts/fable_irae/01_discover_geo.py` queries NCBI GEO DataSets (`gds`) via E-utilities
with four term sets (irAE; ICI+pneumonitis; ICI+colitis; lung-ICI+toxicity), unions the
hits, and pulls `esummary`. 527 unique records → `results/fable_irae/tables/geo_candidates.tsv`
(102 GSE series; 50 ICI+irAE series after keyword scoring →
`geo_candidates_filtered.tsv`). Supplementary-file sizes were then listed from the GEO
FTP mirror to keep downloads tractable.

## Chosen datasets (3 complementary compartments)

| Accession | Compartment | irAE | Platform | n | Processed file used | Size |
|-----------|-------------|------|----------|---|---------------------|------|
| **GSE206300** | Colon mucosa, **epithelial** single-nucleus RNA-seq | ir**Colitis** | 10x 5'/3' (GPL24676) | 26 donors (12 Case / 14 Control), 81,707 nuclei | `ircolitis-tissue-epithelial.h5ad` (raw UMI counts) | 517 MB gz |
| **GSE277136** | Bronchoalveolar lavage fluid scRNA-seq (**lung**) | ICI **pneumonitis** | 10x (GPL24676) | 4 AE + 3 HC, 74,607 cells | `BLFplusNM.h5ad` (`.raw` = log1p, 22,126 genes) | 1.4 GB |
| **GSE319496** | Whole blood **bulk** RNA-seq | irAE Yes/No (mRCC, nivo+ipi) | Illumina (GPL18573) | 51 (29 Yes / 22 No) | `GEO_raw_counts_SampleID.csv.gz` | 0.7 MB |

Rationale: TACSTD2/CLDN4 are epithelial genes, so measurability depends on the compartment.
The three datasets span an epithelial→immune→systemic gradient, which is exactly what the
"if measurable" clause needs to resolve: colon epithelium (target tissue), BAL/lung
(mostly immune with rare epithelium), and peripheral blood (systemic).

## Datasets considered but not used

- GSE206299 (ircolitis tissue immune compartments), GSE206298 (blood) — immune-sorted, so
  epithelial markers not informative; large.
- GSE144469 (Luoma 2020 ICI-colitis) — only a 1 GB RAW tar of raw 10x, no compact
  processed matrix.
- GSE189185 / GSE190564 (colitis/UC spatial) — 8.7 GB / 11 GB, exceed the budget.
- GSE216329 (ICI-pneumonitis T cells) — immune-only, no epithelium.

## Reproducibility

Run `scripts/fable_irae/run.sh`. Large matrices are **not** committed (see
`results/fable_irae/.gitignore`); they are re-downloaded from GEO FTP by
`02_download.py`. Derived tables/figures under `results/fable_irae/` are the committed
processed outputs (< 1 MB total).
