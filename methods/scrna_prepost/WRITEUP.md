# Extra table: public lung-tumor scRNA, TACSTD2 / CLDN4 before vs after ICI

Zhejiang paired IHC (user A7) is taken as given. This extra is the public single-cell catalog and the one testable GEO cohort.

## Search

Public human lung-tumor scRNA with **both** a pre-ICI (or chemo-IO) timepoint and a post timepoint, unmatched patients allowed. GEO searches plus must-try accessions (GSE207422, GSE337519, GSE205335, GSE241934, GSE146100) and 2023–2026 leftovers with Pre/Post in the title. Full inventory: `CATALOG.md`.

**One accession deposits both timepoints with tumor cells:** GSE207422 (Hu et al., *Genome Med* 2023) — 3 unmatched pre-treatment biopsies and 12 unmatched post-treatment resections. GEO does not deposit CopyKAT or barcode-level author labels.

Opened and not pooled:

| Accession | Public deposit | Pre vs post tumor cells |
|---|---|---|
| GSE337519 | 1 unlabeled 10x library (8,612 cells; GSM9856929 “one sample”) | Series text claims 1 pair; no timepoint labels |
| GSE205335 | 33 samples / 26 patients; RECIST, no timepoint field | Multi-sample patients are multi-site, not serial |
| GSE241934 | 11 IIT + 34 real-world resected tumors after neoadjuvant IO+chemo | Post-only |
| GSE146100 | 3 nodules from 1 patient after pembrolizumab | Post-only |
| GSE291670 | 6 post resections (anlotinib + PD-1); pre is HRA001033 | Post-only in GEO |
| GSE179994 and other T-sorted series | Pre and post T cells | No malignant compartment |
| HRA006493 | Paired pre/post chemo-IO or anti-VEGFA | Controlled access |

## GSE207422 (sample unit, unmatched Mann–Whitney U)

Epithelial cells are the GEO-available tumor-cell compartment. A stricter malignant-like call (epithelial minus alveolar / club / ciliated programs) leaves P08 with 6 cells and several post samples with 0–3 cells, so that contrast is n=2 vs 7.

| Compartment | Gene | n pre | n post | Median pre | Median post | Direction | U | p |
|---|---|---|---|---|---|---|---|---|
| Epithelial (all libraries) | TACSTD2 | 3 | 12 | 1.77 | 1.55 | down | 27 | 0.23 |
| Epithelial (all libraries) | CLDN4 | 3 | 12 | 1.59 | 1.53 | down | 20 | 0.84 |
| Epithelial (NMPR post only) | TACSTD2 | 3 | 8 | 1.77 | 1.70 | down | 15 | 0.63 |
| Epithelial (NMPR post only) | CLDN4 | 3 | 8 | 1.59 | 1.58 | down | 13 | 0.92 |
| Malignant-like (≥10 cells) | TACSTD2 | 2 | 7 | 1.81 | 1.81 | down | 7 | 1.00 |
| Malignant-like (≥10 cells) | CLDN4 | 2 | 7 | 1.47 | 1.34 | down | 8 | 0.89 |

Metric: mean log1p(CP10K). Patients are unmatched. P01 is pre-treatment with pathologic response NE.

## Combined estimate

One independent tumor-cell cohort. Direction count is 0 up / 1 down for both genes. Stouffer equals that study’s signed two-sided p: TACSTD2 Z = −1.19, p = 0.23; CLDN4 Z = −0.20, p = 0.84.

Machine tables: `tables/paper_table_scrna_prepost.tsv`, `tables/combined_stouffer.tsv`, `tables/gse207422_per_sample.tsv`.
