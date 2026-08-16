# Repository search summary (Zenodo / figshare / OSF)

**Date:** 2026-08-16
**Scope:** open, processed expression/marker matrices for lung ICI (immune
checkpoint inhibitor) and/or the epithelial markers TACSTD2 (TROP2) and CLDN4.

## What was queried

8 query strings x 3 repositories (Zenodo REST API, figshare v2 search, OSF v2
nodes). Scripts: `scripts/fable_zenodo/search_repositories.py`.

- Total unique records returned: **217** (119 Zenodo, 98 figshare; OSF returned
  no title-matching nodes).
- Access breakdown: 205 open, 11 restricted, 1 embargoed.
- After title-based relevance scoring (`filter_candidates.py`): 119 records with
  relevance >= 3; 40 of those expose a downloadable matrix-like file.

## Selected & downloaded (open, CC-BY, <2GB, DOI-verified)

| DOI | Title | Files used | License |
|-----|-------|-----------|---------|
| 10.5281/zenodo.10731914 | Single-cell RNA-seq profiles of lung adenocarcinoma patients and tumor-bearing mice | `cell_matrix_sce_8LUADs.csv.gz` (58 MB), `cell_annotation_sce_8LUADs.csv` (9.6 MB) | CC-BY-4.0 |
| 10.5281/zenodo.8041882 | Acquired resistance to anti-PD1 therapy in patients with NSCLC associates with immunosuppressive T cell phenotype | `sce-clustering-allsamples-original_final.RDS` (160 MB), `AK_List.csv` | CC-BY-4.0 |

## On-topic candidates verified but not downloaded

- 10.5281/zenodo.11205626 — paired normal/LUAD scRNA (h5, 761 MB); redundant with
  10731914 for the two genes.
- 10.5281/zenodo.13947395 — paired LUAD + healthy IMC/scRNA (processed zip 1.5 GB);
  no ICI arm.

## Restricted Zenodo records skipped (per instruction)

- 10.5281/zenodo.18488822 — Vascular STING / NK cell, SCLC (restricted)
- 10.5281/zenodo.18729311 — Multiomic TME after bronchoscopic ... (restricted)
- 10.5281/zenodo.10911472 — single-cell/spatial inflammatory architecture (restricted)
- 10.5281/zenodo.7990870 — Spatial predictors of immunotherapy response, TNBC (restricted)

All DOIs above resolve (HTTP 200); see `results/fable_zenodo/doi_verification.json`.
