# GSE221322 GeoMx protein DSP — EpCAM/PanCk vs CD8/CD3/PD-1/CD45

Additive spatial-protein readout for **B6** (CLDN4-high epithelium avoids immune niches; taken as given from the Visium/CosMx work). This folder does not re-analyse those assays.

Public only: Monkman et al., *Immunology* 2023 (PMID 37022147); GEO [GSE221322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221322). Immunotherapy-treated NSCLC TMA. Methods: `methods/gse221322_dsp_protein/methods.md`. Script: `scripts/gse221322_dsp/`.

## TL;DR (honest)

**CLDN4 and TROP2 protein are not on the deposited 68-plex.** The closest barrier / epithelial antibodies are **EpCAM (TACSTD1)** and **PanCk**.

On the author-normalized tumour segments, EpCAM is inversely associated with **CD45** (leukocyte) protein. CD8, CD3 and PD-1 have the same sign and smaller |ρ|; they are not significant at this n. PD-1 sits on isotype background.

| Test | n | ρ | p | q (BH, 12 EpCAM tests) |
|---|---:|---:|---:|---:|
| Tumour AOI EpCAM vs CD45 | 47 | **−0.42** | **0.0030** | **0.037** |
| Tumour AOI EpCAM vs CD8 | 47 | −0.19 | 0.21 | 0.63 |
| Tumour AOI EpCAM vs CD3 | 47 | −0.08 | 0.59 | 0.85 |
| Tumour AOI EpCAM vs PD-1 | 47 | −0.09 | 0.54 | 0.85 |
| Patient-mean tumour EpCAM vs CD45 | 40 | **−0.48** | **0.0019** | — |
| Stroma AOI EpCAM vs CD45 | 48 | −0.32 | 0.024 | 0.29 |
| Tumour AOI PanCk vs CD4 | 47 | −0.32 | 0.031 | 0.37 |
| Tumour AOI PanCk vs CD8 | 47 | −0.02 | 0.89 | 0.93 |

Median-split tumour EpCAM: CD45 lower in the high half (Δ = −0.76, Mann–Whitney p = 0.014, n_hi = 24, n_lo = 23). CD8 and CD3 are the same direction and not significant (p = 0.39 and 0.46).

This is compartment-level protein on a TMA, not a single-cell neighborhood. It is extra protein evidence that epithelial/barrier-high tumour AOIs carry less leukocyte protein, using the closest antibodies the panel actually has.

## Panel

68 antibodies including 3 isotype controls (GPL29263: Immune Cell Profiling, IO Drug Target, Immune Cell Typing, Immune Activation, PI3K/AKT, MAPK, Cell Death). Full list: `tables/panel_inventory.csv`.

| Requested | On panel? |
|---|---|
| CLDN4 / Claudin-4 | no |
| TROP2 / TACSTD2 | no |
| any other claudin | no |
| EpCAM (TACSTD1) | yes — closest TACSTD-family / barrier proxy |
| PanCk | yes — epithelial morphology marker |
| CD8, CD3, PD-1 | yes |
| CD45, CD4, GZMB, CD20, FOXP3, CD68, HLA-DR, PD-L1, CD56 | yes |

EpCAM and PanCk are only weakly related in tumour AOIs (ρ = 0.21, p = 0.17, n = 47). They are not interchangeable.

## Cohort (deposited files)

GEO abstract: 42 TMA cores. Deposited tables: **48 paired tumour/stroma AOIs, 41 patients** (16 responder, 24 non-responder, 1 N/A). 35 patients have one core, 5 have two, 1 has three.

QC drop: tumour AOI `026` (0 nuclei, area 54 µm²). After QC: **47 tumour AOIs / 40 patients**, **48 stroma AOIs / 41 patients**, **47 paired cores**.

PD-1 median S/N vs mean isotype = 1.10; **0 / 96 AOIs** exceed 2× mean IgG. Tumour PD-1 tracks Ms IgG1 (ρ = 0.51, p = 2.9×10⁻⁴). CD8 median S/N = 4.4 (98% of AOIs > 2×). Treat PD-1 as near the assay floor.

## Results

### 1. Tumour AOIs — EpCAM vs immune (primary)

Author-normalized matrix. Two-sided Spearman. BH over the 12 prespecified neighborhood proteins.

| Immune protein | ρ | 95% CI (Fisher z) | p | q |
|---|---:|---:|---:|---:|
| CD45 | −0.42 | −0.63, −0.15 | 0.0030 | 0.037 |
| HLA-DR | −0.23 | −0.49, 0.06 | 0.12 | 0.49 |
| CD4 | −0.23 | −0.48, 0.06 | 0.12 | 0.49 |
| CD8 | −0.19 | −0.45, 0.11 | 0.21 | 0.63 |
| PD-1 | −0.09 | −0.37, 0.20 | 0.54 | 0.85 |
| CD3 | −0.08 | −0.36, 0.21 | 0.59 | 0.85 |
| GZMB | −0.08 | −0.36, 0.22 | 0.61 | 0.85 |
| CD20 | −0.06 | −0.34, 0.23 | 0.70 | 0.85 |
| PD-L1 | −0.05 | −0.33, 0.24 | 0.76 | 0.85 |
| FOXP3 | −0.04 | −0.33, 0.25 | 0.77 | 0.85 |
| CD68 | 0.01 | −0.28, 0.30 | 0.95 | 0.95 |
| CD56 | 0.10 | −0.19, 0.37 | 0.51 | 0.85 |

CD8 vs CD3 in the same tumour AOIs: ρ = 0.47, p = 8.7×10⁻⁴ (the T-cell pair behaves).

Partial Spearman, EpCAM vs CD45, residualized on PanCk: ρ = −0.40, p = 0.0058. On nuclei: ρ = −0.42, p = 0.0031. On Histone H3: ρ = −0.45, p = 0.0014. The leukocyte inverse is not explained by PanCk content or nuclei count.

### 2. Patient means and stroma

Averaging cores per patient (tumour n = 40) leaves EpCAM vs CD45 at ρ = −0.48, p = 0.0019. EpCAM vs CD8 remains ρ = −0.16, p = 0.32.

Stroma AOIs (n = 48): EpCAM vs CD45 ρ = −0.32, p = 0.024 (q = 0.29 among 12 stroma EpCAM tests). EpCAM vs CD8 / CD3 / PD-1: ρ = −0.20 / −0.23 / −0.18, all p > 0.12.

### 3. Same-core neighbourhood (tumour barrier vs stroma immune)

47 paired cores. Tumour EpCAM vs stroma CD8 ρ = −0.11, p = 0.48. vs stroma CD3 ρ = −0.06, p = 0.71. vs stroma PD-1 ρ = −0.28, p = 0.052. vs stroma CD45 ρ = −0.12, p = 0.41. Patient-averaged cross-compartment numbers are in `tables/spearman_cross_patient.csv` (n = 40) and stay in the same range.

DSP here is two masks per TMA core, not a radius around a CLDN4-high cell. The within-tumour-mask CD45 inverse is the cleaner protein readout.

### 4. Normalization sensitivity

Simple log2(QC count / area) and log2(QC / nucleus) **do not** reproduce the EpCAM–CD45 inverse (area-scaled ρ = −0.13, p = 0.39; per-nucleus ρ = −0.04, p = 0.77). On those raw-scaled matrices PanCk vs CD8 becomes positive (area-scaled ρ = 0.39, p = 0.0064). The primary claim is therefore tied to the **author-deposited normalized matrix**, which is the file GEO provides for analysis and the scale used in Monkman et al. Both matrices are public; both are reported.

### 5. ICI response (sanity check only)

Not a B6 test. Stromal patient-mean EpCAM is higher in responders than non-responders (median 8.44 vs 6.64, n = 16 vs 24, Mann–Whitney p = 0.035), in the same direction as the published stromal EpCAM result. Tumour EpCAM, PanCk, CD8, CD3, PD-1 and CD45 do not separate response in this re-use of the deposited table (all p > 0.30).

## What this adds

B6 is the spatial-transcriptomic pattern that CLDN4-high epithelium sits away from immune niches. GSE221322 cannot score CLDN4 or TROP2 protein. What it can score is EpCAM, the other TACSTD-family epithelial/barrier protein on a 68-plex IO panel, in tumour vs stroma masks from an ICI-treated NSCLC TMA. On the deposited normalized protein matrix, EpCAM-high tumour AOIs have less CD45 protein (n = 47, ρ = −0.42, p = 0.0030, q = 0.037). CD8/CD3/PD-1 point the same way with this n. That is extra spatial-protein evidence at compartment resolution, not a second single-cell neighborhood test.

## Files

| File | Contents |
|---|---|
| `stats.txt` | Key n / ρ / p |
| `provenance.json` | GEO URLs and key dump |
| `tables/panel_inventory.csv` | 68-plex + requested absences |
| `tables/qc_signal_vs_igg.csv` | S/N vs isotype |
| `tables/aoi_compact.csv` | QC-pass AOIs, barrier + immune |
| `tables/spearman_barrier_vs_immune.csv` | AOI and patient Spearman |
| `tables/spearman_cross_core.csv` / `spearman_cross_patient.csv` | Tumour barrier vs stroma immune |
| `tables/partial_spearman_tumour.csv` | EpCAM vs immune \| PanCk / nuclei / H3 |
| `tables/mannwhitney_highlow.csv` | Median-split EpCAM / PanCk |
| `tables/spearman_qc_area_nuclei_sensitivity.csv` | QC area- and nucleus-scaled |
| `tables/response_mannwhitney.csv` | ICI response sanity check |
| `figures/forest_tumour_barrier_vs_immune.*` | Tumour ρ forest |
| `figures/forest_stroma_barrier_vs_immune.*` | Stroma ρ forest |
| `figures/scatter_tumour_epcam_cd45_cd8.*` | EpCAM vs CD45 and CD8 |
| `figures/box_tumour_epcam_highlow.*` | Median-split boxes |
| `figures/heatmap_rho_compartments.*` | Tumour + stroma ρ |
| `figures/forest_cross_tumour_epcam_vs_stroma.*` | Same-core cross-compartment |

Raw GEO CSVs are not in git. Rebuild with `bash scripts/gse221322_dsp/00_fetch.sh && python3 scripts/gse221322_dsp/analyze.py`.
