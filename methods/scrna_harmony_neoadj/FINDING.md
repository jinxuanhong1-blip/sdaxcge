# Harmony neoadjuvant lung scRNA — FINDING

**Public data only. Additive to User A3. The GSE207422 CopyKAT slide was not re-run.**

**Change of plan (see `FINDING_combos.md`):** do not treat the four-series / three-series object below as the claim. Pairwise and leave-one-out Harmonys were run; extra figures only for combos that hold the User A3 direction. The three-series object (this file) did **not** hold that direction.

## 一句话结论 / TL;DR

Joint Harmony object: **429,593 cells** from **66 patients** (GSE207422 + GSE241934 + GSE291670; embedding used 77,572 sample-stratified cells). GSE205335 has epithelium but was left out (palliative RECIST, not neoadjuvant MPR). After a ≥20 malignant-like-cell floor, malignant-like *TACSTD2* was **not** higher in NMPR than MPR (median 1.48 vs 1.90; n=24 vs 9; p=0.92). *CLDN4* likewise NS (1.34 vs 1.08; p=0.86). Cohort-adjusted OLS: *TACSTD2* β=−0.27 (p=0.28). Unadjusted *TACSTD2* vs joint-embedding T/NK was **positive** (ρ=0.49, n=32, p=0.0044) and is **cohort-driven** (GSE291670 snRNA sits at low score and low T/NK); partial Spearman after residualizing cohort ρ=0.20, p=0.26. LUAD vs LUSC NS (n=21 vs 4; p=0.50). This does not retract or replace User A3.

## n

| | Count |
|---|---|
| Cells in deposited matrices (included series) | 429,593 |
| Cells in Harmony embedding (cap 1,200 / sample) | 77,572 |
| Patients / samples | 66 / 66 |
| Shared gene symbols (3-way intersection) | 16,209 (Harmony used 781) |
| Post-tx samples with MPR label | 63 (MPR 24, NMPR 39) |
| Eligible for malignant tests (≥20 malignant-like) | 33 (MPR 9, NMPR 24) |
| Dropped at the ≥20 floor | 30 (GSE241934 26, GSE207422 4; 2 with zero malignant-like) |
| GSE207422 | 15 samples / 92,330 cells (12 post-tx; 3 pre-tx excluded from MPR tests) |
| GSE241934 | 45 samples / 308,196 cells (IIT 78,691 + Real 229,505) |
| GSE291670 | 6 samples / 29,067 cells (snRNA) |

Harmony epithelial counts track author Epi on GSE241934 (IIT 1,734 vs author 1,699; Real 14,026 vs 12,785). The ≥20 **malignant-like** floor is stricter: epithelial clusters minus a high normal-lung score. Residual tumors after neoadjuvant therapy often fail that filter.

## Left out

- **GSE205335** — author epithelium and malignant cells are present (Epithelial cells=31,519; Malignant cells=28,512; 96,505 identity barcodes). The series is palliative / advanced ICI with **RECIST, not neoadjuvant MPR**. Mixing RECIST into an MPR model is an endpoint swap. Public UMI is an R `dgCMatrix` RDS (not MTX). CellIdentity only was downloaded.

## Primary (pre-specified; patient unit; ≥20 malignant-like)

| Test | n | Result | p |
|---|---|---|---|
| TACSTD2 NMPR vs MPR (median log1p CP10k) | 24 vs 9 | 1.48 vs 1.90 (Δ=−0.42; MPR higher) | 0.92 |
| CLDN4 NMPR vs MPR | 24 vs 9 | 1.34 vs 1.08 (Δ=+0.26) | 0.86 |
| TACSTD2 OLS NMPR vs MPR + cohort | 33 | β=−0.265 | 0.28 |
| CLDN4 OLS NMPR vs MPR + cohort | 33 | β=−0.207 | 0.28 |
| TACSTD2 vs T/NK Spearman | 32 | ρ=+0.49 | 0.0044 |
| CLDN4 vs T/NK Spearman | 32 | ρ=+0.59 | 0.00042 |
| TACSTD2 vs T/NK partial Spearman (cohort) | 32 | ρ=+0.20 | 0.26 |
| CLDN4 vs T/NK partial Spearman (cohort) | 32 | ρ=+0.39 | 0.026 |
| TACSTD2 LUSC vs LUAD | 4 vs 21 | median 1.57 vs 1.87 | 0.50 |
| CLDN4 LUSC vs LUAD | 4 vs 21 | median 1.34 vs 1.47 | 0.50 |

Per-cohort MPR (same floor; often underpowered):

| Cohort | TACSTD2 NMPR vs MPR | p |
|---|---|---|
| GSE207422 | 7 vs 1 (only one MPR cleared the floor; P06/P11/P14 dropped) | not testable |
| GSE241934 | 14 vs 5; median 1.55 vs 2.32 (MPR higher) | 0.16 |
| GSE291670 | 3 vs 3; median 0.062 vs 0.167 (MPR higher; snRNA scale) | 0.20 |

## Sensitivity (same scores; looser floor)

| Floor | TACSTD2 NMPR vs MPR | OLS β (cohort) | TACSTD2 vs T/NK ρ / partial ρ |
|---|---|---|---|
| ≥5 cells | 31 vs 16; 1.60 vs 1.98; p=0.24 | −0.33; p=0.059 | +0.35 p=0.016 / +0.08 p=0.59 |
| ≥10 cells | 29 vs 11; 1.60 vs 1.90; p=0.90 | −0.18; p=0.39 | +0.39 p=0.014 / +0.07 p=0.67 |
| ≥20 cells (primary) | 24 vs 9; 1.48 vs 1.90; p=0.92 | −0.27; p=0.28 | +0.49 p=0.0044 / +0.20 p=0.26 |

No floor recovers NMPR>MPR at p<0.05. The unadjusted positive T/NK ρ shrinks and loses significance for *TACSTD2* after cohort residualization.

## Read this as extra, not as a copy of the given slide

- Cell-level p-values are not reported (pseudoreplication).
- User A3 / GSE207422 CopyKAT IDs are still not on GEO. Malignant-like here is Harmony epithelium minus a high normal-lung score.
- GSE291670 is **nuclear RNA**; raw *TACSTD2*/*CLDN4* CP10k is ~10× lower than the 10x whole-cell series. Harmony batch-corrects the embedding, not the patient-level scores — hence the cohort covariate.
- GSE241934 has no LUSC (LUAD / ASC only). LUAD vs LUSC is almost entirely GSE207422 Adeno vs Squamous (n_LUSC=4).
- Unadjusted ρ>0 vs T/NK is the **opposite** of the User A3 window (−0.40 to −0.50) and is not interpreted as a biological inversion after the partial Spearman.

## Files

- `figures/fig1_umap_harmony_dataset_lineage.png`
- `figures/fig2_malignant_tacstd2_cldn4_by_mpr.png`
- `figures/fig3_malignant_score_vs_tnk.png`
- `figures/fig4_luad_vs_lusc.png`
- `tables/per_sample.tsv`, `stats.json`, `inventory.json`, `sensitivity_floors.json`, `dropouts_lt20_malignant.tsv`
- `paper_snippet.md`
