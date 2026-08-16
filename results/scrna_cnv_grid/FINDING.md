# Malignant-definition grid — which public calls recover the user direction

User A3 direction (taken as given): **malignant TACSTD2 higher in NMPR than MPR**, and **Spearman ρ < 0 vs T/NK fraction**. This slice grids public malignant definitions on GSE207422 and GSE241934. Sample is the unit. Primary metric: mean `log1p(CP10k)` TACSTD2 in the called-malignant compartment. pCR counted as MPR.

Author CopyKAT barcodes are **not** on GEO for GSE207422. GSE241934 deposits `major.cell.type==Epi` (IIT 1,699 cells; paper CopyKAT malignant N=1,669).

---

## Definitions that recover both directions

### GSE207422, 12 post-treatment samples (NMPR 8 / MPR 4)

Two public axes recover both directions **without** third-party T/NK labels:

| Malignant def | Cells | T/NK def | n (NMPR/MPR) | mean NMPR | mean MPR | p | ρ | p(ρ) |
|---|---|---|---|---|---|---|---|---|
| `marker_all_epi` | 10,669 | lineage T/NK | 12 (8/4) | 1.598 | 1.256 | 0.109 | −0.098 | 0.762 |
| `marker_all_epi` | 10,669 | UMI CD3E/CD8A/NKG7 | 12 (8/4) | 1.598 | 1.256 | 0.109 | −0.035 | 0.914 |
| `marker_epi_nnl_p75` | 8,001 | lineage T/NK | 10 (7/3) | 1.950 | 0.971 | 0.183 | −0.394 | 0.260 |
| `marker_epi_nnl_zero` | 4,945 | lineage T/NK | 8 (6/2) | 1.964 | 0.342 | 0.071 | −0.238 | 0.570 |
| `intersect_markerNNL_cnvP95` | 5,033 | lineage T/NK | 8 (7/1) | 2.074 | 0.907 | 0.250 | −0.310 | 0.456 |
| `thirdparty_drmref` | 2,051 | lineage T/NK | 12 (8/4) | 1.571 | 1.120 | 0.154 | −0.126 | 0.697 |

`marker_all_epi` and `thirdparty_drmref` keep all 12 post-treatment samples. Normal-lung filters and intersections drop MPR samples that have almost no residual epithelium (P11/P14/P06), so those n/p values are not the same contrast.

Closest to the user’s ρ band with a complete n=12:

| Malignant def | T/NK def | n | NMPR vs MPR | ρ | p(ρ) |
|---|---|---|---|---|---|
| `thirdparty_drmref` | DRMref T/NK | 12 (8/4) | 1.571 vs 1.120, p=0.154 | **−0.490** | 0.106 |
| `cnv_mu2` | DRMref T/NK | 12 (8/4) | 1.831 vs 1.445, p=0.283 | **−0.427** | 0.167 |
| `marker_epi_nnl_p75` | DRMref T/NK | 10 (7/3) | 1.950 vs 0.971, p=0.183 | **−0.758** | **0.011** |

`marker_epi_nnl_p75` × DRMref T/NK is the only primary row with p(ρ)<0.05. It uses a third-party T/NK fraction and drops two MPR samples.

### GSE207422 CNV-only (CopyKAT-style p95 / p90 / mean+2SD / 2-means)

NMPR > MPR on all four cuts (Δ +0.39 to +0.46; exact p=0.073–0.283, n=12). Versus **public** lineage/UMI T/NK the ρ is ~0 to +0.12. Versus **DRMref T/NK** the ρ is negative (p95 ρ=−0.350, p=0.265; mu2 ρ=−0.427, p=0.167). So CNV-aneuploid epithelium recovers the group direction on its own, and recovers both directions only when T/NK is the third-party label.

### GSE241934 IIT (EGFR-mut trial, 11 patients: NMPR 7 / MPR 4)

Author `Epi`, marker epithelium, and their intersection: **ρ negative** (author Epi vs author T/NK: ρ=−0.073, p=0.832) but **NMPR < MPR** (1.590 vs 1.778, p=0.315). CNV-p95 on author/marker Epi: neither half (NMPR < MPR; ρ positive).

No IIT definition recovers both directions.

### GSE241934 real-world EGFR-WT (34 tumors; 33 with author Epi: NMPR 20 / MPR+pCR 13)

Author `Epi` vs author T/NK: n=33, mean NMPR 1.515 vs MPR 1.524, p=0.927, ρ=−0.106, p=0.557 (ρ-only). Most marker and CNV-p95 rows are NMPR < MPR.

The only Real row that is formally “both” is `cnv_gmm_on_authorEpi` (n=32, 20/12): Δ=+0.016, p=0.892, ρ=−0.144, p=0.430. That is a near-tie, not a usable NMPR>MPR split.

`intersect_markerNNL_cnvP95` on Real is the opposite group direction (NMPR 1.840 vs MPR 2.399, p=0.014, n=29).

---

## Compact lead table (preferred public T/NK)

GSE207422: lineage T/NK. GSE241934: author T/NK. Metric = mean log1p(CP10k). min malignant cells = 0.

| Dataset | Cohort | Malignant def | n | NMPR | MPR | Δ | p | ρ | p(ρ) | direction |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE207422 | post_n12 | `thirdparty_drmref` | 12 | 1.571 | 1.120 | +0.450 | 0.154 | −0.126 | 0.697 | both |
| GSE207422 | post_n12 | `marker_all_epi` | 12 | 1.598 | 1.256 | +0.342 | 0.109 | −0.098 | 0.762 | both |
| GSE207422 | post_n12 | `marker_epi_nnl_p75` | 10 | 1.950 | 0.971 | +0.978 | 0.183 | −0.394 | 0.260 | both |
| GSE207422 | post_n12 | `marker_epi_nnl_zero` | 8 | 1.964 | 0.342 | +1.622 | 0.071 | −0.238 | 0.570 | both |
| GSE207422 | post_n12 | `cnv_p95` | 12 | 1.864 | 1.453 | +0.411 | 0.154 | +0.077 | 0.812 | NMPR>MPR only |
| GSE207422 | post_n12 | `cnv_p90` | 12 | 1.851 | 1.395 | +0.456 | 0.154 | +0.063 | 0.846 | NMPR>MPR only |
| GSE207422 | post_n12 | `cnv_mu2` | 12 | 1.831 | 1.445 | +0.386 | 0.283 | +0.028 | 0.931 | NMPR>MPR only |
| GSE207422 | post_n12 | `cnv_gmm` | 12 | 1.818 | 1.372 | +0.446 | 0.073 | 0.000 | 1.000 | NMPR>MPR only |
| GSE207422 | post_n12 | `intersect_markerNNL_cnvP95` | 8 | 2.074 | 0.907 | +1.167 | 0.250 | −0.310 | 0.456 | both |
| GSE207422 | post_n12 | `intersect_drmref_cnvP95` | 12 | 2.060 | 1.734 | +0.326 | 0.283 | +0.182 | 0.572 | NMPR>MPR only |
| GSE207422 | post_n12 | `intersect_drmref_allEpi` | 12 | 1.639 | 1.350 | +0.288 | 0.214 | +0.112 | 0.729 | NMPR>MPR only |
| GSE241934 | IIT_EGFRmut | `author_epi` | 11 | 1.590 | 1.778 | −0.188 | 0.315 | −0.073 | 0.832 | ρ<0 only |
| GSE241934 | IIT_EGFRmut | `marker_all_epi` | 11 | 1.555 | 1.764 | −0.209 | 0.315 | −0.073 | 0.832 | ρ<0 only |
| GSE241934 | IIT_EGFRmut | `cnv_p95_on_authorEpi` | 11 | 1.985 | 2.056 | −0.071 | 0.648 | +0.255 | 0.450 | neither |
| GSE241934 | REAL_EGFRWT | `author_epi` | 33 | 1.515 | 1.524 | −0.008 | 0.927 | −0.106 | 0.557 | ρ<0 only |
| GSE241934 | REAL_EGFRWT | `marker_all_epi` | 33 | 1.466 | 1.520 | −0.054 | 0.699 | −0.066 | 0.716 | ρ<0 only |
| GSE241934 | REAL_EGFRWT | `cnv_p95_on_authorEpi` | 30 | 1.752 | 2.166 | −0.415 | 0.082 | −0.099 | 0.604 | ρ<0 only |
| GSE241934 | REAL_EGFRWT | `cnv_gmm_on_authorEpi` | 32 | 1.703 | 1.687 | +0.016 | 0.892 | −0.144 | 0.430 | both (near-tie) |

---

## What to take from the grid

On **GSE207422**, the definitions that recover the user direction on a complete n=12 with **public** T/NK are **marker epithelium** (`marker_all_epi`) and **third-party DRMref malignant**. CNV-aneuploid epithelium recovers NMPR>MPR at the same n and the same p≈0.15 as DRMref, and recovers the negative-ρ half when T/NK is DRMref, not when T/NK is marker/UMI.

On **GSE241934**, author epithelium (the deposited CopyKAT-style Epi compartment) recovers a weak negative ρ and does not recover NMPR>MPR in either the 11-patient IIT or the 33-patient real-world set.

n=11–12 is small: a true ρ=−0.45 has two-sided Spearman p≈0.14 at n=12. Direction and p are reported separately. Every T/NK pairing, %positive, and pseudobulk row is in the per-dataset `grid_stats.tsv` files.

## Files

- `GSE207422/` — cell calls, per-sample scores, grid, figures
- `GSE241934/` — IIT + Real cell calls, per-sample scores, grid
- `grid_stats_primary.tsv` / `direction_lead.tsv` — combined primary rows
- `methods/scrna_cnv_grid/METHODS.md` — methods only
