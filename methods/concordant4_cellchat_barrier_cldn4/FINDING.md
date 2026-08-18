# FINDING — Concordant-four CellChat barrier call (CLDN4-only)

ADDITIVE. **Thesis already correct. Ligands stay.**
CLDN4 only. No dual-high. Concordant four only
(GSE123902 + GSE131907 + GSE205335 + GSE189357).
Do **not** add GSE148071 / GSE127465 / CD45-only. This is **not** a full-pool.

Patient/donor/locked-sample is the unit. Honest ligand n = units with
**both compartments** (malignant/epithelium high+low arms and T/NK).

Thesis (already correct; this folder reports the four-set pooled ligand call):

- Barrier/inhibitory outgoing **UP from CLDN4-high** (F11R, NECTIN2–TIGIT, CDH1, LGALS9) is **ON-thesis**.
- IFN/T-recruit outgoing UP from CLDN4-low is the KD-like arm (often **not** detected; CXCL9/10 sparse).
  Do not bury barrier-up-in-high as a recruit-up skip.

Primary split is malignant **Q4 vs Q1**. Extra: %pos and median.

## Honest n

Locked four n=65 (13+21+22+9) is **not** the ligand n.

| gate | n | note |
|---|---:|---|
| Locked four (given PR #459 singles) | 65 | 13+21+22+9; T/NK ρ not re-audited |
| Locked GSE123902 donors | 13 | marker-malignant / epithelium proxy; PRIMARY preferred |
| Locked GSE131907 samples | 21 | author Malignant cells / tS*; T lymphocytes + NK cells |
| Locked GSE205335 patients | 22 | author lineage.sub Malignant cells; lineage.total T/NK |
| Locked GSE189357 patients | 9 | marker-malignant / epithelium proxy; TD1–TD9 |
| Both compartments (n_mal≥10, n_tnk≥20) | 65 | honest inventory |
| Q4 vs Q1 LR | **64** | primary ligand n (n_mal≥40, ≥10/arm, T/NK≥20) |
| %pos LR extra | 62 | CLDN4>0 vs =0 |
| Median-split LR extra | 62 | |
| Tertile extra | 64 | exclusive arms |

GSE123902 / GSE189357 malignant labels are thin → epithelium marker-malignant
(EPCAM\|KRT8\|KRT18\|KRT19 > 0 and PTPRC == 0) → T/NK.
GSE131907 / GSE205335 use author malignant → T/NK.
TACSTD2 is never a gate. PVRL2 is aliased to NECTIN2 on GSE131907.

## Primary — barrier / inhibitory family ΔP (high − low)

ON-thesis. Primary table: `results/family_barrier_inhibitory.tsv` (Q4 vs Q1).

| pair | axis | expect | n | 123902 | 131907 | 205335 | 189357 | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| JAM1_ITGAL_ITGB2 | F11R | high>low | 44 | 6 | 13 | 16 | 9 | +0.188 | 1.14e-13 | high>low | yes |
| NECTIN2_TIGIT | NECTIN2-TIGIT | high>low | 55 | 10 | 19 | 18 | 8 | +0.127 | 1.11e-10 | high>low | yes |
| CDH1_ITGAE_ITGB7 | CDH1 | high>low | 42 | 6 | 15 | 17 | 4 | +0.101 | 4.55e-13 | high>low | yes |
| CDH1_KLRG1 | CDH1 | high>low | 32 | 2 | 10 | 16 | 4 | +0.082 | 4.66e-10 | high>low | yes |
| LGALS9_HAVCR2 | LGALS9 | high>low | 27 | 8 | 9 | 9 | 1 | +0.023 | 0.00542 | high>low | yes |
| LGALS9_CD44 | LGALS9 | high>low | 49 | 11 | 16 | 17 | 5 | +0.097 | 1.02e-06 | high>low | yes |
| LGALS9_CD45 | LGALS9 | high>low | 48 | 11 | 16 | 16 | 5 | +0.124 | 7.19e-07 | high>low | yes |
| **FAMILY** barrier_inhibitory | family | high>low | 61 | 12 | 20 | 20 | 9 | +0.099 | 3.99e-11 | high>low | yes |

Family aggregate: n=61 mean ΔP=+0.099 p_W=3.99e-11 agrees=yes.

Called-out axes (F11R, NECTIN2–TIGIT, CDH1, LGALS9) are the high-side barrier ligands.

## KD-like arm — IFN / T-recruit / MHC-I (expect low > high)

Do not file the barrier result as a recruit-up skip.

CXCL9/10 (KD-like recruit arm): CXCL9_CXCR3 undetected at expr_prop≥0.1 (honest n=1; units with any detect flag=1); CXCL10_CXCR3 n=10 mean ΔP=-0.017 p_W=0.922.

Primary table: `results/family_ifn_recruit_mhci.tsv`.

| pair | axis | expect | n | 123902 | 131907 | 205335 | 189357 | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| CXCL9_CXCR3 | CXCL9/10-CXCR3 | low>high | 1 | 1 | 0 | 0 | 0 | +0.001 | NA | high>low | thin |
| CXCL10_CXCR3 | CXCL9/10-CXCR3 | low>high | 10 | 2 | 5 | 3 | 0 | -0.017 | 0.922 | low>high | direction_only |
| CCL5_CCR5 | CCL5 | low>high | 4 | 0 | 2 | 2 | 0 | -0.005 | 0.875 | low>high | direction_only |
| CCL5_CCR1 | CCL5 | low>high | 3 | 0 | 2 | 1 | 0 | +0.004 | 0.75 | high>low | no |
| HLA-A_CD8A | HLA-CD8 | low>high | 56 | 11 | 19 | 18 | 8 | +0.039 | 6.21e-06 | high>low | opposite |
| HLA-B_CD8A | HLA-CD8 | low>high | 56 | 11 | 19 | 18 | 8 | +0.037 | 4.22e-05 | high>low | opposite |
| HLA-C_CD8A | HLA-CD8 | low>high | 56 | 11 | 19 | 18 | 8 | +0.041 | 1.93e-04 | high>low | opposite |
| **FAMILY** ifn_recruit_mhci | family | low>high | 56 | 11 | 19 | 18 | 8 | +0.036 | 5.99e-05 | high>low | opposite |

Family aggregate: n=56 mean ΔP=+0.036 p_W=5.99e-05 agrees=opposite.

## Extra: median-split barrier family (same units, different gate)

| pair | axis | expect | n | 123902 | 131907 | 205335 | 189357 | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| JAM1_ITGAL_ITGB2 | F11R | high>low | 43 | 5 | 12 | 17 | 9 | +0.163 | 2.27e-13 | high>low | yes |
| NECTIN2_TIGIT | NECTIN2-TIGIT | high>low | 53 | 9 | 17 | 19 | 8 | +0.097 | 2.39e-10 | high>low | yes |
| CDH1_ITGAE_ITGB7 | CDH1 | high>low | 43 | 6 | 15 | 18 | 4 | +0.081 | 2.27e-13 | high>low | yes |
| CDH1_KLRG1 | CDH1 | high>low | 33 | 2 | 10 | 17 | 4 | +0.065 | 2.33e-10 | high>low | yes |
| LGALS9_HAVCR2 | LGALS9 | high>low | 26 | 6 | 9 | 9 | 2 | +0.016 | 1.40e-04 | high>low | yes |
| LGALS9_CD44 | LGALS9 | high>low | 50 | 10 | 16 | 18 | 6 | +0.063 | 2.64e-05 | high>low | yes |
| LGALS9_CD45 | LGALS9 | high>low | 49 | 10 | 16 | 17 | 6 | +0.081 | 3.13e-05 | high>low | yes |
| **FAMILY** barrier_inhibitory | family | high>low | 59 | 11 | 18 | 21 | 9 | +0.078 | 1.55e-10 | high>low | yes |

## Comparison rows (PR #424 / #474 — not re-audited)

The new result is the four-set pooled ligand call. These rows are context only.

| source | combo | row | split | n | ΔP | p | note |
|---|---|---|---|---:|---:|---|---|
| PR #424 | GSE131907+GSE205335 | NECTIN2–TIGIT | median | 39 | +0.094 | 3.6e-12 | given |
| PR #474 | GSE123902+GSE205335 | FAMILY barrier | median | 32 | +0.087 | 2.00e-08 | given |
| this | concordant four | FAMILY barrier | Q4 vs Q1 | 61 | +0.099 | 3.99e-11 | new |

## Extra figures

- `figures/fig_given_singles_rho.png` — PR #459 singles (not re-audited)
- `figures/fig_honest_n.png` — given n vs ligand-eligible n
- `figures/fig_n_per_unit.png` — per-unit malignant / T/NK floors
- `figures/fig_family_barrier_bars.png` / `fig_family_barrier_forest.png`
- `figures/fig_called_out_barrier.png` — F11R, NECTIN2–TIGIT, CDH1, LGALS9
- `figures/fig_family_ifn_bars.png` / `fig_family_ifn_forest.png`
- `figures/fig_family_barrier_units.png` / `fig_family_ifn_units.png`
- `figures/fig_family_by_cohort.png`
- `figures/fig_comparison_424_474.png` — comparison rows vs this call

## What is not claimed

- PR #424 / #474 numbers are comparison rows. They were not recomputed.
- PR #459 T/NK Spearmans are given. They were not re-audited.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071, GSE127465, and CD45-only libraries are not added.
- This is not a CellChat discovery screen and not a 7-pool.
- CellChat R / LIANA R were not run. Probability is Jin et al. 2021
  (10% truncated mean, Kh=0.5, expr_prop≥0.1).

## Reproduce

```bash
python3 methods/concordant4_cellchat_barrier_cldn4/scripts/download.py
python3 methods/concordant4_cellchat_barrier_cldn4/scripts/analyze.py
```

Hill constants: trim=0.1, Kh=0.5, expr_prop=0.1.
Arm floor ≥10; T/NK compartment floor ≥20; Q4 vs Q1 also n_mal≥40.
