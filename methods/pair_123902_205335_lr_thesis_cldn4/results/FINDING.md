# FINDING — CLDN4-only thesis-aligned ligand test on GSE123902 + GSE205335

ADDITIVE. **CLDN4 only. No dual-high.** Patient/donor is the unit
(GSE123902 donor; GSE205335 locked patient from PR #459).
Do **not** add GSE148071. Do **not** re-audit the T/NK Spearman.

Thesis (already correct; this folder tests it on the given cut):

- CLDN4-high → more barrier/inhibitory outgoing to T/NK (F11R, NECTIN2–TIGIT, CDH1, LGALS9)
- CLDN4-low / KD-like → more IFN / T-recruit / MHC-I outgoing (CXCL9/10–CXCR3, CCL5, HLA–CD8)

The pair that **differs** is taken as given from PR #459:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE205335 | %pos | 35 | −0.522 (0.002, 0%) | −0.802 (0.005; 9 vs 9) |

Singles (given, same PR): GSE123902 n=13 ρ=−0.659; GSE205335 n=22 ρ=−0.435.

This folder scores **both pre-specified families** as CellChat-style outgoing
Hill probabilities (Jin et al. 2021; 10% truncated mean, Kh=0.5, expr_prop≥0.1)
from CLDN4-high vs CLDN4-low malignant cells to same-unit T/NK.
CellChat R was not run. Pairs were not discovered.

## Honest n

Given combo n=35 is **not** the ligand n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (median: ≥20
malignant and ≥10/arm plus ≥10 T/NK; Q4 vs Q1:
n_mal≥40 and ≥10/arm).

| gate | n | note |
|---|---:|---|
| Given combo (do not re-audit) | 35 | 13 donors + 22 patients |
| Locked GSE123902 donors | 13 | marker-malignant; PRIMARY preferred over METASTASIS |
| Locked GSE205335 patients | 22 | author Malignant cells; T/NK cells |
| Median-split LR | **34** | primary ligand n |
| Tertile extra | 34 | exclusive arms |
| Q4 vs Q1 extra | 34 | n_mal≥40 |

GSE123902 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0).
GSE205335 malignant = author `lineage.sub` == Malignant cells.
TACSTD2 is never a gate.

## Family 1 — barrier / inhibitory (expect CLDN4-high > low)

Primary table: `results/family_barrier_inhibitory.tsv`.

| pair | axis | expect | n | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---|---|---|
| JAM1_ITGAL_ITGB2 | F11R | high>low | 22 | +0.206 | 4.77e-07 | high>low | yes |
| NECTIN2_TIGIT | NECTIN2-TIGIT | high>low | 28 | +0.125 | 7.45e-09 | high>low | yes |
| CDH1_ITGAE_ITGB7 | CDH1 | high>low | 24 | +0.104 | 1.19e-07 | high>low | yes |
| CDH1_KLRG1 | CDH1 | high>low | 19 | +0.072 | 3.81e-06 | high>low | yes |
| LGALS9_HAVCR2 | LGALS9 | high>low | 15 | +0.022 | 0.00429 | high>low | yes |
| LGALS9_CD44 | LGALS9 | high>low | 28 | +0.055 | 0.0116 | high>low | yes |
| LGALS9_CD45 | LGALS9 | high>low | 27 | +0.074 | 0.0164 | high>low | yes |
| **FAMILY** barrier_inhibitory | family | high>low | 32 | +0.087 | 2.00e-08 | high>low | yes |

Family aggregate: n=32 mean ΔP=+0.087 p_W=2.00e-08 agrees=yes.

## Family 2 — IFN / T-recruit / MHC-I (expect CLDN4-low / KD-like > high)

Primary table: `results/family_ifn_recruit_mhci.tsv`.

| pair | axis | expect | n | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---|---|---|
| CXCL9_CXCR3 | CXCL9/10-CXCR3 | low>high | 1 | +0.006 | NA | high>low | thin |
| CXCL10_CXCR3 | CXCL9/10-CXCR3 | low>high | 6 | -0.003 | 0.625 | low>high | direction_only |
| CCL5_CCR5 | CCL5 | low>high | 2 | -0.003 | 0.5 | low>high | thin |
| CCL5_CCR1 | CCL5 | low>high | 1 | -0.002 | NA | low>high | thin |
| HLA-A_CD8A | HLA-CD8 | low>high | 29 | +0.026 | 0.0148 | high>low | opposite |
| HLA-B_CD8A | HLA-CD8 | low>high | 29 | +0.027 | 0.0243 | high>low | opposite |
| HLA-C_CD8A | HLA-CD8 | low>high | 29 | +0.030 | 0.0386 | high>low | opposite |
| **FAMILY** ifn_recruit_mhci | family | low>high | 30 | +0.024 | 0.0282 | high>low | opposite |

Family aggregate: n=30 mean ΔP=+0.024 p_W=0.0282 agrees=opposite.

## Extra figures

- `figures/fig_given_combo_rho.png` — given %pos ρ (not re-audited)
- `figures/fig_honest_n.png` — given n vs ligand-eligible n
- `figures/fig_n_per_unit.png` — per-unit malignant / T/NK floors
- `figures/fig_family_barrier_bars.png` / `fig_family_barrier_forest.png`
- `figures/fig_family_ifn_bars.png` / `fig_family_ifn_forest.png`
- `figures/fig_family_barrier_units.png` / `fig_family_ifn_units.png`

## What is not claimed

- The n=35 ρ=−0.522 / Q4 r=−0.802 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071 is not added.
- Cell-pooled permutations are not the test. Patient/donor is the unit.
- This is not a CellChat discovery screen. Only the two pre-specified families.

## Reproduce

```bash
python3 methods/pair_123902_205335_lr_thesis_cldn4/scripts/download.py
python3 methods/pair_123902_205335_lr_thesis_cldn4/scripts/analyze.py
```

Hill constants: trim=0.1, Kh=0.5, expr_prop=0.1.
