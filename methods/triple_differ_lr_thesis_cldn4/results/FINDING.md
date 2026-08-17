# FINDING — CLDN4-only thesis-aligned ligand test on the triple that differs

ADDITIVE. **CLDN4 only. No dual-high.** Patient/donor/locked-sample is the unit
(GSE123902 donor; GSE131907 locked sample; GSE205335 patient from PR #459).
Do **not** add GSE148071. Do **not** pile the 7-pool. Do **not** re-audit the T/NK Spearman.

Thesis (already correct; this folder tests it on the given cut):

- CLDN4-high → more barrier/inhibitory outgoing to T/NK (F11R, NECTIN2–TIGIT, CDH1, LGALS9)
- CLDN4-low / KD-like → more IFN / T-recruit / MHC-I outgoing (CXCL9/10–CXCR3, CCL5, HLA–CD8)

The triple that **differs** is taken as given from PR #459:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE131907+GSE205335 | %pos | 56 | −0.522 (7.2e-05, 0%) | −0.735 (0.001; 14 vs 14) |

Singles (given, same PR): GSE123902 n=13 ρ=−0.659; GSE131907 n=21 ρ=−0.522; GSE205335 n=22 ρ=−0.435.

This folder scores **both pre-specified families** as CellChat-style outgoing
Hill probabilities (Jin et al. 2021; 10% truncated mean, Kh=0.5, expr_prop≥0.1)
from CLDN4-high vs CLDN4-low malignant cells to same-unit T/NK.
CellChat R was not run. Pairs were not discovered.

## Honest n

Given combo n=56 is **not** the ligand n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (median: ≥20
malignant and ≥10/arm plus ≥10 T/NK; Q4 vs Q1:
n_mal≥40 and ≥10/arm).

| gate | n | note |
|---|---:|---|
| Given triple (do not re-audit) | 56 | 13 donors + 21 samples + 22 patients |
| Locked GSE123902 donors | 13 | marker-malignant; PRIMARY preferred over METASTASIS |
| Locked GSE131907 samples | 21 | author Malignant cells / tS*; T lymphocytes + NK cells |
| Locked GSE205335 patients | 22 | author lineage.sub Malignant cells; lineage.total T/NK |
| Median-split LR | **53** | primary ligand n |
| Tertile extra | 55 | exclusive arms |
| Q4 vs Q1 extra | 55 | n_mal≥40 |

GSE123902 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0).
GSE131907 malignant = author `Cell_subtype` in {Malignant cells, tS1, tS2, tS3}.
GSE205335 malignant = author `lineage.sub` == Malignant cells.
TACSTD2 is never a gate.

## Family 1 — barrier / inhibitory (expect CLDN4-high > low)

Primary table: `results/family_barrier_inhibitory.tsv`.

| pair | axis | expect | n | 123902 | 131907 | 205335 | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---|---|---|
| JAM1_ITGAL_ITGB2 | F11R | high>low | 34 | 5 | 12 | 17 | +0.166 | 1.16e-10 | high>low | yes |
| NECTIN2_TIGIT | NECTIN2-TIGIT | high>low | 28 | 9 | 0 | 19 | +0.125 | 7.45e-09 | high>low | yes |
| CDH1_ITGAE_ITGB7 | CDH1 | high>low | 39 | 6 | 15 | 18 | +0.084 | 3.64e-12 | high>low | yes |
| CDH1_KLRG1 | CDH1 | high>low | 29 | 2 | 10 | 17 | +0.070 | 3.73e-09 | high>low | yes |
| LGALS9_HAVCR2 | LGALS9 | high>low | 24 | 6 | 9 | 9 | +0.017 | 2.62e-04 | high>low | yes |
| LGALS9_CD44 | LGALS9 | high>low | 44 | 10 | 16 | 18 | +0.056 | 8.70e-05 | high>low | yes |
| LGALS9_CD45 | LGALS9 | high>low | 43 | 10 | 16 | 17 | +0.073 | 1.21e-04 | high>low | yes |
| **FAMILY** barrier_inhibitory | family | high>low | 50 | 11 | 18 | 21 | +0.077 | 5.45e-13 | high>low | yes |

Family aggregate: n=50 mean ΔP=+0.077 p_W=5.45e-13 agrees=yes.

## Family 2 — IFN / T-recruit / MHC-I (expect CLDN4-low / KD-like > high)

Primary table: `results/family_ifn_recruit_mhci.tsv`.

| pair | axis | expect | n | 123902 | 131907 | 205335 | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---|---|---|
| CXCL9_CXCR3 | CXCL9/10-CXCR3 | low>high | 2 | 1 | 1 | 0 | +0.004 | 0.5 | high>low | thin |
| CXCL10_CXCR3 | CXCL9/10-CXCR3 | low>high | 10 | 3 | 4 | 3 | -0.001 | 0.496 | low>high | direction_only |
| CCL5_CCR5 | CCL5 | low>high | 4 | 0 | 2 | 2 | +0.003 | 1 | high>low | no |
| CCL5_CCR1 | CCL5 | low>high | 3 | 0 | 2 | 1 | +0.004 | 0.5 | high>low | no |
| HLA-A_CD8A | HLA-CD8 | low>high | 46 | 10 | 17 | 19 | +0.027 | 3.03e-05 | high>low | opposite |
| HLA-B_CD8A | HLA-CD8 | low>high | 46 | 10 | 17 | 19 | +0.031 | 6.66e-05 | high>low | opposite |
| HLA-C_CD8A | HLA-CD8 | low>high | 46 | 10 | 17 | 19 | +0.036 | 1.63e-04 | high>low | opposite |
| **FAMILY** ifn_recruit_mhci | family | low>high | 47 | 11 | 17 | 19 | +0.028 | 1.34e-04 | high>low | opposite |

Family aggregate: n=47 mean ΔP=+0.028 p_W=1.34e-04 agrees=opposite.

## Extra figures

- `figures/fig_given_combo_rho.png` — given %pos ρ (not re-audited)
- `figures/fig_honest_n.png` — given n vs ligand-eligible n
- `figures/fig_n_per_unit.png` — per-unit malignant / T/NK floors
- `figures/fig_family_barrier_bars.png` / `fig_family_barrier_forest.png`
- `figures/fig_family_ifn_bars.png` / `fig_family_ifn_forest.png`
- `figures/fig_family_barrier_units.png` / `fig_family_ifn_units.png`
- `figures/fig_family_by_cohort.png` — extra: family ΔP split by cohort

## What is not claimed

- The n=56 ρ=−0.522 / Q4 r=−0.735 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071 is not added. This is not the 7-pool.
- Cell-pooled permutations are not the test. Patient/donor is the unit.
- This is not a CellChat discovery screen. Only the two pre-specified families.

## Reproduce

```bash
python3 methods/triple_differ_lr_thesis_cldn4/scripts/download.py
python3 methods/triple_differ_lr_thesis_cldn4/scripts/analyze.py
```

Hill constants: trim=0.1, Kh=0.5, expr_prop=0.1.
MIN_CELLS_COMP=20 is documented but the primary gate is the 10-cell arm.
