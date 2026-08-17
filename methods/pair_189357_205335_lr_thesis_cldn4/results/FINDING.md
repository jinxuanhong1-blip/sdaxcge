# FINDING — CLDN4-only thesis-aligned ligand test on GSE189357 + GSE205335

**Verdict:** On 31 median-paired patients (9 GSE189357 + 22 GSE205335), the barrier/inhibitory family is **7/7 pairs in the thesis direction** (high > low; family n=30, mean ΔP=+0.090, p=3.86e-07). The IFN/T-recruit/MHC-I family is **not**: CXCL9/10 and CCL5 are sparse (honest n=0–3); classical HLA–CD8 is slightly higher from CLDN4-high (n=27, mean ΔP≈+0.02, p>0.05). Barrier-up-in-high is a primary arm, not a recruit leftover. Scores are CellChat-style Hill *P*, not secretion or contact.

ADDITIVE. **CLDN4 only. No dual-high.** Patient is the unit
(GSE189357 TD1–TD9 marker-malignant; GSE205335 locked author-malignant patients).
Do **not** re-audit the T/NK Spearman. p-values are descriptive.

Thesis (already correct; this folder tests it on the given cut):

- CLDN4-high → more barrier/inhibitory outgoing to T/NK (F11R, NECTIN2–TIGIT, CDH1, LGALS9)
- CLDN4-low / KD-like → more IFN / T-recruit / MHC-I outgoing (CXCL9/10–CXCR3, CCL5, HLA–CD8)

The pair that **differs** is taken as given from PR #459:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE189357+GSE205335 | %pos | 31 | −0.478 (0.009, 0%) | −0.750 (0.010; 8 vs 8) |

Singles (given, same PR): GSE189357 n=9 ρ=−0.600; GSE205335 n=22 ρ=−0.435.

This folder scores **both pre-specified families** as CellChat-style outgoing
Hill probabilities (Jin et al. 2021; 10% truncated mean, Kh=0.5, expr_prop≥0.1)
from CLDN4-high vs CLDN4-low malignant cells to same-patient T/NK.
CellChat R was not run. Pairs were not discovered.

## Honest n

Given combo n=31 is **not** the ligand n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (median: ≥20
malignant and ≥10/arm plus ≥10 T/NK; Q4 vs Q1:
n_mal≥40 and ≥10/arm).

| gate | n | note |
|---|---:|---|
| Given combo (do not re-audit) | 31 | 9 + 22 patients |
| Locked GSE189357 patients | 9 | marker-malignant; TD1–TD9 |
| Locked GSE205335 patients | 22 | author Malignant cells; T/NK cells |
| Median-split LR | **31** | primary ligand n |
| Tertile extra | 30 | exclusive arms |
| Q4 vs Q1 extra | 30 | n_mal≥40 |

GSE189357 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0).
GSE205335 malignant = author `lineage.sub` == `Malignant cells`.
TACSTD2 is never a gate.

## Family 1 — barrier / inhibitory (expect CLDN4-high > low)

Primary table: `results/family_barrier_inhibitory.tsv`.

| pair | axis | expect | n | n_189357/n_205335 | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---|---:|---|---|---|
| JAM1_ITGAL_ITGB2 | F11R | high>low | 26 | 9/17 | +0.213 | 2.98e-08 | high>low | yes |
| NECTIN2_TIGIT | NECTIN2-TIGIT | high>low | 27 | 8/19 | +0.115 | 1.49e-08 | high>low | yes |
| CDH1_ITGAE_ITGB7 | CDH1 | high>low | 22 | 4/18 | +0.111 | 4.77e-07 | high>low | yes |
| CDH1_KLRG1 | CDH1 | high>low | 21 | 4/17 | +0.067 | 9.54e-07 | high>low | yes |
| LGALS9_HAVCR2 | LGALS9 | high>low | 11 | 2/9 | +0.030 | 0.00195 | high>low | yes |
| LGALS9_CD44 | LGALS9 | high>low | 24 | 6/18 | +0.058 | 0.034 | high>low | yes |
| LGALS9_CD45 | LGALS9 | high>low | 23 | 6/17 | +0.077 | 0.0384 | high>low | yes |
| **FAMILY** barrier_inhibitory | family | high>low | 30 | 9/21 | +0.090 | 3.86e-07 | high>low | yes |

Family aggregate: n=30 mean ΔP=+0.090 p_W=3.86e-07 agrees=yes.

## Family 2 — IFN / T-recruit / MHC-I (expect CLDN4-low / KD-like > high)

Primary table: `results/family_ifn_recruit_mhci.tsv`.

| pair | axis | expect | n | n_189357/n_205335 | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---|---:|---|---|---|
| CXCL9_CXCR3 | CXCL9/10-CXCR3 | low>high | 0 | 0/0 | NA | NA | flat/NA | thin |
| CXCL10_CXCR3 | CXCL9/10-CXCR3 | low>high | 3 | 0/3 | -0.015 | 1 | low>high | direction_only |
| CCL5_CCR5 | CCL5 | low>high | 2 | 0/2 | -0.003 | 0.5 | low>high | thin |
| CCL5_CCR1 | CCL5 | low>high | 1 | 0/1 | -0.002 | NA | low>high | thin |
| HLA-A_CD8A | HLA-CD8 | low>high | 27 | 8/19 | +0.023 | 0.0954 | high>low | no |
| HLA-B_CD8A | HLA-CD8 | low>high | 27 | 8/19 | +0.021 | 0.17 | high>low | no |
| HLA-C_CD8A | HLA-CD8 | low>high | 27 | 8/19 | +0.021 | 0.279 | high>low | no |
| **FAMILY** ifn_recruit_mhci | family | low>high | 27 | 8/19 | +0.020 | 0.178 | high>low | no |

Family aggregate: n=27 mean ΔP=+0.020 p_W=0.178 agrees=no.

## Extra figures

- `figures/fig_given_combo_rho.png` — given %pos ρ (not re-audited)
- `figures/fig_honest_n.png` — given n vs ligand-eligible n
- `figures/fig_n_per_unit.png` — per-unit malignant / T/NK floors
- `figures/fig_family_barrier_bars.png` / `fig_family_barrier_forest.png`
- `figures/fig_family_ifn_bars.png` / `fig_family_ifn_forest.png`
- `figures/fig_family_barrier_units.png` / `fig_family_ifn_units.png`

## What is not claimed

- The n=31 ρ=−0.478 / Q4 r=−0.750 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- Cell-pooled permutations are not the test. Patient is the unit.
- This is not a CellChat discovery screen. Only the two pre-specified families.
- Malignant definitions differ (marker vs author) and are stated on every row.

## Reproduce

```bash
python3 methods/pair_189357_205335_lr_thesis_cldn4/scripts/download.py
python3 methods/pair_189357_205335_lr_thesis_cldn4/scripts/analyze.py
```

Hill constants: trim=0.1, Kh=0.5, expr_prop=0.1.
