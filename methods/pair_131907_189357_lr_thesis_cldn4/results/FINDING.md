# FINDING — CLDN4-only thesis-aligned ligand test on GSE131907 + GSE189357

ADDITIVE. **CLDN4 only. No dual-high.** Patient is the unit
(GSE131907 locked sample from PR #459; GSE189357 patient TD1–TD9).
Do **not** add GSE148071. Do **not** re-audit the T/NK Spearman.

Thesis (already correct; this folder tests it on the given cut):

- CLDN4-high → more barrier/inhibitory outgoing to T/NK (F11R, NECTIN2–TIGIT, CDH1, LGALS9)
- CLDN4-low / KD-like → more IFN / T-recruit / MHC-I outgoing (CXCL9/10–CXCR3, CCL5, HLA–CD8)

The pair that **differs** is taken as given from PR #459:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE131907+GSE189357 | %pos | 30 | −0.542 (0.00291, 0%) | −0.619 (0.042; 9 vs 7) |

Singles (given, same PR): GSE131907 n=21 ρ=−0.522; GSE189357 n=9 ρ=−0.600.

This folder scores **both pre-specified families** as CellChat-style outgoing
Hill probabilities (Jin et al. 2021; 10% truncated mean, Kh=0.5, expr_prop≥0.1)
from CLDN4-high vs CLDN4-low malignant cells to same-unit T/NK.
CellChat R was not run. Pairs were not discovered.

## Honest n

Given combo n=30 is **not** the ligand n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (median: ≥20
malignant and ≥10/arm plus ≥10 T/NK; Q4 vs Q1:
n_mal≥40 and ≥10/arm).

| gate | n | note |
|---|---:|---|
| Given combo (do not re-audit) | 30 | 21 samples + 9 patients |
| Locked GSE131907 samples | 21 | author Malignant cells / tS*; T lymphocytes + NK cells |
| Locked GSE189357 patients | 9 | marker-malignant; TD1–TD9 |
| Median-split LR | **28** | primary ligand n |
| Tertile extra | 30 | exclusive arms |
| Q4 vs Q1 extra | 30 | n_mal≥40 |

GSE131907 malignant = author `Cell_subtype` in {Malignant cells, tS1, tS2, tS3}.
GSE189357 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0).
TACSTD2 is never a gate.

Median dropped **EBUS_13** and **NS_16** (CLDN4 almost all zero; high arm <10).
Those two still pass tertile/Q4 extras. Ligand n is pair-specific after the
expr_prop≥0.10 detection floor — not 28 for every pair.

## Verdict

On this cut, **family 1 agrees** and **family 2 does not**.

Barrier/inhibitory outgoing is higher from CLDN4-high malignant cells to
same-unit T/NK (family n=27, mean ΔP=+0.066, p_W=1.04e-06). Every
pre-specified pair in that family is high>low.

IFN / T-recruit / MHC-I is **not** higher from CLDN4-low cells. CXCL9/10
and CCL5 are mostly undetected (n=1–4, thin). HLA–CD8 is detected
(n=25) but goes the **opposite** way (CLDN4-high > low). The family
aggregate is therefore opposite, and it is carried by HLA–CD8, not by
the chemokine pairs.

## Family 1 — barrier / inhibitory (expect CLDN4-high > low)

Primary table: `results/family_barrier_inhibitory.tsv`.

| pair | axis | expect | n | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---|---|---|
| JAM1_ITGAL_ITGB2 | F11R | high>low | 21 | +0.118 | 9.54e-07 | high>low | yes |
| NECTIN2_TIGIT | NECTIN2-TIGIT | high>low | 8 | +0.045 | 0.00781 | high>low | yes |
| CDH1_ITGAE_ITGB7 | CDH1 | high>low | 19 | +0.053 | 3.81e-06 | high>low | yes |
| CDH1_KLRG1 | CDH1 | high>low | 14 | +0.056 | 1.22e-04 | high>low | yes |
| LGALS9_HAVCR2 | LGALS9 | high>low | 11 | +0.007 | 0.00684 | high>low | yes |
| LGALS9_CD44 | LGALS9 | high>low | 22 | +0.073 | 1.21e-04 | high>low | yes |
| LGALS9_CD45 | LGALS9 | high>low | 22 | +0.090 | 8.06e-05 | high>low | yes |
| **FAMILY** barrier_inhibitory | family | high>low | 27 | +0.066 | 1.04e-06 | high>low | yes |

Family aggregate: n=27 mean ΔP=+0.066 p_W=1.04e-06 agrees=yes.

## Family 2 — IFN / T-recruit / MHC-I (expect CLDN4-low / KD-like > high)

Primary table: `results/family_ifn_recruit_mhci.tsv`.

| pair | axis | expect | n | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---|---|---|
| CXCL9_CXCR3 | CXCL9/10-CXCR3 | low>high | 1 | +0.002 | NA | high>low | thin |
| CXCL10_CXCR3 | CXCL9/10-CXCR3 | low>high | 4 | +0.002 | 1 | high>low | no |
| CCL5_CCR5 | CCL5 | low>high | 2 | +0.008 | 0.5 | high>low | thin |
| CCL5_CCR1 | CCL5 | low>high | 2 | +0.007 | 0.5 | high>low | thin |
| HLA-A_CD8A | HLA-CD8 | low>high | 25 | +0.023 | 7.15e-04 | high>low | opposite |
| HLA-B_CD8A | HLA-CD8 | low>high | 25 | +0.022 | 0.00418 | high>low | opposite |
| HLA-C_CD8A | HLA-CD8 | low>high | 25 | +0.025 | 0.0173 | high>low | opposite |
| **FAMILY** ifn_recruit_mhci | family | low>high | 25 | +0.021 | 0.00418 | high>low | opposite |

Family aggregate: n=25 mean ΔP=+0.021 p_W=0.00418 agrees=opposite.

## Extra figures

- `figures/fig_given_combo_rho.png` — given %pos ρ (not re-audited)
- `figures/fig_honest_n.png` — given n vs ligand-eligible n
- `figures/fig_n_per_unit.png` — per-unit malignant / T/NK floors
- `figures/fig_family_barrier_bars.png` / `fig_family_barrier_forest.png`
- `figures/fig_family_ifn_bars.png` / `fig_family_ifn_forest.png`
- `figures/fig_family_barrier_units.png` / `fig_family_ifn_units.png`

## What is not claimed

- The n=30 ρ=−0.542 / Q4 r=−0.619 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071 is not added.
- Cell-pooled permutations are not the test. Patient is the unit.
- This is not a CellChat discovery screen. Only the two pre-specified families.

## Reproduce

```bash
python3 methods/pair_131907_189357_lr_thesis_cldn4/scripts/download.py
python3 methods/pair_131907_189357_lr_thesis_cldn4/scripts/analyze.py
```

Hill constants: trim=0.1, Kh=0.5, expr_prop=0.1.
