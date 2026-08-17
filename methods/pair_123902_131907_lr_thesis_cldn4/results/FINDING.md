# FINDING — CLDN4-only thesis-aligned ligand test on GSE123902 + GSE131907

ADDITIVE. **CLDN4 only. No dual-high.** Patient/donor is the unit
(GSE123902 donor; GSE131907 locked sample from PR #459).
Do **not** add GSE148071. Do **not** re-audit the T/NK Spearman.

Thesis (already correct; this folder tests it on the given cut):

- CLDN4-high → more barrier/inhibitory outgoing to T/NK (F11R, NECTIN2–TIGIT, CDH1, LGALS9)
- CLDN4-low / KD-like → more IFN / T-recruit / MHC-I outgoing (CXCL9/10–CXCR3, CCL5, HLA–CD8)

The pair that **differs** is taken as given from PR #459:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE131907 | %pos | 34 | −0.575 (0.00053, 0%) | −0.700 (0.012; 10 vs 8) |

Singles (given, same PR): GSE123902 n=13 ρ=−0.659; GSE131907 n=21 ρ=−0.522.

## Verdict

On this given cut, **family 1 supports the thesis** and **family 2 does not**
as an outgoing malignant → T/NK score.

- Barrier/inhibitory (F11R, NECTIN2–TIGIT, CDH1, LGALS9): every pre-specified
  pair is higher from CLDN4-high malignant cells. Family aggregate n=29,
  mean ΔP=+0.064, Wilcoxon p=2.6×10⁻⁸ (27/29 units positive).
- IFN / T-recruit / MHC-I: CXCL9/10–CXCR3 and CCL5 are **not detected at
  scale** (honest n=2–7). Classical HLA–CD8 is higher in CLDN4-high
  (n=27, opposite of the KD-like arm). Family aggregate is therefore
  opposite (n=28, mean ΔP=+0.026).

This folder scores **both pre-specified families** as CellChat-style outgoing
Hill probabilities (Jin et al. 2021; 10% truncated mean, Kh=0.5, expr_prop≥0.1)
from CLDN4-high vs CLDN4-low malignant cells to same-unit T/NK.
CellChat R was not run. Pairs were not discovered.

## Honest n

Given combo n=34 is **not** the ligand n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (median: ≥20
malignant and ≥10/arm plus ≥10 T/NK; Q4 vs Q1:
n_mal≥40 and ≥10/arm).

| gate | n | note |
|---|---:|---|
| Given combo (do not re-audit) | 34 | 13 donors + 21 samples |
| Locked GSE123902 donors | 13 | marker-malignant; PRIMARY preferred over METASTASIS |
| Locked GSE131907 samples | 21 | author Malignant cells / tS*; T lymphocytes + NK cells |
| Median-split LR | **31** | primary ligand n; LX699, EBUS_13, NS_16 fail (zero-inflated CLDN4) |
| Tertile extra | 34 | zero-inflated units can put almost all cells in both tails |
| Q4 vs Q1 extra | 34 | same caveat; median is the honest primary |

GSE123902 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0).
GSE131907 malignant = author `Cell_subtype` in {Malignant cells, tS1, tS2, tS3}.
TACSTD2 is never a gate.

## Family 1 — barrier / inhibitory (expect CLDN4-high > low)

Primary table: `results/family_barrier_inhibitory.tsv`.

| pair | axis | expect | n | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---|---|---|
| JAM1_ITGAL_ITGB2 | F11R | high>low | 17 | +0.087 | 1.53e-05 | high>low | yes |
| NECTIN2_TIGIT | NECTIN2-TIGIT | high>low | 9 | +0.083 | 0.00391 | high>low | yes |
| CDH1_ITGAE_ITGB7 | CDH1 | high>low | 21 | +0.050 | 9.54e-07 | high>low | yes |
| CDH1_KLRG1 | CDH1 | high>low | 12 | +0.062 | 4.88e-04 | high>low | yes |
| LGALS9_HAVCR2 | LGALS9 | high>low | 15 | +0.006 | 0.0103 | high>low | yes |
| LGALS9_CD44 | LGALS9 | high>low | 26 | +0.067 | 5.76e-05 | high>low | yes |
| LGALS9_CD45 | LGALS9 | high>low | 26 | +0.084 | 5.76e-05 | high>low | yes |
| **FAMILY** barrier_inhibitory | family | high>low | 29 | +0.064 | 2.61e-08 | high>low | yes |

Family aggregate: n=29 mean ΔP=+0.064 p_W=2.61e-08 agrees=yes.
NECTIN2–TIGIT is detected in GSE123902 only (n=9). LGALS9–HAVCR2 is small
(+0.006) but still high>low.

## Family 2 — IFN / T-recruit / MHC-I (expect CLDN4-low / KD-like > high)

Primary table: `results/family_ifn_recruit_mhci.tsv`.

| pair | axis | expect | n | mean ΔP | p_W | observed | agrees |
|---|---|---|---:|---:|---|---|---|
| CXCL9_CXCR3 | CXCL9/10-CXCR3 | low>high | 2 | +0.004 | 0.5 | high>low | thin |
| CXCL10_CXCR3 | CXCL9/10-CXCR3 | low>high | 7 | +0.005 | 0.562 | high>low | no |
| CCL5_CCR5 | CCL5 | low>high | 2 | +0.008 | 0.5 | high>low | thin |
| CCL5_CCR1 | CCL5 | low>high | 2 | +0.007 | 0.5 | high>low | thin |
| HLA-A_CD8A | HLA-CD8 | low>high | 27 | +0.026 | 2.59e-05 | high>low | opposite |
| HLA-B_CD8A | HLA-CD8 | low>high | 27 | +0.028 | 1.25e-04 | high>low | opposite |
| HLA-C_CD8A | HLA-CD8 | low>high | 27 | +0.034 | 2.36e-04 | high>low | opposite |
| **FAMILY** ifn_recruit_mhci | family | low>high | 28 | +0.026 | 2.37e-04 | high>low | opposite |

Family aggregate: n=28 mean ΔP=+0.026 p_W=2.37e-04 agrees=opposite.
CXCL9/10 and CCL5 fail the 10% expr_prop floor in most units. The family
row is carried by HLA–CD8, which is **not** a recruit-down copy of the
given T/NK Spearman.

## Extra figures

- `figures/fig_given_combo_rho.png` — given %pos ρ (not re-audited)
- `figures/fig_honest_n.png` — given n vs ligand-eligible n
- `figures/fig_n_per_unit.png` — per-unit malignant / T/NK floors
- `figures/fig_family_barrier_bars.png` / `fig_family_barrier_forest.png`
- `figures/fig_family_ifn_bars.png` / `fig_family_ifn_forest.png`
- `figures/fig_family_barrier_units.png` / `fig_family_ifn_units.png`

## What is not claimed

- The n=34 ρ=−0.575 / Q4 r=−0.700 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071 is not added.
- Cell-pooled permutations are not the test. Patient/donor is the unit.
- This is not a CellChat discovery screen. Only the two pre-specified families.

## Reproduce

```bash
python3 methods/pair_123902_131907_lr_thesis_cldn4/scripts/download.py
python3 methods/pair_123902_131907_lr_thesis_cldn4/scripts/analyze.py
```

Hill constants: trim=0.1, Kh=0.5, expr_prop=0.1.
