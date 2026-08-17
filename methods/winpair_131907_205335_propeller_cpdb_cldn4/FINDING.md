# FINDING — propeller + CellPhoneDB-style on winning pair GSE131907+GSE205335

ADDITIVE **CLDN4 only**. No dual-high TACSTD2×CLDN4. **No GSE207422.**
Winning pair is the PR #320 author-malignant **%pos** cut vs T/NK
(Spearman n=43 ρ=−0.479; Q4 vs Q1 n=23 r=−0.705). Patient/sample is the unit.
LIANA is not run (that agent is separate). CellChat is not run.

## Honest n

| Cohort | Locked units | Unit | Q1 | Q4 | Compared | Unique patients in tails |
|---|---:|---|---:|---:|---:|---:|
| GSE131907 | 21 | sample (n_mal≥20; tumor origins) | 6 | 5 | 11 | 11 |
| GSE205335 | 22 | patient (≥20 mal + ≥20 T/NK) | 6 | 6 | 12 | 12 |
| Pair | 43 | mixed | 12 | 11 | 23 | 22 |

GSE131907 T/NK extract is **sample-level** (PR #320). The locked n=21
tumor samples with `n_malignant ≥ 20` are metastatic / advanced sites
(mLN, tL/B, mBrain) — no tLung sample meets that malignant-count floor.
Quartiles are cut **within cohort** on malignant CLDN4 %pos (`pd.qcut` on
average ranks). Q4 vs Q1 uses the tails only. Cells are not n.
GSE205335 B is B+plasma. GSE205335 Q4 mixes 3 SCLC with ADC/SQ (PR #362).

## 1. Propeller / speckle-style composition

Transform = empirical logit of the fraction (Phipson et al. 2022 *Bioinformatics*,
`speckle::propeller`). Test = Welch two-sample t on logit (cohort) or OLS
`logit ~ Q4 + cohort` (pair). BH on the pre-specified primary family:
**T/NK + B × 2 cohorts (4 tests)**. CD8 is nested in T/NK and is secondary.
eBayes moderation across 3–4 compartments is a companion (k is too small
for a stable prior; unmoderated logit p is the primary p).

### Primary family (cohort × T/NK, B)

| cohort | compartment | n_Q1/n_Q4 | median frac Q1 | Q4 | Δ | logit Δ | p | BH q | recover q<0.10 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| GSE131907 | TNK | 6/5 | 0.390 | 0.048 | -0.342 | -1.910 | 0.066 | 0.088 | yes |
| GSE131907 | B | 6/5 | 0.084 | 0.005 | -0.079 | -2.924 | 0.00676 | 0.027 | yes |
| GSE205335 | TNK | 6/6 | 0.514 | 0.116 | -0.398 | -1.682 | 0.027 | 0.053 | yes |
| GSE205335 | B | 6/6 | 0.063 | 0.034 | -0.028 | -0.395 | 0.505 | 0.505 | no |

### Winning-pair OLS (companion; BH across T/NK, B, CD8)

| compartment | n_Q1/n_Q4 | median frac Q1 | Q4 | Δ | logit β_Q4 | p | BH q |
|---|---:|---:|---:|---:|---:|---:|---:|
| TNK | 12/11 | 0.409 | 0.105 | -0.304 | -1.790 | 0.00116 | 0.00174 |
| B | 12/11 | 0.077 | 0.023 | -0.054 | -1.599 | 0.00934 | 0.00934 |
| CD8 | 12/11 | 0.142 | 0.015 | -0.128 | -2.100 | 0.000515 | 0.00154 |

**Propeller readout:** 3/4 primary tests have logit Δ < 0 and BH q < 0.10.
Fraction MWU (PR #320 companion) is in `tables/propeller_composition.tsv` and is not the primary p.

## 2. CellPhoneDB-style documented mean score (MHC-I + T-recruit)

Score = mean of partner means on log1p(CP10k); complexes = min of subunit means
(Efremova 2020 *Nat Protoc*; Garcia-Alonso 2022 *Nat Protoc*).
`expr_prop` = 0.10. Outgoing = author-malignant → **same-unit** T/NK.
Contrast = **between-unit CLDN4 Q4 vs Q1** (not the within-patient high/low LIANA split).
A pair enters the table if detected in ≥3 Q1 and ≥3 Q4 units. BH within each scope×cohort focus family.
This is **not** a CellChat probability and **not** LIANA `mt.cellphonedb`.

Winning-pair focus pairs that passed the detect gate: **13**. 11/13 have median Δ < 0 (weaker in Q4). **0/13 reach BH q < 0.05** in the weaker-in-Q4 direction.

T-recruit dropout is the n, not a hidden negative: CXCL9/10/11–CXCR3, CCL5–CCR5, CX3CL1–CX3CR1, CCL19/21–CCR7 and CXCL13–CXCR5 did not reach ≥3 Q1 and ≥3 Q4 `expr_prop` units. The only T-recruit pair in the table is CXCL16–CXCR6 (5/4). KIR / LILRB MHC-I pairs likewise failed the detect gate.

### Winning-pair focus table

| pathway | pair | n_Q1/n_Q4 | median Q1 | Q4 | Δ | r | p | BH q |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| MHC_I | HLA-E–KLRD1 | 10/8 | 0.892 | 0.684 | -0.208 | -0.475 | 0.101 | 0.680 |
| MHC_I | HLA-E–KLRC1 | 5/4 | 0.773 | 0.446 | -0.326 | -0.600 | 0.190 | 0.680 |
| MHC_I | HLA-E–KLRC1+KLRD1 | 5/3 | 0.773 | 0.419 | -0.353 | -0.600 | 0.250 | 0.680 |
| MHC_I | HLA-E–KLRK1 | 6/4 | 0.911 | 0.517 | -0.394 | -0.500 | 0.257 | 0.680 |
| MHC_I | HLA-B–CD8A | 12/7 | 1.614 | 1.376 | -0.238 | -0.333 | 0.261 | 0.680 |
| MHC_I | HLA-B–CD8B | 11/5 | 1.558 | 1.305 | -0.253 | -0.309 | 0.377 | 0.766 |
| MHC_I | HLA-E–KLRC2 | 3/3 | 0.853 | 0.583 | -0.270 | -0.333 | 0.700 | 0.965 |
| MHC_I | HLA-E–KLRC2+KLRD1 | 3/3 | 0.853 | 0.583 | -0.270 | -0.333 | 0.700 | 0.965 |
| MHC_I | HLA-A–CD8B | 11/5 | 1.428 | 1.404 | -0.023 | -0.127 | 0.743 | 0.965 |
| MHC_I | HLA-A–CD8A | 12/7 | 1.470 | 1.473 | +0.002 | -0.024 | 0.967 | 1 |
| MHC_I | HLA-C–CD8A | 12/7 | 1.333 | 1.349 | +0.016 | +0.024 | 0.967 | 1 |
| MHC_I | HLA-C–CD8B | 11/5 | 1.335 | 1.315 | -0.020 | +0.018 | 1 | 1 |
| T_recruit | CXCL16–CXCR6 | 5/4 | 0.423 | 0.308 | -0.115 | -0.400 | 0.413 | 0.766 |

Per-cohort rows: `tables/cpdb_focus.tsv`.

## What was not run

- LIANA `mt.cellphonedb` (separate agent).
- CellChat R.
- Dual-high TACSTD2×CLDN4.
- GSE207422.
- scCODA HMC / Milo / muscat.

## Files

| File | Role |
|---|---|
| `tables/propeller_composition.tsv` | **Table 1** — propeller / speckle-style |
| `tables/cpdb_focus.tsv` | **Table 2** — CellPhoneDB-style MHC-I + T-recruit |
| `tables/units_with_quartiles.tsv` | Locked units + Q labels |
| `figures/fig_propeller_fractions_q4q1.png` | Extra: T/NK and B boxes |
| `figures/fig_propeller_effects.png` | Extra: logit Δ forest |
| `figures/fig_honest_n.png` | Extra: honest n |
| `figures/fig_cpdb_mhci.png` | Extra: MHC-I Δ |
| `figures/fig_cpdb_trecruit.png` | Extra: T-recruit Δ |
| `figures/fig_cpdb_focus_pair.png` | Extra: focus pair bars |
| `figures/fig_extra_pair_boxes.png` | Extra: pair composition boxes |

## Reproduce

```bash
python3 methods/winpair_131907_205335_propeller_cpdb_cldn4/scripts/download.py
python3 methods/winpair_131907_205335_propeller_cpdb_cldn4/scripts/analyze.py
```

