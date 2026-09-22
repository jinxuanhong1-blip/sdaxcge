# FINDING — largest TACSTD2 (TROP2+) barrier-ligand Δ, concordant-4

ADDITIVE. TACSTD2 gate for PPT **TROP2+ barrier face** (Part1).
Does not replace the locked CLDN4 max-effect table (PR 713, +27.8 pp) or CellChat PR 616 (+0.0037).
No dual-high. No CLDN4 gate here. Concordant four only.
GSE148071, GSE127465, GSE154826, GSE200563, and E-MTAB-13526 are not added.
Senders are malignant TACSTD2-high vs TACSTD2-low. Receiver is T/NK.
The patient is the unit. Cell-pooled tests are not reported.

## Headline

Largest eligible family Δ is **expr_prop_pp** on the **decile** split. Per-ligand mean **+26.14** (percentage points of ligand-positive malignant cells). The family sum of F11R + NECTIN2 + CDH1 + LGALS9 is **+104.5**. That sum is not a single proportion.

n=59. Wilcoxon p=5.42e-11. Sign-flip one-sided p=1.00e-04. Patients with family sum > 0: 95%. Cohort means of the sum: +87.72 / +99.63 / +121.5 / +95.28 (n 11+18+21+9). Random-effects I²=0%.

Q4 vs Q1, same score (locked gate style): n=63, per-ligand mean +25.19, family sum +100.8, Wilcoxon p=7.23e-12, cohort means +75.59 / +87.93 / +128.7 / +100.5 (n 13+20+21+9). This Q4 figure matches the TACSTD2 crude row in PR 716.

**PPT one-liner (Part1 TROP2+ barrier face):** paste `results/PPT_SLIDE.md`. Headline number is **+26 pp**, not CellChat Hill +0.074 and not population-scaled +0.004.

LIANA / CellChat on the same decile cells (per-ligand means): CellPhoneDB lr_means **+0.126**, Connectome **+0.123**, CellChat Hill (no pop.) **+0.074**, Hill×pop **+0.0040**, LIANA log2FC **+0.364**. All 4/4 cohorts positive.

Selection tier: 4/4 ligands each 4/4 cohorts (max absolute family sum).
The menu and the gates were fixed in the script before the ranking. The winner is the largest family sum among scores that pass. Units are not interchangeable. A percentage-point sum across four ligands is the absolute gap in how many sender cells express each ligand, gated on a detected T/NK receptor. It is not a CellChat probability. Fold-changes and population-scaled Hill sums are in the menu below and are smaller numbers on their own scales.
Chemokine family is also 4/4 on this scale because CXCL16 alone is +16.4 pp; CXCL9/10/11/CCL5 are not barrier-sized. Barrier per-ligand mean (+26.14) remains larger than chemokine per-ligand mean (+3.47).

## Why CellChat population-scaled Δ stays tiny

Official CellChat with `population.size=TRUE` multiplies every edge by sender × receiver proportions, so a real ligand gap becomes a probability near zero (CLDN4 PR 616 family sum +0.0037). The same seven pairs on the TACSTD2 gate, recomputed here:

| score | split | n | mean pair-sum Δ | cohorts + | p Wilcoxon | p sign-flip |
|---|---|---:|---:|---:|---|---|
| cellchat_hill | decile | 59 | +0.402 | 4/4 | 2.15e-10 | 1.00e-04 |
| cellchat_hill_pop | decile | 59 | +0.022 | 4/4 | 2.28e-10 | 1.00e-04 |
| cellchat_prod | decile | 59 | +0.585 | 4/4 | 3.55e-10 | 1.00e-04 |
| cellchat_hill | q4q1 | 63 | +0.364 | 4/4 | 1.14e-10 | 1.00e-04 |
| cellchat_hill_pop | q4q1 | 63 | +0.027 | 4/4 | 3.45e-11 | 1.00e-04 |
| cellchat_prod | q4q1 | 63 | +0.525 | 4/4 | 2.92e-10 | 1.00e-04 |

## Barrier ligands

| ligand | n | mean Δ | median | frac > 0 | cohorts + | 123902 | 131907 | 205335 | 189357 | p Wilcoxon |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| F11R | 59 | +27.15 | +30.47 | 66% | 4/4 | +20.67 | +22.48 | +30.35 | +36.95 | 3.85e-08 |
| NECTIN2 | 59 | +30.49 | +27.93 | 92% | 4/4 | +28.22 | +29.66 | +35.29 | +23.69 | 9.89e-11 |
| CDH1 | 59 | +33.50 | +35.00 | 78% | 4/4 | +22.46 | +32.93 | +41.09 | +30.44 | 3.52e-09 |
| LGALS9 | 59 | +13.41 | +11.45 | 68% | 4/4 | +16.37 | +14.56 | +14.81 | +4.205 | 2.89e-07 |

## IFN / recruit, separate from HLA–CD8

Chemokine family (CXCL9, CXCL10, CXCL11, CCL5, CXCL16), same score and split: sum **+17.37**, per-ligand mean +3.474, n=59, Wilcoxon p=2.13e-06, cohorts high>low 4/4, cohort means +36.18 / +16.79 / +13.08 / +5.556.

HLA–CD8 (HLA-A, HLA-B, HLA-C) is not recruitment. Same score: sum **+41.24**, per-ligand mean +13.75, Wilcoxon p=7.74e-09, cohorts high>low 4/4.

| ligand | family | mean Δ | cohorts + | p Wilcoxon |
|---|---|---:|---:|---|
| CXCL9 | ifn_recruit | +0.202 | 2/4 | 0.465 |
| CXCL10 | ifn_recruit | +1.327 | 3/4 | 0.128 |
| CXCL11 | ifn_recruit | +0.363 | 1/4 | 0.317 |
| CCL5 | ifn_recruit | -0.894 | 0/4 | 0.068 |
| CXCL16 | ifn_recruit | +16.37 | 4/4 | 1.73e-06 |
| HLA-A | mhc_cd8 | +14.42 | 4/4 | 1.75e-08 |
| HLA-B | mhc_cd8 | +13.58 | 4/4 | 1.13e-08 |
| HLA-C | mhc_cd8 | +13.24 | 4/4 | 6.29e-08 |

## Full menu (barrier family sum, both splits)

Eligible scores had to be positive, Wilcoxon p<0.05, positive in all four cohorts (each n≥3), larger per ligand than the chemokine family, and, for untrimmed CP10k scores, at least half as large after 10% trimming. Each of the four ligands also had to be positive. The winning tier is named above.

| method | split | n | family sum | per ligand | frac>0 | cohorts + | ligands 4/4 | p | separates | trim ok |
|---|---|---:|---:|---:|---:|---:|---|---|---|---|
| expr_prop_pp | decile | 59 | +104.5 | +26.14 | 95% | 4/4 | yes | 5.42e-11 | yes | yes |
| conn_z | decile | 59 | +3.679 | +0.920 | 93% | 4/4 | no | 3.93e-10 | yes | yes |
| cp10k_delta_trim | decile | 59 | +2.161 | +0.540 | 93% | 4/4 | yes | 9.21e-11 | yes | yes |
| cp10k_delta | decile | 59 | +2.147 | +0.537 | 92% | 4/4 | no | 9.71e-11 | yes | yes |
| cp10k_log2fc_trim | decile | 59 | +1.864 | +0.466 | 95% | 4/4 | yes | 7.87e-11 | yes | yes |
| cp10k_log2fc | decile | 59 | +1.510 | +0.377 | 93% | 4/4 | no | 7.87e-11 | yes | yes |
| liana_logfc | decile | 59 | +1.456 | +0.364 | 95% | 4/4 | yes | 7.08e-11 | yes | yes |
| natmi_spec | decile | 59 | +1.031 | +0.258 | 95% | 4/4 | no | 9.21e-11 | yes | yes |
| cpdb_means | decile | 59 | +0.505 | +0.126 | 95% | 4/4 | yes | 7.08e-11 | yes | yes |
| conn_prod | decile | 59 | +0.492 | +0.123 | 93% | 4/4 | yes | 9.21e-11 | yes | yes |
| cellchat_prod | decile | 59 | +0.423 | +0.106 | 95% | 4/4 | yes | 7.87e-11 | yes | yes |
| sca_lrscore | decile | 59 | +0.362 | +0.091 | 95% | 4/4 | no | 9.21e-11 | yes | yes |
| cellchat_hill | decile | 59 | +0.297 | +0.074 | 95% | 4/4 | yes | 8.29e-11 | yes | yes |
| cellchat_hill_pop | decile | 59 | +0.016 | +0.003963 | 95% | 4/4 | yes | 6.71e-11 | yes | yes |
| expr_prop_pp | q4q1 | 63 | +100.8 | +25.19 | 95% | 4/4 | yes | 7.23e-12 | yes | yes |
| conn_z | q4q1 | 63 | +3.464 | +0.866 | 95% | 4/4 | no | 2.15e-11 | yes | yes |
| cp10k_delta_trim | q4q1 | 63 | +1.746 | +0.437 | 94% | 4/4 | yes | 3.76e-11 | yes | yes |
| cp10k_delta | q4q1 | 63 | +1.592 | +0.398 | 90% | 4/4 | no | 1.84e-10 | yes | yes |
| cp10k_log2fc_trim | q4q1 | 63 | +1.570 | +0.392 | 94% | 4/4 | yes | 2.36e-11 | yes | yes |
| liana_logfc | q4q1 | 63 | +1.233 | +0.308 | 94% | 4/4 | yes | 3.43e-11 | yes | yes |
| cp10k_log2fc | q4q1 | 63 | +1.148 | +0.287 | 90% | 4/4 | no | 1.18e-10 | yes | yes |
| natmi_spec | q4q1 | 63 | +0.908 | +0.227 | 95% | 4/4 | yes | 1.55e-11 | yes | yes |
| cpdb_means | q4q1 | 63 | +0.427 | +0.107 | 94% | 4/4 | yes | 3.43e-11 | yes | yes |
| conn_prod | q4q1 | 63 | +0.408 | +0.102 | 92% | 4/4 | yes | 9.83e-11 | yes | yes |
| cellchat_prod | q4q1 | 63 | +0.374 | +0.093 | 94% | 4/4 | yes | 3.43e-11 | yes | yes |
| sca_lrscore | q4q1 | 63 | +0.315 | +0.079 | 94% | 4/4 | no | 2.47e-11 | yes | yes |
| cellchat_hill | q4q1 | 63 | +0.265 | +0.066 | 94% | 4/4 | yes | 1.55e-11 | yes | yes |
| cellchat_hill_pop | q4q1 | 63 | +0.019 | +0.004862 | 94% | 4/4 | yes | 1.70e-11 | yes | yes |

## Honest n

Inventory units: 65. Q4 vs Q1 requires ≥40 malignant cells, ≥10 cells in each arm, and ≥20 T/NK cells. Decile uses the outer 10% and the same arm floor, so patients with fewer than 100 malignant cells drop out.
Units with n_mal≥40 in the inventory built here: 64 (reference expectation 64; P4001 has 27 malignant cells).
EBUS_13 has flat malignant TACSTD2 (all zeros) and is skipped, so TACSTD2 Q4 n=63 (CLDN4 Q4 was 64). Decile n=59.

## What is not claimed

- This is an expression ligand–receptor contrast. It is not spatial exclusion and it does not replace CosMx TACSTD2 short-range cold (PR 726/735) or the locked CLDN4 CosMx 8/8.
- Population-scaled CellChat probabilities are not re-issued as +0.0037. The pair-sum table above is the recomputed Hill score on the TACSTD2 gate.
- A larger number on a CP10k or percentage-point scale is a different unit from a probability. The menu shows both.
- CLDN4 is not the gate here. No dual-high. No mediation claim. GSE148071 is not added.
- Part1 slide claim is TROP2+ malignant cells carry a stronger barrier-ligand face toward T/NK; immune exclusion ρ stays on CLDN4 %pos (PR 539/712).

## Reproduce

```bash
bash methods/concordant4_ccc_maxeffect_tacstd2/scripts/download.sh /tmp/concordant4_raw
python3 methods/concordant4_ccc_maxeffect_tacstd2/scripts/prepare_131907.py
gcc -O3 -o /tmp/extract_131907 methods/concordant4_ccc_maxeffect_tacstd2/scripts/extract_131907.c -lz
/tmp/extract_131907 \
  /tmp/concordant4_raw/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz \
  /tmp/concordant4_cache/gse131907/genes_request.txt \
  /tmp/concordant4_cache/gse131907/keep_idx.i32 \
  208506 \
  /tmp/concordant4_cache/gse131907
python3 methods/concordant4_ccc_maxeffect_tacstd2/scripts/run_maxeffect.py
```

