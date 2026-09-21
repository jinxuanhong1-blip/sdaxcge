# FINDING — largest barrier-ligand communication Δ, concordant-4

ADDITIVE. Does not replace the locked CellChat probability table (PR 616, family sum +0.0037).
CLDN4 only. No dual-high. No TACSTD2 gate. Concordant four only.
GSE148071, GSE127465, GSE154826, GSE200563, and E-MTAB-13526 are not added.
Senders are malignant CLDN4-high vs CLDN4-low. Receiver is T/NK.
The patient is the unit. Cell-pooled tests are not reported.

## Headline

On the decile split, CLDN4-high malignant cells are **+27.8 percentage points** more often positive for a barrier ligand than CLDN4-low cells from the same patient (mean of F11R, NECTIN2, CDH1, LGALS9). The family sum of those four gaps is **+111.1** percentage points. That sum is not a single proportion.

n=60. Wilcoxon p=3.69e-11. Sign-flip one-sided p=1.00e-04 (none of 10,000 flips reached the observed mean). Patients with a positive family sum: 95%. Cohort means of the sum: +85.2 / +112.3 / +132.4 / +90.2 (n 11+19+21+9). Random-effects I²=39%. All four ligands are positive in all four cohorts.

The locked Q4 vs Q1 split is almost the same size and keeps the patient who drops out of the decile: n=64, per-ligand mean **+25.6** percentage points, family sum **+102.4**, Wilcoxon p=4.91e-12, cohort means +85.3 / +87.1 / +129.7 / +99.1 (n 13+21+21+9).

Per ligand, decile: F11R +29.6, NECTIN2 +33.1, CDH1 +33.8, LGALS9 +14.6. Q4 vs Q1: F11R +27.4, NECTIN2 +30.5, CDH1 +32.3, LGALS9 +12.3.

This is a CellPhoneDB expression-proportion contrast. An edge counts only when the T/NK receptor is present in at least 10% of T/NK cells (NECTIN2–CD96/TIGIT, F11R–LFA-1, CDH1–integrin/KLRG1, LGALS9–CD44/CD45/TIM3). The score is how many more sender cells carry the ligand, not a CellChat probability.

The same cells on a fold-change scale are smaller: 10% trimmed CP10k log2FC averages **+0.45 per ligand** on the decile (about 1.4-fold; family sum +1.82). Trimmed CP10k abundance differs by **+0.51 CP10k per ligand**. The large absolute effect is the gain in ligand-positive cells, not a large jump in counts per positive cell.

Population-scaled Hill probabilities stay small. The seven-pair sum (the PR 616 pair list) is +0.026 on Q4 vs Q1 and +0.021 on the decile. Without the population-size weight that seven-pair sum is +0.35 (Q4) and +0.41 (decile). Official CellChat in PR 616 was +0.0037. The recomputed population-scaled sum is the same kind of object and is still far below the percentage-point result. It is not forced to equal +0.0037.

Selection tier: each of the four ligands positive in 4/4 cohorts. The stricter tier (chemokine family not positive in every cohort) did not pass, because CXCL16 lifts the chemokine family sum. CXCL9, CXCL10, CXCL11, and CCL5 do not. Barrier per-ligand mean (+27.8) is larger than the chemokine per-ligand mean (+2.9).

The menu and the gates were fixed in the script before the ranking. The winner is the largest family sum among scores that pass. Units are not interchangeable.

## Why +0.0037 was small

PR 616 summed official CellChat probabilities with `population.size=TRUE`. That multiplies every edge by the sender proportion times the receiver proportion, so a real ligand gap becomes a probability near zero. The same seven pairs, recomputed here:

| score | split | n | mean pair-sum Δ | cohorts + | p Wilcoxon | p sign-flip |
|---|---|---:|---:|---:|---|---|
| cellchat_hill | decile | 60 | +0.409 | 4/4 | 9.89e-11 | 1.00e-04 |
| cellchat_hill_pop | decile | 60 | +0.021 | 4/4 | 2.45e-10 | 1.00e-04 |
| cellchat_prod | decile | 60 | +0.591 | 4/4 | 1.29e-10 | 1.00e-04 |
| cellchat_hill | q4q1 | 64 | +0.354 | 4/4 | 9.95e-11 | 1.00e-04 |
| cellchat_hill_pop | q4q1 | 64 | +0.026 | 4/4 | 2.42e-11 | 1.00e-04 |
| cellchat_prod | q4q1 | 64 | +0.510 | 4/4 | 2.08e-10 | 1.00e-04 |

## Barrier ligands

| ligand | n | mean Δ | median | frac > 0 | cohorts + | 123902 | 131907 | 205335 | 189357 | p Wilcoxon |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| F11R | 60 | +29.64 | +34.09 | 70% | 4/4 | +19.39 | +27.33 | +34.87 | +34.88 | 1.65e-08 |
| NECTIN2 | 60 | +33.06 | +32.99 | 90% | 4/4 | +31.79 | +32.24 | +39.90 | +20.39 | 1.38e-10 |
| CDH1 | 60 | +33.77 | +41.39 | 77% | 4/4 | +24.84 | +34.21 | +39.13 | +31.25 | 3.52e-09 |
| LGALS9 | 60 | +14.58 | +10.11 | 63% | 4/4 | +9.174 | +18.53 | +18.50 | +3.663 | 1.82e-07 |

## IFN / recruit, separate from HLA–CD8

Chemokine family (CXCL9, CXCL10, CXCL11, CCL5, CXCL16), same score and split: sum **+14.73**, per-ligand mean +2.946, n=60, Wilcoxon p=2.74e-05, cohorts high>low 4/4, cohort means +27.71 / +12.90 / +13.93 / +4.587.

HLA–CD8 (HLA-A, HLA-B, HLA-C) is not recruitment. Same score: sum **+49.02**, per-ligand mean +16.34, Wilcoxon p=9.92e-10, cohorts high>low 4/4.

| ligand | family | mean Δ | cohorts + | p Wilcoxon |
|---|---|---:|---:|---|
| CXCL9 | ifn_recruit | +0.046 | 1/4 | 1.000 |
| CXCL10 | ifn_recruit | +0.585 | 3/4 | 0.263 |
| CXCL11 | ifn_recruit | -0.079 | 1/4 | 0.655 |
| CCL5 | ifn_recruit | -0.404 | 1/4 | 0.345 |
| CXCL16 | ifn_recruit | +14.58 | 4/4 | 4.33e-06 |
| HLA-A | mhc_cd8 | +16.69 | 4/4 | 1.70e-09 |
| HLA-B | mhc_cd8 | +16.14 | 4/4 | 8.75e-10 |
| HLA-C | mhc_cd8 | +16.19 | 4/4 | 2.80e-09 |

## Full menu (barrier family sum, both splits)

Eligible scores had to be positive, Wilcoxon p<0.05, positive in all four cohorts (each n≥3), larger per ligand than the chemokine family, and, for untrimmed CP10k scores, at least half as large after 10% trimming. Each of the four ligands also had to be positive. The winning tier is named above.

| method | split | n | family sum | per ligand | frac>0 | cohorts + | ligands 4/4 | p | separates | trim ok |
|---|---|---:|---:|---:|---:|---:|---|---|---|---|
| expr_prop_pp | decile | 60 | +111.1 | +27.76 | 95% | 4/4 | yes | 3.69e-11 | yes | yes |
| conn_z | decile | 60 | +3.930 | +0.983 | 93% | 4/4 | no | 4.80e-11 | yes | yes |
| cp10k_delta_trim | decile | 60 | +2.047 | +0.512 | 97% | 4/4 | yes | 3.51e-11 | yes | yes |
| cp10k_delta | decile | 60 | +2.013 | +0.503 | 95% | 4/4 | no | 7.26e-11 | yes | yes |
| cp10k_log2fc_trim | decile | 60 | +1.819 | +0.455 | 97% | 4/4 | yes | 3.51e-11 | yes | yes |
| cp10k_log2fc | decile | 60 | +1.478 | +0.370 | 95% | 4/4 | no | 5.32e-11 | yes | yes |
| liana_logfc | decile | 60 | +1.468 | +0.367 | 97% | 4/4 | yes | 3.51e-11 | yes | yes |
| natmi_spec | decile | 60 | +1.111 | +0.278 | 97% | 4/4 | no | 3.51e-11 | yes | yes |
| cpdb_means | decile | 60 | +0.509 | +0.127 | 97% | 4/4 | yes | 3.51e-11 | yes | yes |
| conn_prod | decile | 60 | +0.493 | +0.123 | 97% | 4/4 | yes | 3.51e-11 | yes | yes |
| cellchat_prod | decile | 60 | +0.421 | +0.105 | 97% | 4/4 | yes | 3.51e-11 | yes | yes |
| sca_lrscore | decile | 60 | +0.401 | +0.100 | 97% | 4/4 | no | 3.51e-11 | yes | yes |
| cellchat_hill | decile | 60 | +0.300 | +0.075 | 97% | 4/4 | yes | 3.51e-11 | yes | yes |
| cellchat_hill_pop | decile | 60 | +0.015 | +0.003769 | 97% | 4/4 | yes | 3.51e-11 | yes | yes |
| expr_prop_pp | q4q1 | 64 | +102.4 | +25.60 | 97% | 4/4 | yes | 4.91e-12 | yes | yes |
| conn_z | q4q1 | 64 | +3.450 | +0.863 | 92% | 4/4 | yes | 1.57e-11 | yes | yes |
| cp10k_delta_trim | q4q1 | 64 | +1.697 | +0.424 | 97% | 4/4 | yes | 4.26e-12 | yes | yes |
| cp10k_log2fc_trim | q4q1 | 64 | +1.541 | +0.385 | 97% | 4/4 | yes | 4.26e-12 | yes | yes |
| cp10k_delta | q4q1 | 64 | +1.479 | +0.370 | 94% | 4/4 | no | 2.27e-11 | yes | yes |
| liana_logfc | q4q1 | 64 | +1.215 | +0.304 | 97% | 4/4 | no | 6.21e-12 | yes | yes |
| cp10k_log2fc | q4q1 | 64 | +1.086 | +0.271 | 94% | 4/4 | no | 2.85e-11 | yes | yes |
| natmi_spec | q4q1 | 64 | +0.891 | +0.223 | 95% | 4/4 | no | 9.02e-12 | yes | yes |
| cpdb_means | q4q1 | 64 | +0.421 | +0.105 | 97% | 4/4 | no | 6.21e-12 | yes | yes |
| conn_prod | q4q1 | 64 | +0.399 | +0.100 | 94% | 4/4 | no | 2.27e-11 | yes | yes |
| cellchat_prod | q4q1 | 64 | +0.366 | +0.091 | 94% | 4/4 | yes | 1.31e-11 | yes | yes |
| sca_lrscore | q4q1 | 64 | +0.300 | +0.075 | 94% | 4/4 | no | 1.50e-11 | yes | yes |
| cellchat_hill | q4q1 | 64 | +0.261 | +0.065 | 95% | 4/4 | yes | 6.21e-12 | yes | yes |
| cellchat_hill_pop | q4q1 | 64 | +0.019 | +0.004688 | 95% | 4/4 | yes | 5.39e-12 | yes | yes |

## Honest n

Inventory units: 65. Q4 vs Q1 requires ≥40 malignant cells, ≥10 cells in each arm, and ≥20 T/NK cells. Decile uses the outer 10% and the same arm floor, so patients with fewer than 100 malignant cells drop out.
Units with n_mal≥40 in the inventory built here: 64 (reference expectation 64; P4001 has 27 malignant cells).

## What is not claimed

- This is an expression ligand–receptor contrast. It is not spatial exclusion and it does not replace the CosMx 8/8 result.
- Population-scaled CellChat probabilities are not re-issued as +0.0037. The pair-sum table above is the recomputed Hill score.
- A larger number on a CP10k or percentage-point scale is a different unit from a probability. The menu shows both.
- TACSTD2 is not a gate. GSE148071 is not added.

## Reproduce

```bash
bash methods/concordant4_ccc_maxeffect_cldn4/scripts/download.sh /tmp/concordant4_raw
python3 methods/concordant4_ccc_maxeffect_cldn4/scripts/prepare_131907.py
gcc -O3 -o /tmp/extract_131907 methods/concordant4_ccc_maxeffect_cldn4/scripts/extract_131907.c -lz
/tmp/extract_131907 \
  /tmp/concordant4_raw/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz \
  /tmp/concordant4_cache/gse131907/genes_request.txt \
  /tmp/concordant4_cache/gse131907/keep_idx.i32 \
  208506 \
  /tmp/concordant4_cache/gse131907
python3 methods/concordant4_ccc_maxeffect_cldn4/scripts/run_maxeffect.py
```

