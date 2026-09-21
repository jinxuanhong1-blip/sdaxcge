# FINDING — CellChat v2 full network and sender roles, concordant-four

ADDITIVE. **Thesis already correct: CLDN4-high malignant cells are the barrier senders.**
This does not replace the locked 14-pair CellChat table (PR 540).
CLDN4 only. No dual-high. No TACSTD2 gate.
Concordant four only: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Do **not** add GSE148071 / GSE127465 / CD45-only.

Engine: R + Seurat 5.5.1 + CellChat 2.2.0.9001.
Python streams the GSE131907 text matrix only.
Senders = malignant CLDN4 Q4 vs Q1. Receivers in the extracted readout = T/NK.
Honest n = patient / locked sample. Cell-pooled tests are not reported.

## What was added

- Full CellChatDB.human network (secreted, ECM-receptor, cell-cell contact), not only 14 pairs.
- Official overexpressed L-R filter (`thresh.p = 0.05`) for that full network.
- Signaling roles on each pathway and on the two thesis family subnetworks:
  sender = weighted out-degree, receiver = in-degree, mediator = flow betweenness, influencer = information centrality
  (`netAnalysis_computeCentrality` / `computeCentralityLocal`).
- Three permutation procedures on the extracted barrier/inhibitory and IFN/recruit outgoing scores.

Thesis readout:

- Barrier/inhibitory outgoing **UP from CLDN4-high** is ON-thesis.
- IFN/recruit outgoing UP from CLDN4-low is the KD-like arm.
  Report it even when the sign is opposite. Do not bury barrier-up-in-high as a recruit-up skip.

## Honest n

| gate | n | note |
|---|---:|---|
| Locked four | 65 | 13+21+22+9 |
| Inventory units | 65 | GSE123902=13, GSE131907=21, GSE189357=9, GSE205335=22 |
| Both compartments in inventory (n_mal>=10, n_tnk>=20) | 65 | |
| Q4 vs Q1 CellChat units (primary) | **64** | status ok |
| Cells used after cap 200 / group (median high / low / TNK) | | 160 / 160 / 200 |
| Full-network L-R interactions per unit (median) | 642.5 | overexpressed filter |

GSE123902 / GSE189357: epithelium marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0), then T/NK markers.
GSE131907 / GSE205335: author malignant, then author T/NK. Normal-tissue biopsies in GSE205335 are excluded.
CLDN4 rank is log1p(CP10k) using the **full** UMI library size, then Q4 vs Q1.
CellChat input is the same LogNormalize (scale 1e4) with that full library size.
Each group is capped at 200 cells (patient-specific seed) before `computeCommunProb`.

## Primary — extracted outgoing to T/NK

CellChat `computeCommunProb`, truncatedMean trim=0.1, population.size=TRUE, nboot=100 on the pre-specified pairs.
`delta_sum` = (sum of high→T/NK probabilities) − (sum of low→T/NK) over the family, undetected pairs contribute 0.
`mean_detected_delta` keeps only pairs with probability > 0 on at least one arm (closer to the locked ligand table).
`delta_outdeg` is the sender-role contrast on the family subnetwork (all targets, not only T/NK).

Across-unit permutation: sign-flip of the patient deltas, B=10000, seed=3979. One-sided p is in the thesis direction.
Within-unit permutation: shuffle CLDN4-high vs low labels among the malignant cells, T/NK fixed, B=49, then recompute with CellChat. Stouffer combines those one-sided p-values. `frac_units_label_p05` is in the family table.

| family | metric | expect | n | 123902 | 131907 | 205335 | 189357 | mean Δ | p_signflip thesis | p_signflip two | p_W | observed | agrees |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|---|
| barrier_inhibitory | delta_sum | high>low | 64 | 13 | 21 | 21 | 9 | +0.0037 | 1.00e-04 | 1.00e-04 | 8.92e-10 | high>low | yes |
| barrier_inhibitory | mean_detected_delta | high>low | 62 | 13 | 20 | 20 | 9 | +0.0007 | 1.00e-04 | 1.00e-04 | 4.19e-10 | high>low | yes |
| barrier_inhibitory | delta_outdeg | high>low | 64 | 13 | 21 | 21 | 9 | +0.0040 | 1.00e-04 | 1.00e-04 | 1.21e-09 | high>low | yes |
| ifn_recruit | delta_sum | low>high | 64 | 13 | 21 | 21 | 9 | +0.0024 | 1 | 1.00e-04 | 5.50e-06 | high>low | opposite |
| ifn_recruit | mean_detected_delta | low>high | 58 | 11 | 19 | 19 | 9 | +0.0008 | 1 | 1.00e-04 | 6.37e-06 | high>low | opposite |
| ifn_recruit | delta_outdeg | low>high | 64 | 13 | 21 | 21 | 9 | +0.0024 | 1 | 1.00e-04 | 6.14e-06 | high>low | opposite |

### Alignment

- Barrier outgoing (primary `delta_sum`): barrier_inhibitory delta_sum: n=64 mean Δ=+0.0037 observed=high>low agrees=yes signflip_thesis=1.00e-04 wilcox=8.92e-10 stouffer_label=3.19e-24.
- Barrier sender role (`delta_outdeg`): barrier_inhibitory delta_outdeg: n=64 mean Δ=+0.0040 observed=high>low agrees=yes signflip_thesis=1.00e-04 wilcox=1.21e-09 stouffer_label=4.65e-24.
- IFN/recruit outgoing (primary `delta_sum`): ifn_recruit delta_sum: n=64 mean Δ=+0.0024 observed=high>low agrees=opposite signflip_thesis=1 wilcox=5.50e-06 stouffer_label=1.
- IFN/recruit sender role (`delta_outdeg`): ifn_recruit delta_outdeg: n=64 mean Δ=+0.0024 observed=high>low agrees=opposite signflip_thesis=1 wilcox=6.14e-06 stouffer_label=1.

Direction `agrees=yes` means the mean patient Δ has the thesis sign.
A sign-flip p below 0.05 is the across-patient permutation support for that sign.
The IFN/recruit arm is reported on its own line even if it is opposite.

## Pre-specified pairs (detected units only)

Sign-flip is across units with a detected pair. This is not a discovery screen.

| pair | family | expect | n | mean Δ | p_signflip thesis | observed | agrees |
|---|---|---|---:|---:|---|---|---|
| JAM1_ITGAL_ITGB2 | barrier_inhibitory | high>low | 41 | +0.0016 | 1.00e-04 | high>low | yes |
| NECTIN2_TIGIT | barrier_inhibitory | high>low | 54 | +0.0008 | 1.00e-04 | high>low | yes |
| CDH1_ITGAE_ITGB7 | barrier_inhibitory | high>low | 41 | +0.0006 | 1.00e-04 | high>low | yes |
| CDH1_KLRG1 | barrier_inhibitory | high>low | 30 | +0.0005 | 1.00e-04 | high>low | yes |
| LGALS9_HAVCR2 | barrier_inhibitory | high>low | 28 | +0.0001 | 0.134 | high>low | yes |
| LGALS9_CD44 | barrier_inhibitory | high>low | 52 | +0.0006 | 1.00e-04 | high>low | yes |
| LGALS9_CD45 | barrier_inhibitory | high>low | 51 | +0.0010 | 5.00e-04 | high>low | yes |
| CXCL9_CXCR3 | ifn_recruit | low>high | 1 | +0.0000 | NA | high>low | opposite |
| CXCL10_CXCR3 | ifn_recruit | low>high | 10 | -0.0001 | 0.424 | low>high | yes |
| CCL5_CCR5 | ifn_recruit | low>high | 5 | -0.0000 | 0.244 | low>high | yes |
| CCL5_CCR1 | ifn_recruit | low>high | 4 | +0.0000 | 0.684 | high>low | opposite |
| HLA-A_CD8A | ifn_recruit | low>high | 58 | +0.0009 | 1 | high>low | opposite |
| HLA-B_CD8A | ifn_recruit | low>high | 58 | +0.0009 | 1 | high>low | opposite |
| HLA-C_CD8A | ifn_recruit | low>high | 58 | +0.0009 | 1 | high>low | opposite |

## Control — full-network outgoing (not a thesis family)

Sum of every overexpressed L-R probability from the sender to T/NK, and the sender out-degree on that whole network.
This asks whether CLDN4-high is a stronger sender in general, not only on the barrier pairs.

- `delta_sum`: n=64 mean Δ=+0.0279 observed=high>low signflip_two=1.00e-04 signflip_high>low=1.00e-04 wilcox=5.03e-12.
- `delta_outdeg`: n=64 mean Δ=+0.0840 observed=high>low signflip_two=1.00e-04 signflip_high>low=1.00e-04 wilcox=3.61e-12.

## Full-network pathway roles

Pathways with a finite sender-role contrast: 165.
BH q-values across pathways are exploratory. They are not the thesis test.
Pre-specified pathway tags (the pathway that contains a thesis pair; the pathway may also contain other interactions):

| pathway | tag | n | mean Δ outdeg | mean Δ to T/NK | p_signflip two | q_BH | observed |
|---|---|---:|---:|---:|---|---|---|
| MHC-I | ifn_recruit | 61 | +0.0064 | +0.0064 | 1.00e-04 | 3.55e-04 | high>low |
| JAM | barrier_inhibitory | 63 | +0.0019 | +0.0008 | 1.00e-04 | 3.55e-04 | high>low |
| CDH1 | barrier_inhibitory | 56 | +0.0011 | +0.0005 | 1.00e-04 | 3.55e-04 | high>low |
| NECTIN | barrier_inhibitory | 58 | +0.0007 | +0.0006 | 1.00e-04 | 3.55e-04 | high>low |
| GALECTIN | barrier_inhibitory | 45 | +0.0030 | +0.0017 | 4.00e-04 | 0.00119 | high>low |
| CXCL | ifn_recruit | 35 | +0.0001 | +0.0002 | 0.0073 | 0.0165 | high>low |
| CCL | ifn_recruit | 23 | +0.0000 | +0.0000 | 0.451 | 0.605 | high>low |

- Barrier-tagged pathways with mean out-degree high>low: 4 / 4.
- IFN/recruit-tagged pathways with mean out-degree low>high: 0 / 3.

## What is not claimed

- TACSTD2 is not a gate. This is not dual-high.
- GSE148071, GSE127465, and CD45-only libraries are not added.
- The pathway table is not a discovery screen and not a 7-pool.
- Cell-pooled tests are not reported.
- With three nodes, mediator and influencer scores are secondary. The sender role is out-degree.
- PR 540 used a ligand-subset library size and `nboot = 1`. This run uses the full UMI library size and real permutations. Numbers are not expected to match that table digit for digit.

## Reproduce

```bash
Rscript methods/cellchat_v2_network_roles_concordant4/scripts/install_packages.R
bash methods/cellchat_v2_network_roles_concordant4/scripts/download.sh /tmp/concordant4_raw
Rscript methods/cellchat_v2_network_roles_concordant4/scripts/run_cellchat_v2.R --raw=/tmp/concordant4_raw
```

Parameters: max_cells=200, nboot_full=20, nboot_pair=100, label_perm=49, signflip=10000, version=v2.
Q4 vs Q1 requires n_mal >= 40, each arm >= 10, T/NK >= 20.

