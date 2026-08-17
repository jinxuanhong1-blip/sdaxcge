# FINDING — triple GSE131907 + GSE148071 + GSE205335, CLDN4-only CellChat

ADDITIVE. **CLDN4 only.** No dual-high. No GSE207422. The 131907+205335-only CellChat is a different agent and is not re-run. PR #320 pair Q4 *r*=−0.705 and GSE148071 n=25 partial ρ=−0.49 are **given** and are not re-audited. Patient is the unit (GSE131907 is sample-level).

p-values are descriptive.

## Honest n (triple patient table first)

| cohort | unit | malignant | floor | n | CellChat floor pass |
|---|---|---|---|---:|---:|
| GSE131907 | sample | author Malignant cells | PR #320 extract | **21** | 21 |
| GSE148071 | patient | putative epithelium | ≥25 / ≥25 (given partial) | **25** | 25 |
| GSE205335 | patient | author malignant | ≥20 / ≥20 | **22** | 22 |
| **triple** | mixed | — | — | **68** | **68** |

Do not write n=44 (GSE131907 series) or n=42 (GSE148071 deposited) or n=26 (GSE205335 GEO patients). The computable units are above.

## Combo rho table (locked singles → new triple)

Primary family is malignant CLDN4 vs same-unit T/NK. Spearman pool is DerSimonian–Laird on Fisher-z of the **given** singles (not re-audited). The PR #320 pair Q4 row is listed as given.

| analysis | combo | k | N | effect (p) | source |
|---|---|---:|---:|---|---|
| spearman_single_given | GSE131907 | 1 | 21 · Q1/Q4=6/5 | -0.522 (0.0152) | PR #320; not re-audited |
| q4q1_single_given | GSE131907 | 1 | 11 · Q1/Q4=6/5 | -0.600 (0.126) | PR #320; not re-audited |
| spearman_single_given | GSE205335 | 1 | 22 · Q1/Q4=6/6 | -0.435 (0.0429) | PR #320; not re-audited |
| q4q1_single_given | GSE205335 | 1 | 12 · Q1/Q4=6/6 | -0.778 (0.026) | PR #320; not re-audited |
| spearman_single_given | GSE148071 | 1 | 25 | -0.490 (0.0129) | given n=25 partial; not re-audited |
| q4q1_pair_given | GSE131907+GSE205335 | 2 | 23 · Q1/Q4=12/11 | -0.705 (3.01e-04, I²=0%) | PR #320 author/tnk/pct; not re-audited |
| spearman_triple_combo | GSE131907+GSE148071+GSE205335 | 3 | 68 | -0.483 (5.23e-05, I²=0%) | NEW triple Fisher-z of locked singles |
| spearman_pair_given | GSE131907+GSE205335 | 2 | 43 | -0.479 (0.00153, I²=0%) | PR #320 author/tnk/pct reconstructed from locked singles (not a re-audit of Q4 r=−0.705) |
| q4q1_148071_table_extra | GSE148071 | 1 | 13 · Q1/Q4=7/6 | +0.000 (1) | NEW Q4 on the n=25 eligibility table |

**New triple Spearman** (k=3, N=68): ρ=-0.483 p=5.23e-05 I²=0%. Locked pair Spearman GSE131907+GSE205335 is ρ=-0.479 (PR #320 Q4 *r*=-0.705 is given and not re-ranked).

GSE148071 n=25 partial ρ=−0.49 is the locked single. The eligibility table’s T/NK-fraction Spearman is a different cut and is not used to replace that given ρ.

### Locked singles (not re-audited)

| cohort | n | unit | ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) |
|---|---:|---|---|---|
| GSE131907 | 21 | sample | -0.522 (0.0152) | -0.600 (0.126; 6/5) |
| GSE205335 | 22 | patient | -0.435 (0.0429) | -0.778 (0.026; 6/6) |
| GSE148071 | 25 | patient | -0.490 (0.0129) | — |

Full combo table: [`results/combo_rho_table.tsv`](results/combo_rho_table.tsv). Patient table: [`results/triple_patient_table.tsv`](results/triple_patient_table.tsv).

## Ligand table (outgoing CLDN4-high malignant → same-patient T/NK)

Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs. Outgoing = malignant → **same-patient** T/NK. Within-cohort CLDN4 %pos quartiles (scales are not mixed). Test = Mann–Whitney on per-patient *P* (detected arm n≥3). CellChat R was not run. Cell-pooled means are not the test.

Detect-gated outgoing rows: 75. p<0.05: **3**. Only these three are claimed. Q4 vs Q1 stacks **within-cohort** tails (do not read as a global 68-patient quartile).

| pair | class | n_Q1/n_Q4 | cohorts | median P Q1 | median P Q4 | Δ | r | p |
|---|---|---|---:|---:|---:|---:|---:|---|
| MDK–NCL ** | other | 19/17 | 3 | 0.517 | 0.724 | +0.207 | +0.659 | 7.83e-04 |
| APP–CD74 ** | other | 18/17 | 3 | 0.482 | 0.618 | +0.136 | +0.431 | 0.0306 |
| MDK–ITGA4_ITGB1 ** | other | 16/9 | 3 | 0.211 | 0.412 | +0.200 | +0.528 | 0.0338 |
| ICAM1–ITGAL | other | 14/4 | 3 | 0.040 | 0.182 | +0.142 | +0.643 | 0.0614 |
| F11R–ITGAL_ITGB2 | barrier | 14/6 | 3 | 0.135 | 0.312 | +0.177 | +0.548 | 0.0622 |
| SPP1–CD44 | other | 12/6 | 3 | 0.508 | 0.056 | -0.451 | -0.556 | 0.0668 |
| CLEC2B–KLRB1 | other | 15/8 | 3 | 0.039 | 0.076 | +0.037 | +0.467 | 0.0755 |
| ICAM1–ITGAL_ITGB2 | other | 13/4 | 3 | 0.097 | 0.240 | +0.142 | +0.615 | 0.079 |
| ICAM1–SPN | other | 14/9 | 3 | 0.079 | 0.141 | +0.062 | +0.444 | 0.0832 |
| HLA-E–CD94:NKG2C | inhibitory | 3/3 | 2 | 0.521 | 0.134 | -0.386 | -1.000 | 0.1 |
| HLA-E–KLRC2 | inhibitory | 3/3 | 2 | 0.343 | 0.036 | -0.307 | -1.000 | 0.1 |
| JAG1–NOTCH1 | other | 3/3 | 1 | 0.001 | 0.016 | +0.015 | +1.000 | 0.1 |
| MIF–CD74_CD44 | other | 19/16 | 3 | 0.841 | 0.651 | -0.190 | -0.309 | 0.124 |
| LAMA3–CD44 | other | 11/4 | 3 | 0.059 | 0.212 | +0.153 | +0.545 | 0.138 |
| HLA-E–CD94:NKG2A | inhibitory | 6/3 | 3 | 0.309 | 0.038 | -0.271 | -0.667 | 0.167 |
| BAG6–NCR3 | other | 6/3 | 2 | 0.019 | 0.041 | +0.022 | +0.667 | 0.167 |
| HLA-E–KLRK1 | inhibitory | 6/4 | 1 | 0.406 | 0.168 | -0.238 | -0.583 | 0.171 |
| CDH1–ITGAE_ITGB7 | barrier|inhibitory | 10/8 | 3 | 0.065 | 0.130 | +0.065 | +0.400 | 0.173 |

**MDK–NCL** is the lead outgoing pair (higher in CLDN4 Q4; n=19 vs 17; all three cohorts). That is the same pair as the GSE205335-only extra (PR #362); here it is the triple same-patient test, not a cell-pooled redo of 131907+205335. APP–CD74 and MDK–ITGA4/ITGB1 are the only other p<0.05 rows. F11R–ITGAL/ITGB2 (barrier) is p=0.062 and is **not** claimed. SPP1–CD44 is lower in Q4 (p=0.067), not claimed.

CXCL16–CXCR6 is slightly **higher** in Q4 (ΔP +0.016, p=0.91) — not a recruit-down story. NECTIN2–TIGIT is detected (11/9) but p=0.49. **CD274–PDCD1** is detected in 8 units only and does not pass the ≥3 / ≥3 detect gate. CXCL9–CXCR3 is not detect-gated.

GSE205335 Q4 remains SCLC-mixed (given). GSE131907 NS_12 is at the T/NK floor (22 cells) and is kept.

Full LR table: [`results/ligand_table.tsv`](results/ligand_table.tsv) and [`results/lr_q4q1_all.tsv`](results/lr_q4q1_all.tsv).

## What was not done

- No dual-high TACSTD2×CLDN4 score.
- No GSE207422-only CellChat.
- No 131907+205335-only CellChat redo.
- GSE131907 tS1–tS3 / tLung-only CellChat is the other agent’s slice.
- CellChat R and LIANA were not run.

## Files

- `results/combo_rho_table.tsv` — locked singles + new triple pool
- `results/triple_patient_table.tsv` — honest n / quartiles
- `results/ligand_table.tsv` — CellChat-style outgoing table
- `figures/scatter_triple_cldn4_tnk.png`
- `figures/forest_triple_spearman.png`
- `figures/q4q1_tnk_boxes.png`
- `figures/fig_extra_ligand_table.png`
- `figures/n_cells_per_unit.png`

Reproduce: `python3 methods/triple_131907_148071_205335_cellchat_cldn4/analyze.py`

