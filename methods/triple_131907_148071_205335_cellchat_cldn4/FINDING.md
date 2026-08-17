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

## Ligand table (outgoing CLDN4-high malignant → T/NK)

Not scored in this write-up (matrix step skipped or no pair passed the detect gate). Re-run without `--skip-cellchat`.

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

