# FINDING — Pair GSE207422 + GSE131907, CLDN4-only (patient T/NK then CellChat outgoing)

**Additive. CLDN4 only. No dual-high.** GSE207422 is included. Prior single-cohort
CellChat folders (`methods/scrna_cellchat_cldn4`, `methods/gse131907_cellchat_cldn4`)
and the malignant Q4 T/NK extract (`methods/cldn4_malig_q4_tnk`) are taken as given
and were **not** re-run. Matrices were not re-downloaded.

Primary tables: [`results/combo_rho.tsv`](results/combo_rho.tsv) and
[`results/combo_lr_table.tsv`](results/combo_lr_table.tsv)
(also `ligand_table.tsv`). Extra figures under [`figures/`](figures/).

## 1. Patient-level malignant CLDN4 vs T/NK (honest n)

Unit is the **patient** on GSE207422 and the **sample** on GSE131907 (Kim atlas
has one tumor-origin sample per patient in the author-malignant extract).
Quartiles are **within cohort**, then Q4 vs Q1 T/NK is pooled as rank-biserial *r*
on Fisher-z (DerSimonian–Laird). Spearman is the same RE on Fisher-z(ρ).
p-values are descriptive.

| Item | Public? | n | Note |
|---|---|---:|---|
| GSE207422 patients (DRMref malignant) | yes | **12** | locked A3 table; post-treatment resections |
| GSE207422 Q4 vs Q1 tails | yes | **3 vs 3** | n=12 so each tail is 3; poolable but thin |
| GSE131907 author-malignant samples | yes | **21** | `n_malignant ≥ 20`; origins mBrain,mLN,tL/B |
| GSE131907 Q4 vs Q1 tails | yes | **6 vs 5** | sample-level |
| **Pair Spearman N** | yes | **33** | 12 patients + 21 samples |
| **Pair Q4+Q1 compared** | yes | **17** | 9 Q1 + 8 Q4 |
| tLung tS1–tS3 in the primary T/NK row | **no** | 0 | primary tLung epithelium is not `Malignant cells` |
| Dual-high TACSTD2×CLDN4 | **no** | 0 | CLDN4-only |
| Cell-pooled CellChat n (below) | yes | see §2 | not a 44-patient mixed model |

Do not write n=44 (Kim series patients) or n=15 (Hu series samples) for the pair T/NK test.
The computable pair n is **33** for Spearman and **17** compared tails for Q4 vs Q1.

### Singles

| cohort | malig | score | unit | n | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) |
|---|---|---|---|---:|---|---|
| GSE207422 | author_DRMref | mean | patient | 12 | -0.091 (0.779) | -0.333 (0.7; 3/3) |
| GSE131907 | author_malig | mean | sample | 21 | -0.396 (0.0755) | -0.533 (0.177; 6/5) |
| GSE131907 (sensitivity) | author_malig | pct | sample | 21 | -0.522 (0.0152) | -0.600 (0.126; 6/5) |
| GSE131907 tLung (not primary) | author_epi | mean | sample | 11 | -0.045 (0.894) | -0.111 (1; 3/3) |

### Combo (primary = mean / mean)

| analysis | k | N | effect | p | I² | 95% CI | Stouffer p |
|---|---:|---:|---|---|---:|---|---|
| Spearman ρ | 2 | 33 | **-0.300** | 0.108 | 0% | [-0.596, +0.067] | 0.107 |
| Q4 vs Q1 r | 2 | 17 | **-0.483** | 0.0804 | 0% | [-0.807, +0.064] | 0.177 |
| Spearman ρ (131907 %pos sensitivity) | 2 | 33 | -0.373 | 0.0984 | 30% | — | — |
| Q4 vs Q1 r (131907 %pos sensitivity) | 2 | 17 | -0.536 | 0.0471 | 0% | — | — |

GSE207422 alone is near-null (ρ = −0.09, n=12). GSE131907 author-malignant mean is
negative and larger. The pair stays **CLDN4-negative vs T/NK** with I² = 0, but the
Q4 vs Q1 p is above 0.05 on the primary mean/mean cut. The %pos sensitivity on
GSE131907 is the stronger single (see table). That is the honest pair, not a hidden n.

## 2. CellChat-style outgoing CLDN4-high → T/NK (combo LR)

CellChat R was not run. Each cohort already has Jin et al. 2021 Hill probability on
CellChatDB v2 protein pairs + 100 permutations of CLDN4-high/low among malignant
cells (`expr_prop ≥ 0.10`, p < 0.05; smallest p = 1/101 = 0.0099).

| Cohort | Kept split | Mal high / low | T/NK | Patients/samples | Detected outgoing | Sig outgoing |
|---|---|---:|---:|---:|---:|---:|
| GSE207422 | median_post (epithelial proxy) | 4,203 / 4,204 | 33,760 | 12 | 82 | 68 |
| GSE131907 | tumor tertile (author malig + tS*) | 10,379 / 10,379 | 34,741 | 32 | 61 | 60 |

Means are **cell-pooled**. BD_immune07 is ~61% of GSE207422 post epithelium;
EBUS_28 is ~27% of GSE131907 Mal_high. That is not a 12+32 patient mixed model.

### Concordance (outgoing, significant in ≥1 cohort)

| Concordance | n pairs |
|---|---:|
| both significant, same direction | **26** |
| both significant, discordant | 18 |
| only GSE207422 | 24 |
| only GSE131907 | 16 |

### Shared same-direction outgoing (higher in CLDN4-high)

| Pair | Pathway | Class | Concordance | ΔP 207422 | ΔP 131907 | arm 422 | arm 907 |
|---|---|---|---|---:|---:|---|---|
| CD99_CD99 | CD99 | other | both_sig_same | +0.037 | +0.078 | high | high |
| ICAM1_ITGAL_ITGB2 | ICAM | other | both_sig_same | +0.054 | +0.045 | high | high |
| ICAM1_SPN | ICAM | other | both_sig_same | +0.043 | +0.054 | high | high |
| LAMC2_CD44 | LAMININ | other | both_sig_same | +0.056 | +0.035 | high | high |
| CDH1_ITGAE_ITGB7 | CDH1 | barrier|inhibitory | both_sig_same | +0.041 | +0.039 | high | high |
| HLA-E_CD8A | MHC-I | inhibitory | both_sig_same | +0.008 | +0.057 | high | high |
| ICAM1_ITGAL | ICAM | other | both_sig_same | +0.052 | +0.011 | high | high |
| HLA-E_CD8B | MHC-I | inhibitory | both_sig_same | +0.013 | +0.046 | high | high |
| MDK_ITGA4_ITGB1 | MK | other | both_sig_same | +0.027 | +0.028 | high | high |
| MDK_NCL | MK | other | both_sig_same | +0.023 | +0.030 | high | high |
| HLA-F_CD8A | MHC-I | inhibitory | both_sig_same | +0.013 | +0.039 | high | high |
| JAM1_ITGAL_ITGB2 | JAM | barrier | both_sig_same | +0.011 | +0.031 | high | high |
| HLA-DRB5_CD4 | MHC-II | other | both_sig_same | +0.035 | +0.000 | high | high |
| HLA-F_CD8B | MHC-I | inhibitory | both_sig_same | +0.008 | +0.025 | high | high |
| FN1_CD44 | FN1 | other | both_sig_same | +0.002 | +0.030 | high | high |
| CXCL16_CXCR6 | CXCL | recruit | both_sig_same | +0.019 | +0.007 | high | high |
| CDH1_KLRG1 | CDH1 | barrier|inhibitory | both_sig_same | +0.002 | +0.006 | high | high |
| FN1_ITGA4_ITGB1 | FN1 | other | both_sig_same | +0.001 | +0.004 | high | high |
| FN1_ITGA4_ITGB7 | FN1 | other | both_sig_same | +0.001 | +0.003 | high | high |

### Shared same-direction outgoing (higher in CLDN4-low)

| Pair | Pathway | Class | Concordance | ΔP 207422 | ΔP 131907 | arm 422 | arm 907 |
|---|---|---|---|---:|---:|---|---|
| SPP1_CD44 | SPP1 | other | both_sig_same | -0.036 | -0.108 | low | low |
| LAMC1_CD44 | LAMININ | other | both_sig_same | -0.085 | -0.018 | low | low |
| MIF_CD74_CD44 | MIF | other | both_sig_same | -0.006 | -0.035 | low | low |
| SPP1_ITGA4_ITGB1 | SPP1 | other | both_sig_same | -0.026 | -0.014 | low | low |
| MIF_CD74_CXCR4 | MIF | other | both_sig_same | -0.009 | -0.028 | low | low |
| PVR_TIGIT | PVR | barrier|inhibitory | both_sig_same | -0.014 | -0.003 | low | low |
| LAMB2_CD44 | LAMININ | other | both_sig_same | -0.003 | -0.013 | low | low |

### Pre-specified key pairs (always shown, even if undetected)

| Pair | Pathway | Class | Concordance | ΔP 207422 | ΔP 131907 | arm 422 | arm 907 |
|---|---|---|---|---:|---:|---|---|
| CCL5_CCR5 | CCL | recruit | not_differential | +0.000 | +0.000 | ns | ns |
| CD274_PDCD1 | PD-L1 | inhibitory | only_207422_undetected_131907 | -0.006 | +0.000 | low | ns |
| CDH1_ITGAE_ITGB7 | CDH1 | barrier|inhibitory | both_sig_same | +0.041 | +0.039 | high | high |
| CDH1_KLRG1 | CDH1 | barrier|inhibitory | both_sig_same | +0.002 | +0.006 | high | high |
| CXCL10_CXCR3 | CXCL | recruit | not_differential | +0.000 | +0.000 | ns | ns |
| CXCL16_CXCR6 | CXCL | recruit | both_sig_same | +0.019 | +0.007 | high | high |
| CXCL9_CXCR3 | CXCL | recruit | not_differential | +0.000 | +0.000 | ns | ns |
| HLA-E_CD8A | MHC-I | inhibitory | both_sig_same | +0.008 | +0.057 | high | high |
| HLA-E_CD8B | MHC-I | inhibitory | both_sig_same | +0.013 | +0.046 | high | high |
| HLA-G_CD8A | MHC-I | inhibitory | only_131907_undetected_207422 | +0.000 | +0.006 | ns | high |
| HLA-G_CD8B | MHC-I | inhibitory | only_131907_undetected_207422 | +0.000 | +0.004 | ns | high |
| JAM1_ITGAL_ITGB2 | JAM | barrier | both_sig_same | +0.011 | +0.031 | high | high |
| LGALS9_CD44 | GALECTIN | inhibitory | only_131907 | +0.013 | +0.038 | ns | high |
| LGALS9_CD45 | GALECTIN | inhibitory | only_131907 | +0.016 | +0.043 | ns | high |
| MIF_CD74_CD44 | MIF | other | both_sig_same | -0.006 | -0.035 | low | low |
| MIF_CD74_CXCR4 | MIF | other | both_sig_same | -0.009 | -0.028 | low | low |
| NECTIN2_CD226 | NECTIN | barrier|inhibitory | only_207422_undetected_131907 | +0.022 | NA | high | ns |
| NECTIN2_TIGIT | NECTIN | barrier|inhibitory | only_207422_undetected_131907 | +0.096 | NA | high | ns |
| PVR_TIGIT | PVR | barrier|inhibitory | both_sig_same | -0.014 | -0.003 | low | low |
| SPP1_CD44 | SPP1 | other | both_sig_same | -0.036 | -0.108 | low | low |

Shared barrier / inhibitory outgoing that stay **higher in CLDN4-high** on both
cohorts (when detected): CDH1–ITGAE/ITGB7, JAM1–ITGAL/ITGB2, HLA-E–CD8A/B.
CXCL16–CXCR6 is also **higher** in CLDN4-high on both kept splits (recruit-up,
not recruit-down). NECTIN2–TIGIT and CD274–PDCD1 are GSE207422-only (not detected
on GSE131907). Do not invent those edges on the Kim atlas.

Full combo table: [`results/combo_lr_table.tsv`](results/combo_lr_table.tsv).
Significant-in-either outgoing slice: [`results/ligand_table.tsv`](results/ligand_table.tsv).

## Extra figures

- [`figures/fig_honest_n.png`](figures/fig_honest_n.png) — pair n vs compared tails
- [`figures/fig_combo_rho_forest.png`](figures/fig_combo_rho_forest.png)
- [`figures/fig_combo_q4q1_forest.png`](figures/fig_combo_q4q1_forest.png)
- [`figures/fig_scatter_GSE207422.png`](figures/fig_scatter_GSE207422.png)
- [`figures/fig_scatter_GSE131907.png`](figures/fig_scatter_GSE131907.png)
- [`figures/fig_q4q1_box_GSE207422.png`](figures/fig_q4q1_box_GSE207422.png)
- [`figures/fig_q4q1_box_GSE131907.png`](figures/fig_q4q1_box_GSE131907.png)
- [`figures/fig_extra_lr_concordance.png`](figures/fig_extra_lr_concordance.png)
- [`figures/fig_extra_shared_outgoing.png`](figures/fig_extra_shared_outgoing.png)
- [`figures/fig_extra_ligand_table.png`](figures/fig_extra_ligand_table.png)

## What is not claimed

- This is **not** dual-high TACSTD2×CLDN4 and **not** a 6-unit / 4-unit search.
- GSE207422 CellChat uses marker epithelium (CopyKAT IDs are not public); the
  patient T/NK row uses DRMref malignant. Those are different malignant calls.
- GSE131907 CellChat includes tLung tS1–tS3 (32 samples). The patient T/NK row
  uses author `Malignant cells` with n_mal≥20 (21 samples; mets-heavy).
- Permutation tests on cell-pooled truncated means are not a 33-unit mixed model.
- CXCL16–CXCR6 up in CLDN4-high is the opposite of an immune-cold recruit-down story.
- NECTIN2–TIGIT / PD-L1–PD-1 are not detected on GSE131907.

## Reproduce

```bash
python3 methods/pair_207422_131907_cellchat_cldn4/analyze.py
```

Requires the committed `data/` extracts only. GEO UMI matrices are not needed.

