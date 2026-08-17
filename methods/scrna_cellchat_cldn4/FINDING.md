# FINDING — CellChat-style CLDN4-high vs low vs T/NK (GSE207422)

Additive to prior TACSTD2 CellChat (`methods/scrna_cellchat`). That folder is taken as given and was **not** re-run. **CLDN4 only.**

**Verdict: mixed ligands, not a TACSTD2-style immune-cold copy.** On public GSE207422 the kept split is the **all-post CLDN4 median** (4,203 vs 4,204 epithelial cells, **12 patients**; T/NK 33,760). Barrier/checkpoint-adjacent outgoing pairs are higher in CLDN4-high (NECTIN2–TIGIT ΔP +0.096; CDH1–ITGAE/ITGB7 +0.041; HLA-E–CD8; JAM1). The same high arm also has **CXCL16–CXCR6 up** (ΔP +0.019) and **IFNG incoming up** (ΔP +0.015). MHC-II→CD4 outgoing is *higher* in CLDN4-high (opposite the prior TACSTD2 NMPR kept split). CD274–PDCD1 and PVR–TIGIT are higher from CLDN4-low. Recruit-up = 1 on median_post **and** on the NMPR tertile, so the TACSTD2 keep-NMPR rule (recruit-up = 0) **does not fire**.

Public UMI only (Hu et al., *Genome Medicine* 2023, PMID 36869384). CellChat R and LIANA were not run. Pairs below use Hill probability + 100 high/low permutations in `scripts/analyze.py`. A pair is significant if detected (`expr_prop ≥ 0.10` on both sides), `P > 0`, and permutation `p < 0.05`. Smallest possible p with 100 permutations is 1/101 = 0.0099.

No edge is drawn unless it meets that rule. The ligand table is [`results/ligand_table.tsv`](results/ligand_table.tsv) (also `ligand_table_kept.tsv`). Extra figure: [`results/fig_extra_ligand_table.png`](results/fig_extra_ligand_table.png). Honest per-sample n: [`results/per_sample_post.tsv`](results/per_sample_post.tsv) and [`fig_n_per_sample.png`](results/fig_n_per_sample.png).

## English

### Data and n

| Item | n | Note |
| --- | ---: | --- |
| Genes × barcodes in the GEO UMI | 24,292 × 92,330 | one matrix |
| CellChatDB v2 protein pairs with every subunit in the matrix | 1,537 | of 2,239 protein pairs |
| Post-treatment samples used for splits | 12 | MPR+pCR = 4; NMPR = 8 |
| Post-treatment epithelial cells | 8,407 | malignant proxy (CopyKAT IDs not public) |
| Post-treatment T+NK | 33,760 | T and NK merged as TNK |
| Pre-treatment biopsies | excluded | 3 samples |

**BD_immune07** (NMPR) is **5,121 / 8,407 (60.9%)** of post-treatment epithelial cells and **5,121 / 6,914 (74.1%)** of NMPR epithelial cells. Group means are **cell-pooled**, so this sample can dominate NMPR and pooled contrasts. Unit of inference for a patient claim would be n=4 vs 8; that test is **not** what the permutation does.

### Per-sample n (post-treatment)

| Sample | Response | Epithelial | T+NK | Mean CLDN4 (epi) |
| --- | --- | ---: | ---: | ---: |
| BD_immune03 | MPR | 804 | 5,060 | 1.08 |
| BD_immune06 | MPR | 64 | 2,924 | 1.56 |
| BD_immune11 | MPR | 439 | 1,690 | 1.90 |
| BD_immune14 | MPR | 186 | 2,863 | 1.47 |
| BD_immune02 | NMPR | 107 | 2,899 | 1.92 |
| BD_immune04 | NMPR | 145 | 5,616 | 1.50 |
| BD_immune07 | NMPR | 5,121 | 1,087 | 1.68 |
| BD_immune09 | NMPR | 271 | 3,283 | 1.72 |
| BD_immune10 | NMPR | 208 | 4,080 | 1.05 |
| BD_immune12 | NMPR | 482 | 1,386 | 0.87 |
| BD_immune13 | NMPR | 35 | 1,347 | 1.66 |
| BD_immune15 | NMPR | 545 | 1,525 | 1.46 |

### Splits compared (CLDN4 `log1p(CP10k)` on epithelium)

Tertile high/low drops the middle third.

| Split | Mal high | Mal low | TNK | Detected LR | Sig LR (out / in) | Cold/barrier | Recruit-up in high |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| tertile_post | 2803 | 2803 | 33760 | 220 | 93 (63 / 30) | 12 | 1 |
| **median_post** | **4203** | **4204** | **33760** | **224** | **100 (68 / 32)** | **13** | **1** |
| tertile_nmpr | 2305 | 2305 | 21223 | 203 | 83 (58 / 25) | 13 | 1 |
| tertile_mpr | 498 | 498 | 12537 | 233 | 111 (77 / 34) | 10 | 3 |
| combo_mpr | 2408+395 | 2228+575 | 21223+12537 | 438 | 214 (152 / 62) | 10 | 4 |

Combo score adds NMPR and MPR arms and is not comparable to a single contrast.

### Kept split

**median_post** (all 12 post-treatment patients). The TACSTD2 keep rule would have taken NMPR-only tertile if that arm had recruit-up = 0 and MPR did not. Here **both** median_post and tertile_nmpr have recruit-up = 1 (CXCL16–CXCR6 higher in CLDN4-high) and the same cold/barrier score (13). median_post is kept because it uses 12 patients rather than 8. MPR-only tertile has recruit-up = 3 and score 10 (not kept). Combo score 10 sums two arms and is not used.

Kept n: CLDN4-high epithelial **4,203** cells / **12** patients (mean CLDN4 2.240); CLDN4-low **4,204** / **12** (mean 0.866); T/NK **33,760** / **12**. Median cut = 1.599 `log1p(CP10k)`.

Cold/barrier score = 13 (barrier-up 5, inhib-up 8, recruit-up 1, recruit-down 0, attack-in-down 3, barrier-down 2).

Signed pairs that matter for the mixed verdict (all p_min = 0.010 unless noted):

| Pair | Dir | ΔP | Note |
| --- | --- | ---: | --- |
| NECTIN2–TIGIT | out | +0.096 | barrier / inhibitory, higher in high |
| CDH1–ITGAE/ITGB7 | out | +0.041 | barrier |
| JAM1–ITGAL/ITGB2 | out | +0.011 | barrier |
| HLA-E–CD8B / CD8A | out | +0.013 / +0.008 | inhibitory |
| CXCL16–CXCR6 | out | +0.019 | **recruit-up in high** |
| HLA-DRB1–CD4 | out | +0.071 | MHC-II higher in high (opposite TACSTD2 NMPR) |
| CD274–PDCD1 | out | −0.006 | higher from CLDN4-low |
| PVR–TIGIT | out | −0.014 | higher from CLDN4-low (p=0.020) |
| IFNG–IFNGR1/2 | in | +0.015 | attack **into high** |
| FASLG–FAS | in | −0.002 | slightly into low (p=0.050) |
| CD8A–CEACAM5 | in | +0.252 | CEACAM5 on CLDN4-high epithelium |

### Ligand table (kept split, significant ΔP only)

Among 6,148 directed tests: detected 224; **significant 100** (outgoing Mal→TNK 68; incoming TNK→Mal 32). Ligand table rows: **100** differential pairs (68 outgoing, 32 incoming; outgoing higher in CLDN4-high = 36, higher in low = 32).

Full table: [`results/ligand_table.tsv`](results/ligand_table.tsv). Extra figure: [`results/fig_extra_ligand_table.png`](results/fig_extra_ligand_table.png).

#### Outgoing Mal → T/NK, higher in CLDN4-high

| Pair | Pathway | Class | ΔP (high−low) | P high | P low | p min |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| NECTIN2_TIGIT | NECTIN | barrier|inhibitory | +0.096 | 0.382 | 0.286 | 0.010 |
| HLA-DRB1_CD4 | MHC-II | other | +0.071 | 0.295 | 0.224 | 0.010 |
| GDF15_TGFBR2 | GDF | other | +0.065 | 0.150 | 0.085 | 0.010 |
| LAMC2_CD44 | LAMININ | other | +0.056 | 0.425 | 0.370 | 0.010 |
| ICAM1_ITGAL_ITGB2 | ICAM | other | +0.054 | 0.152 | 0.098 | 0.010 |
| ICAM1_ITGAL | ICAM | other | +0.052 | 0.146 | 0.094 | 0.010 |
| HLA-DQB1_CD4 | MHC-II | other | +0.047 | 0.099 | 0.052 | 0.010 |
| HLA-DPA1_CD4 | MHC-II | other | +0.044 | 0.257 | 0.213 | 0.010 |
| ICAM1_SPN | ICAM | other | +0.043 | 0.119 | 0.076 | 0.010 |
| CD55_ADGRE5 | ADGRE | other | +0.042 | 0.336 | 0.293 | 0.010 |
| HLA-DRA_CD4 | MHC-II | other | +0.041 | 0.403 | 0.362 | 0.010 |
| CDH1_ITGAE_ITGB7 | CDH1 | barrier|inhibitory | +0.041 | 0.253 | 0.212 | 0.010 |

#### Outgoing Mal → T/NK, higher in CLDN4-low

| Pair | Pathway | Class | ΔP (high−low) | P high | P low | p min |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| LAMB1_CD44 | LAMININ | other | -0.106 | 0.019 | 0.125 | 0.010 |
| COL4A5_CD44 | COLLAGEN | other | -0.091 | 0.168 | 0.259 | 0.010 |
| LAMC1_CD44 | LAMININ | other | -0.085 | 0.204 | 0.289 | 0.010 |
| LAMA3_CD44 | LAMININ | other | -0.060 | 0.157 | 0.217 | 0.010 |
| COL4A5_ITGA1_ITGB1 | COLLAGEN | other | -0.048 | 0.075 | 0.123 | 0.010 |
| LAMC1_ITGA1_ITGB1 | LAMININ | other | -0.047 | 0.093 | 0.140 | 0.010 |
| LAMB1_ITGA1_ITGB1 | LAMININ | other | -0.047 | 0.008 | 0.054 | 0.010 |
| LAMA5_CD44 | LAMININ | other | -0.038 | 0.256 | 0.294 | 0.010 |
| SPP1_CD44 | SPP1 | other | -0.036 | 0.356 | 0.392 | 0.020 |
| LAMA3_ITGA1_ITGB1 | LAMININ | other | -0.031 | 0.069 | 0.100 | 0.010 |
| SPP1_ITGA4_ITGB1 | SPP1 | other | -0.026 | 0.197 | 0.222 | 0.020 |
| LAMB3_CD44 | LAMININ | other | -0.025 | 0.252 | 0.277 | 0.050 |

#### Incoming T/NK → Mal, higher into CLDN4-high

| Pair | Pathway | Class | ΔP (high−low) | P high | P low | p min |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| CD8A_CEACAM5 | CEACAM | other | +0.252 | 0.286 | 0.034 | 0.010 |
| SEMA4D_PLXNB2 | SEMA4 | other | +0.047 | 0.344 | 0.297 | 0.010 |
| CD99_CD99 | CD99 | other | +0.037 | 0.692 | 0.654 | 0.010 |
| GZMA_F2RL1 | PARs | other | +0.031 | 0.073 | 0.042 | 0.010 |
| GZMA_PARD3 | PARs | other | +0.024 | 0.593 | 0.569 | 0.010 |
| HLA-DPB1_CD4 | MHC-II | other | +0.023 | 0.023 | 0.000 | 0.010 |
| CD6_ALCAM | CD6 | other | +0.022 | 0.387 | 0.365 | 0.010 |
| CD96_NECTIN1 | CD96 | other | +0.022 | 0.550 | 0.528 | 0.010 |

#### Incoming T/NK → Mal, higher into CLDN4-low

| Pair | Pathway | Class | ΔP (high−low) | P high | P low | p min |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| CD96_PVR | CD96 | other | -0.018 | 0.022 | 0.040 | 0.020 |
| LGALS9_CD44 | GALECTIN | inhibitory | -0.004 | 0.037 | 0.041 | 0.010 |
| NAMPT_ITGA5_ITGB1 | VISFATIN | other | -0.003 | 0.000 | 0.003 | 0.010 |
| FASL_FAS | FASLG | attack | -0.002 | 0.014 | 0.017 | 0.050 |
| TNFSF14_LTBR | LIGHT | other | -0.002 | 0.031 | 0.033 | 0.010 |
| TNF_TNFRSF1A | TNF | attack | -0.002 | 0.013 | 0.015 | 0.010 |
| LGALS9_P4HB | GALECTIN | inhibitory | -0.002 | 0.055 | 0.056 | 0.010 |
| TNFSF10_TNFRSF10B | TRAIL | attack | -0.002 | 0.012 | 0.014 | 0.050 |

### What is not claimed

- These are permutation tests on cell-pooled truncated means, not a patient-level mixed model. Patient n = 4 MPR vs 8 NMPR.
- BD_immune07 supplies a large share of post epithelial UMIs.
- CellChat R visualizations were not generated. Figures show only pairs that pass p < 0.05.
- Author cell-type and CopyKAT objects were not used.
- TACSTD2 was not used to define high/low. This is not a dual-high or TACSTD2 re-analysis.
- GSE253013 / GSE241934 were not used.

### Files

`results/ligand_table.tsv`, `ligand_table_kept.tsv`, `key_pairs_kept.tsv`, `per_sample_post.tsv`, `summary.json`, `fig_extra_ligand_table.png`, `fig_n_per_sample.png`, `fig_n_*.png`, `fig_nsig_*.png`, `fig_top_*.png`.

---

## 中文

**加性、只切 CLDN4。结论是混合配体，不是 TACSTD2 NMPR 那种免疫冷拷贝。** 先前 TACSTD2 CellChat 当作已给，不重跑。公开 GSE207422 一张 UMI（24,292 基因 × 92,330 细胞）。术后 **MPR 4（含 pCR）/ NMPR 8**。上皮（恶性代理）8,407，T+NK 33,760。**BD_immune07** 占术后上皮 60.9%。均值是细胞池化，不是患者混合模型。

保留切分 **median_post**（12 例；NMPR 三分位同样 recruit-up=1，不能套 TACSTD2 的“只留 NMPR”规则）：CLDN4-high 4,203，low 4,204，T/NK 33,760。NECTIN2–TIGIT / CDH1 在 high 更高，但 CXCL16–CXCR6 和 IFNG 进入 high 也更高；CD274–PDCD1 反而在 low 更高。显著差异配体–受体 **100** 行，见 `results/ligand_table.tsv`。额外图 `fig_extra_ligand_table.png`。
