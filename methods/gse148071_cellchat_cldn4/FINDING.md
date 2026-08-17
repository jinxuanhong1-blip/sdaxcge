# FINDING — CellChat-style CLDN4-high vs low vs T/NK (GSE148071)

**Verdict: cell-pooled only; do not read as a 42-patient result.** On public GSE148071 (Wu et al. 2021), **25 / 42** biopsies have ≥25 putative malignant (epithelial) cells **and** ≥25 T/NK cells. In that eligible pool, CLDN4-high epithelium has higher outgoing **LGALS9–PTPRC**, **NECTIN2–TIGIT**, **HLA-E–CD8A**, **F11R–ITGAL/ITGB2**, and **MDK** toward T/NK (permutation p = 0.0099). **CD274–PDCD1 is not detected.** **CXCL16–CXCR6 is slightly higher** in CLDN4-high (not lower). T/NK is thin: **4,025** cells; **P7 is 26.3%** of them. Cell-pooled truncated means are not a patient mixed model.

This folder is **CLDN4 only**. TACSTD2 is not used to define groups. Epithelial cells are **putative** malignant (no GEO labels / CopyKAT).

---

## English

### Honest n

| Item | n | Note |
| --- | ---: | --- |
| Patients deposited | **42** | Wu et al., *Nat Commun* 2021, PMID 33953163; stage III/IV NSCLC biopsies |
| Patients eligible (≥25 epi **and** ≥25 T/NK) | **25** | unit that can be scored; **not 42** |
| Patients excluded | **17** | mostly epithelium with almost no T/NK (P17 7,322 epi / 0 TNK; P3 7,107 / 17; P41 5,471 / 3) |
| Cells in GEO tar | 89,887 | 42 `GSM*_P*_exp.txt.gz` matrices |
| Epithelial cells, all / eligible | 51,213 / **19,832** | marker-argmax; putative malignant |
| T/NK cells, all / eligible | 4,161 / **4,025** | T + NK merged |
| CellChatDB v2 protein pairs in matrix | 1,905 | of 2,239 |
| Kept split | **median** | CLDN4 `log1p(CP10k)` cut **1.39** among eligible epithelium |
| Mal_high / Mal_low / TNK (kept) | **9,916 / 9,916 / 4,025** | 25 patients each group |
| Detected / significant directed tests | 112 / 52 | 100 label permutations; smallest p = 1/101 = 0.0099 |
| Significant differential pairs | **52** | 39 outgoing Mal→TNK, 13 incoming TNK→Mal |

Eligible epithelium is dominated by **P1 (16.9%)**, P10 (12.8%), P6 (11.5%), P9 (10.7%), P18 (9.5%) — top five = **61%** of 19,832 cells. Eligible T/NK is dominated by **P7 (1,060 / 4,025 = 26.3%)**, P42 (11.3%), P40 (9.3%). Per-sample counts: `results/per_sample.tsv`.

Tertile (cuts 0.73 / 1.89; 6,614 / 6,611 / 4,025) has the same cold/barrier score (10) and the same recruit-up count (1). Median is kept only because it uses more malignant cells.

### Ligand table — outgoing Mal CLDN4-high → T/NK (kept median)

Pairs with permutation p < 0.05 on at least one arm and ΔP ≠ 0. Full table: `results/ligand_table.tsv`.

**Higher from CLDN4-high (ΔP = P_high − P_low > 0):**

| Pair | Pathway | Class | P high | P low | ΔP | p high |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| LGALS9–PTPRC (CD45) | GALECTIN | inhibitory | 0.095 | 0 | +0.095 | 0.0099 |
| MDK–NCL | MK | other | 0.733 | 0.670 | +0.064 | 0.0099 |
| MDK–ITGA4/ITGB1 | MK | other | 0.285 | 0.228 | +0.058 | 0.0099 |
| CLEC2B–KLRB1 | CLEC | other | 0.070 | 0.034 | +0.036 | 0.0099 |
| HLA-E–CD8A | MHC-I | inhibitory | 0.284 | 0.254 | +0.031 | 0.0099 |
| HLA-C–CD8A | MHC-I | other | 0.349 | 0.319 | +0.030 | 0.0099 |
| LGALS9–CD44 | GALECTIN | inhibitory | 0.022 | 0 | +0.022 | 0.0099 |
| HLA-B–CD8A | MHC-I | other | 0.431 | 0.409 | +0.021 | 0.0099 |
| **NECTIN2–TIGIT** | NECTIN | barrier\|inhibitory | 0.034 | 0.014 | +0.020 | 0.0099 |
| CD55–ADGRE5 | ADGRE | other | 0.077 | 0.058 | +0.019 | 0.0099 |
| HLA-F–CD8A | MHC-I | inhibitory | 0.051 | 0.034 | +0.017 | 0.0099 |
| COL1A1–CD44 | COLLAGEN | other | 0.071 | 0.058 | +0.013 | 0.0099 |
| **F11R–ITGAL/ITGB2** | JAM | barrier | 0.081 | 0.069 | +0.012 | 0.0099 |
| ICAM1–ITGAL/ITGB2 | ICAM | other | 0.058 | 0.048 | +0.009 | 0.0099 |
| CXCL16–CXCR6 | CXCL | recruit | 0.003 | 0.001 | +0.002 | 0.0099 |

**Higher from CLDN4-low (ΔP < 0; significant on the low arm):**

| Pair | Pathway | P high | P low | ΔP | p low |
| --- | --- | ---: | ---: | ---: | ---: |
| SPP1–CD44 | SPP1 | 0.274 | 0.415 | −0.141 | 0.0099 |
| SPP1–ITGA4/ITGB1 | SPP1 | 0.171 | 0.280 | −0.109 | 0.0099 |
| LAMB3–CD44 | LAMININ | 0.237 | 0.311 | −0.074 | 0.0099 |
| MIF–CD74/CD44 | MIF | 0.240 | 0.314 | −0.073 | 0.0099 |
| FN1–CD44 | FN1 | 0.005 | 0.071 | −0.066 | 0.0099 |
| MIF–CD74/CXCR4 | MIF | 0.193 | 0.257 | −0.064 | 0.0099 |

### Incoming T/NK → Mal (kept)

CD8A–CEACAM5 is higher into CLDN4-high (ΔP +0.092; CEACAM5 is the epithelial receptor). TNF–TNFRSF1A is **lower** into high (ΔP −0.009). GZMA–PARD3 is higher into high (ΔP +0.054).

### What is not supported

- **n = 42** as the communication n. Seventeen biopsies cannot be scored.
- **CD274–PDCD1**, PVR–TIGIT, CXCL9/10–CXCR3: not detected (`expr_prop < 0.10` on at least one side).
- A clean “CLDN4-high = chemokine-low” story: CXCL16–CXCR6 is **up** in high.
- Patient-level inference, ICI response, or histology split (GEO has age/sex only).
- Author malignant IDs. Epithelium is marker-argmax.

### Files

`results/ligand_table.tsv`, `n_table.tsv`, `per_sample.tsv`, `fig_extra_ligand_table.png`, `fig_n_per_sample.png`, `summary.json`.

---

## 中文

**结论：只能报告细胞池化结果，不能写成 42 例。** GSE148071 存档 42 例晚期 NSCLC 活检，**只有 25 例**同时有 ≥25 个上皮（代理恶性）和 ≥25 个 T/NK。这 25 例里，CLDN4 高上皮对 T/NK 的外向配体更高的是 **LGALS9–PTPRC、NECTIN2–TIGIT、HLA-E–CD8A、F11R、MDK**（置换 p=0.0099）。**CD274–PDCD1 未检出。** CXCL16–CXCR6 在高组略高，不是更低。T/NK 一共 **4,025** 个细胞，**P7 占 26.3%**。推断单位不是患者。

上皮无 GEO 标签，是 marker-argmax 的**假定恶性**。排除的 17 例里有大批上皮、几乎没有 T/NK（P17 / P3 / P41）。
