# FINDING — GSE131907 LIANA/LR from CLDN4-high malignant to T/NK

**Verdict:** On 29 paired tumor samples (29 patients), CLDN4-high vs CLDN4-low malignant → T/NK focus pairs do **not** support a coordinated T-recruit drop (3/17 median Δ < 0; **0 FDR < 0.05** on the down side). MHC-I outgoing is **higher** from CLDN4-high (9 pairs FDR < 0.05, Δ > 0). CXCL9/10/11–CXCR3 fail the 10% expression filter in malignant cells. The LR table is a ranked co-expression list, not a causal claim.

Kim et al., *Nat Commun* 2020 (PMID 32385277); GEO [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907). Treatment-naive LUAD atlas. **No ICI / MPR labels.**

---

## English

### Honest n

| Item | n | Note |
| --- | ---: | --- |
| Cells in UMI matrix | 208,506 | author barcodes |
| GEO samples / patients | 58 / 44 | 58 / 44 in the series |
| Tumor-origin samples | 37 | tLung, tL/B, mLN, mBrain, PE |
| Author malignant cells | 31,136 | subtype ∈ {Malignant cells, tS1, tS2, tS3} |
|  … tLung tS1/tS2/tS3 | 6,352 | Kim tumor-specific epi, not the 'Malignant cells' label |
|  … author Malignant cells | 24,784 | mets / tL-B / mLN / mBrain |
| CLDN4-high / low malignant | 15,569 / 15,567 | global median log1p(CP10k) = 1.798 |
| T / NK / T+NK (all samples) | 79,676 / 11,551 / 91,227 | author Cell_type |
| T/NK in tumor-origin samples | 48,012 | receivers for the pooled LR table |
| Samples with any malignant + T/NK | 32 | descriptive |
| **Paired samples (unit of test)** | **29** | ≥10 high, ≥10 low, ≥20 T/NK |
| **Unique patients in paired set** | **29** | do not count cells as n |
| Pairs scored / pass expr_prop 0.10 | 2735 / 53 | CellPhoneDB v5 + overlay |
| LIANA | ran_cellphonedb_method | secondary |
| CellChat | not_run_R_unavailable | not run |

PE epithelial cells are **unlabeled** in the author file (0 `Malignant cells`); they are **not** counted as malignant. nLung AT1/AT2/Club/Ciliated are **not** malignant. Patients can contribute more than one tumor site; paired unique-patient n is 29, not 44. Dropped from the paired gate: LUNG_T09 (5 malignant), NS_16 (79 malignant, 0 CLDN4-high), EBUS_13 (376 malignant, 0 CLDN4-high).

Per-sample counts: `results/n_cells_samples.tsv`.

### Method

Documented CellPhoneDB-style score on log1p(CP10k): partner expression = **min of subunit means**; pair score = **mean of the two partner means** (Efremova et al. 2020; Garcia-Alonso et al. 2022). `pass_expr_prop` requires both partners in ≥10% of cells in their group. Pooled scores are descriptive. The test is a **paired Wilcoxon** of high vs low **per tumor sample**. FDR is Benjamini–Hochberg within the outgoing contrast. This is **not** a CellChat communication probability and does **not** observe secretion or spatial contact.

### LR table — CLDN4-high malignant → T/NK (pooled, pass expr_prop)

Full table: `results/lr_table.tsv`. Top 25 by score:

| ligand | receptor | pathway | ligand_frac | receptor_frac | cpdb_mean_score |
|---|---|---|---|---|---|
| PPIA | BSG | other | 0.96 | 0.33 | 1.439 |
| CXCL14 | CXCR4 | T_recruit_partial | 0.31 | 0.79 | 1.401 |
| APP | CD74 | other | 0.70 | 0.77 | 1.313 |
| HLA-A | CD8A | MHC_I | 0.85 | 0.24 | 1.254 |
| HLA-B | CD8A | MHC_I | 0.87 | 0.24 | 1.206 |
| HLA-A | CD8B | MHC_I | 0.85 | 0.20 | 1.204 |
| HLA-B | CD8B | MHC_I | 0.87 | 0.20 | 1.156 |
| HLA-C | CD8A | MHC_I | 0.85 | 0.24 | 1.114 |
| HLA-C | CD8B | MHC_I | 0.85 | 0.20 | 1.063 |
| CD58 | CD2 | other | 0.18 | 0.73 | 0.832 |
| SEMA4D | PTPRC | other | 0.11 | 0.69 | 0.693 |
| PTGES3 | PTGER4 | other | 0.78 | 0.21 | 0.657 |
| HLA-E | KLRD1 | MHC_I | 0.70 | 0.18 | 0.619 |
| IGFBP3 | TMEM219 | other | 0.46 | 0.27 | 0.606 |
| PTGES3 | PTGER2 | other | 0.78 | 0.10 | 0.566 |
| CD55 | ADGRE5 | other | 0.64 | 0.23 | 0.563 |
| APP | SORL1 | other | 0.70 | 0.16 | 0.536 |
| CD47 | SIRPG | other | 0.68 | 0.11 | 0.529 |
| DHCR24 | RORA | other | 0.47 | 0.34 | 0.501 |
| ICAM1 | SPN | other | 0.52 | 0.21 | 0.469 |
| CDH1 | ITGAE+ITGB7 | other | 0.66 | 0.14 | 0.466 |
| CDH1 | KLRG1 | other | 0.66 | 0.11 | 0.456 |
| ALCAM | CD6 | other | 0.45 | 0.30 | 0.444 |
| CD44 | TYROBP | other | 0.40 | 0.22 | 0.434 |
| CEACAM5 | CD8A | MHC_I_partial | 0.27 | 0.24 | 0.405 |

### Focus axes (T-recruit / IFN / MHC-I) that pass expr_prop

| ligand | receptor | pathway | ligand_frac | receptor_frac | cpdb_mean_score |
|---|---|---|---|---|---|
| HLA-A | CD8A | MHC_I | 0.85 | 0.24 | 1.254 |
| HLA-B | CD8A | MHC_I | 0.87 | 0.24 | 1.206 |
| HLA-A | CD8B | MHC_I | 0.85 | 0.20 | 1.204 |
| HLA-B | CD8B | MHC_I | 0.87 | 0.20 | 1.156 |
| HLA-C | CD8A | MHC_I | 0.85 | 0.24 | 1.114 |
| HLA-C | CD8B | MHC_I | 0.85 | 0.20 | 1.063 |
| HLA-E | KLRD1 | MHC_I | 0.70 | 0.18 | 0.619 |
| CXCL16 | CXCR6 | T_recruit | 0.40 | 0.12 | 0.261 |

### Patient/sample-level high vs low (outgoing)

Median Δ = CLDN4-high − CLDN4-low. Negative = weaker from the high state. Paired n = **29 samples / 29 patients**.

| ligand | receptor | pathway | n_samples | median_delta | pval | padj |
|---|---|---|---|---|---|---|
| HLA-B | CD8A | MHC_I | 27 | +0.082 | 0.000536 | 0.00264 |
| HLA-B | CD8B | MHC_I | 23 | +0.082 | 0.00214 | 0.00804 |
| HLA-A | CD8A | MHC_I | 27 | +0.087 | 5.49e-05 | 0.000536 |
| HLA-E | KLRC2 | MHC_I | 4 | +0.089 | 0.125 | 0.215 |
| HLA-E | KLRC2+KLRD1 | MHC_I | 4 | +0.089 | 0.125 | 0.215 |
| HLA-C | CD8A | MHC_I | 27 | +0.091 | 3.02e-05 | 0.000344 |
| HLA-E | KLRD1 | MHC_I | 24 | +0.092 | 2.38e-07 | 9.42e-06 |
| HLA-A | CD8B | MHC_I | 23 | +0.093 | 0.000253 | 0.00154 |
| HLA-C | CD8B | MHC_I | 23 | +0.102 | 0.000253 | 0.00154 |
| HLA-E | KLRC1 | MHC_I | 9 | +0.120 | 0.00391 | 0.0129 |
| HLA-E | KLRC1+KLRD1 | MHC_I | 8 | +0.123 | 0.00781 | 0.0232 |
| CEACAM5 | CD8A | MHC_I_partial | 18 | +0.033 | 0.0539 | 0.106 |
| CXCL10 | CXCR3 | T_recruit | 4 | -0.012 | 0.375 | 0.502 |
| CCL4 | CCR5 | T_recruit | 3 | -0.003 | 0.5 | 0.59 |
| CCL3 | CCR5 | T_recruit | 3 | -0.003 | 1 | 1 |
| CXCL12 | CXCR4 | T_recruit | 3 | +0.017 | 0.25 | 0.373 |
| CX3CL1 | CX3CR1 | T_recruit | 4 | +0.025 | 0.125 | 0.215 |
| CXCL16 | CXCR6 | T_recruit | 18 | +0.026 | 0.0342 | 0.0712 |
| SFTPD | ADGRE5 | other | 17 | -0.046 | 0.0887 | 0.171 |
| PPIA | BSG | other | 29 | -0.024 | 0.0274 | 0.0601 |
| SEMA4D | PTPRC | other | 21 | +0.005 | 0.029 | 0.062 |
| CD58 | CD2 | other | 23 | +0.007 | 0.0112 | 0.0304 |
| LIPA | RORA | other | 26 | +0.007 | 0.00794 | 0.0232 |
| LPAR2 | ADGRE5 | other | 24 | +0.008 | 0.015 | 0.0374 |
| JAG1 | CD46 | other | 9 | +0.012 | 0.00391 | 0.0129 |
| ITGAV+ITGB1 | ADGRE5 | other | 28 | +0.012 | 0.00518 | 0.0164 |
| CLEC2B | KLRF1 | other | 8 | +0.015 | 0.0234 | 0.0545 |
| PTGES | PTGER4 | other | 16 | +0.016 | 0.00058 | 0.00269 |
| PTGES | PTGER2 | other | 9 | +0.016 | 0.0117 | 0.0309 |
| PODXL2 | SELL | other | 15 | +0.017 | 0.00116 | 0.00482 |
| LGALS9 | HAVCR2 | other | 10 | +0.018 | 0.0273 | 0.0601 |
| CD44 | TYROBP | other | 26 | +0.018 | 0.0435 | 0.088 |
| LGALS9 | P4HB | other | 24 | +0.018 | 0.000494 | 0.0026 |
| ALCAM | CD6 | other | 27 | +0.020 | 0.0104 | 0.0293 |
| BAG6 | NCR3 | other | 9 | +0.020 | 0.00391 | 0.0129 |
| TGM2 | ADGRG1 | other | 7 | +0.028 | 0.0156 | 0.0374 |
| ICAM4 | ITGAL+ITGB2 | other | 15 | +0.032 | 6.1e-05 | 0.000536 |
| APP | CD74 | other | 28 | +0.034 | 0.000792 | 0.00348 |
| PLAUR | ITGA4+ITGB1 | other | 16 | +0.056 | 0.0155 | 0.0374 |
| APP | SORL1 | other | 17 | +0.059 | 0.00168 | 0.00663 |
| CD55 | ADGRE5 | other | 28 | +0.067 | 0.000194 | 0.00153 |
| F11R | ITGAL+ITGB2 | other | 17 | +0.072 | 3.05e-05 | 0.000344 |
| ICAM1 | SPN | other | 27 | +0.073 | 1.49e-08 | 1.18e-06 |
| CDH1 | ITGAE+ITGB7 | other | 23 | +0.076 | 7.15e-07 | 1.88e-05 |
| IGFBP3 | TMEM219 | other | 19 | +0.080 | 0.00042 | 0.00237 |
| CDH1 | KLRG1 | other | 14 | +0.083 | 0.000244 | 0.00154 |
| ICAM1 | ITGAL | other | 17 | +0.109 | 1.53e-05 | 0.000241 |
| ICAM1 | ITGAL+ITGB2 | other | 17 | +0.109 | 1.53e-05 | 0.000241 |

### LIANA (`mt.cellphonedb`)

Status: `ran_cellphonedb_method`. If a LIANA table is present, p-values are within-object specificity on a downsampled object (≤2,000 cells/group, 50 permutations), **not** the patient-level test above.

| source | target | ligand_complex | receptor_complex | lr_means | cellphone_pvals |
|---|---|---|---|---|---|
| Malig_CLDN4high | T | APP | CD74 | 1.33 | 0 |
| Malig_CLDN4high | NK | APP | CD74 | 1.31 | 0 |
| Malig_CLDN4high | NK | MDK | SORL1 | 1.12 | 0 |
| Malig_CLDN4high | NK | HLA-B | KIR3DL2 | 1.08 | 1 |
| Malig_CLDN4high | T | MDK | SORL1 | 1.06 | 0 |
| Malig_CLDN4high | NK | HLA-C | KIR2DL3 | 1 | 1 |
| Malig_CLDN4high | T | CD58 | CD2 | 0.861 | 0 |
| Malig_CLDN4high | NK | HLA-E | KLRC1_KLRD1 | 0.745 | 1 |
| Malig_CLDN4high | NK | HLA-E | KLRC1 | 0.745 | 1 |
| Malig_CLDN4high | NK | CLEC2B | KLRF1 | 0.701 | 0 |
| Malig_CLDN4high | T | LGALS9 | CD44 | 0.686 | 0 |
| Malig_CLDN4high | T | HBEGF | CD44 | 0.658 | 0 |
| Malig_CLDN4high | T | SPP1 | CD44 | 0.658 | 0 |
| Malig_CLDN4high | T | FGFR2 | CD44 | 0.639 | 0 |
| Malig_CLDN4high | NK | HLA-E | KLRC2_KLRD1 | 0.598 | 1 |

### What this is not

- Not ICI / MPR (GSE131907 is treatment-naive).
- Not inferCNV/CopyKAT recomputed malignant IDs.
- Not CellChat. LIANA permutation p-values (if present) are within-object specificity, not patient tests.
- Not spatial proximity or protein secretion.
- Cells are not the sample size.

### Files

`results/lr_table.tsv` (all scored pairs), `results/lr_table_pass.tsv` (53 pairs passing expr_prop), `results/lr_table_high_vs_low.tsv`, `results/liana_cldn4high_to_tnk.csv`, `results/n_table.tsv`, `results/n_cells_samples.tsv`, `results/summary.json`, `results/figures/`.

---

## 中文

**结论：** 在 **29 个配对肿瘤样本 / 29 名患者** 上，CLDN4 高 vs 低恶性细胞 → T/NK 的焦点对**不支持**协同招募下降（3/17 中位 Δ < 0；下行 **0 个 FDR < 0.05**）。MHC-I 出站在 CLDN4 高侧更高（9 对 FDR < 0.05）。CXCL9/10/11–CXCR3 未过恶性细胞 10% 表达门槛。LR 表是共表达排序，不是因果机制。

GSE131907（Kim 2020）治疗初治 LUAD。作者注释恶性细胞 31,136 个（tLung 用 tS1/tS2/tS3；转移灶用 Malignant cells）。CLDN4 按恶性细胞全局中位数拆高/低。配对检验单位是肿瘤样本，不是细胞：**29 个样本 / 29 名患者**。
胸水上皮在作者文件里未标恶性，不计入。正常肺 AT1/AT2/Club/Ciliated 不是恶性。
LR 表是共表达排序，不是招募机制。
