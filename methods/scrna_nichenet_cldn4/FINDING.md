# FINDING — NicheNet-style ligands from CLDN4-high malignant cells

**Verdict: the signed claim is not supported.** On public GSE207422 (n=4 MPR including pCR, n=8 NMPR), T/NK cytotoxicity is **not** lower in NMPR (mean 0.786 vs 0.745, U=16, **p=1.0**) and T/NK IFN is **not** significantly different (0.381 vs 0.245, U=23, **p=0.28**; NMPR is higher, not lower). Malignant CLDN4 vs T/NK cytotoxicity is **ρ=−0.11, p=0.73**; vs IFN **ρ=+0.20, p=0.53**.

What the prior *does* return, honestly: among 147 potential ligands from CLDN4-high malignant cells, the a priori **IFN** set is recovered far more cleanly than cytotoxicity. Top ligand **IFITM1** (Pearson 0.533, AUROC 0.94, p=1.1e−21). That is prior structure (ISG / MHC ligands sit next to IFN genes in NicheNet-v2), **not** evidence that CLDN4-high cells induce or repress IFN in this cohort. Combinatorial MPR vs NMPR **activity ranks of shared ligands are identical** (same prior × same gene set). The combinatorial difference is the **expressed ligand set**, not a new ranking.

GSE253013 was **not used** (9.3 GB; no MPR/NMPR). TACSTD2 was **not** used to call senders (companion dual-high n=674).

---

## English

### Data and n

| Item | n | Note |
| --- | --- | --- |
| Post-tx patients | 12 | Hu et al. GSE207422; pCR=P06 counted as MPR |
| MPR / NMPR patients | **4 / 8** | unit of inference |
| DRMref-annotated cells | 30,877 | Liu et al. *NAR* 2024; barcode-matched to GEO |
| Malignant cells | 2,051 | not Hu CopyKAT |
| CLDN4-high malignant | **1,026** | CLDN4 `log1p(CP10k)` ≥ malignant median (1.39) |
| Dual-high companion (not used) | 674 | TACSTD2 and CLDN4 both ≥ median |
| CLDN4-high cells MPR / NMPR | 210 / 816 | **counts only** — do not treat as replicates |
| T/NK cells | 11,083 | CD4+T + CD8+T + NK |
| Potential ligands (pooled) | **147** | ≥10% CLDN4-high + T/NK receptor in NicheNet-v2 LR |
| Potential ligands MPR / NMPR | 155 / 137 | MPR call is noisier (210 sender cells) |
| Background genes | 276 | panel ∩ prior ∩ T/NK-expressed (**not** full transcriptome) |
| IFN / cytotoxicity genes in prior | 16 / 14 | all a priori members present |

Per-patient malignant and T/NK scores: `results/sample_metrics.tsv`. P06 has 15 malignant / 9 CLDN4-high cells and was **kept**.

### E1 / E2 — signed T/NK programs (patient-level)

| Test | n | Result |
| --- | --- | --- |
| NMPR vs MPR, T/NK cytotoxicity | 8 vs 4 | mean 0.786 vs 0.745; median 0.707 vs 0.804; U=16; **p=1.0** |
| NMPR vs MPR, T/NK IFN | 8 vs 4 | mean 0.381 vs 0.245; median 0.325 vs 0.227; U=23; **p=0.28** |
| NMPR vs MPR, IFN − cytotoxicity | 8 vs 4 | −0.406 vs −0.500; **p=0.37** |
| Spearman malignant CLDN4 vs cytotoxicity | 12 | **ρ=−0.112, p=0.73** |
| Spearman malignant CLDN4 vs IFN | 12 | **ρ=+0.203, p=0.53** |
| Spearman CLDN4 vs T/NK fraction | 12 | ρ=−0.09, p=0.78 |

NMPR malignant CLDN4 is higher in direction (1.47 vs 1.18, p=0.21) — same NS direction as other A3 folders. That does **not** carry a cytotoxicity-down T/NK program. IFN is higher in NMPR in direction, not lower.

Gene-wise patient Wilcoxon (p<0.15, unadjusted): NMPR-up includes **IFIT1, OAS1, MX1, ISG15, IFI6, IFI44L, ISG20, MX2, OAS3** (IFN/ISG, **same sign as the IFN score**). NMPR-down includes **BTLA, CD160, TOX2, FGFBP2**. GZMB is higher in NMPR (1.12 vs 0.68, p=0.21). IFNG p=0.81. These p-values are not FDR<0.05.

### E3 — ligand activity (unsigned NicheNet-v2 prior)

Primary sender = CLDN4-high malignant (CLDN4 only); receiver = all post-tx T/NK. Rank by Pearson of prior target scores vs gene-set membership on 276 background genes.

**A priori IFN (16 genes). Top 8 of 147:**

| Rank | Ligand | Pearson | AUROC | pearson p |
| --- | --- | --- | --- | --- |
| 1 | **IFITM1** | 0.533 | 0.941 | 1.1e−21 |
| 2 | HLA-F | 0.444 | 0.894 | 9.8e−15 |
| 3 | VSIG10 | 0.443 | 0.895 | 1.0e−14 |
| 4 | CD40 | 0.430 | 0.910 | 7.5e−14 |
| 5 | HLA-DRB5 | 0.404 | 0.884 | 3.0e−12 |
| 6 | HLA-B | 0.400 | 0.880 | 4.7e−12 |
| 7 | MDK | 0.350 | 0.851 | 2.3e−09 |
| 8 | SDC2 | 0.322 | 0.844 | 4.4e−08 |

IFITM1 / MHC class I ligands sit next to IFN genes in the published prior. AUROC 0.94 is prior recovery of an IFN list, **not** a signed CLDN4→IFN axis in GSE207422.

**A priori cytotoxicity (14 genes). Top 8:** HLA-DRA 0.295; TYROBP 0.294; HLA-DQB1 0.256; CCL4 0.237; CD47 0.215; CD58 0.211; CD274 0.201; HLA-A 0.189. These Pearson values do **not** mean cytotoxicity is repressed. CD274 AUROC on cytotoxicity is 0.64 (best AUROC in that list) with Pearson 0.20.

Patient-level ligand means in CLDN4-high cells (Wilcoxon, unadjusted, 147 ligands): IL7, ALCAM, HMGB1, AGRN are **higher in MPR** (p=0.013–0.048). NECTIN2 is higher in NMPR (p=0.048). None survive a 147-test view.

### E4 — combinatorial MPR vs NMPR receivers

Ligand **activity** for a shared ligand is the same number in MPR and NMPR settings (prior × gene set do not change). Shared potential ligands: **120**.

- NMPR-only potential (17): ADGRE5, CD38, CEACAM5, CLCF1, COL17A1, COL8A1, EDN1, EFNB2, HP, ICAM4, ITGB1, ITGB2, PTDSS1, **SPP1**, TFPI, TIMP3, ULBP2.
- MPR-only potential (35): includes IL15, IL7, IL6, CD40, CD274, PDCD1LG2 — MPR senders are 210 cells, so the 10% rule is still noisy.

This is **not** evidence that NMPR T/NK are a different NicheNet receiver program.

### What this is not

- Not Hu CopyKAT malignant IDs.
- Not a full-transcriptome `nichenetr` run (no R; background is 276 panel genes).
- Not GSE253013.
- Not TACSTD2∩CLDN4 dual-high senders (that folder is `methods/scrna_nichenet/`; companion dual-high here is 674 cells and was not used).
- Not a significant CLDN4 → IFN-down or cytotoxicity-down inductive axis in this cohort.

### Files

`results/n_table.tsv`, `sample_metrics.tsv`, `patient_level_tests.tsv`, `top_ligands.tsv`, `top_ligands_primary.tsv`, `ligand_activity_all.tsv`, `fig1`–`fig4`.

---

## 中文

**结论：符号方向的说法不成立。** GSE207422 术后 **MPR 4 例（含 pCR）/ NMPR 8 例**，患者水平 T/NK 杀伤评分 NMPR 并不低（0.786 vs 0.745，**p=1.0**），IFN 也不显著更低（0.381 vs 0.245，**p=0.28**，方向是 NMPR 更高）。恶性 CLDN4 对杀伤 **ρ=−0.11**、对 IFN **ρ=+0.20**。

先验能诚实报告的是：CLDN4 高恶性细胞 147 个潜在配体里，IFN 基因集第一名是 **IFITM1**（Pearson 0.533，AUROC 0.94）。这是先验结构（ISG/MHC 配体靠近 IFN 基因），**不是**本队列里 CLDN4 高细胞诱导或抑制 IFN 的证据。杀伤基因集前几名是 HLA-DRA、TYROBP、HLA-DQB1、CCL4、CD47、CD58、CD274、HLA-A。MPR/NMPR 组合分析里，**共有配体的活性分数相同**（120 个）；差别只在谁被 10% 规则叫成“表达”（NMPR 独有 17 个，含 SPP1；MPR 独有 35 个，含 IL15/CD40/CD274，MPR 发送端 210 个细胞，不稳定）。

发送端只用 CLDN4（中位数 1.39；1,026 / 2,051 个恶性细胞）。TACSTD2∩CLDN4 双高伴随计数 674，**未用于排序**。GSE253013 未用。推断单位是患者，不是细胞。
