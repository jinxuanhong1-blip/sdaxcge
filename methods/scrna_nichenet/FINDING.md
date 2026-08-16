# FINDING — NicheNet-style ligands from TACSTD2-high / CLDN4-high malignant cells

**Verdict: the signed claim is not supported.** On public GSE207422 (n=4 MPR including pCR, n=8 NMPR), T/NK cytotoxicity is **not** lower in NMPR (mean 0.786 vs 0.745, U=16, **p=1.0**) and exhaustion is **not** higher in a detectable way (0.372 vs 0.325, U=19, **p=0.68**). Malignant TACSTD2/CLDN4 vs T/NK cytotoxicity is **ρ=+0.13, p=0.70** (wrong sign); vs exhaustion **ρ=+0.03, p=0.91**.

What the prior *does* return, honestly: among 146 potential ligands from high-barrier malignant cells, **CD274** is the top ligand for the a priori exhaustion set (Pearson 0.129, AUROC 0.58, p=0.040). That is a weak prior hit. Combinatorial MPR vs NMPR **activity ranks of shared ligands are identical** (same prior × same gene set). The combinatorial difference is the **expressed ligand set**, not a new ranking.

GSE253013 was **not used** (9.3 GB; no MPR/NMPR).

---

## English

### Data and n

| Item | n | Note |
| --- | --- | --- |
| Post-tx patients | 12 | Hu et al. GSE207422; pCR=P06 counted as MPR |
| MPR / NMPR patients | **4 / 8** | unit of inference |
| DRMref-annotated cells | 30,877 | Liu et al. *NAR* 2024; barcode-matched to GEO |
| Malignant cells | 2,051 | not Hu CopyKAT |
| High-barrier malignant | **674** | TACSTD2 **and** CLDN4 `log1p(CP10k)` ≥ malignant median (1.50 / 1.39) |
| High-barrier cells MPR / NMPR | 76 / 598 | **counts only** — do not treat as replicates |
| T/NK cells | 11,083 | CD4+T + CD8+T + NK |
| Potential ligands (pooled) | **146** | ≥10% high-barrier + T/NK receptor in NicheNet-v2 LR |
| Potential ligands MPR / NMPR | 167 / 140 | MPR call is unstable (76 sender cells) |
| Background genes | 256 | panel ∩ prior ∩ T/NK-expressed (**not** full transcriptome) |

Per-patient malignant and T/NK scores: `results/sample_metrics.tsv`. P06 has 15 malignant / 8 high-barrier cells and was **kept**.

### E1 / E2 — signed T/NK programs (patient-level)

| Test | n | Result |
| --- | --- | --- |
| NMPR vs MPR, T/NK cytotoxicity | 8 vs 4 | mean 0.786 vs 0.745; median 0.707 vs 0.804; U=16; **p=1.0** |
| NMPR vs MPR, T/NK exhaustion | 8 vs 4 | mean 0.372 vs 0.325; U=19; **p=0.68** |
| NMPR vs MPR, exhaustion − cytotoxicity | 8 vs 4 | −0.415 vs −0.421; **p=0.81** |
| Spearman malignant barrier vs cytotoxicity | 12 | **ρ=+0.126, p=0.70** |
| Spearman malignant barrier vs exhaustion | 12 | **ρ=+0.035, p=0.91** |
| Spearman barrier vs T/NK fraction | 12 | ρ=−0.37, p=0.24 |

NMPR malignant barrier is higher in direction (1.52 vs 1.15, p=0.11) — same NS direction as other A3 folders. That does **not** carry a cytotoxicity-down / exhaustion-up T/NK program.

Gene-wise patient Wilcoxon (p<0.15, unadjusted): NMPR-up = AIMP1, ISG15, JMJD6, MX1, NAMPT (IFN/stress, **not** the exhaustion list). NMPR-down includes **BTLA, CD160, TOX2, FGFBP2** — several exhaustion/cytotoxicity genes go **down**, opposite the claim. GZMB is higher in NMPR (1.12 vs 0.68, p=0.21). PDCD1 p=1.0.

### E3 — ligand activity (unsigned NicheNet-v2 prior)

Primary sender = TACSTD2-high ∩ CLDN4-high malignant; receiver = all post-tx T/NK. Rank by Pearson of prior target scores vs gene-set membership on 256 background genes.

**A priori exhaustion (14 genes). Top 8 of 146:**

| Rank | Ligand | Pearson | AUROC | pearson p |
| --- | --- | --- | --- | --- |
| 1 | **CD274** | 0.129 | 0.580 | 0.040 |
| 2 | ALCAM | 0.053 | 0.463 | 0.40 |
| 3 | TNFSF9 | 0.051 | 0.481 | 0.41 |
| 4 | COL6A1 | 0.035 | 0.375 | 0.58 |
| 5 | PDGFA | 0.025 | 0.414 | 0.69 |
| 6 | CALR | 0.023 | 0.383 | 0.71 |
| 7 | HLA-DRB1 | 0.015 | 0.390 | 0.81 |
| 8 | CD40 | 0.015 | 0.459 | 0.81 |

Only CD274 is even weakly above noise. AUROC 0.58 is not a strong recovery of the exhaustion list.

**A priori cytotoxicity (14 genes). Top 8:** HLA-DRA 0.300; TYROBP 0.296; CCL4 0.280; HLA-DQB1 0.261; CD47 0.218; CD58 0.214; CD274 0.212; HLA-A 0.201. These Pearson values are larger because MHC/immune ligands sit near immune genes in the prior — they do **not** mean cytotoxicity is repressed.

TGFB1, LGALS9, HLA-E, CD274 expression in high-barrier cells does **not** differ NMPR vs MPR (all p≥0.57). SPP1 is higher in NMPR (0.40 vs 0.02, p=0.084, unadjusted). ALCAM and HMGB1 are **higher in MPR** (p=0.048, unadjusted; 146 ligands tested).

### E4 — combinatorial MPR vs NMPR receivers

Ligand **activity** for a shared ligand is the same number in MPR and NMPR settings (prior × gene set do not change). Shared potential ligands: **129**.  

- NMPR-only potential (11): ADGRE5, CD38, CLCF1, EDN1, EFNB2, HP, ICAM4, ITGB1, PTDSS1, **SPP1**, TFPI.  
- MPR-only potential (38): includes IL15, IL7, CXCL10, CD40, MICB — but MPR senders are only 76 cells, so the 10% rule is noisy.

This is **not** evidence that NMPR T/NK are a different NicheNet receiver program.

### What this is not

- Not Hu CopyKAT malignant IDs.  
- Not a full-transcriptome `nichenetr` run (no R; background is 256 panel genes).  
- Not GSE253013.  
- Not a significant CD274 / TGFB / IL10 inductive axis onto exhausted T cells in this cohort.

### Files

`results/n_table.tsv`, `sample_metrics.tsv`, `patient_level_tests.tsv`, `top_ligands.tsv`, `ligand_activity_all.tsv`, `fig1`–`fig4`.

---

## 中文

**结论：符号方向的说法不成立。** GSE207422 术后 **MPR 4 例（含 pCR）/ NMPR 8 例**，患者水平 T/NK 杀伤评分 NMPR 并不低（0.786 vs 0.745，**p=1.0**），耗竭也不显著更高（0.372 vs 0.325，**p=0.68**）。恶性 TACSTD2/CLDN4 对杀伤 **ρ=+0.13**、对耗竭 **ρ=+0.03**。

先验能诚实报告的是：高屏障恶性细胞 146 个潜在配体里，耗竭基因集第一名是 **CD274**（Pearson 0.129，AUROC 0.58）。很弱。MPR/NMPR 组合分析里，**共有配体的活性分数相同**；差别只在谁被 10% 规则叫成“表达”（NMPR 独有 11 个，含 SPP1；MPR 独有 38 个，但 MPR 高屏障只有 76 个细胞，不稳定）。

GSE253013 未用（9.3 GB，无 MPR/NMPR）。推断单位是患者，不是细胞。
