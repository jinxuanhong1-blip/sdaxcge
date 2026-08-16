# B2 rework · CCLE NSCLC **protein** TACSTD2–CLDN4

**Honest verdict: n=118 and ρ=0.69 are not the same cohort.**

- **n=118** = CCLE **RPPA** 20180123 lung lines that are not SCLC (117 NSCLC + 1 carcinoid). RPPA has **no TACSTD2/TROP2 and no CLDN4** antibody. This n cannot be the correlation sample.
- **ρ=0.69** = Gygi/Nusinow CCLE **mass-spec** proteomics, lung lines with **both** proteins quantified: Spearman **0.693 (n=45)**. Pearson r = 0.687. Bootstrap 95% CI [0.48, 0.82].
- NSCLC-only MS (closer to the user’s histology label) is **n=36, ρ=0.73** — not 0.69 and not n=118.
- The earlier DepMap **RNA** lung result (ρ=0.61, n=214) was the wrong modality for this claim.

We did not tune filters to hit 0.69.

---

## 中文

### 用户更正之后还对不对？

用户说 0.69 是 **CCLE NSCLC 蛋白、n=118**，不是 DepMap RNA。公开蛋白矩阵核对后：

| 声称 | 实际对得上的公开对象 | 能否计算 TACSTD2–CLDN4 |
|---|---|---|
| n=118 | RPPA 20180123：`_LUNG` 且非 SCLC = **118**（NSCLC 117 + 类癌 1） | **不能**。214 个抗体里没有 TROP2/TACSTD2，也没有 CLDN4（只有 Claudin-7） |
| ρ=0.69 | Gygi MS（Nusinow et al., Cell 2020）肺系、两蛋白都有定量 | **能**。n=**45**，Spearman **0.693**，四舍五入就是 0.69 |

**同一张表里不存在“n=118 且 ρ=0.69”。** 把 RPPA 的样本量贴到 MS 的相关系数上，会得到用户那组数字，但那不是一次完整病例相关。

### 各蛋白数据集

**1. Gygi / Nusinow CCLE TMT-MS（主蛋白相关）**

- 文件：`protein_quant_current_normalized.csv.gz`（Gygi lab，2024-07-15 镜像）
- 两基因都在：`TACSTD2` = `sp\|P09758\|TACD2_HUMAN`，`CLDN4` = `sp\|O14493\|CLD4_HUMAN`
- TACSTD2 在全部 378 列都有值；CLDN4 只有 198 列。肺列 77 个，两蛋白齐全 **45**。

| 队列 | n | Spearman | 四舍五入 | vs 0.69 |
|---|---:|---:|---:|---|
| `_LUNG` 完整病例（主） | **45** | **0.693** | **0.69** | 命中点估计 |
| Table S1 Tissue=Lung | 45 | 0.693 | 0.69 | 同上 |
| 22Q2 lineage=lung | 46 | 0.696 | 0.70 | 邻近 |
| NSCLC only | 36 | 0.731 | 0.73 | 不匹配 |
| lung 非 SCLC | 37 | 0.738 | 0.74 | 不匹配 |
| SCLC | 9 | 0.25 | 0.25 | 不匹配 |
| 全谱系完整病例 | 198 | 0.612 | 0.61 | 不匹配 |

图：`fig_gygi_ms_lung_complete.png`。45 系名单：`gygi_ms_lung_complete_n45.csv`。

**2. CCLE RPPA 20180123（n=118 的来源）**

- 899 个细胞系 × 214 抗体。`_LUNG` = 164。
- 与 22Q2 `sample_info` 对齐：NSCLC **117**，SCLC 44，类癌 1，未注释 2。
- **lung 且非 SCLC = 118** ← 与用户 n 完全一致。
- 抗体表无 TACSTD2/TROP2/CLDN4。**不能**在此 n 上算这对蛋白的 ρ。

**3. ProCan-DepMapSanger（Gonçalves 2022，949 系）**

- 8498 蛋白矩阵含 TACSTD2（`P09758;TACD2_HUMAN`），**不含 CLDN4**（有 CLDN1/3/7）。
- 注释 NSCLC = 101，Lung = 187。都不是 118，且算不了这对相关。

### 和上一轮 RNA 的关系

上一轮 24Q4 肺系 RNA ρ=0.61（n=214）仍然正确，只是**答错了题**。用户数字是蛋白，不是转录本。蛋白上肺完整病例确实约 0.69，但样本量是 45 不是 118。

### 复现

```bash
pip install -r scripts/rework/B2_protein_followup/requirements.txt
python3 scripts/rework/B2_protein_followup/download.py
python3 scripts/rework/B2_protein_followup/analyze.py
```

---

## English

### TL;DR

The user-claimed **ρ=0.69, n=118, CCLE NSCLC protein** splits across two public files:

1. **n=118** is the **RPPA** lung-not-SCLC headcount. RPPA cannot measure this gene pair.
2. **ρ=0.69** is **Gygi MS** lung complete cases, **n=45** (Spearman 0.693). NSCLC-only MS is n=36 / ρ=0.73.

Cite **0.69 (n=45, Gygi CCLE MS, lung lines with both proteins)** if the claim is the correlation. Do not cite n=118 as the correlation sample.

### Methods (pre-specified)

- Download Gygi normalized TMT, CCLE RPPA + antibody table, ProCan 8498 averaged matrix.
- Spearman + Pearson on complete cases. No imputation. No expression filter.
- NSCLC / SCLC labels from DepMap 22Q2 `sample_info` (`lineage_subtype`) and Gygi Table S1 tissue.
- Match rules: n exactly 118; ρ rounds to 0.69 at 2 d.p.

### Limits

- Gygi lung n=45 still mixes some SCLC (DMS114, SHP77, NCI-H146, …). Restricting to NSCLC moves ρ to 0.73 and n to 36.
- TMT ratios are relative, not IHC H-scores. Not ADC-ready protein quantification.
- RPPA500 / TCPA CCLE matrix was not used; the 2018 CCLE RPPA file is the one whose lung-not-SCLC count is exactly 118.
- ProCan lacks CLDN4, so a 949-line MS recompute of this pair is impossible from the public 8498 matrix.
