# STK11 / KEAP1 / KRAS genotype vs. TACSTD2 (TROP2) & CLDN4 (Claudin-4) in lung ICI cohorts

> Parallel slice `fable_stk11`. Every output of this slice lives under
> `notes/fable_stk11/`, `scripts/fable_stk11/`, and `results/fable_stk11/`.

---

## English

### 1. Question
Antibody–drug conjugates (ADCs) against **TACSTD2 (TROP2)** (e.g. sacituzumab
govitecan, datopotamab deruxtecan) and **CLDN4 (Claudin-4)** are being developed
for non-small-cell lung cancer (NSCLC), often positioned as options *after*
immune-checkpoint-inhibitor (ICI) failure. **STK11 (LKB1)** and **KEAP1**
mutations — especially co-occurring with **KRAS** — mark an immune-"cold",
ICI-resistant subset of lung adenocarcinoma (LUAD).

This slice asks, using only public data:
1. **Genotype ↔ expression** — do STK11/KEAP1/KRAS mutations change TACSTD2/CLDN4 mRNA?
2. **Genotype ↔ response** — do these mutations change ICI benefit?
3. **Expression ↔ response** — does TACSTD2/CLDN4 expression itself predict ICI benefit?

Answering all three lets us reason about whether TROP2/Claudin-4 ADCs are a
rational salvage strategy in the STK11/KEAP1-mutant ICI-resistant subset.

### 2. Public datasets (no controlled access)
| Cohort | Source | n | Data used | Axis |
|---|---|---|---|---|
| TCGA-LUAD (PanCancer Atlas) | cBioPortal `luad_tcga_pan_can_atlas_2018` | 566 seq / 510 RNA | mutation + RNA-seq | genotype↔expression |
| TCGA-LUSC (PanCancer Atlas) | cBioPortal `lusc_tcga_pan_can_atlas_2018` | 484 | mutation + RNA-seq | secondary check |
| Rizvi 2015 (Science) | cBioPortal `luad_mskcc_2015` | 35 | mutation + DCB | genotype↔response |
| Hellmann 2018 (Cancer Cell) | cBioPortal `nsclc_mskcc_2018` | 75 | mutation + DCB + PFS | genotype↔response |
| Rizvi 2018 (JCO) | cBioPortal `nsclc_pd1_msk_2018` | 240 | mutation + DCB + PFS | genotype↔response |
| Jung 2019 (Nat Commun) | GEO `GSE135222` | 27 | RNA-seq (TPM) + PFS | expression↔response |

All raw pulls are gene-subset (9 genes) or a single 1.6 MB GEO matrix, so the
entire `results/fable_stk11/data/` footprint is **~2.4 MB** (limit: 2 GB).

### 3. Methods (see `scripts/fable_stk11/`)
- **Mutation status**: a sample is "mutant" for a gene if it carries any
  non-synonymous / functional variant (missense, nonsense, frameshift, splice,
  in-frame indel) in that gene. The wild-type universe is the cohort's
  *sequenced* sample list, so absence of a call means true wild-type, not missing.
- **Expression**: cBioPortal RNA-seq V2 RSEM and GEO TPM, both `log2(x+1)`.
- **Genotype↔expression**: Mann-Whitney U + Welch t-test, medians, log2
  fold-change, Cliff's delta, Cohen's d; Benjamini-Hochberg FDR within cohort.
  A KRAS-restricted co-mutation analysis isolates the KRAS/STK11 ("KL") and
  KRAS/KEAP1 subsets.
- **Genotype↔response**: durable clinical benefit (DCB, harmonized YES/NO)
  by Fisher exact; PFS by Kaplan-Meier + log-rank + univariable Cox HR
  (`lifelines`); pooled effect via Mantel-Haenszel odds ratio.
- **Expression↔response**: Mann-Whitney (DCB proxy = PFS ≥ 6 months),
  univariable Cox (HR per log2 unit), median-split KM.
- **Verification** (`05_verify.py`): 28 automated checks — sample counts,
  mutation-frequency plausibility, independent recomputation of headline
  p-values, determinism, response-direction sanity, and the < 2 GB footprint.

### 4. Results

**Axis 1 — Genotype ↔ expression (TCGA-LUAD, n = 510).**
STK11 mutation is associated with **lower** ADC-target expression:
- TACSTD2: median log2 12.20 (mut) vs 12.72 (WT), log2FC −0.51, MWU p = 3.4×10⁻⁵,
  FDR = 1.0×10⁻⁴, Cliff's δ = −0.30.
- CLDN4: 12.75 vs 13.21, log2FC −0.47, MWU p = 7.7×10⁻⁶, FDR ≈ 3×10⁻⁵, δ = −0.33.

KEAP1 mutation lowers CLDN4 (FDR 0.027) but not TACSTD2. KRAS mutation *alone*
shows no significant effect. The effect concentrates in the **KRAS/STK11 "KL"**
subset: within KRAS-mutant tumours, STK11 co-mutation lowers TACSTD2 (log2FC
−0.71, p = 1.0×10⁻⁴, δ = −0.43) and CLDN4 (−0.68, p = 5.7×10⁻⁵, δ = −0.44).
Both targets are weakly *negatively* correlated with CD8A (TACSTD2 ρ = −0.10,
CLDN4 ρ = −0.12; p < 0.05). TCGA-LUSC carries too few STK11 (n = 5) or KRAS
(n = 7) mutants to be informative — reported for completeness only.

**Axis 2 — Genotype ↔ ICI response (MSK cohorts).**
Whole-cohort genotype effects on DCB are modest and cohort-dependent
(STK11 pooled Mantel-Haenszel OR 0.69; KEAP1 1.53; KRAS 1.41). Restricting to
**KRAS-mutant** tumours recovers the expected Skoulidis-2018 biology: STK11
co-mutation reduces durable benefit in both scorable cohorts (Hellmann 2018 DCB
29% vs 63%, OR 0.24, p = 0.19; Rizvi 2018 20% vs 40%, OR 0.38, p = 0.089;
**pooled MH OR 0.36**), whereas KEAP1 co-mutation is null (pooled OR 0.96). For
PFS, STK11 mutation trends toward worse outcome in Hellmann 2018 (median 2.4 vs
6.8 months, Cox HR 1.95, 95% CI 0.91–4.17, p = 0.087). KRAS alone in the small
Rizvi 2015 cohort associates with *higher* DCB (OR 16, p = 0.011), consistent
with KRAS-mutant tumours being smoking-related and TMB-high.

**Axis 3 — Expression ↔ ICI response (GSE135222, n = 27, 21 events).**
TACSTD2 and CLDN4 expression do **not** predict ICI benefit: Mann-Whitney by DCB
proxy p = 0.61 / 0.69; Cox HR per log2 unit 1.03 (p = 0.78) / 1.05 (p = 0.61);
median-split log-rank p = 0.17 / 0.91. As a positive control, the immune markers
behave as expected (CD8A Cox HR 0.77, CD274 HR 0.69; both protective direction),
confirming the survival coding is correct. The cohort is small, so this is an
exploratory null.

### 5. Synthesis
The STK11-mutant (particularly KRAS/STK11 "KL") LUAD subset is simultaneously
(a) **ICI-resistant** and (b) carries **lower** TROP2 and Claudin-4 mRNA. TROP2
and Claudin-4 expression are not themselves ICI-response biomarkers. The
practical implication: TROP2/Claudin-4 ADCs are a plausible salvage idea for
ICI-refractory NSCLC in general, but the specific STK11/KL subset that most needs
a non-immune option tends to have *lower* target antigen — so target-expression
screening (IHC/RNA) is advisable rather than assuming the ICI-resistant genotype
is ADC-target–enriched.

### 6. Limitations
- MSK ICI cohorts are targeted panels without matched tumour RNA, so
  expression↔response is only tested in the small GSE135222 RNA cohort.
- GSE135222 has no matched STK11/KEAP1 genotype, so the full triangle is closed
  across (not within) cohorts.
- RSEM vs TPM units differ across cohorts; comparisons are always within-cohort.
- Observational public data; no correction for PD-L1, TMB, line of therapy, or
  regimen (mono- vs combination ICI) beyond what each cohort provides.

### 7. Reproduce
```bash
pip install -r scripts/fable_stk11/requirements.txt
bash scripts/fable_stk11/run_all.sh
```
Outputs: tables in `results/fable_stk11/tables/`, figures in
`results/fable_stk11/figures/`, processed data in
`results/fable_stk11/data/processed/`, verification in
`notes/fable_stk11/verification.md`.

---

## 中文

### 1. 问题
针对 **TACSTD2（TROP2）**（如 sacituzumab govitecan、datopotamab deruxtecan）和
**CLDN4（Claudin-4）** 的抗体偶联药物（ADC）正在非小细胞肺癌（NSCLC）中开发，
通常定位为免疫检查点抑制剂（ICI）失败**之后**的选择。**STK11（LKB1）** 与
**KEAP1** 突变——尤其是与 **KRAS** 共突变时——标记了肺腺癌（LUAD）中一个免疫
"冷"、对 ICI 耐药的亚群。

本切片仅使用公共数据回答：
1. **基因型 ↔ 表达**：STK11/KEAP1/KRAS 突变是否改变 TACSTD2/CLDN4 的 mRNA 水平？
2. **基因型 ↔ 疗效**：这些突变是否改变 ICI 获益？
3. **表达 ↔ 疗效**：TACSTD2/CLDN4 表达本身能否预测 ICI 获益？

三条轴合起来，可推断 TROP2/Claudin-4 ADC 在 STK11/KEAP1 突变、ICI 耐药亚群中
是否是合理的挽救策略。

### 2. 公共数据集（均无需受控访问）
| 队列 | 来源 | n | 使用数据 | 轴 |
|---|---|---|---|---|
| TCGA-LUAD | cBioPortal `luad_tcga_pan_can_atlas_2018` | 566 测序 / 510 RNA | 突变 + RNA-seq | 基因型↔表达 |
| TCGA-LUSC | cBioPortal `lusc_tcga_pan_can_atlas_2018` | 484 | 突变 + RNA-seq | 次要验证 |
| Rizvi 2015 (Science) | cBioPortal `luad_mskcc_2015` | 35 | 突变 + DCB | 基因型↔疗效 |
| Hellmann 2018 (Cancer Cell) | cBioPortal `nsclc_mskcc_2018` | 75 | 突变 + DCB + PFS | 基因型↔疗效 |
| Rizvi 2018 (JCO) | cBioPortal `nsclc_pd1_msk_2018` | 240 | 突变 + DCB + PFS | 基因型↔疗效 |
| Jung 2019 (Nat Commun) | GEO `GSE135222` | 27 | RNA-seq (TPM) + PFS | 表达↔疗效 |

所有原始下载仅取 9 个基因子集或单个 1.6 MB 的 GEO 矩阵，因此整个
`results/fable_stk11/data/` 体积约 **2.4 MB**（上限 2 GB）。

### 3. 方法（见 `scripts/fable_stk11/`）
- **突变状态**：若样本在某基因上携带任一非同义/功能性变异（错义、无义、移码、
  剪接、框内插入缺失），则记为该基因"突变"。野生型全集取各队列的**测序**样本列表，
  因此"无变异记录"意味着真正的野生型而非缺失。
- **表达**：cBioPortal RNA-seq V2 RSEM 与 GEO TPM，均取 `log2(x+1)`。
- **基因型↔表达**：Mann-Whitney U + Welch t 检验、中位数、log2 倍数变化、
  Cliff's delta、Cohen's d；队列内 Benjamini-Hochberg FDR。KRAS 限定的共突变分析
  单独刻画 KRAS/STK11（"KL"）与 KRAS/KEAP1 亚群。
- **基因型↔疗效**：持久临床获益（DCB，统一为 YES/NO）用 Fisher 精确检验；PFS 用
  Kaplan-Meier + log-rank + 单因素 Cox HR（`lifelines`）；用 Mantel-Haenszel 合并 OR。
- **表达↔疗效**：Mann-Whitney（DCB 代理 = PFS ≥ 6 个月）、单因素 Cox（每 log2 单位 HR）、
  按中位数分组的 KM。
- **验证**（`05_verify.py`）：28 项自动检查——样本数、突变频率合理性、头条 p 值的
  独立重算、确定性、疗效方向合理性，以及 < 2 GB 体积。

### 4. 结果

**轴 1 — 基因型 ↔ 表达（TCGA-LUAD，n = 510）。**
STK11 突变与 ADC 靶点**更低**的表达相关：
- TACSTD2：中位 log2 12.20（突变）vs 12.72（野生型），log2FC −0.51，
  MWU p = 3.4×10⁻⁵，FDR = 1.0×10⁻⁴，Cliff's δ = −0.30。
- CLDN4：12.75 vs 13.21，log2FC −0.47，MWU p = 7.7×10⁻⁶，FDR ≈ 3×10⁻⁵，δ = −0.33。

KEAP1 突变降低 CLDN4（FDR 0.027）但不影响 TACSTD2。**单独** KRAS 突变无显著效应。
效应集中于 **KRAS/STK11（"KL"）** 亚群：在 KRAS 突变肿瘤中，STK11 共突变使
TACSTD2（log2FC −0.71，p = 1.0×10⁻⁴，δ = −0.43）与 CLDN4（−0.68，p = 5.7×10⁻⁵，
δ = −0.44）进一步降低。两个靶点均与 CD8A 呈弱**负**相关（TACSTD2 ρ = −0.10，
CLDN4 ρ = −0.12；p < 0.05）。TCGA-LUSC 中 STK11（n = 5）与 KRAS（n = 7）突变太少，
不具信息量，仅列出以求完整。

**轴 2 — 基因型 ↔ ICI 疗效（MSK 队列）。**
全队列层面基因型对 DCB 的效应较弱且依赖队列（STK11 合并 Mantel-Haenszel OR 0.69；
KEAP1 1.53；KRAS 1.41）。**限定 KRAS 突变**后恢复了预期的 Skoulidis-2018 生物学：
STK11 共突变在两个可评分队列中均降低持久获益（Hellmann 2018 DCB 29% vs 63%，
OR 0.24，p = 0.19；Rizvi 2018 20% vs 40%，OR 0.38，p = 0.089；**合并 MH OR 0.36**），
而 KEAP1 共突变为零效应（合并 OR 0.96）。PFS 方面，STK11 突变在 Hellmann 2018 中
趋向更差（中位 2.4 vs 6.8 个月，Cox HR 1.95，95% CI 0.91–4.17，p = 0.087）。在小样本
Rizvi 2015 中，单独 KRAS 与**更高** DCB 相关（OR 16，p = 0.011），符合 KRAS 突变肿瘤
多为吸烟相关、TMB 高的特点。

**轴 3 — 表达 ↔ ICI 疗效（GSE135222，n = 27，21 事件）。**
TACSTD2 与 CLDN4 表达**不**预测 ICI 获益：按 DCB 代理的 Mann-Whitney
p = 0.61 / 0.69；每 log2 单位 Cox HR 1.03（p = 0.78）/ 1.05（p = 0.61）；
中位数分组 log-rank p = 0.17 / 0.91。作为阳性对照，免疫标志物方向符合预期
（CD8A Cox HR 0.77，CD274 HR 0.69，均为保护方向），证实生存编码正确。样本量小，
故此为探索性阴性结果。

### 5. 综合结论
STK11 突变（尤其 KRAS/STK11 "KL"）LUAD 亚群同时具备：（a）**对 ICI 耐药**，且
（b）TROP2 与 Claudin-4 mRNA **更低**。TROP2/Claudin-4 表达本身并非 ICI 疗效标志物。
实践含义：TROP2/Claudin-4 ADC 作为 ICI 难治性 NSCLC 的挽救思路总体合理，但最需要
非免疫方案的 STK11/KL 亚群，其靶抗原反而**更低**——因此应先做靶点表达筛查（IHC/RNA），
而不能默认 ICI 耐药基因型富集 ADC 靶点。

### 6. 局限
- MSK ICI 队列为靶向 panel，无配对肿瘤 RNA，故表达↔疗效仅在小样本 GSE135222 中检验。
- GSE135222 无配对 STK11/KEAP1 基因型，故完整"三角"是跨队列（而非同队列内）闭合。
- 各队列 RSEM 与 TPM 单位不同，比较一律在队列内进行。
- 观察性公共数据，除各队列已提供者外，未对 PD-L1、TMB、治疗线数或方案（单药 vs 联合 ICI）作校正。

### 7. 复现
```bash
pip install -r scripts/fable_stk11/requirements.txt
bash scripts/fable_stk11/run_all.sh
```
输出：表格见 `results/fable_stk11/tables/`，图见 `results/fable_stk11/figures/`，
处理后数据见 `results/fable_stk11/data/processed/`，验证见
`notes/fable_stk11/verification.md`。
