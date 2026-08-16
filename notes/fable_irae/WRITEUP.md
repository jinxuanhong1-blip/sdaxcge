# TACSTD2 (TROP2) & CLDN4 (claudin-4) in public ICI irAE datasets
# 公共 ICI 免疫相关不良事件（irAE）数据集中的 TACSTD2（TROP2）与 CLDN4（claudin-4）

*Parallel slice `fable_irae`. All outputs live under `notes/fable_irae/`,
`scripts/fable_irae/`, `results/fable_irae/`.*

---

## English

### 1. Objective
Quantify the two epithelial surface/tight-junction markers **TACSTD2 (TROP2)** and
**CLDN4 (claudin-4)** in **public** immune-checkpoint-inhibitor (ICI) datasets that carry
immune-related adverse event (irAE) annotations relevant to lung ICI therapy
(pneumonitis, colitis, and a systemic blood comparison), and ask two questions:
1. **Are the genes measurable** in each compartment?
2. **Do they track irAE status** (case vs control / irAE Yes vs No)?

Because TACSTD2/CLDN4 are epithelial genes, "measurable" is compartment-dependent, so we
deliberately sampled an epithelial→immune→systemic gradient.

### 2. Datasets (all from NCBI GEO, human, processed data only)

| Accession | Compartment | irAE | n | Processed matrix |
|-----------|-------------|------|---|------------------|
| [GSE206300](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206300) | Colon mucosa **epithelial** single-nucleus RNA-seq | **irColitis** | 26 donors (12 Case / 14 Control); 81,707 nuclei | raw UMI counts (`.h5ad`) |
| [GSE277136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277136) | Bronchoalveolar lavage fluid scRNA-seq (**lung**) | **ICI pneumonitis** | 4 pneumonitis (AE) + 3 healthy (HC); 74,607 cells | log1p-normalized (`.raw`) |
| [GSE319496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319496) | Whole-blood **bulk** RNA-seq (mRCC, nivolumab+ipilimumab) | irAE Yes/No | 51 (29 Yes / 22 No) | raw counts (CSV) |

Candidate discovery (527 GEO records screened) and the selection rationale are in
[`dataset_selection.md`](dataset_selection.md). Raw matrices (~2 GB) are **not** committed;
re-download with `scripts/fable_irae/run.sh`. Only small derived tables/figures are versioned.

### 3. Methods
- **GSE206300 (primary).** Donor-level **pseudobulk**: per donor, CP10K =
  1e4 · Σ(gene counts)/Σ(all counts), then log1p; also per-nucleus detection rate.
  Case (irColitis) vs Control compared with Mann–Whitney U, confirmed by Welch t-test and a
  5,000-sample bootstrap 95 % CI for log2FC (seed 20260816). EPCAM (epithelial +control),
  PTPRC (immune −control), CLDN3/CLDN7/ACTB as references.
- **GSE319496.** CPM + log2(CPM+1); irAE Yes vs No by Mann–Whitney U. Measurability =
  presence in matrix + % samples with non-zero counts.
- **GSE277136.** Markers pulled from `.raw` (log1p, 22,126 genes). Expression localized
  across BAL cell types (`predicted_labels`); AE vs HC compared by per-sample pseudobulk
  (Mann–Whitney U). *AE vs HC is fully confounded with data source/batch — descriptive only.*

### 4. Results

**Measurability gradient (see `figures/summary_measurability.png`):**
`colon epithelium (CLDN4 robust, TACSTD2 at floor) > BALF (both at floor) > blood (CLDN4 low, TACSTD2 absent)`.

**4.1 Colon epithelium — irColitis (GSE206300).**

| Gene | mean CP10K Case | mean CP10K Ctrl | detection (Case/Ctrl) | log2FC [boot 95 % CI] | p (MWU / t) |
|------|-----------------|-----------------|-----------------------|-----------------------|-------------|
| **CLDN4** | 1.625 | 1.532 | 25.0 % / 24.2 % | +0.086 [−0.60, 0.80] | 0.817 / 0.746 |
| **TACSTD2** | 0.00166 | 0.00076 | 0.06 % / 0.02 % | +1.12 [−0.65, 2.52] | 0.588 / 0.224 |
| EPCAM (+ctrl) | — | — | 49.8 % / 49.2 % | +0.156 | 0.396 |
| PTPRC (immune) | — | — | 1.34 % / 0.41 % | +1.54 [0.52, 2.39] | **0.017** |

CLDN4 is robustly expressed in colon epithelium; TACSTD2 sits at the detection floor even
in the epithelial compartment (single-nucleus data undercount some surface-protein mRNAs).
Neither differs between irColitis and control. The positive control behaves as expected:
PTPRC rises in the irColitis "epithelial" compartment (bootstrap CI excludes 0), i.e. immune
infiltration bleeding into the epithelial fraction during active colitis.

**4.2 Whole blood — irAE Yes/No (GSE319496).**
TACSTD2 is **absent from the count matrix** (not measurable). CLDN4 is low
(mean 0.51 CPM, detected in 57 % of samples) and shows **no** irAE association
(log2FC −0.078, MWU p = 0.843).

**4.3 BALF / lung — ICI pneumonitis (GSE277136).**
Both genes exist in `.raw` but expression is at the floor: BAL is immune-dominated
(alveolar macrophages, T cells). Genuine epithelial cells are only **102/74,607 (0.14 %)**
and even those read ~0 for TACSTD2/CLDN4/EPCAM — BAL harvests luminal immune cells, not
intact epithelium. Per-sample pseudobulk AE-vs-HC is non-significant and batch-confounded
(TACSTD2 p = 0.057 at floor; CLDN4 p = 0.21). PTPRC is numerically higher in AE than HC
(mean log1p 1.73 vs 1.18) but the n=4 vs 3 contrast is not significant (MWU p = 0.40).

### 5. Verification
`scripts/fable_irae/07_verify.py` re-derives every headline number independently, adds a
parametric test + bootstrap CI, and checks controls — all pass:
- EPCAM high in epithelium ✔  · PTPRC ~0 in epithelium ✔  · PTPRC high in BALF immune ✔
- Machine-readable verdict + all stats: `results/fable_irae/tables/verification_report.json`,
  consolidated in `results/fable_irae/tables/master_summary.tsv`.

### 6. Limitations
Single-nucleus data undercount surface-protein mRNAs (depresses TACSTD2); BAL contains
essentially no epithelium; GSE319496 is metastatic renal-cell carcinoma (systemic irAE),
not lung; sample sizes are modest (26 / 51 / 7); the AE-vs-HC lung contrast is
batch-confounded. Results are hypothesis-level, not clinical.

### 7. Conclusion
Across these public ICI datasets, **neither TACSTD2 nor CLDN4 shows a significant,
reproducible association with irAE status.** Of the two, **CLDN4 is the more measurable**
(well-expressed in colon epithelium), while **TACSTD2 is largely at the detection floor**
except where genuine epithelium is captured, and is absent from the whole-blood matrix.
Practically, neither is a usable blood-based irAE biomarker here; any epithelial-marker
signal for these toxicities must be sought in target-tissue biopsies, not BAL or blood.

---

## 中文

### 1. 目的
在带有免疫相关不良事件（irAE）标注、且与肺癌 ICI 治疗相关（肺炎 pneumonitis、结肠炎
colitis，并加入一个外周血系统性对照）的**公共** ICI 数据集中，定量两个上皮标志物
**TACSTD2（TROP2）** 与 **CLDN4（claudin-4）**，回答两个问题：
1. 在各个组织/细胞区室中这两个基因**是否可测**？
2. 它们**是否与 irAE 状态相关**（病例 vs 对照 / 有无 irAE）？

由于 TACSTD2/CLDN4 是上皮基因，"可测性"取决于区室，因此我们特意覆盖了
上皮→免疫→系统血液 的梯度。

### 2. 数据集（均来自 NCBI GEO，人类，仅使用已处理数据）

| 编号 | 区室 | irAE 类型 | 样本量 | 处理后矩阵 |
|------|------|-----------|--------|------------|
| [GSE206300](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206300) | 结肠黏膜**上皮**单核 RNA-seq | **免疫性结肠炎 irColitis** | 26 例（12 病例 / 14 对照）；81,707 个核 | 原始 UMI 计数（`.h5ad`） |
| [GSE277136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277136) | 支气管肺泡灌洗液 scRNA-seq（**肺**） | **ICI 肺炎** | 4 肺炎(AE) + 3 健康(HC)；74,607 个细胞 | log1p 归一化（`.raw`） |
| [GSE319496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319496) | 全血 **bulk** RNA-seq（转移性肾癌，纳武+伊匹） | 有/无 irAE | 51（29 有 / 22 无） | 原始计数（CSV） |

候选筛选（共筛查 527 条 GEO 记录）与选择理由见
[`dataset_selection.md`](dataset_selection.md)。原始大矩阵（约 2 GB）**不纳入版本控制**，
可用 `scripts/fable_irae/run.sh` 重新下载；仓库仅保存小体积的衍生表格与图。

### 3. 方法
- **GSE206300（主分析）**：按供者做**伪 bulk（pseudobulk）**：每位供者
  CP10K = 1e4 · Σ(基因计数)/Σ(全部计数)，再取 log1p；并计算每个核的检出率。
  病例（irColitis）vs 对照用 Mann–Whitney U 检验，并以 Welch t 检验与 5,000 次自助法
  （bootstrap，种子 20260816）的 log2FC 95% 置信区间进行验证。EPCAM 为上皮阳性对照，
  PTPRC 为免疫阴性对照，另设 CLDN3/CLDN7/ACTB 参考。
- **GSE319496**：CPM + log2(CPM+1)；有无 irAE 用 Mann–Whitney U 检验。可测性 = 是否在
  矩阵中 + 非零样本比例。
- **GSE277136**：从 `.raw`（log1p，22,126 基因）提取标志物；按 BAL 细胞类型
  （`predicted_labels`）定位表达；AE vs HC 用每样本伪 bulk 的 Mann–Whitney U 比较。
  *AE 与 HC 完全与数据来源/批次混杂，仅作描述。*

### 4. 结果

**可测性梯度（见 `figures/summary_measurability.png`）：**
`结肠上皮（CLDN4 稳定可测，TACSTD2 处于检出下限）＞ BALF（两者均处于下限）＞ 血液（CLDN4 低，TACSTD2 缺失）`。

**4.1 结肠上皮 — irColitis（GSE206300）**

| 基因 | 病例 CP10K | 对照 CP10K | 检出率（病例/对照） | log2FC [自助 95%CI] | p（MWU / t） |
|------|-----------|-----------|---------------------|---------------------|--------------|
| **CLDN4** | 1.625 | 1.532 | 25.0% / 24.2% | +0.086 [−0.60, 0.80] | 0.817 / 0.746 |
| **TACSTD2** | 0.00166 | 0.00076 | 0.06% / 0.02% | +1.12 [−0.65, 2.52] | 0.588 / 0.224 |
| EPCAM（阳性对照） | — | — | 49.8% / 49.2% | +0.156 | 0.396 |
| PTPRC（免疫） | — | — | 1.34% / 0.41% | +1.54 [0.52, 2.39] | **0.017** |

CLDN4 在结肠上皮中稳定表达；TACSTD2 即使在上皮区室也处于检出下限（单核数据会低估部分
膜蛋白 mRNA）。两者在 irColitis 与对照间均无显著差异。阳性对照符合预期：PTPRC 在
irColitis 的"上皮"区室升高（自助 CI 不含 0），提示活动性结肠炎时免疫浸润渗入上皮组分。

**4.2 全血 — 有/无 irAE（GSE319496）**
TACSTD2 **未出现在计数矩阵中**（不可测）。CLDN4 表达很低（均值 0.51 CPM，57% 样本可检出），
与 irAE **无**关联（log2FC −0.078，MWU p = 0.843）。

**4.3 BALF / 肺 — ICI 肺炎（GSE277136）**
两个基因都存在于 `.raw`，但表达处于下限：BAL 以免疫细胞为主（肺泡巨噬细胞、T 细胞）。
真正的上皮细胞仅 **102/74,607（0.14%）**，且这些细胞的 TACSTD2/CLDN4/EPCAM 读数几乎为 0
——BAL 采集的是腔内免疫细胞，而非完整上皮。每样本伪 bulk 的 AE vs HC 不显著且受批次混杂
（TACSTD2 p = 0.057，处于下限；CLDN4 p = 0.21）。PTPRC 在肺炎组数值更高
（mean log1p 1.73 vs 1.18），但 n=4 vs 3 的对比不显著（MWU p = 0.40）。

### 5. 验证
`scripts/fable_irae/07_verify.py` 独立复算全部关键数值，追加参数检验与自助置信区间，并核对
对照，均通过：
- 上皮中 EPCAM 高 ✔ · 上皮中 PTPRC 近 0 ✔ · BALF 免疫区 PTPRC 高 ✔
- 机器可读结论与全部统计：`results/fable_irae/tables/verification_report.json`，
  汇总于 `results/fable_irae/tables/master_summary.tsv`。

### 6. 局限
单核数据低估膜蛋白 mRNA（压低 TACSTD2）；BAL 基本不含上皮；GSE319496 为转移性肾细胞癌
（系统性 irAE），非肺；样本量有限（26 / 51 / 7）；肺部 AE-vs-HC 对比存在批次混杂。
结果为假设层面，非临床结论。

### 7. 结论
在这些公共 ICI 数据集中，**TACSTD2 与 CLDN4 均未显示出与 irAE 状态显著、可重复的关联**。
两者相比，**CLDN4 更可测**（在结肠上皮中表达良好），而 **TACSTD2 基本处于检出下限**
（仅在捕获到真正上皮时可见），并且在全血矩阵中缺失。就实际而言，二者都不适合作为基于血液的
irAE 生物标志物；针对这些毒性的上皮标志物信号需在**靶组织活检**中寻找，而非 BAL 或血液。

---

### Files / 文件
- Scripts: `scripts/fable_irae/01_discover_geo.py` … `08_summary_figure.py`, `run.sh`
- Tables: `results/fable_irae/tables/` (`master_summary.tsv`, `verification_report.json`,
  per-dataset DE/measurability tables, `geo_candidates*.tsv`)
- Figures: `results/fable_irae/figures/` (`summary_measurability.png`,
  `colon_GSE206300_boxplots.png`, `blood_GSE319496_CLDN4.png`, `balf_GSE277136_celltype.png`)
