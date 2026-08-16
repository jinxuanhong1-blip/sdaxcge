# Tacstd2 / Cldn4 by treatment across open mouse-lung ICI RNA datasets

*Parallel slice — outputs confined to `notes/fable_mouse_ici/`,
`scripts/fable_mouse_ici/`, `results/fable_mouse_ici/`.*

---

## English

### Objective
Quantify how the epithelial tumour markers **Tacstd2 (Trop2)** and **Cldn4
(claudin-4)** change with treatment across ten publicly available mouse
lung-cancer immune-checkpoint-inhibitor (ICI) RNA datasets, using **processed
data only**, verifying every accession, skipping any single file > 2 GB (none
qualified), and computing **real statistics** where replication allows.

### Accessions (all verified, all *Mus musculus*, lung, ICI-context)
GSE239485, GSE297630, E-MTAB-13704, GSE241978, GSE330658 (replicated → tested);
GSE197260, GSE133604, GSE129297, GSE297632, GSE222158 (single sample per arm →
descriptive). Full evidence and per-dataset notes are in
[`notes/fable_mouse_ici/DATASETS.md`](../../notes/fable_mouse_ici/DATASETS.md)
and the JSON files beside it.

### Methods (brief)
- **Download** (`scripts/fable_mouse_ici/download_data.sh`): processed
  supplementary files from GEO FTP and EBI BioStudies into a temp cache (raw
  bundles are *not* committed; only derived tables/figures are).
- **Extraction** (`analyze.py`): per data type — log2-normalised matrices used
  as-is (GSE239485, GSE297630, GSE241978); TPM → log2(TPM+1) (GSE330658,
  GSE197260); raw counts → CPM → log2(CPM+1) (E-MTAB-13704); 10x single-cell →
  per-sample pseudobulk (sum over cells) → CPM → log2(CPM+1). Genes matched by
  symbol/alias, or by Ensembl ID / Clariom probe where needed.
- **Statistics** (`stats_and_plots.py`): for each treatment-vs-control contrast
  with ≥2 replicates per arm, a **Welch two-sample t-test** and a
  **Mann-Whitney U** test (two-sided) on log-scale expression, plus
  **Benjamini-Hochberg FDR** across all contrasts within a dataset. No p-values
  are produced for single-replicate arms.

### Headline results (replicate datasets, treatment − control on log2 scale)

| Dataset | Gene | Contrast | log2FC | Welch p | BH-FDR | Direction |
|---|---|---|---|---|---|---|
| **GSE239485** | **Tacstd2** | PolyIC+anti-PD1 vs Control | **+2.09** | **6.3e-4** | **0.0019** | ↑ significant |
| **GSE239485** | **Tacstd2** | PolyIC+anti-PD1+anti-C5aR1 vs Control | **+2.58** | **9.7e-4** | **0.0019** | ↑ significant |
| GSE239485 | Cldn4 | PolyIC+anti-PD1 vs Control | −0.77 | 0.24 | 0.32 | ns |
| **GSE297630** | **Cldn4** | anti-PD-1 vs Control | **−0.46** | **6.8e-4** | **0.0014** | ↓ significant (author FC −1.36, p≈9e-4) |
| GSE297630 | Tacstd2 | anti-PD-1 vs Control | −0.15 | 0.52 | 0.52 | ns |
| **E-MTAB-13704** | **Tacstd2** | VEGFRi/anti-PD-L1 vs vehicle | **+0.41** | **7.6e-4** | **0.0076** | ↑ significant (MWU p=0.016) |
| E-MTAB-13704 | Cldn4 | Cisplatin/anti-PD-L1/anti-CTLA4 vs vehicle | +0.52 | 0.083 | 0.42 | trend ns |
| GSE241978 | Cldn4 | AhR-KO vs Control | −0.59 | 0.18 | 0.18 | ns |
| GSE241978 | Tacstd2 | AhR-KO vs Control | 0.00 | — | — | at detection floor (not expressed) |
| GSE330658 | Tacstd2/Cldn4 | any arm vs Control | ±<0.5 | >0.12 | >0.6 | ns (n=2, low power) |

Full numbers: [`stats_results.csv`](stats_results.csv);
single-replicate values: [`descriptive_single_replicate.csv`](descriptive_single_replicate.csv);
per-sample data: [`per_sample_expression.csv`](per_sample_expression.csv);
figures in [`figures/`](figures/) (per-dataset box/strip plots plus
`overview_log2FC.png`).

### Interpretation
- **Tacstd2 (Trop2) rises after anti-PD-1-based therapy** in the strongest
  (n=8) dataset, GSE239485 — a large, FDR-significant induction (~4–6×) with
  Poly I:C + anti-PD-1, further with the triple combination. E-MTAB-13704
  independently shows a smaller but significant Tacstd2 increase specifically in
  the VEGFRi/anti-PD-L1 arm. The single-cell LLC series diverge at n=1
  (GSE297632 Tacstd2 slightly ↓; GSE129297 combo ↓), so the induction is
  supported mainly by the replicated bulk data.
- **Cldn4 tends to fall with anti-PD-1**: significantly in GSE297630 (matching
  the authors' own statistic) and non-significantly in AhR-KO (GSE241978) and
  the GSE239485 anti-PD-1 arm.
- Net picture from the datasets that permit inference: **anti-PD-1/PD-L1
  treatment is associated with up-regulation of Tacstd2/Trop2 and mild
  down-regulation of Cldn4** in mouse lung tumours, though effect sizes and
  significance vary by model and combination partner.

### Caveats
Single-replicate arms (GSE197260 and all four single-cell series) are
descriptive only. GSE129297 uses unfiltered 10x matrices (noisier pseudobulk).
GSE222158 immune-sorted (CD45/CD3) samples largely lack epithelial cells, so low
Tacstd2/Cldn4 there reflects cell composition, not treatment. GSE241978 is a
genotype knockout, not an antibody arm. GSE330658's RNA subset lacks the study's
anti-PD-L1 arm and has n=2.

### Reproduce
```bash
python3 scripts/fable_mouse_ici/verify_accessions.py
python3 scripts/fable_mouse_ici/inventory_sizes.py
python3 scripts/fable_mouse_ici/fetch_sample_metadata.py
bash   scripts/fable_mouse_ici/download_data.sh
python3 scripts/fable_mouse_ici/analyze.py
python3 scripts/fable_mouse_ici/stats_and_plots.py
```

---

## 中文

### 目标
在十个公开的**小鼠肺癌免疫检查点抑制剂（ICI）RNA 数据集**中，量化上皮性肿瘤标志物
**Tacstd2（即 Trop2）** 与 **Cldn4（claudin-4）** 随治疗的表达变化。要求：**仅使用处理后的数据**、
逐一核实每个编号、跳过任何单个 > 2 GB 的文件（实际无一超标）、并在有生物学重复的情况下计算**真实统计量**。

### 数据集（全部核实，均为小鼠、肺、ICI 背景）
可做统计（有重复）：GSE239485、GSE297630、E-MTAB-13704、GSE241978、GSE330658；
仅描述（每组单样本）：GSE197260、GSE133604、GSE129297、GSE297632、GSE222158。
完整证据与逐数据集说明见
[`notes/fable_mouse_ici/DATASETS.md`](../../notes/fable_mouse_ici/DATASETS.md) 及同目录 JSON。

### 方法（简述）
- **下载**：从 GEO FTP 与 EBI BioStudies 获取处理后的补充文件到临时缓存（原始大文件不入库，仅提交衍生表格与图）。
- **提取**：按数据类型处理——已 log2 归一化矩阵直接使用（GSE239485/GSE297630/GSE241978）；
  TPM → log2(TPM+1)（GSE330658/GSE197260）；原始 counts → CPM → log2(CPM+1)（E-MTAB-13704）；
  10x 单细胞 → 每样本 pseudobulk（按基因对所有细胞求和）→ CPM → log2(CPM+1)。基因按符号/别名匹配，
  必要时按 Ensembl ID 或 Clariom 探针匹配。
- **统计**：对每个"治疗 vs 对照"、且两组各 ≥2 重复的对比，做 **Welch 双样本 t 检验** 与
  **Mann-Whitney U 检验**（双侧），并在每个数据集内部做 **Benjamini-Hochberg FDR** 校正；
  单样本组不给出 p 值。

### 主要结果（有重复的数据集，log2 尺度上"治疗−对照"）

| 数据集 | 基因 | 对比 | log2FC | Welch p | BH-FDR | 方向 |
|---|---|---|---|---|---|---|
| **GSE239485** | **Tacstd2** | PolyIC+anti-PD1 vs 对照 | **+2.09** | **6.3e-4** | **0.0019** | ↑ 显著 |
| **GSE239485** | **Tacstd2** | PolyIC+anti-PD1+anti-C5aR1 vs 对照 | **+2.58** | **9.7e-4** | **0.0019** | ↑ 显著 |
| GSE239485 | Cldn4 | PolyIC+anti-PD1 vs 对照 | −0.77 | 0.24 | 0.32 | 不显著 |
| **GSE297630** | **Cldn4** | anti-PD-1 vs 对照 | **−0.46** | **6.8e-4** | **0.0014** | ↓ 显著（作者报告 FC −1.36，p≈9e-4） |
| GSE297630 | Tacstd2 | anti-PD-1 vs 对照 | −0.15 | 0.52 | 0.52 | 不显著 |
| **E-MTAB-13704** | **Tacstd2** | VEGFRi/anti-PD-L1 vs vehicle | **+0.41** | **7.6e-4** | **0.0076** | ↑ 显著（MWU p=0.016） |
| E-MTAB-13704 | Cldn4 | Cisplatin/anti-PD-L1/anti-CTLA4 vs vehicle | +0.52 | 0.083 | 0.42 | 趋势，不显著 |
| GSE241978 | Cldn4 | AhR-KO vs 对照 | −0.59 | 0.18 | 0.18 | 不显著 |
| GSE241978 | Tacstd2 | AhR-KO vs 对照 | 0.00 | — | — | 处于检测下限（未表达） |
| GSE330658 | Tacstd2/Cldn4 | 各治疗组 vs 对照 | 幅度<0.5 | >0.12 | >0.6 | 不显著（n=2，功效低） |

完整数值见 [`stats_results.csv`](stats_results.csv)；单样本组见
[`descriptive_single_replicate.csv`](descriptive_single_replicate.csv)；逐样本数据见
[`per_sample_expression.csv`](per_sample_expression.csv)；图见 [`figures/`](figures/)。

### 解读
- **Tacstd2/Trop2 在以 anti-PD-1 为基础的治疗后升高**：在样本量最大（n=8）的 GSE239485 中，
  Poly I:C + anti-PD-1 使其显著升高（约 4–6 倍，FDR 显著），三联组合进一步升高；
  E-MTAB-13704 独立地在 VEGFRi/anti-PD-L1 组显示较小但显著的升高。单细胞 LLC 系列在 n=1 时方向不一
  （GSE297632 略降；GSE129297 combo 降），因此该诱导主要由有重复的 bulk 数据支持。
- **Cldn4 在 anti-PD-1 下倾向下降**：在 GSE297630 中显著（与作者统计一致），在 AhR-KO（GSE241978）
  与 GSE239485 的 anti-PD-1 组呈非显著下降。
- 综合可推断的数据集：**在小鼠肺肿瘤中，anti-PD-1/PD-L1 治疗与 Tacstd2/Trop2 上调、Cldn4 轻度下调相关**，
  但效应量与显著性因模型与联合用药而异。

### 注意事项
单样本组（GSE197260 及四个单细胞系列）仅作描述。GSE129297 使用未过滤的 10x 矩阵（pseudobulk 更嘈杂）。
GSE222158 为免疫细胞分选（CD45/CD3），基本不含上皮细胞，故 Tacstd2/Cldn4 偏低反映的是细胞组成而非治疗效应。
GSE241978 为基因敲除而非抗体治疗组。GSE330658 的 RNA 子集不含研究中的 anti-PD-L1 组且 n=2。

### 复现
见上方英文"Reproduce"代码块中的脚本执行顺序。
