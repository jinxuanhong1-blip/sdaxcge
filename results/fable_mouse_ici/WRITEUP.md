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
- Honest net from the original ten: the **GSE239485-style Tacstd2 rise on IO
  is real in that dataset** (and the E-MTAB VEGFRi/aPD-L1 combo), but
  anti-PD-1/PD-L1 **monotherapy is not uniformly Tacstd2-up** (GSE297630 and
  E-MTAB aPD-L1 monotherapy are slightly down, ns). See the hunt extension
  below for a pre-specified combined test.

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
python3 scripts/fable_mouse_ici/hunt_new_accessions.py
python3 scripts/fable_mouse_ici/extend_icb_immune.py
```

---

## Extension — more paired ICB sets, combined Tacstd2-up test, TROP2-high / immune-low

### Hunt (processed, <2 GB)
GEO was re-queried for mouse lung + anti-PD-1/PD-L1 RNA. New **tumour** RNA series with a control-vs-ICB arm and a processed file <2 GB:

| Accession | Model | Primary ICB contrast | n | Tacstd2 log2FC | Welch p (two-sided) |
|---|---|---|---|---|---|
| GSE114601 | GEMM NSCLC | anti-PD1 vs Vehicle | 2 vs 2 | **−1.16** | 0.45 |
| GSE157880 | KP lung | PD-1 0 Gy vs IgG 0 Gy | 2 vs 3 | +0.35 | 0.43 |
| GSE309199 | RPM SCLC | aPD1 vs Ctrl | 3 vs 3 | +0.49 | 0.47 |
| GSE169196 | KPM total tumour | A2V+aPD1 vs IgG *(no aPD1-mono arm)* | 3 vs 3 | +0.13 | 0.84 |

Skipped on purpose (see `notes/fable_mouse_ici/hunt_triage.md`): GSE190264 (in-vitro chemo/MEK, no ICB), GSE114300 / GSE277610 (sorted T cells), GSE184000 (irAE whole lung, not tumour), GSE315010 (RAS inhibitors, no ICB antibody), GSE194166 / GSE129298 (scRNA n=1/arm).

TISMO 49/64 Tacstd2-up after ICB (p=5.8e-5) is treated as an **external all-cancer prior**. The current TISMO site is a SPA; the GitHub repo (`zexian/TISMO_data`) has processing scripts only. We **did not independently recompute** that 49/64 count.

### Combined Tacstd2-up test (one primary contrast per independent tumour dataset)

Pre-specified set (n=7): GSE239485 PolyIC+aPD1; GSE297630 aPD1; E-MTAB-13704 aPD-L1 monotherapy; GSE114601 aPD1; GSE157880 PD-1 0 Gy; GSE309199 aPD1; GSE169196 A2V+aPD1.

| Metric | Result |
|---|---|
| Direction | **4 / 7** Tacstd2 log2FC > 0 |
| Sign test (greater) | **p = 0.50** |
| Stouffer Z on one-sided Welch p (up) | Z = 0.85, **p = 0.20** |
| Only two-sided p<0.05 in the primary set | GSE239485 (+2.09, p=6.3e-4) |

So: **GSE239485-style Tacstd2 rise on IO is consistent with that one well-powered LLC experiment**, but **this public mouse-lung slice does not reproduce a TISMO-like 49/64 consensus**. Anti-PD-1/PD-L1 monotherapy in GEMM/SCLC/LLC-array data is mixed and under-powered. Tables: [`primary_icb_tacstd2.csv`](primary_icb_tacstd2.csv), [`combined_icb_tacstd2.json`](combined_icb_tacstd2.json), figure [`figures/primary_icb_tacstd2.png`](figures/primary_icb_tacstd2.png).

### Sensitivity (same primary table, pre-specified subsets)

| Set | n | Tacstd2 up | Sign p (up) | Stouffer p (up) |
|---|---|---|---|---|
| All 7 primary | 7 | 4 | 0.50 | 0.20 |
| Drop GSE239485 (Poly I:C confounder) | 6 | 3 | 0.66 | 0.68 |
| **Monotherapy only** | 5 | **2** | **0.81** | **0.73** |
| Monotherapy with n≥3/arm | 3 | 1 | 0.88 | 0.79 |

Dropping the Poly I:C experiment **removes the only significant Tacstd2-up call**. Monotherapy alone is 2 up / 3 down. [`sensitivity_combined.json`](sensitivity_combined.json), [`figures/sensitivity_tacstd2.png`](figures/sensitivity_tacstd2.png).

**Radiation caveat (GSE157880):** PD-1 + 4 Gy vs IgG 0 Gy looks like Tacstd2 up (+1.42, p=0.003), but IgG + 4 Gy vs IgG 0 Gy is also up (+1.78, p=0.013). That rise is **radiation, not ICB**. The primary contrast correctly uses the 0 Gy pair.

**Cldn4 on the same 7 primary contrasts:** 6 / 7 down; up-direction sign p=0.99; down-direction sign p=**0.062** (trend, not <0.05). Only GSE297630 is individually significant. [`primary_icb_cldn4.csv`](primary_icb_cldn4.csv), [`figures/primary_icb_cldn4.png`](figures/primary_icb_cldn4.png).

TISMO download was retried (current SPA + `/tismo` `/rtismo` `/datadownload` APIs all 404; GitHub has scripts only). **49/64 was not recomputed.**

### TROP2-high / immune-low
Immune score = mean log-expression of Cd8a, Cd3e, Cd3d, Gzmb, Prf1, Ifng, Cd274, Pdcd1, Cxcl9, Cxcl10, Nkg7 (whichever present).

| Dataset | Spearman ρ (Tacstd2 vs immune, all samples) | p | Median-split: Tacstd2-high has lower immune? |
|---|---|---|---|
| GSE239485 | **+0.56** | 0.0048 | no (treatment-driven: IO raises both) |
| GSE239485 control-only | −0.10 | 0.82 | — |
| E-MTAB-13704 | **+0.58** | 0.0014 | no |
| GSE157880 | −0.20 | 0.46 | yes, ns |
| GSE197260 | −0.21 | 0.64 | yes, ns (n=7, 1/arm) |
| GSE114601 / GSE309199 / GSE169196 | ~0 | >0.5 | no |

**No robust TROP2-high / immune-low subset in these bulk lung matrices.** The two significant correlations are **positive** (Tacstd2 tracks with the immune score), largely because IO/combo arms move both. Baseline-only splits are under-powered and non-significant. [`tacstd2_immune_correlation.csv`](tacstd2_immune_correlation.csv), [`trop2_high_immune_split.csv`](trop2_high_immune_split.csv), [`figures/tacstd2_vs_immune.png`](figures/tacstd2_vs_immune.png).

---

## 增补 — 更多配对 ICB 集、Tacstd2 上调合并检验、TROP2 高 / 免疫低

### 检索
在 GEO 中补检小鼠肺 + anti-PD-1/PD-L1 RNA。新增符合“处理后、<2 GB、肿瘤 RNA、对照 vs ICB”的：GSE114601、GSE157880、GSE309199、GSE169196（后者无 aPD1 单药、只有 A2V+aPD1）。刻意跳过的编号见 `notes/fable_mouse_ici/hunt_triage.md`。

TISMO 的 49/64、p=5.8e-5 作为**外部（全癌种）先验**引用；当前 TISMO 站点为前端 SPA，GitHub 仅有处理脚本，**未能独立重算**该 49/64。

### 合并检验（每个独立肿瘤数据集只取一个主对比，共 7 个）
Tacstd2 升高 **4/7**；符号检验 p=**0.50**；Stouffer（单侧上调）Z=0.85，p=**0.20**。主对比里唯一双侧显著的是 GSE239485（+2.09，p=6.3e-4）。

结论：**GSE239485 式的 IO 后 Tacstd2 升高在该实验中成立**，但**本公开小鼠肺切片不能复现 TISMO 式的 49/64 共识**；单药 anti-PD-1/PD-L1 方向混杂且功效不足。

敏感性：去掉 GSE239485（Poly I:C 混杂）后 3/6 升高（符号 p=0.66，Stouffer p=0.68）；**仅单药**为 2/5 升高（p=0.81 / 0.73）。GSE157880 的 4 Gy 升高是放疗而非 ICB（IgG 4 Gy 同样升高）。同一 7 个主对比上 Cldn4 为 6/7 下降（下降方向符号 p=0.062，趋势未过 0.05）。TISMO 接口仍 404，49/64 **未重算**。

### TROP2 高 / 免疫低
免疫评分 = Cd8a/Cd3e/Cd3d/Gzmb/Prf1/Ifng/Cd274/Pdcd1/Cxcl9/Cxcl10/Nkg7 的均值。GSE239485 与 E-MTAB 的 Spearman 为**显著正相关**（治疗同时抬高两者）；对照子集与其余数据集为弱负或近零、均不显著。**这些 bulk 肺矩阵中没有稳健的 TROP2-high / immune-low 亚群证据。**

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
- 诚实综合：GSE239485 式的 IO 后 Tacstd2 升高在该数据集中是真实的（E-MTAB 的 VEGFRi/aPD-L1
  联合组也升高），但 **anti-PD-1/PD-L1 单药并不一律上调 Tacstd2**（GSE297630 与 E-MTAB aPD-L1
  单药略降、不显著）。预指定的合并检验见下方增补。

### 注意事项
单样本组（GSE197260 及四个单细胞系列）仅作描述。GSE129297 使用未过滤的 10x 矩阵（pseudobulk 更嘈杂）。
GSE222158 为免疫细胞分选（CD45/CD3），基本不含上皮细胞，故 Tacstd2/Cldn4 偏低反映的是细胞组成而非治疗效应。
GSE241978 为基因敲除而非抗体治疗组。GSE330658 的 RNA 子集不含研究中的 anti-PD-L1 组且 n=2。

### 复现
见上方英文"Reproduce"代码块中的脚本执行顺序。
