# TACSTD2 / CLDN4 in the blood of NSCLC patients under immune-checkpoint inhibition
# NSCLC 免疫检查点治疗患者外周血中的 TACSTD2 / CLDN4

*Slice: `fable_blood_ici`. All outputs live under `notes/fable_blood_ici/`,
`scripts/fable_blood_ici/`, `results/fable_blood_ici/`.*

---

## English

### Question
TACSTD2 (TROP2) and CLDN4 are epithelial cell-surface antigens and antibody–drug-
conjugate targets. This slice asks, on **real human NSCLC blood transcriptomic
data from GEO**, two honest questions:
1. Are TACSTD2 / CLDN4 detectable at all in the blood compartment (PBMC / whole
   blood)?
2. If detectable, do their levels associate with immune-checkpoint-inhibitor
   (ICI) response or immune-related adverse events (irAE)?

The task flagged in advance that these targets may be low or undetectable in
blood — so the deliverable is an honest, statistically-backed measurement, not a
forced positive.

### Datasets (both processed, both < 2 GB)
| GEO | Assay / compartment | Samples | Labels used |
|-----|---------------------|---------|-------------|
| **GSE285888** (task-named) | scRNA-seq, baseline **PBMC** | 222,144 cells / 33 NSCLC pts | CR / DR / PD response, irAE severity |
| **GSE305086** | Affymetrix HG-U133 Plus 2.0, **whole blood** | 173 (76 baseline + 76 follow-up + 21 control) | sample type, treatment regimen |

Full provenance and the gene panel are in `datasets.md`.

### Method (real statistics)
- **GSE285888.** Streamed the 486 MB genes × cells matrix, extracting a target +
  reference panel. Detection rate = fraction of cells with count > 0;
  magnitude = CP10K. Association tested on **patient-level pseudobulk** with the
  two-sided **Mann-Whitney U** test (unit = patient, not cell).
- **GSE305086.** Detectability = percentile rank of each probe's mean log2
  intensity among all 26,452 probes (left tail = array background). Disease
  effect via **Mann-Whitney U** (baseline vs control); treatment effect via
  paired **Wilcoxon signed-rank** (baseline vs follow-up).

### Results

**1. Detectability — the targets are essentially absent from circulating cells.**

GSE285888 baseline PBMC (fraction of 222,144 cells with any UMI):

| Gene | Class | % cells detected | mean CP10K |
|------|-------|-----------------:|-----------:|
| PTPRC/CD45 | immune | 97.77 % | 22.46 |
| CD3D | immune | 46.24 % | 4.99 |
| NKG7 | immune | 36.77 % | 18.60 |
| GZMB | effector | 27.24 % | 4.30 |
| **CLDN4** | **target** | **0.18 %** | **0.0039** |
| EPCAM | epithelial | 0.08 % | 0.0015 |
| **TACSTD2** | **target** | **0.08 %** | **0.0012** |

TACSTD2 is detected in **179 / 222,144 cells (0.08 %)** and CLDN4 in **404
(0.18 %)** — the same near-zero range as the epithelial control EPCAM (0.08 %).
This is at the level of ambient/contaminating epithelial transcripts, i.e.
**effectively undetectable in PBMC**. See `gse285888_detectability.png`.

GSE305086 whole-blood array (percentile rank of probe mean among 26,452 probes;
low = array background):

| Gene | Probe | percentile | interpretation |
|------|-------|-----------:|----------------|
| PTPRC | 212587_s_at | 99.3 | clearly expressed |
| PRF1 | 214617_at | 97.6 | clearly expressed |
| **CLDN4** | 1569421_at | **70.7** | low-moderate signal |
| **CLDN4** | 201428_at | **26.4** | at background |
| **TACSTD2** | 202286_s_at | **42.0** | near background |
| **TACSTD2** | 202287_s_at | **20.1** | at background |
| EPCAM | 201839_s_at | 18.1 | background floor |

On bulk whole blood, TACSTD2 sits at the array background (its probes rank
20th–42nd percentile, next to the EPCAM background floor at the 18th). CLDN4 is
**probe-discordant**: its canonical probe `201428_at` is at background (26th)
while `1569421_at` reaches the 70th percentile — so CLDN4 may carry a low-level
blood signal, but the evidence is weak and probe-dependent. See
`gse305086_detectability.png`.

**2. Association with ICI response / irAE — none for the targets.**

Patient-level pseudobulk in GSE285888 (Mann-Whitney U, two-sided):

| Gene | Contrast | metric | p |
|------|----------|--------|--:|
| TACSTD2 | CR vs PD | detection | 0.96 |
| CLDN4 | CR vs PD | CP10K | 0.54 |
| TACSTD2 | severe-irAE vs non | CP10K | 0.68 |
| CLDN4 | severe-irAE vs non | CP10K | 0.46 |

No contrast approaches significance. In GSE305086 the targets are statistically
different between baseline and controls (TACSTD2 p = 0.034; CLDN4 p = 0.0015) and
CLDN4 rises slightly at follow-up (p = 0.005), but the effect sizes are ~0.1–0.2
log2 units at/near background — statistically detectable thanks to n = 76 but
biologically marginal.

**3. The pipeline is sound (positive controls fire).** The same pipeline recovers
the source studies' real signals: in GSE285888 IL1B and CXCL8 are significantly
elevated in severe-irAE patients (p = 0.026 each, the study's reported irAE
predictors), and in GSE305086 GATA3 is strongly depleted in patient baseline
blood vs controls (p = 1.3 × 10⁻¹⁰). So the null result for TACSTD2 / CLDN4 is a
real biological negative, not a processing artefact.

### Honest conclusion
In human NSCLC blood, **TACSTD2 is effectively undetectable** (PBMC ≈ 0.08 % of
cells, whole blood at array background) and **CLDN4 is at best very low and
probe-inconsistent**. Neither shows any association with ICI response or irAE.
Blood transcriptomics is therefore **not** a viable readout for these two
epithelial ADC targets in NSCLC; tumour tissue is required to profile them.

### Reproduce / files
`scripts/fable_blood_ici/01_download.sh` → `02_extract_gse285888.py` →
`03_analyze_gse285888.py` → `04_analyze_gse305086.py`. Result tables, figures and
JSON summaries are in `results/fable_blood_ici/`.

---

## 中文

### 问题
TACSTD2（TROP2）与 CLDN4 是上皮细胞表面抗原，也是抗体偶联药物（ADC）的靶点。
本切片基于 **GEO 上真实的人类 NSCLC 血液转录组数据**，诚实地回答两个问题：
1. 在血液区室（PBMC / 全血）中，TACSTD2 / CLDN4 是否能被检出？
2. 若可检出，其水平是否与免疫检查点抑制剂（ICI）疗效或免疫相关不良事件
   （irAE）相关？

任务已预先提示：这两个靶点在血液中可能很低甚至检不到。因此交付物是一个诚实、
有统计支撑的测量结果，而非强行做出的阳性结论。

### 数据集（均为处理后数据，均 < 2 GB）
| GEO | 平台 / 区室 | 样本 | 使用的标签 |
|-----|-------------|------|-----------|
| **GSE285888**（任务指定） | scRNA-seq，基线 **PBMC** | 222,144 细胞 / 33 例 NSCLC | CR / DR / PD 疗效、irAE 分级 |
| **GSE305086** | Affymetrix HG-U133 Plus 2.0，**全血** | 173（76 基线 + 76 随访 + 21 对照） | 样本类型、治疗方案 |

数据来源与基因面板详见 `datasets.md`。

### 方法（真实统计）
- **GSE285888**：流式读取 486 MB 的“基因 × 细胞”矩阵，仅提取靶点+参考基因面板。
  检出率 = 计数 > 0 的细胞比例；表达量 = CP10K。关联分析在**患者层面
  pseudobulk**（每位患者一个数值）上用双侧 **Mann-Whitney U** 检验，统计单元是
  患者而非细胞。
- **GSE305086**：可检出性 = 每个探针平均 log2 强度在全部 26,452 个探针中的
  百分位（越靠左尾 = 越接近芯片背景）。疾病效应用 **Mann-Whitney U**（基线
  vs 对照）；治疗效应用配对 **Wilcoxon 符号秩**（基线 vs 随访）。

### 结果

**1. 可检出性——靶点在循环细胞中基本不存在。**

GSE285888 基线 PBMC（222,144 细胞中有任意 UMI 的比例）：

| 基因 | 类别 | 检出细胞比例 | 平均 CP10K |
|------|------|-----------:|-----------:|
| PTPRC/CD45 | 免疫 | 97.77 % | 22.46 |
| CD3D | 免疫 | 46.24 % | 4.99 |
| GZMB | 效应 | 27.24 % | 4.30 |
| **CLDN4** | **靶点** | **0.18 %** | **0.0039** |
| EPCAM | 上皮 | 0.08 % | 0.0015 |
| **TACSTD2** | **靶点** | **0.08 %** | **0.0012** |

TACSTD2 仅在 **179 / 222,144 个细胞（0.08 %）** 中检出，CLDN4 在 **404 个
（0.18 %）** 中检出——与上皮对照 EPCAM（0.08 %）处于同一近零区间，属于环境/
污染性上皮转录本水平，即 **在 PBMC 中实际不可检出**（见
`gse285888_detectability.png`）。

GSE305086 全血芯片（探针平均值在 26,452 探针中的百分位，越低 = 越接近背景）：

| 基因 | 探针 | 百分位 | 解读 |
|------|------|-----:|------|
| PTPRC | 212587_s_at | 99.3 | 明确表达 |
| PRF1 | 214617_at | 97.6 | 明确表达 |
| **CLDN4** | 1569421_at | **70.7** | 低-中等信号 |
| **CLDN4** | 201428_at | **26.4** | 位于背景 |
| **TACSTD2** | 202286_s_at | **42.0** | 接近背景 |
| **TACSTD2** | 202287_s_at | **20.1** | 位于背景 |
| EPCAM | 201839_s_at | 18.1 | 背景基线 |

在全血 bulk 芯片上，TACSTD2 处于芯片背景（各探针位于第 20–42 百分位，紧邻
EPCAM 背景基线的第 18 位）。CLDN4 **探针不一致**：经典探针 `201428_at` 位于
背景（第 26 位），而 `1569421_at` 达到第 70 位——因此 CLDN4 或许存在微弱血液
信号，但证据薄弱且依赖探针（见 `gse305086_detectability.png`）。

**2. 与 ICI 疗效 / irAE 的关联——靶点均无。**

GSE285888 患者层面 pseudobulk（双侧 Mann-Whitney U）：

| 基因 | 对比 | 指标 | p |
|------|------|------|--:|
| TACSTD2 | CR vs PD | 检出率 | 0.96 |
| CLDN4 | CR vs PD | CP10K | 0.54 |
| TACSTD2 | 重度 irAE vs 无 | CP10K | 0.68 |
| CLDN4 | 重度 irAE vs 无 | CP10K | 0.46 |

无任何对比接近显著。GSE305086 中靶点在基线与对照间存在统计差异
（TACSTD2 p = 0.034；CLDN4 p = 0.0015），CLDN4 在随访时略升（p = 0.005），
但效应量仅约 0.1–0.2 个 log2 单位且处于背景附近——因 n = 76 而可被统计检出，
生物学意义甚微。

**3. 流程可靠（阳性对照成立）。** 同一流程能复现原研究的真实信号：GSE285888 中
IL1B 与 CXCL8 在重度 irAE 患者中显著升高（各 p = 0.026，正是该研究报告的
irAE 预测因子）；GSE305086 中 GATA3 在患者基线全血中相对对照显著降低
（p = 1.3 × 10⁻¹⁰）。因此 TACSTD2 / CLDN4 的阴性结果是真实的生物学阴性，
而非处理伪影。

### 诚实结论
在人类 NSCLC 血液中，**TACSTD2 实际不可检出**（PBMC 约 0.08 % 细胞，全血处于
芯片背景），**CLDN4 至多极低且探针不一致**。两者均与 ICI 疗效或 irAE 无关联。
因此血液转录组 **不适合** 作为 NSCLC 中这两个上皮 ADC 靶点的读出手段，需使用
肿瘤组织进行检测。

### 复现 / 文件
`scripts/fable_blood_ici/01_download.sh` → `02_extract_gse285888.py` →
`03_analyze_gse285888.py` → `04_analyze_gse305086.py`。结果表格、图像与 JSON
摘要位于 `results/fable_blood_ici/`。
