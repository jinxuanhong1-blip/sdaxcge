# WRITEUP — TROP2 / CLDN ADC ± ICI in Lung: Public Data Catalog & Expression-vs-Outcome Analysis

Slice: `fable_adc`. Outputs live only under `notes/fable_adc/`, `scripts/fable_adc/`, `results/fable_adc/`.

---

## English

### 1. Objective
Catalog **public** omics data relevant to TROP2-directed antibody–drug conjugates (ADCs) —
sacituzumab govitecan (SG) and datopotamab deruxtecan (Dato-DXd) — with or without immune-checkpoint
inhibitors (ICI) in lung cancer, plus CLDN18.2-targeted agents where any public omics touch lung.
Then, using only openly downloadable patient-level data, test whether the ADC target genes track
with immunotherapy outcome.

### 2. What is (and is not) public
The pivotal TROP2-ADC biomarker readouts have **no public patient-level omics**:
- **TROPION-Lung01** (Dato-DXd vs docetaxel, [NCT04656652](https://clinicaltrials.gov/study/NCT04656652)):
  the predictive biomarker is a computational-pathology **TROP2 QCS normalized membrane ratio (NMR)**
  (QCS-NMR+ = ≥75% of tumor cells with membrane/cytoplasm ratio ≤0.56), reported in
  [JCO 10.1200/JCO-24-01544](https://ascopubs.org/doi/10.1200/JCO-24-01544) and WCLC24. Image/scoring
  data are proprietary (AstraZeneca QCS).
- **EVOKE-01 / EVOKE-02** (SG ± pembrolizumab, [NCT05089734](https://clinicaltrials.gov/study/NCT05089734) /
  [NCT05186974](https://clinicaltrials.gov/study/NCT05186974)): TROP2 was assessed by IHC **H-score**;
  both AACR 2025 abstracts ([LB260](https://doi.org/10.1158/1538-7445.am2025-lb260),
  [LB399](https://doi.org/10.1158/1538-7445.am2025-lb399)) conclude efficacy is **independent** of TROP2
  expression. No expression matrices released.

So there is **no honest way** to analyze ADC response vs expression on public data. Instead we use
three fully-open GEO NSCLC cohorts treated with ICI (or chemo-ICI) as surrogates to characterize the
ADC-relevant genes against immunotherapy outcome. Full provenance is in
[`notes/fable_adc/data_catalog.md`](./data_catalog.md). **No accession was invented.**

### 3. Datasets analyzed (open)
| Accession | Treatment | n | Outcome |
|---|---|---|---|
| [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044) | anti–PD-1 mono | 16 | responder (5) / non-responder (11) |
| [GSE135222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222) | anti–PD-(L)1 | 27 | PFS (21 events) |
| [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) | neoadjuvant chemo + anti–PD-1 | 24 | MPR (9) / NMPR (15) + residual-tumor % |

Genes: `TACSTD2` (TROP2 target), `CLDN18` (Claudin-18; 18.2 is the therapeutic epitope but bulk RNA-seq
cannot resolve the isoform), `TOP1` (SG/Dato-DXd payload target), `CD274` (PD-L1).

### 4. Methods
Reproduced by `scripts/fable_adc/download_data.sh` then `scripts/fable_adc/analyze.py`.
- GSE126044 raw counts → log2(CPM+1); Mann–Whitney U (responder vs non-responder), rank-biserial effect.
- GSE135222 TPM → log2(TPM+1); median-split Kaplan–Meier + log-rank on PFS, and univariate Cox HR per SD.
- GSE207422 log2-TPM; Mann–Whitney (MPR vs NMPR) and Spearman vs residual-tumor fraction.
- TROP2–PD-L1 (`TACSTD2`–`CD274`) Spearman correlation in each cohort.

### 5. Results (see `results/fable_adc/tables/` and `figures/`)
- **TROP2 (`TACSTD2`) is not associated with ICI outcome** in any cohort: responder-vs-non-responder
  Mann–Whitney p=0.44 (GSE126044); MPR-vs-NMPR p=0.34 (GSE207422); PFS log-rank p=0.43, Cox HR/SD 1.07
  (95% CI 0.68–1.67) (GSE135222). This is directionally consistent with the EVOKE-01/02 conclusion that
  SG/SG+pembro efficacy is independent of TROP2 expression.
- **CLDN18 is essentially unexpressed in unselected NSCLC**: median log2-TPM ≈ 0.07–0.08 in GSE207422,
  no MPR/NMPR difference (p=1.0). This matches the clinical picture that CLDN18.2 in lung is restricted
  largely to invasive mucinous adenocarcinoma; unselected NSCLC is a poor CLDN18.2 indication.
- **TOP1** (payload target) trended higher-expression → worse PFS (Cox HR/SD 1.54, p=0.10) in GSE135222;
  not significant.
- **PD-L1 (`CD274`)** behaved as expected for ICI — higher in responders (GSE126044, rank-biserial
  −0.42) and protective PFS trend (GSE135222 Cox HR/SD 0.68, p=0.15) — but underpowered.
- **TROP2 and PD-L1 are not positively correlated** (Spearman ρ = −0.23 / +0.02 / −0.18 across the three
  cohorts), i.e. TROP2 expression carries information largely orthogonal to PD-L1 — a biological
  rationale for combining a TROP2-ADC with ICI rather than expecting redundancy.

### 6. Interpretation & limitations
Public transcriptomes reinforce that **TROP2 mRNA is not a useful predictive biomarker for
immunotherapy outcome**, echoing why the ADC field moved to spatial/protein metrics (QCS-NMR, H-score).
Cohorts are small (n=16–27), none received a TROP2-ADC, and `CLDN18` gene-level signal cannot separate
the 18.2 isoform. Findings are hypothesis-level and describe target biology, not ADC response.

---

## 中文

### 1. 目标
系统梳理与 TROP2 抗体偶联药物（ADC）——sacituzumab govitecan（SG）、datopotamab deruxtecan（Dato-DXd）——
在肺癌中联合或不联合免疫检查点抑制剂（ICI）相关的**公开**组学数据，并补充有肺癌公开组学的 CLDN18.2 靶向药物。
随后，仅使用可公开下载的患者级数据，检验 ADC 靶基因是否与免疫治疗结局相关。

### 2. 哪些是公开的、哪些不是
关键的 TROP2-ADC 生物标志物结果**均无公开的患者级组学数据**：
- **TROPION-Lung01**（Dato-DXd 对比多西他赛，[NCT04656652](https://clinicaltrials.gov/study/NCT04656652)）：
  预测性标志物是基于计算病理的 **TROP2 QCS 归一化膜比值（NMR）**（QCS-NMR 阳性 = ≥75% 肿瘤细胞膜/胞质比 ≤0.56），
  见 [JCO 10.1200/JCO-24-01544](https://ascopubs.org/doi/10.1200/JCO-24-01544) 与 WCLC24；图像与评分数据为
  AstraZeneca QCS 专有。
- **EVOKE-01 / EVOKE-02**（SG ± 帕博利珠单抗）：以 IHC **H-score** 评估 TROP2；两篇 AACR 2025 摘要
  （[LB260](https://doi.org/10.1158/1538-7445.am2025-lb260)、[LB399](https://doi.org/10.1158/1538-7445.am2025-lb399)）
  均结论疗效**与 TROP2 表达无关**，且未公开表达矩阵。

因此在公开数据上**无法诚实地**直接分析「ADC 疗效 vs 表达」。我们改用三套完全公开、接受 ICI（或化疗+ICI）的
NSCLC GEO 队列作为替代，来刻画 ADC 相关基因与免疫治疗结局的关系。完整来源见
[`notes/fable_adc/data_catalog.md`](./data_catalog.md)。**未虚构任何编号。**

### 3. 所分析的数据集（公开）
| 编号 | 治疗 | n | 结局 |
|---|---|---|---|
| [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044) | 抗 PD-1 单药 | 16 | 应答 (5) / 无应答 (11) |
| [GSE135222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222) | 抗 PD-(L)1 | 27 | PFS（21 事件） |
| [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) | 新辅助化疗 + 抗 PD-1 | 24 | MPR (9) / NMPR (15) + 残留肿瘤比例 |

基因：`TACSTD2`（TROP2 靶点）、`CLDN18`（Claudin-18；18.2 为治疗表位，但 bulk RNA-seq 无法区分亚型）、
`TOP1`（SG/Dato-DXd 载荷靶点）、`CD274`（PD-L1）。

### 4. 方法
由 `scripts/fable_adc/download_data.sh` 与 `scripts/fable_adc/analyze.py` 可完全复现。
- GSE126044 原始 counts → log2(CPM+1)；Mann–Whitney U（应答 vs 无应答）+ rank-biserial 效应量。
- GSE135222 TPM → log2(TPM+1)；按中位数分组的 Kaplan–Meier + log-rank（PFS），及单因素 Cox（每 SD 的 HR）。
- GSE207422 log2-TPM；Mann–Whitney（MPR vs NMPR）与对残留肿瘤比例的 Spearman 相关。
- 各队列中 TROP2–PD-L1（`TACSTD2`–`CD274`）的 Spearman 相关。

### 5. 结果（见 `results/fable_adc/tables/` 与 `figures/`）
- **TROP2（`TACSTD2`）在所有队列均与 ICI 结局无关**：应答对比 p=0.44（GSE126044）；MPR 对比 p=0.34（GSE207422）；
  PFS log-rank p=0.43、Cox HR/SD 1.07（95% CI 0.68–1.67）（GSE135222）。方向上与 EVOKE-01/02「疗效与 TROP2
  表达无关」一致。
- **CLDN18 在未经选择的 NSCLC 中几乎不表达**：GSE207422 中位 log2-TPM ≈ 0.07–0.08，MPR/NMPR 无差异（p=1.0）。
  这与临床事实吻合——肺癌 CLDN18.2 主要局限于浸润性黏液腺癌，未选择的 NSCLC 并非良好适应证。
- **TOP1**（载荷靶点）在 GSE135222 中呈「高表达→更差 PFS」的趋势（Cox HR/SD 1.54，p=0.10），未达显著。
- **PD-L1（`CD274`）** 表现符合 ICI 预期——应答者更高（GSE126044，rank-biserial −0.42）、PFS 保护趋势
  （GSE135222 Cox HR/SD 0.68，p=0.15），但样本量不足。
- **TROP2 与 PD-L1 无正相关**（三队列 Spearman ρ = −0.23 / +0.02 / −0.18），即 TROP2 表达携带的信息与 PD-L1
  大体正交——为 TROP2-ADC 与 ICI 联合（而非彼此冗余）提供生物学依据。

### 6. 解读与局限
公开转录组进一步支持 **TROP2 mRNA 不是免疫治疗结局的有用预测标志物**，这也解释了 ADC 领域为何转向空间/蛋白层面的
指标（QCS-NMR、H-score）。队列偏小（n=16–27），且无一接受 TROP2-ADC；`CLDN18` 基因层面无法区分 18.2 亚型。
结论为假设层面，描述的是靶点生物学而非 ADC 疗效。

---

### Reproduce
```bash
bash scripts/fable_adc/download_data.sh   # fetch open GEO files -> results/fable_adc/raw/
python3 scripts/fable_adc/analyze.py      # writes tables/, figures/, analysis_summary.json
```
Dependencies: `pandas numpy scipy statsmodels lifelines matplotlib openpyxl`.
