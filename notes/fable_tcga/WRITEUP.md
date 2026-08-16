# TCGA LUAD / LUSC：TACSTD2 (TROP2) 与 CLDN4 的免疫微环境与生存分析
# TCGA LUAD / LUSC: TACSTD2 (TROP2) and CLDN4 vs Immune Contexture and Survival

> **重要声明 / KEY DISCLAIMER**
>
> **中文**：TCGA 队列入组于免疫检查点抑制剂（ICI）获批之前，**患者均未接受 ICI 治疗**。本分析中的任何"免疫相关性"都**不能**被解读为 ICI 疗效或获益的证据；它只描述未经免疫治疗的肿瘤中，靶点表达与免疫浸润估计值之间的横断面关联。
>
> **English**: TCGA patients were accrued before immune-checkpoint-inhibitor approvals and were **NOT treated with ICI**. No association reported here can be read as evidence of ICI response or benefit; these are cross-sectional associations between target expression and estimated immune infiltration in immunotherapy-naive tumors.

---

## 中文

### 一句话结论 (TL;DR)

在 **LUSC（肺鳞癌，n=501）** 中，TACSTD2/TROP2 高表达显著标记"免疫冷"肿瘤：与 T 细胞、CD8 T 细胞、细胞毒性淋巴细胞（MCP-counter）、T 细胞炎症 GEP（Ayers 18 基因）及几乎所有检查点基因均呈负相关（Spearman rho 约 −0.22 至 −0.31，BH-FDR 全部 < 0.05），且经 ABSOLUTE 肿瘤纯度校正后**结论不变**；唯一正相关的免疫成分是中性粒细胞（rho = +0.24）。在 **LUAD（肺腺癌，n=516）** 中，两个基因与免疫特征的关联总体很弱（|rho| ≤ 0.23），TACSTD2 甚至与 CD274 (PD-L1) 和髓系树突细胞呈弱正相关。**两个基因在两个队列中均不具有 OS 或 PFI 预后意义**（Cox 每标准差 HR 0.90–1.09，所有 p > 0.2；log-rank 所有 p ≥ 0.097）。两靶点表达中度共表达（LUAD rho = 0.53，LUSC rho = 0.39）。

### 数据来源（全部为公开开放矩阵，无受控数据）

| 文件 | 来源 | md5 |
|---|---|---|
| TCGA-LUAD.star_tpm.tsv.gz（log2(TPM+1)，GENCODE v36） | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_tpm.tsv.gz) | 73e6ce31c2ce12a8b63df11e53087eb4 |
| TCGA-LUSC.star_tpm.tsv.gz | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUSC.star_tpm.tsv.gz) | 8f35ed18f2497fb3242b69cf9f16b0b8 |
| TCGA-CDR 生存终点（OS/PFI 等；Liu et al., Cell 2018） | [Xena Pan-Cancer Atlas hub](https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp) | 50e4e056d3930bb98c51caf5f083d613 |
| gencode.v36 基因 probemap | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap) | 59d24b459af04b543cf1d5d4161a98fc |
| ABSOLUTE 肿瘤纯度（PanCanAtlas 开放补充文件） | [GDC open API](https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5) | — |
| MCP-counter 标记基因（Becht et al., Genome Biol 2016） | [作者 GitHub](https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt) | 944af31812c3f2aa61526074604e1c77 |

### 方法

- **样本**：仅原发肿瘤（样本码 -01），每位患者取一个样本：LUAD 516 例肿瘤 / 59 例癌旁正常；LUSC 501 / 51。全部肿瘤样本均可匹配 TCGA-CDR 临床与生存数据；ABSOLUTE 纯度覆盖 503 / 493 例。
- **免疫反卷积**：MCP-counter 10 个细胞群评分（标记基因 log2(TPM+1) 均值；上游签名文件中一处符号截断按其 Ensembl ID 解析为 CAVIN2，KIR3DS1 与 MGC40069 不在 GENCODE v36 矩阵内，已如实记录）；Ayers 等 18 基因 T 细胞炎症 GEP（每基因 z 分数取均值）；细胞溶解活性 CYT（GZMA、PRF1 的 log2 均值）；以及单基因检查点/效应分子（CD274、PDCD1、CTLA4、LAG3、TIGIT、HAVCR2、CD8A、GZMB、CXCL9 等）。
- **统计**：Spearman 相关 + 按队列×基因家族做 BH-FDR 校正；另做以 ABSOLUTE 纯度为协变量的秩偏相关（结果几乎不变——TACSTD2 与纯度本身几乎不相关：LUAD rho = 0.04，LUSC rho = 0.01）；四分位组间 Mann-Whitney U；肿瘤 vs 癌旁用非配对 Mann-Whitney U 与配对 Wilcoxon 双重检验。
- **生存**：OS 与 PFI（TCGA-CDR）；中位切分 Kaplan-Meier + log-rank；Cox 比例风险（表达按每标准差），未校正与校正（年龄、性别、AJCC III/IV vs I/II 期）两个模型（lifelines）。

### 结果

**1. 肿瘤 vs 癌旁（`tumor_vs_normal.csv`，图 `fig1`）**：两基因在肿瘤中均上调。最强信号：LUAD 中 CLDN4（中位数 8.18 vs 7.02 log2(TPM+1)，非配对 p = 6.6e-18，配对 p = 9.0e-9）；LUSC 中 TACSTD2（9.57 vs 8.70，p = 4.8e-9，配对 p = 4.0e-3）。LUAD 中 TACSTD2（Δ = +0.38，p = 8.0e-3）与 LUSC 中 CLDN4（Δ = +0.40，p = 0.014）上调幅度较小，配对检验未达显著。

**2. 免疫相关性（`immune_correlations.csv`，图 `fig2`、`fig5`）**

- **LUSC × TACSTD2 —— 最清晰的信号**：与 CXCL9（rho = −0.31；纯度校正后 −0.33）、细胞毒性淋巴细胞（−0.30）、CD8A（−0.28）、LAG3（−0.27）、CD8 T 细胞（−0.27）、PDCD1（−0.26）、CTLA4（−0.26）、GEP18（−0.25）、T 细胞（−0.24）、TIGIT（−0.22）、CYT（−0.20）均显著负相关（FDR 均 < 1e-5）；CD274 弱负相关（−0.09，FDR = 0.036）。唯一显著正相关：**中性粒细胞 +0.24**（校正后 +0.25）。四分位对比：TACSTD2 最高 vs 最低四分位的 CD8 T 细胞评分中位数 1.13 vs 1.86（p = 2.1e-7），GEP18 为 −0.21 vs +0.30（p = 2.5e-6）。
- **LUAD**：关联普遍弱。TACSTD2 与 GZMB（−0.17）、NK（−0.16）、B 系（−0.15）、CD8 T（−0.12）、CYT（−0.12）弱负相关，但与髓系树突细胞（+0.18）和 CD274（+0.13）弱**正**相关；GEP18 无关联（rho = 0.02）。CLDN4 方向类似（NK −0.23、GZMB −0.19、细胞毒性 −0.16）。
- **LUSC × CLDN4**：模式同向但显著更弱（中性粒细胞 +0.20；成纤维细胞 −0.18；GEP18 −0.08 无显著性）。
- 纯度校正后所有主要结论不变（见 `fig2` 右图）。

**3. 生存（`survival_cox.csv`、`survival_logrank.csv`，图 `fig3`、`fig4`）**：无任何显著关联。校正后 Cox（OS，每标准差）：LUAD TACSTD2 HR = 1.02（95% CI 0.87–1.20，p = 0.80）、CLDN4 HR = 1.00（0.85–1.17，p = 0.98）；LUSC TACSTD2 HR = 0.95（0.83–1.09，p = 0.45）、CLDN4 HR = 1.06（0.92–1.21，p = 0.43）。PFI 同样全部不显著。log-rank 最小 p 为 LUAD CLDN4 OS 的 0.097（高表达组趋势更差，未达显著）。

### 解读（谨慎版）

- 若考虑 TROP2-ADC 与 PD-(L)1 联合策略：在 LUSC 中，TROP2 mRNA 高的肿瘤平均而言 T 细胞浸润更少、GEP 更低——即靶点高表达人群与"免疫热"人群**并不重合**；在 LUAD 中两个维度近似独立。**但由于 TCGA 无 ICI 治疗，这不能推断联合或单药 ICI 的疗效方向。**
- 两基因 mRNA 均无独立预后价值，提示其表达高低本身在未经 ICI/ADC 治疗的自然病程中不是强预后因子。

### 局限

1. 非 ICI 队列（见顶部声明）；亦无 ADC 治疗数据。
2. Bulk RNA 反卷积是估计值而非计数；MCP-counter 分数只能做样本间比较，不能跨细胞群比较。
3. mRNA ≠ 蛋白：TROP2/CLDN4 作为 ADC 靶点的临床判断依赖 IHC 蛋白表达，本分析不能替代。
4. TCGA 以早中期手术切除标本为主，人群与当代晚期 NSCLC 试验人群不同。
5. 每患者单一样本，无法评估瘤内异质性；多重校正仅在每个队列×基因的特征家族内进行。

### 复现

```bash
pip install pandas numpy scipy statsmodels lifelines matplotlib
python3 scripts/fable_tcga/download_data.py   # 下载至 /tmp/fable_tcga_data（可用 FABLE_TCGA_DATA 覆盖）
python3 scripts/fable_tcga/run_analysis.py    # 输出至 results/fable_tcga/
```

---

## English

### TL;DR

In **LUSC (lung squamous, n=501)** high TACSTD2/TROP2 marks immunologically **cold** tumors: significant negative Spearman correlations with T cells, CD8 T cells, cytotoxic lymphocytes (MCP-counter), the Ayers 18-gene T-cell-inflamed GEP, and essentially every checkpoint gene (rho ≈ −0.22 to −0.31, all BH-FDR < 0.05), robust to adjustment for ABSOLUTE tumor purity; the only positively correlated population is neutrophils (rho = +0.24). In **LUAD (n=516)** associations are weak overall (|rho| ≤ 0.23), and TACSTD2 is even weakly *positively* correlated with CD274 (PD-L1) and myeloid dendritic cells. **Neither gene is prognostic for OS or PFI in either cohort** (Cox HR per SD 0.90–1.09, all p > 0.2; all log-rank p ≥ 0.097). The two targets are moderately co-expressed (rho = 0.53 LUAD, 0.39 LUSC).

### Data (all open-access matrices; no controlled-access data)

Same table as in the Chinese section above: Xena GDC hub STAR TPM (log2(TPM+1), GENCODE v36) for TCGA-LUAD and TCGA-LUSC; TCGA-CDR curated survival endpoints (Liu et al., Cell 2018) from the Xena Pan-Cancer Atlas hub; ABSOLUTE purity from the open PanCanAtlas supplement via the GDC API; MCP-counter marker genes from the authors' GitHub. MD5 checksums are listed in the table and printed by `download_data.py`.

### Methods

- **Samples**: primary tumors only (sample code -01), one per patient: LUAD 516 tumors / 59 adjacent normals; LUSC 501 / 51. All tumors matched to TCGA-CDR clinical data; ABSOLUTE purity available for 503 / 493.
- **Immune deconvolution**: MCP-counter scores for 10 populations (mean log2(TPM+1) of marker genes; one truncated symbol in the upstream signature file was resolved to CAVIN2 via its Ensembl ID; KIR3DS1 and MGC40069 are absent from the GENCODE v36 matrix and recorded as such); Ayers et al. 18-gene T-cell-inflamed GEP (mean of per-gene z-scores); cytolytic activity CYT (mean log2 of GZMA, PRF1); individual checkpoint/effector genes (CD274, PDCD1, CTLA4, LAG3, TIGIT, HAVCR2, CD8A, GZMB, CXCL9, …).
- **Statistics**: Spearman correlations with BH-FDR within each cohort × gene family; rank-based partial correlations adjusting for ABSOLUTE purity (results essentially unchanged — TACSTD2 itself is uncorrelated with purity: rho = 0.04 LUAD, 0.01 LUSC); Q4-vs-Q1 Mann-Whitney U contrasts; tumor-vs-normal tested with both unpaired Mann-Whitney U and paired Wilcoxon.
- **Survival**: OS and PFI (TCGA-CDR); median-split Kaplan-Meier with log-rank; Cox PH per SD of log2 TPM, unadjusted and adjusted for age, sex, and AJCC stage (III/IV vs I/II), via lifelines.

### Results

**1. Tumor vs adjacent normal** (`tumor_vs_normal.csv`, `fig1`): both genes are up in tumor. Strongest: CLDN4 in LUAD (median 8.18 vs 7.02 log2(TPM+1), unpaired p = 6.6e-18, paired p = 9.0e-9) and TACSTD2 in LUSC (9.57 vs 8.70, p = 4.8e-9, paired p = 4.0e-3). TACSTD2 in LUAD (Δ = +0.38, p = 8.0e-3) and CLDN4 in LUSC (Δ = +0.40, p = 0.014) are modest, with non-significant paired tests.

**2. Immune correlations** (`immune_correlations.csv`, `fig2`, `fig5`)

- **LUSC × TACSTD2 — the clearest signal**: negative with CXCL9 (rho = −0.31; −0.33 purity-adjusted), cytotoxic lymphocytes (−0.30), CD8A (−0.28), LAG3 (−0.27), CD8 T cells (−0.27), PDCD1 (−0.26), CTLA4 (−0.26), GEP18 (−0.25), T cells (−0.24), TIGIT (−0.22), CYT (−0.20), all FDR < 1e-5; CD274 weakly negative (−0.09, FDR = 0.036). The only significant positive population is **neutrophils, +0.24** (+0.25 adjusted). Quartile contrast: CD8 T-cell score median 1.13 (Q4) vs 1.86 (Q1), p = 2.1e-7; GEP18 −0.21 vs +0.30, p = 2.5e-6.
- **LUAD**: weak overall. TACSTD2 negative with GZMB (−0.17), NK (−0.16), B lineage (−0.15), CD8 T (−0.12), CYT (−0.12), but weakly **positive** with myeloid dendritic cells (+0.18) and CD274 (+0.13); GEP18 null (rho = 0.02). CLDN4 similar in direction (NK −0.23, GZMB −0.19, cytotoxic −0.16).
- **LUSC × CLDN4**: same direction as TACSTD2 but much weaker (neutrophils +0.20; fibroblasts −0.18; GEP18 −0.08, ns).
- All headline conclusions survive purity adjustment (right panel of `fig2`).

**3. Survival** (`survival_cox.csv`, `survival_logrank.csv`, `fig3`, `fig4`): nothing significant. Adjusted Cox for OS (per SD): LUAD TACSTD2 HR = 1.02 (95% CI 0.87–1.20, p = 0.80), CLDN4 HR = 1.00 (0.85–1.17, p = 0.98); LUSC TACSTD2 HR = 0.95 (0.83–1.09, p = 0.45), CLDN4 HR = 1.06 (0.92–1.21, p = 0.43). All PFI models also null. Smallest log-rank p is 0.097 (LUAD CLDN4, OS; trend toward worse survival in the high group, not significant).

### Interpretation (cautious)

- For TROP2-ADC ± PD-(L)1 combination thinking: in LUSC, TROP2-mRNA-high tumors have on average less T-cell infiltration and lower GEP, i.e. the target-high population does **not** coincide with the "immune-hot" population; in LUAD the two axes are approximately independent. **Because TCGA is not ICI-treated, this says nothing about ICI efficacy in either subgroup.**
- Neither gene's mRNA carries independent prognostic value here, so expression level per se is not a strong natural-history prognostic factor in ICI/ADC-naive resected NSCLC.

### Limitations

1. Not an ICI-treated cohort (see disclaimer); no ADC-treated data either.
2. Bulk-RNA deconvolution yields estimates, not counts; MCP-counter scores compare samples within a population, not across populations.
3. mRNA ≠ protein: clinical ADC target assessment relies on IHC; this analysis is not a substitute.
4. TCGA is dominated by resected, earlier-stage tumors and predates modern therapy; the population differs from contemporary advanced-NSCLC trial populations.
5. One sample per patient (no intratumoral heterogeneity); FDR correction applied within each cohort × gene feature family only.

### Reproduce

```bash
pip install pandas numpy scipy statsmodels lifelines matplotlib
python3 scripts/fable_tcga/download_data.py   # data land in /tmp/fable_tcga_data (override with FABLE_TCGA_DATA)
python3 scripts/fable_tcga/run_analysis.py    # writes results/fable_tcga/
```

### Output files

| Path | Content |
|---|---|
| `results/fable_tcga/tumor_vs_normal.csv` | Tumor vs adjacent-normal tests |
| `results/fable_tcga/immune_correlations.csv` | Spearman + purity-adjusted partial correlations, FDR |
| `results/fable_tcga/quartile_immune_comparison.csv` | Q4 vs Q1 immune-score contrasts |
| `results/fable_tcga/survival_logrank.csv` / `survival_cox.csv` | KM log-rank; Cox per-SD (unadjusted + adjusted) |
| `results/fable_tcga/cohort_summary.json` | Sample counts, purity correlations, TACSTD2–CLDN4 co-expression |
| `results/fable_tcga/fig1–fig5 (*.png)` | Boxplots, correlation heatmaps, KM curves (OS, PFI), scatter plots |
