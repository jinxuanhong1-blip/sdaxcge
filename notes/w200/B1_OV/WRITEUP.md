# w200/B1_OV — TCGA-OV：TACSTD2 (TROP2) 与 CLDN4 的免疫微环境与生存分析
# w200/B1_OV — TCGA-OV: TACSTD2 (TROP2) and CLDN4 vs Immune Contexture and Survival

B1（fable_tcga，TCGA LUAD/LUSC）流程的 TCGA-OV（卵巢浆液性囊腺癌）类比分析。
Analog of the B1 slice (fable_tcga, TCGA LUAD/LUSC) applied to TCGA-OV (ovarian serous cystadenocarcinoma).

> **重要声明 / KEY DISCLAIMER**
>
> **中文**：TCGA-OV 队列入组于免疫检查点抑制剂（ICI）获批之前，**患者均未接受 ICI 治疗**（一线为铂类/紫杉醇时代）。本分析中的任何"免疫相关性"都**不能**被解读为 ICI 疗效或获益的证据；它只描述未经免疫治疗的肿瘤中，靶点 mRNA 与免疫浸润估计值之间的横断面关联。
>
> **English**: TCGA-OV patients were accrued before immune-checkpoint-inhibitor approvals and were **NOT treated with ICI** (platinum/taxane era). No association reported here can be read as evidence of ICI response or benefit; these are cross-sectional associations between target mRNA and estimated immune infiltration in immunotherapy-naive tumors.

---

## 中文

### 一句话结论 (TL;DR)

在 **TCGA-OV（n=421 原发肿瘤）** 中，B1 在 LUSC 观察到的"TROP2 高 = 免疫冷"模式**没有重现**：TACSTD2 与各免疫特征的关联全部很弱（|rho| ≤ 0.23），且方向偏**正**——NK（+0.23）、中性粒细胞（+0.22）、HAVCR2（+0.19）、细胞毒性淋巴细胞（+0.17）、GEP18（+0.16）均为弱正相关（FDR < 0.05）；ABSOLUTE 纯度校正后仅细胞毒性淋巴细胞（+0.13）、HAVCR2（+0.15）、NK（+0.19）、中性粒细胞（+0.19）与 B 系（−0.19）保留显著。CLDN4 与免疫特征基本无关联（仅中性粒细胞 +0.17，校正后 +0.21）。与肺癌不同，**两个基因在 OV 中均显示适度的不良生存信号**：校正年龄与 FIGO 分期后，OS 每标准差 HR：TACSTD2 = 1.22（95% CI 1.06–1.40，p = 0.005）、CLDN4 = 1.20（1.04–1.37，p = 0.010）；中位切分 log-rank OS p 分别为 0.015 与 0.0023（高表达组更差）。PFI 仅 TACSTD2 校正后勉强显著（HR = 1.15，p = 0.031）。效应量均不大，且未跨模型做多重校正，应视为**待独立验证的适度信号**。两靶点共表达弱于肺癌（ρ = 0.32 vs LUAD 0.53 / LUSC 0.39）。在 Bausch-Fluck 2018 表面组宇宙中，CLDN4 是 TACSTD2 的 **#27 / 2,677** 共表达伙伴（n = 421，Spearman ρ = 0.320），**不是第 1**（第 1 为 TGFA，ρ = 0.470）。

### 与 B1（LUAD/LUSC）的关键差异

| 维度 | B1 LUSC | B1 LUAD | B1_OV（本分析） |
|---|---|---|---|
| TACSTD2 × T 细胞/CD8/GEP18 | 明确负相关（−0.22 至 −0.31） | 近似无关联 | 弱**正**相关（GEP18 +0.16，纯度校正后不显著） |
| 唯一稳定正相关成分 | 中性粒细胞 +0.24 | 髓系树突细胞 +0.18 | 中性粒细胞 +0.19–0.22（两基因均有，校正后仍显著） |
| OS/PFI 预后 | 全部无显著 | 全部无显著 | TACSTD2 与 CLDN4 高表达均与更差 OS 适度相关 |
| 肿瘤 vs 癌旁 | 可做（51–59 例癌旁） | 可做 | **不可做**：矩阵中无 -11 癌旁样本（如实记录，未伪造） |
| TACSTD2–CLDN4 表面组秩 | （B1 肺癌切片未做表面组全排名） | （同上） | **n=421，ρ=0.320，rank=#27/2677**（非第 1；#1=TGFA） |

### 数据来源（全部为公开开放矩阵，无受控数据；md5 由 download_data.py 打印）

| 文件 | 来源 | md5 |
|---|---|---|
| TCGA-OV.star_tpm.tsv.gz（log2(TPM+1)，GENCODE v36） | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-OV.star_tpm.tsv.gz) | 339cb343896ae6665ac14c1aea840667 |
| TCGA-CDR 生存终点（OS/PFI；Liu et al., Cell 2018） | [Xena Pan-Cancer Atlas hub](https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp) | 50e4e056d3930bb98c51caf5f083d613 |
| gencode.v36 基因 probemap | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap) | 59d24b459af04b543cf1d5d4161a98fc |
| ABSOLUTE 肿瘤纯度（PanCanAtlas 开放补充文件） | [GDC open API](https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5) | 8ea2ca92c8ae58350538999dfa1174da |
| MCP-counter 标记基因（Becht et al., Genome Biol 2016） | [作者 GitHub](https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt) | 944af31812c3f2aa61526074604e1c77 |
| in-silico human surfaceome table S3（Bausch-Fluck et al., PNAS 2018） | [steveneschrich/surfaceome 镜像](https://raw.githubusercontent.com/steveneschrich/surfaceome/main/data-raw/surfy/table_S3_surfaceome.xlsx) | bf95ed25ffc933fbe0c0383ab19b8652 |

### 方法（与 B1 相同，OV 特异性调整已注明）

- **样本**：仅原发肿瘤（-01 码），每位患者一个样本：421 例（矩阵含 422 个 -01 样本，1 例患者的重复取样只保留其一；7 个复发 -02 样本排除）。**无癌旁正常样本**（-11 计数为 0），故 B1 的肿瘤 vs 癌旁分析在 OV 不可行，已在 `cohort_summary.json` 中如实记录；`fig1` 改为仅展示肿瘤内两基因的表达分布。
- **免疫反卷积**：与 B1 完全一致——MCP-counter 10 个细胞群（标记基因 log2(TPM+1) 均值；KIR3DS1 与 MGC40069 不在 GENCODE v36 矩阵内，已记录）；Ayers 18 基因 GEP（每基因 z 分数均值）；CYT（GZMA、PRF1 log2 均值）；单基因检查点/效应分子。
- **统计**：Spearman + 每基因家族内 BH-FDR；以 ABSOLUTE 纯度（409/421 例可得）为协变量的秩偏相关；Q4 vs Q1 Mann-Whitney U。
- **生存**：OS 与 PFI（TCGA-CDR；n=419 有有效随访时间）；中位切分 KM + log-rank；Cox（每标准差），未校正与校正两个模型。**OV 特异性调整**：TCGA-OV 的 `ajcc_pathologic_tumor_stage` 全部缺失，分期改用 FIGO `clinical_stage`（III/IV vs I/II，剥离 A/B/C 亚期）；队列全部为女性，故去掉 B1 的性别协变量。校正模型 = 年龄 + FIGO III/IV。
- **表面组共表达排名**（与 B1_BRCA 类比一致）：以 TACSTD2 为锚点，对 Bausch-Fluck 2018 table S3「in silico surfaceome only」中、本矩阵可测且方差非零的全部表面基因计算 Spearman，报告 CLDN4 的 n / ρ / rank。未预先圈定有利邻域。Pearson 仅作敏感性。重复符号按最高均值行折叠（PAR / 多映射 GENCODE ID）。

### 结果

**1. 肿瘤 vs 癌旁**：**未执行**——TCGA-OV RNA-seq 矩阵不含任何 -11 癌旁样本。这是数据事实而非分析选择；未用 GTEx 等外部正常组织替代（跨项目批次效应会使该对比不可靠）。

**2. 免疫相关性（`immune_correlations.csv`，图 `fig2`、`fig5`）**

- **TACSTD2**：所有关联均弱（|rho| ≤ 0.23）。边际显著正相关：NK +0.23、中性粒细胞 +0.22、HAVCR2 +0.19、细胞毒性淋巴细胞 +0.17、GEP18 +0.16、CYT +0.15、单核系 +0.14、T 细胞 +0.13、CD274 +0.12（FDR 均 < 0.05）；B 系为唯一负相关（−0.11）。注意 TACSTD2 与 ABSOLUTE 纯度弱负相关（rho = −0.12，p = 0.014），纯度校正后 GEP18/CYT/T 细胞/CD274 均降至不显著，仅细胞毒性淋巴细胞（+0.13）、HAVCR2（+0.15）、NK（+0.19）、中性粒细胞（+0.19）、B 系（−0.19）保留 FDR < 0.05。**结论：OV 中 TROP2 高既不标记"免疫冷"，也谈不上强"免疫热"，只是弱的偏热趋势。**
- **CLDN4**：基本无免疫关联。仅中性粒细胞 +0.17（校正后 +0.21）稳定显著；CD8 T 细胞边际弱负（−0.14，FDR = 0.044），校正后不显著。
- 四分位对比（`quartile_immune_comparison.csv`）：TACSTD2 Q4 vs Q1 的细胞毒性淋巴细胞中位数 1.00 vs 0.86（p = 0.013）、GEP18 +0.12 vs −0.13（p = 0.029）；CLDN4 Q4 vs Q1 的 CD8 T 细胞 0.60 vs 0.80（p = 0.004，方向为负）。效应量都小。

**3. 生存（`survival_cox.csv`、`survival_logrank.csv`，图 `fig3`、`fig4`）**

- **OS**：TACSTD2 log-rank p = 0.015；Cox 未校正 HR = 1.12/SD（p = 0.104），校正年龄+分期后 HR = 1.22（95% CI 1.06–1.40，p = 0.005）。CLDN4 log-rank p = 0.0023；未校正 HR = 1.16（p = 0.027），校正后 HR = 1.20（1.04–1.37，p = 0.010）。两基因高表达均与更差 OS 相关。
- **PFI**：TACSTD2 log-rank p = 0.040，校正 HR = 1.15（1.01–1.30，p = 0.031）；CLDN4 全部不显著（校正 HR = 1.05，p = 0.42）。
- 诚实提示：HR ≈ 1.15–1.22/SD 属适度效应；本切片共 8 个 Cox 模型与 4 个 log-rank 未做跨检验校正；未校正与校正模型间的差异（TACSTD2）提示与年龄/分期存在一定混杂结构，需独立队列验证。

**4. TACSTD2–CLDN4 表面组共表达排名（`surface_rank_headline.csv`、`surface_rank.md`、`coexpression_TACSTD2_surfaceome.csv`，图 `fig6`、`fig7`）**

| | |
|---|---|
| n | **421** 原发肿瘤（-01，每患者 1 样本） |
| ρ（Spearman） | **0.320**（p = 1.64e-11；BH-FDR q = 1.51e-09） |
| rank | **#27 / 2,677**（第 99.029 百分位） |
| 是否第 1 | **否**。第 1 为 TGFA（ρ = 0.470） |
| Pearson 敏感性 | r = 0.407，rank #10 / 2,677 |

宇宙构成：table S3 去重 UniProt 基因 2,799 个；113 个不在本矩阵 GENCODE v36 符号中；8 个零方差无法定义相关（CLDN22、GP1BB、ICAM4、LRRC24、OR14I1、OR4N4、OR5M9、SERINC4），已如实剔除。CLDN4 在表面组中属**高位伙伴（约前 1%）但不是第一**。同一对在 B1_BRCA 类比为 #4 / 2,618（ρ = 0.348）。此处 ρ 与 `cohort_summary.json` 中成对 Spearman（0.32）一致，未另造数字。

**5. 其他**：CLDN4 与纯度 rho = +0.09（p = 0.065）。

### 解读（谨慎版）

- 对 TROP2-ADC（如 sacituzumab govitecan 类）在卵巢癌的联合策略思考：与 LUSC 不同，OV 中 TROP2-mRNA 高的肿瘤**并不**更"冷"；靶点表达与免疫浸润近似正交或弱正相关。**但 TCGA 无 ICI/ADC 治疗，本数据不能推断任何疗法的疗效方向。**
- 两基因 mRNA 高表达与更差 OS 的适度关联（校正后仍在）值得在独立 OV 队列（如 AOCS/GSE9891、TothillOV）中验证后再作解读；以当前证据不足以称其为预后标志物。
- TACSTD2–CLDN4 共表达在 OV 表面组中排第 27（ρ = 0.320），不能表述为"最强表面伙伴"或与 BRCA（#4）等强。

### 局限

1. 非 ICI/ADC 治疗队列（见顶部声明）。
2. **无癌旁正常对照**——肿瘤 vs 正常上调无法在本队列评估。
3. Bulk 反卷积是估计值；MCP-counter 分数只能样本间比较；OV 高纯度/高间质异质性使反卷积更不稳定，纯度校正已给出但只是线性近似。
4. mRNA ≠ 蛋白；ADC 靶点判断依赖 IHC。
5. TCGA-OV 几乎全为晚期（FIGO III/IV 393/418 有分期者）高级别浆液性癌，分期协变量变异度小（I/II 仅 25 例），"校正分期"作用有限；铂类时代人群与当代试验人群不同。
6. 生存分析未跨模型做多重校正；每患者单一样本。
7. 表面组宇宙依赖 2018 年 in-silico 预测名单与 UniProt 基因名 ↔ GENCODE v36 符号的字符串匹配；未做同义词（如 PVRL4/NECTIN4）手工合并，故部分伙伴可能因别名未入宇宙。零方差基因已剔除而非填 0。

### 复现

```bash
pip install pandas numpy scipy statsmodels lifelines matplotlib openpyxl
python3 scripts/w200/B1_OV/download_data.py   # 下载至 /tmp/w200_b1_ov_data（可用 W200_B1_OV_DATA 覆盖）
python3 scripts/w200/B1_OV/run_analysis.py                 # 免疫 + 生存 → results/w200/B1_OV/
python3 scripts/w200/B1_OV/rank_surface_coexpression.py    # 表面组 n/ρ/rank
```

---

## English

### TL;DR

In **TCGA-OV (n=421 primary tumors)** the "TROP2-high = immune-cold" pattern that B1 found in LUSC does **not** replicate: all immune associations of TACSTD2 are weak (|rho| ≤ 0.23) and skew **positive** — NK cells (+0.23), neutrophils (+0.22), HAVCR2 (+0.19), cytotoxic lymphocytes (+0.17), GEP18 (+0.16), all BH-FDR < 0.05; after ABSOLUTE-purity adjustment only cytotoxic lymphocytes (+0.13), HAVCR2 (+0.15), NK (+0.19), neutrophils (+0.19) and B lineage (−0.19) remain significant. CLDN4 is essentially uncorrelated with immune features (only neutrophils, +0.17 marginal / +0.21 adjusted). Unlike lung, **both genes show a modest adverse survival signal in OV**: adjusted for age and FIGO stage, OS HR per SD is 1.22 for TACSTD2 (95% CI 1.06–1.40, p = 0.005) and 1.20 for CLDN4 (1.04–1.37, p = 0.010); median-split log-rank OS p = 0.015 and 0.0023 respectively (high = worse). For PFI only TACSTD2 is marginally significant after adjustment (HR = 1.15, p = 0.031). Effects are modest and not corrected across models — treat as **a moderate signal pending independent validation**. Target co-expression is weaker than in lung (ρ = 0.32 vs 0.53 LUAD / 0.39 LUSC). Among the Bausch-Fluck 2018 surfaceome, CLDN4 is TACSTD2's **#27 / 2,677** co-expression partner (n = 421, Spearman ρ = 0.320) — **not #1** (#1 is TGFA, ρ = 0.470).

### Key differences vs B1 (LUAD/LUSC)

| Axis | B1 LUSC | B1 LUAD | B1_OV (this slice) |
|---|---|---|---|
| TACSTD2 × T cells/CD8/GEP18 | clearly negative (−0.22 to −0.31) | ~null | weakly **positive** (GEP18 +0.16; ns after purity adjustment) |
| Only robust positive population | neutrophils +0.24 | myeloid DC +0.18 | neutrophils +0.19–0.22 (both genes, survives adjustment) |
| OS/PFI prognosis | all null | all null | high TACSTD2 and CLDN4 both modestly associated with worse OS |
| Tumor vs adjacent normal | feasible (51–59 normals) | feasible | **infeasible**: no -11 samples in the matrix (recorded honestly, not fabricated) |
| TACSTD2–CLDN4 surfaceome rank | (B1 lung slice did not rank the full surfaceome) | (same) | **n=421, ρ=0.320, rank=#27/2677** (not #1; #1=TGFA) |

### Data (all open-access matrices; no controlled-access data)

Same table as the Chinese section: Xena GDC hub STAR TPM for TCGA-OV (log2(TPM+1), GENCODE v36); TCGA-CDR curated endpoints (Liu et al., Cell 2018); GENCODE v36 probemap; ABSOLUTE purity from the open PanCanAtlas supplement via the GDC API; MCP-counter marker genes from the authors' GitHub; Bausch-Fluck 2018 in-silico surfaceome table S3 (steveneschrich/surfaceome mirror). MD5 checksums are in the table above and printed by `download_data.py`.

### Methods (identical to B1; OV-specific adaptations flagged)

- **Samples**: primary tumors only (-01), one per patient: 421 (the matrix has 422 -01 samples; one duplicate vial per patient collapsed; 7 recurrent -02 samples excluded). **No adjacent normals exist** (zero -11 samples), so the B1 tumor-vs-normal arm is impossible here; this is recorded in `cohort_summary.json`, and `fig1` instead shows tumor-only expression distributions. No external normal tissue (e.g. GTEx) was substituted, since cross-project batch effects would make that comparison unreliable.
- **Immune deconvolution**: exactly as B1 — MCP-counter 10 populations (mean log2(TPM+1) of markers; KIR3DS1 and MGC40069 absent from the GENCODE v36 matrix, recorded); Ayers 18-gene GEP (mean of per-gene z-scores); CYT (mean log2 of GZMA, PRF1); individual checkpoint/effector genes.
- **Statistics**: Spearman with BH-FDR within each gene family; rank-based partial correlations adjusting for ABSOLUTE purity (available for 409/421); Q4-vs-Q1 Mann-Whitney U.
- **Survival**: OS and PFI (TCGA-CDR; n=419 with valid follow-up time); median-split KM + log-rank; Cox per SD, unadjusted and adjusted. **OV adaptations**: `ajcc_pathologic_tumor_stage` is entirely missing for OV, so FIGO `clinical_stage` is used (III/IV vs I/II, stripping A/B/C substages); the cohort is all-female, so B1's sex covariate is dropped. Adjusted model = age + FIGO III/IV.
- **Surfaceome co-expression rank** (same definition as the B1_BRCA analog): TACSTD2 is the anchor; every gene on the Bausch-Fluck 2018 table S3 "in silico surfaceome only" sheet that is measured and non-zero-variance in this matrix is scored by Spearman; CLDN4's n / ρ / rank is reported. No favourable neighbourhood was pre-selected. Pearson is a sensitivity only. Duplicate symbols collapsed to the highest-mean row (PAR / multi-mapped GENCODE IDs).

### Results

**1. Tumor vs adjacent normal**: **not performed** — the TCGA-OV RNA-seq matrix contains no solid-tissue-normal samples. This is a property of the data, not an analysis choice.

**2. Immune correlations** (`immune_correlations.csv`, `fig2`, `fig5`)

- **TACSTD2**: everything is weak (|rho| ≤ 0.23). Marginally significant positives: NK +0.23, neutrophils +0.22, HAVCR2 +0.19, cytotoxic lymphocytes +0.17, GEP18 +0.16, CYT +0.15, monocytic lineage +0.14, T cells +0.13, CD274 +0.12 (all FDR < 0.05); B lineage is the only negative (−0.11). TACSTD2 is weakly negatively correlated with purity itself (rho = −0.12, p = 0.014); after purity adjustment GEP18/CYT/T cells/CD274 drop to non-significance, and only cytotoxic lymphocytes (+0.13), HAVCR2 (+0.15), NK (+0.19), neutrophils (+0.19) and B lineage (−0.19) keep FDR < 0.05. **Bottom line: in OV, TROP2-high does not mark cold tumors, and the "hot" trend is weak.**
- **CLDN4**: essentially null. Only neutrophils +0.17 (adjusted +0.21) is robust; CD8 T cells weakly negative (−0.14, FDR = 0.044), not significant after adjustment.
- Quartile contrasts (`quartile_immune_comparison.csv`): TACSTD2 Q4 vs Q1 cytotoxic-lymphocyte median 1.00 vs 0.86 (p = 0.013), GEP18 +0.12 vs −0.13 (p = 0.029); CLDN4 Q4 vs Q1 CD8 T cells 0.60 vs 0.80 (p = 0.004, negative direction). All small effects.

**3. Survival** (`survival_cox.csv`, `survival_logrank.csv`, `fig3`, `fig4`)

- **OS**: TACSTD2 log-rank p = 0.015; Cox unadjusted HR = 1.12/SD (p = 0.104), adjusted HR = 1.22 (95% CI 1.06–1.40, p = 0.005). CLDN4 log-rank p = 0.0023; unadjusted HR = 1.16 (p = 0.027), adjusted HR = 1.20 (1.04–1.37, p = 0.010). High expression of either gene tracks with worse OS.
- **PFI**: TACSTD2 log-rank p = 0.040, adjusted HR = 1.15 (1.01–1.30, p = 0.031); CLDN4 null (adjusted HR = 1.05, p = 0.42).
- Honesty notes: HR ≈ 1.15–1.22 per SD is a modest effect; 8 Cox models and 4 log-rank tests were run in this slice with no cross-test correction; the unadjusted-vs-adjusted shift for TACSTD2 indicates some confounding structure with age/stage. Independent validation is required.

**4. TACSTD2–CLDN4 surfaceome rank** (`surface_rank_headline.csv`, `surface_rank.md`, `coexpression_TACSTD2_surfaceome.csv`, `fig6`, `fig7`)

| | |
|---|---|
| n | **421** primary tumors (-01, one per patient) |
| ρ (Spearman) | **0.320** (p = 1.64e-11; BH-FDR q = 1.51e-09) |
| rank | **#27 / 2,677** (99.029th percentile) |
| Is #1? | **No.** #1 is TGFA (ρ = 0.470) |
| Pearson sensitivity | r = 0.407, rank #10 / 2,677 |

Universe construction: 2,799 unique UniProt gene names on table S3; 113 absent from GENCODE v36 symbols in this matrix; 8 zero-variance (rho undefined) dropped and listed (CLDN22, GP1BB, ICAM4, LRRC24, OR14I1, OR4N4, OR5M9, SERINC4). CLDN4 is a **high** partner (~top 1%) but not the top one. The same pair ranked #4 / 2,618 in the B1_BRCA analog (ρ = 0.348). The ρ here is the same pairwise Spearman already stored in `cohort_summary.json` (0.32); no number was invented.

**5. Other**: CLDN4 vs purity rho = +0.09 (p = 0.065).

### Interpretation (cautious)

- For TROP2-ADC ± PD-(L)1 combination thinking in ovarian cancer: unlike LUSC, TROP2-mRNA-high OV tumors are **not** colder; target expression and immune infiltration are roughly orthogonal to weakly positive. **TCGA has no ICI or ADC exposure, so nothing here speaks to efficacy of either modality.**
- The modest adverse-OS association of both genes (surviving covariate adjustment) is a hypothesis for validation in independent OV cohorts (e.g. AOCS/Tothill GSE9891), not an established prognostic marker.
- TACSTD2–CLDN4 co-expression ranks #27 of the OV surfaceome (ρ = 0.320). That is not "the strongest surface partner" and is weaker than the B1_BRCA analog (#4).

### Limitations

1. Not an ICI/ADC-treated cohort (see disclaimer).
2. **No adjacent-normal control** — tumor-vs-normal upregulation cannot be assessed in this cohort.
3. Bulk deconvolution yields estimates; MCP-counter scores compare samples, not populations; OV's high purity/stromal heterogeneity makes deconvolution less stable, and linear purity adjustment is only an approximation.
4. mRNA ≠ protein; clinical ADC target calls rely on IHC.
5. TCGA-OV is almost entirely advanced-stage (FIGO III/IV in 393/418 staged) high-grade serous carcinoma; the stage covariate has little variance (only 25 stage I/II), so "stage-adjusted" is weak; the platinum-era population differs from contemporary trial populations.
6. No multiplicity correction across survival models; one sample per patient.
7. The surfaceome universe is a 2018 in-silico prediction list matched by UniProt gene name ↔ GENCODE v36 symbol string; synonyms (e.g. PVRL4/NECTIN4) were not manually merged, so some partners may be missing from the universe. Zero-variance genes were dropped, not imputed.

### Reproduce

```bash
pip install pandas numpy scipy statsmodels lifelines matplotlib openpyxl
python3 scripts/w200/B1_OV/download_data.py   # data land in /tmp/w200_b1_ov_data (override with W200_B1_OV_DATA)
python3 scripts/w200/B1_OV/run_analysis.py                 # immune + survival → results/w200/B1_OV/
python3 scripts/w200/B1_OV/rank_surface_coexpression.py    # surfaceome n / ρ / rank
```
