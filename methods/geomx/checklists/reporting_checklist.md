# Reporting checklist — GeoMx WTA compartment / mixed-effects / OS
# 报告清单 — GeoMx WTA 分区 / 混合效应 / 总生存

Tick every item before submission. Items map to playbook sections.
投稿前逐条勾选。条目对应手册章节。

## Design and units / 设计与分析单元

- [ ] Independent unit stated as **patient**, not AOI. `n` in every results sentence is patients (or events). / 独立单元写明为**患者**而非 AOI；结果句中的 `n` 是患者数或事件数。
- [ ] Nesting drawn: slide → patient → core → ROI → segment. / 画出嵌套：玻片 → 患者 → core → ROI → 分区。
- [ ] AOIs per patient per segment tabulated (balance). / 按患者×分区列出 AOI 数（平衡性）。
- [ ] Slide × cohort and slide × outcome cross-tabs shown; perfect confounding declared if present. / 给出玻片×队列、玻片×结局交叉表；若完全混杂须声明。
- [ ] Public data accession listed (GSE271689 / GSE292098); no restricted counts. / 列出公开数据登录号；无受限计数。

## QC / 质控

- [ ] Per-segment nuclei and area floors, not one threshold for all compartments. / 分区设置细胞核与面积下限，而非全部分区共用一个阈值。
- [ ] AOIs in → out at each QC step, **split by segment**. / 每步质控的 AOI 进出数，**按分区拆分**。
- [ ] LOQ formula and constants (`n = 2`, `minLOQ = 2`) reported. / 报告 LOQ 公式与常数。
- [ ] Gene filter described as union (within-compartment) vs intersection (cross-compartment). / 写明基因过滤用的是并集（分区内）还是交集（跨分区）。
- [ ] Median genes detected, area and nuclei per segment. / 各分区检出基因、面积、细胞核的中位数。

## Normalization / 归一化

- [ ] Primary method **declared before DE**. / 主方法在看 DE **之前**指定。
- [ ] At least one alternative shown (RLE + PCA by segment / slide / area / detection). / 至少展示一种备选（RLE + 按分区/玻片/面积/检出率着色的 PCA）。
- [ ] Joint vs per-segment normalization stated, and matched to the contrast type. / 写明合并还是分区归一化，并与对比类型匹配。
- [ ] Batch-corrected matrix used only for visualization. / 批次校正矩阵仅用于可视化。

## Mixed models / 混合模型

- [ ] Patient random effect (or patient×segment pseudobulk) in every AOI-level model. / 每一个 AOI 层面模型都有患者随机效应（或患者×分区 pseudobulk）。
- [ ] Random slope used for co-existing segments; random intercept only for between-patient contrasts. / 共存分区用随机斜率；患者间对比只用随机截距。
- [ ] `consensus.correlation` reported; second round if > 0.5. / 报告 `consensus.correlation`；> 0.5 时跑第二轮。
- [ ] Singular-fit rate reported if `lmer` was used genome-wide. / 若全基因组使用 `lmer`，报告奇异拟合比例。
- [ ] Satterthwaite or Kenward–Roger df for small `n`. / 小样本使用 Satterthwaite 或 Kenward–Roger 自由度。
- [ ] BH **within contrast**; number of contrasts stated. / BH 在**对比内部**进行；写明对比个数。
- [ ] Pseudobulk sensitivity run; direction clashes resolved toward pseudobulk. / 做了 pseudobulk 敏感性分析；方向冲突以 pseudobulk 为准。

## TACSTD2 / CLDN4

- [ ] Both genes declared prespecified. / 两基因均声明为事先指定。
- [ ] Q1 Level: % AOIs above LOQ per segment. / Q1 水平：各分区高于 LOQ 的 AOI 比例。
- [ ] Q2 Specificity: mixed-model LFC + 95% CI + % patients with LFC > 0. / Q2 特异性：混合模型 LFC + 95% CI + LFC > 0 的患者比例。
- [ ] Q3 Heterogeneity: ICC, reliability `R`, ROIs needed for `R ≥ 0.8`. / Q3 异质性：ICC、可靠性 `R`、达到 `R ≥ 0.8` 所需 ROI 数。
- [ ] Epithelial spillover index (genes listed) correlated with immune-AOI target expression. / 上皮渗漏指数（列出基因）与免疫 AOI 靶点表达的相关。
- [ ] Estimate repeated in the lowest spillover tertile. / 在最低渗漏三分位中重复估计。
- [ ] Conclusion phrased as a **lower bound** on compartment specificity. / 结论表述为分区特异性的**下界**。
- [ ] No claim that immune cells "express TACSTD2/CLDN4" without a spillover control. / 没有在缺乏渗漏对照的情况下声称免疫细胞"表达 TACSTD2/CLDN4"。

## Survival / 生存

- [ ] Primary model is **one row per patient**. / 主模型为**每位患者一行**。
- [ ] Exposure continuous, HR per 1 SD; no primary dichotomization. / 暴露为连续变量，HR 为每 1 个 SD；主分析不分二分类。
- [ ] Events, patients, events-per-variable stated. / 写明事件数、患者数、每变量事件数。
- [ ] Reliability `R` next to every HR. / 每个 HR 旁有可靠性 `R`。
- [ ] PH (`cox.zph`) and functional form (spline) checked. / 检查了比例风险与函数形式。
- [ ] Cohort/slide stratified, not naively adjusted, when baseline hazards differ. / 基线风险不同时对队列/玻片分层而非简单校正。
- [ ] AOI-level clustered / frailty Cox labelled **sensitivity**. / AOI 层面聚类/脆弱 Cox 标明为**敏感性分析**。
- [ ] If a cut point exists: locked in training, applied unchanged to validation; multiplicity/optimism handled. / 若有切点：在训练集锁定、原样用于验证集；处理了多重性/乐观偏倚。
- [ ] KM has numbers at risk. / KM 曲线带风险人数。
- [ ] Genome-wide OS scan, if any, labelled hypothesis-generating. / 若有全基因组 OS 扫描，标明为假设生成。
- [ ] Apparent C-index of a score fit on the same data is **not** reported as performance. / 未把同一数据上拟合评分的表观 C-index 当作性能。

## Reproducibility / 可重复性

- [ ] Config file (`study_config.yml`) archived with the paper. / 配置文件随论文存档。
- [ ] Session info (`*_sessionInfo.txt`) in the supplement. / 补充材料含会话信息。
- [ ] Software versions listed. / 列出软件版本。
- [ ] Random seed recorded. / 记录随机种子。
- [ ] Public GEO accessions and PKC version listed. / 列出公开 GEO 登录号与 PKC 版本。
