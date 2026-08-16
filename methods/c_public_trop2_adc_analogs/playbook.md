# Methods playbook — public TROP2-ADC analog scoring

**EN.** Methods only. No results, no numbers, no biological claims. This playbook describes how to score CLDN4 / tight-junction and IFN / MHC-I / APM transcripts after a public TROP2-ADC (or a declared TROP2-loss analog) versus its own control. It does not describe private SKB264 experiments.

**中文。** 仅方法。无结果、无数字、无生物学结论。本手册说明如何在公开 TROP2-ADC（或事先声明的 TROP2 缺失对照）相对其自身对照中，对 CLDN4 / 紧密连接与 IFN / MHC-I / APM 转录本打分。不描述私有 SKB264 实验。

---

## Scope / 范围

**EN.** Eligible series are public GEO (or equivalent) expression matrices in which (i) a TROP2-directed ADC is the treatment arm, or (ii) TACSTD2 / Trop2 is genetically lost, and both CLDN4 and at least the C4 IFN/MHC-I panel genes are measured. Label every non-lung series as a non-lung analog. Never rename sacituzumab govitecan (IMMU132 / SG) or datopotamab as SKB264 / sac-TMT.

**中文。** 合格系列须为公开 GEO（或同等）表达矩阵，且 (i) 处理臂为 TROP2 导向 ADC，或 (ii) TACSTD2 / Trop2 遗传缺失，并且同时测到 CLDN4 与至少 C4 IFN/MHC-I 面板基因。所有非肺系列必须标为非肺类似。不得把 sacituzumab govitecan（IMMU132 / SG）或 datopotamab 改称为 SKB264 / sac-TMT。

## Contrast / 对比

**EN.** Pre-specify one primary contrast per series: ADC (or KO) versus the depositor’s vehicle / parental / wild-type control. Combination arms (ADC + second drug) are secondary and are not used to rescue a null primary. Experimental unit is the library / tumor / PDX model as deposited. Do not concatenate series.

**中文。** 每个系列预先指定一个主对比：ADC（或 KO）对提交者的溶剂 / 亲本 / 野生型对照。联合用药臂为次要，不得用来挽救主对比的阴性。实验单位为所提交的文库 / 肿瘤 / PDX。不得拼接不同系列。

## Expression scale / 表达尺度

**EN.** Prefer author-processed counts or FPKM/RPKM as deposited. Analyse `log2(x + 1)`. Do not treat FPKM as DESeq2 counts. If only author DESeq2 tables exist, report those coefficients and also recompute a Welch test on any accompanying normalized counts.

**中文。** 优先使用提交者处理后的 counts 或 FPKM/RPKM。分析 `log2(x + 1)`。不得把 FPKM 当作 DESeq2 counts。若只有作者 DESeq2 表，报告其系数，并在附带的标准化 counts 上另做 Welch 检验。

## Gene sets / 基因集

**EN.** Freeze sets before looking at fold-changes:

1. Targets: TACSTD2, CLDN4.
2. C4 IFN/MHC-I panel: IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A.
3. APM and a broader IFN/ISG list.
4. Junction / barrier list with CLDN4 kept inside the set.

For mouse matrices, use a pre-specified 1:1 Ensembl map. Classical H2 genes stand in for HLA-A/B. Record symbols with no ortholog as absent. Do not build a composite “TJ index” as the primary CLDN4 readout.

**中文。** 看倍数变化前冻结基因集：靶基因 TACSTD2、CLDN4；C4 IFN/MHC-I 六基因面板；APM 与更宽 IFN/ISG；含 CLDN4 的连接/屏障列表。小鼠矩阵使用预先指定的 1:1 Ensembl 对照。经典 H2 基因代替 HLA-A/B。无直系同源的符号记为缺失。不得把复合 “TJ 指数” 当作 CLDN4 的主读出。

## Statistics / 统计

**EN.** Unpaired contrasts: Welch t-test on `log2(x+1)`. Paired PDX: paired t on the same scale; Wilcoxon signed-rank as a sensitivity. Genome-wide BH-FDR is exploratory. Gene-set test: Mann–Whitney of member log2FC versus expressed background, plus a two-sided binomial on the sign of the set. Sample scores are the mean of per-gene z-scores across libraries in that series. Always report n, log2FC, and p. Underpowered arms (n<3) are labeled underpowered, not negative.

**中文。** 非配对：对 `log2(x+1)` 做 Welch t。配对 PDX：同尺度配对 t，Wilcoxon 符号秩作敏感性。全基因组 BH-FDR 为探索性。基因集：成员 log2FC 对表达背景的 Mann–Whitney，外加集合符号的双侧二项。样本分数为该系列内逐基因 z 的均值。必须报告 n、log2FC 与 p。n<3 的臂标为效力不足，而不是阴性。

## Honesty rules / 诚实规则

**EN.** A series whose deposited samples were not treated with a TROP2-ADC cannot be scored as “after TROP2-ADC.” Knockout / shRNA series are TROP2-loss analogs. Papers that mention sacituzumab but deposit a different treatment (CHK1i, osimertinib, untreated PDX) are catalogued as cannot-test. Do not claim public analogs are SKB264.

**中文。** 所提交样本未经 TROP2-ADC 处理的系列，不得打成 “TROP2-ADC 之后”。敲除 / shRNA 系列是 TROP2 缺失类似。论文提到 sacituzumab 但提交的是其他处理（CHK1i、奥希替尼、未处理 PDX）的，目录记为无法检验。不得声称公开类似物就是 SKB264。
