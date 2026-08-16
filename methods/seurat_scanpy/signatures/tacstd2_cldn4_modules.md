# TACSTD2 / CLDN4 module signatures — curation record
# TACSTD2 / CLDN4 模块签名——构建记录

> **This file intentionally ships with EMPTY gene lists.** No signature is invented here.
> Fill in `signatures/modules.template.json` with genes you can justify, and record the
> provenance below. The analysis templates read your curated lists, not any list from this repo.
>
> **本文件有意保留空基因列表**，不编造任何签名。请在 `signatures/modules.template.json` 中填入
> 可辩护的基因，并在下方记录来源。分析模板读取你自己整理的列表，而非本仓库预置的任何列表。

## Why TACSTD2 and CLDN4 / 为什么是 TACSTD2 与 CLDN4

**EN** — `TACSTD2` (TROP2) and `CLDN4` (Claudin-4) are epithelial cell-surface proteins and clinically
relevant antibody-drug-conjugate (ADC) targets in lung / NSCLC. In the Human Lung Cell Atlas (HLCA v2)
annotation panels, `TACSTD2` appears as an epithelial / airway-epithelium marker, and `CLDN4` is a
tight-junction epithelial marker. A per-cell **module score** (rather than a single gene) gives a more
stable readout of the associated epithelial/tumor program.

**中文** — `TACSTD2`（TROP2）与 `CLDN4`（Claudin-4）是上皮细胞表面蛋白，也是肺癌/NSCLC 中临床相关的
抗体偶联药物（ADC）靶点。在人类肺细胞图谱 HLCA v2 的注释面板中，`TACSTD2` 是上皮/气道上皮标记，
`CLDN4` 是紧密连接上皮标记。使用每细胞**模块分数**（而非单一基因）能更稳定地反映相关上皮/肿瘤程序。

## How to build each module (choose a transparent route) / 如何构建模块（选择可追溯路径）

1. **Anchor + curated program / 锚基因 + 精选程序**
   Start from the anchor (`TACSTD2` or `CLDN4`) and add genes from a **named, versioned, citable**
   source — e.g. an MSigDB gene set (Hallmark / Reactome epithelial or tight-junction sets) or a
   published HLCA epithelial marker panel. Record the exact set name and version below.

2. **Data-driven co-expression / 数据驱动共表达**
   Derive genes co-expressed with the anchor **in your own data** (Spearman correlation on the
   epithelial subset, or an `hdWGCNA` / `Hotspot` module containing the anchor). Record the recipe
   (subset, method, threshold, date, object hash) so it is reproducible.

3. **Validate / 验证**
   Confirm the module (a) localizes to epithelial/malignant cells, (b) correlates with the anchor gene,
   (c) is stable across samples (score per-sample), and (d) does not merely track total counts
   (correlate score vs `nCount`/`total_counts`; regress out if needed).

## Provenance log (fill this in) / 来源记录（请填写）

| Module | Source type (1/2/3) | Exact source + version | Date | Curator | Notes |
|---|---|---|---|---|---|
| TACSTD2_module | | | | | |
| CLDN4_module | | | | | |

<!-- Do NOT paste unsourced gene lists here. Every gene must trace to a row above. -->
<!-- 请勿在此粘贴无来源的基因列表；每个基因都应能追溯到上表某一行。 -->
