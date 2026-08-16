# KEYNOTE lung processed expression data / KEYNOTE 肺癌处理后表达数据

**Search date / 检索日期:** 2026-08-16  
**Scope / 范围:** named KEYNOTE lung-cancer trial specimens; RNA-seq, NanoString, or gene-by-sample counts; anonymous public access; files under 2 GB.

## Bottom line / 结论

**EN.** I found **no dataset satisfying all four requirements**: (1) a named KEYNOTE trial, (2) a lung-cancer cohort, (3) a processed patient-level RNA/NanoString/count matrix, and (4) anonymous open download. Therefore, no biological matrix was downloaded. This is a zero-result, not a failed download. The sponsor publications consistently route genetic or exploratory biomarker data through application, scientific review, a data-sharing agreement, and often a secure analysis environment.

The closest literal GEO hit is not a lung KEYNOTE dataset. GSE216069/GSE216055 contain open processed single-nucleus RNA and Slide-seq counts from three KEYNOTE-001 biopsies, but the paper explicitly identifies that patient as **melanoma**. The same deposit contains separate NSCLC method-validation specimens; those are not the KEYNOTE-001 specimens. Mixing the two would create a false “KEYNOTE lung” dataset.

**中文。** 未发现同时满足以下四项的数据集：(1) 明确命名的 KEYNOTE 试验；(2) 肺癌队列；(3) 患者级处理后 RNA/NanoString/计数矩阵；(4) 无需申请即可公开下载。因此本次没有下载生物学矩阵。这是诚实的“零结果”，不是下载失败。相关申办方论文均说明：遗传或探索性生物标志物数据通常需要提交研究方案、接受科学审查、签署数据共享协议，并常常只能在安全分析环境中使用。

最接近的 GEO 字面命中并不是 KEYNOTE 肺癌数据。GSE216069/GSE216055 确实公开了来自 KEYNOTE-001 三次活检的单核 RNA 和 Slide-seq 处理后计数，但论文明确说明该患者患有**黑色素瘤**。同一存档中的 NSCLC 样本只是另外的技术验证样本，并非 KEYNOTE-001 样本；把两者混在一起会错误地制造出“KEYNOTE 肺癌”数据集。

## What exists / 已找到的内容

The machine-readable, evidence-linked catalog is [`results/gpt_keynote/catalog.tsv`](../../results/gpt_keynote/catalog.tsv).

| Trial/resource | Lung | Expression evidence | Open processed patient matrix | Assessment |
|---|---:|---|---:|---|
| KEYNOTE-001 NSCLC | yes | pivotal report is PD-L1 IHC; later sponsor biomarker work | no | controlled/on-request |
| GSE216069/GSE216055 | no for KEYNOTE samples | snRNA-seq + Slide-seq V2 counts | yes | excluded: KEYNOTE patient is melanoma |
| KEYNOTE-010 | yes | samples referenced by sponsor biomarker analyses | no | controlled/on-request |
| KEYNOTE-028 SCLC | yes | NanoString-derived T-cell-inflamed GEP + TMB | no | only aggregate publication results |
| KEYNOTE-042 | yes | PD-L1 IHC + WES/TMB | no RNA matrix | controlled/on-request |
| KEYNOTE-189/407 | yes | WES/TMB + selected mutations | no RNA matrix | controlled/on-request |
| KEYNOTE-495 | yes | tumor RNA-derived 18-gene GEP + TMB | no | aggregate groups/scores only |
| KEYNOTE-782 | yes | bulk tumor RNA-seq, 69 evaluable, 11 signatures | no | abstract only; no accession/matrix |
| SU2C-MARK | yes | open processed bulk RNA-seq TPM | yes | excluded: not a KEYNOTE trial |

中文概述：KEYNOTE-028、495 和 782 明确使用了表达相关检测，但公开资料只给出签名分数、分组或汇总统计，没有患者级基因表达矩阵。KEYNOTE-042、189、407 的已发表探索性分子分析主要是 WES/TMB，不应误标为 RNA 数据。SU2C-MARK 有公开的肺癌 RNA-seq TPM 数据，但它不是 KEYNOTE 试验，因此只作为近似资源列入目录，未下载。

## TACSTD2 and CLDN4 / TACSTD2 与 CLDN4

**EN.**

- Neither **TACSTD2** nor **CLDN4** is part of the published 18-gene T-cell-inflamed GEP used in the KEYNOTE NanoString-derived analyses.
- KEYNOTE-782 used whole-transcriptome RNA-seq, so both genes were in principle measurable, but patient-level values are not public and cannot be checked.
- Both genes should be represented by the whole-transcriptome SU2C-MARK assay, but the 1.18 GB archive was not inspected because that cohort is not KEYNOTE; this is not reported as a verified presence.
- No claim of absence was made when the matrix itself was unavailable; the catalog says “not assessable.”

**中文。**

- **TACSTD2** 和 **CLDN4** 均不属于已发表的 18 基因 T-cell-inflamed GEP，因此不会出现在相应 KEYNOTE NanoString 衍生签名中。
- KEYNOTE-782 使用全转录组 RNA-seq，理论上可测到这两个基因，但患者级数据未公开，无法核查数值。
- SU2C-MARK 使用全转录组检测，理论上应覆盖这两个基因；但由于该队列并非 KEYNOTE，本次未检查其 1.18 GB 存档，也不把它写成“已验证存在”。
- 当原始矩阵不可获得时，本报告不武断地写“缺失”，而是标记为“无法评估”。

Detailed per-resource calls are in [`results/gpt_keynote/tacstd2_cldn4.tsv`](../../results/gpt_keynote/tacstd2_cldn4.tsv).

## Downloads and reproducibility / 下载与复现

`results/gpt_keynote/download_manifest.tsv` is intentionally header-only: there were zero eligible biological downloads. Open files above the limit were not downloaded; controlled files were not accessed; open near-matches outside the KEYNOTE-lung intersection were not downloaded.

Run:

```bash
python3 scripts/gpt_keynote/catalog_keynote.py
```

The script:

1. writes the curated catalog and gene checks;
2. queries NCBI GEO for the literal term `KEYNOTE`;
3. records the live API response and interpretation;
4. refuses automatic download unless a curated row is simultaneously KEYNOTE, lung, processed, and open;
5. writes SHA-256 checksums.

For a network-free deterministic rebuild:

```bash
python3 scripts/gpt_keynote/catalog_keynote.py --offline
```

## Key evidence / 主要证据

- KEYNOTE-001 trial scope: <https://clinicaltrials.gov/study/NCT01295827>
- KEYNOTE-001 NSCLC report: <https://doi.org/10.1056/NEJMoa1501824>
- GEO near-match and explicit melanoma identity: <https://doi.org/10.1038/s41588-022-01268-9>
- KEYNOTE-028 GEP/TMB analysis: <https://doi.org/10.1200/JCO.2018.78.2276>
- KEYNOTE-042 WES/TMB analysis: <https://doi.org/10.1016/j.annonc.2023.01.011>
- KEYNOTE-189/407 WES/TMB and MSD sharing statement: <https://doi.org/10.1016/j.jtocrr.2022.100431>
- KEYNOTE-495 RNA-derived GEP: <https://doi.org/10.1038/s41591-023-02385-6>
- KEYNOTE-782 RNA-seq signatures: <https://doi.org/10.1200/JCO.2024.42.16_suppl.8578>
- MSD data-sharing documentation: <https://engagezone.msd.com/ds_documentation.php>
- Open non-KEYNOTE SU2C-MARK near-match: <https://doi.org/10.5281/zenodo.11179623>

## Limits / 局限

**EN.** “Not found” means not anonymously discoverable through the searched repositories, publication supplements, trial records, and sponsor sharing statements as of the search date; it does not prove that no internal matrix exists. Conference abstracts may precede later deposits. Sponsor-approved access may return derived biomarker covariates rather than raw or normalized gene counts.

**中文。** “未找到”表示截至检索日期，在所查数据库、论文补充材料、试验注册和申办方共享声明中，没有可匿名下载的数据；这并不证明内部不存在矩阵。会议摘要之后仍可能新增存档。即使申请获批，申办方也可能仅提供衍生的生物标志物协变量，而非原始或标准化基因计数。
