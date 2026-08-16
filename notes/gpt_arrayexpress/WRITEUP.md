# ArrayExpress/BioStudies lung ICI / TACSTD2 / CLDN4 search

Search date: 2026-08-16 UTC  
Scope: ArrayExpress studies hosted by EMBL-EBI BioStudies; GEO mirrors excluded.

## 中文

### 结论

在 204 个可复现同义词交叉检索中，共得到 80 个唯一研究：47 个是 `E-GEOD-*`/GSE 镜像，21 个经完整元数据复核后不符合疾病或干预范围。仓库检索前只有 `README.md`，没有既用 GEO series，因此“已使用 GEO”排除集为空；仍然排除了所有 GEO 镜像。

本轮补查的 leftover 肺 ICI 相关、且 processed matrix 可用的非 GEO 研究只有 **E-MTAB-15235**。它的标题和摘要提到 KRAS 突变 NSCLC 以及 HUDSON 试验中 ceralasertib + durvalumab，但**入库矩阵的实验因素是 ATRi（ceralasertib）对 vehicle，不是 ICI 治疗或临床 response**。矩阵 80 列全部能映射到 SDRF；TACSTD2（`ENSG00000184292`）和 CLDN4（`ENSG00000189143`）都在。

| accession | 纳入级别 | 内容 | 开放 processed data |
|---|---|---|---|
| [E-MTAB-13704](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13704) | 核心 | 原位肺肿瘤 GEMM；vehicle、aPD-L1 及 4 个联合治疗组，27 个样本 | 1 个矩阵，13,938,914 B |
| [E-MTAB-15883](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15883) | 核心 | 皮下 Lewis lung carcinoma；NC、CTX、ICI、ICI+CTX；PD-1 fate-mapped 单细胞 | 6 个文件，505,837,273 B |
| [E-MTAB-15235](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15235) | leftover | NSCLC 相关记录；80 个样本的 ceralasertib vs vehicle 矩阵；摘要提及 durvalumab | 1 个矩阵，29,797,813 B |
| [E-MTAB-9451](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9451) | 背景 | 5 位 NSCLC 患者配对肺肿瘤/外周血；无 ICI 治疗标签 | 10 个文件，64,790,171 B |
| [E-MTAB-10633](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10633) | 仅元数据 | 全肺放疗 ± Bintrafusp（TGF-β trap/anti-PD-L1） | 无 |
| [E-MTAB-8867](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8867) | 仅元数据 | ICI 心肌炎；9 人中 3 人有 NSCLC 病史 | 无 |

另有非肺、直接 TROP2/TACSTD2 提交：E-MTAB-11466、E-MTAB-16433、E-MTAB-16843、E-MTAB-16849（有 processed）；E-MTAB-11377、E-MTAB-11382（仅元数据）。未发现直接非 GEO 的 CLDN4 中心提交。未下载超 2 GB 的 CRC TROP2 矩阵。

共下载 34 个文件、614,639,597 B（含 leftover 矩阵和标签元数据）。每个 processed 文件均公开且严格小于 2,000,000,000 B。SHA-256 见 `results/gpt_arrayexpress/download_manifest.json`。

### 标签与分析

- **E-MTAB-13704：标签完整。** SDRF 将 `R1`–`R27` 映射到 6 个 stimulus 组（n=4–5）。两靶点都在矩阵中。相对 vehicle 的描述性均值 log2 ratio：TACSTD2 在 aPD-L1 为 −0.107，CLDN4 为 −0.097。文件名虽为 `raw_counts`，值却是非整数，未做计数检验。
- **E-MTAB-15883：标签完整。** 21,041 个 annotated cells 全部匹配 barcode。TACSTD2 各组非零比例 0.073%–0.340%，CLDN4 为 0.170%–0.363%。这是免疫细胞数据，低表达不能外推为肿瘤细胞不表达。
- **E-MTAB-15235：有 ATRi 标签，没有 ICI/response 标签。** 80/80 列可映射到 SDRF Source Name。CTG1955 中 TACSTD2 组均值约 10.01–10.26，CLDN4 约 8.66–9.18；CTG2026 中 TACSTD2 为负值（约 −2.93 至 −2.19），CLDN4 约 2.07–2.70。这些是入库矩阵的描述性均值，不是检验统计量。SDRF 把全部样本标成 CT26/`Mus musculus`，矩阵却用人类 Ensembl ID；分组按列名中的个体、化合物和天数。
- **E-MTAB-9451：** 有患者/组织/病理标签，无 ICI 标签；397-gene immune panel 不含 TACSTD2/CLDN4。
- **E-MTAB-10633 / E-MTAB-8867：** 有相关标签，无 processed data。

补查排除的邻近记录：E-MTAB-12508（肺肿瘤 DC-therapy，非 ICI，无 processed）；E-MTAB-13713 / E-MTAB-13703（ceralasertib，无肺/ICI 标签）；E-MTAB-13823（CT26 TREX1，非肺 ICI）。

### 检索与复现

```bash
python3 scripts/gpt_arrayexpress/search_biostudies.py
python3 scripts/gpt_arrayexpress/download_biostudies.py
python3 scripts/gpt_arrayexpress/analyze_downloads.py
```

## English

### Bottom line

204 paginated synonym-intersection queries returned 80 unique ArrayExpress records. Full-metadata review classified 47 as GEO mirrors and 21 as disease/intervention mismatches. The repository had no pre-existing used-GEO list.

The only leftover non-GEO lung ICI-adjacent series with a usable processed matrix is **E-MTAB-15235**. The record discusses KRAS-mutant NSCLC and HUDSON ceralasertib plus durvalumab, but the deposited 80-sample matrix is ATRi versus vehicle. Both TACSTD2 and CLDN4 are present as human Ensembl IDs. No leftover series provided ICI-treated samples plus both genes beyond the already analyzed E-MTAB-13704 and E-MTAB-15883.

### Labels and analysis

- **E-MTAB-13704** has usable ICI combination labels. Descriptive mean log2 ratios versus vehicle were −0.107 for TACSTD2 and −0.097 for CLDN4 under aPD-L1 alone. No count-based test was run because the so-called raw-count file is non-integer.
- **E-MTAB-15883** has cell-level ICI labels. TACSTD2 was nonzero in 0.073%–0.340% and CLDN4 in 0.170%–0.363% of 21,041 annotated cells.
- **E-MTAB-15235** has ceralasertib/vehicle and day-7/day-10 labels, not ICI or response labels. All 80 matrix columns map to SDRF source names. Group means were computed from the deposited values only: TACSTD2 stayed near 10.0–10.3 in CTG1955 and negative in CTG2026; CLDN4 stayed near 8.7–9.2 in CTG1955 and 2.1–2.7 in CTG2026. SDRF marks every sample as CT26/`Mus musculus` while the matrix uses human gene IDs.
- **E-MTAB-9451** cannot be used for TACSTD2/CLDN4: both genes are absent from its 397-gene targeted panel.
- Direct non-lung TROP2 submissions exist (E-MTAB-11466/16433/16843/16849). No direct non-GEO CLDN4-centered submission was found.

These are descriptive summaries, not clinical-response claims. See `target_expression_by_group.csv`, `label_audit.csv`, and `analysis_summary.json`.

### Reproducibility and limits

The scripts paginate the ArrayExpress collection API, fetch complete study JSON, reject GEO identifiers, enforce the per-file <2 GB cap, verify byte sizes, and hash downloads. The leftover pass added HUDSON/ceralasertib/bintrafusp/CTLA-4 terms. Text search cannot recover a study whose archived metadata never mentions a searched concept.
