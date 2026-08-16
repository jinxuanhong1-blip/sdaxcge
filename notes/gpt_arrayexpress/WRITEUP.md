# ArrayExpress/BioStudies lung ICI / TACSTD2 / CLDN4 search

Search date: 2026-08-16 UTC  
Scope: ArrayExpress studies hosted by EMBL-EBI BioStudies; GEO mirrors excluded.

## 中文

### 结论

在 151 个可复现的同义词交叉检索中，共得到 70 个唯一研究：39 个是 `E-GEOD-*`/GSE 镜像，26 个经完整元数据复核后不符合疾病或干预范围，保留 5 个非 GEO 研究。仓库检索前只有 `README.md`，没有既用 GEO series，因此“已使用 GEO”排除集为空；仍然排除了所有 GEO 镜像。

| accession | 纳入级别 | 内容 | 开放 processed data |
|---|---|---|---|
| [E-MTAB-13704](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13704) | 核心 | 原位肺肿瘤 GEMM；vehicle、aPD-L1 及 4 个联合治疗组，27 个样本 | 1 个矩阵，13,938,914 B |
| [E-MTAB-15883](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15883) | 核心 | 皮下 Lewis lung carcinoma；NC、CTX、ICI、ICI+CTX；PD-1 fate-mapped 单细胞 | 6 个文件，505,837,273 B |
| [E-MTAB-9451](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9451) | 背景 | 5 位 NSCLC 患者配对肺肿瘤/外周血；研究动机涉及 checkpoint response，但患者未接受 ICI | 10 个文件，64,790,171 B |
| [E-MTAB-10633](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10633) | 仅元数据 | 全肺放疗 ± Bintrafusp（TGF-β trap/anti-PD-L1），7 个 pooled 样本 | 无 |
| [E-MTAB-8867](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8867) | 仅元数据 | ICI 心肌炎心肌活检；9 人中 3 人有 NSCLC 病史 | 无 |

共下载 27 个文件、584,712,106 B，包括 17 个 processed 文件及用于标签审计的 10 个 IDF/SDRF 文件。每个 processed 文件均公开且严格小于 2,000,000,000 B；SHA-256、原始 URL 和精确字节数见 `results/gpt_arrayexpress/download_manifest.json`。大型下载保留在本地 `results/gpt_arrayexpress/downloads/`，可由脚本和 manifest 重建。

### 标签与分析

- **E-MTAB-13704：标签完整。** SDRF 可将矩阵 `R1`–`R27` 一一映射到 6 个 stimulus 组（每组 n=4–5），还提供动物、肿瘤数和采样信息。TACSTD2 与 CLDN4 均存在于 processed matrix。相对 vehicle 的描述性均值 log2 ratio：
  - TACSTD2：aPD-L1 −0.107；ATRi/aPD-L1 +0.108；VEGFRi/aPD-L1 +0.378。
  - CLDN4：aPD-L1 −0.097；ATRi/aPD-L1 −0.242；VEGFRi/aPD-L1 −0.499。
  - 文件名虽为 `raw_counts`，值却是非整数，因此没有把它当原始计数做显著性检验。
- **E-MTAB-15883：标签完整。** 21,041 个 annotated cells 全部匹配 matrix barcode；治疗标签为 NC 5,847、CTX 5,299、ICI 4,113、ICI+CTX 5,782。两个靶点都极少表达：TACSTD2 各组非零比例 0.073%–0.340%，CLDN4 为 0.170%–0.363%。这是以免疫细胞为重点的数据，低表达不能外推为肿瘤细胞不表达。
- **E-MTAB-9451：有患者、组织、病理标签，但无 ICI treatment/response 标签。** 10 个库均使用相同的 397-gene targeted immune panel；TACSTD2 和 CLDN4 都不在 panel 中，故不能分析这两个基因。
- **E-MTAB-10633：** 有 compound 与 irradiation 标签，但无 processed data。
- **E-MTAB-8867：** 有原发癌种、ICI 药物和 High/Low CD8 标签；只有 3 个 NSCLC 病例，且无 processed data，也没有 response 标签。

以上只作描述性分析，不把 treatment group 当作临床 response。完整数值见 `target_expression_by_group.csv`，标签枚举见 `label_audit.csv`。

### 检索与复现

`search_biostudies.py` 对 ArrayExpress/BioStudies collection API 完整分页，组合以下概念：

- lung：`lung`、`"lung cancer"`、`"non-small cell lung"`、`NSCLC`、`LLC`
- ICI：PD-1/PD-L1 各种连字符形式、`ICI`、checkpoint/blockade、immunotherapy，以及 nivolumab、pembrolizumab、atezolizumab、durvalumab、cemiplimab、ipilimumab
- targets：`TACSTD2`、`TROP2`、`TROP-2`、`CLDN4`、`claudin 4`、`claudin-4`

另外单独运行 target-only 查询，以防肺相关文本分词异常。每个候选再通过 `/studies/{accession}` 获取完整元数据；`E-GEOD-*` 或元数据内含 GSE/E-GEOD 标识的研究统一作为 GEO 镜像排除。检索是按已归档元数据穷尽同义词，不代表能发现元数据从未提及这些概念的研究。

```bash
python3 scripts/gpt_arrayexpress/search_biostudies.py
python3 scripts/gpt_arrayexpress/download_biostudies.py
python3 scripts/gpt_arrayexpress/analyze_downloads.py
```

## English

### Bottom line

The 151 reproducible synonym-intersection queries returned 70 unique studies. Full-metadata review classified 39 as `E-GEOD-*`/GSE mirrors and 26 as disease/intervention mismatches, leaving five non-GEO records. The repository contained only `README.md` before this slice, so there was no pre-existing used-GEO list; all GEO mirrors were excluded regardless.

Two studies have directly analyzable processed data: **E-MTAB-13704** (27 in-situ lung GEMM tumors across vehicle, aPD-L1, and four combination groups) and **E-MTAB-15883** (PD-1-fate-mapped single cells from subcutaneous LLC under NC, CTX, ICI, and ICI+CTX). **E-MTAB-9451** was retained as NSCLC checkpoint-context data, but its patients were not ICI-treated. **E-MTAB-10633** and **E-MTAB-8867** are relevant metadata-only records with no deposited processed data.

In total, 27 files (584,712,106 bytes) were downloaded: 17 processed files plus 10 IDF/SDRF metadata files needed to audit labels. Every processed file is open and individually below 2,000,000,000 bytes. Exact sizes, URLs, and SHA-256 checksums are in `download_manifest.json`.

### Labels and analysis

- **E-MTAB-13704 has usable treatment labels.** SDRF maps all 27 matrix columns to six groups (n=4–5). Both targets are present. Descriptive mean log2 ratios versus vehicle were −0.107 for TACSTD2 and −0.097 for CLDN4 under aPD-L1 alone. No count-based significance test was run because the file called `raw_counts` contains non-integer values.
- **E-MTAB-15883 has cell-level treatment labels.** All 21,041 annotated cells matched matrix barcodes. TACSTD2 was nonzero in 0.073%–0.340% and CLDN4 in 0.170%–0.363% of cells across groups. This immune-cell-focused resource cannot establish tumor-cell target absence.
- **E-MTAB-9451 has patient, tissue, and histology labels, but no ICI-treatment or response label.** Its common 397-gene targeted panel contains neither TACSTD2 nor CLDN4.
- **E-MTAB-10633** has compound and irradiation labels but no processed data.
- **E-MTAB-8867** has primary-cancer, ICI-drug, and High/Low-CD8 labels; only three subjects had NSCLC, processed data are absent, and response is not labeled.

These are descriptive treatment-group summaries, not clinical-response claims. See `results/gpt_arrayexpress/target_expression_by_group.csv`, `label_audit.csv`, and `analysis_summary.json` for machine-readable results.

### Reproducibility and limits

The scripts fully paginate the authoritative ArrayExpress collection endpoint in BioStudies, fetch complete study records, reject GEO identifiers, enforce the per-file size cap, verify byte sizes, hash downloads, and analyze labels/targets without third-party Python packages. The search is exhaustive over the documented synonym matrix and archived metadata; no text search can recover a study whose metadata never mentions any searched lung, ICI, or target concept.
