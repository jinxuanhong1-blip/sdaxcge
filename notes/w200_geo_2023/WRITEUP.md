# GEO 2023 leftover lung ICI series: TACSTD2 / CLDN4

Slice `w200` / `results/w200/GEO_2023/`. Independent re-search of GEO for
**2023** human lung immune-checkpoint-inhibitor (ICI) series, minus the one
2023 lung cohort already tested in `fable_geo_2022_2023`.

## 1. What "leftover" means

The earlier slice `fable_geo_2022_2023` searched 2022–2023 together and
**actually tested TACSTD2/CLDN4 vs ICI outcome in only one 2023 lung series**:
**GSE207422** (NSCLC neoadjuvant anti-PD-1 + chemo, n=24 baseline tumours;
directionally lower markers in MPR, all p > 0.2). GSE243238 (2023 acral
melanoma) was a non-lung cross-check.

This slice **does not re-test those two**. Leftover = every other 2023 GSE
returned by an independent human × lung × ICI search. Most leftovers were
listed in the earlier catalog but never probed for gene presence.

Already excluded as analysed:
- GSE207422: already analysed in fable_geo_2022_2023 (NSCLC neoadjuvant anti-PD-1+chemo, TACSTD2/CLDN4 vs MPR)
- GSE243238: already analysed as non-lung melanoma cross-check in fable_geo_2022_2023

## 2. Search

- Database: NCBI GEO via E-utilities `db=gds`, `gse[ETYP]`, *Homo sapiens*,
  PDAT 2023-01-01 … 2023-12-31.
- Union of a broad lung×ICI query plus per-drug and per-histology queries.
- **56** unique GSE series dated 2023.
- **2** already analysed (excluded).
- **21** leftover series pass the text filter
  (human + lung + ICI + expression).

Files: `search_candidates_2023.csv`, `leftover_classification.csv`,
`leftover_probe.csv`, `leftover_suppl_files.csv`, `leftover_catalog.csv`.

## 3. Leftover relevant series

| Accession | n | Category | TACSTD2 | CLDN4 | Outcome keys | Verdict |
|---|---:|---|---|---|---|---|
| GSE150255 | 8 | bulk_or_other | False | True | — | leftover: processed matrix peeked; TACSTD2/CLDN4 not both present or outcome missing |
| GSE185204 | 40 | single_cell;blood_or_pbmc | False | False | — | leftover: single-cell / sorted-immune; not bulk tumour epithelium for TACSTD2/CLDN4 vs ICI |
| GSE185205 | 36 | single_cell;blood_or_pbmc |  |  | — | leftover: single-cell / sorted-immune; not bulk tumour epithelium for TACSTD2/CLDN4 vs ICI |
| GSE185206 | 108 | single_cell;blood_or_pbmc |  |  | — | leftover: single-cell / sorted-immune; not bulk tumour epithelium for TACSTD2/CLDN4 vs ICI |
| GSE186446 | 32 | single_cell;blood_or_pbmc | False | False | — | leftover: single-cell / sorted-immune; not bulk tumour epithelium for TACSTD2/CLDN4 vs ICI |
| GSE189804 | 9 | single_cell;cell_line | False | False | — | leftover: single-cell / sorted-immune; not bulk tumour epithelium for TACSTD2/CLDN4 vs ICI |
| GSE213902 | 6 | single_cell;blood_or_pbmc | False | False | — | leftover: single-cell / sorted-immune; not bulk tumour epithelium for TACSTD2/CLDN4 vs ICI |
| GSE216297 | 286 | platelet | False | False | — | leftover: platelet RNA, not tumour |
| GSE217451 | 9 | bulk_or_other | True | True | — | leftover: H1650 si-hMENA cell-line RNA-seq; both genes present; not a patient ICI cohort |
| GSE218402 | 8 | cell_line | True | False | — | leftover: cell-line / in-vitro, not a patient ICI outcome cohort |
| GSE223779 | 6 | single_cell |  |  | — | leftover: single-cell / sorted-immune; not bulk tumour epithelium for TACSTD2/CLDN4 vs ICI |
| GSE224099 | 12 | single_cell;cell_line;blood_or_pbmc | True | True | — | leftover: melanoma TIL CD4 T-cell subsets (not lung epithelium); not TACSTD2/CLDN4 tumour ICI |
| GSE224216 | 9 | chromatin;targeted_panel | True | True | — | leftover: H2030 si-hMENA cell-line RNA-seq; both genes present; not a patient ICI cohort |
| GSE224246 | 8 | chromatin;targeted_panel |  |  | — | leftover: chromatin assay (ATAC/ChIP), not gene expression of TACSTD2/CLDN4 |
| GSE229353 | 7 | single_cell |  |  | pathology | leftover: single-cell / sorted-immune; not bulk tumour epithelium for TACSTD2/CLDN4 vs ICI |
| GSE235048 | 15 | blood_or_pbmc;flow_or_cytof | True | False | — | leftover: blood/PBMC, not tumour epithelium |
| GSE235500 | 34 | single_cell;has_mouse_or_xeno |  |  | — | leftover: single-cell / sorted-immune; not bulk tumour epithelium for TACSTD2/CLDN4 vs ICI |
| GSE235603 | 86 | bulk_or_other |  |  | — | leftover: SuperSeries of tumour-infiltrating Treg scRNA; no bulk tumour matrix |
| GSE238006 | 8 | has_mouse_or_xeno | False | False | — | leftover: processed matrix peeked; TACSTD2/CLDN4 not both present or outcome missing |
| GSE248378 | 29 | bulk_or_other | True | True | — | leftover ANALYSED: post-durvalumab non-MPR bulk FPKM, both genes present; GEO has arm/histology only; recurrence recovered from Nat Commun source data Fig.5d via unique ITGAE FPKM (9 recur / 20 no recur) |
| GSE248830 | 44 | targeted_panel | False | False | — | leftover: targeted panel — gene presence must be checked; often lacks TACSTD2/CLDN4 |

## 4. Association tests on leftovers

Leftover series flagged `usable_maybe` (both genes seen in an open processed
matrix **and** a per-sample outcome field): **3**.

Tests run: **5** (4 outcome/arm tests + 1 TACSTD2–CLDN4 Spearman). Outcome
tests with p < 0.05: **0**. The Spearman (ρ = 0.77, p < 10⁻⁴) is
co-expression, not an ICI-outcome test. GSE248378 is Twist-exome capture
RNA-seq (~15k genes), post-treatment, non-MPR only.

| Cohort | Gene | Outcome | n1 | median1 | n2 | median2 | p | Cliff δ |
|---|---|---|---:|---:|---:|---:|---:|---:|
| GSE248378 | TACSTD2 | recurrence_in_nonMPR | 20 | 41.98 | 9 | 80.01 | 0.0562 | -0.456 |
| GSE248378 | TACSTD2 | GEO_treatment_arm | 18 | 44.925 | 11 | 122.63 | 0.2517 | -0.263 |
| GSE248378 | CLDN4 | recurrence_in_nonMPR | 20 | 38.97 | 9 | 77.37 | 0.3108 | -0.244 |
| GSE248378 | CLDN4 | GEO_treatment_arm | 18 | 41.04 | 11 | 87.21 | 0.1568 | -0.323 |
| GSE248378 | TACSTD2_vs_CLDN4 | spearman_within_cohort | 29 | 45.92 | 29 | 52.6 | 0.0 |  |

Skipped after a usable_maybe flag:
- GSE217451: sample overlap <4 (n=0)
- GSE224099: expr parse failed: list index out of range
- GSE224216: sample overlap <4 (n=0)

## 5. Honest conclusion

The only leftover 2023 lung series that can be tested is **GSE248378**
(Altorki et al., NCT02904954): post-durvalumab ± SBRT *non-MPR* resected
tumours, n=29, both genes present. GEO itself has no recurrence label.
Recurrence (9 vs 20) was recovered from the paper's public Source Data
Fig. 5d by matching unique ITGAE FPKM values — not by guessing that a
sample title ending in "R" means recurrence (that guess is false: POD23R /
POD26R / POD27R / POD32R did *not* recur).

TACSTD2 is higher in the 9 tumours that later recurred (median FPKM 80.0 vs
42.0; Mann–Whitney p=0.056; Cliff δ −0.46). CLDN4 is in the same direction
and weaker (77.4 vs 39.0; p=0.31). Neither test is p<0.05. The cohort is
**post-treatment and non-MPR only**, so it cannot address baseline
prediction or MPR. It is not a second GSE207422-class baseline
whole-transcriptome ICI-response cohort.

Every other leftover is single-cell T cells, platelets, PBMCs, cell-line
siRNA, ATAC, or a targeted panel without these two genes. Do not treat the
leftover catalog as confirmed biomarker evidence. The TACSTD2 p=0.056
result is hypothesis-generating and selected.

## 6. Limitations

- Text filters can miss a series whose GEO title/summary never says lung or
  ICI; per-drug queries reduce but do not eliminate that risk.
- Gene peeking reads the first column of the smallest processed file < 80 MB.
  A gene present only in a larger or binary (`.rds` / `.RData`) file can be
  missed; those files are recorded in `leftover_suppl_files.csv`.
- No raw-read realignment. Author-processed matrices only.
- Single-cell leftovers were not pseudo-bulked. TACSTD2/CLDN4 are epithelial;
  T-cell / Treg / PBMC series cannot answer the tumour-intrinsic question.

---

# GEO 2023 年肺癌 ICI「剩余」系列：TACSTD2 / CLDN4

切片 `w200` / `results/w200/GEO_2023/`。对 GEO **2023** 年人类肺癌免疫检查点
抑制剂（ICI）系列做独立再检索，并减去 `fable_geo_2022_2023` 中已经检验过的
2023 年肺癌队列。

## 1. 「剩余」的定义

先前切片 `fable_geo_2022_2023` 把 2022–2023 放在一起检索，但 **真正对
TACSTD2/CLDN4 与 ICI 结局做了统计检验的 2023 年肺癌系列只有 GSE207422**
（NSCLC 新辅助抗 PD-1 + 化疗，基线肿瘤 n=24；MPR 组标志物更低，但全部
p > 0.2）。GSE243238（2023 肢端黑色素瘤）是非肺交叉验证。

本切片 **不再重复检验这两套**。剩余 = 独立检索得到的其余全部 2023 年 GSE。
其中多数在先前目录里出现过，但从未核查基因是否存在。

已排除：
- GSE207422: already analysed in fable_geo_2022_2023 (NSCLC neoadjuvant anti-PD-1+chemo, TACSTD2/CLDN4 vs MPR)
- GSE243238: already analysed as non-lung melanoma cross-check in fable_geo_2022_2023

## 2. 检索

- NCBI GEO E-utilities，`gse[ETYP]`，人，PDAT 2023-01-01 至 2023-12-31。
- 宽查询 ∪ 逐药 ∪ 逐组织学。
- **56** 个 2023 年 GSE。
- **2** 个已分析（排除）。
- **21** 个剩余系列通过文本过滤（人 + 肺 + ICI + 表达）。

## 3. 剩余相关系列

见英文表。完整字段在 `leftover_catalog.csv`。

## 4. 剩余系列上的关联检验

`usable_maybe`（开放矩阵中见到两个基因 **且** 有逐样本结局字段）：**3**。
已跑检验：**5**（4 个结局/臂次 + 1 个 TACSTD2–CLDN4 Spearman）。结局检验
p < 0.05：**0**。Spearman（ρ = 0.77）是共表达，不是 ICI 结局。GSE248378
是 Twist 外显子捕获 RNA-seq（约 1.5 万基因），且仅为治疗后非 MPR。

| Cohort | Gene | Outcome | n1 | median1 | n2 | median2 | p | Cliff δ |
|---|---|---|---:|---:|---:|---:|---:|---:|
| GSE248378 | TACSTD2 | recurrence_in_nonMPR | 20 | 41.98 | 9 | 80.01 | 0.0562 | -0.456 |
| GSE248378 | TACSTD2 | GEO_treatment_arm | 18 | 44.925 | 11 | 122.63 | 0.2517 | -0.263 |
| GSE248378 | CLDN4 | recurrence_in_nonMPR | 20 | 38.97 | 9 | 77.37 | 0.3108 | -0.244 |
| GSE248378 | CLDN4 | GEO_treatment_arm | 18 | 41.04 | 11 | 87.21 | 0.1568 | -0.323 |
| GSE248378 | TACSTD2_vs_CLDN4 | spearman_within_cohort | 29 | 45.92 | 29 | 52.6 | 0.0 |  |

## 5. 诚实结论

唯一能做检验的 2023 年剩余肺系列是 **GSE248378**（Altorki 等，NCT02904954）：
度伐利尤单抗 ± SBRT **治疗后、非 MPR** 切除肿瘤，n=29，两个基因都在。
GEO 没有复发标签。复发（9 vs 20）是用论文公开 Source Data 图 5d 的
ITGAE FPKM 一一对上的，**不是**把样本名末尾的 “R” 当成复发
（POD23R/26R/27R/32R 实际未复发）。

TACSTD2 在随后复发的 9 例更高（中位 FPKM 80.0 vs 42.0；Mann–Whitney
p=0.056；Cliff δ −0.46）。CLDN4 同向更弱（77.4 vs 39.0；p=0.31）。
都不是 p<0.05。队列是 **治疗后且仅非 MPR**，不能回答基线预测或 MPR。
它也不是第二套 GSE207422 级别的基线全转录组 ICI 反应队列。

其余剩余系列是单细胞 T 细胞、血小板、PBMC、细胞系 siRNA、ATAC，
或没有这两个基因的靶向 panel。不要把剩余目录当成已证实的标志物证据。
TACSTD2 p=0.056 只是假设生成、且经过选择。

## 6. 限制

- 标题/摘要从未写 lung 或 ICI 的系列可能漏检。
- 基因探测只读 < 80 MB 最小处理文件的第一列；`.rds`/`.RData` 可能漏检。
- 不做原始读段重比对。
- 不对单细胞剩余系列做伪 bulk：TACSTD2/CLDN4 是上皮基因。

