# TACSTD2 / CLDN4 vs lung ICI outcomes — GEO leftover slice, 2017–2018

> All outputs live under `results/w200/GEO_2017_2018/` and
> `scripts/w200/geo_2017_2018/`. Every GSE is a real NCBI accession from
> E-utilities. No IDs, outcomes, or p-values were invented.

## English

### 1. Question
Does any **human lung-cancer immune-checkpoint-inhibitor (ICI)** GEO series
**published 2017–2018** jointly deposit **TACSTD2 (TROP2)** and **CLDN4**
expression plus a per-sample ICI endpoint (RECIST / PFS / OS)?

This is the leftover window after the 2019–2021, 2022–2023, and later GEO
slices already in this repo. Those later slices (GSE126044, GSE135222,
GSE207422, …) are **not** re-analyzed here.

### 2. Search
NCBI `gds` esearch, GSE only, *Homo sapiens*, PDAT 2017-01-01 … 2018-12-31,
lung-cancer terms × ICI / PD-1 / PD-L1 / CTLA-4 / named drugs. Combined query
+ per-drug queries: **22 unique series**. UIDs and esummary records:
`search_uids.json`, `candidates_metadata.json`.

### 3. Honest result
**Zero series in this window can test TACSTD2 or CLDN4 against a lung ICI
outcome.** No Mann–Whitney, Cox, or Kaplan–Meier was run. A “non-significant”
p-value would be fabricated.

Two series are genuine leftover **lung ICI** cohorts. Both fail the gene
requirement, for different deposited-data reasons:

| GSE | What it actually is | n lung | ICI + endpoint on GEO | TACSTD2 / CLDN4 |
|---|---|---|---|---|
| **GSE93157** | Prat *Clin Cancer Res* 2017 (PMID 28487385). FFPE tumor, nCounter PanCancer Immune 730. Nivolumab or pembrolizumab. | 35 (22 non-squamous + 13 squamous); plus 25 melanoma + 5 HNSCC | **Yes.** `best.resp` RECIST (lung: 1 CR / 8 PR / 12 SD / 14 PD) and PFS on all 35. | **Absent.** Grep of `GSE93157_raw_data_values.txt.gz` (770 features): no TACSTD*, TROP2, CLDN*, CLAUDIN*. Immune panel does not include these epithelial genes. |
| **GSE110390** | Higgs *Clin Cancer Res* 2018 (PMID 29716923). Durvalumab trial 1108 / NCT01693562. RNA-seq of 97 NSCLC + 62 UC biopsies. | 97 NSCLC samples in the 21-gene file | **No.** Sample records have tissue + “needle biopsy” only. No RECIST/PFS/OS. | **Absent.** Deposited matrix is **21 IFNG-axis genes** (CXCL9, CD274, IFNG, CD8A, …). Full RSEM transcriptome was computed (STAR/RSEM, hg19) but **not released**. SRA exists; this slice does not reprocess SRA. |

Full triage: `tables/triage_all_candidates.csv`.
Gene check: `tables/gene_presence.csv`.
GSE93157 clinical: `tables/GSE93157_clinical.csv`,
`tables/GSE93157_lung_recist_counts.csv`.
GSE110390 gene list: `tables/GSE110390_21genes.csv`.
Empty test table (intentional): `tables/analyzable_marker_outcome.csv`.

### 4. Do not use the GEO `response` field on GSE93157
Every one of the 35 lung samples is labeled `response: RC_RP_SD`, including
the 14 patients whose `best.resp` is PD. The usable clinical field is
`best.resp` (and `pfs` / `pfse`). This does not matter for TACSTD2/CLDN4
because the genes are not on the panel; it would matter if someone reused
this clinical table for the immune-panel genes.

### 5. Why the other 20 hits are not leftover lung ICI tests
- **GSE91061** (Riaz *Cell* 2017): anti-PD-1 ± CTLA-4 **melanoma**. Abstract
  mentions NSCLC; every GSM is `tissue: melanoma`. Out of scope.
- **T-cell / TIL only, no ICI endpoint:** GSE100860 (nivolumab-bound blood
  CD8), GSE99531 (PD-1-high CD8; paper’s ICI cohort is not on the GSM
  records), GSE90728 / GSE90729 (untreated CD8), GSE99254 (scRNA T cells),
  GSE115305 (Siglec-9 CD8). Epithelial TACSTD2/CLDN4 is the wrong
  compartment even if a count matrix were rebuilt from SRA.
- **Not ICI-treated patients:** GSE101929 / GSE102286 (race comparison),
  GSE124199 (dendritic cells / Wnt1), GSE113972 (n=3 lung TSA discovery),
  GSE106420 (LN/BM TRM).
- **Cell line / in vitro / wrong disease:** GSE109010, GSE109020, GSE87879
  (H1299; “PD 0332991” is palbociclib), GSE111360 (organoids), GSE108819
  (breast DC–T), GSE121682 (prostate + olaparib), GSE54781 (AML),
  GSE120545 (transplant endothelium).

### 6. What this slice does **not** claim
- It does **not** say TACSTD2/CLDN4 are unrelated to ICI response. The
  genes were simply **not measured** (GSE93157) or **not released**
  (GSE110390) in the only 2017–2018 lung ICI GEO deposits.
- It does **not** analyze GSE91061 as a lung surrogate.
- It does **not** quantify TACSTD2/CLDN4 in CD8 series as a negative
  control; those matrices are SRA-only and would not answer the tumor-marker
  question.
- Reprocessing GSE110390 SRA to recover TACSTD2/CLDN4 would still lack
  deposited outcomes, so it would not close the question.

### 7. Bottom line
The 2017–2018 GEO leftover window is **closed and empty** for
TACSTD2/CLDN4 × lung ICI. The first usable open GEO tests of these two
genes against lung ICI endpoints appear in later years (already covered
elsewhere in this repo). Anyone citing “no association in 2017–2018 GEO”
should cite **missing measurements**, not a negative statistical result.

---

## 中文

### 1. 问题
**2017–2018 年**发表、涉及**人类肺癌免疫检查点抑制剂（ICI）**的 GEO 系列中，
是否有任何一套数据**同时**提供肿瘤 **TACSTD2 / CLDN4** 表达和逐例 ICI
终点（RECIST / PFS / OS）？

这是本仓库 2019–2021、2022–2023 等切片之后的**剩余窗口**。后期系列
（GSE126044、GSE135222、GSE207422 等）**不在此重复分析**。

### 2. 检索
NCBI `gds`，仅 GSE，人，PDAT 2017-01-01 至 2018-12-31，肺癌词 × ICI / 药名。
合并查询 + 分药查询共 **22 个真实系列**。编号全部来自 NCBI，无虚构。

### 3. 诚实结论
**本窗口内没有任何系列能检验 TACSTD2 或 CLDN4 与肺癌 ICI 结局的关系。**
未做 Mann–Whitney / Cox / 生存分析。给出“无显著”的 p 值属于编造。

仅有两套是真正的肺癌 ICI 队列，都缺靶基因：

- **GSE93157**（Prat 2017）：35 例肺癌（+ 黑色素瘤 / 头颈）抗 PD-1，
  NanoString 730 免疫面板，**有** RECIST 与 PFS；面板上**没有**
  TACSTD2 / CLDN4（已对 770 个探针名全文检索）。
- **GSE110390**（Higgs 2018）：度伐利尤单抗 97 例 NSCLC RNA-seq，
  GEO 只放了 **21 个 IFNG 轴基因**，样本记录**没有**疗效终点；
  TACSTD2 / CLDN4 不在这 21 个基因中。全文转录组算过但未释放。

### 4. GSE93157 的 `response` 字段不可用
35 例肺癌全部标成 `RC_RP_SD`，其中包括 14 例 `best.resp=PD`。
临床应以 `best.resp` 和 `pfs` 为准。对本问题无影响（基因不在面板上），
但若有人拿这张临床表去分析免疫面板基因，必须避开该字段。

### 5. 其余 20 个命中
黑色素瘤 ICB（GSE91061）、纯 T 细胞 / 单细胞、非 ICI 治疗、细胞系
（含把 CDK4/6 抑制剂 PD 0332991 误中为 PD-1）、前列腺癌、AML、移植内皮等。
详见 `tables/triage_all_candidates.csv`。

### 6. 本切片不声称什么
不声称 TACSTD2/CLDN4 与 ICI 无关；只是 **2017–2018 年 GEO 公开层没有测到
或没有放出这两个基因**。不把黑色素瘤当肺癌替代，不重跑 SRA。

### 7. 一句话
2017–2018 GEO 剩余窗口对 TACSTD2/CLDN4 × 肺癌 ICI **已穷尽且为空**。
应表述为**缺少测量**，而不是阴性统计结果。
