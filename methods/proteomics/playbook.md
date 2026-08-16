# CPTAC/TMT + PRIDE proteomics playbook

**Scope / 范围.** TROP2 (`TACSTD2`, Ensembl `ENSG00000184292`, UniProt P09758) and CLDN4 (`ENSG00000189143`, UniProt O14493) in **CPTAC LUAD/LSCC** TMT gene-abundance tables, public **GDC S3 STAR** gene-counts, and **PRIDE** protein tables. Outputs live only under `methods/proteomics/`.

**Demo / 演示图.** `methods/proteomics/demo/figures/` (PNG+PDF). Rebuild: `python3 methods/proteomics/scripts/run_demo.py`.

---

## 0. IDs and files / 标识与文件

### English

| Entity | IDs you must accept |
| --- | --- |
| TROP2 / TACSTD2 | `ENSG00000184292` (± version `.7`), `TACSTD2`, `TROP2`, P09758 |
| CLDN4 | `ENSG00000189143` (± version `.9`), `CLDN4`, `CLD4`, O14493 |

Two public file families use **Ensembl rows**:

1. **GDC CPTAC-3 open STAR gene-counts** (this playbook’s working S3 source).  
   Bucket: `s3://gdc-cptac-phs001287-2-open/{file_id}/{uuid}.rna_seq.augmented_star_gene_counts.tsv`  
   (`aws s3 cp --no-sign-request …` or `https://api.gdc.cancer.gov/data/{file_id}`).  
   Filter: project `CPTAC-3`, `primary_site = Bronchus and lung`, STAR – Counts, access = open.  
   Columns: `gene_id` (`ENSG00000184292.7`, `ENSG00000189143.9`), `gene_name`, `unstranded`, `tpm_unstranded`, `fpkm_unstranded`.

2. **CPTAC pan-cancer gene-abundance** (historical freeze, Ensembl rows), e.g.  
   `LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt`  
   and the LSCC sibling. Values are **already log2** reference-intensity-normalized TMT. The public HTTP bucket is often **access-gated** (CDS/dbGaP). If you have a local mirror, drop it in and the same parsers apply.

**LinkedOmics CCT** (used for the demo because it is openly downloadable) uses **gene symbols**, TMT **log2-ratio**, and **NArm** (high-missingness genes dropped). That NArm step is why CLDN4 can vanish as a *row*, not as a row of NAs.

**PRIDE** protein tables are heterogeneous: MaxQuant `proteinGroups.txt`, DIA-NN `report.pg_matrix.tsv`, mzTab. Search P09758 / O14493 / `TACSTD2` / `CLDN4`. Absence after a documented search is a result (e.g. ICI plasma PXD042091 often lacks both targets).

### 中文

| 分子 | 必须同时识别的 ID |
| --- | --- |
| TROP2 / TACSTD2 | `ENSG00000184292`（可带版本号 `.7`）、`TACSTD2`、`TROP2`、P09758 |
| CLDN4 | `ENSG00000189143`（可带 `.9`）、`CLDN4`、`CLD4`、O14493 |

两类公开文件以 **Ensembl 为行名**：

1. **GDC CPTAC-3 开放 STAR 基因计数**（本手册可复现的 S3 来源）。  
   `s3://gdc-cptac-phs001287-2-open/{file_id}/…star_gene_counts.tsv`。  
   过滤：`CPTAC-3`、肺、STAR – Counts、open。行名即 `ENSG00000184292.7` / `ENSG00000189143.9`。

2. **CPTAC 泛癌 gene-abundance**（历史 v1.2 冻结，Ensembl 行，**已是 log2**）。公网桶常被 CDS/dbGaP 闸住；有本地镜像即可用同一解析器。

**LinkedOmics CCT**（演示用，可直接下）用 **基因符号**、TMT **log2 比值**，并做了 **NArm**。CLDN4 常常不是“一行全是 NA”，而是 **整行被丢掉**。

**PRIDE** 表格式不统一。按 P09758 / O14493 / 符号检索；检索后确认不存在也是结果，不要编造定量值。

---

## 1. log2 / 对数变换

### English

**Rule.** Detect scale before transforming. Double-logging TMT ratios is the most common silent error.

| Table | Typical values | Action |
| --- | --- | --- |
| CPTAC TMT gene-abundance / LinkedOmics proteome CCT | signed, \|x\| ≲ 15, median near 0 | **already log2-ratio to pool**. Do not log again. |
| GDC STAR `unstranded` | non-negative integers, 10³–10⁶ | raw counts → `log2(CPM + prior)` for limma-trend; or use TPM/FPKM then `log2(x+1)` |
| GDC STAR `tpm` / `fpkm` | ~0–200 for these genes | not log2; `log2(x+1)` if you need log space |
| LinkedOmics LUAD RNA CCT | ~−3 to 8 | already `log2(UQ-FPKM)` |
| LinkedOmics LSCC RNA CCT | often ~5–24 | already log2, **different offset/scale than LUAD**. Do not concatenate LUAD+LSCC RNA without re-normalization. |
| PRIDE LFQ / iBAQ / intensity | 10⁴–10¹⁰, no negatives | `log2(x)` or `log2(x+1)` after 0-handling |
| PRIDE reporter / TMT ratio exports | may already be log2 | use `detect_log_scale()` |

`scripts/lib_io.py:detect_log_scale` / `apply_log2_if_needed` implement the heuristic: negatives + compact range ⇒ already log; wide non-negative ⇒ log2(x+1); ambiguous ⇒ **do not** log (fail safe).

**Do not** mix LUAD and LSCC RNA numeric values in one model. Protein TMT log-ratios are also **within-study** (common pool per plex/study), not pan-cancer comparable without the official freeze.

### 中文

**原则：** 先判断尺度再变换。对已经是 TMT log2 比值的表再取 log，是最常见的静默错误。

| 表 | 典型数值 | 处理 |
| --- | --- | --- |
| CPTAC TMT gene-abundance / LinkedOmics 蛋白 CCT | 有正有负，\|x\| ≲ 15 | **已是相对 pool 的 log2 比**，禁止再 log |
| GDC STAR `unstranded` | 非负整数 | 原始计数 → limma-trend 用 `log2(CPM+prior)` |
| GDC STAR TPM/FPKM | 约 0–200 | 尚未 log；需要对数空间时 `log2(x+1)` |
| LinkedOmics LUAD RNA | 约 −3–8 | 已是 log2(UQ-FPKM) |
| LinkedOmics LSCC RNA | 常为 5–24 | 已是 log2，**与 LUAD 不在同一尺度**，不能直接拼 |
| PRIDE LFQ/iBAQ | 大范围正数 | `log2(x)` / `log2(x+1)` |
| PRIDE TMT 比值导出 | 可能已 log2 | 走检测函数 |

`detect_log_scale`：出现负数且范围紧凑 → 已 log；大范围非负 → log2(x+1)；不确定 → **不变换**。

LUAD 与 LSCC 的 RNA 数值不要放进同一个线性模型。TMT 比值只在同一研究/同一 pool 内可比。

---

## 2. Missingness and CLDN4 NAs / 缺失与 CLDN4

### English

CLDN4 is a small (~22 kDa) tetraspan transmembrane protein with few proteotypic tryptic peptides. In DDA/TMT roll-up it is often **not quantified**. After **NArm**, the gene is **deleted from the matrix**.

**What we observe on the public LinkedOmics NArm TMT tables (demo):**

| Gene | LUAD protein | LSCC protein | LUAD RNA | LSCC RNA | GDC STAR (Ensembl) |
| --- | --- | --- | --- | --- | --- |
| TACSTD2 / `ENSG00000184292` | present, 0% NA (n=110) | present, 0% NA (n=108) | present | present | present (counts/TPM/FPKM) |
| CLDN4 / `ENSG00000189143` | **no row** | **no row** | present, 0% NA | present, 0% NA | present |

Sister claudins are not interchangeable: LUAD protein has `CLDN3`, `CLDN7`, `CLDN18` but not `CLDN1`/`CLDN4`/`CLDN5`. LSCC protein has `CLDN1` but still **no CLDN4**.

**How to report NAs (do not collapse these):**

1. **Absent row** — gene never quantified or dropped by NArm. Concordance and protein DE are **undefined**.  
2. **Row of NAs / high NA fraction** — MNAR (left-censored, low abundance) is the default for TMT, not MCAR. Do not mean-impute a biomarker.  
3. **Wrong ID** — looking up only `ENSG00000189143` in a symbol table, or only `CLDN4` in an Ensembl gene-abundance file, looks like “all NA”. Always query **both** IDs (`lib_io.lookup_gene_rows`).

**What not to do.** Do not impute CLDN4 protein from RNA. Do not treat “CLDN4 NA” as low protein in an ICI model. If you need CLDN4 protein, use targeted MS, PRM, or IHC — not this TMT gene-abundance table.

Script: `scripts/log2_missingness.py`.

### 中文

CLDN4 是约 22 kDa 的四次跨膜蛋白，可定量的胰酶肽很少。DDA/TMT 汇总时经常 **根本没有定量**；再经过 **NArm**，基因从矩阵里 **整行消失**。

**公开 LinkedOmics NArm TMT 表示例（本仓库演示）：**

| 基因 | LUAD 蛋白 | LSCC 蛋白 | LUAD RNA | LSCC RNA | GDC STAR（Ensembl） |
| --- | --- | --- | --- | --- | --- |
| TACSTD2 / `ENSG00000184292` | 有，0% NA | 有，0% NA | 有 | 有 | 有 |
| CLDN4 / `ENSG00000189143` | **无此行** | **无此行** | 有 | 有 | 有 |

不要用其他 claudin 代替：LUAD 蛋白有 CLDN3/7/18，无 CLDN1/4/5；LSCC 有 CLDN1，**仍无 CLDN4**。

**缺失要分三类写，不要混为一谈：**

1. **行不存在** — 未定量或被 NArm 丢掉。蛋白–RNA 相关和蛋白 DE **无法定义**。  
2. **行在但 NA 很多** — TMT 默认按 MNAR（左删失）而不是 MCAR。不要对生物标志物做均值填补。  
3. **ID 用错** — 在符号表里只查 Ensembl、或在 Ensembl 表里只查 `CLDN4`，会看起来像“全是 NA”。必须两种 ID 都查。

**禁止：** 用 RNA 填补 CLDN4 蛋白；把 “CLDN4 NA” 当成低蛋白送进 ICI 模型。需要 CLDN4 蛋白请用靶向质谱 / PRM / IHC。

---

## 3. Protein–RNA concordance / 蛋白–RNA 一致性

### English

Pair on **harmonized case IDs** (`C3L.00001` = `C3L-00001`). Use **Spearman ρ** on samples with both values. Gene-wise RNA–protein ρ in CPTAC LUAD literature is typically ~0.14–0.53; a single gene can sit anywhere in that range.

**Demo (LinkedOmics tumor tables):** TACSTD2 has paired protein+RNA in LUAD and LSCC (see `demo/tables/protein_rna_concordance.tsv`). CLDN4 has RNA only → status `protein_absent`, ρ = NA. That is the correct output.

If you start from S3 STAR files, collapse per-sample TSVs with `extract_star_targets` (Ensembl → symbol), then join to a protein matrix on `sample_id`. Use the same RNA transform on all samples of one cohort (`log2(TPM+1)` **or** the LinkedOmics log2 table — do not mix).

Script: `scripts/protein_rna_concordance.py`.

### 中文

在 **统一后的病例 ID** 上配对（`C3L.00001` = `C3L-00001`）。只对双端都非缺失的样本算 **Spearman ρ**。文献中 CPTAC LUAD 基因水平 RNA–蛋白相关大约 0.14–0.53，单个基因可以落在任何位置。

**演示：** TACSTD2 在 LUAD/LSCC 都能配对；CLDN4 只有 RNA → 状态为 `protein_absent`，ρ 为 NA。这是正确结果，不是程序失败。

若从 S3 STAR 单样本文件出发：用 `extract_star_targets` 抽出两个 Ensembl，再按 `sample_id` 接到蛋白矩阵。同一队列只用一种 RNA 变换。

---

## 4. limma and limma-trend / 差异分析

### English

| Input | Method | Not this |
| --- | --- | --- |
| TMT log2-ratio (CPTAC, LinkedOmics, CDAP gene-abundance) | `limma` + eBayes. `trend=TRUE` optional (mean–variance on log space). | **voom**, raw-count NB models, extra log2 |
| RNA **counts** (GDC STAR `unstranded`) | `log2(CPM + prior)` then **limma-trend** | voom is acceptable; DESeq2 is fine for counts — pick one |
| RNA already log2(FPKM/UQ) | `limma` on the log table | voom (voom expects counts) |

**Design.** Tumor vs NAT is legitimate on CPTAC. Smoking / histology / TMT plex can be covariates if you have them. **ICI response is not a CPTAC column** — do not code a fake responder factor.

**Missing values.** limma wants a complete matrix. For TMT, drop genes below an observation threshold (demo: 70%), then **gene-median impute** and say so. Do not impute CLDN4: it has no row.

**Python vs R.** `scripts/limma_trend.py` is a Smyth (2004) eBayes port for environments without R. `scripts/limma_trend.R` calls Bioconductor `limma` when installed. For a paper, prefer the R implementation; the Python port is for this demo and CI.

**Demo.** LUAD protein Tumor vs NAT, method `limma` (already logged). TACSTD2 is highlighted on the volcano if present. Output: `demo/tables/limma_LUAD_protein_tumor_vs_nat.tsv`.

### 中文

| 输入 | 方法 | 不要用 |
| --- | --- | --- |
| TMT log2 比值 | `limma` + eBayes，可选 trend | **voom**、负二项、再 log2 |
| RNA **原始计数** | `log2(CPM+prior)` 后 **limma-trend** | 与 voom/DESeq2 三选一，不要混报 |
| 已是 log2(FPKM/UQ) 的 RNA | 直接 `limma` | voom（voom 要计数） |

**设计：** CPTAC 上做肿瘤 vs 癌旁合法。吸烟、组织学、TMT 批次可作协变量。**没有 ICI 缓解字段**，不要编造 responder。

**缺失：** 先按观测比例过滤，再基因中位数填补，并在方法里写明。CLDN4 没有行，不要填。

**Python / R：** 无 R 时用 `limma_trend.py`；写论文优先 Bioconductor `limma`（`limma_trend.R`）。

---

## 5. Joining CIBERSORT / xCell / 对接反卷积

### English

Deconvolution is **RNA-derived**. Join it to protein on **case ID**, not aliquot UUID, unless you verified a 1:1 aliquot map.

| File | Orientation | Join key |
| --- | --- | --- |
| CIBERSORT / CIBERSORTx | samples × cell types; first column `Mixture` | normalize `Mixture` → `C3L-#####` |
| xCell | cell types × samples (or transpose) | same ID normalizer |
| LinkedOmics LSCC `molecular_phenotypes` | `Immune.Cluster.rna` = Cold / Warm / Hot Tumor | already case-level |
| ESTIMATE / xCell published supplements | watch `.` vs `-` | `lib_io.normalize_sample_id` |

`scripts/join_deconvolution.py` also computes a few **RNA marker columns** (`CD8A`, `CD274`, `CXCL13`, …) from the RNA matrix so a join exists when you do not have a CIBERSORT run. Label them as marker scores, not CIBERSORT.

**Interpretation limit.** A correlation between TACSTD2 protein and CD8A RNA / “Cold Tumor” in **treatment-naive** LSCC is a **co-expression / microenvironment** statement. It is **not** an ICI-response association (next section).

### 中文

反卷积来自 **RNA**。与蛋白对接请用 **病例 ID**，不要默认用 aliquot UUID。

| 文件 | 方向 | 键 |
| --- | --- | --- |
| CIBERSORT(x) | 样本 × 细胞；首列 `Mixture` | 规范成 `C3L-#####` |
| xCell | 细胞 × 样本（或转置） | 同一 ID 函数 |
| LinkedOmics LSCC 分子表型 | `Immune.Cluster.rna`（Cold/Warm/Hot） | 已是病例级 |

`join_deconvolution.py` 可从 RNA 抽 `CD8A`/`CD274`/`CXCL13` 等，**须标明是 marker，不是 CIBERSORT**。

TACSTD2 蛋白与 CD8A / “Cold Tumor” 在 **未治疗** LSCC 中的相关，只说明微环境共变，**不能写成 ICI 疗效**。

---

## 6. Treatment-naive ≠ ICI response / 未治疗 ≠ 免疫治疗疗效

### English

CPTAC LUAD (Gillette *et al.*, *Cell* 2020) and LSCC (Satpathy *et al.*, *Cell* 2021) are **mostly surgically resected, treatment-naive** tumors. Clinical tables have stage, smoking, mutations, immune clusters — **not** PD-1/PD-L1 exposure, RECIST, ORR, or IO-PFS.

Therefore:

1. You **may** test TACSTD2/CLDN4 protein vs NAT, vs RNA, vs CIBERSORT/xCell/Immune.Cluster.  
2. You **may not** report those tests as “association with ICI response” or as a replication of Bessede *et al.*, *Clin Cancer Res* 2024 (PMID 38048058; OAK/POPLAR atezolizumab, often EGA-controlled).  
3. Immune-hot treatment-naive tumors are not responders; immune-cold tumors are not primary-resistant to ICI. Those labels require **on-trial outcome**.  
4. Stage, smoking, histology, and TMB associate with both epithelial programs and (in other cohorts) IO outcome — confounding if you smuggle CPTAC into an IO meta-analysis.  
5. PRIDE helps only when the project is ICI-annotated **and** the protein table actually contains P09758/O14493.

**Where ICI labels live (outside this slice):** GEO ICI RNA (e.g. GSE135222), trial supplements, controlled EGA. Catalog them; do not invent p-values from CPTAC.

Figure: `demo/figures/06_naive_not_ici.png`.

### 中文

CPTAC LUAD（Gillette, *Cell* 2020）与 LSCC（Satpathy, *Cell* 2021）以 **手术切除、基线未治疗** 为主。临床表有分期、吸烟、突变、免疫聚类，**没有** PD-1/PD-L1 用药、RECIST、ORR、IO-PFS。

因此：

1. **可以** 做 TACSTD2/CLDN4 蛋白相对癌旁、相对 RNA、相对 CIBERSORT/xCell/Immune.Cluster。  
2. **不可以** 写成 “与 ICI 疗效相关”，也不可当作 Bessede 2024（OAK/POPLAR）的重复。  
3. 未治疗的 immune-hot ≠ 缓解；immune-cold ≠ 原发耐药。疗效标签必须来自 **带结局的 ICI 队列**。  
4. 分期、吸烟、组织学、TMB 会同时关联上皮程序与（其他队列中的）IO 结局，不能把 CPTAC 偷偷并进 IO 荟萃。  
5. PRIDE 只有在项目本身标注 ICI **且** 蛋白表含 P09758/O14493 时才有用。

ICI 标签在 GEO / 试验附录 / EGA。只建目录，不要用 CPTAC 编造 p 值。

---

## 7. Scripts / 脚本

All paths relative to repo root. Python 3.10+ with `pandas`, `numpy`, `scipy`, `matplotlib`, `statsmodels`.

```text
methods/proteomics/scripts/lib_io.py              # IDs, log2 detect, STAR/CCT/PRIDE/CIBERSORT/xCell
methods/proteomics/scripts/fetch_public.py        # LinkedOmics + GDC/S3 STAR
methods/proteomics/scripts/log2_missingness.py    # scale + CLDN4 presence
methods/proteomics/scripts/protein_rna_concordance.py
methods/proteomics/scripts/limma_trend.py         # Python eBayes
methods/proteomics/scripts/limma_trend.R          # Bioconductor limma (optional)
methods/proteomics/scripts/join_deconvolution.py
methods/proteomics/scripts/make_demo_figures.py
methods/proteomics/scripts/run_demo.py            # writes demo/figures + demo/tables
```

```bash
# full demo (downloads CCT + a few STAR files into demo/cache, gitignored)
python3 methods/proteomics/scripts/run_demo.py

# S3 / GDC STAR only (Ensembl ENSG00000184292 / ENSG00000189143)
python3 methods/proteomics/scripts/fetch_public.py --cache-dir methods/proteomics/demo/cache

# local pan-cancer gene-abundance (Ensembl rows), once you have the files:
python3 methods/proteomics/scripts/log2_missingness.py \
  /path/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt \
  --out-json methods/proteomics/demo/tables/pancan_missing.json

# PRIDE proteinGroups
python3 methods/proteomics/scripts/log2_missingness.py proteinGroups.txt --kind pride \
  --out-json methods/proteomics/demo/tables/pride_missing.json
```

`fetch_public.py` lists GDC file UUIDs and the matching `s3://gdc-cptac-phs001287-2-open/{file_id}/…` URIs so you can pull LUAD/LSCC gene-abundance (STAR) without AWS keys.

### 中文

脚本均在 `methods/proteomics/scripts/`。一键演示：`python3 methods/proteomics/scripts/run_demo.py`。  
S3 STAR（Ensembl 两基因）用 `fetch_public.py`。本地泛癌 gene-abundance 或 PRIDE `proteinGroups.txt` 丢给 `log2_missingness.py`。完整 CCT 缓存在 `demo/cache/`，不进 git。

---

## 8. Demo figures / 演示图

| File | Point |
| --- | --- |
| `demo/figures/01_cldn4_missingness.png` | CLDN4 protein row absent; RNA present; NArm NA histogram |
| `demo/figures/02_log2_scale_check.png` | TMT already log2-ratio; LUAD vs LSCC RNA scales differ |
| `demo/figures/03_protein_rna_concordance.png` | TACSTD2 paired; CLDN4 cannot pair |
| `demo/figures/04_limma_tumor_vs_nat.png` | limma on logged TMT, Tumor − NAT |
| `demo/figures/05_join_immune_deconv.png` | TACSTD2 protein vs Immune.Cluster + CD8A RNA |
| `demo/figures/06_naive_not_ici.png` | treatment-naive ≠ ICI response |

Numbers for claims go in `demo/tables/*.tsv` (concordance, limma, join, STAR targets). Do not quote CPTAC statistics as ICI effects.

### 中文

六张图分别对应：CLDN4 行缺失、log2 尺度、蛋白–RNA、limma、免疫对接、未治疗≠ICI。统计数字只在 `demo/tables/`，禁止把 CPTAC 统计写成 ICI 效应。

---

## 9. Decision checklist / 检查清单

### English

1. Did I query **both** Ensembl IDs and symbols?  
2. Is the table already log2? If yes, stop.  
3. Is CLDN4 a missing **row** or a missing **value**?  
4. Am I about to impute protein from RNA? If yes, stop.  
5. Is this limma on log-ratios or limma-trend on log-CPM? Not voom-on-TMT.  
6. Are CIBERSORT/xCell IDs normalized to case IDs?  
7. Do I have an ICI outcome column from an ICI cohort? If no, I may not say “response”.

### 中文

1. Ensembl 与基因符号是否都查了？  
2. 是否已经是 log2？是则不要再变。  
3. CLDN4 是缺行还是缺值？  
4. 是否正用 RNA 填蛋白？是则停止。  
5. 是 log 比值上的 limma，还是 log-CPM 上的 limma-trend？不是对 TMT 做 voom。  
6. CIBERSORT/xCell 的 ID 是否已规范到病例？  
7. 有没有来自 ICI 队列的疗效字段？没有就不能写 “response”。
