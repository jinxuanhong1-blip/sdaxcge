# B2 CCLE all cancers: TACSTD2 vs CLDN4 rank among claudins / 全癌种 CCLE 中 TACSTD2 与 CLDN4 在 claudin 家族中的位次

Analysis commit / 分析提交: `50815d24227d053ad11aa613c2e09f4cdae94e88`  
Run UTC / 运行 UTC: `2026-08-16T19:50:34Z`  
Catalog / 数据目录: `results/w200/B2_CCLE_all/catalog.tsv`  
Seed / 随机种子: `20260816`

**Honesty / 诚实声明:** This slice reports the measured ranks. It does not assume CLDN4 is rank 1. Cell lines are not tumors and have no immune microenvironment; this is not an ICI test.

本切片报告实测位次，不预设 CLDN4 为第 1。细胞系不是肿瘤、无免疫微环境；不能检验 ICI。

---

## English

### Question and design
Primary question: in **all CCLE/DepMap cancer cell lines** (no lineage filter), what is **CLDN4's rank among classical claudins** (`CLDN` + integer, e.g. CLDN1–CLDN25 as present) by Spearman correlation with **TACSTD2**, at RNA and at protein?

Prespecified contrasts:
1. Confirmatory: all-cancers Spearman ρ(TACSTD2, each CLDN\d+), rank of CLDN4 (1 = highest ρ).
2. Confirmatory: all-cancers RNA median-abundance rank of CLDN4 among the same claudins (protein TMT is not treated as copies/cell).
3. Confirmatory: pairwise ρ(TACSTD2, CLDN4) at RNA and protein.
4. Descriptive: genome-wide RNA abundance rank of TACSTD2 and of CLDN4.
5. Exploratory: the same coexpression rank inside Oncotree lineages with n≥20.

Independent biological unit: one cell line (DepMap `ModelID` / CCLE name). No technical-replicate inflation after averaging Nusinow TenPx replicates of the same line.

### Data provenance and accession verification
Databases queried over HTTPS before download (raw responses in `results/w200/B2_CCLE_all/repro/accessions/`):

| Record | Registry | Query | Confirmed title |
|---|---|---|---|
| `10.25452/figshare.plus.25880521` | Figshare+ API `/v2/articles/25880521` | exact article id 25880521 | DepMap 24Q2 Public |
| `MSV000085836` | MassIVE PROXI `accession=MSV000085836&resultType=full` | exact accession | Quantitative Proteomics of the Cancer Cell Line Encyclopedia |

Accepted files (see `catalog.tsv`): DepMap `Model.csv`, `README.txt`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (full file SHA-256 computed in stream; only TACSTD2+CLDN columns retained); Gygi-hosted `protein_quant_current_normalized.csv.gz` and `Table_S1_Sample_Information.xlsx` (processed tables pointed to by the MassIVE record / Gygi CCLE page). No FASTQ. No unverified identifiers entered analysis. Refused files: none (no FASTQ candidates).

RNA values are publisher `log2(TPM+1)` protein-coding gene expression (DepMap 24Q2). Protein values are TMT10 SPS-MS3 normalized log2 ratios (Nusinow et al., *Cell* 2020; 375 lines in the paper).

### Sample size / 样本量
| Stage | RNA | Protein |
|---|---|---|
| Matrix lines available | 1517 | 375 (unique CCLE names after TenPx collapse) |
| Mapped to DepMap ModelID | 1512 | 375 |
| Analyzed, all cancers | 1517 | 375 |
| TACSTD2 non-NA | 1517 | 375 |
| CLDN4 non-NA | 1517 | 197 |
| Pairwise TACSTD2–CLDN4 | 1517 | 197 |
| RNA∩protein mapped lines | 369 | 369 |

Exclusions: RNA rows without a ModelID were none (index is ModelID). Protein columns that did not match `_TenPx##` were ignored if TenPx columns existed. Lineages with n<20 were omitted from exploratory lineage ranks only; they remain inside ALL_CANCERS. Biological unit is the cell line.

### Methods and reproducibility
- Claudin family: gene symbols matching `^CLDN\d+$` present in that matrix. CLDND1 and similar non-integer symbols are excluded. TACSTD2 is not a claudin and is not inserted into the claudin rank list. Confirmatory ranks use pairwise n≥20.
- Missing values: pairwise deletion. Protein genes with duplicate rows were averaged; TenPx replicates of one CCLE name were averaged.
- Test: two-sided Spearman. Software: Python pandas/numpy/scipy/statsmodels (see `repro/pip-freeze.txt`). Seed `20260816` is recorded; no stochastic sampler is used.
- Entry point: `python3 scripts/B2_CCLE_all/00_download.py && python3 scripts/B2_CCLE_all/01_analyze.py`.
- Environment: `results/w200/B2_CCLE_all/repro/pip-freeze.txt` and `conda-explicit.txt` (conda may be unavailable).

### Multiple testing / 多重检验
Family: `CLDN\d+` genes with pairwise n≥20 and defined ρ against TACSTD2 **within one layer × ALL_CANCERS** (RNA and protein are separate families). m_RNA = 24, m_protein = 9. Raw p: two-sided Spearman. Correction: Benjamini–Hochberg. Prespecified q threshold: 0.05. Ties in ρ use minimum rank. Claudins with n<20 remain in the coverage table but are outside the confirmatory rank/BH family. Lineage tests are exploratory and are **not** folded into the confirmatory BH family.

Confirmatory RNA family: 18 / 24 claudins with q<0.05.  
Confirmatory protein family: 6 / 9 claudins with q<0.05.

### Results
**Confirmatory, all cancers — coexpression rank (primary B2).**

RNA (DepMap 24Q2, n=1517 lines): TACSTD2 vs CLDN4 ρ = 0.714, p = 7.82e-237, n_pairwise = 1517. CLDN4 rank among claudins by ρ = **2 / 24**. Top claudin = **CLDN7** (ρ = 0.726). CLDN4 is NOT the top TACSTD2-coexpressed claudin at RNA (rank 2 of 24; top is CLDN7). CLDN7 and CLDN4 are close (Δρ ≈ 0.012).

Protein (Nusinow 2020, n=375 lines; CLDN4 quantified in 197/375): TACSTD2 vs CLDN4 ρ = 0.619, p = 3.47e-22, n_pairwise = 197. CLDN4 rank among confirmatory-family claudins by ρ = **1 / 9**. Top claudin = **CLDN4** (ρ = 0.619). CLDN4 is the top TACSTD2-coexpressed claudin at protein among confirmatory-family claudins. Only claudins with pairwise n≥20 enter this rank; several CLDN proteins are mostly missing.

**Confirmatory — RNA abundance rank among claudins (median log2(TPM+1)).**  
CLDN4 median-abundance rank = 4 / 24. Protein TMT medians are relative to a common reference and are **not** copies/cell; they are tabled but not used as an abundance claim.

**Descriptive — genome-wide RNA abundance rank (mean log2(TPM+1)).**  
TACSTD2 rank 8367 / 19193; CLDN4 rank 5906 / 19193. Protein genome-wide TMT-median ranks are in `tables/genomewide_protein_abundance_rank.tsv` and must not be read as molecular abundance.

**RNA–protein coupling on mapped overlap (n=369 lines).**  
TACSTD2 RNA vs protein ρ = 0.733, p = 2.42e-63, n = 369.  
CLDN4 RNA vs protein ρ = 0.798, p = 2.49e-44, n = 195.

These are cell-line expression associations. They do not measure a TACSTD2–CLDN4 junction, immune exclusion, or ICI resistance.

### Limitations
- CCLE/DepMap lines are in vitro, usually monoclonal, and lack stroma, immune cells, and in vivo barrier anatomy.
- Nusinow protein covers ~375 lines, far fewer than DepMap RNA; membrane proteins including some claudins are frequently missing (see missingness figure). TMT ratios are relative, not copies/cell.
- DepMap 24Q2 RNA mixes historical CCLE and later models; we did not batch-correct beyond the publisher matrix.
- Lineage ranks are exploratory (n≥20 gate only).
- No ICI labels exist for this slice. Do not read ρ as a resistance biomarker.

## 中文

### 研究问题与设计
主要问题：在**全部** CCLE/DepMap 癌细胞系（不按谱系过滤）中，经典 claudin（符号为 `CLDN`+整数）里，**CLDN4 与 TACSTD2 的 Spearman 相关位次**在 RNA 与蛋白层分别是多少？

预设比较：
1. 验证性：全癌种 ρ(TACSTD2, 各 CLDN\d+)，报告 CLDN4 位次（1 = ρ 最高）。
2. 验证性：全癌种 RNA 上 CLDN4 在同一家族中的中位丰度位次（蛋白 TMT 不视为每细胞拷贝数）。
3. 验证性：TACSTD2–CLDN4 成对 ρ（RNA 与蛋白）。
4. 描述性：TACSTD2 与 CLDN4 的全基因组 RNA 丰度位次。
5. 探索性：Oncotree 谱系内（n≥20）重复共表达位次。

独立生物学单位：一条细胞系（DepMap `ModelID` / CCLE 名）。Nusinow 同一系的 TenPx 重复先平均，不把技术重复当独立样本。

### 数据来源与登录号核验
下载前经 HTTPS 核验（原始响应见 `results/w200/B2_CCLE_all/repro/accessions/`）：

| 记录 | 数据库 | 查询 | 确认标题 |
|---|---|---|---|
| `10.25452/figshare.plus.25880521` | Figshare+ API `/v2/articles/25880521` | 精确文章 id 25880521 | DepMap 24Q2 Public |
| `MSV000085836` | MassIVE PROXI `accession=MSV000085836&resultType=full` | 精确登录号 | Quantitative Proteomics of the Cancer Cell Line Encyclopedia |

接受文件见 `catalog.tsv`：DepMap `Model.csv`、`README.txt`、`OmicsExpressionProteinCodingGenesTPMLogp1.csv`（流式计算全文 SHA-256，仅保留 TACSTD2+CLDN 列）；Gygi 托管的 `protein_quant_current_normalized.csv.gz` 与 `Table_S1_Sample_Information.xlsx`（MassIVE 记录/Gygi CCLE 页指向的处理后表格）。无 FASTQ。分析未纳入未经核验的标识符。拒绝文件：无。

RNA 为 DepMap 24Q2 发布的蛋白编码基因 `log2(TPM+1)`。蛋白为 TMT10 SPS-MS3 归一化 log2 比值（Nusinow 等，*Cell* 2020；论文中 375 系）。

### Sample size / 样本量
| 阶段 | RNA | 蛋白 |
|---|---|---|
| 矩阵中的系 | 1517 | 375（TenPx 合并后的唯一 CCLE 名） |
| 映射到 DepMap ModelID | 1512 | 375 |
| 全癌种分析 | 1517 | 375 |
| TACSTD2 非缺失 | 1517 | 375 |
| CLDN4 非缺失 | 1517 | 197 |
| TACSTD2–CLDN4 成对 | 1517 | 197 |
| RNA∩蛋白已映射系 | 369 | 369 |

排除：RNA 以 ModelID 为索引，无空 ID。若存在 TenPx 列，则忽略非 TenPx 样本列。n<20 的谱系只从探索性谱系位次中去掉，仍留在 ALL_CANCERS。生物学单位是细胞系。

### 方法与可复现性
- Claudin 家族：该矩阵中匹配 `^CLDN\d+$` 的基因。CLDND1 等非整数符号不计入。TACSTD2 不是 claudin，不进入 claudin 位次表。验证性位次要求成对 n≥20。
- 缺失：成对删除。蛋白重复基因行取平均；同一 CCLE 名的 TenPx 重复取平均。
- 检验：双侧 Spearman。软件：pandas/numpy/scipy/statsmodels（见 `repro/pip-freeze.txt`）。种子 `20260816` 已记录；本分析无随机抽样。
- 入口：`python3 scripts/B2_CCLE_all/00_download.py && python3 scripts/B2_CCLE_all/01_analyze.py`。
- 环境：`results/w200/B2_CCLE_all/repro/pip-freeze.txt` 与 `conda-explicit.txt`（可能无 conda）。

### Multiple testing / 多重检验
检验族：在**单一层次 × ALL_CANCERS** 内，成对 n≥20 且 ρ 可定义的 `CLDN\d+` 对 TACSTD2（RNA 与蛋白分族）。m_RNA = 24，m_protein = 9。原始 p：双侧 Spearman。校正：Benjamini–Hochberg。预设 q 阈值：0.05。ρ 并列时取最小位次。n<20 的 claudin 仍在覆盖表中，但不进入验证性位次/BH 族。谱系检验为探索性，**不**并入验证性 BH 族。

验证性 RNA 族：18 / 24 个 claudin q<0.05。  
验证性蛋白族：6 / 9 个 claudin q<0.05。

### 结果
**验证性，全癌种 — 共表达位次（B2 主终点）。**

RNA（DepMap 24Q2，n=1517 系）：TACSTD2 vs CLDN4 ρ = 0.714，p = 7.82e-237，成对 n = 1517。CLDN4 按 ρ 在 claudin 中的位次 = **2 / 24**。最高 claudin = **CLDN7**（ρ = 0.726）。在全癌种 RNA 中，CLDN4 并不是与 TACSTD2 共表达最高的 claudin（第 2 / 24；最高为 CLDN7）。 CLDN7 与 CLDN4 很接近（Δρ ≈ 0.012）。

蛋白（Nusinow 2020，n=375 系；CLDN4 定量 197/375）：TACSTD2 vs CLDN4 ρ = 0.619，p = 3.47e-22，成对 n = 197。CLDN4 按 ρ 在验证性家族中的位次 = **1 / 9**。最高 claudin = **CLDN4**（ρ = 0.619）。在验证性蛋白家族中，CLDN4 是与 TACSTD2 共表达最高的 claudin。 仅成对 n≥20 的 claudin 进入该位次；多种 CLDN 蛋白大量缺失。

**验证性 — RNA 在 claudin 家族中的中位丰度位次（中位 log2(TPM+1)）。**  
CLDN4 中位丰度位次 = 4 / 24。蛋白 TMT 中位数是相对桥样的比值，**不是**每细胞拷贝数；已列表但不作为丰度结论。

**描述性 — 全基因组 RNA 丰度位次（平均 log2(TPM+1)）。**  
TACSTD2 第 8367 / 19193；CLDN4 第 5906 / 19193。蛋白全基因组 TMT 中位位次见 `tables/genomewide_protein_abundance_rank.tsv`，不得读成分子丰度。

**已映射重叠系上的 RNA–蛋白耦合（n=369）。**  
TACSTD2 RNA vs 蛋白 ρ = 0.733，p = 2.42e-63，n = 369。  
CLDN4 RNA vs 蛋白 ρ = 0.798，p = 2.49e-44，n = 195。

以上均为细胞系表达关联，不能证明 TACSTD2–CLDN4 连接事件、免疫排斥或 ICI 耐药。

### 局限性
- CCLE/DepMap 细胞系为体外、多为单克隆，缺乏间质、免疫细胞和体内屏障结构。
- Nusinow 蛋白约 375 系，远少于 DepMap RNA；包括部分 claudin 在内的膜蛋白经常缺失（见缺失图）。TMT 比值为相对定量，不是每细胞拷贝数。
- DepMap 24Q2 RNA 混合历史 CCLE 与后续模型；除发布矩阵外未再做批次校正。
- 谱系位次仅为探索性（仅 n≥20）。
- 本切片无 ICI 标签。不能把 ρ 读成耐药生物标志物。
