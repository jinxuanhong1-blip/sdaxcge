# CPTAC LSCC / LUSC — TACSTD2 and CLDN4 protein + RNA + immune

Parallel slice. Outputs only under `notes/grok_cptac_lscc/`, `scripts/grok_cptac_lscc/`, `results/grok_cptac_lscc/`.

---

## English

### Question
In treatment-naive CPTAC lung squamous cell carcinoma (LSCC; TCGA synonym LUSC), how do **TACSTD2** (`ENSG00000184292`, TROP2, UniProt P09758) and **CLDN4** (`ENSG00000189143`) relate to immune contexture at protein and RNA level? This cohort has **no ICI response labels**. OS/PFS here are not immunotherapy outcomes.

### Data (open S3 freeze, verified HTTP 200)
CPTAC pan-cancer `data_freeze_v1.2_reorganized` on `cptac-pancancer-data` (us-west-2). Filenames taken from the LinkedOmics CPTAC-pancan-LSCC table, not invented.

| Matrix | n | Features | Units |
|---|---:|---:|---|
| Protein tumor | 108 | 12,760 genes | log2 reference-intensity |
| Protein NAT | 99 | 12,760 | same |
| RNA tumor | 108 | 60,669 genes | log2 RSEM coding UQ 1500 |
| RNA NAT | 94 | 60,669 | same |
| Phenotype | 108 | CIBERSORT / ESTIMATE / xCell / HALLMARK / PROGENy / TMB | freeze-derived |
| Survival | 108 | OS / PFS days + event | treatment-naive follow-up |
| Meta | 108 | age, sex, stage, smoking, mutation flags | skip `data_type` row |

Paper: Satpathy et al., *Cell* 2021, PMID 34358469. Quote: newly diagnosed resected LSCC that **received no prior chemotherapy or radiotherapy**. PDC proteome study: PDC000234. Clinical mix in this freeze: 86 male / 22 female; Stage I 39, II 44, III 21, IV 1, missing 3.

Standalone `LSCC_xcell.txt` / `LSCC_cibersort.txt` return HTTP 403. Immune scores live inside `LSCC_phenotype.txt`.

### Identifiers and missingness
Rows matched by Ensembl **prefix** (same IDs as the LUAD sibling).

| Layer | Gene | Row | Quantified | NA |
|---|---|---|---:|---:|
| Protein tumor | TACSTD2 | `ENSG00000184292.7` | 108 | **0** |
| Protein tumor | CLDN4 | `ENSG00000189143.9` | **78** | **30 (27.8%)** |
| RNA tumor | both | same versioned IDs | 108 | 0 |

CLDN4 protein missingness is real TMT dropout, not a join error. All CLDN4 protein tests use pairwise-complete n=78 unless noted.

### Statistics
Spearman (n&lt;4 → NA); two-sided Mann–Whitney U + rank-biserial; paired Wilcoxon (n&lt;6 → NA); BH-FDR within each table; signature score = mean of per-gene z-scores with `used/total` coverage; KM + log-rank (median split) and Cox on the continuous value. Scripts that actually ran: `scripts/grok_cptac_lscc/00_download_freeze.py`, `01_analyze_lscc.py`.

RNA signatures were complete (GEP 18/18, IFNG 6/6, TLS 12/12, CYT 2/2). Protein signatures were partial (GEP 16/18, TLS 9/12, cytotoxic 4/6).

### Results

**1. Tumor vs NAT.** TACSTD2 is higher in tumor than paired NAT at both layers. CLDN4 RNA is higher; CLDN4 protein is not.

| Gene | Layer | Paired n | Median tumor−NAT | Wilcoxon p | FDR |
|---|---|---:|---:|---:|---:|
| TACSTD2 | protein | 99 | +0.388 | 1.19×10⁻⁵ | 1.58×10⁻⁵ |
| TACSTD2 | RNA | 94 | +1.295 | 6.44×10⁻¹¹ | 2.58×10⁻¹⁰ |
| CLDN4 | protein | 72 | +0.192 | 0.117 | 0.117 |
| CLDN4 | RNA | 94 | +0.575 | 3.82×10⁻⁷ | 7.65×10⁻⁷ |

**2. Protein–RNA concordance and co-expression.**

| Comparison | n | Spearman ρ | p | FDR |
|---|---:|---:|---:|---:|
| TACSTD2 protein vs RNA | 108 | **0.864** | 2.35×10⁻³³ | 9.40×10⁻³³ |
| CLDN4 protein vs RNA | 78 | **0.580** | 2.55×10⁻⁸ | 5.10×10⁻⁸ |
| TACSTD2 vs CLDN4 protein | 78 | 0.081 | 0.483 | 0.483 |
| TACSTD2 vs CLDN4 RNA | 108 | **0.446** | 1.30×10⁻⁶ | 1.73×10⁻⁶ |

TACSTD2 protein tracks RNA tightly. TACSTD2 and CLDN4 co-vary at RNA but **not** at protein in the 78 samples with both proteins.

**3. Immune contexture — protein (strongest signal).**
Among 106 protein×phenotype tests, **26** have FDR&lt;0.05. The direction is inverse immune / stroma, driven mainly by **CLDN4 protein**.

CLDN4 protein (n=78), selected FDR-significant hits:

| Score | ρ | p | FDR |
|---|---:|---:|---:|
| ESTIMATE ImmuneScore | −0.432 | 7.88×10⁻⁵ | 0.00166 |
| ESTIMATE ESTIMATEScore | −0.453 | 3.12×10⁻⁵ | 0.00165 |
| xCell immune_score | −0.400 | 2.90×10⁻⁴ | 0.00304 |
| xCell CD8 | −0.387 | 4.66×10⁻⁴ | 0.00408 |
| xCell B cell | −0.407 | 2.18×10⁻⁴ | 0.00286 |
| HALLMARK IFNγ | −0.353 | 1.52×10⁻³ | 0.00840 |
| HALLMARK inflammatory | −0.366 | 9.89×10⁻⁴ | 0.00577 |
| HALLMARK allograft rejection | −0.436 | 6.64×10⁻⁵ | 0.00166 |

CLDN4 protein vs RNA-derived signatures (all FDR&lt;0.05, n=78): TLS ρ=−0.518; Immune_general ρ=−0.498; IFNG_6gene ρ=−0.463; GEP T-cell inflamed ρ=−0.461; CD8_Tcell ρ=−0.401.

TACSTD2 protein (n=108) is **not** significantly anti-correlated with ESTIMATE ImmuneScore (ρ=−0.071, p=0.465, FDR=0.642) or CIBERSORT CD8 (ρ=−0.055, p=0.574). FDR-significant TACSTD2 protein associations are stromal / TGF-β inverse and a modest JAK–STAT positive:

| Score | ρ | p | FDR |
|---|---:|---:|---:|
| xCell stroma | −0.345 | 2.52×10⁻⁴ | 0.00294 |
| PROGENy TGFβ | −0.340 | 3.23×10⁻⁴ | 0.00308 |
| xCell CAF | −0.319 | 7.59×10⁻⁴ | 0.00531 |
| CIBERSORT M2 | −0.282 | 3.15×10⁻³ | 0.0158 |
| HALLMARK EMT | −0.256 | 7.55×10⁻³ | 0.0345 |
| ESTIMATE StromalScore | −0.251 | 8.86×10⁻³ | 0.0383 |
| PROGENy JAK–STAT | +0.250 | 9.11×10⁻³ | 0.0383 |

TACSTD2 protein vs GEP: ρ=−0.042, p=0.664 (n=108). This does **not** reproduce a Bessede-like “high TROP2 = T-cell excluded” pattern in treatment-naive LSCC protein data.

**4. Immune contexture — RNA.**
106 RNA×phenotype tests: **0** FDR&lt;0.05 (21 nominal p&lt;0.05). Strongest nominal hits are TACSTD2 RNA vs xCell stroma (ρ=−0.311, p=0.00106, FDR=0.058) and CAF (ρ=−0.307, p=0.00125, FDR=0.058). CLDN4 RNA vs TGFB_exclusion signature: ρ=−0.289, p=0.00239, FDR=0.038 (only RNA-signature FDR hit).

Median-split MWU on a pre-specified 10-score key set: **0** FDR&lt;0.05 (continuous Spearman is the better-powered test).

**5. Survival (not ICI).**
Complete-case OS n=94 / 23 events; PFS n=92 / 20 events; CLDN4 protein OS n=66 / 14 events. Every log-rank and Cox test is null after FDR (all log-rank p≥0.276, all Cox p≥0.313). Underpowered. Do not interpret as “no prognostic effect.”

### Interpretation vs Bessede
Bessede et al. *Clin Cancer Res* 2024 (PMID 38048058) reported high TACSTD2 with atezolizumab primary resistance and less T-cell infiltration (OAK/POPLAR; RNA often EGA-controlled). This CPTAC LSCC slice is **treatment-naive** and cannot test ICI resistance. In this dataset the immune-cold association at **protein** is carried by **CLDN4**, not TACSTD2. TACSTD2 protein instead tracks lower stroma / TGF-β / CAF and slightly higher JAK–STAT.

### Caveats
- No ICI labels; do not treat OS/PFS as IO endpoints.
- CLDN4 protein NA=30/108; those tumors are dropped from protein tests.
- Freeze deconvolution is RNA-derived; protein–immune correlations mix layers.
- CIBERSORT CD4 naive is constant among CLDN4-complete samples (Spearman undefined; recorded as NA).
- Multiple testing is table-wise BH, not experiment-wide.
- LSCC ≠ LUAD; do not pool with the LUAD sibling without a pre-specified interaction test.

### How to rerun
```bash
pip install -r scripts/grok_cptac_lscc/requirements.txt
python3 scripts/grok_cptac_lscc/00_download_freeze.py
python3 scripts/grok_cptac_lscc/01_analyze_lscc.py
```
Matrices: `data/grok_cptac_lscc/` (not committed). Tables/figures: `results/grok_cptac_lscc/`.

---

## 中文

### 问题
在 CPTAC **治疗初治**肺鳞癌（冻存目录名 **LSCC**，TCGA 同义名 **LUSC**）中，**TACSTD2**（`ENSG00000184292`，TROP2）和 **CLDN4**（`ENSG00000189143`）的蛋白 / RNA 水平与免疫微环境如何相关？该队列**没有 ICI 疗效标签**。此处 OS/PFS **不是**免疫治疗结局。

### 数据
开放 S3 冻存 `data_freeze_v1.2_reorganized`（`cptac-pancancer-data`，us-west-2）。文件名来自 LinkedOmics CPTAC-pancan-LSCC 下载表，并经 HTTP 200 核实。肿瘤蛋白 108 例 × 12760 基因；癌旁蛋白 99；肿瘤 RNA 108 × 60669；癌旁 RNA 94。免疫分数在 `LSCC_phenotype.txt`（CIBERSORT / ESTIMATE / xCell / HALLMARK / PROGENy）。Satpathy 等 *Cell* 2021（PMID 34358469）：新诊断手术切除、**术前未接受化疗或放疗**。PDC 蛋白组：PDC000234。

### 缺失
TACSTD2 蛋白 108/108（NA=0）。**CLDN4 蛋白仅 78/108（NA=30，27.8%）**，为 TMT 未定量，不是拼接错误。RNA 两基因均 108/108。Ensembl 按前缀匹配（与 LUAD 切片同一 ID）。

### 主要结果（真实统计，未编造）
1. **肿瘤 vs 配对癌旁**：TACSTD2 蛋白配对 n=99，中位差 +0.388，Wilcoxon p=1.19×10⁻⁵；RNA 配对 n=94，中位差 +1.295，p=6.44×10⁻¹¹。CLDN4 RNA 升高（p=3.82×10⁻⁷）；CLDN4 蛋白不显著（p=0.117）。
2. **蛋白–RNA**：TACSTD2 ρ=0.864（n=108，p=2.35×10⁻³³）；CLDN4 ρ=0.580（n=78，p=2.55×10⁻⁸）。TACSTD2 与 CLDN4 在 RNA 共表达（ρ=0.446，p=1.30×10⁻⁶），在蛋白层不共表达（ρ=0.081，p=0.483，n=78）。
3. **免疫（蛋白层最强）**：106 项检验中 26 项 FDR&lt;0.05，方向为免疫/基质负相关，主要由 **CLDN4 蛋白**驱动。CLDN4 蛋白 vs ESTIMATE ImmuneScore ρ=−0.432（p=7.88×10⁻⁵，FDR=0.00166，n=78）；vs xCell CD8 ρ=−0.387；vs GEP ρ=−0.461；vs TLS ρ=−0.518。TACSTD2 蛋白与 ImmuneScore（ρ=−0.071，p=0.465）和 CIBERSORT CD8（ρ=−0.055，p=0.574）均不显著；FDR 显著项为基质 / TGF-β / CAF 负相关，以及 JAK–STAT 弱正相关（ρ=+0.250，FDR=0.038）。
4. **免疫（RNA 层）**：106 项中 **0** 项 FDR&lt;0.05。名义最强为 TACSTD2 RNA vs 基质/CAF（FDR=0.058）。
5. **生存**：OS 23 事件 / PFS 20 事件，检验全部不显著。把握度不足，不能写成“无预后意义”。**不是 ICI 生存。**

### 与 Bessede 的关系
Bessede 等（PMID 38048058）在 ICI 队列报告高 TACSTD2 与阿替利珠单抗原发耐药、T 细胞浸润减少。本切片为治疗初治 LSCC，**不能检验 ICI 耐药**。本数据中蛋白层“免疫冷”信号在 **CLDN4**，不在 TACSTD2。

### 局限
无 ICI 标签；CLDN4 蛋白缺失 30 例；表型免疫分数来自 RNA 反卷积；多重检验为表内 BH；不可与 LUAD 切片未经预设交互检验而合并。

### 复现
见上文 English “How to rerun”。图表：`results/grok_cptac_lscc/figures/`；完整表：`results/grok_cptac_lscc/tables/`。
