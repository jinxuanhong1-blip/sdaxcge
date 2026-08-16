# CPTAC LUAD: TACSTD2 / CLDN4 protein + RNA + immune

**Slice:** `notes|scripts|results/grok_cptac_luad/` only.
**Cohort:** CPTAC lung adenocarcinoma, pan-cancer data freeze v1.2 (open S3).
**ICI labels:** none. Gillette et al. *Cell* 2020 prospectively collected **treatment-naive** tumors. OS/PFS in this freeze are surgical follow-up, not immunotherapy outcomes. This slice cannot test Bessede et al. (atezolizumab resistance).

---

## English

### What was tested

Open freeze matrices (HTTP 200, filenames from the LinkedOmics CPTAC-pancan-LUAD index, same S3 prefix `cptac-pancancer-data / data_freeze_v1.2_reorganized/LUAD/`):

| Layer | File | Size |
|---|---|---|
| RNA tumor / NAT | `LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_{Tumor,Normal}.txt` | 34 / 31 MB |
| Protein tumor / NAT | `LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_{Tumor,Normal}.txt` | 20 / 18 MB |
| Immune / pathway scores | `LUAD_phenotype.txt` (CIBERSORT, ESTIMATE, xCell, PROGENy, HALLMARK, TMB, purity) | 205 KB |
| Clinical + survival | `LUAD_meta.txt`, `LUAD_survival.txt` | 12 / 2.5 KB |

Targets: **TACSTD2** `ENSG00000184292.7`, **CLDN4** `ENSG00000189143.9`.
n = **110 tumors** with complete RNA + protein + phenotype IDs; **101 matched NAT** (RNA and protein).

Tests (pre-specified): (1) paired tumor vs NAT Wilcoxon; (2) RNA–protein Spearman; (3) TACSTD2–CLDN4 co-expression; (4) Spearman vs freeze immune scores, BH-FDR within predictor; (5) the same after residualizing on `WES_purity`; (6) median-split Mann–Whitney vs ImmuneScore / CD8 / IFN-γ / CYT / GEP; (7) clinical MW (smoking, STK11/KEAP1/EGFR/KRAS/TP53, sex, stage); (8) OS/PFS log-rank + Cox per 1 SD. No fabricated accessions or p-values.

### Coverage (honest missingness)

| Feature | Tumor non-NA | NAT non-NA |
|---|---|---|
| TACSTD2 RNA | 110 / 110 | 101 / 101 |
| CLDN4 RNA | 110 / 110 | 101 / 101 |
| TACSTD2 protein | 110 / 110 | 101 / 101 |
| CLDN4 protein | **79 / 110 (31 NA)** | **72 / 101 (29 NA)** |

CLDN4 protein analyses are pairwise-complete (n=79 tumor, n=71 paired NAT). Do not treat CLDN4 protein as fully quantified.

Phenotype immune scores are **RNA-derived** (Li et al. *Cell Syst.* 2023; immunedeconv CIBERSORT + xCell). Protein-vs-score tests are therefore cross-layer, not circular with TACSTD2/CLDN4 protein.

### Tumor vs matched NAT (paired Wilcoxon)

| Gene | Layer | n pairs | median tumor | median NAT | Δ (T−N) | p |
|---|---|---|---|---|---|---|
| TACSTD2 | RNA | 101 | 12.24 | 11.71 | **+0.46** | **1.61×10⁻⁵** |
| CLDN4 | RNA | 101 | 12.80 | 11.48 | **+1.23** | **3.88×10⁻¹⁷** |
| TACSTD2 | protein | 101 | 26.53 | 26.94 | **−0.42** | **2.41×10⁻⁹** |
| CLDN4 | protein | 71 | 23.04 | 23.20 | **−0.18** | **0.020** |

RNA is higher in tumor for both genes. **Protein goes the other way** (tumor lower than NAT), especially TACSTD2. This RNA–protein discordance is a real result, not a plotting artifact.

### RNA–protein coupling and co-expression

- TACSTD2 RNA vs protein: Spearman **ρ = 0.666**, p = **1.93×10⁻¹⁵**, n = 110.
- CLDN4 RNA vs protein: ρ = 0.174, p = 0.126, n = 79 (not significant; missingness + weaker coupling).
- TACSTD2 vs CLDN4 RNA: **ρ = 0.392**, p = **2.30×10⁻⁵**, n = 110.
- TACSTD2 vs CLDN4 protein: **ρ = 0.275**, p = **0.014**, n = 79.

### Immune associations (primary endpoints)

Primary list: ESTIMATE Immune/Stromal, CIBERSORT CD8 / Treg / M2 / neutrophil, xCell CD8 / immune / Treg, HALLMARK IFN-γ / IFN-α / inflammatory / TGF-β / EMT, PROGENy JAK-STAT / TGFb / NFκB, TMB, WES_purity. FDR is BH within each predictor × adjustment.

**TACSTD2 protein (n=110) — supports “high TROP2, less T-cell / immune score” in this treatment-naive cohort:**

| Target | ρ | p | q (BH) |
|---|---|---|---|
| xCell_immune_score | **−0.309** | **0.00100** | **0.024** |
| xCell_T_cell_CD8+ | **−0.289** | **0.00219** | **0.034** |
| ESTIMATE_ImmuneScore | −0.185 | 0.052 | 0.217 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | −0.151 | 0.116 | 0.315 |
| CIBERSORT_T_cell_CD8+ | −0.133 | 0.164 | 0.380 |

After residualizing on WES_purity: xCell immune ρ = −0.269, p = 0.0049, q = 0.072; xCell CD8 ρ = −0.257, p = 0.0073, q = 0.093. Direction holds; FDR is marginal once purity is removed.

Median-split TACSTD2 protein (55 high / 55 low):

- xCell CD8: high median 0.011 vs low 0.028, rank-biserial = −0.335, p = 0.0025, q = 0.059
- ESTIMATE ImmuneScore: 7862 vs 8712, r = −0.264, p = 0.017, q = 0.072
- GEP T-cell inflamed (RNA z-mean, 18/18 genes): r = −0.265, p = 0.017, q = 0.072
- CYT (GZMA+PRF1 RNA): r = −0.235, p = 0.034, q = 0.090

**TACSTD2 RNA (n=110):** same CD8 direction but weaker and not FDR-significant among primary scores. CIBERSORT CD8 ρ = −0.201, p = 0.035, q = 0.176. TMB ρ = −0.233, p = 0.015, q = 0.098 (purity-adjusted ρ = −0.253, p = 0.0086, q = 0.060).

**CLDN4 protein (n=79):** no primary immune test reaches q < 0.05. Nominal trends: HALLMARK IFN-γ ρ = −0.199, p = 0.079; xCell CD8 ρ = −0.195, p = 0.085. Median-split (40/39) vs CYT RNA p = 0.018, q = 0.072; vs GEP p = 0.019, q = 0.072; vs xCell CD8 p = 0.026, q = 0.077.

**CLDN4 RNA (n=110):** PROGENy JAK-STAT ρ = 0.273, p = 0.0039, q = 0.069; xCell Treg ρ = −0.264, p = 0.0053, q = 0.069; CIBERSORT M2 ρ = 0.254, p = 0.0075, q = 0.081; CIBERSORT CD8 ρ = −0.225, p = 0.018, q = 0.157. ESTIMATE ImmuneScore is null (ρ = −0.079, p = 0.41).

### Clinical (not ICI)

STK11-mutant tumors have **lower** TACSTD2 RNA (median 11.69 vs 12.34; MW p = 6.8×10⁻⁴, q = 0.0055, n_mut=18 / n_wt=90) and **lower** CLDN4 RNA (12.40 vs 13.05; p = 5.8×10⁻⁴, q = 0.0055). EGFR-mutant tumors have **higher** CLDN4 RNA (13.24 vs 12.60; p = 1.31×10⁻⁶, q = 4.2×10⁻⁵, n=39/69). TP53-mutant: higher CLDN4 RNA (p = 4.5×10⁻⁴, q = 0.0055). Protein-level mutation tests are mostly null (KEAP1 vs TACSTD2 protein p = 0.027, q = 0.12). Smoking and stage are not associated with TACSTD2 protein.

### Survival (underpowered; not ICI)

OS events = 23 among 105 cases with time; PFS events = 23 / 104. All median-split log-rank p > 0.05. Closest: CLDN4 protein OS p = 0.087 (15 events / 75); TACSTD2 protein PFS Cox HR per SD = 1.51, p = 0.071. **Do not interpret as ICI survival.**

### How this sits next to Bessede 2024

Bessede et al. (*Clin Cancer Res*, PMID 38048058) reported high TACSTD2 with atezolizumab primary resistance and less T-cell infiltration in OAK/POPLAR (often EGA-controlled RNA). CPTAC LUAD **cannot** test ICI response. What it can test is the infiltration half in untreated LUAD: **TACSTD2 protein is inversely associated with xCell immune score and xCell CD8 (q < 0.05 among primary endpoints)**. That is consistent in direction, at protein level, in a different (treatment-naive) setting. It is not a replication of the ICI-resistance claim.

### Caveats

1. No ICI labels. Do not call these tumors “hot/cold on IO.”
2. CLDN4 protein missing in 31/110 tumors — underpowered and possibly MNAR (low-abundance drop-out).
3. Freeze immune scores come from RNA; protein-vs-score is the cleaner test for TACSTD2/CLDN4 protein.
4. Tumor vs NAT protein decrease vs RNA increase means bulk RNA ranking is not a substitute for TROP2 protein.
5. OS/PFS events are few; null survival is inconclusive.
6. Multiple testing: primary-list FDR is the inferential claim; the full Spearman table is exploratory.
7. Files >2 GB were not used; isoform/phospho/methylation/CNV/MAF were out of this slice.

### How to rerun

```bash
python3 -m pip install pandas numpy scipy matplotlib statsmodels
python3 scripts/grok_cptac_luad/00_download.py --outdir data/grok_cptac_luad
python3 scripts/grok_cptac_luad/01_analyze.py --data data/grok_cptac_luad --outdir results/grok_cptac_luad
```

Figures: `results/grok_cptac_luad/figures/fig1_tumor_vs_nat.png` … `fig8_missingness.png`.
Tables: `results/grok_cptac_luad/tables/*.tsv`.

---

## 中文

### 测了什么

开放 S3 freeze v1.2（LinkedOmics CPTAC-pancan-LUAD 页面给出的文件名；逐个 HEAD 确认 HTTP 200）。只做 **蛋白 + RNA + 免疫评分**。队列是 Gillette 等 *Cell* 2020 的 **初治、手术切除** LUAD，**没有 ICI 疗效标签**。这里的 OS/PFS 是术后随访，不是免疫治疗生存。

靶基因：TACSTD2 `ENSG00000184292.7`，CLDN4 `ENSG00000189143.9`。肿瘤 n=110（RNA/蛋白/表型 ID 对齐），配对癌旁 n=101。

### 覆盖与缺失

TACSTD2 的 RNA 与蛋白在 110 例肿瘤中均无缺失。**CLDN4 蛋白肿瘤 31/110 为 NA**（癌旁 29/101 NA）。CLDN4 蛋白相关统计一律 pairwise-complete，不能当成“测全了”。表型里的 CIBERSORT / xCell / ESTIMATE 来自 RNA 反卷积，不是蛋白反卷积。

### 肿瘤 vs 配对癌旁（配对 Wilcoxon）

- TACSTD2 RNA：肿瘤更高，n=101，Δ=+0.46，p=1.61×10⁻⁵。
- CLDN4 RNA：肿瘤更高，Δ=+1.23，p=3.88×10⁻¹⁷。
- TACSTD2 **蛋白：肿瘤更低**，Δ=−0.42，p=2.41×10⁻⁹。
- CLDN4 蛋白：肿瘤略低，n=71，Δ=−0.18，p=0.020。

RNA 上调、蛋白下调，方向相反，需要写进结论，不能只用 RNA 代替 TROP2 蛋白。

### RNA–蛋白 与 共表达

- TACSTD2 RNA–蛋白 ρ=0.666，p=1.93×10⁻¹⁵，n=110。
- CLDN4 RNA–蛋白 ρ=0.174，p=0.126，n=79（不显著）。
- TACSTD2–CLDN4 RNA ρ=0.392，p=2.30×10⁻⁵；蛋白 ρ=0.275，p=0.014，n=79。

### 免疫（预先指定终点）

**TACSTD2 蛋白**与 xCell immune score（ρ=−0.309，p=0.00100，q=0.024）及 xCell CD8（ρ=−0.289，p=0.00219，q=0.034）负相关；ESTIMATE ImmuneScore 同向但不过 FDR（ρ=−0.185，p=0.052，q=0.217）。校正 WES 纯度后方向仍在，q 约 0.07–0.09。中位数分组：高 TACSTD2 蛋白的 xCell CD8 更低（p=0.0025，q=0.059），GEP / CYT RNA 签名也更低（p=0.017 / 0.034）。

**TACSTD2 RNA** 与 CIBERSORT CD8 名义负相关（ρ=−0.201，p=0.035，q=0.176），FDR 不过关。与 TMB 负相关（ρ=−0.233，p=0.015，q=0.098）。

**CLDN4 蛋白**无主键免疫终点 q<0.05；IFN-γ / xCell CD8 仅有名义趋势。**CLDN4 RNA** 与 JAK-STAT 正、xCell Treg 负、M2 正（q≈0.07–0.08），与 ESTIMATE ImmuneScore 无关。

### 临床（非 ICI）

STK11 突变：TACSTD2 RNA 与 CLDN4 RNA 均更低（q=0.0055）。EGFR 突变：CLDN4 RNA 更高（q=4.2×10⁻⁵）。蛋白水平突变关联基本不显著。吸烟、分期与 TACSTD2 蛋白无关。

### 生存

OS 事件仅 23/105，全部 log-rank p>0.05。不能当成 ICI 生存，也不能把阴性写成“无预后价值”。

### 与 Bessede 2024 的关系

Bessede 等报告高 TACSTD2 对应阿替利珠单抗原发耐药和更少 T 细胞浸润（OAK/POPLAR，RNA 常为 EGA 受控）。本队列**无法检验 ICI 疗效**。能检验的是初治 LUAD 里“浸润更少”这一半：**TACSTD2 蛋白与 xCell 免疫评分 / CD8 负相关（主键 FDR q<0.05）**，方向一致，但不是耐药结论的重复。

### 限制

无 ICI 标签；CLDN4 蛋白大量 NA；免疫评分来自 RNA；肿瘤/癌旁 RNA 与蛋白方向相反；生存事件少；全表 Spearman 为探索性，推断以主键 FDR 为准。未使用 >2 GB 文件；本切片不做 isoform / 磷酸化 / 甲基化 / CNV / MAF。

### 复现

见上文 English “How to rerun”，或 `notes/grok_cptac_luad/rerun.md`。
