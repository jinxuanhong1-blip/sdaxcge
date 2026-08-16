# TACSTD2 / CLDN4 vs SCLC immune subtype and ICI outcome

**Slice outputs:** `notes/fable_sclc/`, `scripts/fable_sclc/`, `results/fable_sclc/`  
**Accessions verified:** 2026-08-16 (see `DATASETS.md`)  
**Processed download:** 1.49 GB (scratch `/tmp/sclc_data`; not in the repo)  
**No fabricated statistics.** Every number below is written by `scripts/fable_sclc/` from public processed matrices. Full tables: `results/fable_sclc/tables/`.

---

## English

### Question

Does **TACSTD2** (TROP2) or **CLDN4** associate with the inflamed SCLC transcriptional class (**SCLC-I**) or with immune-checkpoint inhibitor (ICI) response in public human datasets (bulk / scRNA / spatial)?

### What is public vs controlled

Open processed expression with ICI outcome is scarce. The only open expression + per-patient ICI outcome we could verify are two GeoMx DSP series from first-line chemo-immunotherapy trials:

| Cohort | Accession | Drug | n patients / ROIs | Genes |
|---|---|---|---|---|
| CANTABRICO | [GSE261345](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261345) | durvalumab + platinum/etoposide | 26 / 121 | CTA (~1.8k); **TACSTD2 yes, CLDN4 no** |
| IMfirst | [GSE261348](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261348) | atezolizumab + platinum/etoposide | 32 / 175 | same panel |

The randomized trial that defined SCLC-I benefit from atezolizumab — **IMpower133**, n=271 RNA-seq — is **controlled**: EGA [EGAS50000000138](https://ega-archive.org/studies/EGAS50000000138) (DAC `devsci-dac-d@gene.com`). CheckMate-032 RNA-seq (n=286, nivolumab ± ipilimumab) has **no public accession**. Full controlled list: `DATASETS.md`.

Subtype / immune-context analyses therefore use two pre-ICI bulk RNA-seq cohorts (George 2015; GSE60052) plus the Chan 2021 HTAN scRNA atlas (cellxgene collection `62e8f058-9c37-48bc-9200-e767f318a8ec`; raw FASTQ is dbGaP `phs002371`).

### Methods (reproducible)

```
cd scripts/fable_sclc && pip install -r requirements.txt && ./run_all.sh
```

- **Bulk subtypes (Gay 2021-style NMF):** top 1,250 variable genes, NMF k=4; factors labelled by Spearman correlation of factor weights with ASCL1 / NEUROD1 / POU2F3; leftover factor = SCLC-I. A bimodal marker-threshold rule is stored as a sensitivity annotation (`subtype_marker_rule`).
- **Immune scores:** T-effector = mean z of CD8A, GZMA, GZMB, PRF1, IFNG, CXCL9; cytolytic = GZMA+PRF1; antigen presentation = HLA-A/B/C, B2M, TAP1/2, NLRC5.
- **GeoMx:** log2(normalized count + 1); patient-level = median across ROIs; PFS/OS reconstructed from GEO dosing / progression / follow-up dates; Cox PH of per-cohort z-scored TACSTD2 with a cohort covariate.
- **scRNA:** log1p(CP10K) from raw counts (per-cell totals streamed from `.raw.X`); donor pseudobulk = mean in epithelial cells (≥50 cells); subtype via the marker rule (n=20 donors is too small for NMF).
- Tests are two-sided. BH FDR is applied **within each table**.

NMF SCLC-I is immune-hot as expected (sanity check, not a discovery): George 2015 Teff median 0.67 vs −0.31, Mann–Whitney U=923, p=1.14×10⁻⁵ (n_I=17 / n_other=64); GSE60052 0.21 vs −0.16, U=496, p=0.026 (n_I=10 / n_other=69). NMF vs marker-rule agreement: 67% (George), 49% (GSE60052).

### Results

#### 1. Bulk: TACSTD2 is higher in SCLC-I in GSE60052; CLDN4 is not

NMF counts — George 2015 (n=81): A 36 / N 18 / I 17 / P 10. GSE60052 tumors (n=79): A 43 / N 17 / I 10 / P 9.

| Cohort | Gene | median A / N / P / I | Kruskal–Wallis H, p (BH) | I vs rest MWU p (BH) |
|---|---|---|---|---|
| George 2015 | TACSTD2 | 2.34 / 1.84 / 3.33 / **3.34** | 5.70, 0.127 (0.127) | 0.127 (0.127) |
| George 2015 | CLDN4 | 6.49 / 6.34 / 5.36 / 5.81 | 8.99, 0.029 (0.059) | 0.124 (0.127) |
| GSE60052 | TACSTD2 | 1.72 / 4.42 / 3.04 / **5.07** | 25.92, **9.90×10⁻⁶ (3.96×10⁻⁵)** | **4.04×10⁻⁵ (1.61×10⁻⁴)** |
| GSE60052 | CLDN4 | 6.27 / 6.52 / 5.64 / 4.50 | 7.70, 0.053 (0.070) | 0.094 (0.127) |

Tumor vs adjacent-normal in GSE60052 (n=79 vs 7): TACSTD2 median 2.48 vs 6.28, U=61, p=6.40×10⁻⁴ (BH 1.28×10⁻³) — tumors **lower**. CLDN4 6.26 vs 6.27, p=0.817.

#### 2. Bulk: TACSTD2 tracks immune / myeloid genes in GSE60052; CLDN4 does not

Spearman (selected; full table `bulk_immune_correlations.csv`):

| Cohort | Pair | ρ | p (BH) |
|---|---|---|---|
| GSE60052 | TACSTD2–Teff | 0.407 | 2.01×10⁻⁴ (0.0024) |
| GSE60052 | TACSTD2–cytolytic | 0.486 | 5.57×10⁻⁶ (1.00×10⁻⁴) |
| GSE60052 | TACSTD2–APM | 0.330 | 0.0029 (0.017) |
| GSE60052 | TACSTD2–CD8A | 0.318 | 0.0043 (0.020) |
| GSE60052 | TACSTD2–CD274 | 0.396 | 3.06×10⁻⁴ (0.0028) |
| GSE60052 | TACSTD2–CD163 | 0.487 | 5.45×10⁻⁶ (1.00×10⁻⁴) |
| George 2015 | TACSTD2–Teff | 0.108 | 0.338 (0.487) |
| George 2015 | TACSTD2–CD274 | 0.274 | 0.013 (0.048) |
| George 2015 | TACSTD2–STAT1 | 0.287 | 0.0095 (0.038) |
| George 2015 | TACSTD2–CD163 | 0.323 | 0.0033 (0.017) |
| GSE60052 | CLDN4–Teff | 0.036 | 0.751 (0.795) |
| GSE60052 | CLDN4–STAT1 | −0.351 | 0.0015 (0.011) |
| George 2015 | CLDN4–Teff | −0.177 | 0.114 (0.228) |

**Read-out:** TACSTD2–immune is **cohort-dependent** (clear in GSE60052; only myeloid / IFN / PD-L1 genes in George). CLDN4 is not an immune-hot marker; the only BH-significant CLDN4 correlation is an **inverse** association with STAT1 in GSE60052.

#### 3. Spatial ICI (GeoMx): TACSTD2 co-varies with T-effector genes inside tumors, but not with RECIST or survival

ROI-level Spearman TACSTD2 vs T-effector score: CANTABRICO ρ=0.687 (n=121, p=3.51×10⁻¹⁸); IMfirst ρ=0.649 (n=175, p=2.86×10⁻²²). Same pattern for CD8A, GZMB, PRF1, IFNG, CD274 (ρ 0.61–0.76, all BH q < 0.001). CXCL9 and CD68 are not significant. CLDN4 is **absent from the CTA panel**.

Patient-level TACSTD2 (median of ROIs) vs best RECIST (CR/PR vs SD/PD):

| Scope | n CR/PR | n SD/PD | median resp / non-resp | MWU p (BH) |
|---|---|---|---|---|
| pooled | 38 | 18 | 7.45 / 7.38 | 0.629 (0.856) |
| CANTABRICO | 17 | 9 | 7.25 / 7.36 | 0.258 (0.773) |
| IMfirst | 21 | 9 | 7.63 / 7.58 | 0.856 (0.856) |

Cox PH, n=58, TACSTD2 z-score + cohort flag: **PFS** HR per SD = 0.937 (95% CI 0.680–1.291), p=0.690, 50 events; **OS** HR = 1.055 (0.777–1.432), p=0.732, 44 events. Kaplan–Meier PFS by cohort-wise median split: 29 high / 29 low, log-rank statistic 0.817, **p=0.366**.

These two single-arm chemo-IO series **cannot** isolate an ICI-specific effect (no chemotherapy-only arm). They also cannot test CLDN4.

#### 4. scRNA (Chan 2021): CLDN4 is a tumor-cell gene; TACSTD2 is rare; donor tests are underpowered

SCLC samples: 77,143 cells / 20 donors. In epithelial (tumor) cells (n=55,876): TACSTD2 detected in **2.17%** of cells (mean log1p CP10K 0.023); CLDN4 in **61.3%** (mean 0.968). Lymphoid / myeloid / mesenchymal positivity for CLDN4 is 2–4%.

Donor epithelial pseudobulk (n=20; A 10 / I-like 5 / N 4 / P 1): TACSTD2 median 0.014 (A) vs 0.122 (I-like); CLDN4 0.980 vs 0.944. Kruskal–Wallis across subtypes: TACSTD2 H=3.12, p=0.210; CLDN4 H=3.61, p=0.164 (BH 0.224). ICI-exposed (n=9) vs naive (n=11): TACSTD2 medians 0.021 vs 0.011, U=70, p=0.129; CLDN4 0.690 vs 1.117, U=33, p=0.224. **No donor-level test is significant.** Treatment is mixed (not a trial) and ICI-exposed donors are not paired pre/post.

### Interpretation (what we can and cannot claim)

1. **TACSTD2 is compatible with an SCLC-I / immune-infiltrated state** in GSE60052 (subtype + Teff/CYT/CD8A/CD274) and in both GeoMx ICI trials (ROI-level T-cell / PD-L1 genes). George 2015 does not replicate the Teff correlation.
2. **That does not translate into ICI benefit in the two open GeoMx trials** (response, PFS, OS all null). The biologically relevant test — IMpower133 atezolizumab vs chemo, with official SCLC-I labels — is controlled (EGAS50000000138).
3. **CLDN4 is a common SCLC tumor-cell gene**, not an immune-hot marker, and is missing from the only open ICI expression matrices.
4. TACSTD2 is **lower** in SCLC bulk tumors than in the seven GSE60052 normal lungs, and is **sparse** in Chan 2021 tumor cells. Any TROP2-directed therapy hypothesis in SCLC needs protein / spatial confirmation; this slice is RNA only.

### Limitations

NMF subtypes are a local re-implementation, not the published Gay labels. GeoMx is a targeted panel on mixed “Full ROI” segments (tumor + stroma). Survival times are reconstructed from GEO dates. scRNA n=20 donors. No randomized ICI arm is public.

---

## 中文

### 问题

在公开的人源 SCLC 组学数据（bulk / 单细胞 / 空间）中，**TACSTD2（TROP2）** 与 **CLDN4** 是否与炎症亚型 **SCLC-I** 或免疫检查点抑制剂（ICI）疗效相关？

### 公开与受控数据

带 ICI 疗效的开放表达数据极少。我们核实到的、同时具备表达矩阵与逐例患者疗效的开放数据，只有两个一线化疗联合免疫的 GeoMx DSP 队列：

| 队列 | 登录号 | 药物 | 患者 / ROI | 基因 |
|---|---|---|---|---|
| CANTABRICO | GSE261345 | 度伐利尤单抗 + 铂类/依托泊苷 | 26 / 121 | CTA；**有 TACSTD2，无 CLDN4** |
| IMfirst | GSE261348 | 阿替利珠单抗 + 铂类/依托泊苷 | 32 / 175 | 同上 |

定义 SCLC-I 从阿替利珠单抗获益的随机试验 **IMpower133**（n=271 RNA-seq）为**受控数据**：EGA EGAS50000000138（DAC：`devsci-dac-d@gene.com`）。CheckMate-032 RNA-seq（n=286）**没有公开登录号**。完整受控清单见 `DATASETS.md`。

亚型 / 免疫微环境分析因此使用两个 ICI 前 bulk 队列（George 2015；GSE60052）以及 Chan 2021 HTAN 单细胞图谱（cellxgene `62e8f058-9c37-48bc-9200-e767f318a8ec`；原始 FASTQ 为 dbGaP `phs002371`）。

### 方法

```
cd scripts/fable_sclc && pip install -r requirements.txt && ./run_all.sh
```

- **Bulk 亚型：** 按 Gay 2021 思路，取变异最大的 1,250 个基因做 NMF（k=4）；因子按与 ASCL1 / NEUROD1 / POU2F3 的 Spearman 相关标注，剩余因子为 SCLC-I。另存双峰阈值 marker 规则作为敏感性注释。
- **免疫评分：** T-effector = CD8A、GZMA、GZMB、PRF1、IFNG、CXCL9 的 z 均值；细胞毒 = GZMA+PRF1；抗原呈递 = HLA-A/B/C、B2M、TAP1/2、NLRC5。
- **GeoMx：** log2(归一化计数+1)；患者水平取 ROI 中位数；PFS/OS 由 GEO 给药 / 进展 / 随访日期重建；Cox 使用队列内 z 标准化的 TACSTD2，并加入队列协变量。
- **单细胞：** 由 raw counts 计算 log1p(CP10K)；供者伪 bulk = 上皮细胞均值（≥50 个细胞）；亚型用 marker 规则（n=20 不宜做 NMF）。
- 检验均为双侧；BH FDR 在**每张表内部**校正。

NMF 得到的 SCLC-I 确为免疫热（合理性检查，非新发现）：George 2015 Teff 中位数 0.67 vs −0.31，U=923，p=1.14×10⁻⁵；GSE60052 为 0.21 vs −0.16，U=496，p=0.026。NMF 与 marker 规则一致率：George 67%，GSE60052 49%。

### 结果

#### 1. Bulk：GSE60052 中 TACSTD2 在 SCLC-I 更高；CLDN4 不是

NMF 计数 — George 2015（n=81）：A 36 / N 18 / I 17 / P 10。GSE60052 肿瘤（n=79）：A 43 / N 17 / I 10 / P 9。

| 队列 | 基因 | 中位数 A / N / P / I | Kruskal–Wallis H, p（BH） | I vs 其余 MWU p（BH） |
|---|---|---|---|---|
| George 2015 | TACSTD2 | 2.34 / 1.84 / 3.33 / **3.34** | 5.70, 0.127（0.127） | 0.127（0.127） |
| George 2015 | CLDN4 | 6.49 / 6.34 / 5.36 / 5.81 | 8.99, 0.029（0.059） | 0.124（0.127） |
| GSE60052 | TACSTD2 | 1.72 / 4.42 / 3.04 / **5.07** | 25.92, **9.90×10⁻⁶（3.96×10⁻⁵）** | **4.04×10⁻⁵（1.61×10⁻⁴）** |
| GSE60052 | CLDN4 | 6.27 / 6.52 / 5.64 / 4.50 | 7.70, 0.053（0.070） | 0.094（0.127） |

GSE60052 肿瘤 vs 癌旁（79 vs 7）：TACSTD2 中位数 2.48 vs 6.28，U=61，p=6.40×10⁻⁴（BH 1.28×10⁻³）——肿瘤**更低**。CLDN4 6.26 vs 6.27，p=0.817。

#### 2. Bulk：TACSTD2 在 GSE60052 与免疫 / 髓系基因相关；CLDN4 不相关

Spearman（节选；全表见 `bulk_immune_correlations.csv`）：

| 队列 | 配对 | ρ | p（BH） |
|---|---|---|---|
| GSE60052 | TACSTD2–Teff | 0.407 | 2.01×10⁻⁴（0.0024） |
| GSE60052 | TACSTD2–细胞毒 | 0.486 | 5.57×10⁻⁶（1.00×10⁻⁴） |
| GSE60052 | TACSTD2–抗原呈递 | 0.330 | 0.0029（0.017） |
| GSE60052 | TACSTD2–CD8A | 0.318 | 0.0043（0.020） |
| GSE60052 | TACSTD2–CD274 | 0.396 | 3.06×10⁻⁴（0.0028） |
| GSE60052 | TACSTD2–CD163 | 0.487 | 5.45×10⁻⁶（1.00×10⁻⁴） |
| George 2015 | TACSTD2–Teff | 0.108 | 0.338（0.487） |
| George 2015 | TACSTD2–CD274 | 0.274 | 0.013（0.048） |
| George 2015 | TACSTD2–STAT1 | 0.287 | 0.0095（0.038） |
| George 2015 | TACSTD2–CD163 | 0.323 | 0.0033（0.017） |
| GSE60052 | CLDN4–Teff | 0.036 | 0.751（0.795） |
| GSE60052 | CLDN4–STAT1 | −0.351 | 0.0015（0.011） |
| George 2015 | CLDN4–Teff | −0.177 | 0.114（0.228） |

**解读：** TACSTD2–免疫相关**依赖队列**（GSE60052 清楚；George 仅髓系 / IFN / PD-L1）。CLDN4 不是免疫热标志；唯一 BH 显著的 CLDN4 相关是 GSE60052 中与 STAT1 的**负相关**。

#### 3. 空间 ICI（GeoMx）：TACSTD2 与瘤内 T 效应基因共变，但与 RECIST / 生存无关

ROI 水平 TACSTD2 与 T-effector：CANTABRICO ρ=0.687（n=121，p=3.51×10⁻¹⁸）；IMfirst ρ=0.649（n=175，p=2.86×10⁻²²）。CD8A、GZMB、PRF1、IFNG、CD274 同样 ρ=0.61–0.76（BH q < 0.001）。CXCL9、CD68 不显著。CTA 面板**没有 CLDN4**。

患者水平 TACSTD2（ROI 中位数）与最佳 RECIST（CR/PR vs SD/PD）：

| 范围 | n CR/PR | n SD/PD | 中位数 有效 / 无效 | MWU p（BH） |
|---|---|---|---|---|
| 合并 | 38 | 18 | 7.45 / 7.38 | 0.629（0.856） |
| CANTABRICO | 17 | 9 | 7.25 / 7.36 | 0.258（0.773） |
| IMfirst | 21 | 9 | 7.63 / 7.58 | 0.856（0.856） |

Cox（n=58，TACSTD2 z + 队列）：**PFS** 每 SD 的 HR=0.937（95% CI 0.680–1.291），p=0.690，50 个事件；**OS** HR=1.055（0.777–1.432），p=0.732，44 个事件。按队列内中位数分组的 PFS KM：高/低各 29 例，log-rank 统计量 0.817，**p=0.366**。

这两个单臂化疗联合免疫队列**无法**分离 ICI 特异性效应（无单纯化疗臂），也无法检测 CLDN4。

#### 4. 单细胞（Chan 2021）：CLDN4 是肿瘤细胞基因；TACSTD2 稀疏；供者水平检验效力不足

SCLC：77,143 个细胞 / 20 名供者。上皮（肿瘤）细胞（n=55,876）：TACSTD2 检出率 **2.17%**（均值 0.023）；CLDN4 **61.3%**（均值 0.968）。淋巴 / 髓系 / 间质中 CLDN4 阳性率仅 2–4%。

供者上皮伪 bulk（n=20；A 10 / I-like 5 / N 4 / P 1）：TACSTD2 中位数 A 0.014 vs I-like 0.122；CLDN4 0.980 vs 0.944。亚型 Kruskal–Wallis：TACSTD2 H=3.12，p=0.210；CLDN4 H=3.61，p=0.164（BH 0.224）。ICI 暴露（n=9）vs 初治（n=11）：TACSTD2 0.021 vs 0.011，U=70，p=0.129；CLDN4 0.690 vs 1.117，U=33，p=0.224。**供者水平检验均不显著。** 治疗背景混杂，并非临床试验，也无配对治疗前后样本。

### 结论（能说与不能说）

1. **TACSTD2 与 SCLC-I / 免疫浸润状态相容**：GSE60052（亚型 + Teff/细胞毒/CD8A/CD274）以及两个 GeoMx ICI 队列（ROI 水平 T 细胞 / PD-L1 基因）。George 2015 未重复 Teff 相关。
2. **这并未转化为两个开放 GeoMx 试验中的 ICI 获益**（疗效、PFS、OS 均为阴性）。真正相关的检验——IMpower133 阿替利珠单抗 vs 化疗、带官方 SCLC-I 标签——是受控数据（EGAS50000000138）。
3. **CLDN4 是常见的 SCLC 肿瘤细胞基因**，不是免疫热标志，且缺失于目前唯一开放的 ICI 表达矩阵。
4. TACSTD2 在 GSE60052 bulk 肿瘤中**低于** 7 例正常肺，在 Chan 2021 肿瘤细胞中**稀疏**。任何 SCLC TROP2 靶向假说需要蛋白 / 空间验证；本切片仅为 RNA。

### 局限

NMF 亚型是本地重现，不是 Gay 原文标签。GeoMx 为混合 “Full ROI” 的靶向面板。生存时间由 GEO 日期重建。单细胞仅 20 名供者。没有公开的随机 ICI 对照臂。
