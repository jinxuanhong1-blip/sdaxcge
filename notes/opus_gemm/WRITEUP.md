# Tacstd2 / Cldn4 in public KRAS/LKB1 and NSCLC GEMM ICI datasets

Bilingual write-up (English, then 中文). Parallel slice outputs live only under `notes/opus_gemm/`, `scripts/opus_gemm/`, `results/opus_gemm/`.

**Scope.** Public mouse **KRAS/LKB1 (KL / STK11-null)** and **NSCLC GEMM / GEMM-derived** immune-checkpoint (ICI) transcriptomes. Lewis lung carcinoma (LLC)–only series and SCLC GEMMs (RPM / RP) were screened and **excluded** from the primary statistics. Focus genes: **Tacstd2** (Trop-2) and **Cldn4**.

**Nothing here is a fabricated p-value.** Every number below is in `results/opus_gemm/focus_contrasts.tsv` or `signature_correlations.tsv` and was computed from the processed GEO matrices with Welch or paired *t* tests and Benjamini–Hochberg FDR.

---

## English

### 1. What was asked, and what is actually public

The question is whether Tacstd2 / Cldn4 track ICI treatment or the immune state in mouse lung-cancer models that are **not** the LLC-only story — especially **Kras;Lkb1** and other **NSCLC GEMMs**.

A GEO DataSets screen of 983 mouse series (`results/opus_gemm/geo_screen.tsv`) was parsed from MINiML (design + per-sample characteristics). Processed files only (counts / TPM / FPKM / VST / CPM); no FASTQ. Vendor total **0.502 GB** (`results/opus_gemm/data_manifest.tsv`, MD5s recorded). Ensembl gene IDs for the focus list were resolved against Ensembl REST **release 116** (`results/opus_gemm/gene_ids.tsv`): Tacstd2 = `ENSMUSG00000051397`, Cldn4 = `ENSMUSG00000047501`.

**Included (primary panel)**

| Accession | Model | ICI on the deposited samples? | n used | Why it is in |
|---|---|---|---|---|
| GSE182228 | Lkb1-deficient LUAD, s.c. | anti-PD-1 ± palbociclib vs vehicle | 12 | **Primary KRAS/LKB1 ICI series** (PMID 36871040) |
| GSE114601 | KP GEMM lung nodules | anti-PD-1 ± JQ1 | 8 (n=2/arm) | Autochthonous KP ICI (PMID 30087114) |
| GSE169194 | Kras;Trp53;Msh2null, unfractionated tumour | A2V ± anti-PD-1 | 9 | KP NSCLC + aPD-1 on top of anti-VEGFA/ANGPT2 (PMID 34380768) |
| GSE246922 | KP, CD45− tumour cells | ICB-relapsed vs parental | 15 | Acquired ICB resistance (PMID 38215748) |
| GSE260596 | 344SQ (KrasLA1;p53) | all tumours on anti-PD-1 ± anti-LAIR1 | 7 | GEMM-derived syngeneic; no ICI-free arm |
| GSE197260 | Egfr-del19 s.c. | sequential gefitinib then 4H2 (aPD-1) | 7 (n=1/arm) | Descriptive only |
| GSE274351 | K / KP / KL LCM adenomas | no | 18 | KL vs KP vs K genotype (PMID 39186651) |
| GSE137396 | KL vs KP GEMM nodules | no | 10 | Independent KL vs KP (PMID 34142094) |
| GSE175479 | K vs KL vs Kras;AMPK tumours | no | 12 | LKB1 / AMPK immune-evasion genotypes (PMID 34667030) |
| GSE274352 | KP/KL lines + IFN-β or STING-V154M | no | 12 / 8 | Mechanism from an ICI-sensitisation paper |
| GSE295685 | KL cells ± TNG260 | no (TNG260 is an aPD-1 sensitiser) | 4 (n=2) | STK11-null CoREST inhibition |
| GSE236258 | CMT167 / KLA ± trametinib | no (paper is PD-1 sensitisation) | 15 | KRAS-mutant NSCLC syngeneic lines |
| GSE241978 | CMT167 AhR KO | no | 6 | PD-L1 / IDO paper |
| GSE217405 | Egfr-mutant lines ± osimertinib | no | 24 | Adaptive-immunity EGFR GEMM paper |
| GSE330658 | Egfr-mutant tumours in vivo | anti-PD-L1 **not deposited** | 8 | Baseline / VEGF only |

Full include/exclude list with reasons: `results/opus_gemm/dataset_verification.tsv`.

**Verified exclusions (not analysed as primary)**

- **SCLC GEMM / RP / RPM + ICI** with public processed RNA: GSE309199 (RPM + aPD-1), GSE297817 / GSE297818 (RP + aPD-1+aCTLA-4), GSE208614 (RP-48 + aPD-1), GSE262975 (SCLC + aPD-L1). Out of the KRAS/LKB1–NSCLC scope.
- **LLC-only ICI RNA**: GSE274960 (LL/2 + aPD-1).
- **No usable gene-level matrix**: GSE194166 / GSE232730 (CD45+ scRNA), GSE244452 (DEG table only), GSE298051 (ChIP-seq), GSE301822 (CRISPR sgRNA counts), GSE285606 (10x MTX, IgG-only on the four deposited samples).

### 2. Methods (so the numbers can be regenerated)

1. `scripts/opus_gemm/fetch_data.py` — download the processed GEO files listed in the registry; write URL / bytes / MD5.
2. `scripts/opus_gemm/build_matrices.py` — map columns to GEO samples when the mapping is auditable (explicit description, GSM in the filename, or a documented key regex). Otherwise the file’s own column labels are used and marked `column_labels_only`. Counts → log2(CPM+1); TPM/FPKM/CPM → log2(x+1); VST / depositor log2 used as provided.
3. `scripts/opus_gemm/analyze.py`
   - **Contrasts** are declared in `scripts/opus_gemm/datasets.py` (not inferred ad hoc).
   - Test: Welch two-sample *t* on the log2 matrix; paired *t* when a pairing key is complete (GSE274352 cell line; GSE217405 line × day).
   - **n&lt;2** in either arm → descriptive, no *p*.
   - **n=2** → *p* is still computed and labelled `underpowered`.
   - Effect size: Hedges’ *g*.
   - Multiple testing: Benjamini–Hochberg **within family** (`ici`, `genotype`, `resistance`, `mechanism`) across both genes.
   - Immune scores: mean of within-dataset *z*-scored signature members (CD8, IFNγ, IFN-I, MHC-I, MHC-II, inhibitory, myeloid-suppressive, leukocyte). Spearman ρ of Tacstd2 / Cldn4 vs each score in in-vivo / sorted datasets; BH over all such tests.

Reproduce: `bash scripts/opus_gemm/run_all.sh`.

### 3. Results — ICI contrasts (pre-specified family)

**No Tacstd2 or Cldn4 ICI contrast reaches FDR &lt; 0.05.** Raw *p* values are all ≥ 0.14.

| Dataset | Contrast | Gene | n | Δ log2 (test−ref) | Welch *t* *p* | FDR (ICI family) |
|---|---|---|---|---|---|---|
| GSE182228 KL | aPD-1 vs vehicle | Tacstd2 | 3/3 | −0.041 | 0.79 | 0.79 |
| GSE182228 KL | aPD-1 vs vehicle | Cldn4 | 3/3 | −0.024 | 0.60 | 0.79 |
| GSE182228 KL | aPD-1+palbo vs palbo | Tacstd2 | 3/3 | −0.65 | 0.53 | 0.79 |
| GSE182228 KL | aPD-1+palbo vs palbo | Cldn4 | 3/3 | −0.30 | 0.52 | 0.79 |
| GSE114601 KP | aPD-1 vs vehicle | Tacstd2 | 2/2 | −1.11 | 0.41 | 0.79 |
| GSE114601 KP | aPD-1 vs vehicle | Cldn4 | 2/2 | −1.14 | 0.28 | 0.79 |
| GSE169194 KPM | A2V+aPD-1 vs A2V | Tacstd2 | 3/3 | +0.20 | 0.77 | 0.79 |
| GSE169194 KPM | A2V+aPD-1 vs A2V | Cldn4 | 3/3 | −0.77 | 0.14 | 0.79 |

**Expression context that matters.** In the primary KRAS/LKB1 ICI series (GSE182228), both genes are **near the floor**: median log2(FPKM+1) Tacstd2 = 0.17, Cldn4 = 0.02. Several libraries are exactly 0. You cannot read a treatment effect out of a gene that is not expressed. 344SQ (GSE260596) Tacstd2 is likewise ~0; CMT167 (GSE236258) Tacstd2 is 0 in 14/15 libraries.

GSE197260 (Egfr-del19, n=1/arm) is descriptive: Tacstd2 TPM-scale log2 is 4.01 (vehicle), 1.00 after gefitinib→aPD-1, 3.70 after gefitinib→aPD-1+aVEGFR2. No *p*-value is reported.

### 4. Results — KRAS/LKB1 genotype

Independent GEMM nodules (GSE137396, n=5+5): Tacstd2 is **higher in KL than KP** (Δ log2 = +1.08, Welch *p* = 0.068, FDR = 0.45). Cldn4 Δ = +0.51, *p* = 0.52.

Bulk K vs KL tumours (GSE175479, n=4+4): Cldn4 higher in KL than Kras-only (Δ = +1.68, *p* = 0.063, FDR = 0.45). Tacstd2 is unchanged (Δ = +0.18, *p* = 0.63).

LCM adenomas (GSE274351): no KL vs KP or KL vs K contrast is significant (all FDR ≥ 0.45). Cldn4 is undetectable in the four normal-lung samples (mean 0), so KP-vs-normal Cldn4 raw *p* = 0.031 is a detection contrast, not a genotype effect, and does not survive FDR.

**These are trends, not claims.** With current public n, LKB1 loss is compatible with modestly higher Tacstd2 (nodules) or Cldn4 (bulk tumours) versus KP / Kras-only, but the genotype family FDR does not support a positive call.

### 5. Results — mechanism / resistance (not ICI on the sample)

The **only contrast that survives FDR** in the whole slice is **Cldn4 in KLA cells after trametinib** (GSE236258): mean log2(CPM+1) 2.70 → 9.42, Δ = +6.71, Welch *p* = 1.15×10⁻⁸, FDR (mechanism family) = 3.0×10⁻⁷, n=3/3. The same gene stays near 0 in CMT167 ± trametinib. Tacstd2 is off in both lines. This is a **cell-line-specific, ICI-off-the-sample** finding from a PD-1-sensitisation paper; it is not an in-vivo ICI result.

Other mechanism tests that do **not** pass FDR: IFN-β vs empty Tacstd2 paired *p* = 0.051 (n=6 lines); TNG260 vs DMSO Tacstd2 *p* = 0.040 but n=2 (labelled underpowered, FDR = 0.37); osimertinib vs DMSO Tacstd2 paired *p* = 0.070 (n=12 pairs). ICB-relapsed vs parental KP CD45− cells: no Tacstd2 / Cldn4 shift (all *p* ≥ 0.13).

### 6. Results — immune-signature correlations

One Spearman test survives BH over the 200 in-vivo/sorted tests: **Tacstd2 vs MHC-II** in GSE246922 KP CD45− ICB-relapse series, ρ = 0.85, n=15, *p* = 6.9×10⁻⁵, FDR = 0.014. That matrix is tumour-cell VST, so this is co-variation among malignant cells, not infiltrate.

Several Cldn4–immune correlations have raw *p* &lt; 0.05 (e.g. Cldn4 vs leukocyte score in GSE114601, ρ = −0.86, n=8, *p* = 0.0065; Cldn4 vs MHC-I in GSE137396, ρ = −0.82, n=10, *p* = 0.0038) but **FDR ≥ 0.11**. They are listed in `signature_correlations.tsv` and should not be quoted as discoveries.

In GSE182228 (low Tacstd2/Cldn4), Cldn4 vs MHC-I ρ = 0.69, *p* = 0.012, FDR = 0.21 — same caveat.

### 7. What this does *not* show

- It does **not** show that Tacstd2 or Cldn4 are ICI-response biomarkers in KL or KP mice. No deposited series has a responder / non-responder label plus a gene-level matrix we could test.
- It does **not** show a significant ICI-induced change in either gene.
- SCLC GEMM ICI RNA exists and was verified; it was left out of the statistics because this slice is KRAS/LKB1 + NSCLC only.
- Platforms are heterogeneous (FPKM, CPM, VST, TPM). Contrasts are **within-dataset only**. There is no cross-dataset meta-analysis of the raw values.

### 8. Files

| Path | Content |
|---|---|
| `scripts/opus_gemm/` | Search, fetch, build, analyse, plot |
| `results/opus_gemm/data_manifest.tsv` | URL, bytes, MD5 of every downloaded file |
| `results/opus_gemm/data/matrices/*.expr.tsv.gz` | log2 matrices used for the tests |
| `results/opus_gemm/data/matrices/*.design.tsv` | sample → group, GSM, mapping source |
| `results/opus_gemm/focus_contrasts.tsv` | every contrast, n, Δ, *g*, *p*, FDR |
| `results/opus_gemm/signature_correlations.tsv` | Spearman table |
| `results/opus_gemm/sample_values.tsv` | per-sample Tacstd2, Cldn4, signature scores |
| `results/opus_gemm/dataset_verification.tsv` | include / exclude with reasons |
| `results/opus_gemm/figures/` | forest, GSE182228 strip, genotype, correlation heatmap |

---

## 中文

### 1. 问题与公开数据范围

本切片只分析**公开的 KRAS/LKB1（KL / STK11 缺失）以及 NSCLC GEMM / GEMM 来源同源移植**免疫检查点（ICI）转录组，看 **Tacstd2（Trop-2）** 和 **Cldn4** 是否随 ICI 处理或免疫状态变化。**不包括**单纯 LLC 故事，也**不包括** SCLC GEMM（RPM / RP）。

对 GEO DataSets 中 983 个小鼠系列做了 MINiML 筛查（`geo_screen.tsv`）。只下载处理后的表达矩阵，合计 **0.502 GB**，MD5 见 `data_manifest.tsv`。基因 ID 用 Ensembl REST **release 116** 核对：Tacstd2 = `ENSMUSG00000051397`，Cldn4 = `ENSMUSG00000047501`。

**纳入的主面板**见上文英文表。完整纳入/排除理由：`dataset_verification.tsv`。

已核对但未纳入主统计的数据包括：SCLC+ICI（GSE309199、GSE297817/818、GSE208614、GSE262975）、仅 LLC 的 GSE274960、以及没有可用基因水平矩阵的 scRNA / ChIP / CRISPR 计数系列。

### 2. 统计（可复现，未编造）

对比在 `datasets.py` 中预先写死。检验：对数矩阵上的 Welch *t*；配对设计用配对 *t*。任一组 n&lt;2 只做描述、不算 *p*。n=2 仍给出 *p* 并标明 underpowered。效应量 Hedges’ *g*。多重检验：按家族（`ici` / `genotype` / `resistance` / `mechanism`）对两个基因一起做 BH FDR。免疫签名为数据集内 *z* 分数均值，再与 Tacstd2/Cldn4 做 Spearman，全部相关一起 BH。

复现：`bash scripts/opus_gemm/run_all.sh`。

### 3. ICI 对比结果

**没有任何 Tacstd2 / Cldn4 的 ICI 对比达到 FDR &lt; 0.05。** 原始 *p* 全部 ≥ 0.14。具体数字见英文表与 `focus_contrasts.tsv`。

关键背景：主系列 **GSE182228（LKB1 缺失 LUAD + anti-PD-1）** 中两个基因几乎不表达（Tacstd2 中位 log2(FPKM+1)=0.17，Cldn4=0.02）。344SQ 的 Tacstd2、CMT167 的 Tacstd2 同样接近 0。在未表达的基因上谈“ICI 调控”没有依据。

GSE197260（Egfr 19 缺失，每臂 n=1）只做描述：gefitinib 后接 aPD-1 的 Tacstd2 低于 vehicle，但没有推断性 *p* 值。

### 4. KRAS/LKB1 基因型

GSE137396（KL vs KP 结节，各 n=5）：Tacstd2 在 KL 更高（Δ=+1.08，*p*=0.068，FDR=0.45）。  
GSE175479（KL vs 仅 Kras，各 n=4）：Cldn4 在 KL 更高（Δ=+1.68，*p*=0.063，FDR=0.45）。  
GSE274351 LCM 腺瘤：KL vs KP / KL vs K 均不显著。

这是**趋势，不是阳性结论**。现有公开样本量不足以在 FDR 控制下声称 LKB1 缺失上调 Tacstd2 或 Cldn4。

### 5. 机制 / 耐药（样本上没有 ICI）

全切片**唯一通过 FDR 的对比**是 GSE236258 中 **KLA 细胞 trametinib 处理后 Cldn4**：2.70 → 9.42（log2 CPM+1），Δ=+6.71，*p*=1.15×10⁻⁸，FDR=3.0×10⁻⁷。CMT167 上 Cldn4 接近 0，Tacstd2 两条细胞系都关着。这是细胞系、且样本上无 ICI 的结果，不能当成体内 ICI 结论。

IFN-β、TNG260、osimertinib 对 Tacstd2 的原始 *p* 在 0.04–0.07，FDR 均 ≥ 0.37。KP ICB 复发相对亲本：两个基因都无显著变化。

### 6. 与免疫签名的相关

200 次 Spearman 中，BH 后只剩 **一条**：GSE246922 KP CD45− ICB 复发系列里 Tacstd2 与 MHC-II，ρ=0.85，n=15，*p*=6.9×10⁻⁵，FDR=0.014。这是肿瘤细胞 VST，不是浸润评分。

其余 raw *p*&lt;0.05 的 Cldn4–免疫相关（如 GSE114601 白细胞签名 ρ=−0.86）FDR 均 ≥ 0.11，**不能当发现引用**。

### 7. 明确做不到的事

- 没有带应答/无应答标签的公开 KL/KP 基因水平矩阵，因此**不能**声称 Tacstd2/Cldn4 是小鼠 ICI 疗效标志物。
- **不能**声称 ICI 显著改变这两个基因。
- SCLC GEMM+ICI 的 RNA 是公开的，已核实，但不在本切片统计里。
- 平台混杂，只做数据集内对比，没有跨数据集合并原始值。

### 8. 文件

与英文第 8 节同一套路径。图在 `results/opus_gemm/figures/`。
