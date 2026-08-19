# Additive Stereo-seq LUAD, CLDN4-only

Open BGI Stereo-seq / Stereo-XCR-seq lung adenocarcinoma exists and was downloaded: **[GSE328481](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328481)** (11 LUAD cell-bin H5ADs, public 2026-05-29). No Visium. No private 8-KL. CLDN4 only.

---

## Hunt (GEO, SRA, CNGB, STOmics)

| Accession | Where | What | Access |
|---|---|---|---|
| **GSE328481** | GEO | Stereo-XCR-seq, 11 LUAD tumors, cell-bin H5AD (~6.1 GB) | **Public — analyzed** |
| HRA018548 | GSA-Human / PRJCA047766 | Same study, Stereo-XCR-seq FASTQ | Controlled (DAC HDAC007001) |
| HRA013757 | GSA-Human / PRJCA047766 | Same study, Stereo-seq FASTQ | Controlled |
| HRA009729 | GSA-Human | Earlier Stereo-XCR paper (NSCLC + ccRCC + IBD) FASTQ | Controlled |
| OMIX014896 | NGDC OMIX | Paper twin of GSE328481 processed matrices | Duplicate of GEO; not re-pulled |
| STT0000123 | STOmicsDB | 4 sections (kidney cancer, **lung cancer**, IBD, CRC), 2.83 GB | Archive page lists files; HTML table empty without login |
| **HRA004240** | GSA-Human | Stereo-seq of **1 lung adenosquamous carcinoma** + WES | **Controlled** (DAC YCH / Shanghai Chest) |
| **CNP0005129** | CNGBdb | LUAD STAS “single-cell + spatial”, 14 samples, 534.97 GB | **Controlled — Apply for data**. Project page does **not** name Stereo-seq; published STAS papers from other groups are Visium/GeoMx |
| STDS0000164 | STOmicsDB curated | NSCLC + brain mets | GeoMx (GSE200563) — not Stereo-seq; not used |
| STT0000179 / STDS0000062 / STT0000184 | STOmicsDB | Mouse lung Stereo-seq | Not human cancer |
| GSE294759 | GEO | Stereo-seq lung **organoids** | Not LUAD/NSCLC tumor |
| GSE260779 | GEO | Mouse lung + human P63+ progenitors | Not lung cancer |
| GSE285990 | GEO | Stereo-seq liver metastases | Not lung primary |
| SRA PRJNA1455526 | SRA | BioProject of GSE328481 | GEO: “Raw data not provided”; raw is GSA-Human controlled |

Full decision log: `results/stereo_seq_luad/tables/hunt_catalog.tsv`.

---

## Methods (short)

Public cell-bin AnnData from GEO FTP (`uns/bin_type=cell_bins`, 500 nm/pixel). QC: ≥50 genes / ≥20 UMI (cell-bin); ≥100 genes / ≥50 UMI (bin50 sums). log1p(CP10K). Epithelial = author `anno` in `{cancer_cell, epi}` (bin50 if ≥50% of cells). CD8A+ = raw count > 0. Nearest-CD8 reference **excludes** epithelial units. CLDN4-high epithelial = CLDN4>0 (top half of detected if n≥40); CLDN4-low = epithelial CLDN4==0. KRT8 residual = OLS residual of log1p(CLDN4) ~ log1p(KRT8). Sample (n=11) is the inference unit. Playbook: `methods/stereo_seq_luad/playbook.md`.

---

## Results — GSE328481, 11 LUAD sections

### 1. Same-unit Spearman CLDN4 vs CD8A

| Unit | median ρ | + / − | Wilcoxon p |
|---|---:|---:|---:|
| cell-bin, all units | **+0.0058** | 9 / 2 | **0.019** |
| cell-bin, epithelial only | +0.0040 | 9 / 2 | 0.083 |
| bin50, all units | +0.0080 | 9 / 2 | 0.10 |
| bin50, epithelial only | +0.0086 | 8 / 3 | 0.042 |

All-unit ρ is a few thousandths. Several per-section p-values are small because each section has 3×10⁴–1.9×10⁵ units; the **effect size is near zero**, not a CD8-exclusion pattern.

### 2. KRT8 residual

| Unit | median residual ρ vs CD8A | Wilcoxon p |
|---|---:|---:|
| cell-bin | **+0.0015** | **0.58** |
| cell-bin, epithelial | +0.0033 | 0.76 |
| bin50 | +0.0052 | 0.37 |
| bin50, epithelial | +0.0048 | 0.28 |

After removing the KRT8 (epithelial-content) axis, CLDN4 is **not** associated with CD8A.

### 3. Nearest CD8A+ distance from CLDN4-high epithelial units

Reference set = CD8A+ **non-epithelial** units. ∆ = median distance(high) − median distance(low); positive = high farther from CD8.

| Unit | median ∆µm | high closer / farther | Wilcoxon p |
|---|---:|---:|---:|
| cell-bin | **−1.93** | 7 / 4 | **0.46** |
| bin50 | **0.00** | 5 / 3 (3 ties) | 0.84 |

Typical nearest CD8A+ is ~100–250 µm. The high-vs-low gap is 1–20 µm in most sections and is not signed consistently. Outlier: LUAD_P8 cell-bin ∆ = **+26.8 µm** (MW p = 1.6×10⁻¹¹); bin50 ∆ = **+74.3 µm**. LUAD_P11 is the opposite (cell-bin ∆ = −12.3 µm, MW p = 2.1×10⁻⁵).

### 4. Detection exclusivity (not a neighborhood claim)

CLDN4 and CD8A almost never share a cell-bin (expected for an epithelial vs T-cell gene). Among units with CLDN4>0 **or** CD8A>0, median Spearman = **−0.62** (11/11 negative). That is zero-inflation / lineage exclusivity. Co-detection rates are near the background CD8A+ rate, e.g. LUAD_P1: CD8A+ in 0.71% of CLDN4+ bins vs 0.68% of CLDN4− bins.

---

## What these numbers say

On **open** Stereo-seq LUAD (GSE328481), same-bin CLDN4–CD8A correlation is ~0; the KRT8 residual is null; CLDN4-high epithelial bins are **not** systematically farther from CD8A+ cells. This slice does not use Visium, CosMx, GeoMx, or any private 8-KL object.

---

## 中文

公开的人肺腺癌 Stereo-seq 存在：**GSE328481**（11 例 LUAD cell-bin H5AD，已下载）。受控：GSA-Human **HRA018548 / HRA013757 / HRA009729**（同一/姊妹 Stereo-XCR 原始 FASTQ）、**HRA004240**（1 例肺腺鳞癌 Stereo-seq）、CNGBdb **CNP0005129**（LUAD STAS，申请下载；页面未写明是否 Stereo-seq）。STOmics **STT0000123** 含肺癌切片但网页未暴露文件。未改做 Visium。未用私有 8-KL。只做 CLDN4。

**同 bin Spearman** cell-bin 中位 ρ = +0.0058（9/11 为正，p = 0.019），效应接近 0。**KRT8 残差** 中位 ρ = +0.0015（p = 0.58）。**最近 CD8 距离** CLDN4 高 vs 低上皮中位差 −1.93 µm（p = 0.46），无稳定更远/更近。非零并集 ρ = −0.62 是检测互斥，不是邻域排除。

表与图：`results/stereo_seq_luad/`。复现：`methods/stereo_seq_luad/playbook.md`。
