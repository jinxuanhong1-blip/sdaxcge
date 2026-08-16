# Leftover 2026-only GEO lung ICI series: TACSTD2 / CLDN4

> Scope: **PDAT year 2026 only**. Prior 2024/2025 analyses (GSE261345,
> GSE261348, GSE233203) are out of scope here.
> Outputs: `results/w200/GEO_2026/leftover_2026_*` and
> `scripts/fable_geo_2026/08_leftover_2026.py`, `09_analyze_leftover_2026.py`.
> Search refreshed **2026-08-16 UTC**. English first, 中文 second.

---

## English

### What “leftover 2026” means
The 2024–2026 search returned **44 series with PDAT 2026**. Live series-matrix
and FTP scans are in `leftover_2026_catalog.tsv` and `leftover_2026_supp.tsv`.
Most cannot test tumor TACSTD2/CLDN4 vs ICI response (blood, cell-line,
methylation/miRNA, n=1, missing labels, or files > 2 GB with no smaller
processed matrix).

Two leftover tumor series **do** have open processed expression **and**
per-sample ICI labels that were not analyzed in the first pass:

| GSE | Why it was leftover | What was used |
|-----|---------------------|---------------|
| **GSE329813** | First pass said “no per-ROI MPR label” because `!Sample_characteristics` only has `tissue`/`batch`. **MPR/NMPR is in `!Sample_title`.** | GeoMx normalized CSV (1.6 MB); 22 patients |
| **GSE292299** | No single series matrix; 3.1 GB `RAW.tar` exceeds the 2 GB cap. **Per-sample Visium H5 files are each < 2 GB** on the GSM FTP. | 16 NSCLC H5 + `sample_metadata.csv` |

### Results

**GSE329813 — NSCLC neoadjuvant pembrolizumab + platinum, post-surgery GeoMx DSP**
- Labels parsed from titles (`ROI k, Patient n, Primary tumor bed|Lymph node, MPR|NMPR`). Not invented.
- TACSTD2 is on the processed 1,812-gene matrix. **CLDN4 is absent.**
- Primary tumor bed, patient-mean TACSTD2: MPR n=11 median 3.94 vs NMPR n=11 median 4.80; Mann–Whitney p = **0.002**; Cliff's δ = **−0.79** (lower in MPR).
- Draining lymph node: p = 0.83, δ = 0.08 (null).
- **Honest reading:** these are **post-treatment residual** ROIs. Lower TACSTD2 in MPR tumor beds is the expected leftover-epithelial / residual-tumor signal after a major pathologic response. It is **not** a pretreatment predictive biomarker. The null LN result is consistent with that composition story.

**GSE292299 — pretreatment Visium, 16 NSCLC, ICI (pembro / nivo / ipi-nivo / chemo-IO)**
- `Tx_Response` in the deposited metadata: 12 R vs 4 NR. Whole-slide log2(CPM+1).
- TACSTD2: R 8.68 vs NR 7.37; p = 0.26; δ = 0.42 (ns).
- CLDN4: R 8.60 vs NR 6.56; p = **0.042**; δ = 0.71. Spearman vs PFS months ρ = 0.55 (p = 0.027; censoring ignored).
- TACSTD2–CLDN4 concordance ρ = 0.52 (p = 0.040).
- **Honest reading:** the CLDN4 contrast is driven by **two of four** NR slides (pleura P14, lung P15) with very low CLDN4; the other two NR overlap the R range. Eight of twelve responders are **brain metastases**. Whole-slide Visium mixes tumor and stroma. n_NR = 4 and no multiplicity correction. This is **exploratory, not confirmatory**.

### 2026 series checked and not used for this question
- Blood / fluid with or without RECIST in the paper but not usable for epithelial genes: GSE305086, GSE306542, GSE310370, GSE315510, GSE317309 (PBMC CITE-seq; 7.6 GB tar), GSE295601 / GSE267026 / GSE267027 / GSE308745 / GSE312336 (blood + mouse).
- Tumor RNA without ICI response labels: GSE311200 (pretreatment NSCLC RNA-seq; titles are EGFR-mutant IDs only), GSE292700 (BRCA TME, 2 patients), GSE299684 (relapsed SCLC, no ICI label), GSE337519 (n=1, no response), GSE318850 (n=1 series record; 8.1 GB Seurat object).
- EGFR-TKI resistance, not ICI: GSE253718.
- Spatial TLS / mixed cancers, no ICI label: GSE322553, GSE327192 (13 GB).
- Cell-line / mouse / chromatin: GSE282437, GSE287798/799/930, GSE288117, GSE293914, GSE298300, GSE298789, GSE303043, GSE304741, GSE309452/453, GSE313827, GSE319755, GSE319956, GSE320009, GSE320531, GSE327377, GSE328294, GSE330687.
- GSE292421 (2026 PDAT): still only TME immunophenotype + unmappable FPKM IDs; not ICI response.
- GSE307029: methylation array, not TACSTD2/CLDN4 expression.

### Reproduce
```bash
python3 scripts/fable_geo_2026/08_leftover_2026.py
python3 scripts/fable_geo_2026/09_analyze_leftover_2026.py
```
Source SHA-256: `leftover_2026_input_manifest.tsv`. Large files are gitignored.

---

## 中文

### “2026 leftover” 指什么
2024–2026 检索中 **PDAT=2026 的系列共 44 个**。完整目录见
`leftover_2026_catalog.tsv`。多数无法检验肿瘤 TACSTD2/CLDN4 与 ICI 疗效
（血液、细胞系、甲基化/miRNA、n=1、无标签、或唯一表达文件 > 2 GB）。

两个此前未分析、但同时具备开放已处理表达 **和** 样本级 ICI 标签的肿瘤系列：

| GSE | 此前漏掉的原因 | 本次使用 |
|-----|----------------|----------|
| **GSE329813** | 特征字段只有 `tissue`/`batch`；**MPR/NMPR 写在样本 title 里** | GeoMx 归一化 CSV；22 例 |
| **GSE292299** | 无单一 series matrix；3.1 GB RAW.tar 超限；**各样本 Visium H5 均 < 2 GB** | 16 例 NSCLC H5 + 元数据 |

### 结果

**GSE329813（新辅助帕博利珠+铂类，术后 GeoMx）**
- 标签来自 title，未杜撰。TACSTD2 在 1812 基因矩阵中；**CLDN4 不在**。
- 原发瘤床患者均值：MPR 11 例中位 3.94 vs NMPR 11 例 4.80；p = **0.002**；δ = **−0.79**（MPR 更低）。
- 引流淋巴结：p = 0.83（无差异）。
- **如实解读：** 这是**治疗后残留**组织。MPR 瘤床 TACSTD2 更低，符合残留上皮/肿瘤细胞更少，**不是**治疗前预测标志物。

**GSE292299（治疗前 Visium，16 例 NSCLC）**
- 12 R vs 4 NR。TACSTD2 p = 0.26（不显著）。CLDN4 p = **0.042**，与 PFS 的 Spearman ρ = 0.55（p = 0.027，忽略删失）。
- **如实解读：** CLDN4 差异主要由 4 例 NR 中的 2 例极低值驱动；另 2 例与 R 重叠。12 例 R 中 8 例是**脑转移**。全切片 Visium 混合肿瘤与间质。n_NR = 4，未做多重校正。属**探索性，非确证**。

### 2026 年已核查但未用于本问题的系列
血液/体液、无 ICI 标签的肿瘤 RNA、EGFR-TKI 耐药、细胞系/小鼠、超大文件且无更小已处理矩阵者，均见英文表，不在此重复编造可用性。
