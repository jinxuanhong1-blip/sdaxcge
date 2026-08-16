# Malignant-restricted TACSTD2 vs T/NK and ICI response

**Datasets:** GSE207422 (MPR) and GSE205335 (RECIST). Processed GEO matrices only.
**Claim:** TACSTD2 scored only in malignant cells is higher in NMPR than MPR, and per-patient malignant TACSTD2 vs T/NK is negative with ρ ≈ −0.4 to −0.5.
**Verdict:** On **DRMref public labels** (12 post-tx patients; not Hu CopyKAT) the claim’s **direction and ρ magnitude are recovered** (NMPR mean 1.571 vs MPR 1.120, U=25, **p=0.154**; vs T/NK **ρ=−0.490, p=0.106**). Both tests are **not significant**. A marker-only malignant filter plus a <10-cell drop **misses** that ρ (near zero) because it NA’s the MPR samples DRMref still scores. **GSE205335 author malignant labels do not show the same pattern.** No invented numbers.

---

## English

### GEO verification

| Series | Design | Endpoint | Cells | Malignant definition used here |
|---|---|---|---|---|
| [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) | Neoadjuvant PD-1 + chemo NSCLC. Hu et al. *Genome Med* 2023, PMID 36869384. 3 pre-tx + 12 post-tx (MPR n=4 incl. pCR; NMPR n=8). | MPR / NMPR (pCR = MPR) | 92,330 on GEO | **Primary:** DRMref marker-based labels (Liu et al. *NAR* 2024): `GSE207422_Tor` + `GSE207422_Sin`, 30,877 post-tx cells, 2,051 `Malignant cells`. **Not** author CopyKAT (those IDs are not on GEO / supplements / author GitHub). **Sensitivity:** in-house epithelial minus normal-lung markers. |
| [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) | Lung cancer on ICI. 33 samples / 26 patients. ADC/SQ/SCLC/NUT. | RECIST PR=R; SD/PD=NR | 96,505 | Authors’ `lineage.sub = Malignant cells` (28,512). |

Processed files only; largest download 499.5 MB. GSE205335 RDS is double-gzipped. GSM→`orig.ident` 33/33. DRMref barcodes match GEO 30,877 / 30,877.

### GSE207422 — DRMref way (recovers ρ ≈ −0.49)

TACSTD2 = per-sample mean `log1p(CP10K)` in DRMref `Malignant cells`, using GEO UMIs joined on barcode. T/NK fraction = (CD8+ T + CD4+ T + NK) / DRMref-annotated cells. All 12 post-tx samples kept (P06 pCR has 15 malignant cells).

| Test | n | Result |
|---|---|---|
| NMPR vs MPR, malignant TACSTD2 | 8 vs 4 | mean **1.571 vs 1.120**; median 1.534 vs 1.048; U=25; **p=0.154** |
| Spearman vs T/NK fraction | 12 | **ρ=−0.490, p=0.106** |
| Spearman vs T/NK cytotoxicity (GZMB/GZMA/PRF1/IFNG/NKG7) | 12 | ρ=+0.154, p=0.633 |
| Spearman vs T/NK exhaustion (PDCD1/HAVCR2/LAG3/TIGIT/TOX) | 12 | ρ=+0.021, p=0.948 |
| Secondary: % TACSTD2+ malignant vs T/NK | 12 | ρ=−0.643, p=0.024 |
| Post-hoc drop P06 (<20 malignant) | 8 vs 3 | TACSTD2 p=0.048; vs T/NK ρ=−0.409, p=0.212 |

The unpublished slide’s ρ ≈ −0.4 to −0.5 vs **T/NK fraction** is the DRMref 12-patient Spearman. It is **not** significant at 0.05 (scipy two-sided p=0.106). It is **not** vs cytotoxicity or exhaustion (those ρ are ~0). NMPR>MPR is the same direction and also NS. Dropping P06 to get p=0.048 is post-hoc on n=3 MPR.

Per-patient malignant TACSTD2 (log1p CP10K): P02 1.039, P03 0.817, P04 2.070, P06 1.567, P07 2.153, P09 1.477, P10 1.526, P11 1.240, P12 0.929, P13 1.829, P14 0.856, P15 1.542.

### GSE207422 — marker-only sensitivity (why ρ looked like ~0)

Without DRMref, epithelial cells with tumor-epi score > normal-lung score were called malignant-like (9,429 vs 2,743 normal epi; TACSTD2 1.88 vs 1.37; SFTPA2 0.09 vs 2.53). Applying Hu’s <10-malignant-cell drop leaves **one MPR** (P03) vs 6 NMPR. Spearman vs T/NK is then n=9 **ρ=−0.100, p=0.798** (post-tx n=7 ρ=−0.036). That is a **definition mismatch**, not a second dataset: DRMref still assigns 15–87 malignant cells to P06/P11/P14, which the <10 rule deletes. All-epithelial TACSTD2 post-tx: NMPR 1.534 vs MPR 1.251, U=26, p=0.109, n=8 vs 4.

### GSE205335 — same pattern?

Author malignant TACSTD2, patient mean of tumor tissues (normals out).

| Test | n | Result |
|---|---|---|
| NR vs R, malignant TACSTD2 | 10 vs 6 | median 0.911 vs 1.098, U=32, **p=0.875** (opposite) |
| Spearman vs T/NK | 16 | **ρ=+0.447, p=0.083** (opposite) |
| ADC+SQ only vs T/NK | 12 | **ρ=−0.021, p=0.948** |
| NSCLC NR vs R TACSTD2 | 8 vs 4 | median 0.911 vs 1.644, p=0.368 |

**GSE205335 does not replicate** NMPR/NR-high TACSTD2 or a negative TACSTD2–T/NK ρ. The positive ρ is pulled by SCLC/NUT (SCLC PR tumors are TACSTD2-low and T/NK-low). Cell-level p≈0 values in `stats_summary.tsv` are pseudoreplication and are not used.

### What was not done

No CopyKAT (author IDs undeposited; 15 GB RAM). No extra ICI scRNA series beyond the assigned pair. TLS12 vs MPR/RECIST NS (post p=0.81; GSE205335 patient p=0.43).

### Files

- `scripts/fable_scrna_ici/04_drmref_recompute.py` — DRMref join + tests
- `results/fable_scrna_ici/gse207422_drmref_sample_metrics.tsv`, `gse207422_drmref_stats.tsv`
- `results/fable_scrna_ici/fig3_drmref_tacstd2.png`
- `results/fable_scrna_ici/drmref_cell_annotation.tsv.gz` — 30,877 cells (from DRMref Seurat `meta.data`)
- Marker-way: `gse207422_sample_metrics.tsv`, `stats_summary.tsv`, `fig1_*.png`, `fig2_*.png`
- `notes/fable_scrna_ici/GEO_VERIFICATION.md`

---

## 中文

### 结论

用 **DRMref 公开细胞标签**（术后 12 例，不是 Hu 的 CopyKAT）可以复现说法的**方向和 ρ 幅度**：恶性 TACSTD2 NMPR 均值 1.571 vs MPR 1.120，U=25，**p=0.154**；对 T/NK 比例 **ρ=−0.490，p=0.106**。两项都不显著。对 T/NK 杀伤/耗竭评分的 ρ 接近 0，不是 −0.4。只用 marker 再丢掉 <10 个恶性细胞时，MPR 几乎被删光，ρ 变成 −0.10——那是定义不同，不是第二个队列。**GSE205335 作者恶性标签没有同一模式**（对 T/NK ρ=+0.45；未缓解 vs 缓解 p=0.88，方向相反）。

### 数据

- GSE207422：GEO 只有样本表，没有 CopyKAT 细胞 ID。DRMref（*NAR* 2024）有 Tor+Sin 两套 Seurat，30,877 个术后细胞、2,051 个 `Malignant cells`，条形码与 GEO 全部对上。
- GSE205335：作者 `Malignant cells` 28,512；RECIST PR=缓解，SD/PD=未缓解。

### GSE207422（DRMref）

恶性细胞里 TACSTD2 的样本均值 `log1p(CP10K)`，T/NK = (CD8+T + CD4+T + NK) / 注释细胞。12 个术后样本全留（P06 pCR 仍有 15 个恶性细胞）。

- NMPR vs MPR：1.571 vs 1.120，U=25，p=0.154。
- vs T/NK 比例：ρ=−0.490，p=0.106（scipy 双侧）。这就是说法里的 −0.4 到 −0.5，但不显著。
- vs 杀伤 / 耗竭：ρ=+0.15 / +0.02。
- 次要：恶性细胞 TACSTD2 阳性率 vs T/NK，ρ=−0.643，p=0.024。
- 事后丢掉 P06：p=0.048（8 vs 3），不能当主结果。

### GSE205335

患者水平恶性 TACSTD2：未缓解中位 0.911 vs 缓解 1.098，p=0.875。vs T/NK ρ=+0.447，p=0.083。仅 ADC+SQ：ρ=−0.021。**不支持同一模式。**
