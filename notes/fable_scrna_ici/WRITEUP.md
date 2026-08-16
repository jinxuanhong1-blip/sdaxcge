# Malignant-restricted TACSTD2 vs T/NK and ICI response

**Datasets:** GSE207422 (MPR) and GSE205335 (RECIST). Processed GEO matrices only.
**Claim tested:** TACSTD2 scored only in malignant cells is higher in non-responders (NMPR / SD+PD) and correlates negatively with T/NK fraction at the patient level (ρ ≈ −0.4 to −0.5).
**Verdict:** The direction NMPR > MPR is seen in GSE207422 residual malignant cells but is **not testable** (only one post-tx MPR sample has ≥10 malignant-like cells). Patient-level Spearman vs T/NK is **near zero**, not −0.4 to −0.5. GSE205335 (author malignant labels) **does not show the same pattern** (TACSTD2 vs T/NK ρ = +0.45, p = 0.083; NR vs R p = 0.87, opposite median). No numbers below are invented.

---

## English

### GEO verification

| Series | Design | Endpoint | Cells | Malignant definition | Files |
|---|---|---|---|---|---|
| [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) | Neoadjuvant PD-1 + chemo, resectable NSCLC. Hu et al. *Genome Med* 2023, PMID 36869384. 3 pre-tx biopsies + 12 post-tx resections (MPR n=4 incl. pCR; NMPR n=8). | MPR / NMPR (pCR counted as MPR, as in the paper) | 92,330 | Authors: CopyKAT on epithelium (stromal reference). **GEO has no per-cell labels.** This slice: marker-score lineage, then epithelial cells whose tumor-epi program (EPCAM/KRTs/CEACAMs) exceeds a normal-lung program (SFTPA1/2, SFTPB/D, AGER, NAPSA, SCGB1A1/3A2, TPPP3, FOXJ1, CAPS) and normal-lung score < 0.4. | UMI matrix 175.5 MB; sample xlsx 11 KB |
| [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) | Lung cancer on ICI. 33 samples / 26 patients (LN, liver, effusion, lung/bronchus + a few normals). Subtypes ADC/SQ/SCLC/NUT. | RECIST PR = R; SD/PD = NR; NE dropped | 96,505 | Authors’ `lineage.sub = Malignant cells` (28,512 cells) | RDS 499.5 MB (double-gzipped); cell identity 719 KB |

Nothing >2 GB was downloaded. GSM→`orig.ident` for GSE205335 matched 33/33 (`results/fable_scrna_ici/gse205335_sample_table.tsv`). Hu et al. do **not** analyze TACSTD2; CopyKAT was not re-run (no author barcodes; dense 24k×92k matrix). Samples with <10 malignant(-like) cells have TACSTD2 means set to NA (paper dropped one NMPR sample with <10 malignant cells).

Lineage QC (GSE207422): epithelial n=12,172, EPCAM 1.54, TACSTD2 1.77, PTPRC 0.09; T/NK n=38,514, CD3E 2.14, TACSTD2 0.02. Malignant-like n=9,429 vs normal epi n=2,743: TACSTD2 1.88 vs 1.37; SFTPA2 0.09 vs 2.53; SCGB1A1 0.01 vs 1.50.

### GSE207422 — does NMPR > MPR for malignant TACSTD2?

Post-treatment residual malignant-like counts: MPR P03=1,207, P06=4, P11=1, P14=2; NMPR P02=2, P04=64, P07=4,845, P09=289, P10=323, P12=438, P13=4, P15=28. After the <10-cell rule, **MPR n=1** (P03 TACSTD2 mean = 0.806) vs **NMPR n=6** (median 1.746, mean 1.709). Mann–Whitney U is undefined with n=1. The single MPR point is below every NMPR value except P12 (0.806). Direction matches the claim; it is not a statistical test.

All-epithelial TACSTD2 (keeps every post-tx sample): NMPR median 1.534 vs MPR 1.251, U=26, **p=0.109**, n=8 vs 4. Same direction, not significant.

T/NK fraction post-tx: NMPR median 0.389 vs MPR 0.530, U=13, p=0.683, n=8 vs 4.

Cell-level (pseudoreplicated; do not use as evidence): NMPR median 2.202 vs MPR 0.792, 5,993 vs 1,214 cells, p≈0. Dominated by P07 (4,845 NMPR cells) vs P03 (1,207 MPR cells).

### GSE207422 — per-patient malignant TACSTD2 vs T/NK

One scRNA sample per patient. Spearman ρ (two-sided):

| Subset | n | TACSTD2 vs T/NK ρ | p | vs CD8 ρ | p |
|---|---|---|---|---|---|
| All patients with ≥10 malignant-like cells | 9 | **−0.100** | 0.798 | −0.367 | 0.332 |
| Post-tx only | 7 | **−0.036** | 0.939 | −0.393 | 0.383 |
| Post-tx, ≥30 malignant-like | 6 | **+0.029** | 0.957 | −0.371 | 0.468 |
| All 15 patients, all-epithelial TACSTD2 vs T/NK | 15 | −0.204 | 0.467 | — | — |

The claimed ρ ≈ −0.4 to −0.5 vs **T/NK is not reproduced**. The CD8 correlations sit near −0.37 with n≤9 and p>0.3. Selection bias: MPR/high-T/NK samples are exactly those with too few residual malignant cells to score TACSTD2.

### GSE205335 — same pattern?

Author malignant TACSTD2, patient-level mean of tumor-tissue samples (normals excluded).

| Test | n | Result |
|---|---|---|
| NR vs R, malignant TACSTD2 | 10 vs 6 | median NR 0.911 vs R 1.098, U=32, **p=0.875** (opposite of NMPR>MPR) |
| NR vs R, T/NK fraction | 10 vs 6 | median NR 0.292 vs R 0.188, U=39, p=0.368 |
| Spearman TACSTD2 vs T/NK (patients) | 16 | **ρ = +0.447, p = 0.083** (opposite sign) |
| Same, ADC+SQ only | 12 | **ρ = −0.021, p = 0.948** |
| Same, sample-level tumor tissues | 21 | ρ = +0.404, p = 0.069 |
| NSCLC only, NR vs R TACSTD2 | 8 vs 4 | median NR 0.911 vs R 1.644, p=0.368 |

GSE205335 does **not** show NMPR/NR-high TACSTD2 or a negative TACSTD2–T/NK correlation. The positive ρ is pulled by SCLC/NUT (SCLC PR tumors have very low TACSTD2 and low T/NK). NSCLC-only is null. Cell-level NR vs R (median 1.078 vs 0.000, p≈0) is pseudoreplication driven by those SCLC PR libraries and is not a patient-level result.

### What was not done

No other ICI scRNA series were added beyond the two assigned GEO accessions. No CopyKAT/inferCNV. No raw FASTQ. TLS 12-chemokine scores were computed but were not associated with MPR or RECIST (GSE207422 post p=0.81; GSE205335 patient p=0.43).

### Files

- `scripts/fable_scrna_ici/` — extractors + `03_analyze.py`
- `results/fable_scrna_ici/stats_summary.tsv` — every test
- `results/fable_scrna_ici/gse207422_sample_metrics.tsv`, `gse205335_patient_metrics.tsv`
- `results/fable_scrna_ici/fig1_response_boxplots.png`, `fig2_tacstd2_vs_tnk.png`
- `notes/fable_scrna_ici/GEO_VERIFICATION.md`

---

## 中文

### 要验证的说法

只在恶性细胞里算 TACSTD2：NMPR（或 RECIST 未缓解）高于 MPR/缓解，且患者水平恶性 TACSTD2 与 T/NK 比例呈负相关（ρ 约 −0.4 到 −0.5）。

### 数据（已核对 GEO；只用处理后矩阵；没有下 >2 GB 的文件）

- **GSE207422**（Hu 等，*Genome Medicine* 2023，PMID 36869384）：新辅助 PD-1 + 化疗、可切除 NSCLC。15 例各 1 个文库，92,330 细胞。作者用 CopyKAT 从上皮里分恶性（成纤维/内皮作对照）；GEO **没有**细胞标签。本切片用 marker 先分谱系，再把「肿瘤上皮程序 > 正常肺程序、且正常肺分数 < 0.4」的上皮当作恶性样细胞。pCR 并入 MPR（与原文一致）。恶性样细胞 <10 的样本 TACSTD2 均值记为缺失（原文恶性表达图也丢掉过 <10 细胞的 NMPR）。
- **GSE205335**：ICI 治疗中的肺癌，33 个样本 / 26 例患者，96,505 细胞。恶性细胞用作者公布的 `Malignant cells`（28,512）。RECIST：PR=缓解，SD/PD=未缓解。

GSE207422 谱系质控：上皮 12,172，EPCAM 1.54、TACSTD2 1.77、PTPRC 0.09；T/NK 38,514，CD3E 2.14、TACSTD2 0.02。恶性样 9,429 vs 正常上皮 2,743：TACSTD2 1.88 vs 1.37；SFTPA2 0.09 vs 2.53。

### GSE207422：NMPR 的恶性 TACSTD2 是否更高？

术后残留恶性样细胞：MPR 里只有 P03=1,207，其余 P06/P11/P14 为 4/1/2。应用 <10 规则后 **MPR 只剩 1 例**（P03 均值 0.806），NMPR 6 例（中位 1.746，均值 1.709）。**不能做 Mann–Whitney**。方向与「NMPR>MPR」一致，但不是检验。

全部上皮 TACSTD2（12 个术后样本都在）：NMPR 中位 1.534 vs MPR 1.251，U=26，**p=0.109**，8 vs 4。方向相同，不显著。术后 T/NK 比例：NMPR 0.389 vs MPR 0.530，p=0.683。

细胞水平（假重复，不能当证据）：NMPR 中位 2.202 vs MPR 0.792，5,993 vs 1,214 细胞，p≈0，几乎是 P07 对 P03。

### GSE207422：患者水平 ρ 是否约 −0.4 到 −0.5？

每例患者 1 个样本。Spearman：

- 恶性 TACSTD2 vs T/NK：全队列 n=9，**ρ=−0.100，p=0.798**；术后 n=7，**ρ=−0.036，p=0.939**；术后且恶性样≥30，n=6，**ρ=+0.029，p=0.957**。
- 与 CD8：n=9 时 ρ=−0.367，p=0.332；术后 n=7 时 ρ=−0.393，p=0.383。幅度接近说法，但样本少、不显著。
- 15 例全部上皮 TACSTD2 vs T/NK：ρ=−0.204，p=0.467。

**没有复现「对 T/NK 的 ρ ≈ −0.4 到 −0.5」。** MPR / 高 T/NK 样本恰恰是残留恶性细胞不够、TACSTD2 被设为缺失的那些，相关会被选择偏倚拉偏。

### GSE205335：有没有同一模式？

作者恶性标签、患者水平（只平均肿瘤组织，去掉 Normal*）：

- 恶性 TACSTD2，未缓解 vs 缓解：中位 0.911 vs 1.098，10 vs 6，**p=0.875**（方向相反）。
- TACSTD2 vs T/NK：**ρ=+0.447，p=0.083，n=16**（符号相反）。
- 仅 ADC+SQ：ρ=**−0.021**，p=0.948，n=12。
- 样本水平肿瘤组织：ρ=+0.404，p=0.069，n=21。

**GSE205335 不支持同一模式。** 正相关系数被 SCLC/NUT 拉高（SCLC 缓解瘤 TACSTD2 和 T/NK 都低）。细胞水平 NR 中位 1.078 vs R 0.000（p≈0）是假重复，不能当患者结论。

### 未做的事

没有再加这两个 GEO 之外的 ICI scRNA 数据集。没有重跑 CopyKAT。TLS 12 趋化因子评分与 MPR/RECIST 均无关联（GSE207422 术后 p=0.81；GSE205335 患者 p=0.43）。
