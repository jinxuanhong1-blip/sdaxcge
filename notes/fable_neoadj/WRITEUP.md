# TACSTD2 / CLDN4 vs MPR in public neoadjuvant lung PD-1/PD-L1 datasets

**Slice:** `notes/fable_neoadj/`, `scripts/fable_neoadj/`, `results/fable_neoadj/`  
**PR:** https://github.com/jinxuanhong1-blip/sdaxcge/pull/57  
**Question:** Does malignant (or bulk-tumor) expression of the ADC targets **TACSTD2 (TROP2)** and **CLDN4 (Claudin-4)** associate with major pathologic response (MPR / pCR) after neoadjuvant PD-1/PD-L1 ± chemotherapy in NSCLC?

**Bottom line / 结论:** No. The only true baseline neoadjuvant-MPR dataset (GSE207422 bulk, n=24) shows a non-significant trend toward *lower* TACSTD2/CLDN4 in MPR. Post-treatment epithelial/malignant leftovers and the advanced-ICI leftover GSE205335 (RECIST, not MPR) do not agree with each other and are all non-significant. Honest mismatch: public GEO does not currently contain a second processed (<2 GB) neoadjuvant NSCLC ICI cohort with malignant-cell TACSTD2 and MPR labels.

---

## English

### 1. What was asked, and what exists on GEO

We screened GEO for open **neoadjuvant lung (NSCLC) PD-1/PD-L1 ± chemo** series with **MPR/pCR** and processed expression <2 GB, then tested **TACSTD2 / CLDN4 vs MPR**. Priority accessions: GSE207422, GSE243013 (if <2 GB), GSE146100, GSE229353, plus any others we could verify. A follow-up asked to add **GSE205335** and leftover series, focusing on **malignant TACSTD2 vs MPR**, and to keep an honest mismatch if the data do not fit.

Inventory (full table: `results/fable_neoadj/tables/dataset_inventory.csv`):

| Accession | Why it matters | Decision |
|---|---|---|
| **GSE207422** | Resectable NSCLC, neoadjuvant anti-PD-1 + platinum. Bulk log2TPM of **24 pre-treatment biopsies** with MPR/pCR/NMPR + residual %. scRNA of 15 samples (mostly post-tx). | **Primary** (bulk baseline) + leftover (scRNA malignant gate) |
| **GSE241934** | NEOTIDE/CTONG2104 + real-world neoadjuvant PD-1 ± chemo. scRNA of resected tumors, epithelial compartment, MPR labels. Matrices 466 MB / 1.25 GB. | **Validation** (post-tx epithelium, not baseline) |
| **GSE146100** | One multiprimary-LUAD patient, 3 nodules, induction pembrolizumab. W2 responded; W1/W3 did not. 229 MB. | **Descriptive only** (n=1 patient) |
| **GSE205335** | 26 lung-cancer patients on ICI; 28,512 published **malignant** cells; RDS 524 MB. | **Leftover / mismatch**: advanced/metastatic, **RECIST not MPR**, not neoadjuvant |
| GSE243013 | 234-patient post-chemo-IO immune atlas | **Excluded**: counts matrix **7.1 GB > 2 GB cap**; CD45/immune-focused |
| GSE229353 | 7 NSCLC, NAC vs pembrolizumab+chemo | **Excluded**: CD45+ immune-sorted — no epithelium |
| GSE280232 | KRAS-mut resectable NSCLC, neoadjuvant ICB | **Excluded**: TIL/CD8, recurrence/STK11, no epithelial MPR readout |
| GSE248378 | Durvalumab ± SBRT, bulk FFPE | **Excluded**: recurrence endpoint, no MPR/pCR on GEO |
| GSE225620 | Tislelizumab + chemo | **Excluded**: whole blood |
| GSE176021/176022 | Neoadjuvant PD-1 | **Excluded**: neoantigen-specific T cells |

So the **only processed public series that is simultaneously (i) neoadjuvant NSCLC ICI, (ii) MPR-labeled, (iii) baseline tumor, and (iv) <2 GB** is **GSE207422 bulk**. Everything else is either post-treatment residual tumor, n=1, immune-sorted, too large, or a different disease setting.

### 2. Methods (kept the same across datasets)

- Genes: **TACSTD2**, **CLDN4** (EPCAM used only for epithelial gating / tumor-content context).
- Contrast: MPR or pCR vs non-MPR. For GSE205335 leftover only: RECIST PR vs SD+PD.
- Test: two-sided Mann–Whitney U. Effect size = **AUC = P(marker higher in responder)** and Cliff's delta. AUC 0.5 = no association; AUC >0.5 = higher in responders.
- 95% CI on AUC: 2000 bootstrap resamples.
- scRNA: per-patient **epithelial or malignant pseudobulk**. GSE241934 uses the published `Epi` label. GSE207422 scRNA has **no published malignant labels on GEO**, so we gated `EPCAM>0 & PTPRC==0`. GSE205335 uses published `lineage.sub == Malignant cells`. QC: ≥20 epithelial/malignant cells per patient.
- Baseline and post-treatment / advanced-ICI cohorts are **not pooled** into one p-value. They answer different questions.

Scripts: `scripts/fable_neoadj/00_dataset_inventory.py` … `06_gse205335_malignant.py`, plus `04_meta_summary.py`.

### 3. Results

Master table: `results/fable_neoadj/tables/MASTER_auc_summary.csv`. Forest: `results/fable_neoadj/figures/MASTER_forest_auc.png`.

#### 3.1 Primary — GSE207422 baseline bulk (the only clean test)

24 pre-treatment biopsies, neoadjuvant anti-PD-1 + platinum, 9 MPR/pCR vs 15 non-MPR.

| Gene | median MPR | median non-MPR | AUC (95% CI) | p | Spearman vs residual tumor |
|---|---|---|---|---|---|
| TACSTD2 | 5.00 | 6.22 | 0.38 (0.15–0.61) | 0.34 | ρ = +0.25, p = 0.23 |
| CLDN4 | 5.85 | 6.27 | 0.36 (0.13–0.58) | 0.26 | ρ = +0.22, p = 0.31 |

Direction: both markers trend **lower** in eventual MPR (higher residual tumor ↔ higher expression). Neither is significant. n=24 is small; CIs include 0.5.

#### 3.2 Leftover — GSE207422 post-tx scRNA, EPCAM+/PTPRC−

15 samples (12 post-tx passing QC: 4 MPR vs 8 non-MPR). TACSTD2 AUC 0.34, p=0.46 (lower in MPR). CLDN4 AUC 0.66, p=0.46 (higher in MPR). Opposite directions, both noise. Pre-treatment scRNA is only P05/P08 (both non-MPR) plus P01 (NE) — not testable.

#### 3.3 Validation — GSE241934 post-tx epithelial pseudobulk

35 patients with ≥20 Epi cells, 10 MPR/pCR vs 25 non-MPR. **Post-resection**, so this is residual epithelium after therapy, not a baseline predictor.

| Gene | AUC (95% CI) | p | direction |
|---|---|---|---|
| TACSTD2 | 0.65 (0.43–0.84) | 0.17 | higher in MPR |
| CLDN4 | 0.50 (0.28–0.71) | 0.99 | null |

TACSTD2 here points the **opposite** way from GSE207422 baseline. That is a biological mismatch (residual vs baseline), not a replication.

#### 3.4 Leftover / mismatch — GSE205335 malignant cells vs RECIST

Not neoadjuvant. Not MPR. Stage IV / ED ICI, RECIST PR vs SD+PD. 16 patients with ≥20 malignant cells and a RECIST call (6 PR / 10 SD+PD). NSCLC-only (ADC/SQ): 4 PR vs 8 SD+PD.

| Subset | Gene | AUC | p |
|---|---|---|---|
| All histologies | TACSTD2 | 0.42 | 0.64 |
| All histologies | CLDN4 | 0.45 | 0.79 |
| NSCLC ADC/SQ | TACSTD2 | 0.56 | 0.81 |
| NSCLC ADC/SQ | CLDN4 | 0.44 | 0.81 |

SCLC cases have very low TACSTD2 (as expected for a neuroendocrine histology) and inflate the “all histologies” PR group. Even after dropping them, nothing is significant.

#### 3.5 Descriptive — GSE146100 (n=1)

EPCAM+ mean TACSTD2: W1 non-resp 1.08, W2 resp 1.09, W3 non-resp 1.48. The non-responding W3 is highest. Anecdote only.

### 4. Honest mismatch (do not over-read)

1. **Endpoint mismatch.** MPR/pCR is a surgical-pathology endpoint. GSE205335 is RECIST in advanced disease. We report it because it is the only <2 GB public lung-ICI scRNA object with published malignant labels and TACSTD2, not because it answers the neoadjuvant question.
2. **Timing mismatch.** GSE207422 bulk is baseline. GSE207422 scRNA, GSE241934, and GSE146100 are on- or post-treatment. Residual malignant TACSTD2 after a good pathologic response is a different quantity from pre-treatment TACSTD2.
3. **Compartment mismatch.** Bulk mixes stroma/immune; scRNA uses epithelium or published malignant cells. GSE207422 scRNA malignant calls are a marker gate, not author CNV/malignant labels.
4. **Direction mismatch.** Baseline bulk: TACSTD2 lower in MPR (n.s.). Post-tx GSE241934 epithelium: TACSTD2 higher in MPR (n.s.). These should not be meta-analyzed as one effect.
5. **Power.** Largest clean MPR contrast is 9 vs 15. No dataset is powered for a modest effect.
6. **GSE243013** would have been the large neoadjuvant scRNA resource (n=234) but the immune counts matrix is 7.1 GB and is not an epithelial/malignant atlas we can process under the cap.

### 5. Conclusion

Public processed GEO data do **not** support TACSTD2 or CLDN4 as predictors of MPR to neoadjuvant PD-1/PD-L1 ± chemo in NSCLC. The single baseline test is a weak, non-significant trend in the *opposite* direction from “high ADC-target → better pathologic response.” Leftover malignant-cell analyses are small, mixed, and (for GSE205335) the wrong endpoint. A real test needs a second baseline bulk or malignant-cell cohort with MPR labels — that object is not on GEO in processed form under 2 GB as of this slice.

---

## 中文

### 1. 问题与公开数据实际长什么样

本切片检索 GEO 上**开放的新辅助肺癌（NSCLC）PD-1/PD-L1 ± 化疗**、带 **MPR/pCR**、处理后表达矩阵 <2 GB 的系列，检验 ADC 靶点 **TACSTD2（TROP2）/ CLDN4（Claudin-4）与 MPR** 的关系。指定优先：GSE207422、GSE243013（若 <2 GB）、GSE146100、GSE229353，以及核验到的其他系列。后续补充 **GSE205335** 与 leftover，聚焦**恶性细胞 TACSTD2 vs MPR**，并允许如实报告错配。

完整清单见 `results/fable_neoadj/tables/dataset_inventory.csv`。

- **GSE207422**：可切除 NSCLC，新辅助抗 PD-1 + 铂类。**24 例治疗前活检** bulk log2TPM，有 MPR/pCR/NMPR 与残存肿瘤比例；另有 15 例 scRNA（多为术后）。**唯一合格的基线预测队列**；scRNA 作 leftover。
- **GSE241934**：NEOTIDE/CTONG2104 + 真实世界新辅助 PD-1 ± 化疗，术后肿瘤 scRNA，有上皮室与 MPR。矩阵 466 MB / 1.25 GB。**验证集，但是治疗后残存上皮，不是基线。**
- **GSE146100**：1 例多原发肺腺癌、3 个结节、诱导帕博利珠单抗。仅描述。
- **GSE205335**：26 例接受 ICI 的肺癌，28,512 个已注释恶性细胞，RDS 524 MB。**Leftover / 错配：晚期/转移、RECIST 不是 MPR、不是新辅助。**
- **排除**：GSE243013（counts 7.1 GB > 2 GB，且是免疫图谱）；GSE229353（CD45+ 分选，无上皮）；GSE280232（TIL/CD8，复发/STK11）；GSE248378（复发终点，无 MPR）；GSE225620（全血）；GSE176021/176022（新抗原 T 细胞）。

同时满足「新辅助 NSCLC ICI + MPR + 治疗前肿瘤 + <2 GB 已处理」的，只有 **GSE207422 bulk**。

### 2. 方法

- 基因：TACSTD2、CLDN4。
- 比较：MPR/pCR vs 非 MPR。GSE205335 leftover 仅用 RECIST PR vs SD+PD。
- 统计：双侧 Mann–Whitney U；效应量 AUC = P(应答者更高) 与 Cliff's δ。AUC=0.5 为无关联。
- scRNA：按患者做上皮/恶性伪 bulk。GSE241934 用作者 `Epi`；GSE207422 scRNA **GEO 无恶性注释**，用 `EPCAM>0 且 PTPRC==0`；GSE205335 用作者 `Malignant cells`。质控：每例 ≥20 个上皮/恶性细胞。
- **不把基线、治疗后、晚期 RECIST 合成一个 p 值。**

### 3. 结果

主表 `MASTER_auc_summary.csv`，森林图 `MASTER_forest_auc.png`。

**主分析 GSE207422 基线 bulk（24 例，9 MPR vs 15 非 MPR）**  
TACSTD2：AUC 0.38（0.15–0.61），p=0.34；与残存肿瘤 Spearman ρ=+0.25，p=0.23。  
CLDN4：AUC 0.36（0.13–0.58），p=0.26；ρ=+0.22，p=0.31。  
方向：MPR 者基线表达**偏低**（残存越多、表达越高）。均不显著。

**Leftover GSE207422 治疗后 scRNA（EPCAM+/PTPRC−，4 vs 8）**  
TACSTD2 AUC 0.34，p=0.46（MPR 更低）；CLDN4 AUC 0.66，p=0.46（MPR 更高）。方向相反，都是噪声。治疗前 scRNA 只有 2 例有 MPR 标签且均为非 MPR，无法检验。

**验证 GSE241934 治疗后上皮伪 bulk（10 vs 25）**  
TACSTD2 AUC 0.65（0.43–0.84），p=0.17（MPR 更高）；CLDN4 AUC 0.50，p=0.99。  
与 GSE207422 **基线方向相反**——这是「残存肿瘤」对「治疗前」的生物学错配，不是重复验证。

**Leftover / 错配 GSE205335 恶性细胞 vs RECIST**  
不是新辅助，不是 MPR。全组织学 6 PR vs 10 SD+PD：TACSTD2 AUC 0.42，p=0.64。仅 NSCLC（ADC/SQ，4 vs 8）：TACSTD2 AUC 0.56，p=0.81。SCLC 的 TACSTD2 极低（神经内分泌表型预期如此），会扭曲全队列。去掉后仍不显著。

**GSE146100（n=1）**  
应答结节 W2 的上皮 TACSTD2 并不高于非应答 W3。仅轶事。

### 4. 必须写明的错配

1. **终点错配**：MPR/pCR 是手术病理终点；GSE205335 是晚期 RECIST。收录它是因为它是目前唯一 <2 GB、带官方恶性注释且能读出 TACSTD2 的公开肺 ICI scRNA，**不是**因为它回答了新辅助问题。
2. **时点错配**：GSE207422 bulk 是基线；其余多为治疗中/后残存。
3. **区室错配**：bulk 混有间质/免疫；scRNA 用上皮或恶性细胞。GSE207422 scRNA 的「恶性」是标记门控，不是作者 CNV 注释。
4. **方向错配**：基线 bulk TACSTD2 在 MPR 偏低（不显著）；治疗后 GSE241934 上皮 TACSTD2 在 MPR 偏高（不显著）。不能合成一个效应。
5. **效能**：最干净的 MPR 对比是 9 vs 15，检测中等效应的能力不足。
6. **GSE243013** 本应是最大的新辅助 scRNA（n=234），但免疫 counts 7.1 GB，且不是上皮/恶性图谱，本切片无法处理。

### 5. 结论

现有已处理的公开 GEO 数据**不支持** TACSTD2 或 CLDN4 作为新辅助 PD-1/PD-L1 ± 化疗 MPR 的预测标志。唯一的基线检验是弱的、不显著的、且方向与「ADC 靶点高 → 病理缓解更好」相反。恶性细胞 leftover 样本小、方向乱；GSE205335 还是错误终点。要做实这个假设，需要第二个带 MPR 标签的**治疗前** bulk 或恶性细胞队列——该对象目前不以 <2 GB 已处理形式存在于 GEO。

---

## Reproducibility

```bash
# data files are downloaded from NCBI GEO FTP into results/fable_neoadj/data/ (gitignored)
python3 scripts/fable_neoadj/00_dataset_inventory.py
python3 scripts/fable_neoadj/01_gse207422_bulk.py
python3 scripts/fable_neoadj/02_gse241934_scrna.py
python3 scripts/fable_neoadj/03_gse146100_scrna.py
# R Matrix::readRDS -> results/fable_neoadj/tables/GSE205335_target_counts.csv
python3 scripts/fable_neoadj/05_gse207422_scrna_malignant.py
python3 scripts/fable_neoadj/06_gse205335_malignant.py
python3 scripts/fable_neoadj/04_meta_summary.py
```

Derived tables and figures live under `results/fable_neoadj/tables/` and `results/fable_neoadj/figures/`. Raw GEO downloads are gitignored.
