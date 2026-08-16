# TACSTD2 (TROP2) / CLDN4 protein and spatial omics vs immune contexture in lung cancer

**English first, then 中文.** All numbers below were computed from downloaded open processed files. CPTAC LUAD/LSCC are treatment-naive surgical cohorts: **there are no ICI response labels**.

Identifiers used (not invented): TACSTD2 / TROP2 = `ENSG00000184292` (UniProt P09758); CLDN4 = `ENSG00000189143` (UniProt O14493). Mouse symbols Tacstd2 / Cldn4 (UniProt Q8BGV3 / O35114) were searched only in PXD059688.

Scripts: `scripts/protein_spatial/`. Tables and figures: `results/protein_spatial/`. File-level catalog: `notes/protein_spatial_catalog.tsv`.

---

## English

### 1. CPTAC LUAD / LSCC tumor protein (freeze v1.2)

Downloaded the specified LUAD tumor protein matrix and the LSCC sibling, plus same-freeze `phenotype`, `survival`, and `meta` tables.

| Cohort | Tumors | TACSTD2 observed / NA | CLDN4 observed / NA | TACSTD2 median (log2) | CLDN4 median (log2) |
| --- | ---: | ---: | ---: | ---: | ---: |
| LUAD | 110 | 110 / 0 (0%) | 79 / 31 (28.2%) | 26.53 | 23.07 |
| LSCC | 108 | 108 / 0 (0%) | 78 / 30 (27.8%) | 26.87 | 23.34 |

CLDN4 protein missingness is real and large. It was **not** associated with ESTIMATE ImmuneScore, xCell immune score, CIBERSORT CD8, or age (Mann–Whitney, all p ≥ 0.072). Complete-case n is therefore used for all CLDN4 tests.

Immune scores in the freeze phenotype table are RNA-derived (ESTIMATE / xCell / CIBERSORT). No ICI, PD-1, or treatment-response field exists in `meta` or `survival`. OS/PFS are post-resection follow-up only.

**Spearman (BH-FDR within each correlation block):**

- LUAD TACSTD2 vs xCell_immune_score: n=110, ρ=−0.309, p=0.00100, FDR=0.032
- LUAD TACSTD2 vs xCell CD8 T cells: n=110, ρ=−0.289, p=0.00220, FDR=0.035
- LUAD TACSTD2 vs ESTIMATE ImmuneScore: n=110, ρ=−0.185, p=0.052, FDR=0.305 (not significant after FDR)
- LSCC CLDN4 vs ESTIMATE ImmuneScore: n=78, ρ=−0.432, p=7.9×10⁻⁵, FDR=0.0025
- LSCC CLDN4 vs xCell_immune_score: n=78, ρ=−0.400, p=2.9×10⁻⁴, FDR=0.0031
- LSCC CLDN4 vs protein CD8A: n=78, ρ=−0.444, p=4.6×10⁻⁵
- Pooled LUAD+LSCC CLDN4 vs ESTIMATE ImmuneScore: n=157, ρ=−0.325, p=3.3×10⁻⁵, FDR=5.2×10⁻⁴
- Pooled TACSTD2 vs ESTIMATE ImmuneScore: n=218, ρ=−0.104, p=0.126, FDR=0.230 (not significant)

Median-split Mann–Whitney on ESTIMATE ImmuneScore: LUAD TACSTD2 high vs low, n=55/55, median 7862 vs 8712, p=0.017; LSCC CLDN4 high vs low, n=39/39, median 7693 vs 8672, p=0.013. Direction is the same: higher TACSTD2/CLDN4 protein, lower RNA immune score.

**Survival (not ICI OS):** events are few. LUAD TACSTD2 OS log-rank p=0.71 (23 events / 105); LUAD CLDN4 OS log-rank p=0.087 (15 events / 75); LUAD TACSTD2 PFS Cox HR per SD=1.51, p=0.071. LSCC TACSTD2 OS log-rank p=0.30 (23 events / 94). No endpoint met p<0.05.

Figures: `figures/cptac_missingness.*`, `cptac_scatter_immune.*`, `cptac_LUAD_corr_heatmap.*`, `cptac_LSCC_corr_heatmap.*`, `cptac_os_km.*`.

### 2. PXD042091 ICI plasma SWATH (human metastatic NSCLC)

Processed files downloaded: `6_data_tf_response.txt` (6364 proteins; responder/nonresponder), `2_dat_cs_all.txt` (7115 proteins), `NSClibrary.txt` (11,625 protein names). All raw `wiff`/`wiff.scan` skipped.

- Spectral library **contains** TACD2_HUMAN (GN=TACSTD2) and CLD4_HUMAN (GN=CLDN4).
- Neither quantified matrix contains P09758, O14493, TACSTD2/TACD2, or CLDN4/CLD4.
- Therefore TACSTD2 and CLDN4 are **not present as quantified plasma proteins** in the open SWATH tables. No ICI-response test was possible for these two proteins.

### 3. PXD059688 mouse LLC anti-PD1 ± ascorbate

Only `20230904_Tumor_AA.mzTab.gz` (323 KB) is an open processed table under 2 GB. `20230904_Tumor_AA.mgf` is 8.49 GB and `20230904_Tumor_AA.msf` is 37.9 GB — both skipped.

The mzTab has **93 PRT rows**, almost all methyltransferases (Dnmt1, Nsun2, Prmt1, Mettl3, …). Tacstd2 / Cldn4 / Q8BGV3 / O35114 are **absent**. This export is not a full LLC proteome, so absence from the full experiment cannot be claimed.

### 4. GSE271689 GeoMx WTA (ICI NSCLC)

Downloaded series matrix, SOFT, and `GSE271689_RAW.tar` (586 DCC files, 35.1 MB). SRA FASTQ skipped. RTS IDs were mapped with the official NanoString Hs WTA v1.0 PKC (`Hs_R_NGS_WTA_v1.0.pkc`; TACSTD2=RTS0027086, CLDN4=RTS0026256). This PKC is a panel map, not a new lung-cancer accession.

GEO characteristics: spotid, tissue=Lung Cancer, cell type = CK / CD45 / CD68, treatment = Immunotherapy. **No OS, PFS, vital status, or response field** is present in the series matrix or SOFT (25,828 SOFT lines checked). ICI OS analysis was therefore **not performed**.

After dropping NTC/NA AOIs: n=579 (CK 212, CD68 192, CD45 175). Detection (count>0): TACSTD2 94.9–99.5% of AOIs; CLDN4 96.0–99.1%.

Aligned-read CPM (epithelial vs leukocyte):

- TACSTD2 median CPM: CK 11.48 vs CD45 4.58, Mann–Whitney p=6.5×10⁻¹⁶
- CLDN4 median CPM: CK 10.36 vs CD45 4.31, Mann–Whitney p=2.3×10⁻¹⁸

AOI-level Spearman of CPM(TACSTD2 or CLDN4) vs CPM(CD8A/PTPRC/…) remains high (CK TACSTD2–CD8A ρ=0.53, n=212). This is **not** interpreted as a clean immune-contexture association: only aligned-total CPM was available (not Q3 of the full ~18k-gene matrix), and AOI area/nuclei are not in GEO. The trustworthy GeoMx result is panel presence + epithelial enrichment.

### 5. E-MTAB-13530 Visium WTA (NSCLC, no ICI in SDRF)

Downloaded 40 Space Ranger `filtered_feature_bc_matrix.h5` files plus SDRF/IDF. FASTQ, `spatial.tar`, and HTML summaries skipped. Both genes are in the WTA matrix (36,601 features).

| Class | Sections | TACSTD2 median detect % | CLDN4 median detect % |
| --- | ---: | ---: | ---: |
| Tumor | 20 | 48.2 | 47.9 |
| Adjacent non-tumor | 16 | 19.3 | 14.8 |
| Healthy donor | 4 | 39.5 | 17.3 |

Tumor vs adjacent detection: TACSTD2 Mann–Whitney p=1.9×10⁻⁵; CLDN4 p=1.6×10⁻⁴.

Spot-level Spearman inside tumor sections is weak for T-cell genes (median across 20 sections: TACSTD2–CD8A ρ=0.037; CLDN4–CD8A ρ=0.029) and modest for CD68 (TACSTD2–CD68 median ρ=0.220; CLDN4–CD68 median ρ=0.183). SDRF has disease and sampling site only — **no ICI labels**.

### 6. Targeted Xenium-IO / CosMx-1K

Skipped as instructed. No CosMx/Xenium accessions were invented.

### 7. What can and cannot be concluded

Supported: in treatment-naive CPTAC tumors, CLDN4 protein (complete cases) is inversely associated with RNA immune / CD8 scores, strongest in LSCC; TACSTD2 protein is inversely associated with xCell immune/CD8 scores in LUAD. CLDN4 protein is missing in ~28% of tumors. Both genes are present and epithelium-enriched in GeoMx/Visium WTA.

Not supported: any ICI-response claim from CPTAC; plasma TACSTD2/CLDN4 tests in PXD042091; Tacstd2/Cldn4 quantification in PXD059688; GeoMx ICI OS (metadata absent).

---

## 中文

### 1. CPTAC LUAD / LSCC 肿瘤蛋白（freeze v1.2）

已下载指定的 LUAD 肿瘤蛋白矩阵及 LSCC 对应文件，并合并同一 freeze 的 `phenotype` / `survival` / `meta`。

| 队列 | 肿瘤数 | TACSTD2 检出/缺失 | CLDN4 检出/缺失 | TACSTD2 中位数 (log2) | CLDN4 中位数 (log2) |
| --- | ---: | ---: | ---: | ---: | ---: |
| LUAD | 110 | 110 / 0（0%） | 79 / 31（28.2%） | 26.53 | 23.07 |
| LSCC | 108 | 108 / 0（0%） | 78 / 30（27.8%） | 26.87 | 23.34 |

CLDN4 蛋白缺失是真实且大量的。缺失与 ESTIMATE ImmuneScore、xCell immune score、CIBERSORT CD8、年龄均无关联（Mann–Whitney，全部 p≥0.072）。后续 CLDN4 检验均用完全观测样本。

表型文件中的免疫评分来自 RNA（ESTIMATE / xCell / CIBERSORT）。`meta` 与 `survival` **没有 ICI / 免疫治疗疗效字段**。OS/PFS 仅为术后随访。

**Spearman（各相关块内 BH-FDR）：**

- LUAD TACSTD2 vs xCell_immune_score：n=110，ρ=−0.309，p=0.00100，FDR=0.032
- LUAD TACSTD2 vs xCell CD8：n=110，ρ=−0.289，p=0.00220，FDR=0.035
- LUAD TACSTD2 vs ESTIMATE ImmuneScore：n=110，ρ=−0.185，p=0.052，FDR=0.305（FDR 后不显著）
- LSCC CLDN4 vs ESTIMATE ImmuneScore：n=78，ρ=−0.432，p=7.9×10⁻⁵，FDR=0.0025
- LSCC CLDN4 vs xCell_immune_score：n=78，ρ=−0.400，p=2.9×10⁻⁴，FDR=0.0031
- LSCC CLDN4 vs 蛋白 CD8A：n=78，ρ=−0.444，p=4.6×10⁻⁵
- 合并 LUAD+LSCC CLDN4 vs ESTIMATE ImmuneScore：n=157，ρ=−0.325，p=3.3×10⁻⁵，FDR=5.2×10⁻⁴
- 合并 TACSTD2 vs ESTIMATE ImmuneScore：n=218，ρ=−0.104，p=0.126，FDR=0.230（不显著）

按中位数分组的 Mann–Whitney（ESTIMATE ImmuneScore）：LUAD TACSTD2 高/低 n=55/55，中位数 7862 vs 8712，p=0.017；LSCC CLDN4 高/低 n=39/39，中位数 7693 vs 8672，p=0.013。方向一致：蛋白越高，RNA 免疫评分越低。

**生存（不是 ICI OS）：** 事件数少。LUAD TACSTD2 OS log-rank p=0.71（23/105）；LUAD CLDN4 OS log-rank p=0.087（15/75）；LUAD TACSTD2 PFS Cox 每 SD HR=1.51，p=0.071。LSCC TACSTD2 OS log-rank p=0.30（23/94）。无一终点 p<0.05。

### 2. PXD042091 ICI 血浆 SWATH（人转移性 NSCLC）

已下载 `6_data_tf_response.txt`（6364 蛋白）、`2_dat_cs_all.txt`（7115 蛋白）、`NSClibrary.txt`（11625 个蛋白名）。原始 wiff 已跳过。

谱图库**含有** TACD2_HUMAN 与 CLD4_HUMAN；两张定量表均**不含** P09758 / O14493 / TACSTD2 / CLDN4。无法对这两个蛋白做 ICI 疗效检验。

### 3. PXD059688 小鼠 LLC 抗 PD1 ± 抗坏血酸

唯一 <2GB 的开放处理文件是 `20230904_Tumor_AA.mzTab.gz`。mgf（8.49 GB）与 msf（37.9 GB）已按规则跳过。mzTab 仅 93 条 PRT，几乎全是甲基转移酶。Tacstd2 / Cldn4 **不在该表中**。这不是全蛋白组，故不能声称全实验中不存在。

### 4. GSE271689 GeoMx WTA（ICI NSCLC）

已下载 series matrix、SOFT 与 DCC tar。SRA FASTQ 已跳过。用官方 Hs WTA v1.0 PKC 映射 RTS（TACSTD2=RTS0027086，CLDN4=RTS0026256）。

GEO 仅有 spotid / 组织 / 细胞类型（CK、CD45、CD68）/ treatment=Immunotherapy。**没有 OS/PFS/疗效字段**（SOFT 25828 行已检查），因此**未做 ICI OS 分析**。

过滤 NTC 后 n=579。CK 区 TACSTD2/CLDN4 的 aligned-read CPM 高于 CD45（Mann–Whitney p=6.5×10⁻¹⁶ 与 2.3×10⁻¹⁸）。AOI 水平与免疫基因的 CPM 相关仍然偏高，因缺少全转录组 Q3 与面积/核数，**不解释为可靠的免疫相关**。可确认的是：WTA 包含两基因，且在上皮（CK）区更富集。

### 5. E-MTAB-13530 Visium WTA

已下载 40 个 Space Ranger h5 与 SDRF。两基因均在矩阵中。肿瘤切片检出率高于癌旁（TACSTD2 p=1.9×10⁻⁵；CLDN4 p=1.6×10⁻⁴）。肿瘤切片内与 CD8A 的 spot 水平相关很弱（TACSTD2 中位 ρ=0.037），与 CD68 中等（中位 ρ=0.220）。SDRF **无 ICI 标签**。

### 6. 靶向 Xenium-IO / CosMx-1K

按要求跳过，未编造 accession。

### 7. 结论边界

支持：治疗初治 CPTAC 中，CLDN4 蛋白（完全观测）与 RNA 免疫/CD8 评分负相关，LSCC 更强；LUAD 中 TACSTD2 蛋白与 xCell 免疫/CD8 负相关。CLDN4 蛋白约 28% 缺失。空间 WTA 中两基因存在且偏上皮。

不支持：用 CPTAC 推断 ICI 疗效；PXD042091 血浆定量这两蛋白；PXD059688 全蛋白组中 Tacstd2/Cldn4 的有无；GSE271689 的 ICI OS（GEO 无生存字段）。
