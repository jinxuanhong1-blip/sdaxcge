# PRIDE / ProteomeXchange lung ICI proteomes — TACSTD2 / CLDN4

Slice outputs only: `notes/grok_pride/`, `scripts/grok_pride/`, `results/grok_pride/`.  
No raw MS (`.raw`, `.wiff`, `.wiff.scan`, `.d`) was downloaded. All numbers below were computed from processed identification/quantification tables or journal supplements.

---

# English

## Question

Are TACSTD2 (TROP2) and CLDN4 (claudin-4) present as proteins in open lung immune-checkpoint-inhibitor (ICI) proteomes, starting with PXD042091, PXD059688, PXD019061 and any other verified open lung ICI datasets?

## Verified open lung ICI proteomes

| Accession | Host | Design | TACSTD2 | CLDN4 |
|---|---|---|---|---|
| **PXD042091** | PRIDE | Human plasma SWATH-MS; metastatic NSCLC + pembrolizumab (n=64; 171 samples) | **Absent** (0/7115) | **Absent** (0/7115) |
| **PXD059688** | PRIDE | LLC1 mouse tumors ± high-dose AA ± anti-PD1; paper: 6737 proteins | **Not in deposited mzTab** (0/93 methyltransferase subset) | **Not in deposited mzTab** |
| **PXD019061** | PanoramaPublic / MassIVE | FFPE NSCLC; targeted checkpoints (n=46) + global spectral counts (n=26); ICI-naive tissue | **Present 24/26** (median 10 counts) | **Absent** (not in 6819 gene symbols) |
| **PXD039141** | iProX | Serum intact glycoproteomics; advanced lung cancer on anti-PD-1/PD-L1 | **Absent** (0/119 glycoproteins) | **Absent** |
| **PXD019573** (+ PXD028364, PXD037365) | PRIDE | NSCLC proteasome/degradome; Durvalumab association in paper | Not in proteasome-centric tables | Not in proteasome-centric tables |

Additional PRIDE hits from keyword search (PXD044740, PXD019774, PXD034772, PXD020191, PXD035347) were inspected and **excluded** as not verified lung-ICI tumor/plasma proteomes (T-cell cytokine phospho, cell-line immunopeptidome, inflammation immunopeptidome, ICI-untreated NSCLC proteogenomics, or melanoma). See `results/grok_pride/dataset_inventory.json`.

## Real findings

### 1. PXD042091 — circulating SWATH-MS, pembrolizumab NSCLC

Processed files from PRIDE FTP (2025/02/PXD042091), not raw wiff:

- `2_dat_cs_all.txt`: **7115** proteins × patient/timepoint columns; UniProt mnemonic `ACCESSION_GENENAME`.
- `6_data_tf_response.txt`: **6364** proteins with baseline responder (R*) vs non-responder (NR*) intensities, `p_val`, `foldchange`.

**TACSTD2 / CLDN4:** no row matches P09758, O14493, TACSTD2, TACD2, TROP2, CLDN4, or CLD4.

This is not a failed parse. The same table contains the paper’s seven-protein signature, all with p&lt;0.01, and **323** proteins with p&lt;0.05 (paper: 324 DEPs):

| Protein | Gene | median R | median NR | foldchange | p |
|---|---|---|---|---|---|
| O15020_SPTN2 | SPTBN2 | 13.49 | 12.57 | +1.11 | 0.00124 |
| Q92696_PGTA | RABGGTA | 12.93 | 12.12 | +1.03 | 0.00104 |
| Q4L180_FIL1L | FILIP1L | 13.84 | 15.31 | −1.29 | 0.00109 |
| Q9UHG0_DCDC2 | DCDC2 | 11.31 | 12.62 | −1.26 | 0.00254 |
| Q7Z3C6_ATG9A | ATG9A | 13.82 | 14.64 | −0.84 | 0.00307 |
| Q9NQ48_LZTL1 | LZTFL1 | 12.19 | 13.51 | −1.29 | 0.00501 |
| Q9UPZ3_HPS5 | HPS5 | 11.78 | 12.49 | −0.85 | 0.00868 |

Junction / epithelial proteins **are** quantified in plasma, but none separate R vs NR at p&lt;0.05:

| Protein | p | foldchange | median R | median NR |
|---|---|---|---|---|
| P16422_EPCAM | 0.35 | −0.40 | 14.79 | 15.42 |
| Q9NY35_CLDN1 | 0.53 | −0.30 | 11.65 | 12.56 |
| O00501_CLD5 | 0.30 | +0.51 | 13.64 | 13.17 |
| P56856_CLD18 | 0.35 | −0.44 | 12.19 | 12.28 |
| Q16625_OCLN | 0.75 | +0.14 | 12.08 | 12.06 |
| P15941_MUC1 | 0.25 | −0.52 | 11.95 | 12.07 |

**Interpretation:** in this pembrolizumab plasma proteome, TACSTD2 and CLDN4 are below the SWATH quantitation set, while other claudins and EpCAM circulate and do not track response. Plasma is a poor matrix for these two membrane proteins.

### 2. PXD059688 — LLC1 anti-PD1 ± ascorbic acid

PRIDE processed export is `20230904_Tumor_AA.mzTab.gz` (324 KB). The companion `.mgf` (7.9 GB peak lists) and `.msf` (35 GB Proteome Discoverer DB) were **not** downloaded.

The mzTab is a Proteome Discoverer “minimal Summary Quantification report” with **93 PRT** rows. Every accession is a methyltransferase (Dnmt1, Prmt1, Mettl3, …). PEP=665, PSM=6935. **Zero** matches to Tacstd2/Q8BGV3 or Cldn4/O35054 in PRT/PEP/PSM.

Kim et al. report 6737 proteins / 6209 in all four groups and point to Supplementary Table S1. That full matrix is **not** in the PRIDE processed files. Frontiers HTML advertises `DataSheet1.docx` but file URLs were not openly retrievable in this environment.

**Interpretation:** presence of Tacstd2/Cldn4 in the *full* LLC1 tumor proteome is **unconfirmed from repository processed data**. The deposited mzTab cannot answer the question.

### 3. PXD019061 — FFPE NSCLC checkpoint + global proteome

Host is PanoramaPublic (targeted Skyline), not PRIDE. Global MS is also MassIVE `MSV000085049`. We used the Sci Rep supplements (no raw MS):

- MOESM1: targeted fmol/µg for PD-1, PD-L1, PD-L2, IDO1, LAG3, TIM-3, ICOSLG, VISTA, GITR, CD40 only. **Neither TACSTD2 nor CLDN4 is in the PRM panel.**
- MOESM5 Table S6: **6819** gene symbols × **26** tumors, spectral counts (paper: 12,301 protein groups at 2-peptide / 3.43% protein FDR).

**TACSTD2 is present.** Detection 24/26 (92.3%), median 10, mean 11.3, max 29 counts. The two zeros are SR12-96 (adenocarcinoma) and V576 (squamous).

**CLDN4 is absent** from Table S6. Other claudins in the same matrix: CLDN1 3/26, CLDN3 4/26, CLDN11 1/26. CLDN5/7/18 are also missing. EPCAM 13/26 (median 0.5); CDH1 and TJP1 26/26.

TACSTD2 vs targeted PD-L1 protein (overlap n=23): Spearman **ρ = 0.28** (weak). By histology in the 26-tumor global set:

| Histology | n | TACSTD2 detected | median counts |
|---|---|---|---|
| Squamous cell carcinoma | 10 | 9/10 | **15.5** |
| Adenocarcinoma | 11 | 10/11 | **8.0** |
| Other / NA | 5 | 5/5 | 2–10 |

**Interpretation:** in ICI-naive NSCLC *tissue*, TACSTD2 protein is almost ubiquitous and higher in squamous than adenocarcinoma. CLDN4 is not recovered by this global DDA spectral-count experiment (low abundance or poor peptide observability in FFPE), so tissue presence of CLDN4 is **not supported** here. The study is checkpoint-focused and ICI-naive: it confirms protein observability, not ICI-response association.

### 4. Other verified open lung ICI proteomes

**PXD039141** (iProX IPX0003162000): serum intact glycopeptides during anti-PD-1/PD-L1. Processed search files: 119 protein accessions, plasma/Ig dominant (`CO3`, `IGHM`, `ANGT`, …). No TACSTD2/CLDN4.

**PXD019573 / PXD028364 / PXD037365** (Javitt et al., Nat Cancer 2023): NSCLC proteasome footprinting; PA200/PSME4 vs PSMB10 associated with Durvalumab. Nature Cancer source tables are proteasome-subunit / cleavage-pattern tables. TACSTD2/CLDN4 do not appear. Open and ICI-related, but not a TACSTD2/CLDN4 matrix.

## Cross-dataset conclusion

1. **TACSTD2 protein is real in lung tumor proteomes** (PXD019061 global FFPE: 24/26), not in pembrolizumab **plasma** SWATH (PXD042091) and not in serum glyco-ICI (PXD039141).
2. **CLDN4 protein is not detected** in any processed open table we could actually search (plasma 7115; FFPE 6819 genes; serum glyco 119; LLC1 mzTab 93).
3. The only ICI-treatment *response* matrix with a full protein list (PXD042091) has neither target; circulating EpCAM/claudin-1/5/18 do not track pembrolizumab response.
4. PXD059688 cannot be used for Tacstd2/Cldn4 until the 6737-protein Supplementary Table S1 is obtained; the PRIDE mzTab is the wrong protein set.

## Methods (no raw MS)

- PRIDE WS v2 project metadata + ProteomeXchange JSON/XML.
- HTTPS FTP listings; download only SEARCH / mzTab / journal xlsx / iProX search txt.
- Target aliases: P09758 / TACSTD2 / TACD2 / TROP2; O14493 / CLDN4 / CLD4; mouse Q8BGV3 / O35054.
- PXD042091 R vs NR p-values taken from the authors’ table (reproduced 323 vs 324 DEPs).
- PXD019061 Spearman on ranks of TACSTD2 spectral counts vs Table S1 PD-L1 fmol/µg.

Scripts: `scripts/grok_pride/`. Tables: `results/grok_pride/presence_table.tsv`, `PXD019061_TACSTD2_vs_PDL1.tsv`, `PXD042091_signature_and_junction.tsv`.

## Limits

- PXD059688 full proteome not in the repository processed export.
- PXD019061 is ICI-naive; spectral counts are not absolute quantification.
- Plasma SWATH library (`NSClibrary.txt`, 268 MB) was not downloaded; absence in the 7115-protein quant table already means the proteins were not quantified for patients.
- No claim that CLDN4 protein is biologically absent from NSCLC — only that it is missing from these open processed matrices.

---

# 中文

## 问题

在开放的肺免疫检查点抑制剂（ICI）蛋白质组中，TACSTD2（TROP2）与 CLDN4（claudin-4）是否以蛋白质形式被检出？优先核验 PXD042091、PXD059688、PXD019061，并纳入其他已核实的开放肺 ICI 蛋白质组。不下载原始质谱。

## 已核实的开放肺 ICI 蛋白质组

| 登录号 | 托管 | 设计 | TACSTD2 | CLDN4 |
|---|---|---|---|---|
| **PXD042091** | PRIDE | 人血浆 SWATH-MS；转移性 NSCLC + 帕博利珠单抗（64 例，171 份血样） | **未检出**（0/7115） | **未检出**（0/7115） |
| **PXD059688** | PRIDE | LLC1 小鼠肿瘤 ± 高剂量抗坏血酸 ± anti-PD1；原文 6737 蛋白 | **不在已提交 mzTab 中**（0/93，甲基转移酶子集） | **不在已提交 mzTab 中** |
| **PXD019061** | PanoramaPublic / MassIVE | FFPE NSCLC；靶向检查点（n=46）+ 全局谱图计数（n=26）；组织为 ICI 初治 | **检出 24/26**（中位 10 counts） | **未出现**（6819 个基因符号中无） |
| **PXD039141** | iProX | 血清完整糖蛋白组；晚期肺癌抗 PD-1/PD-L1 治疗前后 | **未检出**（0/119） | **未检出** |
| **PXD019573**（及 PXD028364、PXD037365） | PRIDE | NSCLC 蛋白酶体/降解组；论文中与度伐利尤单抗反应相关 | 不在蛋白酶体中心表格中 | 不在蛋白酶体中心表格中 |

关键词检索到的 PXD044740、PXD019774、PXD034772、PXD020191、PXD035347 已核验并**排除**（T 细胞细胞因子磷酸蛋白质组、细胞系免疫肽组、非 ICI 治疗队列或黑色素瘤）。详见 `dataset_inventory.json`。

## 实测结果

### 1. PXD042091：帕博利珠单抗 NSCLC 循环 SWATH-MS

来自 PRIDE FTP 的 SEARCH 表（非 raw wiff）：7115 个定量蛋白；基线 R vs NR 表 6364 行。

**TACSTD2 / CLDN4 均不在表中**（P09758、O14493、TACD2、CLD4 等别名全无）。解析可靠：同一张表复现了原文 7 蛋白签名（均 p&lt;0.01），且 p&lt;0.05 蛋白数为 **323**（原文 324）。

血浆中**能量化到** EPCAM、CLDN1、CLDN5、CLDN18、OCLN、MUC1，但 R vs NR 均不显著（最小 p 为 MUC18 的 0.096）。

**含义：** 在该帕博利珠单抗血浆蛋白质组中，TACSTD2/CLDN4 未进入 SWATH 定量集合；其他紧密连接/上皮蛋白虽可循环，但不随疗效分层。血浆不是这两个膜蛋白的合适基质。

### 2. PXD059688：LLC1 anti-PD1 ± 抗坏血酸

仓库中唯一可公开拉取的加工结果是 324 KB 的 mzTab。`.mgf`（7.9 GB）与 `.msf`（35 GB）未下载。mzTab 仅含 **93** 个甲基转移酶，PRT/PEP/PSM 中无 Tacstd2/Cldn4。原文 6737 蛋白矩阵在补充表 S1，PRIDE 加工文件中没有，Frontiers 补充文件 URL 在本环境无法直接获取。

**含义：** 就仓库已提交的加工数据而言，**无法确认** Tacstd2/Cldn4；不能把 mzTab 当成全蛋白质组。

### 3. PXD019061：FFPE NSCLC 检查点 + 全局蛋白质组

靶向 PRM 面板只有 PD-1/PD-L1/PD-L2/IDO1/LAG3/TIM-3/ICOSLG/VISTA/GITR/CD40，**不含** TACSTD2/CLDN4。

全局 Table S6（26 例肿瘤，6819 基因符号）：

- **TACSTD2：24/26 检出**（92.3%），中位谱图数 10，最大 29。
- **CLDN4：表中不存在。** 同表 CLDN1 仅 3/26、CLDN3 仅 4/26。
- 鳞癌中位 15.5（n=10）高于腺癌 8.0（n=11）。
- 与靶向 PD-L1 蛋白 Spearman **ρ=0.28**（n=23，弱相关）。

**含义：** ICI 初治肺癌**组织**中 TACSTD2 蛋白几乎普遍存在，鳞癌高于腺癌。CLDN4 在该 FFPE DDA 谱图计数实验中未被回收（低丰度或肽段可观测性差）。本研究不能回答 ICI 疗效关联，只能回答组织可观测性。

### 4. 其他已核实开放集

**PXD039141：** 抗 PD-1/PD-L1 期血清完整糖肽，119 个蛋白登录号，以血浆蛋白/免疫球蛋白为主，无 TACSTD2/CLDN4。

**PXD019573 系列：** 与度伐利尤单抗相关的蛋白酶体/抗原加工研究，开放且属肺 ICI 相关，但表格是蛋白酶体亚基/切割谱，不是 TACSTD2/CLDN4 矩阵。

## 跨库结论

1. **TACSTD2 蛋白在肺癌组织蛋白质组中是真实存在的**（PXD019061：24/26），但在帕博利珠单抗**血浆** SWATH（PXD042091）和 ICI 血清糖蛋白组（PXD039141）中未定量到。
2. **CLDN4 蛋白在所有实际可检索的开放加工表中均未出现**（血浆 7115、FFPE 6819 基因、血清糖蛋白 119、LLC1 mzTab 93）。
3. 唯一带疗效标签的全蛋白矩阵（PXD042091）两个靶标都没有；循环 EpCAM/CLDN1/5/18 不区分帕博利珠单抗反应。
4. PXD059688 在拿到 6737 蛋白补充表之前，不能用于 Tacstd2/Cldn4 结论。

## 方法与限制

只下载 SEARCH / mzTab / 期刊 xlsx / iProX 检索结果。别名覆盖人鼠 UniProt 与常用基因名。未下载 268 MB 光谱库；定量表缺席已足以说明患者水平未定量。PXD019061 为 ICI 初治、谱图计数非绝对定量。CLDN4 未检出不等于生物学上不存在，只说明这些开放加工矩阵里没有它。

脚本：`scripts/grok_pride/`。主表：`results/grok_pride/presence_table.tsv`。
