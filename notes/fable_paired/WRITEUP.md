# TACSTD2 (TROP2) & CLDN4 dynamics vs ICI response in lung cancer
# 肺癌免疫检查点治疗中 TACSTD2 (TROP2) 与 CLDN4 的动态变化与疗效关系

*Parallel slice — all artifacts live under `notes/fable_paired/`, `scripts/fable_paired/`,
`results/fable_paired/`. Raw GEO downloads stay outside the repo; committed processed data ≈ 3 MB.*

---

## English

### 1. Question
Do the epithelial/tumor genes **TACSTD2 (TROP2)** and **CLDN4 (Claudin-4)** change across the
**pre → on → post** timeline of immune-checkpoint inhibitor (ICI) therapy in lung cancer, and does
that behaviour track with **response**? These two genes are of translational interest because TROP2
is the target of the antibody–drug conjugate datopotamab deruxtecan and CLDN4 marks a claudin-low/
epithelial program, so their behaviour under ICI is relevant to combination/sequencing strategies.

### 2. Public data (paired pre/on/post ICI, lung + orthogonal)
| Timepoint axis | Cohort | What it contributes |
|---|---|---|
| **PRE→POST (same patient)** | **GSE248249** (NSCLC, 13 pairs, PD-(L)1, post = acquired resistance) | **A7 analog** — only public same-patient lung-tumor ICI RNA |
| **PRE** | GSE207422 bulk (NSCLC, n=24), GSE126044 (NSCLC, n=16), GSE135222 (NSCLC, n=27) | baseline TACSTD2/CLDN4 vs response, three independent cohorts |
| **ON** | GSE91061 (Riaz, melanoma, 43 same-patient Pre/On pairs) | within-patient pre→on change (orthogonal, not lung) |
| **POST** | GSE207422 scRNA (NSCLC neoadjuvant, 3 pre + 12 post, **0 pairs**) | epithelial TACSTD2/CLDN4 after neoadjuvant anti-PD-1+chemo |
| **Mouse lung ICB** | GSE246922 (KP, LLC1) | TISMO stand-in (TISMO UI not downloadable) |

Because TACSTD2/CLDN4 are epithelial genes, the single-cell analysis is restricted to the
**epithelial compartment** (`EPCAM+/PTPRC-`) and summarised per sample (pseudobulk). Full methods:
[`data_provenance.md`](./data_provenance.md).

### 3. Results

**A7 analog — does TACSTD2 rise after ICI in public paired lung RNA? No. (Fig. 5)**

Claim A7 (TISMO + Zhejiang IHC H-score 94→121) predicts a **rise**. The hunt
([`hunt_paired_ici_lung.md`](./hunt_paired_ici_lung.md)) found one public same-patient
lung-tumor ICI transcriptome: **GSE248249** (Memon et al., Cancer Cell 2024; 13 pairs;
post = acquired resistance, not an unselected on-treatment biopsy).

| Subset | Gene | n | median Δ log2 | ↑ / ↓ | p | vs A7 |
|---|---|---|---|---|---|---|
| all pairs | TACSTD2 | 13 | **−0.16** | 6 / 7 | 0.95 | **flat** |
| all pairs | CLDN4 | 13 | **−0.89** | 2 / 11 | **0.017** | falls |
| same anatomic site | TACSTD2 | 4 | −0.02 | 2 / 2 | 1.00 | flat |
| same anatomic site | CLDN4 | 4 | −0.74 | 0 / 4 | 0.13 | down (n small) |
| **lung → lung** | TACSTD2 | **1** | +0.60 | 1 / 0 | n/a | one pair only |
| lung → lung | CLDN4 | 1 | −0.91 | 0 / 1 | n/a | one pair only |

**12 of 13 pairs are not lung-to-lung** (mets: LN, adrenal, liver, bone, brain, spine).
The single lung-to-lung pair (Patient 01) rises +0.60 log2 — descriptive, not a test.
Same-site n=4 is a coin flip for TACSTD2. **Verdict still NOT SUPPORTED.**

Mouse lung ICB lines (GSE246922) also fail to show a Tacstd2 rise (KP ICB-resistant vs
parental log2FC −0.20, p=0.21). TISMO’s Shiny portal has no programmatic matrix dump.
No public paired lung TROP2 IHC with reusable numbers (Inoue 2025: 5 ICI patients, no
TROP2 change). Zhejiang IHC was **not** used.

**PRE — baseline vs response (Fig. 2).** In two of three NSCLC cohorts, baseline epithelial
TACSTD2 and CLDN4 were **lower in responders** (i.e. higher in non-responders), though never
significant:

| Cohort | Gene | log2FC (resp−nonresp) | AUC | p |
|---|---|---|---|---|
| GSE207422 bulk | TACSTD2 | −0.71 | 0.38 | 0.34 |
| GSE207422 bulk | CLDN4 | −0.23 | 0.36 | 0.26 |
| GSE126044 | TACSTD2 | −1.14 | 0.36 | 0.44 |
| GSE126044 | CLDN4 | −1.27 | 0.24 | 0.11 |
| GSE135222 | TACSTD2 | +0.09 | 0.43 | 0.61 |
| GSE135222 | CLDN4 | +0.20 | 0.56 | 0.69 |

GSE135222 (durable benefit proxied by PFS≥6 mo) showed essentially no separation, so the baseline
signal is a weak, non-significant trend rather than a robust biomarker.

**ON — within-patient pre→on change (Fig. 3).** In the paired melanoma cohort, **responders
upregulated TACSTD2 on-treatment** (median Δlog2 = +0.46) whereas non-responders did not (−0.06).
The responder-vs-non-responder difference in the on-treatment change was the **only nominally
significant result** (Mann-Whitney p = 0.019, AUC = 0.78). CLDN4 moved in the same direction but
weakly (Δ +0.21 vs +0.02, p = 0.40). Across all patients the paired change was not significant
(TACSTD2 Wilcoxon p = 0.18; CLDN4 p = 0.66), i.e. the effect is response-specific, not uniform.

**POST — epithelial cells after neoadjuvant therapy (Fig. 1).** In resected tumors, epithelial
TACSTD2 was **lower in responders (MPR) than non-responders (NMPR)** (log2FC −0.34, AUC 0.22,
p = 0.15) and slightly lower post- vs pre-treatment (log2FC −0.26, p = 0.29). CLDN4 was flat
(all p > 0.8). This is consistent with residual tumors of non-responders retaining a
TROP2-high malignant epithelium.

**Synthesis (Fig. 4).** A directionally coherent picture emerges: **high epithelial TACSTD2/CLDN4
associates with poorer ICI response** at baseline and in residual (post) tumors, while responders
show a distinctive **on-treatment induction of TACSTD2 from a lower baseline**. TACSTD2 carries a
stronger and more consistent signal than CLDN4.

### 4. Verification
`05_verify.py` runs biological positive controls and passed **7/7**: TACSTD2/CLDN4 are ~60× higher
in epithelial than immune cells; EPCAM/PTPRC compartment markers are specific; TACSTD2 and CLDN4
co-express across bulk tumors (Spearman r = 0.86, p = 8e-8); group sizes match the metadata; the
pipeline is deterministic. See `results/fable_paired/tables/verification_report.csv`.

### 5. Statistical honesty & limitations
- **A7 direction is not reproduced.** GSE248249 TACSTD2 is flat in all 13 pairs (p=0.95) and in
  4 same-site pairs (p=1.00). Only **one** pair is lung-to-lung. The CLDN4 drop (p=0.017) is at
  *acquired resistance* and must not be read as an on-treatment responder effect.
- The melanoma GSE91061 responder-specific TACSTD2 induction (p=0.019) does **not** survive
  Benjamini-Hochberg FDR across the original 12 contrasts (q=0.23). All response associations
  remain **hypothesis-generating**.
- Cohorts are small; regimens are heterogeneous (multiple anti-PD-1 antibodies ± chemotherapy).
- GSE207422 pre-vs-post is cross-patient (3 pre / 12 post) and post samples are confounded by chemo.
- Epithelial calling is marker-based, not CNV-refined malignant identification.
- TISMO matrices are not programmatically downloadable; GSE246922 is the closest public mouse lung ICB RNA.

### 6. Reproduce
```bash
bash scripts/fable_paired/run_all.sh
# outputs -> results/fable_paired/{tables,figures}
```

---

## 中文（Chinese）

### 1. 问题
在肺癌免疫检查点抑制剂（ICI）治疗的 **治疗前 → 治疗中 → 治疗后（pre → on → post）** 时间轴上，
上皮/肿瘤基因 **TACSTD2（TROP2）** 与 **CLDN4（Claudin-4）** 的表达是否发生动态变化，且这种变化是否与
**疗效**相关？这两个基因具有转化价值：TROP2 是抗体偶联药物 datopotamab deruxtecan 的靶点，CLDN4 标志
一种 claudin-low/上皮程序，因此它们在 ICI 下的行为对联合/序贯治疗策略有意义。

### 2. 公共配对数据（pre/on/post ICI，肺癌 + 正交对照）
| 时间轴 | 队列 | 贡献 |
|---|---|---|
| **PRE→POST（同一患者）** | **GSE248249**（NSCLC，13 对，PD-(L)1，治疗后=获得性耐药） | **A7 公共对照** — 唯一公开的同一患者肺癌肿瘤 ICI RNA |
| **PRE 治疗前** | GSE207422 bulk（NSCLC，n=24）、GSE126044（NSCLC，n=16）、GSE135222（NSCLC，n=27） | 三个独立队列的基线 TACSTD2/CLDN4 与疗效 |
| **ON 治疗中** | GSE91061（Riaz，黑色素瘤，43 对同一患者 Pre/On） | 同一患者 pre→on 变化（正交，非肺） |
| **POST 治疗后** | GSE207422 scRNA（NSCLC 新辅助，3 例治疗前 + 12 例治疗后，**0 对**） | 新辅助 anti-PD-1+化疗后上皮 TACSTD2/CLDN4 |
| **小鼠肺癌 ICB** | GSE246922（KP、LLC1） | TISMO 替代（TISMO 界面无法程序化下载） |

由于 TACSTD2/CLDN4 是上皮基因，单细胞分析仅限 **上皮细胞区室**（`EPCAM+/PTPRC-`），并按样本汇总为
伪散装（pseudobulk）。完整方法见 [`data_provenance.md`](./data_provenance.md)。

### 3. 结果

**A7 公共对照 — 公开配对肺癌 RNA 中 TACSTD2 是否在 ICI 后上升？否。（图 5）**

Claim A7（TISMO + 浙江 IHC H-score 94→121）预测 **上升**。检索
（[`hunt_paired_ici_lung.md`](./hunt_paired_ici_lung.md)）找到的唯一公开同一患者肺癌肿瘤
ICI 转录组是 **GSE248249**（Memon 等，Cancer Cell 2024；13 对；治疗后为获得性耐药，
而非未选择的治疗中活检）。

| 子集 | 基因 | n | 中位 Δ log2 | ↑ / ↓ | p | 相对 A7 |
|---|---|---|---|---|---|---|
| 全部配对 | TACSTD2 | 13 | **−0.16** | 6 / 7 | 0.95 | **持平** |
| 全部配对 | CLDN4 | 13 | **−0.89** | 2 / 11 | **0.017** | 下降 |
| 同一解剖部位 | TACSTD2 | 4 | −0.02 | 2 / 2 | 1.00 | 持平 |
| 同一解剖部位 | CLDN4 | 4 | −0.74 | 0 / 4 | 0.13 | 下降（n 小） |
| **肺→肺** | TACSTD2 | **1** | +0.60 | 1 / 0 | 无 | 仅 1 对 |
| 肺→肺 | CLDN4 | 1 | −0.91 | 0 / 1 | 无 | 仅 1 对 |

**13 对中有 12 对不是肺→肺**（转移灶：淋巴结、肾上腺、肝、骨、脑、脊柱）。唯一的肺→肺对
（患者 01）TACSTD2 上升 +0.60 log2，仅为描述、不能做检验。同一部位 n=4 对 TACSTD2 是正负各半。
**结论仍为不支持。** 未使用浙江 IHC。

**PRE — 基线与疗效（图 2）。** 三个 NSCLC 队列中有两个，基线上皮 TACSTD2 与 CLDN4 在**应答者中更低**
（即在非应答者中更高），但均未达显著：

| 队列 | 基因 | log2FC（应答−非应答） | AUC | p |
|---|---|---|---|---|
| GSE207422 bulk | TACSTD2 | −0.71 | 0.38 | 0.34 |
| GSE207422 bulk | CLDN4 | −0.23 | 0.36 | 0.26 |
| GSE126044 | TACSTD2 | −1.14 | 0.36 | 0.44 |
| GSE126044 | CLDN4 | −1.27 | 0.24 | 0.11 |
| GSE135222 | TACSTD2 | +0.09 | 0.43 | 0.61 |
| GSE135222 | CLDN4 | +0.20 | 0.56 | 0.69 |

GSE135222（以 PFS≥6 个月代表持久获益）基本无区分度，因此基线信号是一个**弱且不显著**的趋势，
而非稳健的生物标志物。

**ON — 同一患者 pre→on 变化（图 3）。** 在配对黑色素瘤队列中，**应答者在治疗中上调了 TACSTD2**
（中位 Δlog2 = +0.46），而非应答者未上调（−0.06）。治疗中变化在应答者与非应答者间的差异是**唯一
达到名义显著的结果**（Mann-Whitney p = 0.019，AUC = 0.78）。CLDN4 方向一致但较弱（Δ +0.21 vs +0.02，
p = 0.40）。若不分组，全体患者的配对变化不显著（TACSTD2 Wilcoxon p = 0.18；CLDN4 p = 0.66），
说明该效应是**疗效特异**的，而非普遍的。

**POST — 新辅助治疗后的上皮细胞（图 1）。** 在切除肿瘤中，上皮 TACSTD2 在**应答者（MPR）中低于
非应答者（NMPR）**（log2FC −0.34，AUC 0.22，p = 0.15），且治疗后略低于治疗前（log2FC −0.26，
p = 0.29）。CLDN4 无变化（所有 p > 0.8）。这与非应答者残余肿瘤保留 TROP2 高表达恶性上皮相一致。

**综合（图 4）。** 呈现出方向一致的图景：**上皮 TACSTD2/CLDN4 高表达与更差的 ICI 疗效相关**
（无论基线还是治疗后残余肿瘤）；而应答者表现出独特的**从较低基线出发、在治疗中上调 TACSTD2** 的动态。
TACSTD2 的信号比 CLDN4 更强、更一致。

### 4. 验证
`05_verify.py` 运行了生物学阳性对照，**7/7 全部通过**：TACSTD2/CLDN4 在上皮细胞中比免疫细胞高约 60 倍；
EPCAM/PTPRC 区室标记具特异性；TACSTD2 与 CLDN4 在散装肿瘤中共表达（Spearman r = 0.86，p = 8e-8）；
分组样本量与元数据一致；流程可确定性复现。见 `results/fable_paired/tables/verification_report.csv`。

### 5. 统计诚实性与局限
- **A7 方向未被复现。** GSE248249 的 TACSTD2 在全部 13 对（p=0.95）和 4 对同部位（p=1.00）均为持平。
  仅 **1** 对为肺→肺。CLDN4 下降（p=0.017）是获得性耐药时间点，不能外推到应答者的治疗中活检。
- 黑色素瘤 GSE91061 中应答者 TACSTD2 治疗中上调（p=0.019）在 12 个比较的 BH-FDR 后不再显著
  （q=0.23）。所有疗效相关结论均为**假设生成性**。
- 队列样本量小；方案异质（多种 anti-PD-1 抗体 ± 化疗）。
- GSE207422 的 pre-vs-post 为跨患者（3 例治疗前 / 12 例治疗后），且治疗后样本受化疗混杂。
- 上皮细胞判定基于标记基因，未做 CNV 精细的恶性细胞鉴定。
- TISMO 矩阵无法程序化下载；GSE246922 是最接近的公开小鼠肺癌 ICB RNA。

### 6. 复现
```bash
bash scripts/fable_paired/run_all.sh
# 输出 -> results/fable_paired/{tables,figures}
```

---

### Figures / 图
- `results/fable_paired/figures/fig5_gse248249_paired_pre_post.png` — **A7 analog**: 13 NSCLC pairs, pre vs acquired resistance
- `results/fable_paired/figures/fig6_gse246922_mouse_lung_icb.png` — mouse KP/LLC1 ICB (TISMO stand-in)
- `results/fable_paired/figures/fig1_gse207422_scrna_epithelial.png` — POST: epithelial pre/post & MPR/NMPR
- `results/fable_paired/figures/fig2_baseline_bulk_cohorts.png` — PRE: three baseline NSCLC cohorts
- `results/fable_paired/figures/fig3_gse91061_paired_pre_on.png` — ON: melanoma within-patient pre→on
- `results/fable_paired/figures/fig4_effectsize_forest.png` — effect-size synthesis (response contrasts)
