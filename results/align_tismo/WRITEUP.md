# TISMO Tacstd2: baseline vs ICB — replication

**Claim under test.** In the TISMO in vivo Gene module, Tacstd2 is higher after immune-checkpoint blockade than in the matched control arm: 64 tumor models, 49/64 up, mean ~1.05 → 1.3, paired *p* = 5.8×10⁻⁵. Lung models are to be reported separately. Counts must be honest.

**Verdict.** The headline numbers replicate exactly on the public TISMO tables. They describe **64 ICB *cohorts*** (study × cell line × condition), not 64 independent tumor models. The claimed *p* is the Wilcoxon signed-rank test on those 64 paired cohort means. Lung is two LLC cohorts from one study. Clustering-aware and genome-wide checks keep a directional Tacstd2 shift but shrink what the 5.8×10⁻⁵ figure can be taken to mean.

All scripts live in `scripts/align_tismo/`. All outputs live in `results/align_tismo/`.

---

## English

### 1. What was downloaded

`tismo.cistrome.org` 301-redirects to `tismo.pku-genomics.org`. The site is a Vue SPA. The in vivo Gene module POSTs form-urlencoded requests to two public backends (no login):

| backend | URL | used for |
|---|---|---|
| `/tismo` | `https://tismo.pku-genomics.org/tismo` | vocabularies (`getVivoTreatment`, `getVivoCohort`, `getGene`) and metadata (`vivoMeta`, `cellLineMeta`) |
| `/rtismo` | `https://tismo.pku-genomics.org/rtismo` | per-gene sample table (`/gene/downVivoExprn`) |

The six ICB regimens the module exposes are `antiCTLA4`, `antiCTLA4&antiPD1`, `antiCTLA4&antiPDL1`, `antiPD1`, `antiPDL1`, `antiPDL2`. Asking for all six plus every listed tumor model is the “All ICB / All models” query the website itself runs.

Raw payloads (unmodified) are in `results/align_tismo/data/`:

- `vivo_expression_Tacstd2.csv` — 619 rows, 605 unique SRX/sample IDs
- `cellLineMeta.json`, `vivoMeta.json` — cancer-type lookup
- `reference_genes/` — Cd274, Pdcd1, Ifng, Cd8a, Epcam, Krt8, Actb
- `null_genes.tar.gz` — 500 genes drawn with seed `20240816` from the 21,235-gene TISMO list (loose CSVs are a working cache and are gitignored)

Expression values are the ones TISMO plots: Salmon TPM, then quantile-normalised and ComBat-corrected for visualisation. They are **not** raw TPM. TISMO’s own documentation says within-cohort DESeq2 is computed on counts; the `pvalue` column on the download is that statistic.

### 2. How a “tumor model” is defined here

TISMO’s `cell_line` field is a **cohort label**, e.g. `LLC_GSE155972_Setdb1_KO_antiCTLA4&antiPD1(n=17)`. The `(n=N)` suffix is recomputed per gene, so the same cohort can appear as `(n=10)` for Actb and `(n=11)` for Tacstd2. Scripts strip it (`cohort_key`) before comparing genes.

`Baseline=1` is the study’s own control arm (isotype / no_treatment / vehicle / a non-ICB co-treatment). `Baseline=0` is the ICB-containing arm. Every treated arm in this download carries at least one of PD-1 / PD-L1 / PD-L2 / CTLA-4.

The paired unit is therefore **one TISMO cohort**: mean(ICB arm) − mean(control arm) of that study’s own design. That is the analysis the website invites and the analysis that produces 64 pairs.

### 3. Headline replication (all models)

| quantity | claim | this download |
|---|---|---|
| *n* | 64 | **64** paired cohorts |
| up / down / tie | 49 / — / — | **49 / 15 / 0** |
| mean baseline → ICB | ~1.05 → 1.3 | **1.051 → 1.317** (Δ = +0.265) |
| median baseline → ICB | not stated | 0.245 → 0.601 |
| paired Wilcoxon *p* (two-sided) | 5.8×10⁻⁵ | **5.84×10⁻⁵** |
| paired *t* *p* | — | 3.48×10⁻³ |
| sign test *p* | — | 2.44×10⁻⁵ |
| Cohen’s *d<sub>z</sub>* | — | 0.38 |

The claimed *p* is Wilcoxon, not the *t*-test. Means match after ordinary rounding (1.317 → 1.3). Figure: `figures/Tacstd2_paired_slopes.png`.

### 4. Honest denominators

The 64 is real. It is not 64 independent models.

| unit | count |
|---|---|
| TISMO ICB cohorts with both arms, Tacstd2 | **64** |
| ICB cohorts TISMO serves for other genes | **65** (the extra one is `MOC22_RU31562203_antiPD1`; absent from the Tacstd2 payload) |
| unique cell lines | **17** (TISMO lists 18 including MOC22) |
| unique studies | **22** |
| cancer types | **9** |
| unique samples | 605 (294 control, 311 ICB) |
| control samples reused in two cohorts | **14** (shared isotype arms, e.g. 4T1 GSE130472 old/young CTLA-4 vs PD-L1) |
| samples with Tacstd2 value exactly 0 | 14.2% |
| cohorts whose control mean is < 0.5 | 61% |

Calling these “64 tumor models” is TISMO’s own language. A reader who hears “64 independent syngeneic models” will over-count. The 17 cell lines and their cohort counts:

| cell line | cancer type | cohorts | studies |
|---|---|---|---|
| B16 | Melanoma | 11 | 5 |
| T11 | Mammary cancer, NOS | 9 | 1 (GSE124821) |
| 4T1 | Mammary carcinoma | 7 | 3 |
| CT26 | Colorectal carcinoma | 7 | 3 |
| KPB25L | Mammary cancer, NOS | 6 | 1 (GSE124821) |
| YTN16 | Gastric adenocarcinoma | 5 | 1 (GSE146027) |
| MC38 | Colorectal carcinoma | 3 | 2 |
| 402230 | Sarcoma | 2 | 1 |
| BNL-MEA | Hepatocellular carcinoma | 2 | 1 |
| EMT6 | Mammary carcinoma | 2 | 1 |
| LLC | **Lung carcinoma** | **2** | **1 (GSE155972)** |
| p53-2225L | Mammary cancer, NOS | 2 | 1 (GSE124821) |
| p53-2336R | Mammary cancer, NOS | 2 | 1 (GSE124821) |
| D3UV2, D4M.3A.3, E0771, YUMM1.7 | (various) | 1 each | 1 each |

GSE124821 alone contributes 19 of the 64 pairs (T11, KPB25L, p53-2225L, p53-2336R).

### 5. Lung models

TISMO’s cell-line table has three lung lines (LLC, CMT-167, MLE12). Only **LLC** has an ICB in vivo cohort in this download.

| LLC cohort | control *n* | ICB *n* | mean control | mean ICB | Δ | TISMO DESeq2 *p* |
|---|---|---|---|---|---|---|
| `LLC_GSE155972_antiCTLA4&antiPD1` | 10 | 6 | 0.157 | 0.290 | +0.133 | 0.47 |
| `LLC_GSE155972_Setdb1_KO_antiCTLA4&antiPD1` | 10 | 7 | 0.280 | 1.005 | +0.725 | 0.97 |

Both pairs go up. Mean 0.219 → 0.647. Wilcoxon *p* = 0.5, paired *t* *p* = 0.38. That is the entire lung evidence. It is consistent in direction with the pan-cancer pattern and is **not** an independent lung result.

Dropping LLC leaves 47/62 up, mean Δ = +0.260, Wilcoxon *p* = 1.3×10⁻⁴. The headline is not a lung-driven artefact.

### 6. What happens if the pairs are not treated as independent

The 64-pair Wilcoxon assumes exchangeable pairs. Cohorts that share a cell line or a GEO study are not.

| analysis | *n* | up | mean Δ | Wilcoxon *p* |
|---|---|---|---|---|
| headline (one pair per cohort) | 64 | 49 | +0.265 | 5.84×10⁻⁵ |
| one pair per cell line (average its cohorts first) | 17 | 14 | +0.293 | 3.8×10⁻³ |
| one pair per study | 22 | 18 | +0.225 | 4.2×10⁻³ |
| cluster bootstrap, resample cell lines | 17 clusters | — | +0.265 (95% CI 0.078–0.526) | *p* = 0.024 |
| cluster bootstrap, resample studies | 22 clusters | — | +0.265 (95% CI 0.088–0.534) | *p* = 0.027 |

Leave-one-cell-line-out: the weakest remaining Wilcoxon *p* is 6.4×10⁻⁴ after dropping B16 (11 cohorts). No single cell line creates the sign pattern. YTN16 (gastric, Δ = +1.41) and B16 inflate the *mean*; 4T1 actually pulls the mean *down* (cell-line Δ = −0.17).

### 7. Cancer-type slices (honest *n*)

| cancer type | cohorts | cell lines | studies | up / down | mean Δ | Wilcoxon *p* |
|---|---|---|---|---|---|---|
| Melanoma | 14 | 4 | 7 | 13 / 1 | +0.322 | 8.5×10⁻⁴ |
| Colorectal carcinoma | 10 | 2 | 5 | 8 / 2 | +0.160 | 0.049 |
| Gastric adenocarcinoma | 5 | 1 | 1 | 5 / 0 | +1.41 | 0.062 |
| Mammary cancer, NOS | 19 | 4 | 1 | 12 / 7 | +0.161 | 0.10 |
| Mammary carcinoma | 9 | 2 | 4 | 5 / 4 | −0.055 | 0.82 |
| Lung carcinoma | 2 | 1 | 1 | 2 / 0 | +0.429 | 0.50 |
| Sarcoma | 2 | 1 | 1 | 2 / 0 | +0.212 | 0.50 |
| Hepatocellular carcinoma | 2 | 1 | 1 | 1 / 1 | −0.010 | 1.00 |
| Mammary adenocarcinoma | 1 | 1 | 1 | 1 / 0 | +0.010 | — |

Melanoma is the only slice with both *n* > 10 and *p* < 0.01, and many of those B16 values sit near zero (low Tacstd2). Mammary carcinoma (4T1 + EMT6) does not go up. A “Tacstd2 rises after ICB in every cancer type” reading is not supported.

### 8. TISMO’s own within-cohort DESeq2

50 of 64 cohorts carry a precomputed DESeq2 *p*. **7 / 50** are *p* < 0.05, and all 7 are up. The other 14 have no *p* (TISMO leaves it blank when the test is not computed). So the 49/64 claim is a **meta-sign count on visualised means**, not “49 cohorts are individually significant.” Several of the 7 significant DESeq2 calls have a tiny mean-Δ on the plotted scale (e.g. B16 GSE149825 parental combo, Δ = +0.003, *p* = 1.0×10⁻⁴) because DESeq2 is run on counts, not on the quantile-normalised values.

### 9. Is 49/64 unusual among genes?

The same 64-pair Wilcoxon was run on 500 randomly chosen TISMO genes (seed 20240816) and on seven reference genes.

**Null panel (*n* = 500)**

| null quantity | value |
|---|---|
| median *n* up | 25 / ~64 |
| genes with majority up | 22.8% |
| median mean Δ | −0.003 |
| genes with mean Δ > 0 | 36.4% |
| genes with Wilcoxon *p* < 0.05 | **43.6%** |
| genes with Wilcoxon *p* ≤ Tacstd2’s 5.84×10⁻⁵ | 4.8% |
| genes with *n* up ≥ 49 | 1.4% |
| genes with mean Δ ≥ +0.265 | 1.6% |
| Tacstd2 Wilcoxon *p* × 21,235 genes | 1.24 |

Two facts sit next to each other. First, **this design is anti-conservative**: almost half of random genes reach *p* < 0.05, so 5.8×10⁻⁵ is not a genome-wide discovery *p*. Tacstd2 would not survive Bonferroni. Second, Tacstd2 is still in the **right tail** of the same null on the quantities the claim actually uses (49 ups; Δ ≈ 0.27). The pattern is not “every gene goes up after ICB.” Figure: `figures/Tacstd2_null_calibration.png`.

**Reference genes (same statistic)**

| gene | role | cohorts | up | mean Δ | Wilcoxon *p* |
|---|---|---|---|---|---|
| Cd8a | ICB pharmacodynamic | 65 | 54 | +0.595 | 1.2×10⁻⁷ |
| Cd274 (PD-L1) | ICB pharmacodynamic | 65 | 52 | +0.589 | 2.9×10⁻⁷ |
| Ifng | ICB pharmacodynamic | 63 | 48 | +0.539 | 6.5×10⁻⁷ |
| Pdcd1 (PD-1) | ICB pharmacodynamic | 64 | 49 | +0.431 | 6.8×10⁻⁶ |
| **Tacstd2** | **claim** | **64** | **49** | **+0.265** | **5.8×10⁻⁵** |
| Actb | housekeeper | 65 | 41 | +0.197 | 2.3×10⁻³ |
| Epcam | epithelial content | 65 | 43 | +0.113 | 0.017 |
| Krt8 | epithelial content | 65 | 31 | −0.078 | 0.51 |

ICB pharmacodynamic markers move more than Tacstd2, as they should if the ICB arms are biologically different. Epcam / Krt8 do not copy the Tacstd2 sign pattern, so a simple “more tumour cells after ICB” story is not enough. Actb also drifts up, which is one reason 44% of random genes look “significant.” Figure: `figures/Tacstd2_reference_genes.png`.

Cohort *n* is 63–65 across genes because TISMO silently drops a cohort when a gene has no usable rows. Tacstd2’s 64 is 64 of the 65 ICB cohorts the database can serve.

### 10. What this does and does not show

**Shown.** On the public TISMO in vivo Gene-module tables, Tacstd2’s ICB-arm mean exceeds the matched control-arm mean in 49 of 64 cohorts. The means are 1.05 and 1.32. The Wilcoxon *p* on those 64 pairs is 5.84×10⁻⁵. Those three numbers are not a misremembering of the website.

**Not shown.**

1. **64 independent models.** Honest *n* for independence is closer to 17 cell lines or 22 studies. Clustering-aware *p* is ~0.02–0.03, not 5.8×10⁻⁵.
2. **A lung result.** Lung = 2 LLC cohorts, one study, *p* = 0.5.
3. **Per-cohort differential expression.** TISMO’s own DESeq2 is significant in 7 of 50 tested cohorts.
4. **A Tacstd2-specific genome-wide hit.** Under this design 44% of random genes have *p* < 0.05; Bonferroni over the transcriptome exceeds 1.
5. **A uniform cancer-type effect.** Melanoma carries most of the well-powered sign test; 4T1 / mammary carcinoma does not rise.
6. **Causality or a treatment effect on Tacstd2 transcription.** ICB arms differ from control arms in immune content, tumour size, time point, and co-treatments (radiation, birinapant, Setdb1-KO, diet, …). The contrast is “ICB-containing arm vs that study’s control,” not a clean Tacstd2 induction assay.

**Practical reading.** The user-supplied 49/64, 1.05→1.3, *p* = 5.8×10⁻⁵ line is an accurate description of TISMO’s cohort-level Tacstd2 table. It is a weak, directionally consistent, pan-cancer sign pattern on a visualisation scale, dominated by a few cell lines (YTN16, B16, p53-2336R) and not independently demonstrated in lung.

### 11. How to rerun

```bash
# needs Python 3.12+, pandas, numpy, scipy, matplotlib
python3 scripts/align_tismo/01_download.py --null-genes 500 --workers 6
python3 scripts/align_tismo/02_replicate.py
python3 scripts/align_tismo/03_null_calibration.py
python3 scripts/align_tismo/04_figures.py
# or: bash scripts/align_tismo/run_all.sh
```

`01_download.py` hits the live TISMO API. Everything after that is offline from `results/align_tismo/data/`.

---

## 中文

### 1. 下载了什么

`tismo.cistrome.org` 现 301 到 `tismo.pku-genomics.org`。站点是 Vue 单页。体内 Gene 模块用表单 POST 打两个公开后端（无需登录）：`/tismo` 取词表和元数据，`/rtismo` 取单基因样本表。六个 ICB 方案为 anti-CTLA4 / PD-1 / PD-L1 / PD-L2 及其组合。对全部方案 + 全部瘤株的查询，就是网页上 “All ICB / All models”。

原始、未改动的表在 `results/align_tismo/data/`。Tacstd2 共 619 行、605 个独立样本。表达值是 TISMO 作图用的：Salmon TPM，再分位数标准化 + ComBat，**不是**原始 TPM。表里的 `pvalue` 是 TISMO 按 counts 算的队列内 DESeq2。

### 2. 这里的 “tumor model” 是什么

TISMO 的 `cell_line` 字段是**队列标签**，例如 `LLC_GSE155972_Setdb1_KO_antiCTLA4&antiPD1(n=17)`。末尾 `(n=N)` 按基因重算，同一队列在不同基因上 N 会变，比较前先去掉。`Baseline=1` 是该研究自己的对照臂，`Baseline=0` 是含 ICB 的处理臂。配对单位是：**同一 TISMO 队列里，ICB 臂均值 − 对照臂均值**。这正是网页在做、并且能做出 64 对的分析。

### 3. 全模型主结果（与用户数字对照）

| 量 | 用户说法 | 本次下载 |
|---|---|---|
| *n* | 64 | **64** 个成对队列 |
| 上升 / 下降 / 平 | 49 / — | **49 / 15 / 0** |
| 均值 对照 → ICB | ~1.05 → 1.3 | **1.051 → 1.317**（Δ = +0.265） |
| 配对 Wilcoxon 双侧 *p* | 5.8×10⁻⁵ | **5.84×10⁻⁵** |
| 配对 *t* *p* | — | 3.48×10⁻³ |

用户给的 *p* 是 Wilcoxon，不是 *t* 检验。均值按常规四舍五入即可对上（1.317 → 1.3）。图：`figures/Tacstd2_paired_slopes.png`。

### 4. 诚实的分母

64 这个数是真的，但**不是** 64 个独立瘤株。

| 单位 | 计数 |
|---|---|
| Tacstd2 上双臂齐全的 ICB 队列 | **64** |
| TISMO 对其他基因能给出的 ICB 队列 | **65**（多出来的是 `MOC22_RU31562203_antiPD1`，Tacstd2 载荷里没有） |
| 独立细胞系 | **17**（词表里含 MOC22 共 18） |
| 独立研究 | **22** |
| 癌种 | **9** |
| 独立样本 | 605（对照 294，ICB 311） |
| 被两个队列共用的对照样本 | **14**（例如 4T1 GSE130472 老年/年轻的共用 isotype） |

GSE124821 一家就占 19/64 对（T11、KPB25L、p53-2225L、p53-2336R）。把 64 听成 “64 个独立同源模型” 会多算。

### 5. 肺模型（单独报）

细胞系表里有三条肺系（LLC、CMT-167、MLE12）。本次 ICB 下载里**只有 LLC**。两条队列都来自 GSE155972：亲本 Δ = +0.13，Setdb1-KO Δ = +0.73；都上升；均值 0.22 → 0.65；Wilcoxon *p* = 0.5，配对 *t* *p* = 0.38。这就是全部肺证据。方向与全癌种一致，**不能**当成独立的肺结论。去掉 LLC 后仍是 47/62 上升、Wilcoxon *p* = 1.3×10⁻⁴，主结果不是肺在撑。

### 6. 若不把 64 对当成独立

| 分析 | *n* | 上升 | 均值 Δ | Wilcoxon *p* |
|---|---|---|---|---|
| 主分析（每队列一对） | 64 | 49 | +0.265 | 5.84×10⁻⁵ |
| 先按细胞系平均再配对 | 17 | 14 | +0.293 | 3.8×10⁻³ |
| 先按研究平均再配对 | 22 | 18 | +0.225 | 4.2×10⁻³ |
| 按细胞系整簇重抽样 | 17 簇 | — | +0.265（95% CI 0.078–0.526） | *p* = 0.024 |
| 按研究整簇重抽样 | 22 簇 | — | +0.265（95% CI 0.088–0.534） | *p* = 0.027 |

逐个去掉细胞系：最弱的剩余 Wilcoxon *p* 是去掉 B16（11 个队列）后的 6.4×10⁻⁴。没有哪一株单独制造符号。YTN16（胃，Δ = +1.41）和 B16 抬高均值；4T1 反而把均值往下拉（细胞系 Δ = −0.17）。

### 7. 按癌种切开（诚实 *n*）

黑色素瘤 14 队列、4 株、7 个研究，13/14 上升，*p* = 8.5×10⁻⁴，是唯一既有 *n* > 10 又有 *p* < 0.01 的切片，但不少 B16 值接近 0。乳腺瘤（4T1 + EMT6）9 队列 5 升 4 降，均值 Δ 为负。胃癌 YTN16 5/5 上升但只有 1 株 1 个研究。肺、肉瘤、肝都是 *n* = 2。**不能**读成 “每个癌种 ICB 后 Tacstd2 都升”。

### 8. TISMO 自己的队列内 DESeq2

64 个队列里 50 个带预计算 *p*。**7 / 50** 达到 *p* < 0.05，且这 7 个全是上升。49/64 是**可视化均值上的符号计数**，不是 “49 个队列各自差异显著”。有的 DESeq2 显著队列在作图尺度上 Δ 极小（例如 B16 GSE149825 亲本联合治疗 Δ = +0.003），因为检验走的是 counts，不是分位数标准化值。

### 9. 在基因背景里，49/64 特不特殊？

同一套 64 对 Wilcoxon 跑了 500 个随机基因（种子 20240816）和 7 个参照基因。

随机面板：中位上升数 25；22.8% 的基因过半数上升；**43.6% 的基因 Wilcoxon *p* < 0.05**；只有 4.8% 的 *p* 不大于 Tacstd2；*n* 上升 ≥ 49 的占 1.4%；均值 Δ ≥ +0.265 的占 1.6%；Tacstd2 的 *p* × 21,235 ≈ 1.24。

两句话要一起说。第一，**这个设计偏自由**：近一半随机基因都能 *p* < 0.05，5.8×10⁻⁵ 不能当全基因组发现 *p*，过不了 Bonferroni。第二，就用户用的两个量（49 个上升、Δ ≈ 0.27）而言，Tacstd2 仍在同一零分布的**右尾**。并不是 “ICB 后每个基因都升”。

参照：Cd8a / Cd274 / Ifng / Pdcd1 的 Δ 和显著性都强于 Tacstd2（ICB 药效标志本应如此）。Epcam / Krt8 并不复制 Tacstd2 的符号，单用 “ICB 后瘤细胞比例更高” 解释不够。Actb 也在漂高，这是 44% 随机基因 “显著” 的原因之一。不同基因队列数在 63–65 之间，因为 TISMO 在某基因没有可用行时会默默丢掉该队列。Tacstd2 的 64 是数据库能提供的 65 个 ICB 队列中的 64 个。

### 10. 能说什么，不能说什么

**能说。** 在 TISMO 公开的体内 Gene 模块表上，Tacstd2 的 ICB 臂均值在 64 个队列里有 49 个高于配对对照；均值 1.05 → 1.32；这 64 对的 Wilcoxon *p* = 5.84×10⁻⁵。这三个数字不是记错网页。

**不能说。**

1. **64 个独立模型。** 按独立性更接近 17 株或 22 个研究；考虑聚类后 *p* 约 0.02–0.03。
2. **肺的结论。** 肺 = 2 条 LLC 队列、1 个研究，*p* = 0.5。
3. **逐队列差异表达。** TISMO 自己的 DESeq2 只在 50 个已测队列里的 7 个显著。
4. **Tacstd2 的全基因组命中。** 此设计下 44% 随机基因 *p* < 0.05；乘转录组后 > 1。
5. **各癌种一致。** 有把握的符号检验主要在黑色素瘤；4T1 / 乳腺瘤不升。
6. **因果或 Tacstd2 被 ICB 诱导。** ICB 臂与对照臂在免疫成分、瘤体积、时间点、联合处理（放疗、birinapant、Setdb1-KO、饮食等）上都不干净。这是 “含 ICB 的臂 vs 该研究对照”，不是 Tacstd2 诱导实验。

**读法。** 用户给的 49/64、1.05→1.3、*p* = 5.8×10⁻⁵ 是对 TISMO 队列水平 Tacstd2 表的准确描述。它是可视化尺度上、方向一致、全癌种的弱符号模式，主要由少数细胞系（YTN16、B16、p53-2336R）带动，肺里没有独立验证。

### 11. 复现命令

```bash
python3 scripts/align_tismo/01_download.py --null-genes 500 --workers 6
python3 scripts/align_tismo/02_replicate.py
python3 scripts/align_tismo/03_null_calibration.py
python3 scripts/align_tismo/04_figures.py
```

`01_download.py` 访问线上 TISMO；其后步骤只读 `results/align_tismo/data/`。

---

## File map

| path | contents |
|---|---|
| `scripts/align_tismo/tismo_client.py` | public API client + CSV loader |
| `scripts/align_tismo/01_download.py` | fetch vocabularies, metadata, Tacstd2, references, null panel |
| `scripts/align_tismo/02_replicate.py` | paired stats, lung slice, clustering, leave-one-line-out |
| `scripts/align_tismo/03_null_calibration.py` | 500-gene null + reference genes |
| `scripts/align_tismo/04_figures.py` | four PNGs |
| `results/align_tismo/tables/claim_vs_result.csv` | claim vs observed, one row per item |
| `results/align_tismo/tables/cohort_level_Tacstd2.csv` | the 64 pairs |
| `results/align_tismo/tables/paired_summaries_Tacstd2.csv` | every subset in §3–§7 |
| `results/align_tismo/figures/` | slopes, cancer type, null, references |
