# TISMO 肺模型：Tacstd2 / Cldn4，基线 vs ICB
# TISMO lung-only: Tacstd2 / Cldn4, baseline vs ICB paired

> 全部数字来自 TISMO 官方接口与官方矩阵（2026-08-16 下载），脚本在 `scripts/hunt_tismo_lung/`，表在 `results/hunt_tismo_lung/`。未编造登录号或 p 值。
> All numbers are from the official TISMO API/matrices (retrieved 2026-08-16). Scripts: `scripts/hunt_tismo_lung/`. Tables: `results/hunt_tismo_lung/`. No accessions or p-values were invented.

---

## 1. 诚实结论 / Honest verdict

| 问题 / Question | 结论 / Answer |
|---|---|
| 全癌种 Tacstd2 49/64、p=5.8e-5 能否复现？ / Does all-cancer Tacstd2 49/64, p=5.8e-5 reproduce? | **能 / Yes.** 49/64 组均值上升；Wilcoxon 符号秩 **W=439, p=5.84×10⁻⁵**。这就是用户的 5.8e-5。二项检验 49/64 vs 0.5 的双侧 p 是 **2.44×10⁻⁵**，**不是** 5.8e-5。 / 49/64 group-means up; Wilcoxon signed-rank **W=439, p=5.84e-5** (this is the user’s 5.8e-5). The two-sided binomial p is **2.44e-5**, not 5.8e-5. |
| 肺模型能否检验 49/64？ / Can lung models test 49/64? | **不能 / No.** TISMO 里带 ICB 配对的肺模型只有 **2 个官方比较组**，都是 **LLC / GSE155972**。2/2「上升」不是对 49/64 的检验。 / Only **2 official ICB groups**, both **LLC / GSE155972**. 2/2 “up” is not a test of 49/64. |
| 肺里 Tacstd2 是否被 ICB 诱导？ / Is Tacstd2 induced by ICB in lung? | **证据不足，且被一个离群点驱动。** TISMO 自己的 DESeq2：WT p=0.47，Setdb1-KO p=0.97。KO 组均值 0.28→1.00 几乎全来自 **SRX8918393 = 5.41**；去掉后 Welch p=0.94。WT 组在表达地板上 MWU p=0.045、Welch p=0.15。 / **Not supported, and driven by one outlier.** TISMO DESeq2: WT p=0.47, KO p=0.97. The KO mean 0.28→1.00 is almost entirely **SRX8918393 = 5.41**; leave-one-out Welch p=0.94. WT sits on the floor (MWU 0.045, Welch 0.15). |
| Cldn4 是否重复 Tacstd2 的 49/64？ / Does Cldn4 repeat the Tacstd2 49/64 pattern? | **否 / No.** 全癌种 **34/65 up**，Wilcoxon **p=0.12**。肺 2/2 组均值上升，但样本水平不显著（Welch p=0.17 / 0.48），同样被 SRX8918393（Cldn4=3.23）拉高。 / All-cancer **34/65 up**, Wilcoxon **p=0.12**. Lung 2/2 group-means up, sample-level NS, same outlier. |
| TISMO 还有没有更多肺 ICB？ / Are there more TISMO lung ICB models? | **没有。** 肺细胞系只有 LLC、CMT-167、MLE12。CMT-167 仅 GSE100412 原位、无 ICB（n=3）。MLE12 仅体外未处理。KP / 344SQ 等不在 TISMO。 / **No.** Only LLC, CMT-167, MLE12. CMT-167 is orthotopic GSE100412, no ICB (n=3). MLE12 is untreated vitro only. KP / 344SQ are not in TISMO. |

**一句话 / One line:** 全癌种 Tacstd2 49/64 是真的（作为组均值 Wilcoxon，不是二项）；**它不是肺结果。** 肺 ICB 在 TISMO 里就是 LLC 皮下瘤，Tacstd2/Cldn4 都在地板上。

---

## 2. 方法 / Methods

- **数据 / Data:** `POST https://tismo.pku-genomics.org/rtismo/gene/downVivoExprn`（`tumorList=["All"]`，`type=3` CSV）。注释与全矩阵来自 TISMO Data Download 的 Aliyun 份额（`TISMO_vivosample_annotations.csv` 1,518 样本；表达 RDS 21,729×1,518）。
- **切勿** 在服务器端把 `tumorList` 设成 `["LLC"]`：会静默丢掉 **SRX8918393**。
- **配对定义（对齐 49/64）/ Pairing (matches A4):** 每个官方 `cell_line` 组（去掉 `(n=N)`）内，`Baseline==1` 的均值 vs `Baseline==0` 的均值。R 与 NR 合并。方向按均值，不按中位数（中位数方向是 40/64，对不上用户数字）。
- **检验 / Tests:** 组水平 Wilcoxon 符号秩 + 配对 t + 二项/符号检验；肺样本水平 Welch + Mann-Whitney；TISMO 自带 DESeq2 Wald p 原样报告。
- **值 / Values:** 官网用于作图的 ComBat 校正 log2(TPM+1)。TISMO 的组 p 值来自 counts 上的 DESeq2，不是这些作图值。
- **对照基因 / Controls:** Actb, Epcam, Krt8, Cldn3, Cldn7, Cd8a, Ifng, Gzmb, Cd274。

---

## 3. 全癌种 vs 肺：组水平方向 / All-cancer vs lung group-level direction

| 基因 / Gene | 范围 / Scope | up / n | 均值基线→ICB / mean base→ICB | Wilcoxon p | 配对 t p | 二项 p |
|---|---|---|---|---|---|---|
| **Tacstd2** | 全癌种 / all | **49/64** | 1.051 → 1.307 | **5.84e-5** | 0.0036 | 2.44e-5 |
| **Tacstd2** | 肺 / lung | 2/2 | 0.219 → 0.647 | 0.50 | 0.38 | 0.50 |
| **Cldn4** | 全癌种 / all | **34/65** | 1.985 → 2.098 | **0.123** | 0.16 | 0.61 |
| **Cldn4** | 肺 / lung | 2/2 | 0.198 → 0.405 | 0.50 | 0.35 | 0.50 |
| Actb | 全癌种 / all | 41/65 | 11.64 → 11.83 | 0.0023 | 6.4e-4 | 0.033 |
| Cd8a | 全癌种 / all | 54/65 | 1.54 → 2.15 | 1.2e-7 | 3.6e-7 | 6.0e-8 |
| Ifng | 全癌种 / all | 50/65 | 1.09 → 1.70 | 2.3e-7 | 1.2e-7 | 7.1e-6 |
| Gzmb | 全癌种 / all | 52/65 | 2.76 → 3.64 | 9.1e-9 | 2.5e-9 | 1.2e-6 |
| Cd274 | 全癌种 / all | 52/65 | 3.09 → 3.67 | 2.9e-7 | 8.2e-8 | 1.2e-6 |

**解读 / Read this as:**

1. 用户的 Tacstd2 49/64 + 5.8e-5 **复现了**，单位是 **64 个 TISMO ICB 比较组**（约 17 条细胞系 × 22 个研究），不是 64 个独立肿瘤模型。B16 单独占 11 组。存在伪重复。
2. **Cldn4 没有这个方向一致性。** 34/65 与抛硬币无异。
3. Actb 也「上升」41/65（p≈0.002），说明作图值上有轻微全局偏移；但 Cd8a/Ifng/Gzmb/Cd274 的效应强得多，流水线本身是通的。
4. 肺 n=2 的 Wilcoxon/t **只是诊断**，不能当证据。

图：`figures/all_cancer_delta_waterfall.png`（肺组标红）、`figures/sign_counts.png`。

---

## 4. 唯一的肺 ICB 研究：GSE155972 LLC / The only lung ICB study

Griffin et al., LLC ± Setdb1-KO，抗 PD-1 + 抗 CTLA-4。皮下 flank，不是原位肺。TISMO 把 WT+ICB 标成 NR、Setdb1-KO+ICB 标成 R——**应答与基因型完全混淆**。

### 4.1 组均值（对齐 49/64 的定义）/ Group means

| 基因 | 组 | n base/ICB | mean base | mean ICB | Δ | TISMO DESeq2 p |
|---|---|---|---|---|---|---|
| Tacstd2 | LLC WT | 10/6 | 0.157 | 0.290 | +0.133 | **0.467** |
| Tacstd2 | LLC Setdb1-KO | 10/7 | 0.280 | 1.005 | +0.725 | **0.970** |
| Cldn4 | LLC WT | 10/6 | 0.128 | 0.209 | +0.081 | **0.627** |
| Cldn4 | LLC Setdb1-KO | 10/7 | 0.267 | 0.601 | +0.334 | **0.618** |
| Cd8a | LLC WT | 10/6 | 1.83 | 4.16 | +2.33 | **8.0e-4** |
| Cd8a | LLC Setdb1-KO | 10/7 | 3.03 | 4.49 | +1.46 | **1.5e-3** |
| Cd274 | LLC WT | 10/6 | 3.75 | 5.30 | +1.55 | **6.7e-4** |
| Cd274 | LLC Setdb1-KO | 10/7 | 4.56 | 5.81 | +1.25 | **9e-6** |

ICB 在这个模型里**确实**提高了 Cd8a / Cd274（TISMO DESeq2 显著）。Tacstd2 / Cldn4 **不显著**。

### 4.2 样本水平 + 离群点 / Sample-level + outlier

| 基因 | 组 | Welch p | MWU p | 去掉离群点后 Welch p |
|---|---|---|---|---|
| Tacstd2 | WT | 0.15 | **0.045** | （无离群点） |
| Tacstd2 | Setdb1-KO | 0.37 | 0.31 | **0.94**（丢掉 SRX8918393） |
| Cldn4 | WT | 0.17 | 0.17 | （无离群点） |
| Cldn4 | Setdb1-KO | 0.48 | 0.89 | 0.29（丢掉 SRX8918393） |

SRX8918393（Setdb1-KO ICB responder）Tacstd2=**5.41**，其余 6 个 KO ICB 值 ≤0.52。同一只老鼠 Cldn4=**3.23**，其余 ≤0.33。组均值「上升」是这一只样本。

图：`figures/lung_gse155972_strips.png`。

---

## 5. 本轮多挖到的肺内容 / What “hunt more” actually added

Previous TISMO agents stopped at Tacstd2 ICB CSVs. This run also did Cldn4, every TISMO lung study, vitro cytokines, and baseline ranks.

### 5.1 Cldn4 全癌种（此前没算过）/ All-cancer Cldn4 (not computed before)

34 up / 29 down / 2 tie，共 **65** 组（含 Tacstd2 缺失的 MOC22 组）。Wilcoxon p=0.12。**Cldn4 不是 Tacstd2 的 49/64 同类信号。**

### 5.2 其余 TISMO 肺 in-vivo（无 ICB 配对）/ Other TISMO lung in-vivo (no ICB pair)

71 个肺样本、9 个研究。除 GSE155972 外都不能做 baseline-vs-ICB。

- **CMT-167 原位（GSE100412, n=3）:** Tacstd2 均值 0.52，Cldn4 均值 1.73。在 TISMO 基线细胞系排名里 Tacstd2 第 32/65、Cldn4 第 29/65。有表达，但**没有 ICB 臂**。
- **LLC 基线** 在 9 个研究里都低：Tacstd2 均值 0.25（第 41/65），Cldn4 均值 0.17（第 49/65）。LLC 是间质/低分化模型，不是 Tacstd2-high 上皮模型。
- 原位 LLC（GSE100412 / GSE131271）Tacstd2 同样在地板上，部分样本精确为 0。

### 5.3 体外细胞因子 / In-vitro cytokines

Manguso LLC（RTM28723893）亲本 + clone7/8/9：Tacstd2 与 Cldn4 在 IFNβ / IFNγ / TNFα / 未处理下几乎全是 **0**。同一批 Cd274 被 IFN 强烈诱导（Δ≈+3.6 到 +4.5，Welch p<1e-3）。**不是检测失败，是 LLC 不表达这两个基因。**

MLE12（GSE103548，仅未处理 n=4）：Tacstd2 均值 **5.88**（TISMO 肺里唯一的 Tacstd2-high 模型），Cldn4 ≈0。无细胞因子、无 ICB、无 in-vivo。

### 5.4 免疫相关 / Immune correlation (descriptive)

全肺 71 样本：Tacstd2 vs CD8 T (mMCP-counter) Spearman ρ=0.12，p=0.32。Cldn4 ρ=0.26，p=0.028（未校正；跨研究混合，不能当生物学结论）。GSE155972 内两者都不显著。

---

## 6. 不能声称的东西 / What this does **not** show

- 不能说「肺也是 49/64」。n=2。
- 不能说「LLC 上 Tacstd2/Cldn4 被 ICB 诱导」。官方 DESeq2 不显著；KO 均值被一只离群样本拖动。
- 不能把 GSE155972 的 R vs NR 当成治疗应答：R=Setdb1-KO，NR=WT。
- 不能把皮下 LLC 当成原位肺癌。TISMO 里真正的原位肺（GSE100412、GSE131271）没有 ICB。
- 不能把全癌种 Tacstd2 方向一致性外推到 Cldn4。

---

## 7. 复现 / Reproduce

```bash
python3 scripts/hunt_tismo_lung/download_tismo.py          # official gene CSVs + metadata
python3 scripts/hunt_tismo_lung/download_tismo.py --rds    # optional ~250 MB Aliyun matrices
python3 scripts/hunt_tismo_lung/analyze.py
python3 scripts/hunt_tismo_lung/plot.py
```

已提交的 `notes/hunt_tismo_lung/raw/*_vivo.csv` 足够复现 49/64 与肺 ICB 样本统计。描述性全肺 / 体外表需要 Aliyun RDS（或已写出的 `extracted_*_gene_matrix.csv`）。

---

## 8. 文件 / Outputs

| 文件 | 内容 |
|---|---|
| `icb_paired_sign_summary.csv` | 全癌种 vs 肺，Tacstd2/Cldn4/对照的 up/down 与 p |
| `icb_paired_group_stats.csv` | 每个官方 ICB 组的均值、Δ、方向、DESeq2 p |
| `lung_gse155972_sample_level.csv` | 肺唯一 ICB 研究的 Welch/MWU + 离群点敏感性 |
| `lung_icb_per_sample.csv` | GSE155972 逐样本表达 |
| `lung_vivo_catalog.csv` | TISMO 全部 71 个肺 in-vivo 样本的设计 |
| `lung_all_studies_descriptive.csv` | 9 个肺研究的 Tacstd2/Cldn4 均值 |
| `baseline_cellline_rank_*.csv` | 基线细胞系排名 |
| `lung_vitro_*.csv` | LLC/MLE12 体外 |
| `verdict.json` | 机器可读结论 |
| `figures/` | waterfall / strip / sign-count |
