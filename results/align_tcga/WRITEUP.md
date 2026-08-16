# USER-ALIGN writeup — TACSTD2 vs immune signatures (TCGA NSCLC)

**Do the numbers go the same direction as the user?**
**Yes.** In public TCGA LUAD+LUSC primary tumours, TACSTD2 (TROP-2) is
negatively correlated with immune / cytotoxic / exhaustion signatures;
the sign survives purity control by **ABSOLUTE** and **CPE**, and stays
negative (sometimes losing *p* < 0.05) after **ESTIMATE**. Tight-junction
intersection genes CLDN1 / CLDN4 / CLDN7 / F11R / PARD3 are **positive**.
No statistic below is invented: every ρ / *p* / *n* is written by the
scripts from the downloaded open matrices.

---

## 中文

### 结论（先看这个）

用户的 bulk 发现是：TCGA LUAD+LUSC（NSCLC，n≈1031）里 TACSTD2 与免疫 /
细胞毒 / 耗竭签名的 Spearman ρ **几乎全为负**，用肿瘤纯度（ESTIMATE 或
ABSOLUTE/CPE）做偏相关后 **仍然显著**。

我们用同一方法、只使用公开矩阵复算：

| 核对项 | 用户 | 本复算 | 方向是否一致 |
|--------|------|--------|--------------|
| 队列 | TCGA LUAD+LUSC，n≈1031 | 原发瘤 `-01`，**n=1017**（LUAD 515 + LUSC 502） | 是（n 差 14，见下） |
| Spearman 符号 | 几乎全负 | 19/19 免疫特征为负；18/19 *p*<0.05（仅 CD274/PD-L1 不显著） | **是** |
| 偏相关后仍负 | 是 | ABSOLUTE 19/19 负、17/19 仍 *p*<0.05；CPE 18/19 负且 18/19 仍显著；ESTIMATE 17/19 负、12/19 仍显著 | **是（ABSOLUTE/CPE 更贴用户表述；ESTIMATE 更保守）** |
| TJ 交集基因 | CLDN1/4/7、F11R、PARD3 | 合并队列 5/5 为正，偏相关后仍正 | **是** |

n=1017 vs ≈1031：Xena HiSeqV2 原发瘤 `-01` 为 515+502。差约 14 例，最可能
来自 GDC STAR 多出来的样本、或把少数非 `-01` / 重复 vial 算进去。我们没有
补造样本。

### 方法（对齐用户）

1. 表达：UCSC Xena `HiSeqV2`（RSEM，log2(norm+1)），只留 barcode 以 `-01`
   结尾的原发瘤。
2. 纯度：Aran *Nat Commun* 2015 Supplementary Data 1 的 **ESTIMATE、
   ABSOLUTE、CPE**（CPE 是四法共识，不是第三种算法）。
3. 签名：均值 z-score。CYT = GZMA+PRF1（Rooney 2015）；CD8 = CD8A+CD8B；
   细胞毒 9 基因；耗竭 7 基因（PDCD1/CTLA4/LAG3/HAVCR2/TIGIT/BTLA/VSIR）；
   Ayers 2017 的 IFNγ-6 与 18 基因 T-cell-inflamed GEP。基因表见
   `tcga_signature_genes.txt`。
4. 统计：Spearman；一阶偏 Spearman = 秩变量上的偏 Pearson，df = n−3。
   免疫特征内做 BH-FDR。**没有一条数字是手填的。**

注意：ESTIMATE 纯度本身来自免疫/基质表达，用它去偏「TACSTD2 vs 免疫签名」
会部分共线、偏保守。ABSOLUTE 来自拷贝数，与 RNA 签名独立，更适合回答
「是不是纯度混杂」。

### TCGA 合并 NSCLC：签名（用户方法的核心表）

n_Spearman = 1017；ESTIMATE n=1012；ABSOLUTE n=682；CPE n=1016。

| 签名 | Spearman ρ | *p* | 偏ρ ESTIMATE | *p* | 偏ρ ABSOLUTE | *p* | 偏ρ CPE | *p* |
|------|-----------:|----:|-------------:|----:|-------------:|----:|--------:|----:|
| Cytolytic_CYT | −0.1325 | 2.2×10⁻⁵ | −0.0768 | 0.015 | −0.1331 | 5.0×10⁻⁴ | −0.1279 | 4.3×10⁻⁵ |
| CD8_Tcell | −0.1801 | 7.3×10⁻⁹ | −0.1397 | 8.3×10⁻⁶ | −0.1677 | 1.1×10⁻⁵ | −0.1809 | 6.5×10⁻⁹ |
| Cytotoxic_effector | −0.1410 | 6.4×10⁻⁶ | −0.0869 | 0.0057 | −0.1422 | 2.0×10⁻⁴ | −0.1398 | 7.8×10⁻⁶ |
| Exhaustion | −0.1234 | 8.0×10⁻⁵ | −0.0469 | 0.14 | −0.1449 | 1.5×10⁻⁴ | −0.1300 | 3.3×10⁻⁵ |
| IFNgamma_6gene | −0.1211 | 1.1×10⁻⁴ | −0.0592 | 0.060 | −0.1190 | 0.0019 | −0.1141 | 2.7×10⁻⁴ |
| Tcell_inflamed_GEP | −0.1205 | 1.2×10⁻⁴ | −0.0461 | 0.14 | −0.1264 | 9.5×10⁻⁴ | −0.1211 | 1.1×10⁻⁴ |

同一批同时有 ESTIMATE **和** ABSOLUTE 的样本（n=682）上，六条签名的
ABSOLUTE 偏ρ 全部为负且 *p*<0.05；ESTIMATE 偏ρ 仍全负，但只有 CD8
（−0.123，*p*=0.0013）和 Cytotoxic（−0.081，*p*=0.035）过 0.05
（CYT *p*=0.058）。见 `tcga_estimate_vs_absolute_completen.csv`。

免疫单基因 13 个 + 6 个签名 = 19 个特征的符号计数
（`tcga_direction_tally.csv`）：

| 队列 | Spearman 负 / 负且 *p*<0.05 | ESTIMATE | ABSOLUTE | CPE |
|------|------------------------------|----------|----------|-----|
| NSCLC 合并 | 19 / 18 | 17 / 12 | 19 / 17 | 18 / 18 |
| LUAD | 17 / 10 | 17 / 12 | 16 / 9 | 17 / 10 |
| LUSC | 19 / 19 | 19 / 13 | 19 / 19 | 19 / 19 |

LUSC 更强；LUAD 的耗竭 / IFNγ / GEP 在未校正时就不显著。合并信号不是
LUAD 单独能带起来的。CD274（PD-L1）在合并队列与 TACSTD2 基本无关
（ρ=−0.014，*p*=0.65）。

### TACSTD2 自己和纯度

| 队列 | vs ESTIMATE | vs ABSOLUTE | vs CPE |
|------|-------------|-------------|--------|
| NSCLC 合并 | +0.122（*p*=1.0×10⁻⁴，n=1012） | −0.014（*p*=0.72，n=682） | +0.053（*p*=0.093，n=1016） |
| LUAD | −0.006（NS） | +0.021（NS） | +0.021（NS） |
| LUSC | +0.120（*p*=0.007） | −0.060（NS） | +0.009（NS） |

TACSTD2 **并不**随 ABSOLUTE 纯度升高。因此「TACSTD2 高 = 瘤细胞多 = 免疫
基因被稀释」解释不了 ABSOLUTE 偏相关之后仍然为负的结果。ESTIMATE 上的
弱正相关，正是 ESTIMATE 吃进了免疫表达、偏相关会被削弱的原因。

### TJ 交集基因（CLDN1/4/7、F11R、PARD3）

合并 NSCLC：CLDN1 +0.436、CLDN4 +0.251、F11R +0.296、PARD3 +0.164、
CLDN7 +0.105（均 *p*<0.001），三种纯度偏相关后符号不变。
LUAD 里 CLDN7、PARD3 单独不显著；LUSC 五个都显著为正。

### OncoSG LUAD 2020（公开表达 + PURITY）

cBioPortal `luad_oncosg_2020`，有 mRNA 的 n=169，全部带临床 `PURITY`。
六条签名 Spearman ρ 在 −0.32 到 −0.43，偏相关后仍为 −0.26 到 −0.37，
全部 *p*<0.001。方向与 TCGA 相同，|ρ| 更大。TJ：CLDN1/4/7 为正；
F11R、PARD3 为正但不显著。

### Durvalumab 公开 RNA（用户引用的 ρ=−0.65 / −0.46）

PACIFIC / MYSTIC / COAST / NeoCOAST / OAK-POPLAR 都没有可下载的基因×样本
矩阵。找到的公开 durvalumab NSCLC bulk RNA 只有同一试验的两个 GEO 文件
（Altorki *Nat Commun* 2024）：

- **GSE253564 治疗前 FPKM，n=32**：Exhaustion ρ=**−0.642**，
  Cytolytic_CYT ρ=**−0.461**。这是公开数据里最接近用户 −0.65 / −0.46
  的一对。ESTIMATE 偏相关后耗竭仍为负（−0.357，*p*=0.049），CYT 不再显著。
- **GSE248378 治疗后 FPKM，n=29**：同向、|ρ| 更大（CYT −0.81，CD8 −0.76，
  GEP −0.63，IFNγ-6 −0.51）；多数在 ESTIMATE 之后仍显著。该矩阵没有 IFNG
  基因。

不能在没有用户方法细节的前提下断言「这就是他们用的队列/签名」，但
**方向一致，治疗前矩阵的量级也对上了**。详见 `durva_hunt.md`。

### 文件

- `tcga_tacstd2_correlations.csv` — 全特征，三种纯度
- `tcga_signatures_estimate_absolute.csv` — 签名精简表
- `tcga_estimate_vs_absolute_completen.csv` — 同时有 ESTIMATE+ABSOLUTE 的同一 n
- `tcga_direction_tally.csv` / `tcga_tacstd2_vs_purity.csv`
- `oncosg_tacstd2_correlations.csv`
- `durva_tacstd2_correlations.csv` / `durva_hunt.md`
- 脚本：`scripts/align_tcga/`

---

## English

### Answer first

The user’s bulk claim is: in TCGA LUAD+LUSC (NSCLC, n≈1031), Spearman ρ of
TACSTD2 vs immune / cytotoxic / exhaustion signatures is **almost all
negative**, and **remains significant after tumour-purity partial
correlation** (ESTIMATE or ABSOLUTE/CPE).

Recomputed with that method on public open matrices only:

| Check | User | This run | Same direction? |
|-------|------|----------|-----------------|
| Cohort | TCGA LUAD+LUSC, n≈1031 | Primary `-01` tumours, **n=1017** (LUAD 515 + LUSC 502) | Yes (n differs by 14; see below) |
| Spearman sign | almost all negative | 19/19 immune features negative; 18/19 *p*<0.05 (only CD274/PD-L1 NS) | **Yes** |
| Still negative after purity | yes | ABSOLUTE: 19/19 neg, 17/19 still *p*<0.05; CPE: 18/19 neg and 18/19 still sig; ESTIMATE: 17/19 neg, 12/19 still sig | **Yes for ABSOLUTE/CPE; ESTIMATE is more conservative** |
| TJ intersection | CLDN1/4/7, F11R, PARD3 | 5/5 positive in the pool; stay positive after all three purity methods | **Yes** |

n=1017 vs ≈1031: Xena HiSeqV2 primaries are 515+502. The ~14-sample gap is
most likely extra GDC STAR cases or inclusion of non-`-01` / extra vials.
We did not impute samples.

### Method (matched to the user)

1. Expression: UCSC Xena `HiSeqV2` (RSEM, log2(norm+1)), primary tumours
   only (barcode suffix `-01`).
2. Purity: Aran *Nat Commun* 2015 Supplementary Data 1 — **ESTIMATE,
   ABSOLUTE, CPE** (CPE = consensus of four methods).
3. Signatures: mean of per-gene z-scores. CYT = GZMA+PRF1 (Rooney 2015);
   CD8 = CD8A+CD8B; 9-gene cytotoxic; 7-gene exhaustion
   (PDCD1/CTLA4/LAG3/HAVCR2/TIGIT/BTLA/VSIR); Ayers 2017 IFNγ-6 and
   18-gene T-cell-inflamed GEP. Membership: `tcga_signature_genes.txt`.
4. Stats: Spearman; first-order partial Spearman = partial Pearson on
   ranks, df = n−3. BH-FDR within immune features.
   **No number is hand-entered.**

ESTIMATE purity is itself a function of immune/stromal expression, so
partialling it out of TACSTD2-vs-immune is partly circular and
conservative. ABSOLUTE is copy-number-based and independent of the RNA
signatures.

### Pooled TCGA NSCLC — signatures (core table)

n_Spearman = 1017; ESTIMATE n=1012; ABSOLUTE n=682; CPE n=1016.

| Signature | Spearman ρ | *p* | partial ρ ESTIMATE | *p* | partial ρ ABSOLUTE | *p* | partial ρ CPE | *p* |
|-----------|-----------:|----:|-------------------:|----:|-------------------:|----:|--------------:|----:|
| Cytolytic_CYT | −0.1325 | 2.2×10⁻⁵ | −0.0768 | 0.015 | −0.1331 | 5.0×10⁻⁴ | −0.1279 | 4.3×10⁻⁵ |
| CD8_Tcell | −0.1801 | 7.3×10⁻⁹ | −0.1397 | 8.3×10⁻⁶ | −0.1677 | 1.1×10⁻⁵ | −0.1809 | 6.5×10⁻⁹ |
| Cytotoxic_effector | −0.1410 | 6.4×10⁻⁶ | −0.0869 | 0.0057 | −0.1422 | 2.0×10⁻⁴ | −0.1398 | 7.8×10⁻⁶ |
| Exhaustion | −0.1234 | 8.0×10⁻⁵ | −0.0469 | 0.14 | −0.1449 | 1.5×10⁻⁴ | −0.1300 | 3.3×10⁻⁵ |
| IFNgamma_6gene | −0.1211 | 1.1×10⁻⁴ | −0.0592 | 0.060 | −0.1190 | 0.0019 | −0.1141 | 2.7×10⁻⁴ |
| Tcell_inflamed_GEP | −0.1205 | 1.2×10⁻⁴ | −0.0461 | 0.14 | −0.1264 | 9.5×10⁻⁴ | −0.1211 | 1.1×10⁻⁴ |

On the **same n=682** samples that have both ESTIMATE and ABSOLUTE, all
six ABSOLUTE partial ρ stay negative and *p*<0.05; ESTIMATE partial ρ
stay negative but only CD8 (−0.123, *p*=0.0013) and Cytotoxic (−0.081,
*p*=0.035) remain *p*<0.05 (CYT *p*=0.058). File:
`tcga_estimate_vs_absolute_completen.csv`.

Sign counts over 19 immune features (`tcga_direction_tally.csv`):

| Cohort | Spearman neg / neg & *p*<0.05 | ESTIMATE | ABSOLUTE | CPE |
|--------|-------------------------------|----------|----------|-----|
| NSCLC pooled | 19 / 18 | 17 / 12 | 19 / 17 | 18 / 18 |
| LUAD | 17 / 10 | 17 / 12 | 16 / 9 | 17 / 10 |
| LUSC | 19 / 19 | 19 / 13 | 19 / 19 | 19 / 19 |

LUSC is the stronger histology. LUAD exhaustion / IFNγ / GEP are already
non-significant before purity control. CD274 (PD-L1) is uncorrelated with
TACSTD2 in the pool (ρ=−0.014, *p*=0.65).

### TACSTD2 vs purity itself

| Cohort | vs ESTIMATE | vs ABSOLUTE | vs CPE |
|--------|-------------|-------------|--------|
| NSCLC pooled | +0.122 (*p*=1.0×10⁻⁴, n=1012) | −0.014 (*p*=0.72, n=682) | +0.053 (*p*=0.093, n=1016) |
| LUAD | −0.006 (NS) | +0.021 (NS) | +0.021 (NS) |
| LUSC | +0.120 (*p*=0.007) | −0.060 (NS) | +0.009 (NS) |

TACSTD2 does **not** track ABSOLUTE purity. A “high TACSTD2 = more tumour
cells = diluted immune RNA” story therefore does not explain the
ABSOLUTE-adjusted negatives. The weak positive vs ESTIMATE is expected
if ESTIMATE has already swallowed immune expression.

### TJ intersection genes (CLDN1/4/7, F11R, PARD3)

Pooled NSCLC: CLDN1 +0.436, CLDN4 +0.251, F11R +0.296, PARD3 +0.164,
CLDN7 +0.105 (all *p*<0.001); signs unchanged after ESTIMATE / ABSOLUTE /
CPE. In LUAD alone, CLDN7 and PARD3 are NS; in LUSC all five are
significantly positive.

### OncoSG LUAD 2020 (public expression + PURITY)

cBioPortal `luad_oncosg_2020`, n=169 with mRNA; all have clinical
`PURITY`. The six signatures have Spearman ρ −0.32 to −0.43; after
PURITY, −0.26 to −0.37; all *p*<0.001. Same direction as TCGA, larger
|ρ|. TJ: CLDN1/4/7 positive; F11R and PARD3 positive but NS.

### Durvalumab public RNA (cited ρ=−0.65 / −0.46)

PACIFIC / MYSTIC / COAST / NeoCOAST / OAK-POPLAR have no downloadable
gene × sample matrix. The only open durvalumab NSCLC bulk-RNA matrices
are two GEO files from the same trial (Altorki *Nat Commun* 2024):

- **GSE253564 pre-treatment FPKM, n=32:** Exhaustion ρ=**−0.642**,
  Cytolytic_CYT ρ=**−0.461**. Closest public match to the cited
  −0.65 / −0.46. After ESTIMATE, exhaustion stays negative (−0.357,
  *p*=0.049); CYT does not.
- **GSE248378 post-treatment FPKM, n=29:** same sign, larger |ρ|
  (CYT −0.81, CD8 −0.76, GEP −0.63, IFNγ-6 −0.51); most remain
  significant after ESTIMATE. IFNG is missing from this matrix.

We cannot assert these are the user’s exact cohort/signature pair without
their methods, but **the direction matches and the pre-treatment
magnitudes line up**. Details: `durva_hunt.md`.

### Files

- `tcga_tacstd2_correlations.csv` — all features, three purity methods
- `tcga_signatures_estimate_absolute.csv` — signature-only table
- `tcga_estimate_vs_absolute_completen.csv` — same-n ESTIMATE vs ABSOLUTE
- `tcga_direction_tally.csv` / `tcga_tacstd2_vs_purity.csv`
- `oncosg_tacstd2_correlations.csv`
- `durva_tacstd2_correlations.csv` / `durva_hunt.md`
- Scripts: `scripts/align_tcga/`
