# TACSTD2-high = TLS / B / T-low in public lung RNA? Honest slice

Primary claim tested (Bessede + TROP2/ADC literature): **TACSTD2-high lung
tumours are TLS-low, B-cell-low and T-cell-low, and this survives tumour-purity
correction.** CLDN4 is the paired epithelial comparator.

**Short answer.** After DNA-based ABSOLUTE purity correction the full package
(TLS + B + CD8 T) holds in **one** large public set: **TCGA-LUSC** (n = 493).
It does **not** hold as a TLS effect in TCGA-LUAD (TLS_Cabrita even flips
positive). An independent LUAD microarray atlas (GSE72094, n = 442) supports
B / plasma / CD8 after an epithelial-score correction, but not the Cabrita TLS
signature. Every open lung ICI RNA matrix we could download is n ≤ 27 and is
underpowered; TACSTD2 does not predict PFS, DCB, RECIST or MPR in those files.
Patil / OAK / POPLAR / IMpower RNA is EGA-controlled and was **listed only**.

All numbers below are from `results/opus_tls/tables/`. No statistic was typed
by hand from memory.

---

## English

### What was tested

| Item | Choice |
|---|---|
| Genes | TACSTD2 (TROP2), CLDN4; epithelial comparators EPCAM, KRTs, CDH1, MUC1 |
| TLS / B / T scores | Cabrita 8-gene TLS; 12-chemokine TLS; Meylan imprint; B-cell; plasma; Tfh; CD8; Ayers IFN-γ |
| Purity | TCGA: ABSOLUTE (DNA) ± methylation leukocyte fraction. Others: pan-epithelial score |
| HOLDS rule (pre-specified) | n ≥ 40, partial Spearman ρ < 0, p < 0.05 |
| Open ICI | GSE135222 (n=27, PFS), GSE126044 (n=16, RECIST), GSE207422 (n=24, MPR) |
| Open atlases | TCGA-LUAD 528, TCGA-LUSC 501, GSE72094 442, GSE81089 197, GSE190265 26 |
| Not used | Patil 2022 / OAK / POPLAR / IMpower (EGA). See `DATA_ACCESS.md` |

Signatures are transcriptional TLS *surrogates*. None of these open matrices
carry a pathologist TLS call.

### Where the hypothesis holds after purity

From `where_it_holds.csv` and bootstrap CIs in `key_partial_ci.csv`
(2 000 resamples of the partial Spearman).

**1. TCGA-LUSC — holds as a package (ABSOLUTE, n = 493)**

| Signature | ρ_adj | 95% CI | p |
|---|---:|---|---:|
| TLS_Cabrita | −0.186 | −0.269, −0.100 | 3.3 × 10⁻⁵ |
| TLS_12chemokine | −0.313 | −0.393, −0.227 | 1.3 × 10⁻¹² |
| TLS_imprint | −0.243 | −0.322, −0.162 | 4.7 × 10⁻⁸ |
| B_cell | −0.223 | −0.304, −0.133 | 5.5 × 10⁻⁷ |
| Plasma_cell | −0.144 | −0.230, −0.049 | 1.3 × 10⁻³ |
| Tfh | −0.337 | −0.414, −0.257 | 1.6 × 10⁻¹⁴ |
| T_cell_CD8 | −0.286 | −0.368, −0.202 | 9.7 × 10⁻¹¹ |
| IFNg_Ayers | −0.263 | −0.345, −0.178 | 3.2 × 10⁻⁹ |

Q4 vs Q1 TACSTD2 (n = 125 vs 126): CD8 rank-biserial r = −0.37, p = 4.8 × 10⁻⁷;
B-cell r = −0.25, p = 5.1 × 10⁻⁴; TLS_Cabrita r = −0.22, p = 2.9 × 10⁻³.
Stouffer combination across purity tertiles stays significant for TLS, B and
CD8 (`quartile_contrasts.csv`).

OLS of inverse-normal TLS_Cabrita on TACSTD2 + ABSOLUTE + age/sex/stage:
β = −0.132, p = 2.3 × 10⁻³, n = 479 (`multivariable_models.csv`).

Transcriptome-wide rank after ABSOLUTE (`transcriptome_wide_null.csv`):
TACSTD2 sits at percentile 1.3% of 32 628 genes vs CD8 (rank 408), 2.1% vs
B-cell (rank 701), 2.4% vs TLS_Cabrita (rank 795). Empirical two-sided
p-values against the whole transcriptome are 0.078 / 0.134 / 0.192 — so
TACSTD2 is in the anti-immune tail but is **not** a unique genome-wide
outlier. EPCAM is *less* anti-CD8 than TACSTD2 after purity in LUSC
(EPCAM ρ_adj = −0.057 vs TACSTD2 −0.286), which argues this is not “any
epithelial gene”.

CIBERSORT (orthogonal, relative fractions): TACSTD2 vs CD8 ρ_adj = −0.166,
p = 2.2 × 10⁻⁴. Plasma CIBERSORT goes the **other** way in LUSC
(ρ_adj = +0.122, p = 6.6 × 10⁻³). Signature-plasma and CIBERSORT-plasma
therefore disagree in squamous tumours; that disagreement is reported, not
smoothed over.

**2. GSE72094 LUAD microarray — B / plasma / T hold; Cabrita TLS does not**
(epithelial-score partial, n = 442; no ABSOLUTE)

| Signature | ρ_adj | 95% CI | p | call |
|---|---:|---|---:|---|
| TLS_Cabrita | −0.019 | −0.117, +0.074 | 0.69 | NO_EVIDENCE |
| TLS_12chemokine | −0.073 | −0.168, +0.028 | 0.13 | NO_EVIDENCE |
| TLS_imprint | −0.119 | −0.210, −0.023 | 0.012 | HOLDS |
| B_cell | −0.115 | −0.210, −0.024 | 0.016 | HOLDS |
| Plasma_cell | −0.220 | −0.303, −0.132 | 2.9 × 10⁻⁶ | HOLDS |
| T_cell_CD8 | −0.146 | −0.239, −0.050 | 0.002 | HOLDS |

OLS + epithelial score shrinks the B-cell coefficient to ns
(β = −0.074, p = 0.19) while plasma stays (β = −0.202, p = 1.3 × 10⁻⁴).
Two models, two answers for B-cell; plasma is the robust LUAD signal.

**3. TCGA-LUAD — plasma yes; TLS no (sign flip); CD8 signature no**
(ABSOLUTE, n = 515)

| Signature | ρ_adj | 95% CI | p | call |
|---|---:|---|---:|---|
| TLS_Cabrita | **+0.092** | +0.004, +0.178 | 0.037 | OPPOSITE |
| B_cell | −0.104 | −0.187, −0.021 | 0.018 | HOLDS (weak) |
| Plasma_cell | −0.242 | −0.322, −0.156 | 2.9 × 10⁻⁸ | HOLDS |
| T_cell_CD8 | −0.063 | −0.153, +0.024 | 0.15 | NO_EVIDENCE |

CIBERSORT plasma matches the signature (ρ_adj = −0.248, p = 1.2 × 10⁻⁸).
CIBERSORT CD8 is negative (ρ_adj = −0.134, p = 0.002) even though the
CD8 *signature* is not. LUAD therefore supports “TACSTD2-high = plasma-low”
and, depending on the T-cell readout, a weak CD8-low; it does **not**
support “TACSTD2-high = TLS-low”.

**4. GSE81089 Uppsala mixed NSCLC (n = 197) — almost null**

Only TLS_imprint holds after epithelial correction (ρ_adj = −0.184,
CI −0.328 to −0.047, p = 0.0099). Cabrita TLS, B-cell, plasma and CD8
are all ns. Splitting by histology does not rescue a LUSC-like package:
SqCC (n = 106) keeps imprint and plasma after epithelial correction
(ρ ≈ −0.23, p ≈ 0.02) but loses Cabrita TLS and CD8; AC (n = 67) keeps
only Tfh (ρ_adj = −0.26, p = 0.038). This is an independent resected
cohort and it does **not** replicate TCGA-LUSC.

**5. Open ICI RNA — underpowered, no TACSTD2–outcome link**

| Cohort | n | Endpoint | TACSTD2 result |
|---|---:|---|---|
| GSE135222 anti-PD-1/L1 | 27 | PFS Cox | HR 1.03 (0.82–1.31), p = 0.78; +epithelial HR 0.95, p = 0.73 |
| GSE135222 | 27 | DCB ≥ 6 mo | MW p = 0.61 (7 DCB vs 20 NDB) |
| GSE126044 anti-PD-1 | 16 | RECIST | MW p = 0.44 (5 R vs 11 NR) |
| GSE207422 nivo/tori + chemo | 24 | MPR, pre | MW p = 0.34 (9 MPR vs 15 NMPR) |

GSE126044 CD8 signature *does* separate responders (rank-biserial r = +1.0,
p = 4.6 × 10⁻⁴). The immune scores are not dead; the TACSTD2 test is
simply too small. Crude TACSTD2–immune ρ in GSE207422 is strongly
negative (imprint −0.71, p = 1.2 × 10⁻⁴; CD8 −0.58, p = 0.003; B −0.52,
p = 0.009) and shrinks to ns after epithelial correction except imprint
(ρ_adj = −0.48, p = 0.022). With n = 24 that residual is a hint, not a
replication.

Pooled ICI partial correlations (n_total = 67) are all ns with I² up to
81% (`meta_by_histology_group.csv`).

### CLDN4

CLDN4 tracks the epithelial programme more tightly than TACSTD2
(ρ vs epithelial 0.47–0.85 across cohorts) and is a worse TLS-low marker
after purity. In TCGA-LUSC, CLDN4 vs Cabrita TLS is ns after ABSOLUTE
(ρ_adj = +0.014, p = 0.76); vs CD8 ρ_adj = −0.095, p = 0.034 (weak).
Do not treat CLDN4 as interchangeable with TACSTD2 for this hypothesis.

### Confounding (why crude ρ is not enough)

In TCGA, TLS / B / plasma all fall with ABSOLUTE purity
(LUAD ρ = −0.53 / −0.52 / −0.36; LUSC ρ = −0.49 / −0.50 / −0.52).
TACSTD2 itself is almost orthogonal to ABSOLUTE (LUAD ρ = 0.042;
LUSC ρ = 0.006) but correlates with the RNA epithelial score
(LUAD 0.48, LUSC 0.19). So: (i) crude TACSTD2–immune ρ is *not* a
purity artefact via ABSOLUTE in LUSC — TACSTD2 is not a purity gene;
(ii) epithelial-score adjustment in non-TCGA cohorts is conservative
and can eat a real signal because TACSTD2 is part of the epithelial
programme. That is why GSE72094 TLS_Cabrita dies after epithelial
correction while LUSC TLS survives ABSOLUTE.

### Meta-analysis (do not over-pool)

Random-effects meta of the two TCGA sets for TACSTD2 vs Cabrita TLS
after ABSOLUTE: ρ_RE = −0.048, p = 0.73, I² = 95%. Pooling LUAD + LUSC
hides a real histology interaction. The honest summary is
**histology-stratified**, not a single NSCLC ρ.

### What this does *not* show

- It does not show that TACSTD2-high tumours lack *histologic* TLS.
  We have no HE / CD20 / CD23 labels in these files.
- It does not show that TACSTD2 causes TLS exclusion.
- It does not validate the claim on atezolizumab / OAK / POPLAR /
  IMpower RNA (EGA). Those are the datasets that could actually test
  whether TACSTD2-high / TLS-low predicts ICI outcome.
- It does not survive as a general LUAD TLS effect.
- Open ICI n is too small to claim a treatment interaction.

### Practical takeaway for the Bessede / TROP2-ADC hypothesis

Use **TCGA-LUSC** as the public lung set where TACSTD2-high = TLS-low =
B-low = T-low after DNA purity correction. Use **GSE72094** as a LUAD
set where B / plasma / CD8 are low after an epithelial correction but
Cabrita TLS is not. Do **not** cite TCGA-LUAD as TLS-low. For ICI
outcome, wait for EGA (OAK / Patil) or a new open cohort with n ≫ 50.

Figures: `results/opus_tls/figures/fig1_heatmap_tacstd2_adj.png` through
`fig6_verdict_grid.png`.

### Follow-up: independent LUSC + formal histology interaction

The first pass left a hole: TCGA-LUSC was the only large set where the
full TLS+B+CD8 package survived DNA purity. Five extra open matrices
were scored the same way (`replication_lusc_histology.csv`).

**Formal interaction in pooled TCGA (n = 1008, ABSOLUTE + HC3).**
Signature ~ TACSTD2 + LUSC + TACSTD2×LUSC + purity. A negative
interaction means a steeper anti-immune TACSTD2 slope in LUSC.

| Signature | β_TACSTD2 in LUAD | β_interaction | p_interaction |
|---|---:|---:|---:|
| TLS_Cabrita | +0.068 (p=0.080) | **−0.215** | **3.2 × 10⁻⁴** |
| TLS_12chemokine | +0.010 (p=0.79) | **−0.285** | **6.6 × 10⁻⁷** |
| T_cell_CD8 | −0.063 (p=0.097) | **−0.191** | **8.0 × 10⁻⁴** |
| Tfh | −0.092 (p=0.014) | **−0.212** | **1.2 × 10⁻⁴** |
| IFNg_Ayers | +0.020 (p=0.61) | **−0.244** | **3.3 × 10⁻⁵** |
| B_cell | −0.100 (p=0.009) | −0.074 | 0.19 |
| Plasma_cell | −0.228 (p=1.5 × 10⁻⁷) | **+0.125** | 0.042 |

So the LUAD/LUSC split for TLS and CD8 is a real interaction, not two
unrelated p-values. B-cell is negative in both histologies (no
interaction). Plasma is *more* negative in LUAD.

**Independent LUSC / SCC after epithelial correction**

| Cohort | n | Cabrita TLS | B-cell | CD8 | call |
|---|---:|---|---|---|---|
| GSE4573 Raponi LUSC Affy | 130 | −0.124 p=0.16 | −0.165 p=0.061 | **−0.250 p=0.004** | CD8 + 12-chemokine (−0.263 p=0.0026) + IFN-γ hold; Cabrita/B do not |
| GSE17710 Wilkerson LUSC | 56 | −0.212 p=0.12 | −0.205 p=0.13 | −0.240 p=0.077 | all negative; only Tfh holds (p=0.045). Pathologist `tumor_percent` residual: all ns |
| GSE103584 SCC RNA-seq | 31 | −0.501 p=0.0048 | −0.428 p=0.018 | −0.382 p=0.037 | same direction as TCGA-LUSC; **UNDERPOWERED** (n<40), not a HOLDS call |
| GSE50081 SCC Affy | 42 | +0.113 p=0.48 | **+0.444 p=0.0037** | +0.227 p=0.15 | B / plasma / Tfh **OPPOSITE** |
| GSE19188 SCC Affy | 27 | +0.21 p=0.30 | +0.44 p=0.026 | +0.04 p=0.85 | underpowered; B opposite |

**Independent ADC arrays often go the other way.** GSE19188 ADC (n=45):
Cabrita TLS +0.34, B +0.44, CD8 +0.40 (all OPPOSITE). GSE50081 ADC
plasma +0.26 p=0.0037 (OPPOSITE).

**OS in TCGA is null for TACSTD2.** LUAD and LUSC Cox HR for TACSTD2
≈ 1.06 / 0.97, p > 0.4, with or without purity, TLS, or the
TACSTD2×TLS product. B-cell signature is protective in LUAD
(HR 0.78, p=0.007) — a positive control that the immune scores move OS
when they should. This slice is about *composition*, not prognosis.

**Revised takeaway.** The public set where TACSTD2-high = TLS+B+T-low
*after DNA purity* remains **TCGA-LUSC**. That histology interaction is
now a single-model p = 3×10⁻⁴ (Cabrita) / 8×10⁻⁴ (CD8). Independent
LUSC arrays do **not** reproduce the full package: GSE4573 supports
CD8 / 12-chemokine / IFN-γ only; GSE17710 is directional and ns;
GSE50081 SCC is opposite for B-lineage. Do not generalise “squamous =
Bessede pattern” beyond TCGA-LUSC without EGA ICI RNA.

Figures: `fig7_tcga_interaction.png`, `fig8_replication_grid.png`.

---

## 中文

### 问题

Bessede 及 TROP2 ADC 文献的核心表述：**肺肿瘤 TACSTD2 高 = TLS / B 细胞 /
T 细胞低，且不是纯度伪相关。** 本切片在所有可公开下载、单文件 < 2 GB
的肺 ICI 或图谱 RNA 上做检验。Patil / OAK / POPLAR / IMpower 在 EGA，
**只列清单、不下载**（见 `DATA_ACCESS.md`）。

预先规定的 HOLDS：n ≥ 40，纯度校正后偏 Spearman ρ < 0 且 p < 0.05。
TCGA 用 ABSOLUTE（DNA）± 甲基化白细胞比例；其余队列用泛上皮评分。

### 校正后哪里成立

**TCGA-LUSC（n = 493，ABSOLUTE）——整包成立，也是唯一整包成立的大队列。**
Cabrita TLS ρ = −0.186（CI −0.269 至 −0.100，p = 3.3×10⁻⁵）；
12-chemokine −0.313；imprint −0.243；B 细胞 −0.223；浆细胞 −0.144；
Tfh −0.337；CD8 −0.286；IFN-γ −0.263。四分位 Q4 vs Q1 与纯度三分位
Stouffer 仍显著。全转录组排序：TACSTD2 落在 CD8 / B / TLS 偏相关的
1.3% / 2.1% / 2.4% 分位（32 628 基因），经验 p 为 0.078 / 0.134 / 0.192，
故它在抗免疫尾部，但**不是**全基因组唯一离群点。LUSC 里 EPCAM 对 CD8
的校正 ρ 仅 −0.057，说明不是“任意上皮基因”。CIBERSORT CD8 同向
（−0.166，p = 2.2×10⁻⁴）；CIBERSORT 浆细胞在 LUSC **反向**
（+0.122，p = 0.0066），与签名浆细胞不一致，如实报告。

**GSE72094 LUAD 芯片（n = 442，上皮评分校正）——B / 浆 / CD8 成立，
Cabrita TLS 不成立。** TLS_Cabrita ρ = −0.019（CI 跨 0，p = 0.69）；
B −0.115（p = 0.016）；浆 −0.220（p = 2.9×10⁻⁶）；CD8 −0.146（p = 0.002）。
OLS 加上皮评分后 B 细胞系数不再显著（p = 0.19），浆细胞仍在。

**TCGA-LUAD（n = 515，ABSOLUTE）——浆细胞成立；TLS 反向；CD8 签名不成立。**
TLS_Cabrita ρ = **+0.092**（CI +0.004 至 +0.178，p = 0.037，OPPOSITE）；
B −0.104（p = 0.018，弱）；浆 −0.242（p = 2.9×10⁻⁸）；CD8 签名 −0.063
（p = 0.15）。CIBERSORT 浆细胞一致（−0.248）；CIBERSORT CD8 为负
（−0.134，p = 0.002）。故 LUAD 支持“TACSTD2 高 = 浆细胞低”，
**不支持**“= TLS 低”。

**GSE81089 乌普萨拉混合 NSCLC（n = 197）——几乎全空。** 仅 imprint
成立（ρ = −0.184，p = 0.0099）。按腺癌 / 鳞癌切开也组不成 LUSC 那样的
TLS+B+CD8 包。这是独立手术队列，**不能**当作 TCGA-LUSC 的重复。

**公开 ICI RNA 全部 n ≤ 27，TACSTD2 与疗效无关联。**
GSE135222 PFS HR 1.03（0.82–1.31），p = 0.78；DCB p = 0.61。
GSE126044 RECIST p = 0.44（但 CD8 签名本身能分开缓解者，p = 4.6×10⁻⁴，
说明签名没死）。GSE207422 MPR p = 0.34；粗相关很负，上皮校正后只剩
imprint（p = 0.022），n = 24 只能当提示。

### CLDN4

比 TACSTD2 更贴上皮程序，纯度校正后不是可靠的 TLS-low 标记。
LUSC 中 CLDN4 vs Cabrita TLS 校正后 ns（ρ = +0.014，p = 0.76）。
不要把 CLDN4 和 TACSTD2 当成同一假设。

### 为什么必须做纯度校正

TCGA 中 TLS / B / 浆细胞都随 ABSOLUTE 纯度下降（LUAD 约 −0.5）。
TACSTD2 **几乎不随 ABSOLUTE 变**（LUAD 0.042，LUSC 0.006），但与 RNA
上皮分相关。因此：LUSC 的 TACSTD2–免疫负相关不是“纯度基因”伪影；
非 TCGA 队列的上皮评分校正偏保守，可能吃掉真信号（GSE72094 的
Cabrita TLS 就是这样死的）。

两套 TCGA 随机效应合并 Cabrita TLS：ρ_RE = −0.048，p = 0.73，I² = 95%。
**不能**把 LUAD+LUSC 捏成一个 NSCLC 相关系数。

### 不能声称的事

没有病理 TLS 金标准；没有因果；没有 OAK / Patil 的 ICI 结局验证；
不是普遍的 LUAD TLS 效应；公开 ICI 样本不够谈治疗交互。

### 给后续切片的用法

要引用“公开数据里 TACSTD2 高 = TLS/B/T 低且过了纯度校正”，请写
**TCGA-LUSC**。要写 LUAD，请写 GSE72094 的 B/浆/CD8，并写明 Cabrita
TLS 不成立。不要引用 TCGA-LUAD 为 TLS-low。ICI 结局等 EGA（OAK/Patil）
或 n ≫ 50 的新公开队列。

### 续做：独立鳞癌队列 + 组织学交互

合并 TCGA LUAD+LUSC（n=1008，ABSOLUTE，HC3）做
签名 ~ TACSTD2 + LUSC + TACSTD2×LUSC + 纯度。负交互 = 鳞癌里
TACSTD2 的抗免疫斜率更陡：Cabrita TLS β_int = −0.215，p = 3.2×10⁻⁴；
CD8 −0.191，p = 8.0×10⁻⁴；12-chemokine p = 6.6×10⁻⁷。B 细胞两组织学
都负、交互不显著（p=0.19）。浆细胞在腺癌更负（交互 +0.125，p=0.042）。

独立公开鳞癌：GSE4573（n=130）只对 CD8 / 12-chemokine / IFN-γ 成立，
Cabrita TLS 与 B 细胞不成立。GSE17710（n=56）方向一致但除 Tfh 外不显著，
病理 `tumor_percent` 校正后全空。GSE103584 鳞癌 n=31 方向同 TCGA-LUSC
但按规则算 UNDERPOWERED。GSE50081 鳞癌 B/浆细胞反向。若干腺癌芯片
（GSE19188 ADC）是 TACSTD2 高 = 免疫高。TCGA 里 TACSTD2 与 OS 无关
（HR≈1，p>0.4）；B 细胞签名在 LUAD 保护（HR 0.78，p=0.007）。

结论不变且更硬：能写进“纯度校正后整包成立”的公开集仍是
**TCGA-LUSC**；鳞癌不是普遍规律；ICI 结局仍要等 EGA。
