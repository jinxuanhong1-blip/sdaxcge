# w200/C4_GSE22493_ifn — GSE22493 CLDN4 siRNA：IFN / MHC-I / APM 基因列表
# w200/C4_GSE22493_ifn — GSE22493 CLDN4 siRNA: IFN / MHC-I / APM gene list

C4 切片：只分析 **GSE22493**（SKOV-3-IP-Luc，CLDN4 慢病毒 siRNA vs CLDN4 过表达对照，双色芯片，n=3）。不问其他 CLDN4-loss 数据集，不把本切片写成“跨队列重复”。

C4 slice: **GSE22493 only** (SKOV-3-IP-Luc; lentiviral CLDN4 siRNA vs CLDN4-overexpression control; two-color spotted array; n=3). Other CLDN4-loss accessions are out of scope. This is not a replication meta-analysis.

> **一句话 / Bottom line**
>
> **中文**：这张芯片**不支持**“CLDN4 敲低打开 IFN / MHC-I / APM”。优先 6 基因在作者沉积 log2 比值上是 3 上 / 3 下，且无一达到名义 p < 0.05；MHC-I/APM 12 个可测基因里 9 个中位数向下。ScanArray 原始 Cy5/Cy3 更糟：优先基因全部向下。CLDN4 本身在阵列上**不能干净确认敲低**，对照是过表达而不是 scramble/WT。低置信度、偏负或无效的证据，不是机制成立的证据。
>
> **English**: This array does **not** support “CLDN4 knockdown opens IFN / MHC-I / APM”. The six priority genes are 3 up / 3 down on the deposited log2 ratios, none at nominal p < 0.05; 9/12 measurable MHC-I/APM genes have a negative median. Raw ScanArray Cy5/Cy3 is worse: every priority gene is down. CLDN4 itself is **not cleanly confirmed** as knocked down on the array, and the control is overexpression, not scramble/WT. Low-confidence negative-to-null evidence — not a mechanistic confirmation.

---

## 中文

### 数据是什么（先把设计说清楚）

| 项 | 事实 |
|---|---|
| 登录号 | [GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493)（无配套论文；提交者 Zhijian Gao, BWH） |
| 平台 | [GPL10555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10555) BWH Human Release 3.0（Operon 60-mer 点制芯片，~36k 探针） |
| 细胞 | SKOV-3-IP-Luc（人卵巢癌），**不是肺** |
| 样本 | GSM558700 / 701 / 702，各一张双色芯片 |
| 通道 | ch1 **Cy3 = CLDN4 过表达（对照）**；ch2 **Cy5 = CLDN4 慢病毒 siRNA** |
| 沉积值 | series matrix `VALUE` = **log2(敲低/对照)**（作者已归一化） |
| 未做的事 | 没有 scramble/WT 对照；GEO 文本写了 dye-swap，但三张芯片染料方向相同，仓库里没有独立的 swap 样本 |

这是 **过表达 vs 敲低**，不是 **野生型 vs 敲低**。即使 IFN 基因真的变了，也不能单独归因于“CLDN4 丢失打开抗原呈递”。

### 基因列表（预先固定，不是看完结果再挑）

- **优先 6 基因**：IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A（用户点名；TAP1/TAP2 放在 APM 面板里一起报）
- **MHC-I / APM（16）**：HLA-A/B/C, B2M, NLRC5, TAP1, TAP2, TAPBP, PSMB8, PSMB9, PSMB10, ERAP1, ERAP2, CALR, CANX, PDIA3
- **IFN_IMMUNE（73）**：I 型 ISG + IFN 信号 + 部分趋化因子/配体（与平行 CLDN4-loss 切片同一列表）

**别名（写明，避免假装“ISG15 不在芯片上”）**：GPL10555 没有 `ISG15` ORF，但有历史符号 **G1P2**（IFI-15K）。本切片把它映射为 ISG15，并在表里标注 `G1P2->ISG15`。`IFI6` 同理用了 `G1P3`。

芯片上**没有注释**的 APM/IFN 基因：NLRC5, ERAP1, ERAP2, PDIA3, MX2, RSAD2, IFITM1, IFITM3, IFI44L, DDX60, XAF1, CMPK2, EIF2AK2, IRF9, IL32, IFNL1。CXCL9、CCL5 只有 1 张芯片有值，不进入基因集检验。

### 方法

1. **主分析**：GEO series matrix 沉积 `VALUE`。基因在每张芯片上取该基因全部探针的中位数，再在阵列间取中位数/均值。n≥2 时做对 0 的单样本 t（**描述性**；n=3 检验力极低）。面板内 BH-FDR。
2. **基因集**：每个可测基因（≥2 阵列）的阵列均值 log2，对集合 vs 背景做双侧 Mann-Whitney U。
3. **敏感性**：ScanArray `Ch2 Median−B` / `Ch1 Median−B` 的 log2（两通道背景减后均须 >0）。这不是作者的正式比值；只用来看沉积值会不会被原始强度推翻。
4. **不制造 p 值**：缺失基因保持 NA。不做 limma/全基因组 DE 再回填。

### CLDN4 敲低体检（先过这一关）

`cldn4_diagnostic.tsv`，探针 17169：

| 阵列 | 沉积 log2(KD/ctrl) | ScanArray Cy3 / Cy5（Median−B） | ScanArray log2 |
|---|---:|---|---:|
| GSM558700 | 缺失 | 89 / **−39** | 不可用（Cy5 背景减后非正） |
| GSM558701 | **−1.74** | 251 / 584 | **+1.22（与沉积反号）** |
| GSM558702 | **−0.71** | 1031 / 330 | −1.64 |

沉积的两张可用芯片方向为负（敲低 < 过表达对照），但第三张缺失，且 GSM558701 的原始强度与沉积值**反号**。对照还是过表达。**不能声称阵列水平确认了干净的 CLDN4 敲低。** 下游 IFN 解读必须带着这个前提。

### 优先 6 基因（主结果）

来自 `priority_genes.tsv`（沉积值；基因 = 每阵列探针中位数）：

| 基因 | 映射 | 阵列 log2 | 中位数 | 方向 | p（单样本 t） |
|---|---|---|---:|---|---:|
| IFI27 | ORF | −2.64 / **+2.84** / −2.32 | −2.32 | DOWN | 0.73 |
| OAS2 | ORF | −0.29 / +2.81 / +0.11 | +0.11 | UP | 0.46 |
| IFIT1 | ORF | −2.84 / −1.12 / −2.12 | −2.12 | DOWN | 0.055 |
| MX1 | ORF | −0.07 / NA / +0.79 | +0.36 | UP | 0.56 |
| ISG15 | **G1P2→ISG15** | +0.71 / +1.79 / −1.06 | +0.71 | UP | 0.62 |
| HLA-A | 11 探针，DESCRIPTION 前缀 | −1.60 / +0.80 / −1.74 | −1.60 | DOWN | 0.41 |

- 3 上 / 3 下。**没有任何优先基因在 BH-FDR 或甚至名义 p < 0.05 下显著。**
- 最接近一致的是 **IFIT1 向下**（三张芯片全负，p = 0.055，仍不是显著上调）。
- IFI27 / OAS2 / HLA-A 在 GSM558701 上经常反号——单张芯片在拖均值。
- HLA-A 的 11 个探针彼此打架（见 `probe_level.tsv`：例如 7509 为 −3.84 / +1.35 / −6.64，8068 为 +2.31 / −0.86 / +1.34）。**不能说 MHC-I 经典分子上升。**

### MHC-I / APM

`apm_genes.tsv`：12/16 可测。中位数向上的只有 **PSMB8、PSMB9、CALR**。向下的包括 **HLA-A/B/C、B2M、TAP1、TAP2、TAPBP、PSMB10、CANX**。NLRC5 / ERAP1 / ERAP2 / PDIA3 不在平台上。

这不是“抗原呈递机器被打开”。免疫蛋白酶体亚基（PSMB8/9）略偏正，但 TAP 和经典 MHC-I / B2M 偏负，而且同样被 GSM558701 的反号拖着。没有一个 APM 基因的面板内 q < 0.05。

### 更宽的 IFN 集合

`geneset_stats.tsv`（对背景、双侧 MW）：

| 集合 | 可测/列表 | 均值向上/向下 | 集合中位 log2 | MW p | 相对背景 |
|---|---|---|---:|---:|---|
| PRIORITY | 6/6 | 3 / 3 | −0.17 | 0.67 | DOWN |
| APM | 12/16 | 3 / 9 | −0.45 | 0.16 | DOWN |
| IFN_IMMUNE | 58/73 | 25 / 33 | −0.23 | 0.29 | DOWN |

集合水平**没有**显著的 IFN 上移。背景基因中位 log2 ≈ −0.03，IFN 集合还略更负。个别基因（JAK2 沉积中位 +1.80，p = 0.033；PSMB9 +1.71，p = 0.091）不能被说成通路打开——未过面板 FDR，且与 HLA/B2M/TAP/IFIT1 的方向相反。

### 敏感性：沉积值 vs ScanArray 原始比值

`sensitivity_deposited_vs_scanarray.tsv`：

| 基因 | 沉积中位 | ScanArray 中位 | 符号一致？ |
|---|---|---|---|
| IFI27 | DOWN (−2.32) | DOWN (−2.41） | 是（原始三张全负；沉积的 +2.84 是 701 的反号） |
| OAS2 | UP (+0.11) | **DOWN (−1.53)** | **否** |
| IFIT1 | DOWN (−2.12) | DOWN (−3.69) | 是 |
| MX1 | UP (+0.36) | **DOWN (−0.56)** | **否** |
| ISG15 | UP (+0.71) | **DOWN (−1.69)** | **否** |
| HLA-A | DOWN (−1.60) | DOWN (−1.20) | 是 |
| CLDN4 | DOWN (−1.23) | DOWN (−0.21，两张可用且反号） | 弱 |

作者沉积比值里仅有的三个“向上”优先基因（OAS2、MX1、ISG15），在原始 Cy5/Cy3 上都变成向下。**如果有人只引用沉积中位数说“ISG15/OAS2 升了”，那是选择性报告。** 反过来，也不能把 ScanArray 的“全向下”写成“CLDN4 敲低抑制 IFN”——敲低未确认、对照是过表达、芯片又旧又吵。

### 诚实结论

1. **GSE22493 不能当作“CLDN4 丢失 → IFN/MHC-I/APM 打开”的公开阳性证据。**
2. 它也不是干净的阴性：设计是过表达 vs 敲低，CLDN4 探针自相矛盾，沉积值与原始强度在 GSM558701 上多次反号。
3. 能说的上限：在这三张 2010 年卵巢癌双色芯片上，优先 IFN/MHC 列表**没有**可辩护的上调；经典 MHC-I / B2M / TAP 偏下；IFIT1 是最一致的下调。
4. 不外推到肺、不外推到 ICI、不与未公开的私有 RNA-seq 强行对齐。

### 局限（请原样保留）

- n=3；单样本 t 和 MW 都只是描述。
- 自定义旧芯片；HLA 多探针不一致；部分 ISG 用历史符号。
- 文本声称 dye-swap，沉积样本没有。
- 提交者写敲低用 Western / qRT-PCR / GFP 验证过，但这些读数不在 GEO 文件里。
- 无论文、无作者官方 DE 表。
- 卵巢 SKOV-3，不是 NSCLC。

### 复现

```bash
pip install pandas numpy scipy statsmodels matplotlib
python3 scripts/w200/C4_GSE22493_ifn/download_data.py   # → /tmp/w200_c4_gse22493_ifn
python3 scripts/w200/C4_GSE22493_ifn/run_analysis.py    # → results/w200/C4_GSE22493_ifn/
```

公开文件 md5（`download_data.py` 打印）：series matrix `8f8a07cb8f3f3a90396c3feb1eb3ea0f`；family SOFT `5d88cef4c873370970e6ba6be82d7b48`；RAW tar `9bae4611ed5e75f28b387b32f8f793db`。

---

## English

### What this file is

GEO [GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493) / [GPL10555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10555): three two-color Operon 60-mer arrays (GSM558700–702) on **SKOV-3-IP-Luc ovarian** cells. Channel 1 is Cy3 **CLDN4 overexpression (control)**; channel 2 is Cy5 **CLDN4 lentiviral siRNA**. The deposited `VALUE` is the author-normalized **log2(knockdown/control)**. There is no scramble/WT arm. The record mentions dye-swaps; all three GSMs have the same dye assignment.

This contrast is **overexpression versus knockdown**, not wild-type versus knockdown, and it is **not lung**.

### Fixed lists

Priority: IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A. APM (16) and IFN_IMMUNE (73) as in `scripts/w200/C4_GSE22493_ifn/gene_sets.py`. ISG15 is present only as historical **G1P2** (IFI-15K) and is labelled `G1P2->ISG15`. Absent on the platform: NLRC5, ERAP1, ERAP2, PDIA3 and several ISGs listed in `geneset_stats.tsv`.

### Methods

Primary: deposited series-matrix log2 ratios; per-array median across probes, then median/mean across arrays; one-sample t-test vs 0 when n≥2 (descriptive). Gene-set shift: two-sided Mann-Whitney of gene-mean log2 vs background. Sensitivity: ScanArray log2((Ch2 Median−B)/(Ch1 Median−B)) when both channels are positive. No p-values are invented for missing genes.

### CLDN4 check

Probe 17169 is missing on GSM558700 (Cy5 background-subtracted median = −39). The two deposited values are negative (−1.74, −0.71), but GSM558701 raw Cy5/Cy3 is **+1.22** — opposite sign. Knockdown is **not confirmed** on the array. The control is overexpression.

### Priority panel

3 up / 3 down on deposited medians. **None** reach nominal p < 0.05 (IFIT1 is the most consistent: all three arrays negative, p = 0.055, direction **down**). IFI27, OAS2 and HLA-A flip on GSM558701. HLA-A’s 11 probes disagree with each other (`probe_level.tsv`). This is not an MHC-I-up call.

### APM and IFN set

12/16 APM genes measured: only PSMB8, PSMB9, CALR have a positive median; HLA-A/B/C, B2M, TAP1, TAP2, TAPBP, PSMB10, CANX are negative. No APM gene has within-panel q < 0.05. IFN_IMMUNE (58 measured): median mean-log2 −0.23 vs background −0.03, Mann-Whitney p = 0.29. No set-level IFN up-shift.

### Sensitivity

The only deposited “up” priority genes (OAS2, MX1, ISG15) are **down** on raw ScanArray Cy5/Cy3. IFI27 raw intensities are negative on all three arrays; the deposited +2.84 on GSM558701 is a sign flip relative to the TIFF quantitation. Citing deposited OAS2/ISG15 as “up after CLDN4 siRNA” is selective. Citing ScanArray “all down” as “CLDN4 loss represses IFN” is equally unsupported: knockdown is unconfirmed, the control is overexpression, and the array is noisy.

### Honest ceiling

GSE22493 is **not** public positive evidence that CLDN4 siRNA opens IFN / MHC-I / APM. It is also not a clean negative, because the design and the CLDN4 probe are compromised. The most that can be said: on these three 2010 ovarian two-color arrays, the pre-specified IFN/MHC list has **no defensible up-signature**; classical MHC-I / B2M / TAP trend down; IFIT1 is the most consistent down gene. Do not export this to lung, ICI, or a private RNA-seq story.

### Reproduce

Same commands as the Chinese section. Outputs live in `results/w200/C4_GSE22493_ifn/`.
