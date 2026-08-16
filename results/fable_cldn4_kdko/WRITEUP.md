# CLDN4 loss-of-function meta-analysis: TACSTD2, junctions, and IFN/immune genes

*Parallel slice `fable_cldn4_kdko`. All code in `scripts/fable_cldn4_kdko/`,
data and tables in `results/fable_cldn4_kdko/`, notes in
`notes/fable_cldn4_kdko/`.*

---

# English

## Question
Across every public dataset in which CLDN4 (claudin-4) is genuinely
knocked down / knocked out (shRNA / siRNA / CRISPR / genetic KO), does the
loss of CLDN4 change **TACSTD2 (Trop-2)**, **cell-junction genes**, and
**interferon / immune genes**? Lung tissue prioritised; other tissues kept
only when the perturbation is a *bona fide* CLDN4 perturbation.

## ★ Directional focus (user priority): does CLDN4 loss OPEN IFN/MHC-I/APM?
The user's private RNA-seq shows CLDN4 KD/KO **up-regulating** the interferon /
MHC-class-I / antigen-presentation machinery (IFI27, OAS2, IFIT1, MX1, ISG15,
HLA-A, TAP1/2). I tested one-sidedly which public CLDN4-loss datasets reproduce
that **upward** direction (`scripts/fable_cldn4_kdko/panel_ifn_apm.py`,
`results/fable_cldn4_kdko/panel_ifn_apm_results.json`,
`panel_core_by_dataset.csv`).

**Verdict — one public dataset reproduces it, and it is the in-vivo epithelial (lung) KO:**

| Dataset | Model | Extended IFN+MHC-I+APM panel | median log2FC | one-sided test (panel UP) | Reproduces "opening"? |
|---|---|---|---|---|---|
| **GSE50927** | mouse **lung** Cldn4 KO | **34/39 genes up** | **+0.43** | Mann-Whitney vs bg **p = 7.2e-11**; Wilcoxon **p = 3.1e-7**; sign test **p = 1.2e-6** | **✅ YES, strongly** |
| GSE207704 | human breast CRISPR CLDN4-/- | 7/20 up | -0.31 | p(up) = 0.94 (i.e. it *closes*) | ❌ opposite |
| GSE22493 | human ovarian CLDN4 siRNA (low conf.) | 16/34 up | -0.26 | p(up) = 0.96 | ❌ opposite |

In **GSE50927** the user's exact core panel behaves as their private data predicts:
**OAS2 +1.49 (FDR 1.5e-3)**, **ISG15 +1.08 (FDR 3.1e-3)**, **IFIT1 +0.58 (FDR 2.8e-2)**,
**IFI27/Ifi27l2a +0.62 (FDR 0.22)**, **MX1 +0.63 (p 0.12)**, **TAP1 +0.18** — 6/8 up.
The MHC-I/APM arm opens through **Nlrc5 +0.75 (FDR 0.09, the MHC-I transactivator)**,
**B2m +0.66 (FDR 3.4e-2)**, **Psmb9 +0.90 (FDR 3.7e-2)**, **Psmb8 +0.43 (FDR 0.29)**,
**Stat1 +0.39**, **Irf1 +0.38**, **H2-D1 +0.35** (an HLA-A/B ortholog); only H2-K1
(-0.07) and Tap2 (-0.01) are flat. Broader ISGs corroborate: **Oas3 +1.82 (FDR 1.6e-3)**,
**Rsad2 +1.20**, **Oas1a +0.76**. This is a coherent, multi-gene interferon+antigen-presentation
"opening," exactly matching the private-data direction.

**Honest bottom line for the priority question:** the public evidence for
CLDN4-loss *opening* IFN/MHC-I/APM is real but **specific to the in-vivo,
non-transformed lung epithelium (GSE50927)**. In the two human *cancer cell line*
perturbations the module moves the **opposite** way (down). So the private-data
direction is corroborated by the biologically closest public model (epithelial,
in vivo), while cancer-cell-line KO/KD is not a faithful surrogate for it.
Caveat: GSE50927's core WT-vs-KO contrast has limited biological replication, but
the signal rests on 34/39 panel genes moving up concordantly (several at FDR<0.05),
which is not explicable by single-gene noise.

## What was searched and what qualified
Exhaustive NCBI GEO (`gds`) + SRA search (see
`notes/fable_cldn4_kdko/dataset_inventory.md` for the full term list and the
complete accepted/rejected table). Every candidate accession was verified at
the **sample level** through the GEO `acc.cgi` interface, not just by title.

Three datasets are genuine CLDN4 loss-of-function experiments with usable
processed data (all files < 14 MB; nothing hit the 2 GB skip rule):

| Accession | Organism / tissue | Perturbation | Data used |
|---|---|---|---|
| **GSE207704** | Human breast cancer (T47D, MCF7) | CRISPR **CLDN4-/-** vs WT | RNA-seq FPKM (per group) |
| **GSE50927** | **Mouse lung** | genetic **Cldn4 KO** vs WT | author edgeR DE table |
| **GSE22493** | Human ovarian (SKOV-3) | **CLDN4 siRNA** vs CLDN4-high control | 2-colour array (low confidence) |

Perturbation confirmed by the CLDN4/Cldn4 signal itself: FPKM drops in both
breast lines (MCF7 86.6→52.3, T47D 43.4→20.6); mouse lung `Cldn4`
logFC = **-6.06, FDR = 4.1e-26**. GSE22493's own CLDN4 probe is inconsistent
across arrays (log2 KD/ctrl = -1.64 / +1.22 / NA), so it is reported as
supporting-only. Roughly a dozen other CLDN4-mentioning series (ceRNA tumour
panels, other-claudin KOs — CLDN7 GSE26055, CLDN18 GSE48443 —, lung marker
studies, a 4C-seq locus study SRP263109) were rejected; reasons are tabulated
in the notes.

## Statistics
- **GSE50927**: author edgeR `logFC / PValue / FDR` per gene (KO vs WT).
- **GSE207704**: log2FC(KO/WT) per cell line from FPKM (+1 pseudocount),
  averaged over the two lines; genes filtered to FPKM ≥ 1.
- **GSE22493**: per-array log2(Cy5 KD / Cy3 control) on background-subtracted
  medians, median-centred per array, one-sample t-test across the 3 arrays.
- **Gene-set level**: two-sided Mann-Whitney U comparing each marker set's
  fold-changes against all other genes (`geneset_stats.csv`).

## Results

### 1) TACSTD2 / Trop-2 — context dependent
- **Breast CRISPR CLDN4-/- (GSE207704): TACSTD2 goes DOWN**, consistently in
  both lines — log2FC = **-1.00 (MCF7)** and **-0.62 (T47D)**, mean **-0.81**
  (~1.4–2×). So in breast cancer epithelium, losing CLDN4 lowers Trop-2.
- **Mouse lung Cldn4 KO (GSE50927): TACSTD2 unchanged** — `Tacstd2`
  logFC = -0.14, p = 0.49, FDR = 1.0.
- **Ovarian (GSE22493): not testable** — TACSTD2/Trop-2 is not on that array.

There is a real CLDN4→Trop-2 coupling in breast cancer cells, but it does not
generalise to normal mouse lung.

### 2) Junction genes — no collapse, no compensation
- **Breast**: a modest but statistically significant coordinated **down**
  shift of the junction set (median -0.088 vs +0.012 background,
  Mann-Whitney **p = 0.049**, n = 30). Individual claudins/TJ genes move only
  slightly (CLDN3 -0.30, CLDN7 -0.14, OCLN -0.18, CDH1 -0.05); no compensatory
  claudin is switched on.
- **Lung**: junction set shift is trivial in size (median -0.047) though
  formally significant (**p = 0.019**, n = 45); every individual junction gene
  has FDR = 1.0 (`Cldn3` +0.09, `Cldn7` +0.11, `Ocln` -0.05, `Tjp1` -0.10,
  `Cdh1` -0.14). Biologically, the tight-junction program is intact.
- **Ovarian**: no junction-set shift (p = 0.91).

Takeaway: CLDN4 loss does **not** trigger tight-junction collapse or a
compensatory claudin-switch at the transcript level in any of the three
systems.

### 3) IFN / immune genes — perturbed, and direction flips by context
This is the strongest and most interesting signal.
- **Mouse lung Cldn4 KO: interferon program is UP**, highly significant
  (median +0.39 vs 0, Mann-Whitney **p = 4.8e-13**, n = 69). Core ISGs rise
  with real FDRs: **Isg15 +1.08 (FDR 3.1e-3)**, **B2m +0.66 (FDR 3.4e-2)**,
  **Stat1 +0.39 (p = 0.026)**. CLDN4-null lung is in a heightened
  interferon/antigen-presentation state.
- **Breast CRISPR CLDN4-/-: interferon program is DOWN**, also highly
  significant (median -0.20 vs +0.012, **p = 4.7e-5**, n = 30); e.g. ISG15
  -1.87 in T47D.
- **Ovarian siRNA (low confidence): IFN set DOWN** (median -0.35,
  **p = 0.025**), directionally matching breast, but the shaky CLDN4 probe
  means this only weakly corroborates.

So CLDN4 loss reproducibly *engages* the interferon/immune module, but the
sign is context-dependent: **up** in mouse lung in vivo, **down** in human
breast/ovarian cancer cells in vitro.

## Bottom line
1. **TACSTD2/Trop-2**: down ~1.4–2× when CLDN4 is knocked out in breast
   cancer cells; unchanged in normal mouse lung; untestable in the ovarian
   array. A CLDN4–Trop-2 link exists but is cancer-epithelium-specific.
2. **Junctions**: no transcriptional collapse and no compensatory claudin
   induction anywhere; at most a small coordinated dip.
3. **IFN/immune**: the clearest effect — CLDN4 loss perturbs the interferon
   program strongly (up in lung KO, down in breast/ovarian cancer).

## Caveats
- GSE207704 provides only per-group FPKM (two cell lines act as the
  replicates), so cell-line-level fold-changes are solid but within-group
  variance is unavailable; the two lines are treated as biological replicates.
- GSE50927's core WT-vs-KO contrast has low biological replication; the ISG
  result rests on many genes moving concordantly (multi-gene FDR < 0.05),
  which argues against single-gene noise, but replication is limited.
- GSE22493 is an old two-colour Operon array whose control channel is
  CLDN4-*overexpressing* (not scramble) and whose CLDN4 probe is inconsistent,
  so it is treated as supporting evidence only.

## Reproduce
```
pip install pandas numpy scipy statsmodels
cd scripts/fable_cldn4_kdko
python3 search_geo.py          # GEO candidate sweep
python3 search_broaden.py      # SRA + broad nets
python3 verify_accessions.py   # sample-level verification
# download processed files listed in dataset_inventory.md into results/.../raw/
python3 analyze.py             # DE + gene-set statistics -> results/*.csv, *.json
```

---

# 中文

## 问题
在所有公开的、真正对 CLDN4（claudin-4）进行敲低/敲除的数据集中
（shRNA / siRNA / CRISPR / 基因敲除），CLDN4 缺失是否会改变
**TACSTD2（Trop-2）**、**细胞连接基因** 和 **干扰素/免疫基因**？
优先肺组织；其他组织仅在确属真实 CLDN4 扰动时纳入。

## ★ 方向性重点（用户优先级）：CLDN4 缺失是否"打开"IFN/MHC-I/APM？
用户的私有 RNA-seq 显示 CLDN4 敲低/敲除会**上调**干扰素 / MHC-I 类 /
抗原提呈机器（IFI27、OAS2、IFIT1、MX1、ISG15、HLA-A、TAP1/2）。我用单侧检验
考察了哪些公开的 CLDN4 缺失数据集能重现这一**上调**方向
（`scripts/fable_cldn4_kdko/panel_ifn_apm.py`、
`results/fable_cldn4_kdko/panel_ifn_apm_results.json`、`panel_core_by_dataset.csv`）。

**结论——只有一个公开数据集重现了它，而且正是体内上皮（肺）敲除：**

| 数据集 | 模型 | 扩展 IFN+MHC-I+APM 面板 | 中位 log2FC | 单侧检验（面板上调） | 重现"打开"？ |
|---|---|---|---|---|---|
| **GSE50927** | 小鼠**肺** Cldn4 敲除 | **34/39 基因上调** | **+0.43** | Mann-Whitney vs 背景 **p = 7.2e-11**；Wilcoxon **p = 3.1e-7**；符号检验 **p = 1.2e-6** | **✅ 是，且很强** |
| GSE207704 | 人乳腺 CRISPR CLDN4-/- | 7/20 上调 | -0.31 | p(上调) = 0.94（即被*关闭*） | ❌ 相反 |
| GSE22493 | 人卵巢 CLDN4 siRNA（低置信度） | 16/34 上调 | -0.26 | p(上调) = 0.96 | ❌ 相反 |

在 **GSE50927** 中，用户点名的核心面板正如其私有数据所预测：
**OAS2 +1.49（FDR 1.5e-3）**、**ISG15 +1.08（FDR 3.1e-3）**、**IFIT1 +0.58（FDR 2.8e-2）**、
**IFI27/Ifi27l2a +0.62（FDR 0.22）**、**MX1 +0.63（p 0.12）**、**TAP1 +0.18**——8 个里 6 个上调。
MHC-I/APM 分支通过 **Nlrc5 +0.75（FDR 0.09，MHC-I 主转录激活因子）**、
**B2m +0.66（FDR 3.4e-2）**、**Psmb9 +0.90（FDR 3.7e-2）**、**Psmb8 +0.43**、
**Stat1 +0.39**、**Irf1 +0.38**、**H2-D1 +0.35**（HLA-A/B 的直系同源）而"打开"；
只有 H2-K1（-0.07）和 Tap2（-0.01）持平。更广泛的 ISG 佐证：
**Oas3 +1.82（FDR 1.6e-3）**、**Rsad2 +1.20**、**Oas1a +0.76**。这是一个连贯的、
多基因的干扰素+抗原提呈"打开"，方向与私有数据完全一致。

**对该优先问题的诚实结论：** 公开数据中支持 CLDN4 缺失"打开"IFN/MHC-I/APM 的
证据是真实的，但**特定于体内、非转化的肺上皮（GSE50927）**。在两个人*癌细胞系*
的扰动中，该模块方向**相反**（下调）。因此私有数据方向被生物学上最接近的公开模型
（体内、上皮）所印证，而癌细胞系 KO/KD 并非其忠实替身。
注意：GSE50927 的核心 WT-vs-KO 对比生物学重复有限，但结论建立在 34/39 个面板基因
一致上调（多个 FDR<0.05）之上，无法用单基因噪声解释。

## 检索与纳入
在 NCBI GEO（`gds`）与 SRA 上进行了穷尽式检索（完整检索词与全部
纳入/排除表见 `notes/fable_cldn4_kdko/dataset_inventory.md`）。每个候选号
都通过 GEO `acc.cgi` 在**样本层面**逐一核验，而非仅凭标题。

三个数据集属于真正的 CLDN4 功能缺失实验且有可用的处理后数据
（所有文件 < 14 MB，无一触及 2 GB 跳过阈值）：

| 编号 | 物种/组织 | 扰动方式 | 使用数据 |
|---|---|---|---|
| **GSE207704** | 人乳腺癌（T47D, MCF7） | CRISPR **CLDN4-/-** vs 野生型 | RNA-seq FPKM（按组） |
| **GSE50927** | **小鼠肺** | 基因 **Cldn4 敲除** vs 野生型 | 作者 edgeR 差异表 |
| **GSE22493** | 人卵巢癌（SKOV-3） | **CLDN4 siRNA** vs CLDN4 高表达对照 | 双色芯片（低置信度） |

扰动本身已由 CLDN4/Cldn4 信号证实：两株乳腺系 FPKM 下降
（MCF7 86.6→52.3，T47D 43.4→20.6）；小鼠肺 `Cldn4`
logFC = **-6.06，FDR = 4.1e-26**。GSE22493 自身的 CLDN4 探针在各芯片间
方向不一致（log2 KD/对照 = -1.64 / +1.22 / NA），故仅作辅助证据。
另有十余个仅"提及"CLDN4 的系列（ceRNA 肿瘤组织、其他 claudin 敲除——
CLDN7 的 GSE26055、CLDN18 的 GSE48443——、肺 marker 研究、位点 4C-seq
SRP263109）被排除，理由见 notes。

## 统计方法
- **GSE50927**：直接采用作者 edgeR 的 `logFC / PValue / FDR`（KO vs WT）。
- **GSE207704**：由 FPKM（+1 伪计数）计算每株细胞的 log2FC(KO/WT)，
  再对两株取平均；仅保留 FPKM ≥ 1 的基因。
- **GSE22493**：对背景扣除中位强度计算每张芯片
  log2(Cy5 敲低 / Cy3 对照)，按芯片做中位数归一，3 张芯片做单样本 t 检验。
- **基因集层面**：用双侧 Mann-Whitney U 检验，将每个 marker 基因集的
  变化幅度与其余全部基因比较（见 `geneset_stats.csv`）。

## 结果

### 1）TACSTD2 / Trop-2 —— 依赖具体情境
- **乳腺 CRISPR CLDN4-/-（GSE207704）：TACSTD2 下调**，两株一致——
  log2FC = **-1.00（MCF7）**、**-0.62（T47D）**，均值 **-0.81**（约 1.4–2 倍）。
  即乳腺癌上皮中 CLDN4 缺失会降低 Trop-2。
- **小鼠肺 Cldn4 敲除（GSE50927）：TACSTD2 无变化**——
  `Tacstd2` logFC = -0.14，p = 0.49，FDR = 1.0。
- **卵巢（GSE22493）：无法检验**——该芯片上没有 TACSTD2/Trop-2 探针。

乳腺癌细胞中确实存在 CLDN4→Trop-2 的耦合，但不能推广到正常小鼠肺。

### 2）连接基因 —— 既不崩塌，也无代偿
- **乳腺**：连接基因集出现幅度不大但具统计学意义的协同**下移**
  （中位数 -0.088 vs 背景 +0.012，Mann-Whitney **p = 0.049**，n = 30）。
  单个 claudin/紧密连接基因变化很小（CLDN3 -0.30，CLDN7 -0.14，
  OCLN -0.18，CDH1 -0.05），没有任何 claudin 被代偿性开启。
- **肺**：连接基因集位移幅度极小（中位数 -0.047），虽形式上显著
  （**p = 0.019**，n = 45），但每个连接基因的 FDR = 1.0
  （`Cldn3` +0.09，`Cldn7` +0.11，`Ocln` -0.05，`Tjp1` -0.10，`Cdh1` -0.14），
  紧密连接程序在生物学上保持完整。
- **卵巢**：连接基因集无位移（p = 0.91）。

结论：在三个系统中，CLDN4 缺失在转录层面**都没有**导致紧密连接崩塌，
也没有代偿性的 claudin 切换。

### 3）干扰素/免疫基因 —— 被扰动，且方向随情境反转
这是最强、也最有意思的信号。
- **小鼠肺 Cldn4 敲除：干扰素程序上调**，极显著
  （中位数 +0.39 vs 0，Mann-Whitney **p = 4.8e-13**，n = 69）。核心 ISG
  显著升高：**Isg15 +1.08（FDR 3.1e-3）**、**B2m +0.66（FDR 3.4e-2）**、
  **Stat1 +0.39（p = 0.026）**。CLDN4 缺失的肺处于干扰素/抗原提呈增强态。
- **乳腺 CRISPR CLDN4-/-：干扰素程序下调**，同样极显著
  （中位数 -0.20 vs +0.012，**p = 4.7e-5**，n = 30）；例如 T47D 中 ISG15 -1.87。
- **卵巢 siRNA（低置信度）：干扰素集下调**（中位数 -0.35，**p = 0.025**），
  方向与乳腺一致，但因 CLDN4 探针不稳定，仅弱支持。

即 CLDN4 缺失可重复地**激活/扰动**干扰素-免疫模块，但方向依情境而定：
体内小鼠肺**上调**，体外人乳腺/卵巢癌细胞**下调**。

## 总体结论
1. **TACSTD2/Trop-2**：乳腺癌细胞敲除 CLDN4 后下调约 1.4–2 倍；正常小鼠肺
   无变化；卵巢芯片无法检验。CLDN4–Trop-2 关联存在，但为癌上皮特异性。
2. **连接基因**：任何系统都未见转录层面崩塌，也无代偿性 claudin 诱导，
   最多是小幅协同下移。
3. **干扰素/免疫**：效应最明确——CLDN4 缺失强烈扰动干扰素程序
   （肺敲除上调，乳腺/卵巢癌下调）。

## 注意事项（局限）
- GSE207704 仅提供按组 FPKM（两株细胞充当重复），故细胞株层面的倍数变化
  可靠，但缺少组内方差；分析中将两株视为生物学重复。
- GSE50927 的核心 WT-vs-KO 对比生物学重复很少；ISG 结论依赖多个基因
  一致同向变化（多基因 FDR < 0.05），这可排除单基因噪声，但重复数有限。
- GSE22493 是较老的双色 Operon 芯片，其对照通道是 CLDN4**过表达**
  （而非 scramble），且 CLDN4 探针方向不一致，故仅作辅助证据。

## 复现
```
pip install pandas numpy scipy statsmodels
cd scripts/fable_cldn4_kdko
python3 search_geo.py          # GEO 候选检索
python3 search_broaden.py      # SRA + 宽泛检索
python3 verify_accessions.py   # 样本层面核验
# 按 dataset_inventory.md 下载处理后文件到 results/.../raw/
python3 analyze.py             # 差异表达 + 基因集统计 -> results/*.csv, *.json
```

## Output files / 输出文件
- `results/fable_cldn4_kdko/panel_ifn_apm_results.json` — **directional IFN/MHC-I/APM test（用户优先级）**
- `results/fable_cldn4_kdko/panel_core_by_dataset.csv` — 用户核心面板逐基因跨数据集
- `results/fable_cldn4_kdko/analysis_results.json` — full结果
- `results/fable_cldn4_kdko/geneset_stats.csv` — gene-set Mann-Whitney 统计
- `results/fable_cldn4_kdko/markers_GSE207704.csv` / `_GSE50927.csv` / `_GSE22493.csv`
- `results/fable_cldn4_kdko/GSE22493_CLDN4_kd_diagnostic.csv` — 敲低效率诊断
- `results/fable_cldn4_kdko/geo_verified.json` — 样本级核验记录
- `notes/fable_cldn4_kdko/dataset_inventory.md` — 完整纳入/排除表
