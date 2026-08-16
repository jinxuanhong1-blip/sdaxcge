# C4 analog — GSE245459 SKOV3 shTACSTD2: CLDN4 and IFN/APM

**Slice:** `notes|scripts|results/w200/C4_GSE245459/`
**Dataset:** [GSE245459](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245459) (Han et al., *J Cancer* 2024; SKOV3 bulk RNA-seq FPKM).
**Claim C4 (private):** CLDN4 KD **opens** IFN/MHC-I (`IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A`).
**This slice:** TACSTD2 (TROP2) shRNA analog in **ovarian SKOV3**, not a CLDN4-KD test and not lung. Primary contrast is **shTACSTD2 vs shNC, no cisplatin**.

---

## English

### What was tested

GEO supplementary `GSE245459_fpkm.anno.txt.gz` (HTTP 200, 14 MB). Twelve SKOV3 libraries, n = 3 per group:

| Group | GEO titles | Perturbation |
|---|---|---|
| shNC | shNC1–3 | scramble |
| shTACSTD2 | sh1–3 | shTACSTD2 |
| shNC+DDP | shNCDDP1–3 | scramble + cisplatin 10 µg/ml, 48 h |
| shTACSTD2+DDP | shDDP1–3 | shTACSTD2 + cisplatin |

Author design (GEO + paper): TACSTD2 in platinum resistance via Rap1/PI3K/AKT. They did **not** claim an IFN/APM phenotype. We do not re-analyse SRA counts.

Pre-specified (frozen in `scripts/w200/C4_GSE245459/gene_sets.py` before looking at fold-changes):

1. Knockdown QC: TACSTD2.
2. C4 junction gene: CLDN4.
3. C4 IFN/MHC-I panel (6 genes).
4. APM (19 genes) and a broader IFN/ISG set (55).
5. Junction set (context) and Rap1/PI3K/AKT transcripts (paper axis; not C4).

Statistics on `log2(FPKM+1)`: Welch t-test (n=3 vs n=3), BH-FDR genome-wide; gene-set Mann–Whitney of log2FC vs expressed background (mean FPKM ≥ 1); sample scores = mean of per-gene z across the 12 libraries. FPKM is not counts. n=3 is small. Do not treat q-values as RNA-seq limma/edgeR.

Primary contrast **KD_noDDP**. Cisplatin arms are secondary. C4 support in this analog would require **CLDN4 down and IFN/APM up**. Anything else is reported as such.

### Knockdown QC — TACSTD2 is down

| Group | mean FPKM |
|---|---|
| shNC | 5.40 |
| shTACSTD2 | **0.034** |
| shNC+DDP | **19.88** |
| shTACSTD2+DDP | 0.276 |

sh vs shNC: log2FC = **−2.63**, p = **5.3×10⁻⁵**, q = **0.004**. On cisplatin: log2FC = **−4.03**, p = **4.4×10⁻⁵**, q = **0.007**. The shRNA is visible at RNA. Cisplatin **induces** TACSTD2 in scramble cells (log2FC = **+1.69**, p = 0.0034, q = 0.022).

### CLDN4 — down after TACSTD2 KD (no cisplatin)

| Group | mean FPKM |
|---|---|
| shNC | 3.15 |
| shTACSTD2 | **0.101** |
| shNC+DDP | 0.365 |
| shTACSTD2+DDP | 0.239 |

sh vs shNC: log2FC = **−1.92**, p = **3.3×10⁻⁵**, q = **0.0037**. In this line, TACSTD2 loss takes CLDN4 RNA with it (~30-fold). That is the only part of the analog that matches a TACSTD2→CLDN4 story.

It does **not** isolate CLDN4: CLDN3 (−1.48, q=0.005) and CLDN7 (−1.39, q=0.014) fall too. EPCAM also falls (−1.25, q=0.006). TJP1 and OCLN do not. Junction-set median log2FC = −0.23 vs background −0.09 (MW p = 0.007) — a coordinated dip, not a scaffold collapse.

On cisplatin, CLDN4 is already low in shNC+DDP (log2FC vs shNC = **−1.60**, p = 2.2×10⁻⁵). Adding TACSTD2 KD does not move it further (log2FC = −0.14, p = 0.099). Cisplatin **dissociates** TACSTD2 (up) from CLDN4 (down).

### IFN / APM — closed, not opened (primary contrast)

This is the C4 test. Direction is the **opposite** of “opening.”

**C4 panel, shTACSTD2 vs shNC (all 6 genes down):**

| Gene | FPKM shNC | FPKM sh | log2FC | p | q |
|---|---|---|---|---|---|
| IFI27 | 0.45 | 0.084 | −0.42 | 0.022 | 0.084 |
| OAS2 | 0.47 | 0.021 | −0.52 | 0.017 | 0.072 |
| IFIT1 | 15.8 | 3.50 | **−1.89** | 0.0029 | **0.024** |
| MX1 | 0.37 | 0.16 | −0.25 | 0.032 | 0.110 |
| ISG15 | 62.3 | 40.3 | −0.62 | 0.014 | 0.064 |
| HLA-A | 18.6 | 0.090 | **−4.17** | 1.0×10⁻⁷ | **6.7×10⁻⁴** |

Panel vs background: median log2FC −0.57 vs −0.09, MW **p = 0.0029**; binomial 6/6 down, p = 0.031. Sample z-mean: shNC +1.18 vs sh −0.78, Welch p = **6.9×10⁻⁴**.

**APM:** 15/19 genes down; median log2FC −0.40 vs −0.09, MW **p = 0.006**. B2M −0.51 (q=0.023), TAP1 −0.40 (q=0.040), TAPBP −1.28 (q=0.006), PSMB8 −1.54 (q=0.010), HLA-C −0.40 (q=0.028). TAP2 is the clear exception (**+0.70**, q=0.040). HLA-B −0.35 (q=0.050).

**IFN/ISG:** 50/54 present genes down; median −0.44 vs −0.09, MW **p = 5.1×10⁻⁹**. STAT1 −0.49 (q=0.037). CXCL10 is low-expressed and slightly down.

HLA-A (ENSG00000206503) collapses from ~19 to ~0.09 FPKM in all three sh replicates. MHC class I annotation is messy; the replicate concordance is real in this matrix. Do not silently drop it, and do not treat a single HLA locus as the whole APM story — the rest of the APM/ISG set moves the same way.

**Verdict on C4:** TACSTD2 KD in untreated SKOV3 **does not open** IFN/MHC-I or APM. It **closes** them while CLDN4 falls. That is evidence **against** a simple “lose CLDN4 → open IFN/APM” reading in this system. It is consistent with the public SKOV-3 CLDN4-siRNA array (GSE22493) in the parallel CLDN4-KD sweep, where the IFN set also went down (low-confidence array). Breast CRISPR CLDN4−/− (GSE207704) was also IFN-down. The only public CLDN4-loss IFN-up call in that sweep was mouse lung Cldn4 KO (GSE50927). SKOV3 matches the cancer-cell-line sign, not the mouse-lung sign.

### Cisplatin interaction — do not swap this in as the C4 result

Cisplatin in scramble cells **also** closes IFN/APM (C4 panel 6/6 down, MW p=0.003; APM MW p=4.8×10⁻⁶; IFN MW p=7.6×10⁻⁸) while raising TACSTD2 and lowering CLDN4.

Relative to that suppressed shNC+DDP baseline, shTACSTD2+DDP scores sit higher (C4 panel z-mean delta +0.92, p=0.003; APM +0.84, p=0.019; IFN +0.79, p=0.014). Gene-set MW on KD_DDP: C4 panel **null** (p=0.39), APM up p=0.039, IFN up p=2.0×10⁻⁴. CLDN4 is not further reduced.

This is a **drug × knockdown interaction**, not C4. Possible readings: (i) TACSTD2-high cells under cisplatin keep IFN/APM lower; (ii) scoring artifact from a low DDP-in-NC baseline; (iii) death / cell-state mix after 48 h cisplatin. The paper’s WGCNA used shNC / shNC+DDP / sh+DDP and emphasised platinum resistance, not IFN. None of this rescues the primary no-drug C4 analog.

Rap1/PI3K/AKT transcripts are a **null** gene-set shift in KD_noDDP (MW p=0.53). The paper’s mechanism is signalling/protein, not this mRNA module. We do not confirm or refute Rap1/AKT here.

### Bottom line

| Pre-specified question | Honest call |
|---|---|
| Is TACSTD2 knocked down? | **Yes** (RNA). |
| Does CLDN4 fall? | **Yes**, no-drug arm (q=0.0037). Other claudins fall too. |
| Does IFN/MHC-I / APM open? | **No. They close.** Opposite of C4. |
| Does this replicate C4? | **No.** Analog support would need CLDN4 down **and** IFN/APM up. We have the first and the opposite of the second. |
| Cisplatin arm | Interaction exists; **not** the C4 test. |

### Caveats

- Ovarian SKOV3, not lung ICI tissue. Analog of a mechanism, not a cohort replication.
- Perturbation is TACSTD2 shRNA, not CLDN4 KD. CLDN4 and TACSTD2 are confounded.
- Author FPKM, not re-quantified counts. No limma/voom. n=3.
- CLDN4 basal FPKM is modest (~3). HLA-A drop is extreme; MHC mapping caveat stands.
- No protein, TEER, MHC surface, or T-cell killing in this GEO deposit.
- Paper question is cisplatin IC50 / Rap1-PI3K-AKT, not IFN.

### Reproduce

```
pip install pandas numpy scipy matplotlib statsmodels
python3 scripts/w200/C4_GSE245459/download.py
python3 scripts/w200/C4_GSE245459/analyze.py
```

Outputs: `results/w200/C4_GSE245459/tables/` and `figures/`.

---

## 中文

### 测了什么

公开 GEO 补充文件 `GSE245459_fpkm.anno.txt.gz`。SKOV3 十二个文库，每组 n=3：shNC、shTACSTD2（GEO 标题 sh1–3）、shNC+顺铂、shTACSTD2+顺铂（10 µg/ml，48 h）。作者论文做的是 TACSTD2 经 Rap1/PI3K/AKT 介导铂耐药，**没有**声称 IFN/APM 表型。

C4 私有主张：CLDN4 敲低**打开** IFN/MHC-I。本切片是卵巢 SKOV3 的 **TACSTD2 shRNA 类比**，不是 CLDN4 敲低，也不是肺。主对比：**无顺铂的 shTACSTD2 vs shNC**。支持 C4 类比需要 **CLDN4 下降且 IFN/APM 上升**。

方法：`log2(FPKM+1)`，Welch t，全基因组 BH-FDR；基因集对表达背景做 Mann–Whitney；样本分为基因 z 均值。FPKM 不是 counts，n=3，不能当成 limma。

### 敲低质控

TACSTD2：shNC 5.40 → sh 0.034 FPKM，log2FC = **−2.63**，q = **0.004**。顺铂背景下同样被打掉（log2FC −4.03）。乱码对照里顺铂反而**诱导** TACSTD2（log2FC +1.69）。

### CLDN4

无药：3.15 → 0.101 FPKM，log2FC = **−1.92**，q = **0.0037**。这条线上 TACSTD2 丢失会带走 CLDN4 RNA。CLDN3/CLDN7/EPCAM 也降，TJP1/OCLN 不明显，不能说“只动了 CLDN4”。顺铂本身就把 CLDN4 压低（log2FC −1.60），再敲 TACSTD2 几乎不再降（p=0.099）。顺铂下 TACSTD2↑、CLDN4↓，二者解耦。

### IFN / APM（主对比）——关上，不是打开

C4 六基因**全部下调**（IFI27/OAS2/IFIT1/MX1/ISG15/HLA-A）。IFIT1 log2FC −1.89（q=0.024）；HLA-A −4.17（q=6.7×10⁻⁴）。面板相对背景 MW **p=0.0029**。APM 15/19 下调，MW **p=0.006**。IFN/ISG 50/54 下调，MW **p=5.1×10⁻⁹**。

**结论：在未加药的 SKOV3 里，TACSTD2 敲低并不打开 IFN/MHC-I 或 APM，而是关闭。** 这与“丢掉 CLDN4 → 打开 IFN/APM”的简单读法相反。方向与平行切片里 SKOV-3 CLDN4 siRNA（GSE22493，低置信芯片）以及乳腺 CRISPR CLDN4−/−（GSE207704）的 IFN 下调一致；公开数据里 IFN 上调主要出现在小鼠肺 Cldn4 KO（GSE50927）。

### 顺铂交互 —— 不能替换成 C4 结果

乱码+顺铂本身也关闭 IFN/APM。相对这条被压低的基线，shTACSTD2+DDP 的 z 分更高（C4 面板 p=0.003），但 C4 六基因 MW 为 **null**，且 CLDN4 不再降。这是药×敲低交互，不是 C4。Rap1/PI3K/AKT 转录模块在主对比上为 null（MW p=0.53），不能用来证实或否定作者的蛋白信号轴。

### 一句话

敲低成立；无药时 CLDN4 下降；IFN/APM **下降而不是上升**。C4 类比**不成立**。顺铂臂有交互，不能拿来救 C4。

### 限制

卵巢不是肺；扰动是 TACSTD2 不是 CLDN4；作者 FPKM；n=3；无蛋白/TEER/杀伤；HLA 注释需保留但 MHC 位点本身乱。
