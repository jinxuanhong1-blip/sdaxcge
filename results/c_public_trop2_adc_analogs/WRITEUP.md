# EXTRA analog figures — public TROP2-ADC / TROP2-loss expression series

**Slice:** `notes|scripts|methods|results/c_public_trop2_adc_analogs/`
**Private C3–C5 SKB264 mechanism (CLDN4 internalization, shared IFN/MHC-I with CLDN4-KD, TJ down in PDX) is taken as given and is not re-derived.**
**These series are public analogs. None of them is SKB264 / sac-TMT / MK-2870.**

---

## English

### What was tested

GEO search (2026-08-16) for sacituzumab, govitecan, IMMU132, Trodelvy, SKB264, sac-TMT, tirumotecan, datopotamab. **Zero public GSE records for sac-TMT / SKB264 / datopotamab.** IMMU132 hits that actually deposit an ADC-versus-control matrix are GSE312098, GSE311016, and GSE304294 (sacituzumab govitecan). Candidate accessions named in the request were opened and classified before scoring.

| Accession | Deposited contrast | ADC vs control? | Tissue | Call |
|---|---|---|---|---|
| **GSE312098** | CX-1 ± IMMU132, 2 d (also GSK / combo) | **yes, n=3 vs 3** | CRC cell line | primary SG analog |
| **GSE311016** | 5 CRC PDX ± IMMU132, day 29 | **yes, 5 pairs** | CRC PDX | primary SG analog |
| **GSE304294** | KYSE30 ± IMMU132, 1 d (also IACS / combo) | **yes, n=2 vs 3** | ESCC | SG analog, underpowered |
| GSE334497 | 4T1 Trop2 KO vs WT tumors | no | mouse TNBC | TROP2-loss analog |
| GSE289287 | T-47D Trop2 KO xenografts | no | luminal breast | TROP2-loss analog (not TNBC) |
| GSE235812 | untreated breast P0 vs matching PDX P1 | no | breast | untreated correlation only |
| GSE245459 | SKOV3 shTACSTD2 ± cisplatin | no | ovarian | already scored in `results/w200/C4_GSE245459/` |
| GSE302284 | residual NSCLC / osimertinib vs vehicle | no | lung | not an SG treatment contrast |
| GSE278664 | HGSOC biopsies, treatment = prexasertib | no | ovarian | SG paper, wrong treatment |

Pre-specified sets (frozen in `scripts/c_public_trop2_adc_analogs/gene_sets.py`): TACSTD2, CLDN4; C4 IFN/MHC-I panel `IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A`; APM; broader IFN; junction/barrier. Statistics on `log2(x+1)`: Welch t (unpaired) or paired t (PDX); BH-FDR genome-wide is exploratory; gene-set Mann–Whitney of log2FC vs expressed background plus a sign binomial. FPKM is not counts. Combination arms were not used to rescue a null primary.

Support for a **public ADC analog of the given SKB264 direction** would be CLDN4 / junction down and IFN/MHC-I/APM up after TROP2-ADC versus control. Anything else is reported as such.

### EXTRA Fig 1–2 — GSE312098 CX-1 IMMU132 (best public ADC analog)

IMMU132 = sacituzumab govitecan. **Not SKB264.** Non-lung CRC line, 2 days, n=3 vs 3.

| Gene | mean FPKM ctrl | IMMU132 | log2FC | p | q |
|---|---:|---:|---:|---:|---:|
| TACSTD2 | 41.5 | 67.2 | **+0.68** | 5.7×10⁻⁴ | 0.015 |
| **CLDN4** | 91.8 | 50.1 | **−0.86** | **1.7×10⁻⁵** | **0.0043** |
| CLDN1 | 6.03 | 10.1 | +0.66 | 9.3×10⁻⁵ | 0.0076 |
| CLDN7 | 47.2 | 42.7 | −0.14 | 0.064 | 0.21 |
| IFI27 | 29.8 | 45.6 | +0.61 | 0.011 | 0.069 |
| OAS2 | 4.27 | 8.98 | +0.94 | 0.025 | 0.12 |
| IFIT1 | 5.70 | 8.42 | +0.49 | 0.0065 | 0.051 |
| MX1 | 0.20 | 0.17 | −0.038 | 0.50 | 0.63 |
| ISG15 | 115 | 179 | +0.63 | 0.0033 | 0.034 |
| HLA-A | 146 | 154 | +0.078 | 0.14 | 0.37 |
| B2M | 74.1 | 133 | +0.83 | 4.6×10⁻⁵ | 0.0061 |
| TAP1 | 3.25 | 4.70 | +0.43 | 0.070 | 0.23 |

C4 panel: 5/6 up (MX1 unexpressed), median log2FC +0.55 vs background, MW **p = 0.012**. APM 15/19 up, MW **p = 0.0056**. IFN 44/10 up, MW **p = 2.6×10⁻⁸**. Junction set as a whole is **null** (14 up / 13 down, MW p = 0.99) because CLDN4 down is offset by CLDN1/CDH1-class genes up. Sample z-scores: C4 panel Δ = +1.26, p = 0.0025; APM Δ = +0.85, p = 0.0068; IFN Δ = +0.98, p = 0.020; junction Δ = −0.016, p = 0.91.

**Honest read:** in this SG-treated CRC line, CLDN4 RNA falls and IFN/APM rise. That is a public analog in the same *direction* as the given SKB264 claim. It is still SG, still CX-1, still 2-day FPKM, still not internalization or protein.

### EXTRA Fig 3 — GSE311016 CRC PDX IMMU132 (paired, day 29)

Five PDX models, control vs IMMU132, paired t. **Not SKB264.** Non-lung.

CLDN4 mean FPKM 129 → 98, log2FC **−0.44**, paired p = **0.060**, q = 0.60. CLDN7 −0.26, p = 0.033, q = 0.60. TACSTD2 −0.42, p = 0.27. C4 panel 5/6 up but median +0.18, MW p = 0.18 (null). Broader IFN set MW p = 0.011 (up); APM and junction sets null. All four z-scores are NS (C4 p = 0.70; junction p = 0.14). ISG15 looks up in some models and down in others (paired log2FC +1.26, p = 0.32).

**Honest read:** CLDN4 trends down after 29-day SG in CRC PDX and does not reach p<0.05. IFN is not a clean on-treatment opening at the C4 panel. Do not pool with GSE312098.

### EXTRA Fig 4 — GSE304294 KYSE30 IMMU132 (n=2 vs 3, 1 day)

**Underpowered. Not SKB264.** Non-lung ESCC.

CLDN4 **rises** (74 → 139 FPKM, log2FC **+0.91**, p = 3.9×10⁻⁵, q = 0.045). TACSTD2 +1.16, p = 0.0099. Junction set is up (MW p = 0.0013), the opposite of TJ-down. IFIT1 is down (−0.61, p = 2.9×10⁻⁴); ISG15 +0.21, p = 0.036; HLA-A +0.62, p = 0.011. APM set 17/19 up, MW p = 3.4×10⁻⁶. C4 panel MW p = 0.11 (null).

**Honest read:** this 1-day ESCC SG analog does **not** reproduce CLDN4/TJ down. APM goes up. n=2 in the IMMU arm. Do not average it with CX-1.

### EXTRA Fig 5 — TROP2-loss analogs (not ADC)

**GSE334497** 4T1 mouse TNBC, Trop2 KO vs WT, n=5 vs 5. Tacstd2 log2FC **−3.82**, p = 0.0011 (KO is real). Cldn4 −0.82, p = 0.25. Cldn1 −2.18, p = 0.035. Junction set 3 up / 20 down, MW **p = 1.8×10⁻⁵** (down). C4 IFN panel 3/3, MW p = 0.56 (null). Broader IFN MW p = 0.022 (up); APM MW p = 0.037 (up). Junction z-score Δ = −0.64, p = 0.023; IFN z-score NS.

**GSE289287** T-47D luminal xenografts, Trop2 KO n=4 vs WT n=3. Author DESeq2 TACSTD2 log2FC **−3.26**, p = 1.6×10⁻⁶⁵. CLDN4 author log2FC **+0.28**, p = 0.13; Welch on normCounts +0.29, p = 0.041, q = 0.51. C4 panel **6/6 up**, median +0.99, MW **p = 2.1×10⁻⁴**, binomial p = 0.031. Broader IFN MW p = 5.3×10⁻⁹. Junction set null. This is **not TNBC**.

**Honest read:** TROP2 loss is not a substitute for TROP2-ADC. 4T1 KO lowers a junction program without a clean Cldn4 or IFN call. T-47D KO opens IFN/ISG and does **not** lower CLDN4.

### EXTRA Fig 6 — GSE235812 untreated breast (no ADC contrast)

Deposited RNA is untreated patient tumor (P0) and matching PDX (P1). The companion paper used SG in other experiments; **this matrix has treatment = no treatment**. 34 / 36 htseq libraries parsed from `GSE235812_RAW.tar`.

TACSTD2 vs CLDN4 Spearman ρ = **0.57**, p = 4.6×10⁻⁴, q = 0.0018, n = 34. CLDN7 ρ = 0.67; CLDN1 ρ = 0.64; OCLN ρ = 0.65; F11R ρ = 0.64 (all q<0.001). C4 IFN genes vs TACSTD2 are all NS (ISG15 ρ = 0.29, p = 0.097; HLA-A ρ = 0.32, p = 0.067).

**Cannot** be scored as after TROP2-ADC.

### EXTRA Fig 7 and leftovers

Gene-set forest (Fig 7) is ADC-only. GSE245459 is the existing SKOV3 shTACSTD2 C4 analog (IFN/APM **closed**, not opened) and was not re-run. GSE302284 is osimertinib / residual tissue, not SG vs control (left to the dedicated GSE302284 slice). GSE278664 is a prexasertib HGSOC trial. No public sac-TMT RNA.

### Verdict

| Question | Public analog answer |
|---|---|
| Is there public SKB264 / sac-TMT RNA vs control? | **No.** |
| After SG (IMMU132), does CLDN4 fall? | **Yes in CX-1 (GSE312098).** Trend in CRC PDX (p=0.060). **Opposite** in KYSE30 n=2. |
| After SG, does IFN/MHC-I/APM rise? | **Yes in CX-1.** Mixed/NS in CRC PDX. APM up, C4 panel NS in KYSE30. |
| Do these data equal SKB264 C3–C5? | **No.** Different drug, different tissues, RNA not internalization, no lung ADC-vs-control matrix. |

---

## 中文

### 测了什么

2026-08-16 检索 GEO：sacituzumab / govitecan / IMMU132 / Trodelvy / SKB264 / sac-TMT / tirumotecan / datopotamab。**没有** sac-TMT / SKB264 / datopotamab 的公开 GSE。真正存了 ADC 对对照矩阵的 IMMU132 系列是 GSE312098、GSE311016、GSE304294（均为 sacituzumab govitecan）。用户点名的候选号均先打开、分类，再打分。

私有 C3–C5 SKB264 机制（CLDN4 内化、与 CLDN4-KD 共享的 IFN/MHC-I、PDX 中 TJ 下降）**当作已知**，此处不重做。本切片只提供**额外公开类似图**。IMMU132 ≠ SKB264。

预设基因集与统计见 `scripts/c_public_trop2_adc_analogs/gene_sets.py` 与 `methods/c_public_trop2_adc_analogs/playbook.md`。报告 n / log2FC / p。FPKM 不是 counts。n=2 标为效力不足。

### 主类似：GSE312098 CX-1 IMMU132（非肺）

CLDN4 log2FC **−0.86**（p=1.7×10⁻⁵，q=0.0043）。C4 面板 5/6 上调，MW p=0.012；APM MW p=0.0056；IFN MW p=2.6×10⁻⁸。连接基因集整体为阴性（CLDN4 降、CLDN1 升）。这是公开 SG 类似，方向与给定的 SKB264 说法一致，**不是** SKB264，也不是内化/蛋白。

### GSE311016 CRC PDX（5 对，第 29 天）

CLDN4 log2FC **−0.44**，配对 p=**0.060**。C4 面板 MW 阴性。宽 IFN 集 MW p=0.011，z 分数不显著。不能与 CX-1 合并。

### GSE304294 KYSE30（n=2 对 3，1 天）

CLDN4 **升高** log2FC **+0.91**（p=3.9×10⁻⁵）。连接集上调。APM 上调。效力不足。**不支持** TJ 下降。

### TROP2 缺失（不是 ADC）

GSE334497（4T1 TNBC）：Tacstd2 确降；连接集下降（MW p=1.8×10⁻⁵）；Cldn4 本身 p=0.25；IFN 面板阴性。GSE289287（T-47D，非 TNBC）：CLDN4 不降（作者 DESeq2 +0.28，p=0.13）；C4 面板 6/6 上调。

### GSE235812

未处理乳腺 P0/P1。TACSTD2–CLDN4 ρ=0.57（n=34，p=4.6×10⁻⁴）。IFN 基因与 TACSTD2 无显著相关。**不能**当作 ADC 之后。

### 结论

公开世界里**没有** SKB264 / sac-TMT 的处理对对照 RNA。SG（IMMU132）在 CX-1 上出现 CLDN4 下降 + IFN/APM 上升；CRC PDX 为趋势；ESCC n=2 方向相反。这些是额外类似，不是 C3–C5。

---

## Reproduce

```bash
python3 scripts/c_public_trop2_adc_analogs/download.py
python3 scripts/c_public_trop2_adc_analogs/analyze.py
```

Figures: `results/c_public_trop2_adc_analogs/figures/fig_extra1_*.png` … `fig_extra7_*.png`.
Numbers: `key_stats.json`, `tables/key_genes_all_contrasts.tsv`, `tables/geneset_stats.tsv`.
