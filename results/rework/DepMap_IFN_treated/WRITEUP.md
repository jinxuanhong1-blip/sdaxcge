# DepMap IFN-treated rework: PRISM, CRISPR IFN-gene effects, protein ISGs, MHC-I protein

**Slice:** `notes|scripts|results/rework/DepMap_IFN_treated/` only.
**Why this slice exists:** the unstimulated test is a **mismatch**, not a confirmation. In DepMap 24Q4 NSCLC culture RNA (n = 148), TACSTD2 vs Hallmark IFN-α was Spearman **ρ = +0.377**, p = 2.38×10⁻⁶ (hunt: `results/hunt_depmap_ifn`). That is **positive**, not negative. That RNA test is **closed**. This folder does not recompute Hallmark IFN-α/γ RNA scores.

**Reworked question:** if TROP2-high NSCLC cells “suppress IFN,” the signal should appear in a layer the basal transcriptome cannot see — IFN-pathway **treatment**, IFN-pathway **dependency**, or IFN / MHC-I **protein**.

**Data**

| Layer | Source | What it actually is |
|---|---|---|
| Models + TACSTD2 RNA | DepMap Public 24Q4, Figshare+ `10.25452/figshare.plus.27993248.v1` | Unstimulated `log2(TPM+1)`. Used only as the *predictor*. |
| PRISM viability | Repurposing Public 23Q2, `10.6084/m9.figshare.23600310.v4` | 5-day small-molecule LFC at ~2.5 µM. |
| CRISPR gene effect | DepMap 24Q4 Chronos | Knockout fitness, not expression. |
| Protein | Nusinow et al., *Cell* 2020 CCLE MS (Gygi) | TMT abundance, not surface MHC-I FACS. |

**Primary cohort:** NSCLC with TACSTD2 RNA, n = **148** (LUAD 80 / LUSC 27 / other NSCLC 41). Same type rules as the unstimulated hunt. Overlap: PRISM n = **99**, CRISPR n = **97**, protein n = **65**. SCLC (n = 59 RNA; n = 25 CRISPR) is sensitivity only.

**Honest gap — say it first:** there is **no** public DepMap IFN-α/γ–treated transcriptome, and recombinant IFN biologics are **not** in PRISM 23Q2 (or the extended primary list, 6,658 IDs). “IFN-treated DepMap/PRISM” in public data means **IFN-pathway drugs** (JAK inhibitors; IFN inducers / IFNR agonists), not IFN protein on the cells.

**Verdict:** **Does not rescue the claim.** PRISM JAK / IFN-inducer viability vs TACSTD2 is **null**. CRISPR IFN type-I and APM gene-effect means vs TACSTD2 are **null** (nominal single-gene hits do not survive BH). MHC-I **protein** vs TACSTD2 RNA or TACSTD2 protein is **null**. The only protein association that is even nominally consistent is **positive** (TACSTD2 protein vs ISG protein ρ = +0.363, p = 0.003) — the **same sign** as the +0.377 RNA mismatch, opposite “TROP2 suppresses IFN.”

---

## English

### What changed conceptually

The original cell-intrinsic claim predicted **negative** TACSTD2–IFN coupling. Unstimulated Hallmark IFN-α RNA went the other way (+0.377). Three interpretations were still open:

1. The relevant IFN is **induced**, not basal — need IFN-treated viability.
2. The relevant phenotype is **dependency** on IFN sensing / APM, not ISG transcription.
3. The relevant MHC-I is **protein**, because RNA MHC-I was already weak (ρ = +0.149, p = 0.071 in the closed hunt).

This slice tests those three. It does **not** re-open Hallmark RNA.

### 1. PRISM: no recombinant IFN; pathway drugs are null

PRISM is a small-molecule multiplexed viability screen (Corsello et al., *Nat Cancer* 2020). Interferon alfa / gamma proteins are not library members. The public proxies, pre-specified:

- **JAK-inhibitor core** (ruxolitinib, tofacitinib, baricitinib, fedratinib, itacitinib, upadacitinib, AZD1480, CYT387, plus tool compounds in `01_analyze.py`).
- **IFN inducers / IFNR agonists** (imiquimod, tilorone, bropirimine, MIW-815, pidotimod, propagermanium).

NSCLC with both TACSTD2 and PRISM LFC, n = 99. More negative LFC = more killing.

| Test | n | ρ | p |
|---|---|---|---|
| Median LFC, JAK-inhibitor core vs TACSTD2 | 99 | **+0.031** | 0.76 |
| Median LFC, IFN inducer / IFNR agonist vs TACSTD2 | 99 | **−0.050** | 0.62 |

Every single focus compound is null after BH within the focus set (all q ≈ 0.99). Strongest raw JAK hits are WHI-P154 (ρ = −0.137, p = 0.20) and NVP-BSK805 (ρ = −0.136, p = 0.20). Ruxolitinib ρ = +0.054, p = 0.61. Imiquimod ρ = +0.018, p = 0.86. None of these sit in the extreme tail of the 6,415-drug background.

**Read:** TROP2-high NSCLC lines are not more resistant (or more sensitive) to public PRISM IFN-pathway drugs. This is **not** an IFN-α/γ treatment experiment. Do not write “IFN-treated PRISM confirmed suppression.”

### 2. CRISPR: IFN / APM gene effects vs TACSTD2 are null

Chronos more negative = more dependent. If TROP2-high cells suppressed IFN sensing, they should be **less** dependent on IFNAR/JAK/STAT/APM (ρ vs TACSTD2 **positive**). If they relied on a basal IFN program, they should be **more** dependent (ρ **negative**). Either coherent signed pattern would have been interesting. We do not see one.

Primary NSCLC CRISPR overlap, n = 97. BH-FDR among IFN type-I / type-II / negative-regulator / APM genes.

| Target (Chronos) | ρ vs TACSTD2 RNA | p | q |
|---|---|---|---|
| IFN type-I mean (IFNAR1/2, TYK2, JAK1, STAT1/2, IRF9) | **+0.077** | 0.45 | — (pre-specified mean) |
| APM mean (B2M, TAP1/2, TAPBP, NLRC5, HLA-A/B/C) | **+0.017** | 0.87 | — |
| JAK1 | −0.212 | 0.037 | **0.38** |
| TYK2 | +0.223 | 0.028 | **0.38** |
| TAP1 | +0.200 | 0.050 | **0.38** |
| STAT1 | −0.046 | 0.66 | 0.92 |
| IFNAR1 | +0.131 | 0.20 | 0.63 |
| B2M | +0.040 | 0.70 | 0.92 |
| HLA-A / HLA-B / HLA-C | +0.009 / −0.125 / +0.089 | all p ≥ 0.22 | all q ≥ 0.63 |

JAK1 and TYK2 are **opposite signs** at the same nominal threshold. That is not an IFN-dependency axis. After BH, nothing is significant. TACSTD2 itself is not a lung essential (median Chronos = −0.051). SCLC (n = 25 CRISPR) is also null for STAT1 / IFNAR1 / IFNGR1 / B2M / JAK1 (all p > 0.31).

**Read:** TROP2-high NSCLC lines do not systematically depend more or less on IFN sensing or APM. Expression correlation is not CRISPR suppression.

### 3. Protein IFN / ISG scores: same sign as the RNA mismatch, weaker

CCLE MS (Nusinow 2020), NSCLC with protein, n = 65. ISG protein score = within-cohort z-mean of detected ISGs (STAT1/2, IRF9, MX1, ISG15, IFIT1/3, OAS1/2/3, HLA-A/B/C, B2M, TAP1/2, PSMB8/9/10, and others listed in `tables/key_stats.json`). TACSTD2 RNA vs TACSTD2 protein ρ = **0.836** (n = 65) — the protein assay sees TROP2.

| Test | n | ρ | p | q (within MHC-I + ISG proteins) |
|---|---|---|---|---|
| ISG protein score vs TACSTD2 **RNA** | 65 | **+0.202** | 0.11 | — |
| ISG protein score vs TACSTD2 **protein** | 65 | **+0.363** | **0.0030** | — (score, not in the per-protein FDR) |
| MX1 protein vs TACSTD2 RNA | 65 | +0.309 | 0.012 | 0.32 |
| STAT2 protein vs TACSTD2 RNA | 65 | +0.250 | 0.045 | 0.44 |
| STAT1 protein vs TACSTD2 RNA | 65 | +0.009 | 0.95 | 0.99 |
| ISG15 protein vs TACSTD2 RNA | 65 | −0.037 | 0.77 | 0.99 |

The protein–protein ISG association is the strongest new number in this folder. It is **positive**. It agrees with unstimulated Hallmark IFN-α RNA (+0.377), not with “TROP2 suppresses IFN.” Per-protein FDR does not call MX1 or STAT2. Do not upgrade +0.363 into a mechanism.

### 4. MHC-I protein, not RNA: null

This was the cleanest leftover test. The closed hunt already showed MHC-I **RNA** was only marginal (ρ = +0.149, p = 0.071). Protein could have gone negative. It does not.

| Test | n | ρ | p |
|---|---|---|---|
| MHC-I protein score (HLA-A, HLA-B, HLA-C, B2M) vs TACSTD2 **RNA** | 65 | **+0.073** | 0.57 |
| MHC-I protein score vs TACSTD2 **protein** | 65 | **+0.172** | 0.17 |
| HLA-A protein vs TACSTD2 RNA | 65 | +0.036 | 0.78 |
| HLA-B protein vs TACSTD2 RNA | 65 | +0.028 | 0.83 |
| HLA-C protein vs TACSTD2 RNA | 65 | +0.119 | 0.34 |
| B2M protein vs TACSTD2 RNA | 65 | −0.003 | 0.98 |
| Median-split TACSTD2 RNA → MHC-I protein (MW) | 32 / 33 | — | **0.98** |

**Read:** TROP2-high NSCLC lines do not have less MHC-I protein. They also do not have more. The protein layer does not flip the RNA mismatch into the claimed negative.

CLDN4 protein still tracks TACSTD2 RNA (ρ = +0.592, n = 37, p = 1.2×10⁻⁴). The epithelial co-expression from the closed hunt is still there. It is not an IFN gate.

### How this sits next to the +0.377 mismatch

| Layer | Sign vs TACSTD2 | Supports “TROP2 suppresses IFN”? |
|---|---|---|
| Unstimulated Hallmark IFN-α RNA (closed) | **+0.377** | No — opposite |
| PRISM JAK / IFN-inducer viability | ~0 | No |
| CRISPR IFN type-I / APM Chronos | ~0 | No |
| ISG **protein** score | **+** (protein–protein 0.36; RNA–protein 0.20) | No — same sign as RNA |
| MHC-I **protein** | ~0 | No |

A tumor-level “TROP2 high, T-cell / IFN low” pattern (Bessede 2024; CPTAC infiltration slices in this repo) is still a **microenvironment / lineage / STK11** problem. It is not a cell-intrinsic TROP2 → IFN transcriptional block, not an IFN-pathway drug-resistance trait in PRISM, not an IFN-gene CRISPR dependency, and not low MHC-I protein in CCLE MS.

### Caveats

1. **No recombinant IFN treatment.** PRISM proxies are small molecules. Imiquimod is a TLR7 agonist, not IFN-α. JAK inhibitors hit other kinases. A true IFN-α/γ dose–response on these 148 lines does not exist in DepMap public files.
2. Chronos is fitness, not IFN-reporter activity. A line can fail to induce MHC-I after IFN-γ and still be Chronos-null for IFNGR1 in rich media.
3. Nusinow MS is whole-cell TMT, n = 65 NSCLC, missingness is not random. It is not surface HLA-I.
4. Hallmark IFN-γ RNA includes immune genes a pure line may not use. That is why this rework left Hallmark RNA alone and moved to protein / CRISPR / PRISM.
5. Nominal JAK1 / TYK2 / MX1 / STAT2 p-values are **exploratory**. Primary inference is the pre-specified means / scores and the MHC-I protein tests. Those are null or positive.
6. No CRISPR of TACSTD2 measuring IFN output. No CLDN4-loss × IFN-treatment interaction. The closed hunt already showed the RNA interaction was null and TROP2hi/CLDN4lo was n = 12.
7. SCLC CRISPR n = 25 is underpowered; reported only as a sign check.

### How to rerun

```bash
python3 -m pip install pandas numpy scipy matplotlib statsmodels
python3 scripts/rework/DepMap_IFN_treated/00_download.py --outdir data/rework/DepMap_IFN_treated
python3 scripts/rework/DepMap_IFN_treated/01_analyze.py \
  --data data/rework/DepMap_IFN_treated \
  --outdir results/rework/DepMap_IFN_treated
```

Figures: `figures/fig1_tacstd2_vs_IFN_type1_Chronos.png` … `fig8_protein_forest.png`.
Tables: `tables/*.tsv`, `tables/key_stats.json`.

---

## 中文

### 为什么重做

未刺激 DepMap NSCLC 里 TACSTD2 对 Hallmark IFN-α **不是负相关**，而是 **ρ = +0.377**。那个 RNA 检验已经结束，这里不再算 Hallmark。剩下还能问的是：IFN **处理**、IFN 通路 **依赖**、IFN/MHC-I **蛋白**。

### 公开数据缺口（先说清楚）

DepMap 没有 IFN-α/γ 处理过的转录组。PRISM 23Q2 扩展库里也没有干扰素蛋白，只有 JAK 抑制剂和 IFN 诱导剂 / IFNR 激动剂。不能写成“IFN 处理的 DepMap/PRISM 证实了抑制”。

### 结论（如实）

**救不回“TROP2 高抑制 IFN”**。

- PRISM：JAK 核心中位 LFC vs TACSTD2，ρ = **+0.031**，p = 0.76（n=99）；IFN 诱导剂 ρ = **−0.050**，p = 0.62。单个化合物 FDR 全空。
- CRISPR（n=97）：I 型 IFN 基因效应均值 ρ = **+0.077**，p = 0.45；APM 均值 ρ = **+0.017**，p = 0.87。JAK1（ρ=−0.21）和 TYK2（ρ=+0.22）名义显著但方向相反，BH 后 q=0.38。
- ISG **蛋白** vs TACSTD2 RNA：ρ = **+0.202**，p = 0.11；vs TACSTD2 **蛋白**：ρ = **+0.363**，p = 0.003。方向和 +0.377 RNA **相同**，不是负的。
- MHC-I **蛋白**（HLA-A/B/C + B2M）vs TACSTD2 RNA：ρ = **+0.073**，p = 0.57；vs TACSTD2 蛋白：ρ = +0.172，p = 0.17。中位数分割 p = 0.98。四个蛋白单独都是空。

肿瘤里 TROP2 高、T 细胞/IFN 低，更可能是微环境、谱系或 STK11，不是这株癌细胞自己的 IFN 开关。
