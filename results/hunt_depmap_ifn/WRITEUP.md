# DepMap / CCLE lung lines: TACSTD2 vs IFN / MHC-I / APM and vs CLDN4

**Slice:** `notes|scripts|results/hunt_depmap_ifn/` only.
**Data:** DepMap Public 24Q4 (Figshare+ `10.25452/figshare.plus.27993248.v1`), protein-coding RNA `log2(TPM+1)`.
**Primary cohort:** NSCLC cell lines with expression, n = **148** (LUAD 80 / LUSC 27 / other NSCLC 41). SCLC (n=59) is a pre-specified sensitivity cohort, not mixed into the primary test.
**ICI labels:** none. These are cultured lines. There is no T-cell infiltrate and no IFN-γ treatment in this matrix.

**User claim tested:** “TROP2 high suppresses IFN until CLDN4 lost.”

**Verdict (honest):** **Not supported.** In NSCLC lines, TACSTD2 is **positively** associated with basal Hallmark IFN-α/γ scores. TACSTD2 and CLDN4 are tightly co-expressed. The TACSTD2 × CLDN4 interaction on IFN is **null**. TROP2-high / CLDN4-low is rare (n=12) and does **not** have higher IFN than TROP2-high / CLDN4-high. MHC-I is only marginal. This cannot prove or refute immune-cell IFN suppression in tumors.

---

## English

### What this slice can and cannot test

DepMap RNA is **cancer-cell-intrinsic basal transcription**. That is useful: bulk-tumor IFN scores are dominated by immune cells. It is also a hard limit:

- No “hot/cold tumor,” no ICI response, no IFN stimulation time course.
- “CLDN4 lost” here means **low CLDN4 RNA**, not a validated homozygous deletion.
- “Suppresses” is causal language. All tests below are associations.

Pre-specified primary tests (NSCLC, n=148): (1) Spearman TACSTD2 vs Hallmark IFN-γ / IFN-α / IFN-signaling cassette / MHC-I / APM; (2) TACSTD2 vs CLDN4; (3) OLS interaction TACSTD2 × CLDN4 on those signatures (HC3, subtype dummies); (4) TACSTD2–IFN Spearman in CLDN4-high vs CLDN4-low (median split) + Fisher z; (5) median-split quadrants, key contrast TROP2hi/CLDN4hi vs TROP2hi/CLDN4lo; (6) partial Spearman adjusting CLDN4 / EMT / STK11. BH-FDR within predictor among the five primary signatures (+ CLDN4 for TACSTD2).

Signatures (within-cohort z-mean; coverage in `tables/signature_coverage.tsv`):

| Signature | Genes used |
|---|---|
| MHC-I | HLA-A, HLA-B, HLA-C, B2M (4/4) |
| APM | TAP1/2, TAPBP/TAPBPL, PSMB8/9/10, NLRC5, ERAP1/2, CALR, CANX, PDIA3, PSME1/2/3 (16/16) |
| APMS8 | Thompson *JITC* 2020 8-gene set (8/8) |
| IFN signal | STAT1/2, IRF1/9, JAK1/2, IFNGR1/2, IFNAR1/2 (10/10) |
| Hallmark IFN-γ / IFN-α | MSigDB 2024.1 Hs (198/200 and 97/97) |

### TACSTD2 vs CLDN4

NSCLC: Spearman **ρ = 0.682**, p = **1.46×10⁻²¹**, n = 148. LUAD even tighter (ρ = 0.756, n = 80). Both genes anti-correlate with a compact EMT score (TACSTD2 ρ = −0.664; CLDN4 ρ = −0.747). They mark the same epithelial state. Discordant TROP2-high / CLDN4-low lines are uncommon.

### TACSTD2 vs IFN / MHC-I / APM (NSCLC n=148)

Direction is the **opposite** of “TROP2 high suppresses IFN.”

| Target | ρ | p | q (BH) |
|---|---|---|---|
| Hallmark IFN-α | **+0.377** | **2.38×10⁻⁶** | **7.1×10⁻⁶** |
| Hallmark IFN-γ | **+0.344** | **1.81×10⁻⁵** | **3.6×10⁻⁵** |
| APM (16) | +0.185 | 0.024 | 0.036 |
| MHC-I (4) | +0.149 | 0.071 | 0.086 |
| IFN signal cassette | +0.107 | 0.194 | 0.194 |
| APMS8 (Thompson) | −0.026 | 0.758 | — (exploratory) |

Gene-level (exploratory FDR within TACSTD2 tests): IRF1 ρ = +0.243, q = 0.022; HLA-C ρ = +0.216, q = 0.036. HLA-A, B2M, TAP1, STAT1, NLRC5 are not FDR-significant. The Hallmark hit is **ISG-like**, not a clean MHC-I / TAP / STAT1 cassette.

CLDN4 vs the same signatures is also **positive** for IFN-α (ρ = +0.277, q = 0.0033) and IFN-γ (ρ = +0.242, q = 0.0076), weaker than TACSTD2, and null for MHC-I / APM / IFN-signal.

Median-split Mann–Whitney is weaker than Spearman (IFN-γ p = 0.099; IFN-α p = 0.018). Continuous rank correlation is the pre-specified test; do not treat the median-split as a failure of the IFN-α/γ association.

### Does CLDN4 gate a TROP2–IFN suppression?

**No evidence for the claimed gate, and the sign is wrong.**

OLS (z-scored TACSTD2 + CLDN4 + interaction + subtype): IFN-γ interaction β = −0.031, p = 0.51, q = 0.64. IFN-α, MHC-I, APM interactions all p > 0.29. TACSTD2 main effect on IFN-γ is **positive** (β = +0.115, p = 0.016).

CLDN4-median strata (74 / 74): TACSTD2 vs IFN-γ ρ = **+0.426** (p = 1.5×10⁻⁴) in CLDN4-high and **+0.357** (p = 0.0018) in CLDN4-low. Fisher z p = **0.63**. Both strata are positive. The claim predicted a negative association that disappears after CLDN4 loss.

Quadrants (median × median):

| Quadrant | n | IFN-γ median |
|---|---|---|
| TROP2hi CLDN4hi | 62 | +0.016 |
| TROP2hi CLDN4lo | **12** | −0.049 |
| TROP2lo CLDN4hi | 12 | −0.035 |
| TROP2lo CLDN4lo | 62 | −0.099 |

Kruskal p = 0.34. Key contrast TROP2hi/CLDN4hi vs TROP2hi/CLDN4lo: MW p = 0.48. The off-diagonal cells are small **because the two genes co-vary**. The claim needs a populated TROP2-high / CLDN4-lost group; DepMap NSCLC does not provide one.

Partial Spearman, IFN-γ: TACSTD2 | CLDN4 ρ = +0.182, p = 0.027 (still positive). CLDN4 | TACSTD2 ρ = +0.106, p = 0.20 (null). Adjusting EMT **strengthens** TACSTD2–IFN-γ (ρ = +0.451). Adjusting CLDN4+STK11 remains positive (ρ = +0.169, p = 0.040).

### Subtype sensitivity

| Cohort | n | TACSTD2 vs IFN-γ ρ | p | TACSTD2 vs CLDN4 ρ |
|---|---|---|---|---|
| LUAD | 80 | +0.268 | 0.016 | +0.756 |
| LUSC | 27 | +0.274 | 0.167 | +0.532 |
| SCLC | 59 | +0.374 | 0.0035 | +0.348 |

LUSC is underpowered. SCLC still goes **positive**, not negative. Mixing SCLC into “all malignant lung” (n=208) inflates TACSTD2–IFN-γ (ρ = +0.51 range after EMT residualization) because SCLC sits at the low-TACSTD2 / low-IFN end. That is why NSCLC was primary.

**LUAD-only exploratory (not the primary test):** CLDN4-high (n=40) TACSTD2–IFN-γ ρ = +0.540, p = 3.2×10⁻⁴; CLDN4-low ρ = +0.045, p = 0.78; Fisher p = 0.016. After residualizing CLDN4, TACSTD2–IFN-γ in LUAD is null (ρ = +0.030, p = 0.79). This is a CLDN4-shared epithelial/IFN coupling, **still positive**, and the TROP2hi/CLDN4lo cell in LUAD is n = 5. Do not rewrite the claim around this stratum.

### STK11 (positive-control confounder)

STK11 hotspot-or-damaging (28 / 120): lower TACSTD2 (p = 0.019, q = 0.022), lower CLDN4 (p = 0.011, q = 0.016), lower IFN-γ / IFN-α / MHC-I / APM (all q ≤ 0.012). That is the expected STK11–IFN-low axis. It tracks with **lower**, not higher, TACSTD2. KEAP1 (n=13) trends the same way for APM (q = 0.046).

### How this sits next to Bessede 2024 and tumor slices

Bessede et al. (*Clin Cancer Res*, PMID 38048058) reported high TACSTD2 with atezolizumab resistance and fewer T cells in OAK/POPLAR. CPTAC LUAD protein in this repo is consistent with the **infiltration** half (high TACSTD2 protein, lower xCell immune / CD8). DepMap answers a different question: do TROP2-high **cancer cells** themselves transcribe less IFN / MHC-I / APM at baseline?

**They do not.** If anything, TROP2-high / CLDN4-high epithelial NSCLC lines have **higher** basal Hallmark IFN scores. A tumor-level “TROP2 high, IFN/T-cell low” pattern is therefore more likely **microenvironment / lineage mix / STK11**, not a cell-intrinsic TROP2→IFN transcriptional block that CLDN4 loss releases.

### Caveats

1. Basal culture RNA, not IFN-treated, not protein, not surface MHC.
2. Hallmark IFN-γ includes genes a pure cancer line may not use (e.g. GZMA, XCL1). The IFN-α set is cleaner and agrees in sign.
3. MHC-I and TAP/NLRC5 are weak or null; do not say “APM is up with TROP2” as a strong claim (q = 0.036 for the 16-gene mean; APMS8 is null).
4. TROP2hi/CLDN4lo n=12 (LUAD n=5): the gating contrast is underpowered by construction.
5. Mutation calls are hotspot OR damaging; not copy-number loss of CLDN4/TACSTD2.
6. Multiple testing: primary-list FDR is the inferential claim; gene-level and LUAD-stratum tables are exploratory.
7. No CRISPR effect of TACSTD2 or CLDN4 on IFN was tested. Expression correlation ≠ suppression.

### How to rerun

```bash
python3 -m pip install pandas numpy scipy matplotlib statsmodels
python3 scripts/hunt_depmap_ifn/00_download.py --outdir data/hunt_depmap_ifn
python3 scripts/hunt_depmap_ifn/01_analyze.py --data data/hunt_depmap_ifn --outdir results/hunt_depmap_ifn
```

Figures: `results/hunt_depmap_ifn/figures/fig1_tacstd2_vs_cldn4.png` … `fig8_distributions.png`.
Tables: `results/hunt_depmap_ifn/tables/*.tsv` and `key_stats.json`.

---

## 中文

### 能测什么、不能测什么

DepMap 是**癌细胞本身的基线转录**（无免疫浸润、无 IFN 刺激、无 ICI 标签）。适合问“TROP2 高的肺癌细胞自己有没有更低的 IFN/MHC-I/APM”，**不能**问肿瘤热/冷或免疫治疗耐药。所谓 “CLDN4 lost” 在这里只是 **CLDN4 RNA 低**，不是经验证的纯合缺失。“抑制”是因果措辞，下面全是相关。

预设主队列：NSCLC n=148（LUAD 80 / LUSC 27 / 其他 NSCLC 41）。SCLC n=59 只做敏感性，不混进主检验。

### 结论（如实）

**用户假设不成立。** NSCLC 里 TACSTD2 与 Hallmark IFN-α/γ **正相关**，不是负相关。TACSTD2 与 CLDN4 共表达很强（ρ=0.682）。交互项不显著。TROP2 高 / CLDN4 低的细胞系很少（n=12），IFN 并不比 TROP2 高 / CLDN4 高更高。MHC-I 只有边缘相关。

### 主结果（NSCLC n=148）

- TACSTD2 vs CLDN4：ρ = **0.682**，p = **1.46×10⁻²¹**
- TACSTD2 vs Hallmark IFN-α：ρ = **+0.377**，q = **7.1×10⁻⁶**
- TACSTD2 vs Hallmark IFN-γ：ρ = **+0.344**，q = **3.6×10⁻⁵**
- TACSTD2 vs APM（16 基因）：ρ = +0.185，q = 0.036（弱）
- TACSTD2 vs MHC-I：ρ = +0.149，p = 0.071，q = 0.086（不显著）
- IFN 信号盒（STAT/JAK/IRF/受体）：ρ = +0.107，p = 0.19（无）
- Thompson 8 基因 APMS：ρ = −0.026（无）
- 交互 TACSTD2×CLDN4 → IFN-γ：β = −0.031，p = 0.51
- CLDN4 中位数分层：高/低两组 TACSTD2–IFN-γ 都是正（ρ=+0.43 / +0.36），Fisher p=0.63
- 四象限 Kruskal p=0.34；关键对比 TROP2hi/CLDN4hi vs TROP2hi/CLDN4lo p=0.48
- 偏相关：TACSTD2 | CLDN4 对 IFN-γ 仍为正（ρ=+0.18，p=0.027）；CLDN4 | TACSTD2 不显著

STK11 突变（28/120）伴随 **更低** 的 TACSTD2、CLDN4 和 IFN/MHC/APM（阳性对照方向正确）。TROP2 低、而不是 TROP2 高，和这条 IFN 低轴走在一起。

LUAD 单独看，CLDN4 高组里 TACSTD2–IFN-γ 更正（ρ=+0.54 vs +0.05，Fisher p=0.016），校正 CLDN4 后 TACSTD2 效应消失。这仍是**同向共变**，不是“TROP2 抑制 IFN、CLDN4 丢失后解除”。且 LUAD 的 TROP2hi/CLDN4lo 只有 5 株。

与 Bessede 2024 / CPTAC 蛋白切片的关系：那些说的是肿瘤里 T 细胞少，不是癌细胞基线 IFN 转录低。DepMap 说明 **细胞本底并不支持 “TROP2 高抑制 IFN、直到 CLDN4 丢失”**。肿瘤水平的 IFN/T 细胞低，更可能来自微环境、谱系混合或 STK11，而不是这条细胞自主开关。
