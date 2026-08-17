# FINDING — GSE68465 CLDN4 Q4 vs Q1 prerank GSEA (IFN / MHC / TJ)

**Additive only. CLDN4 only.** Public Director's Challenge LUAD microarray (Shedden et al., *Nat Med* 2008, PMID [18641660](https://pubmed.ncbi.nlm.nih.gov/18641660/); GEO [GSE68465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465)). This folder does **not** audit or retract `methods/gse68465_cldn4` (CLDN4 vs CD8A after ESTIMATE purity). TACSTD2 is **not** an anchor here. No slide was re-scored. Unit is the **tumour array**.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `methods/gse285029_cldn4_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

GEO deposits **462** arrays (**443** `Lung Adenocarcinoma` + **19** `Normal`). The paper text says 442 LUAD. Honest tumour n here is the GEO filter: **443**. The GSEA contrast is the quartile arms (**111 vs 111**), not n=443 and not n=462. No ICI labels. This is bulk U133A — an IFN / MHC NES can be infiltrate, not a tumour-cell-intrinsic call.

---

## 一句话 / TL;DR

| Contrast | Honest n | Split | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE68465 CLDN4 Q4 vs Q1 | **111 vs 111** of 443 LUAD | Welch *t* on log2(MAS5) | Hallmark IFN-γ NES -2.223 FDR 0.001; MHC-I / APM NES -1.840 FDR 0.001; KEGG tight junction NES +1.795 FDR 0.001 |

Primary rank is Welch *t* (Q4 − Q1) on log2(MAS5). CLDN4 is a member of KEGG tight junction, so the primary TJ NES includes the stratification gene. After dropping CLDN4 from the rank: NES +1.649 FDR 0.001 nom p 0.001 (n=143). U133A coverage is incomplete; n in rank is the intersection.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public GEO series matrix only. No CEL / MAS5 reprocess. |
| File | `GSE68465_series_matrix.txt.gz` |
| Platform | GPL96 Affymetrix HG-U133A |
| Arrays | **n = 462** (22,283 probes) |
| Tumours used | **n = 443** `disease_state = Lung Adenocarcinoma` |
| Held out | **19** `Normal` (Stratagene) |
| Paper vs GEO | Shedden text = 442 LUAD; GEO filter = **443**. This page uses 443. |
| Genes after max-mean | 13101 (first `///` HUGO from official GPL96) |
| Transform | `log2(MAS5)` (deposited values are MAS5 intensity, min>0) |
| Anchor | **CLDN4 only** (`201428_at`; named probe `201428_at`) |
| Split | Q4 vs Q1 on CLDN4 (≤P25 vs ≥P75). Ties stay in the tail. |
| Honest GSEA n | Q4 = 111, Q1 = 111 (not 443, not 462) |
| CLDN4 log2 cuts | P25 = 9.1168; P75 = 10.1974 |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (all 443 LUAD); Welch *t* after dropping CLDN4 from the rank |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| U133A aliases | `WARS1→WARS`, `MARCHF1→MARCH1`, `PALS1→MPP5` (documented HGNC updates only) |
| FDR | BH inside the 3 headline sets, **per rank** |
| Out of scope | TACSTD2 ranks; ICI response; CEL reprocess; sample-permutation GSEA; CD8 / ESTIMATE re-run |

Series: [GSE68465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465). Title on GEO: *caArray_jacob-00182: gene expression–based survival prediction in lung adenocarcinoma*. Paper: Shedden et al., *Nat Med* 2008; PMID [18641660](https://pubmed.ncbi.nlm.nih.gov/18641660/).

U133A intersection (after aliases): Hallmark IFN-γ 179/200 on U133A (177 direct, 2 alias, 21 absent); MHC-I / APM 20/21 on U133A (20 direct, 0 alias, 1 absent); KEGG tight junction 144/169 on U133A (143 direct, 1 alias, 25 absent); Hallmark IFN-α 81/97 on U133A (80 direct, 1 alias, 16 absent). `CD274` and `NLRC5` are **absent** on this platform.

CLDN4 itself is a KEGG TJ member and is the split gene. Its Welch *t* on this contrast is **+28.42**. TJ NES on the full rank is therefore partly circular. The drop-CLDN4 row is the non-circular companion.

---

## Headline NES (Welch *t*, Q4 − Q1)

13101 genes ranked. **111 Q4 vs 111 Q1** of 443 LUAD. CLDN4 itself is the stratification gene (Welch *t* = +28.421; by construction the top of the rank).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -2.223 | 0.001 | 0.001 | 179 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -1.840 | 0.001 | 0.001 | 20 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.795 | 0.001 | 0.001 | 144 |

Hallmark IFN-α (secondary, not in BH): NES -1.763 nom p 0.001 (n=81; not in headline FDR).

Leading-edge (first 25, Welch rank):

- IFN-γ: `CFH,TRIM25,CIITA,FPR1,CMKLR1,SAMHD1,PSMB10,IRF2,MT2A,RAPGEF6,UBE2L6,IFIH1,NMI,CD69,IL7,ST8SIA4,APOL6,B2M,MARCH1,SLC25A28,ARID5B,IL15RA,P2RY14,BANK1,IL2RB`
- MHC-I / APM: `HLA-E,HLA-F,CANX,TAP2,PSMB8,ERAP1,PSMB10,B2M,PSMB9,TAPBPL,ERAP2,PDIA3,IRF1,TAP1`
- KEGG TJ: `CLDN4,CLDN3,ERBB2,OCLN,EZR,CLDN7,F11R,JUN,TJP3,RAC1,ARPC1B,ARPC1A,DLG3,LLGL2,EPB41L4B,PRKCZ,CLDN9,STK11,SCRIB,PARD3,ACTN4,PRKAB1,MYH9,CDK4,CLDN15`

---

## Sensitivity: drop CLDN4 from the rank

Same Q4 vs Q1 samples and the same three headline sets. CLDN4 is removed from the ranked list so KEGG TJ is not scored on the gene used to define the split.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -2.237 | 0.001 | 0.001 | 179 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -1.810 | 0.002 | 0.002 | 20 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.649 | 0.001 | 0.001 | 143 |

---

## Sensitivity: Spearman ρ vs continuous CLDN4

Uses all **n = 443** LUAD arrays. Rank = per-gene Spearman ρ vs CLDN4. BH-FDR is again inside the three headline sets. This is **not** the GSEA n — the quartile contrast remains 111 vs 111.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -2.199 | 0.001 | 0.001 | 179 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -1.994 | 0.001 | 0.001 | 20 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.758 | 0.001 | 0.001 | 144 |

---

## Extra figures

- `methods/gse68465_cldn4_gsea/figures/fig_headline_nes.png`

Tables: `methods/gse68465_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `inventory.tsv`, `quartile_split.tsv`, `geneset_coverage.tsv`.

---

## What this is not

- Not a TACSTD2 analysis. CLDN4 only.
- Not a re-run or retraction of `methods/gse68465_cldn4` (CD8A / ESTIMATE purity).
- Not n=462 and not the paper's 442. Tumour n = 443; GSEA n = 111 vs 111.
- Not CEL / RMA / fRMA from RAW.tar. The deposited MAS5 series matrix is used as deposited.
- Not sample-permutation GSEA. The engine is gene-set permutation on a prerank.
- Not a cell-intrinsic IFN/MHC call. This is bulk U133A surgical LUAD.
- Not an ICI / RECIST / PD-L1 analysis. This series has no ICI labels.
- Not a claim that CLDN4 *induces* IFN or MHC. Association only.

---

## 中文摘要

只补公开 **GSE68465** Director's Challenge LUAD 芯片上 **CLDN4 Q4 vs Q1** 的 prerank GSEA（IFN / MHC / TJ），不审不撤已有 `gse68465_cldn4`（CD8 / 纯度）页。只做 CLDN4，不做 TACSTD2。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- GEO 462 张芯片：443 LUAD + 19 Normal。论文写 442；本页用 GEO 过滤 **n=443**。GSEA 对比是四分位臂 **111 vs 111**，不要把 443 写成 GSEA 的 n。
- Headline（Welch *t*，log2 MAS5）：Hallmark IFN-γ NES -2.223 FDR 0.001; MHC-I / APM NES -1.840 FDR 0.001; KEGG tight junction NES +1.795 FDR 0.001。
- CLDN4 是 KEGG tight junction 成员；去掉 CLDN4 后再跑：NES +1.649 FDR 0.001 nom p 0.001 (n=143)。
- U133A 覆盖不全（Hallmark IFN-γ 179/200 on U133A (177 direct, 2 alias, 21 absent); MHC-I / APM 20/21 on U133A (20 direct, 0 alias, 1 absent); KEGG tight junction 144/169 on U133A (143 direct, 1 alias, 25 absent); Hallmark IFN-α 81/97 on U133A (80 direct, 1 alias, 16 absent)）。`CD274` / `NLRC5` 不在此平台。
