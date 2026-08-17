# FINDING — GSE166449 CLDN4 Q4 vs Q1 prerank GSEA (IFN / MHC / TJ and vs CD274)

**Additive only. CLDN4 only. Honest n = 22.** Public pretreatment advanced LUAD pembrolizumab bulk (Lee / Hwang et al., *JITC* 2021, PMID 33857424; GEO [GSE166449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166449)). This folder does **not** audit or retract `results/w200/GSE166449`. TACSTD2 is **not** an anchor here.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `methods/gse285029_cldn4_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

**Honest n.** The deposited matrix has **22** pretreatment tumours (7 responder / 15 non-responder from GEO titles). CLDN4 quartiles are **Q4 = 6 vs Q1 = 6**. That is the smallest quartile contrast this repo has previously allowed. Welch *t* on 6 vs 6 is underpowered. The Spearman-vs-CLDN4 rank uses all 22 samples and is the complementary read, not a rescue.

CD274 is a locked **comparator**, not a second claim gene. CLDN4 vs CD274 Spearman on n=22 is ρ = -0.020 (p = 0.930). The same three sets are also ranked by CD274 Q4 vs Q1 (6 vs 6).

Do not headline response OR. The sibling GSE166449 page already found no persuasive CLDN4–response association (n=7 / 15).

---

## 一句话 / TL;DR

| Contrast | Honest n | Split | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE166449 CLDN4 Q4 vs Q1 | 22 (Q4=6, Q1=6) | CLDN4 quartiles on deposited log2(TPM+1) | Hallmark IFN-γ NES +1.292 FDR 0.024; MHC-I / APM NES +1.087 FDR 0.345; KEGG tight junction NES +1.703 FDR 0.003 |

Primary rank is Welch *t* (Q4 − Q1). Hallmark IFN-γ is FDR 0.024 on this 6 vs 6 split and **NS** on the Spearman-vs-CLDN4 rank that uses all 22 (NES +0.747 FDR 0.949). Do not quote the quartile IFN-γ FDR as a continuous-CLDN4 fact. MHC-I / APM is NS on both ranks. CLDN4 is a member of KEGG tight junction, so the primary TJ NES includes the stratification gene. After dropping CLDN4 from the rank: NES +1.686 FDR 0.003 nom p 0.001 (n=164). TJ stays up on Spearman (NES +1.742 FDR 0.003).

CLDN4 vs CD274 (all 22): Spearman ρ = -0.020, p = 0.930. CD274 Welch *t* on the CLDN4 Q4 vs Q1 split = +0.231. CD274 Q4 vs Q1 GSEA (comparator): Hallmark IFN-γ NES +2.368 FDR 0.001; MHC-I / APM NES +2.208 FDR 0.001; KEGG tight junction NES +1.307 FDR 0.007. CD274 is itself a Hallmark IFN-γ member, so the CD274-split IFN-γ NES is partly circular. The point of that row is the contrast with CLDN4, not a new PD-L1 GSEA claim.

This is bulk ICI RNA on n=22. An IFN / MHC NES can be infiltrate. Do not write a cell-intrinsic call.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public author TPM only. No FASTQ / SRA. |
| File | `GSE166449_Raw_gene_TPM_matrix.txt.gz` |
| Samples | **n = 22** pretreatment LUAD (GEO immunotherapy; paper = pembrolizumab) |
| Response labels | 7 Responder / 15 Non-responder from GEO titles (inventory only; not a headline) |
| Genes | 20345 symbols on the deposited matrix |
| Transform | **none** — deposited values are already log2(TPM+1) |
| Anchor | **CLDN4 only** |
| Split | Q4 vs Q1 on CLDN4 (≤P25 vs ≥P75). Ties stay in the tail. |
| n Q4 / Q1 | 6 / 6 |
| CLDN4 log2 cuts | P25 = 1.0059; P75 = 1.8684 |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (all 22); Welch *t* after dropping CLDN4 |
| Comparator | CD274: Spearman vs CLDN4 (n=22); CD274 Q4 vs Q1 GSEA on the same three sets |
| TJ honesty | Drop CLDN4 from the rank and re-run the three headline sets |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets, **per rank** |
| Out of scope | TACSTD2 ranks; headline response OR; FASTQ; sample-permutation GSEA |

Series: [GSE166449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166449). Paper: Lee / Hwang et al., *J Immunother Cancer* 2021; PMID [33857424](https://pubmed.ncbi.nlm.nih.gov/33857424/).

---

## Headline NES (Welch *t*, CLDN4 Q4 − Q1)

19669 genes ranked. Honest n = **6 vs 6**. CLDN4 itself is the stratification gene (Welch *t* = +7.489; by construction the top of the rank).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +1.292 | 0.024 | 0.016 | 198 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +1.087 | 0.345 | 0.345 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.703 | 0.003 | 0.001 | 165 |

Hallmark IFN-α (secondary, not in BH): NES +1.527 nom p 0.004 (n=95; not in headline FDR).

Leading-edge (first 25, Welch rank):

- IFN-γ: `PLA2G4A,CFB,PTGS2,VAMP8,BATF2,LGALS3BP,LY6E,NFKBIA,ICAM1,UPP1,METTL7B,USP18,TNFSF10,CASP7,SRI,SECTM1,IFITM3,PLSCR1,STAT3,ST3GAL5,NAMPT,IRF7,IRF9,RIPK1,ZNFX1`
- MHC-I / APM: `CALR,CANX,PDIA3,TAPBP,HLA-G,ERAP1,HLA-E,TAPBPL,B2M,HLA-C,HLA-B,HLA-A`
- KEGG TJ: `CLDN4,MARVELD3,CRB3,CTTN,F11R,ERBB2,CCND1,CGN,LLGL2,TJP3,SCRIB,MYH14,MARVELD2,CLDN3,OCLN,PARD6B,EZR,SLC9A3R1,CLDN7,TJP1,CLDN10,NEDD4L,RAC1,ARPC1A,AMOTL2`

---

## Sensitivity: drop CLDN4 from the rank

Same Q4 vs Q1 samples and the same three headline sets. CLDN4 is removed from the ranked list so KEGG TJ is not scored on the gene used to define the split.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +1.292 | 0.016 | 0.011 | 198 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +1.068 | 0.358 | 0.358 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.686 | 0.003 | 0.001 | 164 |

---

## Sensitivity: Spearman ρ vs continuous CLDN4

Uses all **n = 22** samples. Rank = per-gene Spearman ρ vs CLDN4. BH-FDR is again inside the three headline sets. This does not throw away the middle 10 tumours.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +0.747 | 0.949 | 0.949 | 198 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +0.814 | 0.926 | 0.617 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.742 | 0.003 | 0.001 | 165 |

---

## Comparator: vs CD274

CD274 is present on the deposited matrix. This is **not** a PD-L1 IHC result.

| Item | Value |
|---|---|
| n | 22 |
| CLDN4 vs CD274 Spearman ρ | -0.020 |
| two-sided p | 0.930 |
| CD274 Welch *t* on CLDN4 Q4 vs Q1 | +0.231 |
| CD274 Q4 / Q1 n | 6 / 6 |
| CD274 log2 cuts | P25 = 1.7648; P75 = 3.4404 |

Honest read: on n=22 the CLDN4–CD274 correlation is negative and not significant. Do not treat CD274 as a surrogate for CLDN4 on this series. CD274 is a Hallmark IFN-γ member, so the CD274 Q4 IFN-γ NES is partly circular; MHC-I / APM is the cleaner CD274-high comparator (NES +2.208 FDR 0.001).

**CD274 Q4 vs Q1 GSEA** (same three headline sets; BH inside the three; positive NES = up in CD274 Q4):

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +2.368 | 0.001 | 0.001 | 198 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +2.208 | 0.001 | 0.001 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.307 | 0.007 | 0.007 | 161 |

---

## Extra figures

- `methods/gse166449_cldn4_gsea/figures/fig_headline_nes.png`
- `methods/gse166449_cldn4_gsea/figures/fig_cldn4_vs_cd274.png`

Tables: `methods/gse166449_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `inventory.tsv`, `quartile_split.tsv`, `cd274_comparator.tsv`.

---

## What this is not

- Not a TACSTD2 analysis. CLDN4 only.
- Not a re-run or retraction of `results/w200/GSE166449`.
- Not a response-OR headline. That test is already NS on n=7 / 15.
- Not FASTQ / salmon / DESeq2 from SRA (author TPM is used as deposited).
- Not sample-permutation GSEA. The engine is gene-set permutation on a prerank.
- Not a cell-intrinsic IFN / MHC / PD-L1 call. This is bulk RNA on n=22.
- Not a large-n result. Q4 vs Q1 is 6 vs 6. Quote the n with the NES.

---

## 中文摘要

只补公开 **GSE166449**（诚实 n=22）上 **CLDN4 Q4 vs Q1** 的 prerank GSEA（IFN-γ / MHC-I / KEGG TJ），并报告与 **CD274** 的对照。不审不撤已有 `results/w200/GSE166449` 页。只做 CLDN4，不做 TACSTD2。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- 22 例 ICI 前 LUAD 肿瘤 bulk。Q4=6，Q1=6。6 vs 6 是本仓库四分位 GSEA 的下限，功效不足。
- Headline（Welch *t*，6 vs 6）：Hallmark IFN-γ NES +1.292 FDR 0.024; MHC-I / APM NES +1.087 FDR 0.345; KEGG tight junction NES +1.703 FDR 0.003。
- IFN-γ 在 Spearman（n=22）上是 NS（NES +0.747 FDR 0.949），不要把四分位 IFN-γ FDR 写成连续 CLDN4 事实。
- CLDN4 是 KEGG tight junction 成员；去掉 CLDN4 后再跑：NES +1.686 FDR 0.003 nom p 0.001 (n=164)。
- CLDN4 vs CD274（n=22）：ρ = -0.020，p = 0.930。CD274 高组 IFN/MHC 强，CLDN4 不是它的替代。
