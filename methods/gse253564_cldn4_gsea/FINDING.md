# FINDING — GSE253564 leftover CLDN4 Q4 vs Q1 GSEA (IFN / MHC / TJ) and vs CD274

**Additive only.** Public pre-treatment FPKM from the leftover neoadjuvant **durvalumab ± SBRT** series. This folder does **not** audit or retract the TACSTD2 keratin/TJ leftover (`a8_ici_gsea`) or the MPR/PFS leftover (`opus_geo_leftover`). **CLDN4 only** — TACSTD2 is not used to split samples.

Prerank engine is the same as `methods/cldn4_ko_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (high minus low). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv` and `tables/correlations.tsv`.

---

## 一句话 / TL;DR

| Contrast | Honest n | CLDN4 Q4 vs Q1 | IFN-γ / MHC-I / KEGG TJ NES (FDR) | CLDN4 vs CD274 |
|---|---|---|---|---|
| GSE253564 pre-treatment leftover | **8 vs 8** of 32 | Welch *t* +6.47; mean log2FC +3.352 | Hallmark IFN-γ NES -1.892 FDR 0.001; MHC-I / APM NES -1.264 FDR 0.108; KEGG tight junction NES +2.265 FDR 0.001 | ρ = -0.228 p = 0.209 (n=32) |

GEO deposits **32** pre-treatment tumours. The GSEA contrast is the quartile arms (**8 vs 8**), not n=32. The CD274 Spearman uses every sample with both genes (n=32).

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public author FPKM only. No FASTQ / SRA. |
| File | `GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz` |
| Scale | log2(FPKM+1); duplicate symbols collapsed to the highest-mean locus |
| Split | CLDN4 quartiles: Q1 ≤ 5.044, Q4 ≥ 6.262 (log2(FPKM+1)) |
| Honest n | GEO n=32; Q4 vs Q1 = **8 vs 8** |
| Sign | positive = up in CLDN4 Q4 |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (n=32); Welch *t* after dropping CLDN4 from the rank |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets |
| CD274 | Spearman on all complete pairs; MWU of CD274 in Q4 vs Q1 |
| Out of scope | TACSTD2 splits; GSE248378 post-treatment; MPR/PFS re-test; Hallmark EMT / keratin |

Series: [GSE253564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564). Title on GEO: *Features associated with enhanced proliferation define lung tumors that respond to dual immune checkpoint and radiation therapy*. Platform GPL24676. PMID 38401548.

CLDN4 is a **KEGG tight-junction member**. Q4 vs Q1 is defined on CLDN4, so the TJ NES is partly circular. The drop-CLDN4 sensitivity is the non-circular TJ number.

---

## Headline NES (Welch *t*, Q4 vs Q1)

20128 genes ranked. CLDN4 itself: Welch *t* **+6.47**, mean log2(FPKM+1) difference Q4−Q1 **+3.352** (the split worked).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -1.892 | 0.001 | 0.001 | 197 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -1.264 | 0.108 | 0.108 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +2.265 | 0.001 | 0.001 | 164 |

Hallmark IFN-α (secondary, not in BH): NES -1.231 nom p 0.096 (n=94; not in headline FDR).

Leading-edge (first 25, Welch rank):

- IFN-γ: `PSME2,PIM1,IL15,CASP3,NCOA3,IL15RA,PSMB10,CXCL9,TAP1,IL18BP,CD38,PARP14,IL7,SELP,B2M,FGL2,SLAMF7,RNF31,HLA-DMA,CD274,FAS,ST8SIA4,PARP12,SOCS1,BTG1`
- MHC-I / APM: `PSMB10,ERAP1,TAP1,ERAP2,B2M,HLA-F,PSMB9,IRF1,TAPBPL`
- KEGG TJ: `CLDN4,CRB3,RAC1,MARVELD2,CTTN,EPB41L4B,MYL12B,MYH14,OCLN,F11R,MYL6,CCND1,LLGL2,MAGI1,MYH11,NEDD4,PRKAB2,ROCK2,PARD6B,AMOTL2,PARD3,SCRIB,CLDN7,TJP1,ACTG1`

---

## Sensitivity

Same three headline sets.

| Rank | n used | Result |
|---|---|---|
| Spearman ρ vs continuous CLDN4 | 32 | Hallmark IFN-γ NES -1.839 FDR 0.001; MHC-I / APM NES -1.729 FDR 0.004; KEGG tight junction NES +2.096 FDR 0.001 |
| Welch *t* after dropping CLDN4 from the rank | 8 vs 8 | Hallmark IFN-γ NES -1.891 FDR 0.001; MHC-I / APM NES -1.285 FDR 0.101; KEGG tight junction NES +2.169 FDR 0.001 |

---

## CLDN4 vs CD274 (and companion scores)

Spearman, two-sided, complete cases. IFN-γ and MHC-I scores are the mean of within-cohort gene-wise z-scores (Hallmark IFN-γ genes present; MHC-I = HLA-A/B/C + B2M).

| Pair | n | ρ | p | Note |
|---|---:|---:|---:|---|
| CLDN4 vs CD274 | 32 | -0.228 | 0.209 | primary leftover correlation |
| CLDN4 vs Hallmark IFN-γ score | 32 | -0.213 | 0.241 | companion, not a GSEA substitute |
| CLDN4 vs MHC-I score (HLA-A/B/C + B2M) | 32 | -0.117 | 0.525 | companion |
| CD274 Q4 vs Q1 (MWU) | 8 vs 8 | med 2.902 vs 3.388 | 0.161 | same quartile arms as GSEA |

---

## Extra figures

- `methods/gse253564_cldn4_gsea/figures/fig_headline_nes.png`
- `methods/gse253564_cldn4_gsea/figures/fig_cd274_and_nes.png`

Tables: `methods/gse253564_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `correlations.tsv`, `inventory.tsv`.

---

## What this is not

- Not a TACSTD2 analysis and not a dual-high gate.
- Not a re-run of A8 leftover keratin/TJ/EMT (that split was TACSTD2).
- Not an MPR / PFS / recurrence test (already in `opus_geo_leftover`; CLDN4 vs MPR was null).
- Not GSE248378 (post-treatment residual tumours; 0 MPR in that matrix).
- Not sample-permutation GSEA. Quartile n = 8 vs 8; the engine is gene-set permutation on a prerank.
- Not a claim that CLDN4 *induces* IFN or PD-L1. Association only. Bulk FFPE mixes epithelium and infiltrate.

---

## 中文摘要

只补公开 **GSE253564** 治疗前 leftover durvalumab ± SBRT bulk 的 **CLDN4 Q4 vs Q1** prerank GSEA（IFN / MHC / TJ）以及 CLDN4 vs CD274，不审不撤已有页。不用 TACSTD2 分层。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- GEO n=32。四分位对比是 **8 vs 8**，不要把 n=32 写成 GSEA 的 n。
- Headline（Welch *t*）：Hallmark IFN-γ NES -1.892 FDR 0.001; MHC-I / APM NES -1.264 FDR 0.108; KEGG tight junction NES +2.265 FDR 0.001。
- CLDN4 vs CD274：ρ = -0.228，p = 0.209，n=32。
- CLDN4 在 KEGG TJ 里；去掉 CLDN4 后再跑的 TJ 才是非循环数字。
