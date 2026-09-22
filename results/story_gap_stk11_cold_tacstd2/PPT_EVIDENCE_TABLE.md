# STORY GAP FILL — Disease context slide

**Claim to support (not to invent):** KL / STK11 tumors are immunologically cold, and Tacstd2/Cldn4 are high in that disease context.

**Method rule:** Numbers copied from existing public-analysis PRs + published Skoulidis 2018 only. No re-fit. No fabrication. Bulk vs scRNA labeled on every row.

---

## Slide headline (honest)

> **Mouse cell-line bulk locks KL > KP Tacstd2/Cldn4 (GSE137244).**  
> **STK11/KL is cold / ICI-resistant in human literature + public MSK.**  
> **Human TCGA bulk does *not* show STK11 → high TACSTD2/CLDN4 — it shows the reverse.**  
> Do not paste human bulk as the Tacstd2/Cldn4-high limb.

---

## PPT paste table

| # | Limb | Assay class | Dataset | Number to show | Supports story? | One-line caveat |
|---|---|---|---|---|---|---|
| 1 | **Cold / ICI-resistant** | Literature (human) | Skoulidis 2018 SU2C KRAS-mut LUAC | KL ORR **7.4%** vs KP **35.7%** vs K-only **28.6%** (P<0.001) | **YES** | Cite paper; not recomputed here |
| 2 | Cold TME | Literature (human IHC) | Skoulidis PROSPECT | STK11mut lower CD3 (P=0.0019) and CD8 (P=0.0072) | **YES** | Published IHC; no local matrix |
| 3 | Cold / ICI (public reanalysis) | Bulk genotype→DCB | MSK Hellmann 2018 + Rizvi 2018 (KRAS-restricted) | STK11 co-mut DCB OR **0.24** / **0.38**; pooled MH OR **0.36** | **YES (direction)** | Small n; Fisher p 0.19 / 0.089; PR #49 |
| 4 | **Tacstd2/Cldn4 high** | **Mouse cell-line bulk** | **GSE137244** KL n=5 vs KP n=5 | Tacstd2 **Δ=+3.238**; Cldn4 **Δ=+5.570**; exact MW **p=0.00794**; complete KL>KP | **YES — LOCK** | Cell-line RNA-seq, not tumor TME; handoff lock |
| 5 | Tacstd2 high (support) | Mouse tumor **bulk** | GSE179500 Lkb1XTR (p53-WT) | LKB1-off > restored: author padj **0.032**; reanalysis MW **p=0.003**, n=10 vs 7 | **YES (modest)** | Restore design; KL-vs-never-lost KT ns (p=0.45); PR #113 |
| 6 | Cldn4 high (support) | Mouse tumor **scRNA** | GSE179502 sorted neoplastic | Cldn4 mouse mean 0.300 → 0.143 on Lkb1 restore; complete rank sep; Welch p=0.036 | **YES (direction)** | Unit = mouse n=3 vs 3; MW floor p=0.1; PR #530 |
| 7 | Tacstd2 in KL epi | Mouse **scRNA** | GSE179502 | 54.7% Tacstd2+; ~14% high subset; NonRestored vs Restored mouse MW **p=0.10** | **PARTIAL** | Underpowered genotype; immune depleted — cannot test cold; PR #113 |
| 8 | Tacstd2 compartment | Mouse **scRNA** | GSE180963 K vs KL | KL epi Tacstd2+ **62.8–64.7%** vs immune ~4%; immune frac ~82% both | **SELECTIVITY only** | **n=1 mouse/genotype**; cold **not** reproduced as fewer immune cells; PR #137/#113 |
| 9 | **Human STK11 → TACSTD2/CLDN4** | Human tumor **bulk** | TCGA-LUAD | STK11mut **lower** TACSTD2 log2FC **−0.51** (FDR 1e-4); CLDN4 **−0.47** (FDR ~3e-5); KL subset −0.71 / −0.68 | **NO — reverse** | Do not claim human bulk STK11 = high targets; PR #49 |
| 10 | Human STK11 TACSTD2 | Human **scRNA** | GSE280232 neoadjuvant ICB | TACSTD2 epithelial (72–85% in EPCAM+); STK11mut vs wt epi **not powered** (n≈2 vs 1 rich) | **NO genotype claim** | Mostly sorted T cells (~0% TACSTD2); PR #113 |

---

## Two-column slide layout (recommended)

**Left — Cold disease context**

- Skoulidis: KL ORR 7.4%; PROSPECT CD8↓  
- Public MSK KRAS-restricted STK11 DCB OR ~0.24–0.38 (pooled MH 0.36)

**Right — Tacstd2/Cldn4 high (where it actually locks)**

- **GSE137244 cell-line bulk:** Tacstd2 +3.24, Cldn4 +5.57, MW p=0.00794  
- Support: GSE179500 bulk LKB1-off Tacstd2↑; GSE179502 Cldn4↓ on restore  
- Footer: **Human TCGA bulk = STK11 → lower TACSTD2/CLDN4** (honest gap)

---

## Bulk vs scRNA cheat sheet (do not blur)

| Class | What it can say | What it cannot say |
|---|---|---|
| GSE137244 cell-line bulk | KL > KP Tacstd2/Cldn4 at library n=5 vs 5 | Tumor immune coldness; human antigen level |
| GSE179500 mouse bulk | LKB1 restore lowers Tacstd2 (sample-level) | Immune exclusion; human STK11 |
| Mouse scRNA (179502/180963) | Epithelial expression / restore direction (small n) | Powered genotype→immune or genotype→antigen without caveats |
| TCGA human bulk | STK11 → **lower** TACSTD2/CLDN4; weak TACSTD2–CD8A ρ | STK11 = Tacstd2/Cldn4-high disease |
| GSE280232 human scRNA | TACSTD2 is epithelial in KRAS tumors | STK11mut > wt TACSTD2 |

---

## Verdict line for slide footer

> Disease context = STK11/KL cold (Skoulidis + public MSK). Tacstd2/Cldn4-high locks on **mouse KL cell-line bulk GSE137244** (+3.24 / +5.57, p=0.00794), with mouse Lkb1-loss support. **Human bulk TCGA is reverse for antigen.** Keep bulk ≠ scRNA.

---

## Do not say

- Do not write “human STK11 tumors are Tacstd2/Cldn4-high” from TCGA — they are lower.
- Do not use GSE180963 n=1 as KL>K Tacstd2 or as immune-fraction cold.
- Do not merge private 8KL scRNA with public mouse.
- Do not call TISMO LLC a KL/KP lung model.
- Do not upgrade GSE280232 mixed epi % into a genotype effect.
- Do not treat GSE137244 IFN soft trends as the locked cold claim (cold = STK11 clinical/TME literature; antigen = GSE137244).
