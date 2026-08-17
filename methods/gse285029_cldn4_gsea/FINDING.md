# FINDING — GSE285029 CLDN4 Q4 vs Q1 prerank GSEA

**Additive only. CLDN4 only.** Public pre-ICI NSCLC WTS, **n = 234** (Koh et al., *JITC* 2025; GEO [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029)). This folder does **not** audit or retract any slide. TACSTD2 is **not** an anchor here.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `methods/gse289287_trop2ko_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

GEO has no RECIST, histology, PD-L1 IHC, TMB, or purity. This is bulk tumour WTS, not a cell-intrinsic call.

---

## 一句话 / TL;DR

| Contrast | n | Split | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE285029 CLDN4 Q4 vs Q1 | 234 (Q4=59, Q1=59) | CLDN4 quartiles on log2(clip0+1) | Hallmark IFN-γ NES +1.377 FDR 0.001; MHC-I / APM NES +1.368 FDR 0.049; KEGG tight junction NES +1.373 FDR 0.001 |

Primary rank is Welch *t* (Q4 − Q1) on the author matrix. MHC-I / APM is FDR **0.049** on this rank (nom p 0.049) and FDR **0.097** on the Spearman sensitivity — reported as computed, not rounded into a stronger call. CLDN4 is a member of KEGG tight junction, so the primary TJ NES includes the stratification gene. After dropping CLDN4 from the rank: NES +1.350 FDR 0.001 nom p 0.001 (n=164).

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public author WTS only. No FASTQ / SRA. |
| File | `GSE285029_WTS_expr_count_235_032820.txt.gz` |
| Samples | **n = 234** (Case1–Case234). Author ICI–RNA-seq cohort. |
| Genes | 20067 symbols on the deposited matrix |
| Transform | `log2(pmax(x,0)+1)` (same clip as the sibling GSE285029 score analysis) |
| Anchor | **CLDN4 only** |
| Split | Q4 vs Q1 on CLDN4 (≤P25 vs ≥P75). Ties stay in the tail. |
| n Q4 / Q1 | 59 / 59 |
| CLDN4 log2 cuts | P25 = 5.5013; P75 = 6.1950 |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (all 234 samples) |
| TJ honesty | Drop CLDN4 from the rank and re-run the three headline sets |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets, **per rank** |
| Out of scope | TACSTD2 ranks; response labels; FASTQ; sample-permutation GSEA |

Series: [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029). Title on GEO: pre-ICI NSCLC tumour WTS (PD-1 or PD-L1 blockade). Paper: Koh et al., *J Immunother Cancer* 2025; PMID [40050048](https://pubmed.ncbi.nlm.nih.gov/40050048/).

The filename says “count”. Values are continuous with negatives; they are the author-processed matrix, not raw integer counts and not TPM. Rank correlations / Welch *t* on the clipped log matrix are the claim.

---

## Headline NES (Welch *t*, Q4 − Q1)

20062 genes ranked. CLDN4 itself is the stratification gene (Welch *t* = +11.368; by construction the top of the rank).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +1.377 | 0.001 | 0.001 | 198 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +1.368 | 0.049 | 0.049 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.373 | 0.001 | 0.001 | 165 |

Hallmark IFN-α (secondary, not in BH): NES +1.516 nom p 0.001 (n=95; not in headline FDR).

Leading-edge (first 25, Welch rank):

- IFN-γ: `PARP14,LY6E,TRIM14,TDRD7,HERC6,TRIM26,HELZ2,CMPK2,TRAFD1,PLA2G4A,USP18,MX2,MVP,OAS3,DDX60,CASP7,IFIH1,ZNFX1,PELI1,CASP3,CASP8,IRF2,PARP12,UPP1,IDO1`
- MHC-I / APM: `TAP2,TAPBP,HLA-A,TAPBPL,ERAP1,HLA-F,HLA-E,IRF1,CANX,HLA-C,B2M,HLA-B,PDIA3,TAP1,NLRC5,PSMB9,CALR,PSMB8`
- KEGG TJ: `CLDN4,CGN,MARVELD2,LLGL2,TJP3,CRB3,CLDN3,DLG3,ERBB2,CLDN7,OCLN,EPB41L4B,PARD6B,TJP2,F11R,CLDN23,PRKCZ,PARD3,CTTN,NEDD4L,EZR,MYH14,MAPK9,CCND1,MAP3K5`

---

## Sensitivity: drop CLDN4 from the rank

Same Q4 vs Q1 samples and the same three headline sets. CLDN4 is removed from the ranked list so KEGG TJ is not scored on the gene used to define the split.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +1.391 | 0.001 | 0.001 | 198 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +1.343 | 0.049 | 0.049 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.350 | 0.001 | 0.001 | 164 |

---

## Sensitivity: Spearman ρ vs continuous CLDN4

Uses all **n = 234** samples. Rank = per-gene Spearman ρ vs CLDN4. BH-FDR is again inside the three headline sets.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +1.227 | 0.010 | 0.007 | 198 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +1.283 | 0.097 | 0.097 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.391 | 0.003 | 0.001 | 165 |

---

## Extra figures

- `methods/gse285029_cldn4_gsea/figures/fig_headline_nes.png`

Tables: `methods/gse285029_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `inventory.tsv`, `quartile_split.tsv`.

---

## What this is not

- Not a TACSTD2 analysis. CLDN4 only.
- Not a re-run of the sibling GSE285029 score / ρ tables (`results/w200/A11_GSE285029`). Those stay as written.
- Not FASTQ / salmon / DESeq2 from SRA (author matrix is used as deposited).
- Not sample-permutation GSEA. The engine is gene-set permutation on a prerank.
- Not a cell-intrinsic IFN call. This is bulk pre-ICI NSCLC WTS.
- Not a response / PD-L1 IHC / purity analysis. GEO does not release those labels.

---

## 中文摘要

只补公开 **GSE285029**（n=234）上 **CLDN4 Q4 vs Q1** 的 prerank GSEA，不审不撤已有页。只做 CLDN4，不做 TACSTD2。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- 234 例 ICI 前 NSCLC 肿瘤 WTS。Q4=59，Q1=59。
- Headline（Welch *t*）：Hallmark IFN-γ NES +1.377 FDR 0.001; MHC-I / APM NES +1.368 FDR 0.049; KEGG tight junction NES +1.373 FDR 0.001。
- CLDN4 是 KEGG tight junction 成员；去掉 CLDN4 后再跑：NES +1.350 FDR 0.001 nom p 0.001 (n=164)。
