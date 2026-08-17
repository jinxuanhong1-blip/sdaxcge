# FINDING — GSE126044 CLDN4 Q4 vs Q1 prerank GSEA (IFN / MHC / TJ)

**Additive only.** Public Cho 2020 anti-PD-1 NSCLC bulk RNA-seq. This folder does **not** audit or retract any slide. Single-gene CLDN4 vs ESTIMATE ImmuneScore is **already known** (Spearman r = **−0.524**, p = **0.037**, n = **16**; `methods/bulk_immune`) and is recorded, not re-tested. The extra is prerank GSEA on **CLDN4 Q4 vs Q1**.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `scrna_pseudobulk_gsea_meta` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

**Honest n.** Cohort n = **16**. After the quartile cut, CLDN4 Q4 n = **4**, Q1 n = **4**. That is the contrast. Gene-set permutation on a prerank does **not** become a 16-sample test.

---

## 一句话 / TL;DR

| Contrast | Honest n | Known CLDN4 vs ImmuneScore | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE126044 CLDN4 Q4 vs Q1 | 4 vs 4 (cohort 16) | Spearman r=−0.524, p=0.037, n=16 (given; not re-tested) | Hallmark IFN-γ NES -2.670 FDR 0.001; MHC-I / APM NES -2.216 FDR 0.001; KEGG tight junction NES +1.023 FDR 0.326 |

Primary rank is Welch *t* on log2(CPM+1). Q4 is **4/4 non-responder** and includes **2 FFPE**; Q1 is **2/4 responder**, all fresh. Bulk biopsy, so IFN / MHC NES can be infiltrate.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public GEO counts only. No FASTQ / SRA. |
| File | `GSE126044_counts.txt.gz` |
| Cohort | 16 pre-treatment NSCLC biopsies, anti-PD-1 (Cho 2020, PMID 32879421) |
| Labels | 5 responder / 11 non-responder; 11 fresh / 5 FFPE (all FFPE are NR) |
| Target | CLDN4 log2(CPM+1) |
| Contrast | Q4 vs Q1 via `gsea_core.quartile_split` (≤P25 vs ≥P75) |
| Honest n | Q4 = 4 (Dis_01, Dis_12, Dis_06, Dis_08); Q1 = 4 (Dis_04, Dis_03, Dis_07, Dis_10) |
| Sign | t>0 / log2FC>0 = **up in CLDN4 Q4** |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | mean log2(CPM+1) difference Q4−Q1; Spearman ρ vs continuous CLDN4 (n=16) |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets |
| Out of scope | ImmuneScore re-test; TACSTD2 median A8 GSEA rerun; FASTQ |

Series: [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044). Title on GEO: *Genome-wide identification of differentially methylated promoters and enhancers associated with response to anti-PD-1 therapy in non-small cell lung cancer*.

---

## Known single-gene vs ImmuneScore (not re-tested)

From `methods/bulk_immune/results/demo/GSE126044_correlation_spearman.tsv` (same 16 samples, same CLDN4):

| Target | Score | n | Spearman r | p | FDR (within family) |
|---|---|---:|---:|---:|---:|
| CLDN4 | `estimate:ImmuneScore` | 16 | −0.524 | 0.037 | 0.142 |

This extra does not recompute ESTIMATE / MCP-counter / xCell.

---

## Headline NES (Welch *t* rank, Q4 minus Q1)

18558 genes ranked. CLDN4 itself: median log2(CPM+1) Q4 = **4.429**, Q1 = **1.107** (the split gene; not a GSEA claim).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -2.670 | 0.001 | 0.001 | 198 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -2.216 | 0.001 | 0.001 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.023 | 0.326 | 0.326 | 164 |

Hallmark IFN-α (secondary, not in BH): NES -2.101 nom p 0.001 (n=95; not in headline FDR).

Leading-edge (first 25, Welch *t* rank):

- IFN-γ: `CASP8,SECTM1,OASL,IL15RA,HLA-DQA1,ST3GAL5,MYD88,PELI1,SELP,TAPBP,C1R,PIM1,STAT1,PSMB8,BPGM,LAP3,BATF2,HLA-DRB1,IFIH1,PTPN1,HLA-A,ARID5B,BST2,FGL2,HLA-DMA`
- MHC-I / APM: `TAP1,CANX,CALR,TAPBP,PSMB8,HLA-A,HLA-F,PSMB10,HLA-E,HLA-C,HLA-B,B2M,IRF1,PSMB9,HLA-G`
- KEGG TJ: `CLDN4,CLDN20,LLGL2,SCRIB,RAC1,ERBB2,AMOT,PRKAA2,PARD6B,TJAP1,MAP2K7,LLGL1,TJP2,MARVELD2,ARPC1A,CGN,ACTN4,OCLN,AMOTL2,CLDN7,PRKCZ,PPP2R1B,NEDD4L,SRC,EPB41L4B`

---

## Sensitivity: log2FC rank (Q4 − Q1 mean)

Same three sets, rank = mean log2(CPM+1) difference instead of Welch *t*. Still n = 4 vs 4.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -2.949 | 0.001 | 0.001 | 198 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -2.170 | 0.001 | 0.001 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.418 | 0.005 | 0.005 | 165 |

---

## Sensitivity: Spearman vs continuous CLDN4 (uses all 16)

Rank = per-gene Spearman ρ vs CLDN4. This uses the full cohort (n=16) and is **not** a Q4 vs Q1 test. Reported so the extreme-quartile NES is not the only cut.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -2.605 | 0.001 | 0.001 | 198 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -2.123 | 0.001 | 0.001 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.130 | 0.078 | 0.078 | 164 |

---

## Extra figures

- `methods/gse126044_cldn4_gsea/figures/fig_headline_nes.png`
- `methods/gse126044_cldn4_gsea/figures/fig_qc_and_nes.png`

Tables: `methods/gse126044_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `sample_table.tsv`, `inventory.tsv`.

---

## What this is not

- Not a re-test of CLDN4 vs ImmuneScore / CD8 / MCP NK (already in `methods/bulk_immune`).
- Not the A8 TACSTD2 **median** GSEA on this same series (that was 8 vs 8, TJ/keratin/EMT).
- Not FASTQ / salmon / DESeq2 from SRA.
- Not sample-permutation GSEA. Honest contrast n = **4 vs 4**; the engine is gene-set permutation on a prerank.
- Not a response GSEA. Q4 happens to be 4/4 NR and includes 2 FFPE; that confounding is recorded, not adjusted.

---

## 中文摘要

只补公开 **GSE126044**（Cho 2020，抗 PD-1 NSCLC，n=16）上 **CLDN4 Q4 vs Q1** 的 prerank GSEA，不审不撤已有页。CLDN4 与 ESTIMATE ImmuneScore 的单基因相关已经知道（Spearman r=−0.524，p=0.037，n=16），这里不重测。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- 诚实 n：队列 16；四分位切割后 Q4 = 4，Q1 = 4。这不是 16 对 16 的样本置换检验。
- Q4 为 4/4 非应答，含 2 例 FFPE；Q1 为 2/4 应答、全部新鲜组织。
- Headline（Welch *t* 秩）：Hallmark IFN-γ NES -2.670 FDR 0.001; MHC-I / APM NES -2.216 FDR 0.001; KEGG tight junction NES +1.023 FDR 0.326。
