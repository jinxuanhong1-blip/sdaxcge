# FINDING — GSE31210 CLDN4 Q4 vs Q1 prerank GSEA (IFN / MHC / TJ)

**Additive extra. CLDN4 only.** Public Okayama / Kohno Japanese stage I–II LUAD (Okayama et al., *Cancer Res* 2012, PMID 22261853; GEO [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210); GPL570). This folder does **not** audit or retract any slide. TACSTD2 is **not** an anchor here.

**CD8 ρ is already known** in `methods/gse31210_cldn4_immune` (CLDN4 vs CD8A ρ=−0.341, n=226) and is **not** re-tested or re-headlined. This folder only adds prerank GSEA of IFN-γ / MHC-I / KEGG TJ on the same tumor set.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `methods/gse285029_cldn4_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

No ICI labels. Bulk MAS5 microarray, not a cell-intrinsic call.

---

## 一句话 / TL;DR

| Contrast | n | Split | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE31210 CLDN4 Q4 vs Q1 | 226 (Q4=57, Q1=57) | CLDN4 quartiles on log2(MAS5+1) | Hallmark IFN-γ NES -1.614 FDR 0.001; MHC-I / APM NES -1.184 FDR 0.125; KEGG tight junction NES +1.807 FDR 0.001 |

Primary rank is Welch *t* (Q4 − Q1) on unique-mapped max-mean collapse.

**Honest read:** Hallmark IFN-γ is **down** in CLDN4 Q4 (NES -1.614 FDR 0.001). That is the same direction as the already-known CLDN4–CD8A anti-correlation, but it is a program NES, not a CD8 re-test. MHC-I / APM is the same sign and **NS** (FDR 0.125; Spearman sensitivity FDR 0.083) — reported as computed, not rounded into a call. KEGG TJ is up (NES +1.807 FDR 0.001). CLDN4 is a member of that set, so the primary TJ NES includes the stratification gene. After dropping CLDN4 from the rank: NES +1.700 FDR 0.001 nom p 0.001 (n=157).

---

## Honest n

| item | n | rule |
|---|---:|---|
| arrays on the series matrix | 246 | GEO `GSE31210_series_matrix.txt.gz` |
| primary lung tumors | **226** | `tissue: primary lung tumor` |
| adjacent / non-tumor | 20 | dropped |
| CLDN4 Q4 / Q1 | 57 / 57 | ≤P25 vs ≥P75; ties stay in the tail |
| genes after unique-mapped max-mean | 20848 | multi-mapped `///` probes dropped |
| genes ranked (Welch) | 20848 | finite Welch *t*, Q4 minus Q1 |

Same tumor rule as `methods/gse31210_cldn4_immune`. Not a new cohort hunt.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public GEO series matrix only. No CEL / FASTQ / SRA. |
| File | `GSE31210_series_matrix.txt.gz` + `GPL570.annot.gz` |
| Samples | **n = 226** primary lung tumors (of 246 arrays) |
| Platform | GPL570 Affymetrix U133 Plus 2.0 |
| Collapse | unique-mapped max-mean (same as sibling immune folder) |
| CLDN4 probe | `201428_at` |
| Transform | `log2(MAS5+1)` |
| Anchor | **CLDN4 only** |
| Split | Q4 vs Q1 on CLDN4 (≤P25 vs ≥P75). Ties stay in the tail. |
| n Q4 / Q1 | 57 / 57 |
| CLDN4 log2 cuts | P25 = 9.1062; P75 = 10.4280 |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (all 226 tumors) |
| TJ honesty | Drop CLDN4 from the rank and re-run the three headline sets |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets, **per rank** |
| Out of scope | CD8 ρ (already known); TACSTD2 ranks; ICI labels; sample-permutation GSEA |

Series: [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210). Okayama et al., *Cancer Res* 2012; PMID [22261853](https://pubmed.ncbi.nlm.nih.gov/22261853/). Treatment-naive surgical LUAD. No ICI arm.

---

## Headline NES (Welch *t*, Q4 − Q1)

20848 genes ranked. CLDN4 itself is the stratification gene (Welch *t* = +27.535; by construction the top of the rank).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -1.614 | 0.001 | 0.001 | 195 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -1.184 | 0.125 | 0.125 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.807 | 0.001 | 0.001 | 158 |

Hallmark IFN-α (secondary, not in BH): NES -0.936 nom p 0.345 (n=92; not in headline FDR).

Leading-edge (first 25, Welch rank):

- IFN-γ: `IFIT2,SP110,CD40,IRF2,IRF1,LYSMD2,VCAM1,SOCS3,HIF1A,FAS,HLA-DQA1,SOCS1,CFH,TXNIP,CD38,PTPN2,C1R,LAP3,ARID5B,CASP1,APOL6,ST8SIA4,RAPGEF6,SSPN,B2M`
- MHC-I / APM: `PSMB9,TAP1,TAPBPL,HLA-B,HLA-A,ERAP2,ERAP1,IRF1,HLA-C,B2M,TAP2,HLA-E,NLRC5`
- KEGG TJ: `CLDN4,CLDN3,CLDN7,CGN,ARPC1A,EPB41L4B,DLG3,CRB3,F11R,MARVELD3,OCLN,LLGL2,ERBB2,MARVELD2,TJP3,NEDD4L,ARPC1B,MYH14,PRKCI,ARPC5L,TUBA4A,EZR,AMOT,ARPC3,MYL12A`

---

## Sensitivity: drop CLDN4 from the rank

Same Q4 vs Q1 samples and the same three headline sets. CLDN4 is removed from the ranked list so KEGG TJ is not scored on the gene used to define the split.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -1.623 | 0.001 | 0.001 | 195 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -1.185 | 0.136 | 0.136 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.700 | 0.001 | 0.001 | 157 |

---

## Sensitivity: Spearman ρ vs continuous CLDN4

Uses all **n = 226** tumors. Rank = per-gene Spearman ρ vs CLDN4. BH-FDR is again inside the three headline sets.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -1.637 | 0.003 | 0.001 | 195 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -1.300 | 0.083 | 0.083 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.699 | 0.003 | 0.002 | 158 |

---

## Extra figures

- `methods/gse31210_cldn4_gsea/figures/fig_headline_nes.png`

Tables: `methods/gse31210_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `inventory.tsv`, `quartile_split.tsv`, `probe_used.tsv`.

---

## What this is not

- Not a re-run of the sibling CD8 / ImmuneScore / TACSTD2 ρ tables (`methods/gse31210_cldn4_immune`). Those stay as written. CD8 ρ is already known.
- Not a TACSTD2 analysis. CLDN4 only.
- Not CEL / RMA / FASTQ from SRA (series-matrix MAS5 is used as deposited).
- Not sample-permutation GSEA. The engine is gene-set permutation on a prerank.
- Not a cell-intrinsic IFN call. This is bulk LUAD microarray.
- Not an ICI-response analysis. Labels are not deposited.

---

## 中文摘要

只补公开 **GSE31210**（诚实 n=226 原发肺腺癌）上 **CLDN4 Q4 vs Q1** 的 prerank GSEA，不审不撤已有页。CD8 ρ 已在 `gse31210_cldn4_immune` 给出，这里不再当主结果。只做 CLDN4，不做 TACSTD2。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- 246 张芯片，**226** 例 `primary lung tumor`，20 例非肿瘤丢掉。Q4=57，Q1=57。
- Headline（Welch *t*）：Hallmark IFN-γ NES -1.614 FDR 0.001; MHC-I / APM NES -1.184 FDR 0.125; KEGG tight junction NES +1.807 FDR 0.001。
- CLDN4 是 KEGG tight junction 成员；去掉 CLDN4 后再跑：NES +1.700 FDR 0.001 nom p 0.001 (n=157)。
