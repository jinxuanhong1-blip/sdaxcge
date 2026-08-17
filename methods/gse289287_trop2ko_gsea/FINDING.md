# FINDING — GSE289287 T-47D Trop-2 KO xenograft prerank GSEA

**Additive only.** Public author DESeq2 on T-47D **Trop-2 KO vs WT xenografts**. This folder does **not** audit or retract any slide. It is **not** lung and it is **not** SKB264.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `scrna_pseudobulk_gsea_meta` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **Trop-2 KO** end of the rank (KO minus WT). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

CLDN4 single-gene padj **0.466** is the known result (log2FC +0.284); it is recorded, not treated as a new claim.

---

## 一句话 / TL;DR

| Contrast | n (KO vs WT) | TACSTD2 log2FC (padj) | CLDN4 log2FC (padj) | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|---|
| GSE289287 T-47D Trop-2 KO xenograft | 4 vs 3 | -3.263 (1.95e-61) | +0.284 (0.466) | Hallmark IFN-γ NES +1.873 FDR 0.003; MHC-I / APM NES +1.081 FDR 0.187; KEGG tight junction NES -1.093 FDR 0.165 |

Primary rank is the author DESeq2 **Wald statistic**. Hosts are **NRG** females (mammary fat pad), so any IFN / APM movement is **tumour-cell-intrinsic**, not adaptive infiltrate.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public author DESeq2 only. No FASTQ / SRA. |
| File | `GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz` |
| Contrast | Trop-2 KO xenografts (GSM8788422/23/24/25; animals 2815/2817/2818/2807) vs WT (GSM8788419/20/21; animals 2812/2808/2810) |
| Sign | log2FC>0 and Wald stat>0 = **up in Trop-2 KO** |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | DESeq2 Wald `stat` |
| Rank (sensitivity) | DESeq2 `log2FoldChange` |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α (reported, not in BH) |
| FDR | BH inside the 3 headline sets |
| Out of scope | DSG2 KO tables; in-vitro T-47D (no Trop-2 KO cell DESeq2 on GEO); lung; SKB264 |

Series: [GSE289287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289287). Title on GEO: *Trop-2 governs anti-metastatic desmosomal integrity* (Vacek / Souček Lab). Platform GPL24676, human T-47D in NRG mammary fat pad.

---

## Headline NES (Wald-stat rank)

26253 genes ranked. TACSTD2 itself: log2FC **-3.263**, Wald **-17.095**, padj **1.95e-61** (KO worked). CLDN4: log2FC **+0.284**, padj **0.466** (known, ns).

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +1.873 | 0.003 | 0.001 | 186 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +1.081 | 0.187 | 0.187 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | -1.093 | 0.165 | 0.110 | 153 |

Hallmark IFN-α (secondary, not in BH): NES +2.211 nom p 0.001 (n=95; not in headline FDR).

Leading-edge (first 25, Wald rank):

- IFN-γ: `IFI44L,ISG15,IFI44,IFITM3,STAT2,MYD88,UBE2L6,DDX58,IFIT3,OAS3,OAS2,MT2A,MX1,MX2,CMPK2,STAT1,BST2,IFIT1,IFI35,IL4R,ISG20,DDX60,TDRD7,PARP12,XAF1`
- MHC-I / APM: `HLA-F,HLA-C,B2M,HLA-B,HLA-A,NLRC5,TAP1,HLA-E,PSMB8,TAPBPL`
- KEGG TJ: `SRC,ARPC1A,IGSF5,PRKAA1,CGNL1,ROCK2,ARPC3,HSPA4,MYH11,PPP2CB,PRKAA2,ERBB2,RAC1,RAP2C,PRKAG2,NF2,ACTR3B,PPP2R2B,ROCK1,ACTR2,CD1D,NEDD4,RHOA,PARD3,PRKACB`

---

## Sensitivity: log2FC rank

Same three sets, rank = author `log2FoldChange` instead of Wald stat.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +1.641 | 0.003 | 0.001 | 186 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +1.121 | 0.250 | 0.167 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | -0.847 | 0.365 | 0.365 | 153 |

---

## Extra figures

- `methods/gse289287_trop2ko_gsea/figures/fig_headline_nes.png`
- `methods/gse289287_trop2ko_gsea/figures/fig_qc_and_nes.png`

Tables: `methods/gse289287_trop2ko_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `key_genes.tsv`, `inventory.tsv`.

---

## What this is not

- Not a re-run of any existing slide or of SKB264.
- Not FASTQ / salmon / DESeq2 from SRA (author table is used as deposited).
- Not lung. T-47D is luminal breast; hosts are NRG mice.
- Not the DSG2-KO arms of the same series.
- Not sample-permutation GSEA. n = 4 vs 3; the engine is gene-set permutation on a prerank.

---

## 中文摘要

只补公开 **GSE289287** T-47D Trop-2 KO 移植瘤的 author DESeq2 prerank GSEA，不审不撤已有页。不是肺癌，不是 SKB264。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = Trop-2 KO 端富集。FDR 只在三个 headline set 内做 BH。

- 4 KO vs 3 WT（NRG 乳腺脂肪垫）。TACSTD2 log2FC -3.263，padj 1.95e-61（敲除成立）。
- CLDN4 log2FC +0.284，padj 0.466（已知，单基因不显著）。
- Headline（Wald 秩）：Hallmark IFN-γ NES +1.873 FDR 0.003; MHC-I / APM NES +1.081 FDR 0.187; KEGG tight junction NES -1.093 FDR 0.165。
