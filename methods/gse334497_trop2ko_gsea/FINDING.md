# FINDING — GSE334497 4T1 Trop2 KO prerank GSEA

**Additive only.** This folder does **not** audit or retract any slide. It is **not** lung and it is **not** SKB264. It is **not** a CLDN4 knockout: the perturbation is CRISPR **Trop2 (Tacstd2)** in 4T1 mammary tumors.

Prerank engine is the same as `methods/scrna_pseudobulk_gsea_meta` (`scripts/gse334497_trop2ko_gsea/gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **Trop2-KO** end of the rank (Welch *t*, KO minus WT). BH-FDR is within the **four headline sets**. Numbers below are written from `tables/gsea_headline.tsv` and `tables/cldn4_log2fc.tsv`.

**Cldn4 gene-level *p* = 0.32 is given and is not re-tested.** Only the log2FC effect size is reported.

---

## 一句话 / TL;DR

| Contrast | n | Tacstd2 log2FC (QC) | Cldn4 log2FC | Cldn4 *p* | IFN-γ | MHC-I/APM | KEGG TJ | keratin |
|---|---|---|---|---|---|---|---|---|
| GSE334497 4T1 Trop2 KO vs WT | **5 vs 5** | -3.821 | -0.817 | **0.32 (given)** | UP NES +2.259 (FDR<0.05) | UP NES +1.688 (FDR<0.05) | DOWN NES -2.025 (FDR<0.05) | DOWN NES -2.539 (FDR<0.05) |

RNA is from **frozen whole-tumor sections** in immunocompetent BALB/c hosts. An IFN / MHC NES on this matrix can be stroma or infiltrate; it is not a tumor-cell-intrinsic call.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public processed file only: `GSE334497_normalized_counts.csv.gz` (GRCm38 Ensembl). No FASTQ / SRA. |
| Model | Mouse **4T1** mammary tumors, CRISPR Trop2 KO vs WT, 3 weeks in BALB/c (Wu *et al.*, *JITC* 2026). Breast, not lung. |
| Libraries | KO: KO162, KO164, KO165, KO172, RESUB-KO163R. WT: control170, RESUB-171R, RESUB-170R, RESUB-169R, RESUB-168R. |
| Transform | log2(normalized count + 1); Ensembl → official symbol, collapse by max mean. |
| Rank | Signed Welch *t*, KO minus WT. Positive = up in Trop2 KO. |
| Engine | Same prerank as `gsea_core.py` (weighted KS *p*=1, 1000 gene-set perms, seed=42, min size 8) |
| Headline sets | Hallmark IFN-γ (MSigDB `mh.all.v2023.2.Mm`); custom MHC-I/APM (human → mouse, HLA→H2); KEGG tight junction; GO keratinization |
| FDR | BH inside the **4** headline sets |
| Cldn4 gene test | **Not re-run.** Given *p* = 0.32. This folder reports log2FC only. |

Primary genes: mouse Cldn4 `ENSMUSG00000047501`; Tacstd2 `ENSMUSG00000051397` (perturbation QC only; not a GSEA-set member).

---

## Headline NES

GEO: [GSE334497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334497). 15073 genes ranked.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +2.259 | 0.001 | 0.001 | 184 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +1.688 | 0.010 | 0.010 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | -2.025 | 0.001 | 0.001 | 135 |
| GO keratinization | `GOBP_KERATINIZATION` | -2.539 | 0.001 | 0.001 | 38 |

| Gene | log2FC (KO − WT) | Note |
|---|---:|---|
| Tacstd2 | -3.821 | Perturbation QC. Knockout worked. |
| Cldn4 | -0.817 | Effect size only. Gene-level *p* = **0.32 (given; not re-tested)**. |

Leading-edge (first 25, from the walk):

| Set | Lead genes |
|---|---|
| Hallmark IFN-γ | `Trim21,H2-Aa,Ripk1,Il18bp,Cxcl9,Stat4,Lcp2,Icam1,Casp1,Nod1,B2m,Cd86,Ifnar2,Tnfaip2,Cd69,Herc6,Nampt,Fgl2,St8sia4,Fcgr1,Fpr1,Psme1,Il10ra,Cxcl10,Cd274` |
| MHC-I / APM | `B2m,H2-Q10,Irf1,Erap1,Nlrc5,Tap1,Tapbpl,Psmb9,H2-K1,Tap2,H2-T23,Psmb8,Psmb10` |
| KEGG tight junction | `Jam2,Erbb2,Llgl1,Rhoa,Ybx3,Tjap1,Actb,Arpc1a,Myl12b,Arhgef18,Myl6b,Amotl1,Actg1,Msn,Myl2,Hspa4,Nedd4l,Ppp2cb,Arpc5l,Cldn4,Cldn15,Tuba1c,Amot,Prkaa2,Tuba4a` |
| GO keratinization | `Krt6b,Krt6a,Cdh3,Tgm3,Loricrin,Evpl,Ivl,Krt16,Casp14,Cdsn,Abca12,Hrnr,Lce1d,Lce1c,Lce1b,Krt1,Krt77,Tmem79,Kazn,Sprr1a,Sfn,Krt5,Lipm,Krt79,Krt7` |

---

## What this is and is not

**Is**

- An additive prerank GSEA on the public 4T1 Trop2 KO 5-vs-5 tumor matrix.
- A report of NES / FDR for IFN-γ, MHC-I/APM, KEGG TJ, and keratin, plus Cldn4 log2FC.

**Is not**

- A CLDN4 knockdown or knockout.
- A lung model or SKB264.
- A re-audit of the given Cldn4 gene-level *p* = 0.32.
- Tumor-cell-intrinsic IFN/APM (bulk immunocompetent tumor RNA).
- Sample-permutation GSEA. The null is gene-set permutation on a prerank.

---

## Figures

- `methods/gse334497_trop2ko_gsea/figures/fig_headline_nes.png`
- `methods/gse334497_trop2ko_gsea/figures/fig_nes_bars.png`

Tables: `methods/gse334497_trop2ko_gsea/tables/gsea_headline.tsv`, `gsea_prerank_all.tsv`, `cldn4_log2fc.tsv`, `inventory.tsv`.

---

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/gse334497_trop2ko_gsea/analyze.py
```

---

## 中文摘要

只补 GSE334497 **4T1 Trop2 KO 肿瘤 5 vs 5** 的 prerank GSEA，不审不撤已有页。不是肺，不是 SKB264，也不是 CLDN4 KO。引擎与 `scrna_pseudobulk_gsea_meta` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = Trop2 KO 端富集。Cldn4 单基因 *p*=0.32 **按给定，不重测**，只报 log2FC。

- Tacstd2 log2FC -3.821（敲除成立）。
- Cldn4 log2FC -0.817（*p*=0.32 给定）。
- IFN-γ：UP NES +2.259 (FDR<0.05)。
- MHC-I/APM：UP NES +1.688 (FDR<0.05)。
- KEGG 紧密连接：DOWN NES -2.025 (FDR<0.05)。
- GO 角化：DOWN NES -2.539 (FDR<0.05)。
- 全瘤冰冻切片、免疫健全宿主，IFN/MHC NES 可以来自间质/浸润，不能写成肿瘤内在。
