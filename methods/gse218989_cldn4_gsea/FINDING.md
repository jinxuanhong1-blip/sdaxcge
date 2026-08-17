# FINDING — GSE218989 CLDN4 Q4 vs Q1 prerank GSEA (IFN-γ / MHC-I / TJ)

**Additive only.** Public SMC–KAIST PD-1/PD-L1 NSCLC bulk TPM (Kang et al., *Nat Commun* 2024, PMID 38744958; GEO [GSE218989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218989)). This folder does **not** audit or retract `methods/gse218989_cldn4_ici`. No slide was re-scored. Unit is the **patient**.

Prerank engine is the same as `methods/cldn4_ko_gsea` / `gse289287_trop2ko_gsea` (`gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4 Q4** end of the rank (Q4 minus Q1). BH-FDR is within the **three headline sets**. Numbers below are written from `tables/gsea_headline.tsv`.

**Do not headline response OR.** CLDN4 vs GEO ICI response is already known NS in `methods/gse218989_cldn4_ici` (MWU p=0.40, AUC=0.47; CD8A / Ayers IFN-γ do separate responders). The Q4 vs Q1 2×2 is recorded in `tables/response_or_note.tsv` only.

---

## 一句话 / TL;DR

| Contrast | n (Q4 vs Q1) | Rank | IFN-γ / MHC-I / KEGG TJ NES (FDR) |
|---|---|---|---|
| GSE218989 CLDN4 Q4 vs Q1 | 89 vs 89 | Welch *t* on log2(TPM+1) | Hallmark IFN-γ NES +0.905 FDR 0.886; MHC-I / APM NES +0.952 FDR 0.884; KEGG tight junction NES +1.191 FDR 0.039 |

Honest read: Hallmark IFN-γ is NS (NES +0.905, FDR 0.886); MHC-I / APM is NS (NES +0.952, FDR 0.884); KEGG tight junction is up in Q4 (NES +1.191, FDR 0.039). Histology is **not deposited**. This is bulk ICI RNA — an IFN / MHC NES can be infiltrate, not a tumour-cell-intrinsic call.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrix | Public author TPM only. No FASTQ / SRA. |
| File | `GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz` |
| n | 355 patients × 19916 protein-coding genes |
| Split | CLDN4 `log2(TPM+1)` quartiles: Q1 ≤ 1.1631, Q4 ≥ 2.9718 |
| Contrast | CLDN4 Q4 (n=89) vs Q1 (n=89) |
| Sign | Welch *t* > 0 = **higher in CLDN4 Q4** |
| Engine | Same prerank as `gsea_core.py` in `cldn4_ko_gsea` |
| Rank (primary) | Welch *t*, Q4 minus Q1 |
| Rank (sensitivity) | Spearman ρ vs continuous CLDN4 (all 355); Welch *t* after dropping CLDN4 from the rank |
| Headline sets | Hallmark IFN-γ, custom MHC-I/APM (21 genes), KEGG tight junction |
| Secondary | Hallmark IFN-α; KEGG TJ without CLDN4 (reported, not in BH) |
| FDR | BH inside the 3 headline sets |
| Out of scope | Headline response OR; LUAD vs LUSC (not deposited); ADC+ICI |

CLDN4 itself is a KEGG TJ member and is the split gene. Its Welch *t* on this contrast is **+26.11**. TJ NES on the full rank is therefore partly circular. The drop-CLDN4 / TJ-minus-CLDN4 rows are the non-circular companions.

---

## Headline NES (Welch *t*, Q4 minus Q1)

19559 genes ranked. 89 Q4 vs 89 Q1.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +0.905 | 0.886 | 0.886 | 200 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +0.952 | 0.884 | 0.589 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.191 | 0.039 | 0.013 | 168 |

Hallmark IFN-α (secondary, not in BH): NES +1.097 nom p 0.200 (n=97; not in headline FDR).

Leading-edge (first 25, Welch rank):

- IFN-γ: `VAMP8,LGALS3BP,LY6E,IFNAR2,RBCK1,TRIM25,PSMA2,ISG15,UPP1,IRF7,DHX58,PSMB2,MVP,ARL4A,METTL7B,CASP3,PSMB10,CASP7,TRIM14,IRF9,SECTM1,TOR1B,IFI35,ISOC1,PSMA3`
- MHC-I / APM: `PSMB10,TAPBPL,CANX,B2M,PSMB8,PSMB9,PDIA3,CALR,IRF1,HLA-G,TAP2,TAPBP`
- KEGG TJ: `CLDN4,CRB3,CLDN7,CLDN3,F11R,CLDN23,LLGL2,ARPC1A,SRC,MARVELD3,DLG3,CGN,CLDN9,RAC1,TJP3,ERBB2,ARPC5L,ARPC1B,ARPC4,RAP2C,PRKAB1,ARPC5,MAPK9,EPB41L4B,ARHGEF18`

---

## Sensitivity

Same three headline sets.

**Spearman ρ vs continuous CLDN4** (all 355 patients; not a quartile cut):

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +0.863 | 0.950 | 0.950 | 200 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +0.876 | 0.950 | 0.724 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.225 | 0.018 | 0.006 | 168 |

**Welch *t* after dropping CLDN4 from the rank** (same Q4/Q1 patients):

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +0.905 | 0.869 | 0.869 | 200 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +0.960 | 0.869 | 0.590 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | +1.173 | 0.069 | 0.023 | 167 |

KEGG TJ minus CLDN4 on the primary rank (not in BH): NES +1.168 nom p 0.026 (n=167; not in headline FDR).

KEGG TJ on the drop-CLDN4 rank (not in BH, same as the TJ row above): NES +1.173 nom p 0.023 (n=167; not in headline FDR).

---

## Response OR — already known NS; not a headline

From `methods/gse218989_cldn4_ici`: CLDN4 vs GEO response MWU p=0.40, AUC=0.47 (168 R / 187 NR). CD8A p=0.0027 and Ayers IFN-γ p=0.00070 do separate responders, so the endpoint is not dead.

Q4 vs Q1 2×2 on the same GEO labels (note only): Q4 37 R / 52 NR vs Q1 47 R / 42 NR; OR 0.64 (95% CI 0.35–1.15), Fisher p=0.176. That is the already-known null, not a new claim. Table: `tables/response_or_note.tsv`.

---

## Extra figures

- `methods/gse218989_cldn4_gsea/figures/fig_headline_nes.png`
- `methods/gse218989_cldn4_gsea/figures/fig_nes_primary_vs_spearman.png`

Tables: `methods/gse218989_cldn4_gsea/tables/gsea_headline.tsv` (the NES table), `gsea_prerank_all.tsv`, `inventory.tsv`, `quartile_n.tsv`, `response_or_note.tsv`.

---

## What this is not

- Not a re-run or retraction of `methods/gse218989_cldn4_ici`.
- Not a response-OR headline. That test is already NS.
- Not FASTQ / salmon / DESeq2 from SRA (author TPM is used as deposited).
- Not a LUAD vs LUSC split — histology is not on GEO or Supplementary Data 8.
- Not an ADC+ICI trial. ICI-only cohort.
- Not sample-permutation GSEA. The engine is gene-set permutation on a prerank.
- Not a tumour-cell-intrinsic IFN/MHC call. Bulk ICI RNA.

---

## 中文摘要

只补公开 **GSE218989** CLDN4 四分位 Q4 vs Q1 的 prerank GSEA（IFN-γ / MHC-I / KEGG TJ），不审不撤已有 `gse218989_cldn4_ici` 页。**不要把缓解 OR 当标题**——那边已经是 NS（GEO MWU p=0.40）。引擎与 `cldn4_ko_gsea` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 Q4 端富集。FDR 只在三个 headline set 内做 BH。

- 89 Q4 vs 89 Q1；秩 = Welch *t*（log2(TPM+1)）。
- Headline：Hallmark IFN-γ NES +0.905 FDR 0.886; MHC-I / APM NES +0.952 FDR 0.884; KEGG tight junction NES +1.191 FDR 0.039。
- CLDN4 在 KEGG TJ 里，又是分组基因，TJ 全秩 NES 有循环成分；去 CLDN4 的伴随结果见上表。
