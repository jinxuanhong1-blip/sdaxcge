# FINDING — public CLDN4 / Cldn4-loss prerank GSEA

**Additive only.** User thesis is taken as given: CLDN4 loss shares the IFN / MHC-I axis; CLDN4-high is the tight-junction barrier. This folder does **not** audit or retract any slide. It is **not** SKB264 and it is **not** a lung-cancer KO (GSE50927 is mouse **whole lung**).

Prerank engine is the same as `methods/scrna_pseudobulk_gsea_meta` (`scripts/cldn4_ko_gsea/gsea_core.py`): weighted KS *p*=1, **1000** gene-set permutations, seed=**42**. Positive NES = the set is enriched at the **CLDN4-loss** end of the rank (KO / siRNA minus WT / CLDN4-high). BH-FDR is within the **five headline sets** per contrast. Numbers below are written from `tables/gsea_headline.tsv` and `tables/cldn4_log2fc.tsv`.

---

## 一句话 / TL;DR

| Contrast | n (loss vs WT) | CLDN4 / Cldn4 log2FC | Do IFN / MHC-I go UP after CLDN4 loss? |
|---|---|---|---|
| GSE207704 T47D CLDN4-/- | 2 vs 2 (group-mean FPKM) | −1.039 | No on this rank: IFN-γ NES −1.567 FDR 0.020; IFN-α NES −1.829 FDR 0.005; MHC-I NES −1.232 FDR 0.209 (n_set=8) |
| GSE207704 MCF7 CLDN4-/- | 2 vs 2 (group-mean FPKM) | −0.750 | Mixed / NS: IFN-γ NES −1.110 FDR 0.185; IFN-α NES −1.008 FDR 0.251; MHC-I NES **+0.868** FDR 0.335 (n_set=8) |
| GSE207704 mean of both lines | descriptive mean rank | −0.894 | No on this rank: IFN-γ NES −1.542 FDR 0.017; IFN-α NES −1.637 FDR 0.010; MHC-I NES −0.819 FDR 0.336 |
| GSE50927 Cldn4 KO naive **whole lung** | **1 vs 1** | −6.061 | **Yes:** IFN-γ NES **+1.508** FDR 0.005; IFN-α NES **+1.587** FDR 0.010; MHC-I NES **+1.423** FDR 0.050 |
| GSE22493 SKOV-3 CLDN4 siRNA vs OE | 3 arrays | −1.225 | Not at FDR<0.05: IFN-γ NES −1.084 FDR 0.323; IFN-α NES −0.611 FDR 0.622; MHC-I NES −0.852 FDR 0.468 |

The only **in-host** CLDN4-loss matrix here is GSE50927 (mouse whole lung, mixed-cell bulk, n=1 vs 1). On that rank, Hallmark IFN-γ, Hallmark IFN-α, and MHC-I/APM are all **UP** after Cldn4 loss (FDR ≤ 0.05). GSE207704 and GSE22493 are in-vitro ranks (collapsed FPKM; OE-control array).

GSE22493 **was included**: GPL10555 ORF (else DESCRIPTION prefix) maps **27878 / 36284** probes to **18682** gene symbols. TACSTD2 is **not** on that platform. On GSE207704, TACSTD2 log2FC is −0.623 (T47D) and −1.008 (MCF7). On GSE50927, Tacstd2 author logFC is −0.145 (submitter FDR = 1). See `tables/tacstd2_log2fc.tsv`.

---

## Locked design (before NES)

| Item | Choice |
|---|---|
| Matrices | Public processed files only. No FASTQ / SRA. |
| Engine | Same prerank as `gsea_core.py` in `scrna_pseudobulk_gsea_meta` |
| Rank | log2FC of CLDN4-loss minus WT (or CLDN4-OE). Positive = up after loss. |
| Headline sets | Hallmark IFN-γ, Hallmark IFN-α, custom MHC-I/APM (21 genes), KEGG tight junction, GO keratinization |
| Mouse Hallmark IFN | MSigDB `mh.all.v2023.2.Mm` (not title-cased human lists) |
| Mouse MHC-I / TJ / keratin | human set → mouse symbols (HLA → classical H2; ERAP2 → Lnpep) |
| FDR | BH inside the 5 headline sets, per contrast |
| GSE207704 | cuffdiff **group-mean FPKM** (replicates already collapsed on GEO) |
| GSE50927 | author edgeR table, **naive whole lung**, not a tumour |
| GSE22493 | series-matrix VALUE; control arm is CLDN4-**overexpressing**, not parental |

Primary genes: CLDN4 `ENSG00000189143`, mouse Cldn4 `ENSMUSG00000047501`. TACSTD2 `ENSG00000184292` is recorded when present; it is not a GSEA-set member.

---

## 1. GSE207704 — CLDN4 KO T47D and MCF7 (breast RNA-seq)

GEO: [GSE207704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207704). File: `GSE207704_CLDN4_RNAseq.txt.gz` (FPKM). Design on the record: T47D WT GSM6310640–41, T47D CLDN4-/- GSM6310642–43, MCF7 WT GSM6310644–45, MCF7 CLDN4-/- GSM6310646–47 (2 vs 2 per line). The open processed table has **one FPKM column per group**, so replicate-level *t* statistics cannot be computed.

CLDN4 itself (highest-FPKM locus after symbol collapse; pseudocount 0.5):

| Line | FPKM WT | FPKM KO | log2FC KO vs WT | n |
|---|---:|---:|---:|---|
| T47D | 42.506 | 20.436 | -1.039 | 2 vs 2, collapsed |
| MCF7 | 85.855 | 50.849 | -0.750 | 2 vs 2, collapsed |

### T47D CLDN4-/- vs WT

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -1.567 | 0.020 | 0.008 | 115 |
| Hallmark IFN-α | `HALLMARK_INTERFERON_ALPHA_RESPONSE` | -1.829 | 0.005 | 0.001 | 60 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -1.232 | 0.209 | 0.144 | 8 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | -0.963 | 0.375 | 0.375 | 121 |
| GO keratinization | `GOBP_KERATINIZATION` | +1.046 | 0.209 | 0.167 | 12 |

IFN / MHC-I after CLDN4 loss on this rank: IFN-γ and IFN-α **DOWN** (FDR 0.020 and 0.005); MHC-I **DOWN** (FDR 0.209; only 8 of 21 set genes are named in the cuffdiff table).

### MCF7 CLDN4-/- vs WT

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -1.110 | 0.185 | 0.111 | 115 |
| Hallmark IFN-α | `HALLMARK_INTERFERON_ALPHA_RESPONSE` | -1.008 | 0.251 | 0.201 | 60 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +0.868 | 0.335 | 0.335 | 8 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | -1.239 | 0.132 | 0.053 | 121 |
| GO keratinization | `GOBP_KERATINIZATION` | -1.556 | 0.090 | 0.018 | 12 |

IFN / MHC-I after CLDN4 loss on this rank: MHC-I **UP** (NES +0.868, FDR 0.335); IFN-γ and IFN-α **DOWN** (FDR 0.185 and 0.251). None of the three are FDR<0.05.

### Mean of both lines (descriptive)

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -1.542 | 0.017 | 0.007 | 115 |
| Hallmark IFN-α | `HALLMARK_INTERFERON_ALPHA_RESPONSE` | -1.637 | 0.010 | 0.002 | 60 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -0.819 | 0.336 | 0.336 | 8 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | -1.223 | 0.102 | 0.082 | 121 |
| GO keratinization | `GOBP_KERATINIZATION` | -1.364 | 0.102 | 0.066 | 12 |

MHC-I / APM on this cuffdiff table is sparse (classical HLA-A/B, B2M, TAP1/2, PSMB8/9 are absent as `gene_short_name`). `n_set_in_rank` is the number that were actually ranked.

---

## 2. GSE50927 — Cldn4 KO mouse **whole lung** (not lung cancer)

GEO: [GSE50927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50927). Author edgeR table `GSE50927_Cldn4lungWTvsKOgenes.csv.gz`. Naive (no ventilator) Cldn4 KO vs WT. **n = 1 vs 1**. Mixed-cell whole lung on a mixed 129S6/C57BL/6/BALB/c background. This is **not** a lung-tumour series.

Cldn4 author logFC = **-6.061** (logCPM 2.251; submitter edgeR FDR 4.07e-26). Sign is KO minus WT.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | +1.508 | 0.005 | 0.001 | 181 |
| Hallmark IFN-α | `HALLMARK_INTERFERON_ALPHA_RESPONSE` | +1.587 | 0.010 | 0.004 | 89 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | +1.423 | 0.050 | 0.030 | 21 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | -0.963 | 0.283 | 0.283 | 154 |
| GO keratinization | `GOBP_KERATINIZATION` | -1.060 | 0.202 | 0.162 | 69 |

IFN / MHC-I after Cldn4 loss: **UP** for all three headline immune sets (IFN-γ FDR 0.005, IFN-α FDR 0.010, MHC-I FDR 0.050). Tight junction and keratinization NES are negative and not FDR<0.05.

Author edgeR *P* / FDR on an unreplicated design assume a dispersion; they are **not** used as the GSEA null. The rank is the deposited logFC.

---

## 3. GSE22493 — SKOV-3 CLDN4 siRNA (mapped, included)

GEO: [GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493), platform [GPL10555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10555). Series-matrix VALUE columns GSM558700–702.

Mapping: **GPL10555 ORF, else DESCRIPTION prefix before '--'**. **27878 / 36284** probes → **18682** symbols. TACSTD2 on platform: **False**.

CLDN4 mean log2(KD / OE) = **-1.225** (per array: GSM558700=NA, GSM558701=-1.737, GSM558702=-0.713). Control arm is CLDN4-overexpressing SKOV-3, not parental WT. n = 3 arrays.

| Set | Term | NES | FDR | nom p | n in rank |
|---|---|---:|---:|---:|---:|
| Hallmark IFN-γ | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | -1.084 | 0.323 | 0.181 | 155 |
| Hallmark IFN-α | `HALLMARK_INTERFERON_ALPHA_RESPONSE` | -0.611 | 0.622 | 0.622 | 70 |
| MHC-I / APM | `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` | -0.852 | 0.468 | 0.375 | 17 |
| KEGG tight junction | `KEGG_TIGHT_JUNCTION` | -1.084 | 0.323 | 0.194 | 139 |
| GO keratinization | `GOBP_KERATINIZATION` | -1.548 | 0.080 | 0.016 | 33 |

IFN / MHC-I after CLDN4 loss on this rank: all three NES are negative; none reach FDR<0.05.

---

## Extra figures

- `methods/cldn4_ko_gsea/figures/fig_headline_nes.png`
- `methods/cldn4_ko_gsea/figures/fig_cldn4_log2fc_and_nes_heatmap.png`

Tables: `methods/cldn4_ko_gsea/tables/gsea_headline.tsv`, `gsea_prerank_all.tsv`, `cldn4_log2fc.tsv`, `tacstd2_log2fc.tsv`, `inventory.tsv`.

---

## What this is not

- Not a re-run of any existing slide or of SKB264.
- Not FASTQ / salmon / DESeq2 from SRA.
- Not a lung-cancer Cldn4 KO. GSE50927 is whole lung.
- Not sample-permutation GSEA. n on GSE207704 and GSE50927 is too small for that; the engine is gene-set permutation on a prerank, same as the scRNA pseudobulk extra.

---

## 中文摘要

只补公开 CLDN4 / Cldn4 缺失的 prerank GSEA，不审不撤已有页。引擎与 `scrna_pseudobulk_gsea_meta` 相同（加权 KS，1000 次基因集置换，seed=42）。正 NES = CLDN4 缺失端富集。

- **GSE207704** 乳腺 T47D / MCF7 CRISPR KO，GEO 只有组均 FPKM（2 vs 2 已合并）。CLDN4 log2FC T47D −1.039、MCF7 −0.750。T47D 上 IFN-γ / IFN-α NES 为负且 FDR<0.05；MCF7 上 MHC-I NES +0.868（FDR 0.335），IFN NES 为负且未过 FDR 0.05。
- **GSE50927** 是小鼠**全肺** Cldn4 KO，不是肺癌。n=1 vs 1。Cldn4 logFC −6.061。IFN-γ / IFN-α / MHC-I 在 Cldn4 缺失端 **UP**（FDR 0.005 / 0.010 / 0.050）。
- **GSE22493** SKOV-3 siRNA 可用 GPL10555 映射到基因符号（18682 个），已纳入。对照是 CLDN4 过表达而非亲本。CLDN4 log2(KD/OE) −1.225。IFN/MHC NES 为负，FDR≥0.05。TACSTD2 不在芯片上。
