# TACSTD2 and CLDN4 with immune/IFN: co-perturbation and dual readout

No public factorial knockout of both TACSTD2 and CLDN4 was found. Where one gene is perturbed and the other is measured in the same matrix, or where an IFN or immune stimulus is applied and both genes are measured, the pair rule below is met in one contrast: SKOV3 shTACSTD2 (GSE245459, no drug, n=3 vs 3). TACSTD2, CLDN4, and the IFN score all fall. The IFN score goes down (mean log2 difference −0.792, Welch p=0.00228). That is co-movement with a lower IFN score, not an IFN increase.

Every other replicated contrast is `null`, or the targeted transcript was not knocked down (`perturbation_failed`). No contrast met `opposite`.

This does not change the locked CosMx, concordant-4, GSE137244, TCGA, or TISMO numbers.

## Rule

A row is a dual readout only when TACSTD2 and CLDN4 (or the mouse symbols) are both in the same contrast, and an IFN score can be computed from at least 8 genes of a fixed panel (`ISG15`, `IFIT1/2/3`, `MX1/2`, `OAS1/2/3`, `IFI6`, `IFI27`, `IFI44`, `IFI44L`, `STAT1`, `IRF1`, `IRF7`, `B2M`, HLA-A/B/C or H2-K1/H2-D1, `TAP1`, `PSMB8/9`, `CD274`, `CXCL10`, `CXCL9`; mouse orthologs as in the script). The score is the mean of log2 expression of the genes present. The p-value is a two-sided Welch test on the per-sample scores, unless the deposit is one value per group.

`supports`: the perturbed gene is down (log2FC < −0.5), the partner is down (p<0.05 and |log2FC|≥0.25), and the IFN score moves (p<0.05 and |delta|≥0.25). For an IFN or dsRNA stimulus, the score must rise and both genes must move in the same direction at that threshold.

`null`: the perturbation or the stimulus is real, but the partner does not clear the threshold, or the IFN score does not.

`opposite`: the partner moves against the perturbed gene while the IFN score moves, or, under a successful IFN stimulus, the two genes move in opposite directions at the threshold.

`perturbation_failed`: the transcript that was supposed to be lost is not down by 0.5 log2. `unreplicated` and `not_evaluable` are reported and are not counted as the three calls.

A second column, `residual_partial_call`, is filled only when the contrast has at least 8 samples. It is a Spearman partial correlation of each gene with the IFN score given the other gene (ranks, linear residual, then Pearson). `supports` there means at least one partial correlation has p<0.05 and every significant partial correlation is positive. It does not replace the pair call.

## Pair calls

| Accession | Contrast | n | TACSTD2 log2FC (p) | CLDN4 log2FC (p) | IFN-score delta (p) | Call |
|---|---|---:|---|---|---|---|
| GSE245459 | SKOV3 shTACSTD2 vs shNC, no drug | 3 vs 3 | −2.628 (5.25×10⁻⁵) | −1.915 (3.29×10⁻⁵) | −0.792 (0.00228) | supports |
| GSE245459 | same knockdown, both arms on cisplatin | 3 vs 3 | −4.033 (4.38×10⁻⁵) | −0.139 (0.0989) | +0.302 (0.0247) | null |
| GSE334497 | 4T1 Trop2 KO vs WT tumor | 5 vs 5 | −3.821 (1.08×10⁻³) | −0.817 (0.253) | +0.176 (0.404) | null |
| GSE289287 | T-47D Trop-2 KO vs WT xenograft | 4 vs 3 | −3.465 (1.13×10⁻⁴) | +0.287 (0.0414) | +0.710 (0.207) | null |
| GSE15212 | SW480 TACSTD2 siRNA vs neg, 72 h | 6 vs 9 | −3.162 (2.15×10⁻¹²) | −0.008 (0.968) | +0.093 (0.0287) | null |
| GSE207704 | CLDN4−/− vs WT, MCF7 and T47D | 2 lines | −0.807 (0.147) | −0.857 (0.111) | −0.453 (0.204) | null |
| GSE50927 | Cldn4 KO vs WT lung, no VILI | 1 vs 1 | −0.145 (author 0.486) | −6.061 (author 1.77×10⁻²⁹) | mean gene logFC +0.713 | unreplicated |
| GSE304294 | KYSE30 IMMU132 vs control, 1 day | 2 vs 3 | +1.159 (0.00987) | +0.911 (3.90×10⁻⁵) | +0.298 (0.00138) | perturbation_failed |
| GSE311016 | CRC PDX IMMU132 vs paired control, day 29 | 5 pairs | −0.421 (0.271) | −0.437 (0.0597) | +0.128 (0.809) | perturbation_failed |
| GSE274940 | EpH4 pan-claudin null vs WT | 3 vs 3 | −0.611 (0.0223) | +0.321 (0.361) | −0.278 (0.282) | perturbation_failed |

Log2FC values are rounded to 3 decimals and p-values to 3 significant figures from `evidence_table.tsv`.

### The one supports row

GSE245459 untreated knockdown: 25 of 26 IFN-panel genes have a negative log2FC. The largest drops are HLA-A −4.166, IFIT1 −1.890, IFIT2 −1.622, IFI6 −1.605, PSMB8 −1.544. CXCL9 is 0 because the gene is absent in these FPKM values. The call is co-movement, and the IFN score falls.

On cisplatin, TACSTD2 is still down (−4.033) but CLDN4 is not (−0.139, p=0.0989), so the pair call is null. The IFN score on that arm is up (+0.302, p=0.0247).

### Null loss-of-function readouts

GSE15212 is the clearest null. Both TACSTD2 siRNAs at 72 h versus the matched negative siRNA drop TACSTD2 by 3.162 log2. CLDN4 changes by −0.008 (p=0.968). The IFN-score delta is +0.093. Its Welch p is 0.0287, and it stays under the 0.25 floor, so it is not an IFN shift.

GSE334497 knocks Trop2 down (−3.821) in immunocompetent 4T1 tumors. Cldn4 is −0.817 with p=0.253, and the IFN score is +0.176 with p=0.404. Columns whose names contain KO are the knockout arm; Tacstd2 is lower in those columns.

GSE289287 knocks Trop-2 down in T-47D xenografts. The sample-level Welch test gives CLDN4 +0.287, p=0.0414. The author DESeq2 table in the same file gives CLDN4 log2FC +0.284, padj=0.466, and TACSTD2 log2FC −3.263, padj=1.95×10⁻⁶¹. The IFN score is +0.710, p=0.207. The pair call is null because the IFN score does not clear p<0.05. The Welch and author tests disagree on CLDN4, and the call does not depend on treating CLDN4 as significant.

GSE207704 deposits one FPKM per genotype per line (the two GEO replicates are already collapsed). Line-level log2FC: TACSTD2 −0.996 and −0.617; CLDN4 −0.706 and −1.008; IFN score −0.302 and −0.604. The two-line t-test is not p<0.05. Call null.

GSE50927 is one WT library and one KO library without ventilation injury. Author edgeR logFC is reported and is not used as a replicate p-value. Cldn4 logFC −6.061. Tacstd2 logFC −0.145, author P=0.486. The mean of 23 deposited IFN-gene logFCs is +0.713. Call unreplicated.

### TROP2 antibody-drug conjugate and the pan-claudin line

IMMU132 does not lower TACSTD2 mRNA in the KYSE30 matrix (GSE304294, +1.159). CLDN4 and the IFN score rise with it. The loss rule therefore returns `perturbation_failed`. OX2 is the IMMU132 block because it is the only block with two replicates. OX1 is the three-replicate block placed before it. As a mapping check that was not part of the call, CDKN1A mean FPKM is 43.3 in OX1 and 128.8 in OX2.

The day-29 PDX contrast (GSE311016, paired by PDX id) moves TACSTD2 by −0.421 (p=0.271), which misses the −0.5 loss cutoff. CLDN4 is −0.437 (p=0.0597). The IFN score is flat.

GSE274940 is an EpH4 pan-claudin-null line. Cldn4 RNA is not down (+0.321, p=0.361), so it is not a CLDN4-loss contrast.

## IFN and immune stimuli (both genes as readouts)

The IFN score rises in every stimulus row below. Neither gene, or only one gene, clears p<0.05 and |log2FC|≥0.25. Pair call is null for all of them. No stimulus row is `opposite`.

| Accession | Stimulus | n | TACSTD2 | CLDN4 | IFN score |
|---|---|---:|---|---|---|
| GSE156295 | A549 type I IFN vs mock | 2 vs 2 | −0.035 (0.500) | −0.901 (0.0650) | +4.061 (0.00310) |
| GSE156295 | HTBE type I IFN vs mock | 2 vs 2 | +0.536 (0.145) | −0.145 (0.428) | +5.036 (5.17×10⁻⁴) |
| GSE215771 | A549 WT IFN-γ vs control | 3 vs 3 | −0.050 (0.423) | −0.439 (0.0790) | +3.345 (3.09×10⁻⁴) |
| GSE306855 | HSC-2 IFN-γ vs PBS | 3 vs 3 | +0.199 (0.0744) | −0.589 (0.00435) | +1.962 (7.46×10⁻⁵) |
| GSE184456 | gingival keratinocyte IFN-β | 3 vs 3 | +0.186 (0.0257) | −0.194 (0.00872) | +2.151 (0.00107) |
| GSE184456 | same cells, IFN-λ | 3 vs 3 | +0.042 (0.368) | −0.124 (0.159) | +3.318 (1.99×10⁻⁴) |
| GSE241852 | MCF10A Poly(I:C), NEO | 3 vs 3 | +0.027 (0.0725) | +0.363 (1.20×10⁻⁴) | +0.838 (3.46×10⁻⁶) |
| GSE241852 | MCF10A Poly(I:C), GRHL2KO | 3 vs 3 | +0.528 (6.42×10⁻⁴) | +0.208 (0.0382) | +0.465 (3.30×10⁻⁵) |
| GSE190899 | IEC organoids IFN-β | 3 lines | +0.236 (0.202) | −0.047 (0.544) | +3.979 (0.00105) |
| GSE190899 | IEC organoids IFN-γ | 2 lines | −0.102 (0.527) | −0.067 (0.544) | +2.389 (2.89×10⁻⁴) |
| GSE190899 | IEC organoids IFN-λ | 3 lines | +0.191 (0.253) | +0.043 (0.370) | +2.640 (8.56×10⁻⁵) |
| GSE239485 | Poly(I:C)+anti-PD-1 vs vehicle | 8 vs 8 | +2.093 (6.31×10⁻⁴) | −0.771 (0.238) | +1.043 (5.93×10⁻⁴) |
| GSE239485 | plus anti-C5aR1 vs vehicle | 8 vs 8 | +2.576 (9.74×10⁻⁴) | +0.586 (0.357) | +1.410 (6.17×10⁻⁴) |

GSE306855 and GSE184456 IFN-β move the two genes in opposite directions, but at least one gene misses the 0.25 floor, so the call stays null. GSE241852 GRHL2KO moves both genes up; CLDN4 log2FC is +0.208, under the floor. GSE190899 tests line-level differences (005 ileum, 005 sigmoid, 007 sigmoid). Letter codes B/G/L/N were accepted only after ISG15 rose on each IFN arm (all checked arms, ISG15 log2FC > 0.5).

GSE239485 is the immune-treatment dual readout. Tacstd2 rises with the IFN score. Cldn4 does not (p=0.238 and p=0.357). Pair call null.

## Residual correlation where n≥8

| Contrast | n | partial ρ TACSTD2 vs IFN given CLDN4 (p) | partial ρ CLDN4 vs IFN given TACSTD2 (p) | Residual call |
|---|---:|---|---|---|
| GSE239485 Poly(I:C)+anti-PD-1 | 16 | +0.672 (0.00435) | −0.053 (0.846) | supports |
| GSE239485 plus anti-C5aR1 | 16 | +0.633 (0.00849) | −0.182 (0.501) | supports |
| GSE334497 Trop2 KO and WT | 10 | −0.294 (0.410) | +0.262 (0.464) | null |
| GSE15212 siRNA and controls | 15 | −0.511 (0.0516) | +0.424 (0.116) | null |
| GSE190899 IFN-β and untreated | 21 | +0.225 (0.327) | −0.231 (0.313) | null |
| GSE190899 IFN-γ and untreated | 12 | +0.329 (0.296) | −0.165 (0.609) | null |
| GSE190899 IFN-λ and untreated | 20 | +0.268 (0.253) | −0.195 (0.409) | null |

In GSE239485 the IFN score stays positively associated with Tacstd2 after Cldn4 is residualized. Cldn4 has no residual association with the score. That is a Tacstd2 residual, and it is not a claim that the two genes form a pair. The GSE15212 TACSTD2 partial correlation is −0.511 with p=0.0516, which misses 0.05.

## What was searched and not scored

Counts and exclusions are in `search_audit.tsv` and `gaps.tsv`. There is no deposited double-knockout RNA-seq of the two genes. Protein papers (PMID 20651236) discuss CLDN4 localization without an IFN transcriptome and without a GEO series. PMID 17651017 measures CLDN4 siRNA and IFN biology without TACSTD2, and it has no GEO link. PMID 35688908 has no GEO link.

## Reproduce

Processed matrices are read from `/tmp/sweep/raw` (GEO FTP supplementary files and GPL4133 annotation). They are not committed.

```bash
python3 scripts/partial_dependence/analyze_partial_dependence.py
python3 scripts/partial_dependence/plot_evidence.py
```

Read `evidence_table.tsv` with `keep_default_na=False`. The call `null` is the word null.

## 中文

公开数据里没有 TACSTD2 与 CLDN4 的双敲转录组。同一对比里两个基因和 IFN 分数都能量到时，按预定规则只有 GSE245459 未加药的 shTACSTD2 达到 supports：TACSTD2、CLDN4 和 IFN 分数一起下降（IFN 分数 −0.792，p=0.00228）。其余有重复的对比是 null，或靶基因转录本没有下降。没有 opposite。GSE239485 里 Tacstd2 在校正 Cldn4 之后仍与 IFN 分数正相关（partial ρ 0.672 和 0.633）；Cldn4 的偏相关不显著。这是 Tacstd2 的剩余相关，不是两个基因成对。
