# ADDITIVE A8 extra — keratin / TJ / EMT GSEA in public ICI-treated lung tumors

**Public data only. Additive to the TCGA A8 slide. That TCGA audit was not re-run.**

User A8 (TROP2-high enrich keratin/TJ, EMT down, conserved) is **taken as given** for TCGA-LUAD/LUSC. This file asks the same locked gene-set question in leftover open **pre-treatment ICI-treated** lung tumor bulk.

## 一句话结论 / TL;DR

TACSTD2-high **keratin / TJ-up is common** in leftover ICI bulk (same direction as A8 TCGA). Hallmark EMT-down is **not** conserved: 4/6 median-split cohorts have Hallmark EMT significantly **up**; the other two are null. **No** median-split leftover cohort is `supportive`.

Per-cohort median verdicts: GSE126044=keratin_TJ_up_Hallmark_EMT_null (n_high=8, n_low=8); GSE135222=keratin_TJ_up_Hallmark_EMT_opposite (n_high=13, n_low=13); GSE166449=keratin_TJ_up_Hallmark_EMT_opposite (n_high=11, n_low=11); GSE253564=keratin_TJ_up_Hallmark_EMT_opposite (n_high=16, n_low=16); GSE190265=keratin_TJ_up_Hallmark_EMT_opposite (n_high=21, n_low=21); GSE283829=keratin_TJ_up_Hallmark_EMT_null (n_high=13, n_low=13).

n=16–43 is **hypothesis-generating**. This does not retract the TCGA A8 slide. Do not quote a pooled ICI NES. Do not write “A8 conserved in ICI tumors.” Hallmark EMT, not GOBP EMT, decides the EMT arm. The single `supportive` call is GSE283829 **quartile** 7 vs 7 — not the locked median split.

## Why this extra exists

The original A8 slide is TCGA-only (ICI-naive surgical resections). Conservation of a keratin/TJ-high, EMT-low TROP2-high state in **tumors that actually received PD-(L)1** is a different question. This leftover slice uses open GEO ICI lung bulk that was not that TCGA matrix.

Excluded on purpose: TCGA-LUAD/LUSC (A8, not re-run); GSE207422 (already GSEA'd in `hunt_tj_gsea`); GSE248378 (post-treatment FPKM); immune-only panels without TACSTD2 (GSE136961, GSE93157); FASTQ/SRA and controlled OAK/POPLAR.

## Cohorts

| Cohort | n | R / NR | Scale | Usable | Note |
|---|---:|---|---|---|---|
| GSE126044 | 16 | 5 / 11 | log2(CPM+1) from deposited counts | True | core open ICI bulk; not in the TCGA A8 slide |
| GSE135222 | 27 | 7 / 20 | log2(TPM+1) from deposited TPM; Ensembl mapped via HGNC | True | core open ICI bulk; DCB is GEO PFS>=180 d |
| GSE166449 | 22 | 7 / 15 | deposited matrix kept if already log2, else log2(TPM+1) | True | core open ICI bulk; response from GEO sample titles |
| GSE253564 | 32 | 0 / 0 | log2(FPKM+1) from deposited pre-treatment FPKM | True | leftover whole-transcriptome ICI tumor matrix; MPR not deposited on GEO |
| GSE190265 | 43 | 14 / 29 | log2(TPM+1) from deposited France3 TPM | True | France4 (GSE190266) omitted: TACSTD2 absent from the deposited TPM cap |
| GSE283829 | 27 | 7 / 10 | log2(CPM+1) from supplementary raw counts; Ensembl mapped via HGNC | True | 2025 leftover; R vs NR is CR vs PD (SD held out) |

## Pre-specified design

| Piece | Choice | Honest limitation |
|---|---|---|
| Split (primary) | TACSTD2 above vs below median | Quartiles on n=16 are 4 vs 4. Median is the locked small-n split. |
| Split (sensitivity) | Quartile only if both arms ≥ 6 | Still discards the middle; underpowered. |
| Ranking | Welch t, high − low | Two-group statistic on n≈8–16/arm. |
| Complementary | Spearman gene vs TACSTD2 prerank | Uses every sample; not a two-group t. |
| Response | R vs NR Welch t if both arms ≥ 5 | Exploratory. Deposited labels only. |
| GSEA | Preranked weighted KS, p=1, 1000 gene-set permutations, seed=42 | Same engine as A8. Gene-set permutation, not sample permutation. |
| FDR | BH within the 12 primary sets per contrast | Not nested GSEA FDR. |
| EMT rule | Hallmark EMT decides the claim | GOBP EMT is recorded and must not be quoted as Hallmark. |
| Sets | Identical `data/genesets/a8_sets.json` primary 12 | TACSTD2 is not a member of the primary sets. |

Positive NES = enriched in TACSTD2-high (or in responders). NES is **NA** when no same-sign null ES exists (the ratio is undefined); still report ES and nominal p from `gsea_prerank_all.tsv`.

## Primary NES — TACSTD2 median split

| Gene set | Want | GSE126044 NES | GSE126044 FDR | GSE135222 NES | GSE135222 FDR | GSE166449 NES | GSE166449 FDR | GSE253564 NES | GSE253564 FDR | GSE190265 NES | GSE190265 FDR | GSE283829 NES | GSE283829 FDR |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | DOWN | +0.923 | 0.473 | +1.512 | 0.004 | +1.416 | 0.006 | +2.776 | 0.001 | +1.311 | 0.007 | -1.059 | 0.276 |
| GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION | DOWN | -0.905 | 0.252 | +1.592 | 0.016 | +1.024 | 0.497 | +1.594 | 0.008 | +1.593 | 0.001 | -1.013 | 0.297 |
| HALLMARK_APICAL_JUNCTION | UP | +1.354 | 0.030 | +1.437 | 0.016 | +1.530 | 0.006 | +2.127 | 0.001 | +1.483 | 0.001 | +2.061 | 0.001 |
| KEGG_TIGHT_JUNCTION | UP | +1.291 | 0.052 | +1.349 | 0.034 | +1.491 | 0.008 | +2.043 | 0.001 | +1.638 | 0.001 | +1.794 | 0.001 |
| GOBP_KERATINIZATION | UP | +1.996 | 0.002 | +2.346 | 0.002 | +0.757 | 0.914 | +2.813 | 0.001 | +1.124 | 0.277 | +2.862 | 0.001 |
| KRT_EPITHELIAL | UP | +2.164 | 0.002 | +2.262 | 0.002 | +1.541 | 0.043 | +3.201 | 0.001 | +1.777 | 0.001 | +2.645 | 0.001 |

Full 12-set tables: `tables/gsea_prerank_all.tsv`.

| Cohort | n_high | n_low | Verdict |
|---|---:|---:|---|
| GSE126044 | 8 | 8 | keratin_TJ_up_Hallmark_EMT_null |
| GSE135222 | 13 | 13 | keratin_TJ_up_Hallmark_EMT_opposite |
| GSE166449 | 11 | 11 | keratin_TJ_up_Hallmark_EMT_opposite |
| GSE253564 | 16 | 16 | keratin_TJ_up_Hallmark_EMT_opposite |
| GSE190265 | 21 | 21 | keratin_TJ_up_Hallmark_EMT_opposite |
| GSE283829 | 13 | 13 | keratin_TJ_up_Hallmark_EMT_null |

## Sensitivity — TACSTD2 quartile (only if both arms ≥ 6)

| Gene set | Want | GSE135222 NES | GSE135222 FDR | GSE166449 NES | GSE166449 FDR | GSE253564 NES | GSE253564 FDR | GSE190265 NES | GSE190265 FDR | GSE283829 NES | GSE283829 FDR |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | DOWN | +0.843 | 0.776 | +1.188 | 0.141 | +2.645 | 0.001 | +0.694 | 0.988 | -1.354 | 0.015 |
| GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION | DOWN | +0.847 | 0.666 | +0.493 | 0.915 | +1.812 | 0.002 | +1.127 | 0.313 | -1.201 | 0.168 |
| HALLMARK_APICAL_JUNCTION | UP | +1.131 | 0.294 | +1.544 | 0.002 | +2.015 | 0.001 | +1.436 | 0.001 | +1.946 | 0.001 |
| KEGG_TIGHT_JUNCTION | UP | +1.022 | 0.505 | +1.828 | 0.002 | +1.938 | 0.001 | +1.616 | 0.001 | +2.212 | 0.001 |
| GOBP_KERATINIZATION | UP | +2.141 | 0.006 | +1.601 | 0.003 | +2.595 | 0.001 | +2.005 | 0.001 | +3.366 | 0.001 |
| KRT_EPITHELIAL | UP | +2.281 | 0.006 | +2.216 | 0.002 | +2.760 | 0.001 | +2.284 | 0.001 | +3.146 | 0.001 |

Full 12-set tables: `tables/gsea_prerank_all.tsv`.

| Cohort | n_high | n_low | Verdict |
|---|---:|---:|---|
| GSE126044 | 4 | 4 | n_too_small |
| GSE135222 | 7 | 7 | mixed |
| GSE166449 | 6 | 6 | keratin_TJ_up_Hallmark_EMT_null |
| GSE253564 | 8 | 8 | keratin_TJ_up_Hallmark_EMT_opposite |
| GSE190265 | 11 | 11 | keratin_TJ_up_Hallmark_EMT_null |
| GSE283829 | 7 | 7 | supportive |

## Complementary — Spearman prerank vs continuous TACSTD2

| Gene set | Want | GSE126044 NES | GSE126044 FDR | GSE135222 NES | GSE135222 FDR | GSE166449 NES | GSE166449 FDR | GSE253564 NES | GSE253564 FDR | GSE190265 NES | GSE190265 FDR | GSE283829 NES | GSE283829 FDR |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | DOWN | +1.002 | 0.426 | +0.770 | 0.874 | +1.441 | 0.003 | +2.190 | 0.001 | +0.737 | 0.993 | -1.174 | 0.137 |
| GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION | DOWN | +0.910 | 0.437 | +0.819 | 0.676 | +0.897 | 0.604 | +1.529 | 0.010 | +1.069 | 0.401 | -1.249 | 0.122 |
| HALLMARK_APICAL_JUNCTION | UP | +1.446 | 0.024 | +1.174 | 0.192 | +1.727 | 0.001 | +1.750 | 0.001 | +1.302 | 0.010 | +2.031 | 0.001 |
| KEGG_TIGHT_JUNCTION | UP | +1.286 | 0.090 | +1.142 | 0.261 | +1.692 | 0.001 | +1.948 | 0.001 | +1.506 | 0.002 | +1.772 | 0.001 |
| GOBP_KERATINIZATION | UP | +1.194 | 0.198 | +2.654 | 0.004 | +1.593 | 0.001 | +2.764 | 0.001 | +1.625 | 0.002 | +3.025 | 0.001 |
| KRT_EPITHELIAL | UP | +2.098 | 0.003 | +2.261 | 0.004 | +2.271 | 0.001 | +2.827 | 0.001 | +1.831 | 0.002 | +2.849 | 0.001 |

Full 12-set tables: `tables/gsea_prerank_all.tsv`.

| Cohort | n_high | n_low | Verdict |
|---|---:|---:|---|
| GSE126044 | 16 | 16 | keratin_TJ_up_Hallmark_EMT_null |
| GSE135222 | 27 | 27 | mixed |
| GSE166449 | 22 | 22 | keratin_TJ_up_Hallmark_EMT_opposite |
| GSE253564 | 32 | 32 | keratin_TJ_up_Hallmark_EMT_opposite |
| GSE190265 | 43 | 43 | keratin_TJ_up_Hallmark_EMT_null |
| GSE283829 | 27 | 27 | keratin_TJ_up_Hallmark_EMT_null |

## Exploratory — responder vs NR (only if both arms ≥ 5)

| Gene set | Want | GSE126044 NES | GSE126044 FDR | GSE135222 NES | GSE135222 FDR | GSE166449 NES | GSE166449 FDR | GSE190265 NES | GSE190265 FDR | GSE283829 NES | GSE283829 FDR |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | DOWN | +2.720 | 0.004 | -1.725 | 0.002 | +0.962 | 0.258 | +2.172 | 0.004 | NA | 0.003 |
| GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION | DOWN | +1.100 | 0.170 | -2.069 | 0.002 | -0.770 | 0.832 | +1.532 | 0.030 | -0.822 | 0.937 |
| HALLMARK_APICAL_JUNCTION | UP | +1.774 | 0.004 | -1.507 | 0.003 | -0.663 | 0.864 | +1.575 | 0.004 | NA | 0.003 |
| KEGG_TIGHT_JUNCTION | UP | -1.126 | 0.192 | -1.446 | 0.002 | +0.873 | 0.319 | +1.464 | 0.015 | NA | 0.003 |
| GOBP_KERATINIZATION | UP | -1.819 | 0.006 | -2.063 | 0.002 | -1.745 | 0.012 | -2.630 | 0.004 | -0.737 | 0.937 |
| KRT_EPITHELIAL | UP | -1.205 | 0.181 | -1.837 | 0.002 | -1.427 | 0.258 | -1.708 | 0.022 | +1.208 | 0.030 |

Full 12-set tables: `tables/gsea_prerank_all.tsv`.

| Cohort | n_high | n_low | Verdict |
|---|---:|---:|---|
| GSE126044 | 5 | 11 | mixed |
| GSE135222 | 7 | 20 | contradicts_TJ_KRT |
| GSE166449 | 7 | 15 | null |
| GSE253564 | 0 | 0 | no_deposited_response_label |
| GSE190265 | 14 | 29 | mixed |
| GSE283829 | 7 | 10 | mixed |

## Complementary signature Spearman (no split)

| Gene set | GSE126044 ρ | GSE135222 ρ | GSE166449 ρ | GSE253564 ρ | GSE190265 ρ | GSE283829 ρ |
|---|---:|---:|---:|---:|---:|---:
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | +0.221 | +0.057 | +0.247 | +0.132 | +0.181 | -0.143 |
| HALLMARK_APICAL_JUNCTION | +0.441 | +0.118 | +0.373 | +0.112 | +0.420 | +0.024 |
| KEGG_TIGHT_JUNCTION | +0.265 | +0.272 | +0.351 | +0.260 | +0.543 | +0.154 |
| GOBP_KERATINIZATION | +0.271 | +0.640 | +0.442 | +0.429 | +0.510 | +0.281 |
| KRT_EPITHELIAL | +0.638 | +0.677 | +0.794 | +0.794 | +0.706 | +0.607 |

## What this does **not** say

- It does not re-compute or revise the TCGA A8 NES table.
- It does not claim a predictive TROP2 × ICI interaction. Association of a program with TACSTD2 inside a treated cohort is not a treatment-effect modifier.
- Small n cannot rule out a modest true effect. Report the observed NES/FDR/n.
- Do not collapse Hallmark EMT and GOBP EMT.
- Do not write “conserved in ICI tumors” unless the median-split Hallmark EMT, a TJ set, and a keratin set all pass FDR<0.05 in the same leftover cohort.

## Reproduce

```bash
pip install -r requirements.txt
bash scripts/a8_ici_gsea/download.sh /tmp/a8_ici_gsea_data
python3 scripts/a8_ici_gsea/analyze.py
```

Methods: `methods/a8_ici_gsea/playbook.md`.
