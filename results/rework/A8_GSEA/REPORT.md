# REWORK A8 — TACSTD2-high keratin/TJ GSEA and EMT-down in TCGA lung

**Self-contained. Public data only. Written to be read without the rest of the repo.**

**Verdict: keratin_TJ_up_Hallmark_EMT_opposite in LUAD; supportive in LUSC (quartile split, BH-FDR<0.05 within the 12 primary sets).**

Keratin / tight-junction programs are enriched in TACSTD2-high tumors in **both** histologies (keratinization LUAD YES, LUSC YES; KEGG TJ LUAD YES, LUSC YES). Hallmark EMT-down is **not** a shared fact: LUSC YES (NES −1.35), LUAD no — Hallmark EMT is significantly **up** in TACSTD2-high LUAD (NES +1.55). The smaller GOBP EMT process set is down in both. Do not write “TROP2-high = EMT-down” for LUAD. Do not quote a pooled NSCLC NES.

## Why this rework exists

User claim A8 said GSEA of TROP2-high vs TROP2-low tumors shows
**keratin / tight-junction programs up** and **EMT down**. That is an
epithelial-identity claim. It can be true as differentiation, true as a
tight-junction module, or true only in squamous tumors (LUSC) where
keratin programs are already the default. This file tests TCGA-LUAD and
TCGA-LUSC separately, reports the actual NES, and does not hide nulls.

## Analysis set

| Cohort | Primary tumors | TACSTD2 Q4 / Q1 | TACSTD2 median (log2 RSEM+1) | TACSTD2 vs ABSOLUTE ρ |
|---|---:|---:|---:|---|
| TCGA-LUAD | 515 | 129 / 129 | 12.63 | +0.007 (p=0.869, n=502) |
| TCGA-LUSC | 502 | 126 / 126 | 13.11 | -0.074 (p=0.102, n=493) |

Primary tumors only (`-01`). One row per 15-character barcode.
Xena HiSeqV2 is log2(RSEM normalized_count + 1).

## Pre-specified design

| Piece | Choice | Honest limitation |
|---|---|---|
| Split | TACSTD2 top vs bottom quartile (primary); median (sensitivity) | Quartiles discard the middle half. That is the usual high-vs-low GSEA design, not a continuous test. |
| Ranking | Welch t, high − low | Two-group statistic. Complementary Spearman of signature scores uses the full cohort. |
| GSEA | Preranked weighted KS, p=1, 1000 gene-set permutations, seed=42 | Gene-set permutation, not sample permutation. NES is comparable within this run, not to Broad GSEA GUI output. |
| FDR | BH within the 12 primary sets (`fdr_primary`); also BH across all ~60 sets (`fdr_all`) | Classical GSEA nested FDR is not used. BH is more transparent. |
| Hallmark | Enrichr MSigDB_Hallmark_2020 = Liberzon 2015 Hallmark 50 | Same collection as MSigDB Hallmark; gene lists can differ by a few symbols from a later MSigDB freeze. |
| KRT / TJ | KEGG 2021 Tight junction; GO BP 2023 keratin/TJ/barrier; MSigDB v2023.2.Hs keratinization / cornification / TJ organization; compact KRT panel | The compact KRT panel is a fixed cytokeratin list, not an MSigDB set. |
| Purity | ABSOLUTE reported vs TACSTD2; **not** used to split | If TACSTD2 were just purity, EMT-down would be “more tumor / less stroma”. |

Positive NES = enriched in TACSTD2-high. Negative NES = depleted in TACSTD2-high
(enriched in TACSTD2-low). TACSTD2 is not a member of the primary sets.

## Direct answer

- **LUAD verdict:** `keratin_TJ_up_Hallmark_EMT_opposite` at FDR<0.05 on primary sets.
- **LUSC verdict:** `supportive` at FDR<0.05 on primary sets.
- **Hallmark EMT (the claim set):** LUAD up (NES=+1.546 p=0.001998 FDR=0.002664); LUSC down (NES=-1.352 p=0.01598 FDR=0.01598).
- **GOBP EMT (secondary process set):** LUAD down (NES=-1.383 p=0.02697 FDR=0.02697); LUSC down (NES=-1.693 p=0.004995 FDR=0.005994).
- **GOBP keratinization up:** LUAD YES; LUSC YES.
- **KEGG tight junction up:** LUAD YES; LUSC YES.
- **Do not quote a pooled NSCLC NES.** LUAD and LUSC are different diseases.
- **Do not collapse Hallmark EMT and GOBP EMT.** They disagree in LUAD.

## Primary NES (quartile split)

| Gene set | Bucket | Want | LUAD NES | LUAD FDR | LUSC NES | LUSC FDR |
|---|---|---|---:|---:|---:|---:|
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | EMT | DOWN | +1.546 | 0.003 | -1.352 | 0.016 |
| GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION | EMT | DOWN | -1.383 | 0.027 | -1.693 | 0.006 |
| HALLMARK_APICAL_JUNCTION | TJ | UP | +2.193 | 0.001 | +1.600 | 0.001 |
| KEGG_TIGHT_JUNCTION | TJ | UP | +2.036 | 0.001 | +1.698 | 0.001 |
| GOBP_TIGHT_JUNCTION_ORGANIZATION | TJ | UP | +1.735 | 0.001 | +1.780 | 0.001 |
| GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY | TJ | UP | +1.712 | 0.006 | +1.800 | 0.001 |
| GOBP_KERATINIZATION | KERATIN_BARRIER | UP | +2.699 | 0.001 | +3.182 | 0.001 |
| GOBP_CORNIFICATION | KERATIN_BARRIER | UP | +1.641 | 0.010 | +1.741 | 0.008 |
| GOBP_ESTABLISHMENT_OF_SKIN_BARRIER | KERATIN_BARRIER | UP | +2.293 | 0.001 | +2.290 | 0.001 |
| GOBP_KERATINOCYTE_DIFFERENTIATION | KERATIN_BARRIER | UP | +2.553 | 0.001 | +2.680 | 0.001 |
| GOBP_EPIDERMAL_CELL_DIFFERENTIATION | KERATIN_BARRIER | UP | +2.348 | 0.001 | +2.587 | 0.001 |
| KRT_EPITHELIAL | KERATIN_BARRIER | UP | +2.337 | 0.001 | +2.774 | 0.001 |

Nominal p, ES, set size, mean Welch t, and leading-edge genes are in
`gsea_prerank_all.tsv`. The four headline rows in prose:

- **HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION**
  - LUAD: NES=+1.546  ES=+0.327  nom p=0.002  FDR_primary=0.003  n=194  mean_t=+0.53
  - LUSC: NES=-1.352  ES=-0.298  nom p=0.016  FDR_primary=0.016  n=194  mean_t=-1.28
- **HALLMARK_APICAL_JUNCTION**
  - LUAD: NES=+2.193  ES=+0.460  nom p=0.001  FDR_primary=0.001  n=194  mean_t=+1.13
  - LUSC: NES=+1.600  ES=+0.312  nom p=0.001  FDR_primary=0.001  n=194  mean_t=+0.01
- **KEGG_TIGHT_JUNCTION**
  - LUAD: NES=+2.036  ES=+0.447  nom p=0.001  FDR_primary=0.001  n=163  mean_t=+1.05
  - LUSC: NES=+1.698  ES=+0.345  nom p=0.001  FDR_primary=0.001  n=163  mean_t=+0.49
- **GOBP_KERATINIZATION**
  - LUAD: NES=+2.699  ES=+0.674  nom p=0.001  FDR_primary=0.001  n=74  mean_t=+2.32
  - LUSC: NES=+3.182  ES=+0.727  nom p=0.001  FDR_primary=0.001  n=78  mean_t=+4.12

### How to read the verdict labels

- `supportive`: Hallmark EMT down **and** a primary TJ set up **and** a primary KRT/barrier set up (FDR<0.05).
- `keratin_TJ_up_Hallmark_EMT_opposite`: TJ+KRT up, but Hallmark EMT is significantly **up** (claim EMT arm fails).
- `keratin_TJ_up_Hallmark_EMT_null`: TJ+KRT up, Hallmark EMT not significant.
- `partial`: Hallmark EMT-down plus TJ-up **or** KRT-up, not both.
- `mixed`: only one of the three arms.
- `null`: none of the three arms at FDR<0.05.
- `contradicts_TJ_KRT`: TJ and KRT/barrier are significantly **down** and none is up.
- GOBP EMT is recorded but **does not** decide the verdict. Hallmark EMT is the user-claim set.

## Complementary: signature Spearman (no split)

z-mean of set members vs continuous TACSTD2. This does not throw away the middle half.

| Gene set | LUAD ρ | LUAD q | LUSC ρ | LUSC q |
|---|---:|---:|---:|---:|
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | +0.058 | 0.256 | -0.104 | 0.034 |
| GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION | -0.138 | 0.004 | -0.313 | 4.92e-12 |
| HALLMARK_APICAL_JUNCTION | +0.188 | 6.76e-05 | +0.015 | 0.797 |
| KEGG_TIGHT_JUNCTION | +0.297 | 5.19e-11 | +0.144 | 0.003 |
| GOBP_TIGHT_JUNCTION_ORGANIZATION | +0.185 | 8.79e-05 | +0.181 | 0.000 |
| GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY | +0.223 | 1.48e-06 | +0.223 | 1.58e-06 |
| GOBP_KERATINIZATION | +0.375 | 1.62e-17 | +0.396 | 2.22e-19 |
| GOBP_CORNIFICATION | +0.289 | 1.87e-10 | +0.239 | 2.35e-07 |
| GOBP_ESTABLISHMENT_OF_SKIN_BARRIER | +0.464 | 4.53e-27 | +0.491 | 4.28e-30 |
| GOBP_KERATINOCYTE_DIFFERENTIATION | +0.439 | 3.25e-24 | +0.489 | 5.36e-30 |
| GOBP_EPIDERMAL_CELL_DIFFERENTIATION | +0.395 | 2.32e-19 | +0.458 | 3.46e-26 |
| KRT_EPITHELIAL | +0.388 | 8.02e-19 | +0.477 | 1.47e-28 |

## Focal genes vs TACSTD2 (Spearman)

| Gene | Class | LUAD ρ | LUAD q | LUSC ρ | LUSC q |
|---|---|---:|---:|---:|---:|
| CLDN1 | TJ | +0.374 | 9.43e-18 | +0.424 | 4.36e-22 |
| CLDN4 | TJ | +0.460 | 2.28e-27 | +0.410 | 5.06e-21 |
| CLDN7 | TJ | +0.054 | 0.291 | +0.401 | 3.16e-20 |
| F11R | TJ | +0.182 | 6.88e-05 | +0.248 | 4.41e-08 |
| PARD3 | TJ | -0.020 | 0.692 | +0.151 | 0.001 |
| OCLN | TJ | +0.047 | 0.346 | +0.175 | 0.000 |
| TJP1 | TJ | +0.036 | 0.470 | +0.162 | 0.000 |
| KRT5 | squamous KRT | +0.175 | 0.000 | +0.368 | 5.07e-17 |
| KRT17 | squamous KRT | +0.293 | 3.61e-11 | +0.229 | 4.56e-07 |
| KRT7 | simple KRT | +0.348 | 1.70e-15 | +0.097 | 0.033 |
| KRT8 | simple KRT | +0.193 | 2.56e-05 | +0.041 | 0.358 |
| KRT18 | simple KRT | +0.118 | 0.013 | +0.043 | 0.353 |
| KRT19 | simple KRT | +0.475 | 3.68e-29 | +0.419 | 7.63e-22 |
| CDH1 | epithelial | +0.345 | 2.72e-15 | +0.206 | 5.73e-06 |
| VIM | EMT | +0.103 | 0.030 | -0.224 | 7.72e-07 |
| ZEB1 | EMT | -0.066 | 0.189 | -0.278 | 6.81e-10 |
| SNAI2 | EMT | +0.008 | 0.863 | +0.098 | 0.033 |

## Sensitivity: median split NES

Same ranking and GSEA, TACSTD2 above vs below median. If the quartile NES
is a tail artifact, the median NES should shrink or flip.

| Gene set | LUAD NES (median) | LUAD FDR | LUSC NES (median) | LUSC FDR |
|---|---:|---:|---:|---:|
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | +1.948 | 0.001 | -1.290 | 0.032 |
| HALLMARK_APICAL_JUNCTION | +2.163 | 0.001 | +1.638 | 0.001 |
| KEGG_TIGHT_JUNCTION | +2.028 | 0.001 | +1.625 | 0.001 |
| GOBP_KERATINIZATION | +2.643 | 0.001 | +3.061 | 0.001 |
| KRT_EPITHELIAL | +2.361 | 0.001 | +2.741 | 0.001 |

## Hallmark context (not the claim)

If TACSTD2-high is just “more epithelial / less stromal”, many Hallmark
sets will move, not only EMT and apical junction. The full Hallmark 50
NES heatmap is `figures/nes_hallmark_heatmap.png`. Top and bottom 5 by
mean NES (LUAD+LUSC, quartile):

| Rank | Hallmark set | mean NES |
|---|---|---:|
| top 1 | HALLMARK_P53_PATHWAY | +2.520 |
| top 2 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | +2.165 |
| top 3 | HALLMARK_APOPTOSIS | +2.107 |
| top 4 | HALLMARK_INTERFERON_ALPHA_RESPONSE | +2.035 |
| top 5 | HALLMARK_ESTROGEN_RESPONSE_LATE | +2.022 |
| bottom 1 | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | -1.499 |
| bottom 2 | HALLMARK_SPERMATOGENESIS | -1.665 |
| bottom 3 | HALLMARK_G2M_CHECKPOINT | -1.967 |
| bottom 4 | HALLMARK_PANCREAS_BETA_CELLS | -1.992 |
| bottom 5 | HALLMARK_E2F_TARGETS | -2.137 |

TACSTD2-high is not a quiet barrier-only state. Across both histologies
the strongest Hallmark enrichments include p53, TNF-α, apoptosis, and
interferon; the strongest depletions include E2F / G2M / MYC targets.
That is stress / interferon / less-proliferative, not “EMT off, done”.

## Honest interpretation

1. **Headline.** Keratin/TJ-up holds in both LUAD and LUSC. Hallmark
   EMT-down holds in LUSC only. In LUAD, Hallmark EMT is significantly
   **enriched** in TACSTD2-high (NES +1.55, FDR 0.003). The continuous
   Hallmark-EMT z-mean Spearman in LUAD is near zero (see table). GSEA
   “up” there is a leading-edge / tail effect (LAMC2, SDC1, TNC, TGFBI),
   not a cohort-wide VIM/FN1/COL mesenchymal program. Still: the claim
   “EMT down” is **false for Hallmark EMT in LUAD**.
2. **Two EMT sets, two answers.** Hallmark EMT is an ECM/stromal module.
   GOBP EMT is a smaller TF/signaling process set and is down in both
   histologies. Quoting only GOBP EMT would manufacture support for A8
   in LUAD. The user claim named Hallmark EMT. That is the number to quote.
3. **EMT-down is still partly tautological where it occurs.** TACSTD2
   is an epithelial surface gene. LUSC Hallmark EMT-down (NES −1.35) is
   modest. VIM/ZEB1 track TACSTD2 negatively in LUSC, not in LUAD.
4. **Hallmark apical junction is a mixed set.** It contains claudins and
   also mesenchymal/immune junction genes (VCAN, VCAM1, THY1, PTPRC).
   A weak or null apical-junction NES is not a failed tight-junction test.
   KEGG tight junction and the GO TJ sets are the cleaner TJ readouts.
5. **Keratinization is a squamous/cornified set, but it is up in LUAD too.**
   GOBP_KERATINIZATION is skin-barrier keratins (KRT1/5/6/16/17, SPRRs,
   LCEs), not only the simple keratins of LUAD (KRT7/8/18/19). NES is
   larger in LUSC (+3.18) than LUAD (+2.70), as expected for squamous
   tissue, but LUAD is not null. Focal genes: LUAD tracks simple KRTs
   (KRT7/19) more than KRT5; LUSC tracks KRT5. That is histology-shaped
   keratin, not one barrier state.
6. **Purity.** TACSTD2 vs ABSOLUTE is reported above. A near-zero rho
   means the high/low split is not a purity split. It does not make the
   GSEA tumor-cell-intrinsic.
7. **Bulk RNA.** These are mixed-tissue tumors. GSEA cannot say TACSTD2
   *causes* tight junctions or blocks EMT.
8. **Not protein, not ADC, not ICI.** TROP2 protein (the ADC target) was
   not measured. This is not a response analysis.
9. **NES implementation.** This is a documented prerank GSEA, not the
   Broad desktop GUI. Do not compare NES magnitudes to a paper that used
   a different ranking or a different MSigDB freeze without re-running.

## Reproduce

```
pip install -r requirements.txt
python scripts/rework_A8_gsea.py
```

Downloads ~30 MB of public tables into `data/` (gitignored) on first run.
Runtime is dominated by 1000-permutation GSEA × 2 cohorts × 2 splits.

## Files

- `gsea_prerank_all.tsv` — NES / ES / nom p / BH-FDR / leading edge, both splits
- `gsea_prerank_quartile.tsv` — quartile rows only
- `gsea_primary_quartile.tsv` — the 12 primary sets, quartile
- `signature_scores_vs_tacstd2.tsv` — z-mean Spearman
- `focal_gene_correlations.tsv`
- `verdicts.tsv`
- `sample_table.tsv` — TACSTD2, quartile/median labels, ABSOLUTE
- `summary.json` / `provenance.json`
- `figures/nes_primary.png`
- `figures/nes_hallmark_heatmap.png`
- `figures/focal_rho.png`

## Data

- Expression: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` and `TCGA.LUSC.sampleMap/HiSeqV2`.
- Purity: GDC `4f277128-f793-4354-a13d-30cc7fe9f6b5` (PanCanAtlas ABSOLUTE).
- Gene sets: `data/genesets/a8_sets.json` (Enrichr Hallmark 2020 + KEGG 2021 + GO BP 2023 + MSigDB v2023.2.Hs keratin/TJ terms).

