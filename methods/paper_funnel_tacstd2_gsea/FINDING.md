# FINDING — Tacstd2-high vs low malignant enrichment (PAPER FUNNEL)

**Verdict: tight junction / cell adhesion / Claudin family do not rank #1–3 in any cohort.**

Public-only re-run. Numbers below are from `python3 analyze.py` (seed 49901, 2000 gene-set permutations). Nothing fabricated.

## Design

| Cohort | What it is | n | Tacstd2 split |
|---|---|---:|---|
| GSE137244 | GEMM nodule–derived **epithelial cell lines** (KP+KL) | 10 | 5 vs 5 |
| GSE164758 | Public GEMM **primary tumors**, untreated KL+KP | 17 | 8 vs 9 |
| GSE137396 | Public GEMM **lung nodules**, KL+KP | 10 | 5 vs 5 |

- Rank metric: Welch *t* (Tacstd2-high − Tacstd2-low); Tacstd2 removed from the ranked list.
- Ranking universe: 589 sets (Hallmark Mm + KEGG_2019_Mouse + WikiPathways_2019_Mouse + GO TJ/adhesion focus + curated Claudin family).
- Rank for the highlight test = position among **positive NES** sets, sorted by NES descending.

### Genotype confound (not hidden)

| Cohort | Tacstd2-high arms | Tacstd2-low arms |
|---|---|---|
| GSE137244 | 5/5 KL | 5/5 KP |
| GSE164758 | 8 KL + 0 KP | 1 KL + 8 KP |
| GSE137396 | 4 KL + 1 KP | 1 KL + 4 KP |

On GSE137244 the Tacstd2 median split is **identical** to the locked KL-vs-KP contrast.

## Highlight test: do TJ / adhesion / Claudin rank #1–3?

**No. `any_cohort_highlight_in_top3 = false`.**

Best positive-NES term per family (from `tables/highlight_rank_summary.tsv`):

| Cohort | Family | Best term | Rank (pos NES) | NES | FDR q |
|---|---|---|---:|---:|---:|
| GSE137244 | tight junction | GOBP_TIGHT_JUNCTION_ORGANIZATION | **21** | +1.990 | 0.086 |
| GSE137244 | cell adhesion | REACTOME_APOPTOTIC_CLEAVAGE_OF_CELL_ADHESION_PROTEINS | **24** | +1.937 | 0.045 |
| GSE137244 | Claudin family | CURATED__CLAUDIN_FAMILY | **123** | +1.231 | 0.185 |
| GSE164758 | tight junction | GOBP_REGULATION_OF_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY | *(not positive)* | **−1.078** | 0.939 |
| GSE164758 | cell adhesion | KEGG Cell adhesion molecules (CAMs) | **10** | +2.063 | 0.00057 |
| GSE164758 | Claudin family | CURATED__CLAUDIN_FAMILY | **296** | +0.687 | 0.905 |
| GSE137396 | tight junction | GOBP_TIGHT_JUNCTION_ORGANIZATION | **23** | +2.025 | 0.025 |
| GSE137396 | cell adhesion | GOBP_CELL_CELL_ADHESION_MEDIATED_BY_INTEGRIN | **51** | +1.768 | 0.055 |
| GSE137396 | Claudin family | CURATED__CLAUDIN_FAMILY | **222** | +1.100 | 0.325 |

Closest call: **KEGG CAMs ranks #10** in GSE164758 (FDR 5.7×10⁻⁴) — still outside #1–3. Claudin family never approaches the top.

## Forced headline NES/FDR (named sets)

From `tables/forced_headline_NES_FDR.tsv`:

| Cohort | Term | NES | NOM p | FDR q | Rank (pos) |
|---|---|---:|---:|---:|---:|
| GSE137244 | KEGG Tight junction | **+1.793** | 0.015 | 0.090 | 31 |
| GSE137244 | KEGG CAMs | **−1.974** | 0.005 | 0.059 | — |
| GSE137244 | CURATED Claudin family | **+1.231** | 0.185 | 0.185 | 123 |
| GSE137244 | GOBP tight-junction organization | **+1.990** | 0.009 | 0.086 | 21 |
| GSE137244 | Hallmark apical junction | **−1.338** | 0.045 | 0.091 | — |
| GSE164758 | KEGG Tight junction | **−1.148** | 0.113 | 0.396 | — |
| GSE164758 | KEGG CAMs | **+2.063** | 0.000 | 0.00057 | 10 |
| GSE164758 | CURATED Claudin family | **+0.687** | 0.905 | 0.905 | 296 |
| GSE164758 | GOBP tight-junction organization | **−1.079** | 0.281 | 1.000 | — |
| GSE164758 | Hallmark apical junction | **+1.494** | 0.002 | 0.028 | 61 |
| GSE137396 | KEGG Tight junction | **+1.194** | 0.111 | 0.297 | 190 |
| GSE137396 | KEGG CAMs | **−1.262** | 0.088 | 0.499 | — |
| GSE137396 | CURATED Claudin family | **+1.100** | 0.325 | 0.325 | 222 |
| GSE137396 | GOBP tight-junction organization | **+2.025** | 0.000 | 0.025 | 23 |
| GSE137396 | Hallmark apical junction | **−1.244** | 0.086 | 0.290 | — |

## What actually ranks #1–3 (positive NES)

| Cohort | #1 | NES / FDR | #2 | NES / FDR | #3 | NES / FDR |
|---|---|---|---|---|---|---|
| GSE137244 | Wiki Cholesterol biosynthesis | 2.59 / 0.33 | KEGG Proteasome | 2.54 / 1.00 | KEGG Aminoacyl-tRNA biosynthesis | 2.47 / 0.54 |
| GSE164758 | KEGG Graft-versus-host disease | 2.44 / 0 | KEGG Allograft rejection | 2.43 / 0 | KEGG Autoimmune thyroid disease | 2.40 / 0 |
| GSE137396 | Wiki Cholesterol biosynthesis | 2.62 / 0 | Hallmark cholesterol homeostasis | 2.51 / 0 | Wiki cholesterol metabolism | 2.49 / 0 |

Full ranked tables: `tables/*__gsea_positive_ranked.tsv` and `tables/gsea_all_cohorts.tsv`.

## Paper sentence (honest)

> Across GSE137244 epithelial cell lines and two public GEMM bulk cohorts, Tacstd2-high vs low preranked GSEA does **not** place tight-junction, cell-adhesion, or Claudin-family programs in ranks #1–3. The strongest named hits are GOBP tight-junction organization (NES +1.99 / +2.03, ranks 21 / 23) and KEGG CAMs in GSE164758 (NES +2.06, rank 10, FDR 5.7×10⁻⁴). Curated Claudin family stays weak (ranks 123–296). Top pathways are metabolic / immune, not barrier.

## Not claimed

- Tacstd2-high as a TJ/Claudin-rank-#1 program in these public mouse matrices.
- Independence from KL genotype on GSE137244 (complete confound).
- LLC / TISMO as KL substitutes.
