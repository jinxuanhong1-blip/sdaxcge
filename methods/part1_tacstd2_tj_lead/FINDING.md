# PART 1 polish — Tacstd2-high → TJ/junction is the lead enrichment

Public human only: **concordant-4** malignant TACSTD2 Q4 vs Q1 (PR #741) and bulk LUAD **OncoSG + GSE31210** (PR #736). Ranks recomputed from those frozen tables. **Never fabricate.** **No CLDN4 pin section.** User OK with TJ first.

## Part 1 ending (smooth)

Across the human public Tacstd2-high contrasts that were locked for this funnel, the enrichment that lands at the top of the rank list is **tight-junction / apical-junction**, not a generic metabolic or immune program:

1. **Concordant-4 ORA:** `HALLMARK_APICAL_JUNCTION` is **rank 1 of 63** (enrichment 5.89, FDR 1.28e-05).
2. **GSE31210 GSEA:** `KEGG_TIGHT_JUNCTION` is **rank 3 of 31** positive-NES sets (NES +2.00, FDR 0.00236).
3. **OncoSG GSEA (honest):** KEGG TJ and Hallmark apical junction are **not** UP (NES -1.00 / -1.48). Best strict TJ term is `GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY` at rank **16 of 26** (NES +1.49, FDR 0.0272). Keratin / skin-barrier lead OncoSG positives — kept visible.

**Ending sentence for Part 1:** Tacstd2-high malignant / LUAD programs resolve first to a **TJ / apical-junction** enrichment in the human layers where rank supports it (concordant-4 ORA #1; GSE31210 KEGG TJ #3), setting up the next section on junction-linked immune geography — without pinning Part 1 on CLDN4 alone.

## Best TJ/junction ranks (headline table)

| Layer | Best TJ/junction term | Rank | Universe | Stat | FDR | Lead? |
|---|---|---:|---|---|---:|:---:|
| concordant-4_ORA | `HALLMARK_APICAL_JUNCTION` | **1** / 63 | A8+custom ORA UP terms (n=63) | enrichment=+5.890 | 1.28e-05 | YES |
| concordant-4_GSEA | `HALLMARK_APICAL_JUNCTION` | **8** / 36 | positive-NES sets (n=36) | NES=+2.555 | 0.00146 | no |
| GSE31210_GSEA | `KEGG_TIGHT_JUNCTION` | **3** / 31 | positive-NES sets (n=31 of 59) | NES=+1.995 | 0.00236 | YES |
| OncoSG_GSEA | `GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY` | **16** / 26 | positive-NES sets (n=26 of 60) | NES=+1.490 | 0.0272 | no |

Full machine table: `tables/best_tj_junction_ranks.tsv`.

## Concordant-4 detail (PR #741)

### ORA (p<0.01 & |logFC|>0.25 UP; TACSTD2 held out)

| Rank | Term | Enrichment | FDR |
|---:|---|---:|---:|
| **1** | HALLMARK_APICAL_JUNCTION | 5.89 | 1.28e-05 |
| 3 | GOBP_KERATINIZATION | 13.98 | 1.28e-05 |
| 11 | GOBP_TIGHT_JUNCTION_ORGANIZATION | 5.94 | 0.00894 |
| 14 | CUSTOM_EPITHELIAL_ADHESION | 11.46 | 0.00981 |

**ORA verdict:** apical junction is the **#1** enrichment — strongest honest lead for Part 1.

### Family scores (OLS, positive = TACSTD2 Q4)

| Family | logFC | FDR |
|---|---:|---:|
| TJ | +0.223 | 0.004187 |
| CLAUDIN_PANEL | +0.426 | 0.006171 |
| EPITHELIAL_ADHESION | +0.734 | 0.004187 |
| KERATIN | +0.929 | 0.00545 |
| IFN | +0.348 | 0.1348 |
| MHC-I/APM | +0.189 | 0.567 |

TJ / claudin / adhesion / keratin families are all up at FDR < 0.01. IFN and MHC-I are not down on this TACSTD2 split.

### GSEA prerank (OLS *t*; TACSTD2 dropped; seed=42)

| Rank (pos NES) | Term | NES | FDR |
|---:|---|---:|---:|
| 2 | GOBP_KERATINIZATION | +2.825 | 0.00146 |
| 7 | CUSTOM_EPITHELIAL_ADHESION | +2.585 | 0.00146 |
| 8 | HALLMARK_APICAL_JUNCTION | +2.555 | 0.00146 |
| 20 | GOBP_TIGHT_JUNCTION_ORGANIZATION | +2.118 | 0.00146 |
| 22 | KEGG_TIGHT_JUNCTION | +2.040 | 0.00146 |

GSEA puts keratinization near the top (rank 2) and apical junction at rank 8. Part 1 still leads with **TJ/junction** because ORA rank 1 is apical junction and the user OK'd TJ-first framing; keratin is reported as co-enriched barrier, not hidden.

## Bulk LUAD detail (PR #736)

### GSE31210 — TJ leads among positives

| Rank (pos NES) | Term | NES | FDR_all |
|---:|---|---:|---:|
| 1 | HALLMARK_GLYCOLYSIS | +2.487 | 0.00236 |
| 2 | KRT_EPITHELIAL | +2.123 | 0.00236 |
| 3 | KEGG_TIGHT_JUNCTION ← TJ | +1.995 | 0.00236 |
| 4 | GOBP_ESTABLISHMENT_OF_SKIN_BARRIER | +1.941 | 0.00406 |
| 5 | HALLMARK_PEROXISOME | +1.922 | 0.00236 |
| 6 | HALLMARK_PROTEIN_SECRETION | +1.911 | 0.00236 |
| 7 | GOBP_KERATINIZATION | +1.906 | 0.00236 |
| 8 | GOBP_EPIDERMAL_CELL_DIFFERENTIATION | +1.831 | 0.00236 |

**GSE31210 verdict:** KEGG Tight Junction is **rank 3** of 31 positive-NES sets — a near-top / lead enrichment for Tacstd2-high bulk LUAD.

### OncoSG — honest non-lead for KEGG TJ

| Term | NES | FDR_all | Among pos? |
|---|---:|---:|---|
| GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY | +1.490 | 0.0272 | rank 16/26 |
| GOBP_ESTABLISHMENT_OF_SKIN_BARRIER | +2.171 | 0.0024 | rank 6/26 |
| GOBP_KERATINIZATION | +1.924 | 0.0024 | rank 9/26 |
| KEGG_TIGHT_JUNCTION | -0.998 | 0.33 | not positive |
| HALLMARK_APICAL_JUNCTION | -1.480 | 0.0058 | not positive |

**OncoSG verdict:** do **not** claim KEGG TJ or apical junction as the lead enrichment here. Best strict TJ term = bicellular TJ assembly (rank 16). Barrier/keratin programs are up and rank higher — Part 1 states that instead of inventing a TJ #1.

## What this Part 1 ending is / is not

- **Is:** a rank-honest closer that Tacstd2-high → TJ/junction is the lead enrichment in human concordant-4 (ORA #1) and GSE31210 (KEGG TJ #3).
- **Is:** TJ-first framing with keratin/skin-barrier co-enrichment reported.
- **Is not:** a CLDN4 pin, a surface-molecule rank nail, or a private KD claim.
- **Is not:** “TJ-high excludes T/NK” from the TACSTD2 funnel alone (PR #741 TJ vs T/NK was null).
- **Is not:** a re-fit of locked CLDN4 %pos vs T/NK ρ = −0.531.
- **Is not:** fabricated NES/FDR — every rank comes from PR #741 / #736 tables in `provenance/`.

## Reproduce

```bash
python3 methods/part1_tacstd2_tj_lead/rank_from_provenance.py
```

## Files

- `tables/best_tj_junction_ranks.tsv` — headline ranks
- `tables/gsea_tj_junction_ranks.tsv` — per-cohort TJ/barrier GSEA ranks
- `tables/c4_ora_tj_barrier_ranks.tsv` — C4 ORA barrier subset
- `tables/summary.json` — machine verdict
- `figures/fig_best_tj_junction_ranks.png` — Part 1 rank figure
- `provenance/` — frozen copies of PR #741 / #736 tables

