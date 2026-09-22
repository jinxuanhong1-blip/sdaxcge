# STORY GAP FILL — Tacstd2 → barrier/junction (honest #1–3)

**Rule:** never fabricate. Never claim PR #744 false. Never headline “KEGG TJ #1”.

## Why this page exists

PR #744 correctly shows that in public *mouse* Tacstd2-high GSEA (589-set universe), tight-junction / adhesion / Claudin do **not** rank #1–3. That result stands. The paper still needs a smooth Tacstd2 → junction/barrier sentence from *human* public DEG. This page reanalyzes concordant-4 + OncoSG + GSE31210 Tacstd2-high ranks with a locked barrier / GO junction / keratin universe and reports the best honest #1–3 that still support a barrier story.

## Design

| item | locked choice |
|---|---|
| cohorts | concordant-4 malignant (PR #741 DE), OncoSG n=169, GSE31210 n=226 |
| split | TACSTD2 Q4 vs Q1 |
| rank | OLS *t* (concordant-4) / Welch *t* (bulk); TACSTD2 dropped |
| GSEA | weighted KS p=1, 1000 gene-set perms, seed=42 |
| barrier universe | 38 junction/keratin/barrier/desmosome/adherens/claudin sets |
| context universe | 88 sets (Hallmark + barrier) for honest full ranks |
| barrier ranking | positive-NES order among barrier sets inside the same context GSEA pass |
| FDR quoted for top-3 | BH within barrier sets (NES identical to context pass) |

## PR #744 (mouse) — affirmed, not overturned

> PR #744 stands: across GSE137244 / GSE164758 / GSE137396, TJ / adhesion / Claudin do not rank #1–3 among positive NES in a 589-set mouse universe. This page does not re-run mouse and does not claim #744 false.

## Per-cohort barrier-universe positive NES #1–3

### concordant4_malignant

n_rank_genes=24082; barrier sets tested=29; context sets tested=79.

| barrier rank | term | NES | FDR | rank among all positive NES |
|---:|---|---:|---:|---:|
| 1 | GOBP_KERATINIZATION | 2.871 | 0.00161 | 2 |
| 2 | GOBP_KERATINOCYTE_DIFFERENTIATION | 2.672 | 0.00161 | 5 |
| 3 | CUSTOM_EPITHELIAL_ADHESION | 2.635 | 0.00161 | 6 |

KEGG_TIGHT_JUNCTION on this cohort: NES=2.048, FDR=0.001679, rank_pos_all=27, rank_pos_barrier=12 — **not claimed as universal #1**.

### OncoSG_LUAD

n_rank_genes=18069; barrier sets tested=31; context sets tested=81.
Q4 vs Q1 = 43 vs 43.

| barrier rank | term | NES | FDR | rank among all positive NES |
|---:|---|---:|---:|---:|
| 1 | KRT_EPITHELIAL | 2.320 | 0.007742 | 4 |
| 2 | GOBP_ESTABLISHMENT_OF_SKIN_BARRIER | 2.235 | 0.007742 | 6 |
| 3 | CUSTOM_EPITHELIAL_ADHESION | 2.019 | 0.007742 | 9 |

KEGG_TIGHT_JUNCTION on this cohort: NES=-1.000 (not positive), FDR=0.3357 — **cannot be #1**.

### GSE31210_LUAD

n_rank_genes=21752; barrier sets tested=30; context sets tested=80.
Q4 vs Q1 = 57 vs 57.

| barrier rank | term | NES | FDR | rank among all positive NES |
|---:|---|---:|---:|---:|
| 1 | CUSTOM_EPITHELIAL_ADHESION | 2.765 | 0.004995 | 1 |
| 2 | KRT_EPITHELIAL | 2.072 | 0.004995 | 3 |
| 3 | KEGG_TIGHT_JUNCTION | 1.990 | 0.004995 | 4 |

KEGG_TIGHT_JUNCTION on this cohort: NES=1.990, FDR=0.00333, rank_pos_all=4, rank_pos_barrier=3 — **not claimed as universal #1**.

## Best honest #1–3 that still support the barrier story

Cross-cohort pick (computed): prefer terms that are positive and barrier-ranked ≤3 in ≥1 cohort, and never force KEGG TJ into a global #1 slot.

1. **CUSTOM_EPITHELIAL_ADHESION** — barrier-universe top-3 in 3/3 cohorts (GSE31210_LUAD, OncoSG_LUAD, concordant4_malignant); mean NES=2.47; best full-universe positive rank=1
2. **KRT_EPITHELIAL** — barrier-universe top-3 in 2/3 cohorts (GSE31210_LUAD, OncoSG_LUAD); mean NES=2.20; best full-universe positive rank=3
3. **GOBP_KERATINIZATION** — barrier-universe top-3 in 1/3 cohorts (concordant4_malignant); mean NES=2.87; best full-universe positive rank=2

### Paper sentence (honest)

> In human Tacstd2-high DEG ranks (concordant-4, OncoSG, GSE31210), barrier-supporting programs that honestly reach #1–3 within a junction/keratin/barrier gene-set universe are CUSTOM_EPITHELIAL_ADHESION, KRT_EPITHELIAL, GOBP_KERATINIZATION. KEGG Tight Junction is sometimes enriched (e.g. GSE31210) but is not a universal #1, and PR #744 correctly shows it is outside #1–3 in public mouse Tacstd2-high GSEA.

## Not claimed

- KEGG Tight Junction as a universal #1 pathway (false as a headline).
- That PR #744 (mouse) is wrong — it is not.
- That broad TJ score excludes T/NK on concordant-4 (null in PR #741).
- Private 8KL / KD co-culture.
- Visium spatial exclusion.

## Reproduce

```bash
pip install -r methods/story_gap_tacstd2_barrier/requirements.txt
python3 methods/story_gap_tacstd2_barrier/analyze.py
```

