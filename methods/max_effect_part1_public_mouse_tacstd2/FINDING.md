# MAX EFFECT Part1: public mouse + GSE137244 Tacstd2

Public data only. Private 8 KL single-cell matrices were not read and were not merged. Every result row is labeled **bulk** or **scRNA**. No values were imputed. LLC is not called KL.

Two maximizations on fixed grids:

1. **Bulk** — Tacstd2 KL > KP Δ and Welch |t|
2. **scRNA** — Tacstd2-high → TJ module concordance (mice up / mice scored)

## 1. Bulk: Tacstd2 KL > KP Δ / |t|

Assay: **bulk** (cell line, primary tumor, nodule, LCM, or array). Endpoint: Tacstd2. Δ = mean(KL) − mean(KP). Primary test stays two-sided exact Mann–Whitney; Welch t is reported beside it. Filters: none, minus Epcam, residual on Epcam, residual on Epcam+Krt8, minus Krt8. Leave-one sample rows are stored with reduced n and are not primary.

### Locked reference (prespecified)

**GSE137244** cultured cells, filter `none`, n = 5 vs 5, log2(FPKM+1):

| metric | value |
|---|---:|
| Tacstd2 Δ | **+3.238** |
| Welch t | **+4.648** |
| exact MW p | **0.00794** |

Complete separation. This is the locked handoff number.

### Unadjusted Tacstd2 by accession (fair full-n)

| accession | setting | assay | n KL vs KP | Δ | \|t\| | MW p | MW floor |
|---|---|---|---:|---:|---:|---:|---:|
| GSE244452 | syngeneic_bulk_tumor | bulk | 3 vs 3 | +7.228 | 15.097 | 0.100 | 0.100 |
| GSE137244 | cultured_cells | bulk | 5 vs 5 | +3.238 | 4.648 | 0.00794 | 0.00794 |
| GSE6135 | mouse_all_histology | bulk | 7 vs 5 | +1.795 | 4.212 | 0.00253 | 0.00253 |
| GSE6135 | mouse_adeno_plus_mixed | bulk | 5 vs 5 | +1.337 | 3.255 | 0.00794 | 0.00794 |
| GSE137396 | in_vivo_nodule | bulk | 5 vs 5 | +1.144 | 2.191 | 0.095 | 0.00794 |
| GSE6135 | mouse_adeno_only | bulk | 4 vs 5 | +0.989 | 3.839 | 0.0159 | 0.0159 |
| GSE164758 | primary_bulk_tumor | bulk | 9 vs 8 | +0.858 | 4.710 | 8.2×10⁻⁵ | 8.2×10⁻⁵ |
| GSE274351 | lcm_early_adenoma | bulk | 5 vs 5 | −0.860 | 0.665 | 0.548 | 0.00794 |
| GSE274352 | cultured_cells_empty_vector | bulk | 3 vs 3 | −1.543 | 4.837 | 0.100 | 0.100 |

### Grid maxima (honest strata)

| selection | accession | filter | n | Δ | \|t\| | MW p | note |
|---|---|---|---:|---:|---:|---:|---|
| max Δ / \|t\| among all fair full-n | GSE244452 | none | 3 vs 3 | +7.228 | 15.097 | 0.100 | **n=3 floor**; cannot reach p<0.05 |
| max Δ among MW-floor < 0.05 | GSE137244 | minus_Krt8 | 5 vs 5 | **+3.478** | **4.882** | 0.00794 | grid extreme on these labels |
| max \|t\| among MW-floor < 0.05 | GSE137244 | minus_Krt8 | 5 vs 5 | +3.478 | 4.882 | 0.00794 | same row |
| locked unadjusted | GSE137244 | none | 5 vs 5 | +3.238 | 4.648 | 0.00794 | prespecified |
| smallest MW p (unadjusted) | GSE164758 | none | 9 vs 8 | +0.858 | 4.710 | 8.2×10⁻⁵ | complete separation |

Epcam adjustment on GSE137244 **shrinks** Tacstd2 (minus Epcam Δ=+0.610; residual Δ=+0.868). Minus Krt8 raises both Δ and |t| relative to locked none; that row is a search result, not a replacement for the locked unadjusted contrast.

Leave-one on GSE137244: dropping KP B6AL10-3 gives |t| = 10.518 (Δ = +2.584) at **5 vs 4**. That raises |t|, lowers Δ, and is not 5 vs 5.

Null / reverse series (GSE274351, GSE274352) stay in the table.

Tables: `tables/bulk_kl_kp/`. Figure: `figures/bulk_tacstd2_delta_vs_abs_t.png`.

## 2. scRNA: Tacstd2-high → TJ concordance (mice up/down)

Assay: **scRNA**. Public integrate mice only: GSE154977, GSE180963, GSE165641. Epithelium is gated without Tacstd2. Within each mouse, Tacstd2-high vs Tacstd2-low epithelial cells; frozen TJ module mean log1p score; one-sided MW that high > low. Concordance = mice with Δ > 0 / mice scored. Eligible specs need ≥4 mice.

Grid: gate ∈ {locked, epcam_only, broad_lung} × Tacstd2 rule ∈ {count>0, ≥2, ≥5, median of positives} × floor ∈ {20/20, 30/30, 50/50} × module ∈ {TJ_TISMO, TJ_EPITHELIAL, CLDN4_TJ_EDGE, TJ_CLAUDIN4}.

### Locked baseline (prespecified)

Locked gate, Tacstd2 count > 0, floor 20/20 (and 30/30 identical on this pool), **TJ_TISMO**:

| metric | value |
|---|---:|
| mice scored | 7 |
| mice up (Δ>0) | **7/7** |
| mice up and p<0.05 | **7/7** |
| mean Δ | +0.261 |
| one-sided binomial p | 0.007812 |

GSE180963_K stays below the high/low floor and is inventoried, not scored. Same 7 mice as the prior paper-funnel DEG→TJ arm.

### Concordance is already at the ceiling

Every eligible grid row (144/144) has **frac_up = 1.0** and **frac_up_sig = 1.0**. The mice-up / mice-scored ratio cannot be raised further on this public pool. Maximization therefore moves to **mean module Δ among perfect concordance at the largest n (7)**.

| selection | module | gate | Tacstd2 rule | floor | n | up | mean Δ |
|---|---|---|---|---:|---:|---:|---:|
| max mean Δ among perfect, n=7 | TJ_CLAUDIN4 | broad_lung | count ≥ 5 | 20/20 | 7 | 7/7 | **+0.620** |
| max mean Δ, TJ_TISMO only, n=7 | TJ_TISMO | broad_lung | count ≥ 5 | 20/20 | 7 | 7/7 | **+0.523** |
| locked TJ_TISMO | TJ_TISMO | locked | count > 0 | 20/20 | 7 | 7/7 | +0.261 |

The broad_lung / count≥5 rows are grid extremes. They do not replace the locked count>0 TJ_TISMO baseline. Cell-level p-values inside a mouse are descriptive; the mouse is the unit for the binomial.

Tables: `tables/scrna_tj/`. Figures: `figures/scrna_tj_concordance.png`, `figures/scrna_winner_mouse_deltas.png`.

## What this does not change

- Locked bulk statement: GSE137244 Tacstd2 Δ = +3.24, n = 5 vs 5, MW p = 0.00794.
- Locked scRNA statement: TJ_TISMO 7/7 Tacstd2-high > low on locked-gate epithelium.
- No private 8KL. No bulk↔scRNA merge. No TISMO LLC as KL.

Private 8 KL mice used: **0**.
