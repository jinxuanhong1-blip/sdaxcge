# Hunt: Tacstd2 (Trop2) in GSE76628 — public data only

**Update:** GSE76628 is not Kras;Lkb1 lung. The A5 recompute that answers
Tacstd2-high vs CD8/NK on a real public KL lung series is
`results/w200/A5_GSE165641/` (GSE165641). GSE76628 rejection:
`results/w200/A5_GSE76628/`.

**Request:** "GSE76628 public only (no private KL mice): Tacstd2-high tumor subset vs CD8/NK."

**Verdict (honest): the requested contrast is infeasible in this dataset.** GSE76628 is not a
tumor single-cell or sorted-population dataset. The closest feasible bulk-level proxies were run
instead and are reported below with their limitations stated plainly. Nothing here supports (or
refutes) a claim that a Tacstd2-high tumor subset is distinguishable from CD8/NK cells.

## What GSE76628 actually is

| Field | Value |
|---|---|
| Title | Stromal-Based Signatures for the Classification of Gastric Cancer [part II] |
| Publication | Uhlik et al., Cancer Res 2016 (PMID [27197264](https://pubmed.ncbi.nlm.nih.gov/27197264/)) |
| Platform | GPL1261, Affymetrix Mouse Genome 430 2.0 — **bulk microarray** |
| Samples | 78 whole flank-skin tissue samples |
| Model | **Ad-VEGF-A164 tumor-surrogate angiogenesis model** (adenoviral VEGF injection into normal flank skin) — there are **no tumors and no tumor cells** |
| Mice | **Athymic nude (Foxn1nu)** females, 4–6 weeks — **no mature T cells exist in these animals** |
| Design | Normal (n=8) + Ad-VEGF at day 5/20/60, each untreated (NT) or treated with DC101 (anti-VEGFR2) or G6-31 (anti-VEGF); n=8 per group (60d DC101 n=6) |

## Why the requested contrast cannot be run

1. **No tumor cells.** The "tumor surrogate" is a VEGF-driven vascular/stromal lesion in normal
   skin. There is no malignant epithelial compartment, hence no "Tacstd2-high tumor subset."
2. **No cell-level resolution.** Bulk arrays of whole tissue; no scRNA-seq, no FACS-sorted CD8 or
   NK populations anywhere in the series or its superseries.
3. **No CD8 T cells in the animals at all.** Nude mice are T-cell deficient. Any bulk `Cd8a`
   signal most plausibly reflects CD8α+ dendritic cells or noise, not CD8 T cells. (NK cells are
   present in nude mice, so NK markers are at least biologically meaningful.)
4. **No KL mice.** The series contains no Kras;Lkb1 animals; the "public only, no private KL
   mice" constraint is satisfied trivially — only public GEO files were used (see
   `provenance.json`).

## What was run instead (closest honest proxies)

Script: `scripts/hunt_gse76628.py`. Processing: GEO series-matrix linear intensities,
log2(x+1); gene level = highest-mean probe per gene; all 78 samples.

### 1. Is Tacstd2 expressed, and does it vary by group?

Tacstd2 has one probe on this array (`1423323_at`). It is **highly expressed in all groups**
(group means 11.0–13.0 log2; 75th–94th within-array percentile) — expected, because Trop2 is
constitutively expressed in normal mouse skin epithelium, which dominates these bulk samples.
The very strong correlation with `Epcam` (Spearman rho = 0.82, p ≈ 6e-20) supports an epithelial
origin of the signal. Group differences are modest but non-random (Kruskal–Wallis p = 2.0e-4),
highest in 60-day G6-treated lesions. See `tacstd2_group_stats.csv` and `tacstd2_by_group.png`.

**This is skin-epithelium Trop2, not tumor Trop2.** It says nothing about tumor-specific
expression.

### 2. Bulk co-abundance of Tacstd2 with CD8/NK marker scores

Across the 78 bulk samples (`tacstd2_vs_immune_correlations.csv`,
`tacstd2_vs_cd8_nk_scores.png`):

| Comparator | Spearman rho | p |
|---|---|---|
| CD8_T score (Cd8a, Cd8b1, Cd3e, Cd3d, Cd3g) | +0.31 | 5.3e-3 |
| NK score (Ncr1, Klrk1, Eomes, Prf1, Gzmb, Nkg7, Il2rb) | -0.37 | 8.9e-4 |
| Epcam (epithelial control) | +0.82 | 5.6e-20 |

Interpretation limits, stated bluntly:

- These are **tissue-level compositional correlations**, not comparisons of expression between
  cell types. A negative Tacstd2–NK correlation here means samples with more NK-like signal have
  proportionally less epithelium (or vice versa) — nothing more.
- The "CD8" score is **not interpretable as CD8 T cells** in T-cell-deficient nude mice.
- `Klrb1c` (NK1.1) has no annotated probe in the GPL1261 annotation used; it is absent from the
  NK score (recorded in `provenance.json`).

## What would be needed to actually answer the question

A dataset with cell-level resolution of tumors and immune cells, e.g. scRNA-seq of mouse lung
tumors (public KL/KP model series exist on GEO) or sorted tumor vs CD8/NK populations. GSE76628
is the wrong substrate for this question regardless of analysis choices.

## Files

| File | Contents |
|---|---|
| `tacstd2_group_stats.csv` | Tacstd2 probe stats per experimental group + Kruskal–Wallis |
| `tacstd2_vs_immune_correlations.csv` | Spearman correlations vs CD8/NK/epithelial markers and panel scores |
| `marker_expression_per_sample.csv` | Per-sample log2 intensity and within-array percentile for all probes of all panel genes |
| `tacstd2_by_group.png` | Tacstd2 by group, boxplots |
| `tacstd2_vs_cd8_nk_scores.png` | Tacstd2 vs CD8/NK panel scores, per-sample scatter |
| `provenance.json` | Source URLs, MD5 checksums, processing notes, missing genes |

## Reproduce

```bash
pip install pandas numpy scipy matplotlib
python3 scripts/hunt_gse76628.py   # downloads public GEO files into data/ if absent
```
