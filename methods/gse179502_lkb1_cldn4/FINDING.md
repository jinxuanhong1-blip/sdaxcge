# Finding — GSE179502 Lkb1 restore, Cldn4-only (mouse)

ADDITIVE. **Cldn4 only.** No dual-high Tacstd2×Cldn4. No human. No GSE179501
(unsorted total cells). Public Cell Ranger aggr only:
`GSE179502_XTR_sorted_scRNAseq_{matrix.mtx,features.tsv,barcodes.tsv}.gz`
(Murray / Winslow *Nat Commun* 2022, PMID 35228570).

FACS-sorted neoplastic epithelium from KT;Lkb1<sup>XTR</sup> KRAS lung tumors.
Author GEO cohort is **Restored** (restorable + tamoxifen) vs **NonRestored**
(non-restorable ± vehicle/tamoxifen, or restorable + vehicle). Honest unit =
**mouse (n=3 vs 3)**. Cell-level tests are descriptive (barcodes are not
independent). MWU at 3 vs 3 cannot go below **p=0.1**.

Thesis already taken as given. This page only asks: when Lkb1 is restored,
does **Cldn4 go down** and/or **IFN / AT2 go up** in neoplastic epithelium?
And in Cldn4-high vs low cells: IFN/MHC down, TJ/barrier up?

## Decision on this matrix

**Cldn4 goes down and AT2 goes up when Lkb1 is restored.** IFN / MHC-I do
**not** rise with restore (flat on n=3 vs 3).

Cldn4-high cells are **TJ/barrier-high and AT2-low**. They are **not**
IFN/MHC-cold: IFN/MHC is flat-to-up, and the IFN Spearman vs Cldn4 collapses
to null after UMI.

## Honest n

| Unit | Count |
|---|---:|
| GEO libraries / mice | **6** |
| Restored (author cohort) | **3** (CM0879, CM0884, ZR1969) |
| NonRestored (author cohort) | **3** (CM0875, ZR1932, ZR1966) |
| Barcodes on the aggr mtx | 16,017 |
| QC cells (UMI≥500, genes≥200, mito≤20%) | **12,643** |
| Smallest two-sided MWU p at 3 vs 3 | **0.1** |

Do not quote n=16,017 cells as the restore n. Do not quote GSE179500’s 17/11
bulk mice as this scRNA n.

## Restore contrast (mouse means; Restored − NonRestored)

| feature | NonRestored mean | Restored mean | Δ | Welch p | MWU p | rank-biserial |
|---|---:|---:|---:|---:|---:|---:|
| **Cldn4** | 0.300 | 0.143 | **−0.157** | **0.036** | **0.1** | **−1.0** |
| frac Cldn4+ | 0.283 | 0.139 | −0.144 | 0.068 | 0.1 | −1.0 |
| **AT2** (14/14) | 2.155 | 2.665 | **+0.510** | **0.022** | **0.1** | **+1.0** |
| IFN A8 (220 genes) | 0.243 | 0.245 | +0.002 | 0.82 | 1.0 | +0.11 |
| IFN CORE (10/10) | 0.048 | 0.049 | +0.001 | 0.76 | 1.0 | +0.11 |
| MHC-I/APM A8 (20) | 0.588 | 0.583 | −0.005 | 0.94 | 1.0 | −0.11 |
| MHC CORE (9/9) | 0.666 | 0.710 | +0.044 | 0.63 | 0.7 | +0.33 |
| TJ (Cldn4 held out) | 0.388 | 0.344 | −0.043 | 0.003 | 0.1 | −1.0 |
| **Stk11** (restore check) | 0.007 | 0.305 | **+0.298** | **2.0×10⁻⁵** | **0.1** | **+1.0** |
| Tacstd2 (audit only) | 0.834 | 0.628 | −0.207 | 0.067 | 0.1 | −1.0 |

Cldn4 and AT2 have **complete rank separation** across the six mice (every
Restored mouse is below every NonRestored mouse for Cldn4; the reverse for
AT2 and Stk11). IFN/MHC do not separate.

Cell-level Cldn4 Δ is the same sign (−0.164) but that p-value is not a
biological n.

## Per-mouse scores (QC cells; log1p CP10k means)

| mouse | cohort | treatment | restorable | n_QC | Cldn4+ | Cldn4 | Stk11 | AT2 | IFN A8 | MHC A8 | TJ |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CM0875 | NonRestored | Tamoxifen | no | 2522 | 0.249 | 0.261 | 0.005 | 2.007 | 0.237 | 0.537 | 0.392 |
| ZR1932 | NonRestored | Vehicle | no | 972 | 0.302 | 0.340 | 0.013 | 2.240 | 0.257 | 0.689 | 0.388 |
| ZR1966 | NonRestored | Vehicle | yes | 2889 | 0.298 | 0.300 | 0.004 | 2.219 | 0.234 | 0.539 | 0.382 |
| CM0879 | Restored | Tamoxifen | yes | 3322 | 0.114 | 0.134 | 0.294 | 2.455 | 0.234 | 0.530 | 0.342 |
| CM0884 | Restored | Tamoxifen | yes | 836 | 0.225 | 0.215 | 0.311 | 2.731 | 0.254 | 0.639 | 0.353 |
| ZR1969 | Restored | Tamoxifen | yes | 2102 | 0.078 | 0.080 | 0.309 | 2.810 | 0.247 | 0.581 | 0.337 |

ZR1966 is restorable but vehicle-treated — author-labeled **NonRestored**,
and Stk11 stays on the floor. That is the correct Lkb1-off control, not a
fourth Restored mouse.

## Cldn4-high vs low (IFN / MHC / TJ / AT2)

Cldn4 is sparse (8–30% cells >0). Q4 vs Q1 uses rank tails (Q1 is the zero
mass). Also Cldn4+ vs Cldn4−. Positive Δ = higher in Cldn4-high.

### All QC cells (descriptive)

| family | gate | n_low / n_high | Δ median | rank-biserial | Spearman ρ | ρ after UMI |
|---|---|---|---:|---:|---:|---:|
| IFN A8 | Q4 vs Q1 | 3160 / 3160 | **+0.014** | +0.21 | +0.15 | **−0.014** (p=0.11) |
| MHC-I A8 | Q4 vs Q1 | 3160 / 3160 | +0.065 | +0.23 | +0.18 | +0.044 |
| **TJ A8 (Cldn4 out)** | Q4 vs Q1 | 3160 / 3160 | **+0.080** | +0.44 | **+0.50** | **+0.35** |
| **AT2** | Q4 vs Q1 | 3160 / 3160 | **−0.935** | −0.26 | **−0.42** | **−0.46** |
| IFN CORE | Q4 vs Q1 | 3160 / 3160 | +0.043 | +0.21 | +0.15 | +0.10 |
| MHC CORE | Q4 vs Q1 | 3160 / 3160 | +0.048 | +0.13 | +0.01 | −0.097 |
| TJ CORE | Q4 vs Q1 | 3160 / 3160 | +0.187 | +0.42 | +0.47 | +0.33 |

Cldn4+ vs Cldn4− is the same sign (TJ +0.13, AT2 −1.43, IFN +0.013).

### Within-mouse Q4 vs Q1, then Wilcoxon across n=6 mice

| family | median Δ (Q4−Q1) | mice Δ<0 / Δ>0 | Wilcoxon p |
|---|---:|---|---:|
| IFN A8 | +0.008 | 0 / **6** | 0.031 |
| MHC-I A8 | +0.030 | 0 / **6** | 0.031 |
| **TJ A8** | **+0.093** | 0 / **6** | 0.031 |
| **AT2** | **−0.976** | **6** / 0 | 0.031 |
| MHC CORE | −0.006 | 3 / 3 | 0.69 |

Wilcoxon two-sided floor at 6/6 same sign is **p=0.031**. IFN/MHC are
slightly **up** in Cldn4-high cells in every mouse — the opposite of
“IFN/MHC down.” The IFN effect is tiny and UMI-sensitive. TJ up and AT2
down are large and stable.

Mouse-mean Spearman (n=6): Cldn4 vs TJ ρ=+0.83 (p=0.042); vs AT2 ρ=−0.71
(p=0.11); vs IFN ρ=+0.14 (p=0.79); vs MHC-I ρ=+0.37 (p=0.47).

## What this does / does not say

**Holds**

- Lkb1 restore (Stk11 back on) lowers neoplastic **Cldn4** and raises the
  **AT2** program. Same direction as the paper’s ATII shift and as the
  sibling bulk series GSE179500 (author Cldn4 log2FC −0.77, not re-scored
  here).
- Cldn4-high neoplastic cells are a **TJ/barrier-high, AT2-low** state.

**Does not hold**

- Restore does **not** open IFN / MHC-I on these six mice.
- Cldn4-high cells are **not** IFN/MHC-low on this mtx.
- Tacstd2 is audit-only (same restore sign as Cldn4). It is not a gate.

## Methods (locked)

- Public GEO mtx + features + barcodes only. Barcode prefix = mouse
  (CM0875, CM0879, CM0884, ZR1932, ZR1966, ZR1969). Cohort / treatment /
  genotype from the GSE179502 SOFT record.
- QC: UMI ≥500, genes ≥200, mito (`mt-*`) ≤20% on the cellranger-filtered
  aggr. Score = mean log1p(CP10k) of genes present on mm10 features.
- AT2 = paper Fig. 4c + standard ATII markers (14/14 present). IFN/MHC/TJ
  = leftover CORE lists plus A8 Hallmark / MHC-I / KEGG+GOBP TJ mapped to
  mouse (Cldn4 and Stk11 held out of TJ).
- Restore tests: Welch *t* and two-sided MWU on **mouse means**. Cell-level
  MWU is written but not the n.
- Cldn4-high vs low: rank-tail Q4 vs Q1 (zeros stay in Q1), Cldn4+ vs −,
  within-mouse deltas + Wilcoxon, and Spearman ± ranking out log UMI.
- Cldn4 only. No dual-high. No human. No GSE179501. No CellChat.
- Reproduce: `python3 methods/gse179502_lkb1_cldn4/analyze.py`

## Files

- `analyze.py` / `gene_sets.py` — download-once GEO cache + scores
- `tables/samples.tsv` — six mice, GEO labels, honest n
- `tables/mouse_scores.tsv` — per-mouse Cldn4 / Stk11 / AT2 / IFN / MHC / TJ
- `tables/restore_contrasts.tsv` — Cldn4 by Lkb1-restore condition
- `tables/cldn4_highlow_ifn.tsv` — Cldn4-high vs low IFN / MHC / TJ / AT2
- `tables/cldn4_highlow.tsv` / `cldn4_highlow_by_mouse.tsv`
- `tables/gene_coverage.tsv` / `one_row.tsv` / `summary.json`
- `figures/cldn4_by_restore.{png,pdf}` — Cldn4 and IFN/AT2/MHC by cohort
- `figures/cldn4_highlow_ifn_tj.{png,pdf}` — Q4−Q1 deltas

## How to read this

- n=3 vs 3 mice. Complete Cldn4 / AT2 / Stk11 separation; MWU still 0.1.
- This is sorted neoplastic epithelium. There is no T/NK fraction here.
- Do not fold this into a human concordant pool. Do not call it dual-high.
