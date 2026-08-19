# RESULTS — leftover OPEN lung ST, CLDN4 vs CD8A

Additive public spatial. CLDN4-only. Not the CosMx demo. No claim-failed page.

Question, per runnable section: does CLDN4 sit away from CD8A? Four numbers: epithelial-spot Spearman, KRT8 residual Spearman, nearest CD8-high distance from CLDN4-high vs CLDN4-low epithelium, and local CD8-high neighbor counts.

Epithelium = KRT8+ (or KRT8 ≥ median if KRT8 is common). CLDN4-high = top quartile of CLDN4 on epithelium, or CLDN4+ if that quartile is 0. CD8-high = top quartile of CD8A among CD8A+ cells/spots. Neighborhood radius = 3× median nearest-neighbor spacing (Stereo-seq / HD cells) or 2 Visium rings.

Reproduce: `python3 methods/leftover_open_st/analyze.py --only stereo` (and `--only gse301973`). Zenodo Visium: load `filtered_feature_bc_matrix.h5` + `tissue_positions_list.csv` as in `analyze.py`.

---

## Runnable sets

### GSE328481 — Stereo-XCR-seq LUAD (11 sections)

Public H5AD with `x/y` and `obsm['spatial']`. CLDN4 and CD8A present. 1,098,344 QC cells.

| Test (section-level) | n | median | n neg / n pos | Wilcoxon p |
|---|---:|---:|---:|---:|
| same-spot ρ(CLDN4, CD8A) on epithelium | 11 | +0.002 | 4 / 7 | 0.32 |
| KRT8 residual ρ(CLDN4, CD8A) | 11 | +0.765 | 1 / 10 | 0.054 |
| nearest CD8 distance, CLDN4-hi − lo | 11 | +20.8 (pixel) | 3 / 8 | 0.17 |

Same-spot ρ is near zero. Residual ρ is **positive** on most sections (not exclusion after KRT8). Nearest-CD8 is farther from CLDN4-high epithelium on 8/11 sections, but the signed-rank test is NS. CD8A is sparse (0.3–1.7% cells). Neighbor-count medians are 0 on every section.

Per-section table: `tables/gse328481_per_section.csv`.

### GSE301973 — leftover Visium HD EGFR NSCLC (2 slides)

StarDist-integrated cells (not the 11 GB 2 µm bin dump). CLDN4 present. Honest n = **2 slides** (GEO: four patients on two slides; public `sample` key is one label per slide).

| Slide | n cells | n epi | same-spot ρ | KRT8 residual ρ | nearest CD8 hi vs lo |
|---|---:|---:|---|---|---|
| A | 45,661 | 1,523 | +0.051 (p=0.047) | +0.395 | 870 vs 931, MW p=0.66 |
| C | 165,630 | 8,459 | −0.010 (p=0.35) | **−0.648** | 1153 vs 1291, MW p=2.1×10⁻⁵ |

Slide C is the only leftover HD contrast with a negative KRT8 residual and *closer* CD8 to CLDN4-high epithelium (opposite of exclusion). n=2: no section-level Wilcoxon.

### Zenodo 7306132 — LUAD Visium with coordinates (11 sections)

STopover / Na–Choi LUAD deposit. Filtered H5 + `tissue_positions_list.csv`. CLDN4 well detected (27–99% spots).

| Test (section-level) | n | median | n neg / n pos | Wilcoxon p |
|---|---:|---:|---:|---:|
| same-spot ρ(CLDN4, CD8A) on epithelium | 11 | +0.005 | 5 / 6 | 0.90 |
| KRT8 residual ρ | 11 | +0.009 | 3 / 8 | 0.37 |
| nearest CD8 distance, CLDN4-hi − lo | 11 | +0.23 hex | 1 / 6 (4 ties) | 0.20 |

This is the clean leftover Visium (coordinates + CLDN4). Same-spot and residual are both **null** at the section level. Strongest single-section negatives: spa10ca01 ρ=−0.073 (p=0.012); spa17ca01 ρ=−0.088 (p=0.0027). Neighbor counts: spa18ca01/02 have median 0 CD8 neighbors on CLDN4-high vs 1 on CLDN4-low (MW p&lt;10⁻⁶) — local, not a cohort effect.

---

## Skipped (real checks)

- **GSE189357** is scRNA (barcodes/features/matrix only). Spatial member is GSE189487 (already used elsewhere).
- **10x Visium HD lung cancer demo**: CLDN4 almost absent (4 UMIs on the 2 µm feature slice). Skip.
- **Figshare 32384337**: no spot matrix / no `tissue_positions`.
- **Tavernari Zenodo 14624390**: restricted.
- **CNP0005129**: controlled CNGBdb apply.
- **GSE200916**: Visium barcodes, no coordinates deposited.

---

## What this does not claim

- It does not claim a leftover meta-ρ or a uniform CLDN4–CD8 exclusion.
- It does not re-run CosMx, E-MTAB-13530, GSE189487, or prior leftover Visium/GeoMx.
- Stereo-seq residual ρ is inflated by zero-inflation; the Visium STopover residual (~0) is the honest leftover Visium number.
- GSE301973 n=2 is too small to call.

Files: `tables/gse328481_per_section.csv`, `tables/gse301973_per_section.csv`, `tables/zenodo_7306132_per_section.csv`, `tables/cross_section_summary.json`, `inventory.md`.
