# PPT paste — TROP2+ barrier face (Part1)

**One line:** TACSTD2-high malignant cells are **+26 percentage points** more often positive for barrier ligands (F11R / NECTIN2 / CDH1 / LGALS9) toward T/NK than TACSTD2-low cells from the same patient.

| | |
|---|---|
| Gate | TACSTD2 (human TROP2), malignant Q outer 10% |
| Score | CellPhoneDB expression proportion (pp) |
| Per-ligand mean | **+26.14 pp** |
| Family sum (4 ligands) | **+104.5 pp** (not a single proportion) |
| n | 59 patients (11+18+21+9) |
| Positive patients | 95% |
| Cohorts | 4/4 (sums +88 / +100 / +122 / +95) |
| Wilcoxon p | 5.4×10⁻¹¹ |
| Sign-flip p | 1.0×10⁻⁴ |
| I² | 0% |

**Ligands (all 4/4 cohorts +):** F11R +27.2 · NECTIN2 +30.5 · CDH1 +33.5 · LGALS9 +13.4

**Q4 vs Q1 (locks to PR 716 crude):** +25.19 pp, n=63, p=7.2×10⁻¹²

**LIANA / CellChat companions (same cells, smaller units):**

| method | per-ligand Δ | family sum | note |
|---|---:|---:|---|
| CellPhoneDB `lr_means` | +0.126 | +0.505 | LIANA CPDB |
| Connectome `expr_prod` | +0.123 | +0.492 | |
| CellChat Hill (no pop.) | +0.074 | +0.297 | |
| CellChat Hill × pop. | +0.0040 | +0.016 | why official CellChat looks tiny |
| LIANA log2FC | +0.364 | +1.46 | |

**Do not say on this slide**

- Mediation / “through CLDN4” (partial shrink is PR 716; CosMx short-range TACSTD2 cold stays after CLDN4).
- Immune-exclusion ρ = −0.53 (that is **CLDN4 %pos**, not TACSTD2).
- Population-scaled CellChat +0.0037 (CLDN4 PR 616) as if it were this TACSTD2 run.
- Dual-high gate. Extra cohorts (GSE148071 etc.). Spatial exclusion from this LR table.

**Compare:** CLDN4 max-effect (PR 713) = **+27.8 pp** on the same score. TACSTD2 is the parallel Part1 face, not a replacement.
