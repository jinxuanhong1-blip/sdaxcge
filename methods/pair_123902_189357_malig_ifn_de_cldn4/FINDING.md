# Pair GSE123902+GSE189357: CLDN4-only malignant IFN/MHC/TJ DE

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. **Not CellChat.**
Tumor-cell-intrinsic program only (marker-malignant / epithelial cells).

Given cut (PR #459, **not re-audited**): GSE123902 + GSE189357 malignant
CLDN4 %pos vs same-unit T/NK, n=22, Spearman ρ=−0.638 (p=0.003).
This extra does **not** re-audit that T/NK ρ.

**Thesis (already correct):** CLDN4 KD / low raises the malignant cell's own
IFN / MHC-I; CLDN4-high should be IFN/MHC down, TJ up.

**Numbers not yet written.** Run `build_malignant_pseudobulk.py` then
`analyze.py` to replace this file with the family DE table.

Q4 vs Q1 tails are **7/5 — thin.** GSE189357 Q4 n=2; single-cohort binary
DE on that cohort is skipped. Patient/donor is the unit.
