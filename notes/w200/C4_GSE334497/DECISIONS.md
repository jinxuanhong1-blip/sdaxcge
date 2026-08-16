# Decisions — C4 analog GSE334497

- Contrast is Trop2 KO minus WT on author-normalized counts. Counts are already size-normalized; do not run DESeq2 size factors again.
- Pre-specified focal genes are Cldn4 and Cxcl9. CORE6 / ISG / APM / T_CYT / OXPHOS / MYC are pre-specified sets, not data-driven.
- Exact permutation (252 splits) is the primary set-level p-value. Competitive MWU is reported but not used to claim IFN opening, because a +0.1 median vs a slightly negative background is not a sample-level program.
- Genome-wide BH-FDR is computed and reported as underpowered (0 genes at FDR 0.05). Calls use pre-specified tests.
- Residual Tacstd2 in KO is not treated as KO failure (bulk tumor).
- This matrix is the Wu et al. depositing series. It cannot independently validate that paper; it can only test whether Claim C4’s IFN/APM opening is visible here.
