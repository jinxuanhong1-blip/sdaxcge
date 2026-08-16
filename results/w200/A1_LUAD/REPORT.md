# CLAIM A1 (LUAD-only): TACSTD2 vs immune features after purity adjustment

**Cohort:** TCGA-LUAD primary tumors (`-01` samples). 515 with RNA-seq; **n = 502** with RNA-seq + ABSOLUTE purity (analysis set). 444 also have Wolf signature scores, 500 have methylation leukocyte fraction.

**Question tested:** Is TACSTD2 (TROP2) mRNA expression associated with immune features in TCGA-LUAD after adjusting for tumor purity?

## Verdict: WEAK / PARTIAL SUPPORT

TACSTD2 shows statistically detectable but **small** associations with a cytotoxicity-poor, mildly immunosuppressive profile in LUAD. These associations survive purity adjustment — but purity adjustment is nearly a no-op here, because **TACSTD2 is uncorrelated with purity** (Spearman rho = 0.007, p = 0.87). There is **no association with overall immune infiltration** (methylation-based leukocyte fraction: purity-adjusted rho = 0.03, FDR = 0.55). All effect sizes are |rho| ≤ 0.22, i.e. explaining <5% of rank variance. Anyone describing TACSTD2-high LUAD as "immune cold" based on this data would be overstating it.

## Key results (partial Spearman | ABSOLUTE purity; BH-FDR across 32 features)

Significant after purity adjustment (FDR < 0.05):

| Feature | Unadjusted rho | Purity-adjusted rho | Adjusted FDR |
|---|---|---|---|
| Wolf TGF-beta response signature | 0.180 | **0.216** | 1.4e-04 |
| GZMB | -0.143 | **-0.158** | 6.4e-03 |
| Wolf IFN-gamma signature | 0.119 | **0.148** | 1.9e-02 |
| GZMA | -0.124 | **-0.135** | 1.9e-02 |
| LAG3 | -0.113 | **-0.122** | 3.7e-02 |
| CIBERSORT Tregs (relative) | 0.117 | **0.119** | 3.7e-02 |
| CIBERSORT CD8 T cells (relative) | -0.118 | **-0.118** | 3.7e-02 |
| Cytolytic activity (mean GZMA/PRF1) | -0.106 | **-0.116** | 3.7e-02 |

Notably **not** significant: leukocyte fraction (rho 0.03), CD8A (rho -0.104, FDR 0.059), CD274/PD-L1 (rho 0.05), PDCD1, CTLA4 (FDR 0.057), FOXP3, CXCL9/10, IDO1, all macrophage features.

Full table: `a1_luad_correlations.tsv`. Per-sample data: `a1_luad_analysis_table.tsv`. Figures: `figures/`.

## Honest interpretation

1. **Direction is coherent but effect sizes are small.** TACSTD2-high LUAD trends toward lower cytotoxic effector expression (GZMA/GZMB/CYT, CIBERSORT CD8) and higher Tregs and TGF-beta response — a mild immunosuppressive skew — with no change in total immune infiltration. The strongest single association in the whole scan is rho = 0.216 (TGF-beta signature).
2. **Purity adjustment does not rescue or kill anything.** Because TACSTD2 is orthogonal to purity, adjusted and unadjusted estimates are nearly identical. The adjustment machinery works (purity vs leukocyte fraction rho = -0.77, as expected), it just has nothing to remove for this gene.
3. **An apparent inconsistency, reported as-is:** the Wolf IFN-gamma *signature* is weakly **positive** (rho = 0.148) while the IFNG *gene* itself is weakly negative (rho = -0.104, FDR = 0.059). These are different measurements (multi-gene module on pan-cancer-normalized expression vs single gene on the LUAD Xena matrix) and both are small; do not cherry-pick either one.
4. **Multiple features are borderline** (CTLA4, IFNG, CD8A, PRF1 at FDR 0.057–0.06). With a slightly different feature list or FDR threshold the "significant" set changes; the small-effect-size picture does not.

## Methods

- Expression: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2`, log2(RSEM normalized count + 1); primary tumors only.
- Purity: PanCanAtlas ABSOLUTE consensus calls (`TCGA_mastercalls.abs_tables_JSedit.fixed.txt`).
- Immune features (32): methylation leukocyte fraction (Thorsson/PanImmune); 5 CIBERSORT relative fractions; 5 Wolf signature scores; 20 canonical immune marker genes + cytolytic score (mean log2 GZMA/PRF1).
- Unadjusted: Spearman. Adjusted: partial Spearman (Pearson on rank residuals after regressing out ranked purity), t-test with n-3 df. BH-FDR across all 32 features, computed separately per method. Replicate aliquots averaged per patient. No sampling/randomness; fully reproducible via `scripts/a1_luad_tacstd2_immune.py`.

## Caveats

- Bulk-RNA observational association; nothing here is causal.
- mRNA only; TROP2 protein (the drug target of sacituzumab govitecan / datopotamab deruxtecan) was not measured.
- CIBERSORT fractions are *relative to leukocyte content*, not absolute; used unfiltered by CIBERSORT p-value.
- TACSTD2 and marker genes come from the same RNA-seq libraries (shared normalization); the Wolf signatures come from a different pan-cancer expression freeze — barcode intersection reduces n to 444 for those rows.
- Purity adjustment via a single covariate cannot remove all compositional confounding, but given TACSTD2 ⟂ purity, residual purity confounding is unlikely to be the story here.

Data provenance (URLs, md5 checksums): `provenance.json`.
