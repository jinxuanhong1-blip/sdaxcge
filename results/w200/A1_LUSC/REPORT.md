# CLAIM A1 (LUSC-only): TACSTD2 vs immune features after purity adjustment

**Cohort:** TCGA-LUSC primary tumors (`-01` samples). 502 with RNA-seq; **n = 493** with RNA-seq + ABSOLUTE purity (analysis set). 478 also have Wolf signature scores, 492 have methylation leukocyte fraction.

**Question tested:** Is TACSTD2 (TROP2) mRNA expression associated with immune features in **TCGA-LUSC only** (not LUAD, not pooled NSCLC) after adjusting for tumor purity?

**Estimator:** partial Spearman ρ of TACSTD2 vs each feature given ABSOLUTE purity (Pearson on rank residuals; t-test with n−3 df). BH-FDR across 32 features.

## Verdict: PARTIAL SUPPORT — RNA T-cell / cytotoxic markers yes; total leukocytes no

TACSTD2-high LUSC has a **modestly T-cell / cytotoxic / checkpoint-low mRNA profile** that survives and slightly strengthens after ABSOLUTE purity adjustment. The strongest purity-adjusted |ρ| is 0.29 (CXCL9). Most T-cell and exhaustion genes sit at partial ρ ≈ −0.20 to −0.25.

That is **not** the same as “immune cold.” The only endpoint that is orthogonal to RNA — DNA-methylation leukocyte fraction — is **null** (partial ρ = +0.012, FDR = 0.82). Wolf IFN-gamma and TGF-beta signatures are also null. Anyone calling TACSTD2-high LUSC an immune desert from these data would be overstating it.

## Key results (partial Spearman | ABSOLUTE purity; BH-FDR across 32 features)

Significant after purity adjustment (FDR < 0.05), strongest first:

| Feature | Unadjusted ρ | Purity-adjusted ρ | Adjusted FDR | n |
|---|---:|---:|---:|---:|
| CXCL9 | −0.267 | **−0.292** | 1.2e-09 | 493 |
| LAG3 | −0.221 | **−0.257** | 1.1e-07 | 493 |
| PDCD1 | −0.195 | **−0.249** | 2.4e-07 | 493 |
| CD8A | −0.205 | **−0.244** | 3.4e-07 | 493 |
| CTLA4 | −0.181 | **−0.240** | 4.4e-07 | 493 |
| TIGIT | −0.173 | **−0.227** | 2.0e-06 | 493 |
| FOXP3 | −0.158 | **−0.218** | 4.6e-06 | 493 |
| CD2 | −0.161 | **−0.217** | 4.6e-06 | 493 |
| CIBERSORT M1 macrophages (rel.) | −0.215 | **−0.215** | 5.5e-06 | 493 |
| NKG7 | −0.168 | **−0.209** | 8.4e-06 | 493 |
| PRF1 | −0.166 | **−0.209** | 8.4e-06 | 493 |
| IFNG | −0.188 | **−0.208** | 8.4e-06 | 493 |
| CD3E | −0.143 | **−0.205** | 1.1e-05 | 493 |
| CXCL10 | −0.176 | **−0.204** | 1.2e-05 | 493 |
| HAVCR2 | −0.139 | **−0.198** | 1.9e-05 | 493 |
| Wolf macrophage CSF1 response | −0.123 | **−0.201** | 1.9e-05 | 478 |
| Wolf lymphocyte infiltration | −0.132 | **−0.201** | 1.9e-05 | 478 |
| GZMB | −0.148 | **−0.185** | 6.7e-05 | 493 |
| Cytolytic activity (mean GZMA/PRF1) | −0.134 | **−0.172** | 2.1e-04 | 493 |
| GZMA | −0.109 | **−0.141** | 2.8e-03 | 493 |
| CIBERSORT CD8 T cells (rel.) | −0.139 | **−0.131** | 5.6e-03 | 493 |
| CIBERSORT activated NK (rel.) | −0.127 | **−0.121** | 1.1e-02 | 493 |
| CD274 (PD-L1) | −0.098 | **−0.110** | 2.1e-02 | 493 |

Notably **not** significant after purity: methylation leukocyte fraction (ρ = +0.012), Wolf IFN-gamma (ρ = −0.066), Wolf TGF-beta (ρ = −0.025), CIBERSORT Tregs (ρ = +0.050), CIBERSORT M2, IDO1, CD68, MS4A1 (FDR 0.078), wound-healing CSR.

Full table: `a1_lusc_correlations.tsv`. Per-sample data: `a1_lusc_analysis_table.tsv`. Figures: `figures/`.

## Honest interpretation

1. **RNA T-cell / cytotoxic / checkpoint genes are consistently negative.** Direction is coherent across CD8A, CD3E, CD2, CYT, GZMA/GZMB/PRF1/NKG7, CXCL9/10, IFNG, and the exhaustion/checkpoint set (PDCD1, CTLA4, LAG3, TIGIT, HAVCR2). Purity adjustment does not create this: every one of these was already negative unadjusted, and partialling makes them slightly more negative.
2. **Effect sizes are modest.** A significant ρ of −0.24 is not “immune desert because of TROP2.” It is a weak-to-moderate rank association in mixed bulk tissue. CXCL9 at −0.29 is the ceiling in this scan.
3. **Total leukocyte content is null.** Methylation leukocyte fraction is the cleanest test here (DNA, not the same RNA library as TACSTD2). Partial ρ = +0.012, FDR = 0.82. The purity-vs-leukocyte check works (ρ = −0.73), so the adjustment is not broken — there is just no TACSTD2–leukocyte association to find.
4. **Purity is not the story, but it is not a no-op either.** TACSTD2 vs ABSOLUTE is weak and non-significant (ρ = −0.074, p = 0.10). Immune mRNA vs purity is strongly negative (as expected). Residual-rank partialling therefore removes a small *positive* confounder and the immune associations become a bit more negative. That is the opposite of “purity invented a cold phenotype.”
5. **Wolf IFN-gamma and TGF-beta are null in LUSC.** Do not import the LUAD-only pattern (there, TGF-beta was the strongest *positive* hit and Wolf IFN-gamma was weakly positive). Those two signatures do not replicate as LUSC facts.
6. **CIBERSORT CD8 is weaker than CD8A mRNA** (partial ρ −0.13 vs −0.24). Relative fractions are compositional (they sum to 1 among leukocytes) and are not absolute CD8 abundance. Use them as a supporting row, not the headline.
7. **This is LUSC only.** These numbers must not be quoted as LUAD or as pooled NSCLC. The matched LUAD-only scan (same 32 features, same estimator) is weaker and more mixed: CD8A was FDR-borderline, leukocyte fraction was also null, and TGF-beta was the top hit. LUSC is where the RNA T-cell / cytotoxic anti-correlation is actually visible.

## Methods

- Expression: UCSC Xena `TCGA.LUSC.sampleMap/HiSeqV2`, log2(RSEM normalized count + 1); primary tumors only.
- Purity: PanCanAtlas ABSOLUTE consensus calls (`TCGA_mastercalls.abs_tables_JSedit.fixed.txt`).
- Immune features (32): methylation leukocyte fraction (Thorsson/PanImmune); 5 CIBERSORT relative fractions (unfiltered by CIBERSORT p-value, matching the LUAD-only A1 script); 5 Wolf signature scores; 20 canonical immune marker genes + cytolytic score (mean log2 GZMA/PRF1).
- Unadjusted: Spearman. Adjusted: partial Spearman (Pearson on rank residuals after regressing out ranked purity), t-test with n−3 df. BH-FDR across all 32 features, computed separately per method. Replicate aliquots averaged per patient. No sampling/randomness; fully reproducible via `scripts/a1_lusc_tacstd2_immune.py`.

## Caveats

- Bulk-RNA observational association; nothing here is causal.
- mRNA only; TROP2 protein (the ADC target) was not measured.
- TACSTD2 and marker genes come from the same RNA-seq libraries (shared normalization / composition). The methylation leukocyte fraction does not, and it is the null result.
- CIBERSORT fractions are *relative to leukocyte content*, not absolute.
- Wolf signatures come from a different pan-cancer expression freeze; barcode intersection reduces n to 478 for those rows.
- Purity adjustment via a single covariate cannot remove all compositional confounding.

## Reproduce

```
pip install -r requirements.txt
python3 scripts/a1_lusc_tacstd2_immune.py
```

Downloads public tables into `data/` (gitignored) on first run. URLs and md5 checksums: `provenance.json`.
