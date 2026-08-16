# A11 methods (bulk)

Full decision rules: `00_analysis_plan.md`. This is the implementation, not a second plan.

## Samples

- TCGA-LUAD / TCGA-LUSC STAR `log2(TPM+1)` from the UCSC Xena GDC hub.
- Restricted to primary tumour barcodes (`XXXX-01`). Adjacent normal (`-11`) used only for the tumour-vs-normal context table.
- ABSOLUTE purity/ploidy from the GDC PanCanAtlas table, joined on the 15-character sample barcode. Samples without a call are kept for naive statistics and dropped (pairwise) for anything that uses purity.

## Genes

HGNC symbols via GENCODE v36 probemap. When two Ensembl IDs share a symbol, the higher-mean locus is kept. The primary panel is the 30 genes on the CD47, Galectin, Nectin and TGF-β axes. FDR is computed inside that panel, separately from the 10 housekeeping genes, so the controls neither dilute nor inflate the correction.

Ligand-only module scores are the mean of per-gene z-scores. F-TBRS is the Mariathasan 2018 fibroblast TGF-β response signature, used as an activity read-out.

## Statistics

- **A1.** Spearman ρ (Fisher-z 95% CI) of each gene with TACSTD2. Tertile split at 1/3 and 2/3 of TACSTD2; high vs low compared with Mann-Whitney / Cliff's δ.
- **A2.** Partial Spearman: rank both variables and the covariate, residualize the covariate by OLS, correlate the residuals. Primary covariate is ABSOLUTE purity. RNA epithelial score is a sensitivity bound (TACSTD2 is itself epithelial, so this can over-adjust).
- **A3.** A1 repeated inside tertiles of ABSOLUTE purity.
- **A4.** Identical pipeline on ACTB, GAPDH, RPL13A, TBP, PPIA, B2M, UBC, SDHA, HPRT1, PGK1.
- Multiple testing: Benjamini-Hochberg inside the primary panel, and separately inside the controls.
- LUAD and LUSC are never pooled. LUSC is replication.

## Post-hoc (labelled as such)

Housekeeping-score residualization of CD47 vs TACSTD2, added after A4 showed that in LUAD, ACTB and UBC correlate with TACSTD2 at least as strongly as CD47 does. Not used for the pre-registered call.

## Single cell (GSE131907)

Author-provided log2TPM matrix and cell annotation. Panel genes streamed out of the 2.9 GB matrix; the full matrix is never loaded.

- **S1.** Mean and percent-expressing per `Cell_type` in tumour-origin cells (`tLung`, `tL/B`, `mLN`, `PE`).
- **S2.** Patient-level mean inside malignant epithelium, then Spearman of that mean against TACSTD2. Primary test is tLung tS1/tS2/tS3 only. LN/bronchus mets, brain mets, and a mixed-site slice are reported separately because both genes drop at brain mets and pooling sites inflates ρ. Patients with <20 (site slices) or <30 (mixed) malignant cells are dropped.
- **S3.** tLung cell-type fractions vs malignant TACSTD2. n = 10; descriptive.

GSE127465 was listed in the plan and was not run.

## Software

Python 3.12, pandas, numpy, scipy. No permutation p-values; CIs are Fisher-z on the reported ρ, which is slightly anti-conservative for partial correlations (df is n−3, CI uses n−1−k). At n ≈ 500 this does not change any call.
