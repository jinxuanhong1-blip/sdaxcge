# Methods

1. **Cohorts.** Public processed scRNA:
   - GSE131907 (Kim et al., *Nat Commun* 2020). Author `Cell_type` / `Cell_subtype`. Tumor-origin samples only (`tLung`, `tL/B`, `mLN`, `mBrain`, `PE`).
   - GSE205335 (Hu et al. lung IO atlas). Author `lineage.sub` / `lineage.total` from the GEO CellIdentity table.
   - **GSE207422 is not analyzed here.** PR #334 already ran CLDN4-only NicheNet-style scoring on that n=12 cohort; signed patient tests were NS.

2. **Unit.** Patient. GSE205335 is already patient-level. GSE131907 GEO files are sample-level (58 samples / 44 patients); eligible samples are collapsed to unique `patient_id` from the series matrix. Sample n is reported separately. Cells are never the independent unit.

3. **Eligibility.** ≥20 author-malignant cells and ≥20 T/NK cells after the patient collapse. Same gate as the GSE205335 lock in PR #320.

4. **Sender / receiver.**
   - Malignant: GSE131907 = tumor-origin ∩ `{Malignant cells, tS1, tS2, tS3}`; GSE205335 = `lineage.sub == Malignant cells`.
   - T/NK: GSE131907 = `T lymphocytes` + `NK cells`; GSE205335 = `lineage.total == T/NK cells`.
   - CLDN4-high = malignant `log1p(CP10k)` ≥ the **cohort** malignant-cell median. TACSTD2 is not part of the call. A companion dual-high count is reported and unused.

5. **Normalization.** `log1p(1e4 * UMI / library size)`.

6. **Prior.** NicheNet-v2 human `ligand_target_matrix_nsga2r_final.rds` and `lr_network_human_21122021.rds` (Zenodo 10.5281/zenodo.7074291; Browaeys et al.). Converted with Python `rdata`. **R / `nichenetr` is not run.** Ligand activity = Pearson, AUROC, and AUPR of the ligand’s prior target scores versus gene-set membership on background genes. Rank by Pearson.

7. **Potential ligand.** Detected in ≥10% of CLDN4-high malignant cells, ≥1 receptor detected in ≥10% of T/NK, and a column in the v2 ligand–target matrix.

8. **Background.** T/NK-expressed extracted genes ∩ prior target rows (plus the a priori program genes). This is narrower than a full-transcriptome `nichenetr` run and is reported as such.

9. **Gene sets.** A priori IFN and cytotoxicity (`scripts/gene_sets.py`). Exhaustion is sensitivity only. Empirical CLDN4-associated T/NK genes: per-patient T/NK gene means, Spearman vs malignant CLDN4 %pos within cohort, DerSimonian–Laird on Fisher-z, p < 0.15 split by sign. That cutoff is honest; it is not FDR<0.05.

10. **Signed patient tests.** Spearman of malignant CLDN4 %pos vs same-patient T/NK IFN / cytotoxicity / fraction. Q4 vs Q1 is Mann–Whitney on the quartile tails with rank-biserial *r*; thin tails (n<8 or a tail <3) are not Fisher-z pooled. Two-cohort merge = DerSimonian–Laird on Fisher-z. Per-patient ligand means in CLDN4-high cells are correlated with the same T/NK programs and merged the same way.

11. **Software.** Python 3; pandas / numpy / scipy / scikit-learn / matplotlib / rdata. Scripts: `00_download.py` → `03_analyze.py`.
