# UCSC Xena TCGA-LUAD: TACSTD2 / CLDN4 vs immune after ESTIMATE and ABSOLUTE purity

**Status:** pipeline written; numbers will be filled after the first run.

**Question:** In TCGA-LUAD primary tumors from UCSC Xena, are TACSTD2 (TROP2) and CLDN4 associated with immune features after adjusting for tumor purity, using both ESTIMATE (RNA) and ABSOLUTE (DNA) purity?

This is not a re-run of `results/w200/A1_LUAD/` (TACSTD2 only, ABSOLUTE only). It adds CLDN4, ESTIMATE vs ABSOLUTE side-by-side, and a GDC STAR TPM sensitivity matrix.

## Methods (locked before looking at results)

- **Primary expression:** UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2`, log2(RSEM normalized count + 1), `-01` primary tumors. This is the RNAseqV2 freeze MD Anderson used for published ESTIMATE scores.
- **Sensitivity expression:** UCSC Xena GDC hub `TCGA-LUAD.star_tpm.tsv.gz`, log2(TPM+1), GENCODE v36 symbols, collapsed to patient-level `-01`.
- **ESTIMATE purity:** MD Anderson LUAD RNAseqV2 scores; purity = `cos(0.6049872018 + 0.0001467884 · ESTIMATE_score)` (Yoshihara 2013). Immune_score and Stromal_score are also kept as features.
- **ABSOLUTE purity:** PanCanAtlas `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` (SNP-array / DNA). Aran 2015 ESTIMATE/ABSOLUTE/CPE used only as a cross-check.
- **Immune features:** methylation leukocyte fraction (Thorsson/PanImmune; RNA-independent); ESTIMATE Immune/Stromal scores; 5 Wolf signatures; 5 CIBERSORT relative fractions; CYT (mean GZMA/PRF1), CD8 score (mean CD8A/CD8B), Ayers GEP18; 21 marker genes.
- **Stats:** Spearman; partial Spearman (Pearson on rank residuals). BH-FDR within each (matrix × predictor × purity method) block. Primary n = samples with HiSeqV2 + both purities.
- **Circularity (pre-declared):** ESTIMATE Immune/Stromal scores are terms in ESTIMATE_score, so partialling ESTIMATE purity out of those scores is circular. ESTIMATE purity is itself RNA-derived, so adjusting RNA immune genes for ESTIMATE purity is only partly independent. ABSOLUTE (DNA) and methylation leukocyte fraction are the non-circular readouts.

Reproduce: `bash scripts/xena_luad_download.sh && python3 scripts/xena_luad_tacstd2_cldn4_immune.py`

## Results

*Pending first run.*

## Verdict

*Pending first run. Will be stated honestly from the numbers, not from the hypothesis.*
