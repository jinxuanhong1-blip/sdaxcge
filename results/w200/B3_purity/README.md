# B3 — TCGA-LUAD: tight-junction (TJ) score vs CD8 score, adjusting for ESTIMATE purity

**Question.** In TCGA lung adenocarcinoma, is the epithelial tight-junction (TJ)
signature associated with CD8 T-cell abundance, and does any association survive
adjustment for tumor purity? TJ and CD8 are both expected to track purity (more
tumor epithelium → higher TJ, lower relative immune signal), so the honest test is
the **partial** correlation controlling for ESTIMATE purity, reported alongside the
raw one.

Reproduce with: `python3 scripts/b3_luad_tj_cd8_purity.py`

## Data

- **Expression:** UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` — gene-level
  `log2(norm_count+1)`, 20,530 genes × 576 samples.
- **Purity:** MD Anderson precomputed **ESTIMATE** scores for TCGA-LUAD RNAseqV2
  (`Stromal_score`, `Immune_score`, `ESTIMATE_score`), 517 samples.
- **Matched primary tumors (sample code `-01`) in both:** **515** (all with a
  defined purity → 515 complete cases).

## Methods

- **TJ score** = mean of per-gene z-scores (across samples) of a curated set of
  **structural** tight-junction genes (all 15 found in the matrix): claudins
  `CLDN1/3/4/7`, occludin `OCLN`, ZO scaffolds `TJP1/2/3`, JAMs `F11R/JAM2/JAM3`,
  MARVEL proteins `MARVELD2/MARVELD3`, cingulins `CGN/CGNL1`. The set is
  deliberately structural (no broad signaling genes). See `tj_genes_used.csv`.
- **CD8 score** = mean of per-gene z-scores of `CD8A`, `CD8B` (both found). See
  `cd8_genes_used.csv`.
- **ESTIMATE purity** = `cos(0.6049872018 + 0.0001467884 · ESTIMATE_score)`
  (Yoshihara et al. 2013).
- **Correlation** = Spearman. **Partial** Spearman = rank-transform TJ, CD8 and
  purity, regress the TJ and CD8 ranks on the purity rank, and Pearson-correlate
  the residuals. **95% CIs** = 2,000× percentile bootstrap (seed 20260816).

## Results (n = 515)

| Test | rho | 95% CI | p |
|---|---|---|---|
| **Unadjusted** TJ vs CD8 | **−0.291** | (−0.370, −0.207) | 1.6e−11 |
| **Partial** TJ vs CD8 \| purity | **−0.257** | (−0.333, −0.173) | 3.4e−09 |
| Context: TJ vs purity | +0.144 | — | 1.1e−03 |
| Context: CD8 vs purity | −0.563 | — | 2.0e−44 |
| Sensitivity: partial \| Immune_score | −0.225 | — | 2.5e−07 |
| Sensitivity: partial \| Stromal_score | −0.288 | — | 2.6e−11 |

**Interpretation (honest).** TJ and CD8 are **modestly, significantly negatively
correlated** (rho ≈ −0.29). Purity confounding is real but **not** the explanation:
CD8 tracks purity strongly (rho = −0.56) while TJ tracks it only weakly (rho =
+0.14), and after removing purity the TJ–CD8 association is only mildly attenuated
(rho = −0.26, still highly significant). The negative sign is robust when instead
adjusting for the ESTIMATE Immune or Stromal score. Effect sizes are modest (each
signature explains only a few percent of the other's rank variance), so this is a
consistent weak inverse relationship, not a strong one.

## Caveats

- The ESTIMATE purity formula was calibrated on Affymetrix data; applied to RNAseqV2
  ESTIMATE scores it is a widely-used but approximate purity proxy. The Immune/Stromal
  sensitivity analyses give the same qualitative answer, so the conclusion does not
  hinge on the exact purity transform.
- Signature scores are simple mean-z-scores, not ssGSEA; the CD8 score uses only the
  lineage coreceptor chains and the TJ score is a specific structural gene set —
  different (equally defensible) gene sets could shift magnitudes but the direction is
  stable across the adjustments tested here.
- Bulk RNA-seq: these are tissue-level associations, not single-cell/spatial
  relationships. No causal claim is made.

## Files

- `sample_scores.csv` — per-sample TJ, CD8, purity, ESTIMATE scores (n=515).
- `correlation_results.json` / `.csv` — all statistics above.
- `tj_genes_used.csv`, `cd8_genes_used.csv` — signature membership + found flags.
- `scatter_tj_cd8.png` — TJ vs CD8, colored by purity.
- `scatter_partial_residuals.png` — purity-residualized TJ vs CD8.
- `purity_associations.png` — TJ~purity and CD8~purity.
