# Claim A2 verification: TACSTD2 (TROP2) vs immune score in durvalumab-treated NSCLC

**User's claim:** Spearman ρ = **-0.65** (raw), **-0.46** (tumor-purity-adjusted), **p = 2e-4**.

## Verdict: substantially reproduced (with caveats on the p-value)

The negative TACSTD2–immune correlation reproduces in both public durvalumab NSCLC
cohorts, and the post-treatment cohort matches the claimed coefficients almost exactly:

| Dataset | n | Raw Spearman ρ | Raw p | Purity-adjusted ρ | Adjusted p |
|---|---|---|---|---|---|
| **Claimed** | ? | **-0.65** | p = 2e-4 (unspecified which) | **-0.46** | — |
| GSE248378 (post-treatment resections) | 29 | **-0.665** | 8.3e-5 | **-0.474** | 1.1e-2 |
| GSE253564 (pre-treatment biopsies) | 32 | -0.716 | 4.0e-6 | -0.530 | 2.2e-3 |

- **Raw ρ:** claimed -0.65; observed -0.665 (post-treatment) and -0.716 (pre-treatment). Match.
- **Purity-adjusted ρ:** claimed -0.46; observed -0.474 (post) and -0.530 (pre). Match
  (post-treatment within 0.014 of the claim).
- **p = 2e-4:** ambiguous in the claim. It is consistent in order of magnitude with the
  *raw* correlation (we get 8.3e-5 post, 4.0e-6 pre). It does **not** describe the
  *purity-adjusted* correlation, which is weaker (p = 0.011 post, p = 0.0022 pre).
  If the claim attaches p = 2e-4 to the adjusted ρ = -0.46, that specific pairing does
  not reproduce at these sample sizes (ρ = -0.46 would need n ≈ 55+ for p = 2e-4).

The effect is also robust to the immune definition: PTPRC (CD45), CD8A, and cytolytic
score (mean log2 of GZMA/PRF1) all correlate negatively with TACSTD2 (see
`results.json`). In the post-treatment cohort every alternative remains significant
after purity adjustment; in the pre-treatment cohort CD8A and CYT attenuate to
non-significance after adjustment (ρ ≈ -0.24 and -0.17), so some of the baseline
signal there is attributable to tumor purity.

## Data provenance

The only public human bulk-RNA-seq datasets of durvalumab-treated NSCLC we could
locate are two GEO series from the same randomized phase II neoadjuvant trial
(durvalumab ± stereotactic radiation, Altorki et al., NCT02904954):

- [GSE253564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564) — pre-treatment samples, FPKM, n = 32
- [GSE248378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378) — post-treatment resected tumors, FPKM, n = 29

We cannot confirm which dataset/definition the user used. The closest published
analogue (Bessede et al., Clin Cancer Res 2024, "TROP2 Is Associated with Primary
Resistance to Immune Checkpoint Inhibition…") used **atezolizumab** (OAK/POPLAR)
data under controlled access at EGA (EGAS00001005013), which is not publicly
downloadable, so it cannot be the source of a reproducible durvalumab claim.

## Methods

- **Immune score / purity:** exact Python port of the ESTIMATE R package v1.0.13
  (rank-normalized ssGSEA, weight exponent 0.25; official `SI_geneset.gmt`
  signatures; 138/141 immune genes present in both matrices).
- **Raw correlation:** Spearman of log2(FPKM+1) TACSTD2 vs ESTIMATE ImmuneScore.
- **Purity adjustment:** partial Spearman (Pearson on ranks, p from t with n-3 df)
  controlling for ESTIMATEScore ranks. ESTIMATE TumorPurity =
  cos(0.6049872018 + 0.0001467884·ESTIMATEScore) is strictly monotone decreasing,
  so rank-based adjustment on ESTIMATEScore is identical to adjusting on purity.
  This matters because the Affymetrix-calibrated purity formula goes out of bounds
  (returns NaN) for 23/32 pre-treatment and 6/29 post-treatment FPKM samples;
  the rank-equivalent adjustment keeps all samples.

## Caveats

1. ESTIMATE's absolute purity values are unreliable on these FPKM data (many out of
   bounds); only purity *ranks* are used, which is all a partial Spearman needs.
2. Both cohorts are small (n = 29–32) and come from one early-stage neoadjuvant
   trial; the claim's provenance (cohort, immune definition, purity method) was
   not specified, so this is a reproduction of the claim's substance, not
   necessarily its exact pipeline.
3. TACSTD2 (TROP2) is epithelial-restricted, so a negative correlation with immune
   infiltration partly reflects tumor-cell content; that is precisely why the
   purity-adjusted estimate (≈ -0.47 to -0.53, still significant) is the more
   meaningful number, consistent with the user's -0.46.

## Files

- `results.json` — all correlations (primary + robustness immune definitions)
- `scores_GSE253564_pretreatment.csv`, `scores_GSE248378_posttreatment.csv` —
  per-sample TACSTD2 expression, ESTIMATE scores, purity
- `tacstd2_vs_immune.png` — scatter plots, both cohorts
- Analysis code: `scripts/claim_A2_analysis.py` (repo root)

## Reproduce

```bash
pip install pandas scipy numpy matplotlib
# data/: GEO supplementary FPKM files + ESTIMATE v1.0.13 gene sets (see scripts)
python3 scripts/claim_A2_analysis.py
```
