# GSE135222: TACSTD2 (TROP2) and CLDN4 vs DCB

**Bottom line: TACSTD2 and CLDN4 do not associate with DCB or PFS in this
cohort. The null is more informative than a simple underpowered shrug:
canonical ICI genes (IFNG, CD8A, immune composite) *do* separate DCB from
NDB here, so the endpoint and sample size are sufficient to detect a known
large effect. TACSTD2/CLDN4 show none of that.**

## Cohort and labels

- GEO [GSE135222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222):
  RNA-seq (RSEM TPM, hg19) of 27 advanced NSCLC patients treated with
  anti-PD-1/PD-L1 (Jung et al., *Nat Commun* 2019, PMID 31537801). This is
  the RNA-seq subset of the Samsung Medical Center n=60 cohort also reported
  in Kim et al., *Clin Epigenetics* 2020, PMID 32762727.
- **Primary DCB** = author-assigned `benefit` in Kim 2020 Supplementary
  Table 1 (Y/N), matched by patient ID. Jung/Kim define DCB as RECIST v1.1
  PR or SD lasting >6 months — not a pure PFS-time cutoff.
  Result: **8 DCB / 19 NDB**.
- GEO itself only provides PFS event and `pfs.time` (days). A PFS ≥ 183 d
  rule gives 7 DCB / 20 NDB and disagrees with the authors on **one**
  RNA-seq patient: NSCLC1708 (PFS 168 d, event, author DCB=Y). That
  discrepancy is in `author_vs_pfs183_discrepancy.csv`. Sensitivity across
  PFS cutoffs (90 / 168 / 183 / 365 d) does not change the TACSTD2/CLDN4
  null.

## Primary result (author DCB)

| Feature | Role | Median DCB | Median NDB | MWU p | BH q (2 genes) | AUC (DCB) | AUC 95% CI |
|---|---|---|---|---|---|---|---|
| TACSTD2 log2(TPM+1) | primary | 7.46 | 7.83 | 1.00 | 1.00 | 0.50 | 0.29–0.72 |
| CLDN4 log2(TPM+1) | primary | 7.45 | 7.23 | 0.90 | 1.00 | 0.52 | 0.30–0.75 |
| Combined z | exploratory | 0.35 | 0.42 | 0.77 | — | 0.46 | 0.25–0.68 |
| CD8A | +control | 3.69 | 2.41 | 0.022 | — | 0.78 | 0.53–0.97 |
| CXCL9 | +control | 4.96 | 2.91 | 0.051 | — | 0.74 | 0.55–0.91 |
| GZMB | +control | 4.52 | 3.65 | 0.095 | — | 0.71 | 0.51–0.89 |
| CD274 | +control | 3.18 | 2.62 | 0.13 | — | 0.69 | 0.47–0.88 |
| IFNG | +control | 1.28 | 0.26 | 0.011 | — | 0.82 | 0.63–0.97 |
| Immune z (CD8A/CXCL9/GZMB) | +control | 0.76 | −0.24 | 0.0065 | — | 0.83 | 0.62–0.99 |

TACSTD2 AUC is exactly 0.50. CLDN4 is noise. The immune composite is the
strongest separator in the file (AUC 0.83).

## Survival (does not use the DCB cutoff)

Univariate Cox on z-scored log2(TPM+1); log-rank on a median split. 21 PFS
events / 27 patients.

| Feature | Cox HR per 1 SD | Cox 95% CI | Cox p | Log-rank p (median) |
|---|---|---|---|---|
| TACSTD2 | 1.07 | 0.68–1.69 | 0.78 | 0.43 |
| CLDN4 | 1.12 | 0.72–1.75 | 0.61 | 0.70 |
| Combined z | 1.11 | 0.70–1.75 | 0.67 | 0.21 |
| IFNG | 0.58 | 0.35–0.97 | 0.039 | 0.039 |
| Immune z | 0.68 | 0.45–1.03 | 0.072 | 0.043 |
| CXCL9 | 0.69 | 0.44–1.08 | 0.10 | 0.034 |
| GZMB | 0.68 | 0.43–1.07 | 0.096 | 0.044 |
| CD8A | 0.73 | 0.49–1.08 | 0.12 | 0.055 |

Higher TACSTD2/CLDN4, if anything, trends toward *worse* PFS (HR > 1) and
is nowhere near significant. IFNG and the immune composite trend the
expected direction (higher expression, longer PFS).

Spearman vs PFS time (ignores censoring): TACSTD2 r = −0.16 (p = 0.44),
CLDN4 r = −0.14 (p = 0.48), IFNG r = +0.42 (p = 0.028). TACSTD2 and CLDN4
are only modestly coexpressed (ρ = 0.35, p = 0.073).

## Sensitivity and influence

- PFS ≥ 90 / 168 / 183 / 365 d: TACSTD2 and CLDN4 remain null (all MWU
  p > 0.2; all AUC CIs include 0.5). The 365 d split is 2 vs 25 and is
  not a serious test.
- Leave-one-out: dropping any single patient, including the long-PFS
  low-expresser NSCLC947 or the label-discrepant NSCLC1708, never moves
  TACSTD2 or CLDN4 to p < 0.5. No one sample is propping up the null.

## Power (honest)

With 8 vs 19, two-sided MWU, α = 0.05, equal-variance normal simulation
(2000 draws):

| True AUC | Power |
|---|---|
| 0.60 | 0.12 |
| 0.70 | 0.34 |
| 0.80 | 0.72 |
| 0.85 | 0.88 |

80% power requires a true AUC around 0.83 — exactly the size of the
immune-z effect that *is* detected. A modest TACSTD2/CLDN4 effect
(AUC 0.65–0.70) would usually be missed. What this cohort *can* rule
out is a CD8A/IFNG-sized association for these two genes.

## Honest caveats

- n = 27 is small. The positive-control success means we are not
  dismissing a dead assay, but it does not license a claim of “no
  association in NSCLC ICI.” It licenses: **no large association in this
  public RNA-seq slice.**
- Author DCB is RECIST-based; we do not have raw CR/PR/SD/PD categories
  for the RNA-seq patients, only the Y/N benefit flag plus GEO PFS.
- Bulk TPM. No histology, PD-L1 IHC, TMB, or purity in the GEO metadata.
  A purity-driven epithelial signal could in principle hide a tumor-cell
  effect; it would not explain why IFNG/CD8A still work.
- Controls were chosen a priori as standard ICI genes, not mined from
  this matrix. They are not part of the two-gene primary family (BH is
  across TACSTD2 and CLDN4 only).
- No multiple-testing correction is applied to the control or survival
  secondary tests; those p-values are descriptive.

## Files

- `sample_data.csv` — per-patient clinical (sex, age, PFS, author DCB,
  PFS≥183 DCB) and expression.
- `stats.csv` — primary MWU / AUC / BH q plus controls.
- `sensitivity_dcb_cutoffs.csv` — author vs PFS 90/168/183/365 d.
- `author_vs_pfs183_discrepancy.csv` — NSCLC1708 only.
- `cox_km.csv` — Cox HR per SD and median-split log-rank.
- `spearman_pfs.csv` — Spearman vs PFS time; TACSTD2–CLDN4 coexpression.
- `leave_one_out.csv` — influence of each sample on primary AUCs.
- `power.csv` — MWU power at 8 vs 19.
- `boxplots_dcb.png`, `boxplots_controls.png`, `roc_dcb.png`,
  `km_curves.png`, `scatter_tacstd2_cldn4.png`.
- Code: `scripts/gse135222_tacstd2_cldn4_dcb.py` (downloads GEO + Kim 2020
  Table 1; hardcoded author labels as fallback; seed 20260816).
