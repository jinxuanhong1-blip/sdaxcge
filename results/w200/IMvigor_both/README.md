# IMvigor210: TACSTD2 and CLDN4 vs ORR and OS

**Bottom line: no significant association for TACSTD2, CLDN4, or the
both-high (AND) gate.** In this public atezolizumab urothelial-carcinoma
cohort, neither gene (nor their mean z-score, nor TACSTD2-high ∩ CLDN4-high)
separated confirmed responders from non-responders, and neither associated
with overall survival. Immune positive controls (CXCL9, CD8A) do associate
in the expected direction, so the endpoints are not inert.

一句话结论：公开 IMvigor210（尿路上皮癌，atezolizumab）里 TACSTD2 / CLDN4
以及双高（AND）与 ORR、OS 均无显著关联；CXCL9 / CD8A 对照方向正确。

## Cohort and data

- Source: IMvigor210CoreBiologies v1.0.0 (Mariathasan et al., *Nature* 2018,
  PMID 29443960), Creative Commons 3.0. Official tarball:
  `http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz`
- Disease / drug: **metastatic urothelial carcinoma**, atezolizumab (anti-PD-L1).
  This is **not** a lung ICI cohort. Biopsy site is recorded (`Tissue`); 10 of
  348 libraries are from lung metastases of urothelial cancer.
- RNA: 31,286 genes × 348 libraries (IGIS Entrez annotation). Expression is
  **log2(TPM+1)** using package gene lengths. TACSTD2 = Entrez 4070,
  CLDN4 = 1364.
- ORR: author `binaryResponse` = confirmed **CR/PR** (n=68) vs
  **SD/PD** (n=230). **NE (n=50) excluded** from ORR tests. This is
  the package's own confirmed-ORR coding, not a recode.
- OS: `os` (months) and `censOS` (1 = death). All 348 libraries have OS;
  one patient contributed two NE libraries with identical OS — patient-level
  OS uses n=347 (larger sizeFactor library kept).
- Primary tests: TACSTD2 and CLDN4 only. BH-FDR is within endpoint
  (2 genes for ORR; 2 genes for continuous Cox OS). Combined z, both-high,
  and CXCL9/CD8A are labeled exploratory / control.

## ORR results (n=298 evaluable)

| Feature | Median CR/PR | Median SD/PD | MWU p | BH q (2 genes) | AUC (CR/PR higher) | AUC 95% CI | rank-biserial |
|---|---:|---:|---:|---:|---:|---|---:|
| TACSTD2 log2(TPM+1) | 7.65 | 7.46 | 0.227 | 0.239 | 0.55 | 0.47–0.63 | 0.10 |
| CLDN4 log2(TPM+1) | 4.73 | 4.41 | 0.239 | 0.239 | 0.55 | 0.47–0.63 | 0.09 |
| Combined z (exploratory) | 0.22 | 0.06 | 0.189 | — | 0.55 | 0.47–0.63 | 0.10 |
| CXCL9 (control) | 4.62 | 3.39 | 2.84e-04 | — | 0.64 | 0.57–0.72 | 0.29 |
| CD8A (control) | 3.07 | 2.77 | 0.033 | — | 0.58 | 0.51–0.66 | 0.17 |

AUC > 0.5 means **higher** expression in confirmed responders. TACSTD2 and
CLDN4 point very slightly that way (AUC 0.55 and
0.55); both 95% CIs include 0.5 and both p-values are
~0.2. That is **not** support for a TROP2/CLDN4-high → worse-ORR claim, and
it is also **not** a significant responder-high finding.

Both-high (each gene ≥ its cohort median) vs rest: ORR
26.0% (104 pts) vs
21.1% (194 pts), Fisher OR
1.31, p=0.386.
Both-high vs both-low: ORR 26.0% vs
18.0%, Fisher p=0.181.

TACSTD2 and CLDN4 co-express (Spearman ρ=0.56,
p=3.29e-26 on ORR-evaluable samples).

## OS results (n=347 patients)

Cox PH on **per-1-SD** log2(TPM+1); HR > 1 = higher expression, shorter OS.

| Feature | n / events | HR per 1 SD | 95% CI | p | BH q (2 genes) | median log-rank p |
|---|---|---:|---|---:|---:|---:|
| TACSTD2 | 347 / 231 | 0.97 | 0.86–1.09 | 0.640 | 0.640 | 0.103 |
| CLDN4 | 347 / 231 | 0.96 | 0.85–1.09 | 0.531 | 0.640 | 0.349 |
| Combined z (exploratory) | 347 / 231 | 0.96 | 0.85–1.09 | 0.542 | — | 0.014 |
| CXCL9 (control) | 347 / 231 | 0.75 | 0.66–0.85 | 1.21e-05 | — | 0.002 |
| CD8A (control) | 347 / 231 | 0.83 | 0.73–0.95 | 0.006 | — | 0.480 |

Both-high vs rest log-rank p=0.116
(n_high=120). IC/Sex/platinum-adjusted Cox HRs for
TACSTD2 and CLDN4 remain ~1 (see `os_cox_adjusted.csv`).

The exploratory **combined-z median** log-rank is p=0.014
(high combined z has *longer* observed OS). That split is **not** a primary
result: the pre-specified continuous Cox for the same score is HR
0.96 (p=0.542). Do not
quote the median p-value as an OS hit.

## Honest caveats

- **Wrong disease for a lung-TROP2 claim.** IMvigor210 is urothelial
  carcinoma. A null here does not prove a null in NSCLC, and a hit here
  would not have transferred either.
- **Null is not proof of no effect.** Effects smaller than about AUC 0.60
  or OS HR ~0.85/1.18 per SD are compatible with these confidence intervals.
  What is *not* supported is a large, consistent TACSTD2/CLDN4–ORR/OS link
  in this public matrix.
- **Direction vs the usual TROP2-cold story.** Point estimates are slightly
  *higher* expression in responders and HR ≈ 1 for OS — the opposite of
  “TROP2/CLDN4-high → worse ICI,” and still non-significant.
- Bulk TPM. No purity residual in the primary test. Epithelial genes are
  diluted by stroma; that can mask a tumor-intrinsic effect (see other
  slices) but does not license reading this as a positive ICI biomarker.
- Association on an all-treated cohort. There is no chemotherapy/SOC arm
  in these public tables, so this is **not** a predictive interaction test.
- Follow-up is capped near 24.5 months; 50 NE patients are out of ORR.
- CXCL9/CD8A are controls that the assay/endpoints can move, not a license
  to hunt other genes post hoc.

## Reproduction

```bash
python3 -m pip install -r scripts/w200/IMvigor_both/requirements.txt
python3 scripts/w200/IMvigor_both/download.py
python3 scripts/w200/IMvigor_both/analyze.py
```

`download.py` pulls the official 122 MB tarball (or reuses a cached copy)
and uses `extract_cds.R` to dump `pData`, `fData`, and counts. Raw files
stay under `data/IMvigor210/` (gitignored). Seed 20260816 for bootstrap.

## Files

- `README.md` — this write-up (numbers filled by `analyze.py`)
- `summary.json` — machine-readable verdict and key stats
- `sample_data.csv` — per-library clinical + expression (no raw counts)
- `orr_stats.csv`, `os_cox.csv`, `os_cox_adjusted.csv`, `os_logrank.csv`
- `both_high_orr.csv`, `quadrant_orr.csv`, `spearman_tacstd2_cldn4.csv`
- `provenance.tsv` — source URL, expected size, sha256 if available
- `orr_boxplots.png/.pdf`, `os_km.png/.pdf`, `tacstd2_cldn4_scatter.png/.pdf`
