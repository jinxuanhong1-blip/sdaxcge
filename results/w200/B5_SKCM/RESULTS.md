# B5_SKCM results — CLDN4 vs ICI response in melanoma

**Verdict: null / inconclusive.** Pre-treatment tumor `CLDN4` mRNA is not associated with ICI response in the three open melanoma RNA-seq cohorts that carry a response label. The result does not meet the positive-result rule written in `ANALYSIS_PLAN.md` (1 of 4 criteria only, and that one is weak).

This is not a well-powered negative. The same pipeline does not recover the known ICI-benefit genes (`CD8A`, `CD274`, `IFNG`, Ayers IFN-γ score) at p < 0.05 either. A null for `CLDN4` here is therefore **uninformative as evidence of no effect**, and should not be cited as a replication failure or as a discovery.

`CLDN4` is also a low-information exposure in this disease: median expression is well below 1 FPKM/TPM in Hugo, Riaz, and TCGA-SKCM, as expected for an epithelial tight-junction gene in a neural-crest tumor.

---

## 1. What was tested

Question, cohorts, sample rules, transforms, and the four-part positive-result rule are in `ANALYSIS_PLAN.md`, committed in this repository *before* any response-associated test was run.

Primary analysis set (pre-treatment, one patient, binary R vs NR):

| Cohort | Therapy | Unit | Patients | R | NR |
|---|---|---|---:|---:|---:|
| Hugo GSE78220 | pembrolizumab | author FPKM | 26 | 14 | 12 |
| Riaz GSE91061 | nivolumab | author FPKM | 33 | 10 | 23 |
| MGH GSE115821 | anti-PD-1 / anti-CTLA-4 | TPM from counts | 8 | 2 | 6 |
| **Total** | | | **67** | **26** | **41** |

Riaz `SD`/`UNK` are excluded from the primary set (plan §4). MGH contributes 8 first-episode patients only; 2 responders. Hugo `Pt27A`/`Pt27B` were averaged. Hugo `Pt16` is on-treatment and was excluded.

Controlled-access series (Gide, Liu, Van Allen) were not used. That is a coverage gap, not a silent exclusion after seeing results.

## 2. `CLDN4` is low in melanoma

`log2(x + 0.1)` medians on the primary set:

| Cohort | Median log2 | Approx. linear | Fraction with log2 < 0 (i.e. x < 0.9) | Low-information flag |
|---|---:|---:|---:|---|
| Hugo | −1.19 | ~0.34 FPKM | 73% | yes |
| Riaz | −2.15 | ~0.13 FPKM | 79% | yes |
| MGH | +0.23 | ~1.1 TPM | 50% | no |
| TCGA-SKCM tumors (n=471, context only) | −0.61 | ~0.56 TPM | 60% | — |

Reads are not zero — almost every sample has a detectable count — but the dynamic range is a spike-and-floor: most tumors sit near the offset, a few are high. That is a poor continuous biomarker. The high tail tracks epithelial / keratinocyte genes, not melanocyte genes (section 6).

## 3. Primary association: no signal

Per-cohort `log2 CLDN4`, R vs NR:

| Cohort | n R / NR | MWU p | AUC (P[R>NR]) | Logistic OR per 1 SD (95% CI) | Logistic p |
|---|---|---:|---:|---|---:|
| Hugo | 14 / 12 | 0.74 | 0.54 | 1.20 (0.53–2.71) | 0.66 |
| Riaz | 10 / 23 | 0.49 | 0.58 | 0.69 (0.29–1.63) | 0.39 |
| MGH | 2 / 6 | 0.64 | 0.67 | 1.37 (0.25–7.50) | 0.72 |

Fixed-effect inverse-variance meta-analysis of the logistic log-odds ratios:

- **OR per 1 SD = 0.96 (95% CI 0.55–1.69), p = 0.90**
- Random-effects identical (I² = 0%, Q p = 0.60, τ² = 0)
- Direction: Hugo and MGH OR > 1, Riaz OR < 1

The rank tests and the logistic model do not even agree on direction in Riaz (AUC 0.58 says responders sit slightly higher; the logistic OR 0.69 says the opposite). That is what a right-skewed, near-zero gene does in n = 33. It is not a finding.

MGH is two responders. Its interval runs from 0.25 to 7.5 and should not be read as “consistent with Hugo.”

One-stage IPD logistic with cohort as a fixed effect (sensitivity S5) is the same answer: OR 0.95 (0.56–1.60), p = 0.85, n = 67.

## 4. Pre-specified positive-result rule

All four were required. Observed:

| Criterion | Required | Observed |
|---|---|---|
| Meta p < 0.05 | yes | **no** (p = 0.90) |
| Same direction in ≥ 2 of 3 cohorts | yes | yes, technically (Hugo + MGH) |
| `CLDN4` in the top 5% of the genome-wide meta ranking | yes | **no** (rank 15,098 / 17,680, 85th percentile) |
| Association survives keratinocyte-score adjustment | yes | **no** (adjusted meta OR 0.88, p = 0.74; MGH model did not converge) |

**4/4 not met. Not a positive result.**

The one criterion that is “yes” is not evidence. Two of the three point estimates are > 1, one of those two is n = 8, and the meta-estimate is 0.96.

## 5. Genome-wide context

The identical per-cohort logistic + IVW meta was run on 17,680 genes present in all three cohorts.

- 235 genes have nominal meta p < 0.05 (expected ~884 under a uniform null; the deficit is consistent with low power, not with a large true-positive set)
- **0 genes** have Benjamini–Hochberg FDR < 0.05
- `CLDN4` FDR = 1.00; it is in the *least* associated fifth of the transcriptome
- The top of the list (WTIP, MMP2, ANKRD55, MFAP2, CLDN16, …) is not an immune-effector module. Nothing in that tail is treated as a hit. `CLDN16` is listed only so it is not later “discovered” from these files; it was not pre-specified and it does not survive FDR.

## 6. Positive controls failed. The assay is underpowered.

Same pipeline, same patients:

| Gene / score | Meta OR per SD | Meta p | Direction + / − |
|---|---:|---:|---|
| CD274 (PD-L1) | 1.43 | 0.21 | 3 / 0 |
| IFNG | 1.34 | 0.28 | 3 / 0 |
| Ayers IFN-γ 6-gene | 1.34 | 0.30 | 2 / 1 |
| CXCL9 | 1.33 | 0.32 | 2 / 1 |
| CD8A | 1.22 | 0.48 | 2 / 1 |
| GZMB | 1.13 | 0.66 | 2 / 1 |
| **CLDN4** | **0.96** | **0.90** | **2 / 1** |

The immune genes sit on the expected side of 1 and still do not reach 0.05. That is the power the plan warned about: n ≈ 67 with a ~40% response rate cannot detect a single-gene OR of ~1.3–1.4. `CLDN4` is weaker than all of them. A powered test of `CLDN4` would need either a much larger open cohort or a larger true effect than anything seen here.

## 7. `CLDN4` tracks epithelium, not melanocytes

Spearman ρ of `CLDN4` vs the pre-specified keratinocyte score (`KRT1`, `KRT5`, `KRT10`, `KRT14`, `SFN`, `DSP`, `KRT6A`) and melanocyte score (`MLANA`, `PMEL`, `TYR`, `DCT`, `SOX10`, `MITF`):

| Cohort | vs keratinocyte score | vs melanocyte score |
|---|---|---|
| Hugo | ρ = 0.09, p = 0.67 | ρ = −0.51, p = 0.008 |
| Riaz | ρ = 0.48, p = 0.005 | ρ = −0.15, p = 0.41 |
| MGH | ρ = 0.79, p = 0.021 | ρ = 0.02, p = 0.96 |
| TCGA-SKCM (n=471) | ρ = 0.49, p = 8×10⁻³⁰ | ρ = −0.27, p = 3×10⁻⁹ |

In TCGA, `CLDN4` also correlates with `EPCAM` (ρ = 0.50), `SFN` (0.50), `KRT5` (0.40), `CLDN3` (0.38) and anti-correlates with `MLANA` (−0.30) and `MITF` (−0.26). The high-`CLDN4` tail in a melanoma biopsy is the expected signature of contaminating epidermis / adnexa, not of a melanocytic program.

Keratinocyte-adjusted logistic (Hugo + Riaz; MGH did not converge): meta OR 0.88 (0.41–1.89), p = 0.74. There is no residual tumor-intrinsic association to interpret.

## 8. Exploratory overall survival (Hugo only)

n = 25, 11 deaths.

- Continuous `z(CLDN4)`: HR 1.10 (0.63–1.92), p = 0.74
- Median split, high vs low: HR 0.35 (0.09–1.33), p = 0.12; log-rank p = 0.11

The two models point in opposite directions (continuous slightly worse, split slightly better). That is a small-n, low-event, skewed-gene artifact. It is not a survival association. OS is not available for Riaz or MGH in the GEO metadata used here.

## 9. Sensitivities (all of them, including the ones that look slightly less null)

| Analysis | What changed | Meta OR (95% CI) | p |
|---|---|---|---:|
| Primary | — | 0.96 (0.55–1.69) | 0.90 |
| S1 | Riaz `SD` coded as NR (10 R / 39 NR in Riaz) | 0.93 (0.53–1.62) | 0.79 |
| S2 | Median-split `CLDN4` high vs low | 1.55 (0.88–2.73) | 0.13 |
| S3 | Riaz from UQ-normalized log-CPM instead of author FPKM | 0.99 (0.57–1.71) | 0.97 |
| S4 | MGH all baseline episodes, still one row per patient | 0.99 (0.56–1.74) | 0.97 |
| S5 | One-stage IPD + cohort FE | 0.95 (0.56–1.60) | 0.85 |

S2 is the least null and is still not significant. The MGH median-split logistic is numerically unstable (complete separation; OR ~10⁴, SE ~6000) and receives ~zero weight in the IVW meta, so S2 is effectively Hugo + Riaz. Fisher exact p-values in those two cohorts are 0.70 and 0.26. S2 is not a finding, and it is not substituted for the primary analysis.

S4 does not add patients: averaging later baseline episodes still collapses to 8 MGH patients.

## 10. What this does *not* show

- It does not show that `CLDN4` is a melanoma ICI biomarker.
- It does not show that `CLDN4` is *not* a melanoma ICI biomarker. The experiment cannot see `CD274`.
- It does not license a claim about CLDN4 protein, IHC, or spatial expression. This is bulk mRNA.
- It does not license a claim about other histologies. Claudin-4 is an epithelial gene; melanoma was a hard place to look.
- It does not evaluate combination ICI, adjuvant ICI, or controlled-access series.

## 11. Files

| Path | Content |
|---|---|
| `ANALYSIS_PLAN.md` | Pre-specified plan (committed first) |
| `summary.json` | Machine-readable primary numbers and verdict |
| `tables/` | Inventory, detectability, per-cohort tests, meta, controls, genome-wide, sensitivities, TCGA, Hugo OS |
| `figures/cldn4_boxplot_by_response.png` | `CLDN4` by R/NR |
| `figures/cldn4_forest.png` | Per-cohort and meta ORs |
| `figures/genomewide_volcano.png` | Genome-wide meta; `CLDN4` in red, ICI controls in blue |
| `figures/cldn4_epithelial_correlation.png` | Spearman ρ vs epithelial / melanocyte genes |
| `figures/tcga_skcm_cldn4_hist.png` | TCGA-SKCM distribution |
| `figures/hugo_os_km.png` | Exploratory KM, median split |
| `scripts/b5_skcm/` | Download, loaders, analysis |

Reproduce:

```
python3 -m venv .venv && .venv/bin/pip install -r scripts/b5_skcm/requirements.txt
.venv/bin/python scripts/b5_skcm/download_data.py
.venv/bin/python scripts/b5_skcm/run_analysis.py
```

Raw GEO / Xena files are not committed (`data/raw/` is gitignored).
