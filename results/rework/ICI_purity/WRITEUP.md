# Rework: ESTIMATE-purity residualized TACSTD2 / CLDN4 / TJ score vs ICI DCB/ORR

**Slice:** `scripts/rework/ICI_purity/` + `results/rework/ICI_purity/` only.
**Why this exists:** the open lung ICI bulk test of **single-gene** TACSTD2 and CLDN4 vs response was **not significant**. This rework asks whether that null is a purity confound (epithelial genes diluted by infiltrate) and whether a **TJ score** (not a single gene) changes the answer.

**Bottom line:** all five requested series have usable bulk RNA. After residualizing TACSTD2, CLDN4, and a locked 5-gene TJ score on ESTIMATEScore, **no DCB/ORR association reaches p < 0.05**. BH q = 1.0 across the 15 primary residual tests. The original single-gene NS result is **not rescued** by purity residualization or by switching to a TJ score.

Reproduce:

```bash
python3 scripts/rework/ICI_purity/analyze.py
```

---

## Requested cohorts (all reported)

| Dataset | Bulk? | Genes | Endpoint used | n (pos / neg) | Status |
|---|---|---|---|---|---|
| **GSE126044** | yes, counts → log2CPM | TACSTD2, CLDN4, TJ_5 5/5 | author responder / non-responder (ORR-like) | 16 (5 / 11) | analyzed |
| **GSE135222** | yes, TPM → log2(TPM+1) | all present (Ensembl) | **DCB** = PFS ≥ 183 d | 27 (7 / 20) | analyzed |
| **GSE166449** | yes, TPM → log2(TPM+1) | all present | author responder / non-responder (ORR-like) | 22 (7 / 15) | analyzed |
| **GSE190265** | yes, TPM → log2(TPM+1) | all present | **DCB** = PFS ≥ 6 months | 43 (14 / 29) | analyzed |
| **GSE207422** | yes, bulk log2TPM (scRNA ignored) | all present | **ORR** = RECIST CR/PR vs SD | 24 (17 / 7) | analyzed |

No requested series was skipped. GSE207422 also has MPR/NMPR (9 / 15); that native pathologic endpoint is reported as secondary and matches the prior raw NS numbers.

---

## Methods

### ESTIMATE purity

Python port of `estimateScore()` from the ESTIMATE R package v1.0.13 (Yoshihara et al., *Nat Commun* 2013):

- Per-sample rank normalization, scaled to `10000 × rank / Ngenes`.
- ssGSEA-like enrichment with weight exponent **0.25** on the official `StromalSignature` and `ImmuneSignature` (141 genes each; `SI_geneset.gmt`).
- `ESTIMATEScore = StromalScore + ImmuneScore`.
- `TumorPurity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)`.

The Affymetrix-calibrated purity cosine is **out of bounds on most RNA-seq samples** here (6/16, 27/27, 12/22, 42/43, 24/24). Primary residualization therefore uses **ESTIMATEScore**, which is defined for every sample and is a strictly monotone transform of TumorPurity wherever the formula is valid (same rank adjustment as claim A2). TumorPurity residuals are in `response_tests.csv` as a sensitivity and are not the primary test.

GMT overlap after HGNC alias/previous-symbol mapping: stromal 139–141 / 141; immune 133–141 / 141 (GSE190265 immune 133/141 is the lowest; still usable).

### Residualization

For each feature *y* (log-expression or TJ score):

```
y ~ 1 + ESTIMATEScore     (OLS)
residual = y − ŷ
```

This is the quantity tested against DCB/ORR. Raw (unadjusted) tests are reported alongside so the rework can be compared to the original NS slice.

### TJ score (not a single gene)

**Primary composite `TJ_5`:** mean of per-gene z-scores of **CLDN1, CLDN4, CLDN7, F11R, PARD3** (the locked TJ programme used in the TCGA/B3–B4 slices). All five genes are present in every cohort.

**Sensitivity `TJ_structural`:** 15-gene structural set (CLDN1/3/4/7, OCLN, TJP1/2/3, F11R, JAM2/3, MARVELD2/3, CGN, CGNL1). All 15 found in every cohort. Full tests are in `response_tests.csv`; they do not change the verdict.

### Endpoints (pre-specified, not tuned)

| Cohort | Rule | Ambiguous cases |
|---|---|---|
| GSE126044 | GEO `patient response` = responder vs non-responder | none |
| GSE135222 | DCB if `pfs.time` ≥ 183 days | 0 patients censored before 183 d |
| GSE166449 | `responder` / `nonresponder` in GEO title | none |
| GSE190265 | DCB if `time_PFS` ≥ 6 months (`samples_info_France3`) | 0 early-censored; PFS in months |
| GSE207422 | ORR if RECIST ∈ {CR, PR}; non-ORR if SD/PD | 4 metadata rows lack RECIST (excluded); MPR stored separately |

GSE190265 supplementary TPM + PFS have **43** samples; the series matrix lists 34 GSM titles, of which 26 string-match the supplementary IDs. Analysis uses the official supplementary 43 (complete TPM+PFS). This is larger than the GEO “35 patients” blurb; we do not drop samples to match the blurb.

### Statistics

- Two-sided Mann–Whitney U; AUC = P(value_benefit > value_no-benefit); rank-biserial = 2·AUC − 1.
- AUC 95% CI: 2000× stratified bootstrap, seed 20260816.
- BH-FDR across the 15 primary residual tests (5 cohorts × TACSTD2, CLDN4, TJ_5) and within each cohort (3 tests).
- Cohorts are **not pooled** (different assays and endpoints).

CD8 score (mean z of CD8A/CD8B) is a **control**, not a primary target. Raw CD8 is expected to track response in immune-hot tumors; residualizing it on ESTIMATEScore (which is mostly immune + stroma) removes that signal by construction.

---

## Primary results (residual on ESTIMATEScore)

| Dataset | Feature | Endpoint | n | AUC (benefit > no-benefit) | 95% CI | p | q (within) | q (all 15) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | TACSTD2 | author ORR | 5 / 11 | 0.47 | 0.13–0.84 | 0.91 | 1.00 | 1.00 |
| GSE126044 | CLDN4 | author ORR | 5 / 11 | 0.42 | 0.11–0.75 | 0.66 | 1.00 | 1.00 |
| GSE126044 | TJ_5 | author ORR | 5 / 11 | 0.49 | 0.15–0.82 | 1.00 | 1.00 | 1.00 |
| GSE135222 | TACSTD2 | DCB | 7 / 20 | 0.41 | 0.20–0.63 | 0.50 | 0.57 | 1.00 |
| GSE135222 | CLDN4 | DCB | 7 / 20 | 0.58 | 0.29–0.84 | 0.57 | 0.57 | 1.00 |
| GSE135222 | TJ_5 | DCB | 7 / 20 | 0.38 | 0.17–0.60 | 0.37 | 0.57 | 1.00 |
| GSE166449 | TACSTD2 | author ORR | 7 / 15 | 0.55 | 0.27–0.82 | 0.73 | 1.00 | 1.00 |
| GSE166449 | CLDN4 | author ORR | 7 / 15 | 0.53 | 0.25–0.81 | 0.84 | 1.00 | 1.00 |
| GSE166449 | TJ_5 | author ORR | 7 / 15 | 0.50 | 0.23–0.75 | 1.00 | 1.00 | 1.00 |
| GSE190265 | TACSTD2 | DCB | 14 / 29 | 0.47 | 0.29–0.66 | 0.79 | 0.79 | 1.00 |
| GSE190265 | CLDN4 | DCB | 14 / 29 | 0.47 | 0.29–0.66 | 0.77 | 0.79 | 1.00 |
| GSE190265 | TJ_5 | DCB | 14 / 29 | 0.46 | 0.29–0.64 | 0.71 | 0.79 | 1.00 |
| GSE207422 | TACSTD2 | RECIST ORR | 17 / 7 | 0.44 | 0.15–0.71 | 0.66 | 0.66 | 1.00 |
| GSE207422 | CLDN4 | RECIST ORR | 17 / 7 | 0.29 | 0.02–0.60 | 0.13 | 0.39 | 1.00 |
| GSE207422 | TJ_5 | RECIST ORR | 17 / 7 | 0.35 | 0.07–0.67 | 0.29 | 0.43 | 1.00 |

Every 95% CI includes 0.5. The smallest residual p is GSE207422 CLDN4 vs RECIST (p = 0.13, AUC = 0.29: residual CLDN4 tends **lower** in ORR+). That is a non-significant trend, not a hit.

### Raw (unadjusted) — confirms the original NS mismatch

| Dataset | TACSTD2 AUC (p) | CLDN4 AUC (p) | TJ_5 AUC (p) |
|---|---|---|---|
| GSE126044 | 0.36 (0.44) | 0.24 (0.11) | 0.33 (0.32) |
| GSE135222 | 0.43 (0.61) | 0.56 (0.69) | 0.44 (0.65) |
| GSE166449 | 0.62 (0.41) | 0.51 (0.95) | 0.48 (0.89) |
| GSE190265 | 0.46 (0.67) | 0.49 (0.93) | 0.47 (0.77) |
| GSE207422 ORR | 0.55 (0.76) | 0.36 (0.32) | 0.46 (0.80) |
| GSE207422 MPR (secondary) | 0.38 (0.34) | 0.36 (0.26) | 0.33 (0.19) |

GSE126044 and GSE166449 raw TACSTD2/CLDN4 p-values match the prior `fable_ici_bulk` table to reported precision (0.4409 / 0.1149 and 0.4069 / 0.9452). GSE207422 raw MPR matches that table (0.340 / 0.257). Residualizing MPR on ESTIMATEScore moves those AUCs to ~0.47–0.51 (p > 0.85).

### Do the genes even track purity?

Spearman of the feature vs ESTIMATEScore (context, not a response test):

- TACSTD2 vs ESTIMATEScore is weak except GSE207422 (ρ = −0.53, p = 0.007).
- TJ_5 vs ESTIMATEScore is negative in GSE135222 (ρ = −0.51, p = 0.007) and GSE207422 (ρ = −0.51, p = 0.011).
- CD8 vs ESTIMATEScore is strongly positive in every cohort (ρ = 0.45 to 0.80), as expected.

Purity confounding of TACSTD2/CLDN4 is **real in some cohorts and modest in others**. Removing it does not create a response association.

### Control: raw CD8 still sees some endpoints

Raw CD8 vs benefit: GSE126044 AUC = 1.00, p = 4.6×10⁻⁴; GSE190265 AUC = 0.71, p = 0.028; GSE135222 / GSE166449 AUC ~0.72–0.74, p = 0.08–0.09; GSE207422 RECIST null. So several binary labels are not inert. After residualizing CD8 on ESTIMATEScore the CD8–response tests become NS (the immune content was the covariate). That is a check on the residualization, not evidence against CD8 as a bulk ICI correlate.

---

## Honest caveats

1. **Small n.** 5 vs 11, 7 vs 20, 7 vs 15, 14 vs 29, 17 vs 7. Only large effects (AUC ≳ 0.85 in the smallest sets) would be detectable. NS here is **inconclusive**, not proof of no association.
2. **Endpoints are heterogeneous.** Author ORR, derived DCB, RECIST ORR, and MPR are not the same estimand. They are not meta-analyzed.
3. **ESTIMATE on RNA-seq.** The purity cosine is Affymetrix-calibrated and mostly OOB; we residualize on ESTIMATEScore. GMT symbols are 2013-era; HGNC aliases recover nearly all genes.
4. **GSE190265 n.** Supplementary France3 TPM/PFS = 43; series matrix = 34 GSM. We analyze 43 and state the mismatch.
5. **GSE166449 TPM values are low** for TACSTD2 (median log2(TPM+1) ~2). Treated as the GEO file claims (raw TPM), same as the prior slice.
6. **Bulk only.** Tumor-intrinsic TROP2/CLDN4 protein or malignant-cell scRNA is not tested.
7. **No claim that high TROP2/TJ predicts ICI resistance in these public bulks.** The data do not support that statement after or before purity residualization.

---

## Files

| Path | Contents |
|---|---|
| `tables/cohort_inventory.csv` | every requested series, n, endpoint, ESTIMATE overlap, gene presence |
| `tables/response_tests.csv` | raw + residual_ESTIMATEScore + residual_TumorPurity; TJ_structural; CD8; GSE207422 MPR |
| `tables/feature_vs_purity.csv` | Spearman of each feature vs ESTIMATEScore / TumorPurity / ImmuneScore |
| `tables/sample_scores.csv` | per-sample expression, TJ scores, ESTIMATE, residuals, labels |
| `tables/gene_coverage.csv` | which TJ/target/CD8 genes mapped |
| `tables/download_manifest.csv` | URL, bytes, sha256 |
| `tables/key_stats.json` | machine-readable summary |
| `figures/*_raw_vs_residual_box.png` | raw vs residual boxplots |
| `figures/*_feature_vs_ESTIMATEScore.png` | feature vs purity proxy, colored by label |
| `figures/primary_residual_auc_forest.png` | 15 primary AUCs with bootstrap CIs |
| `scripts/rework/ICI_purity/analyze.py` | self-contained runner |
| `scripts/rework/ICI_purity/SI_geneset.gmt` | official ESTIMATE signatures |

Data cache: `/tmp/ici_purity_data` (not committed). Override with `ICI_PURITY_DATA` / `ICI_PURITY_OUT`.
