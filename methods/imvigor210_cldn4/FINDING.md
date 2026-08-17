# IMvigor210 urothelial analog — CLDN4 vs CD274 / response

**Analog, not lung.** Public pretreatment RNA from IMvigor210 (metastatic urothelial carcinoma, atezolizumab; Mariathasan et al., *Nature* 2018). This is **bladder / urothelial**, not NSCLC. Do not quote these numbers as a lung ICI result.

**Additive only.** Prior TACSTD2 / CLDN4 vs ORR / OS on this same package (PR 154, 207, 215) is taken as given. This folder only adds **CLDN4 vs CD274** and restates **CLDN4 vs response** on the same matrix, with CD274 vs response beside it.

Official processed CountDataSet: `IMvigor210CoreBiologies` **1.0.0** (CC BY 3.0). **No GEO / GSE accession is claimed.** Raw RNA is EGA `EGAD00001003977` (controlled; not used).

Reproduce: `python3 methods/imvigor210_cldn4/analyze.py`

---

## Honest n

| item | n | rule |
|---|---:|---|
| RNA libraries in `cds` | **348** | CountDataSet columns |
| unique `ANONPT_ID` | **347** | package patient id |
| patients after sizeFactor dedup | **347** | one extra library for patient 10285 (both NE); keep larger sizeFactor |
| CLDN4 / CD274 / CD8A / CXCL9 finite | **347** | one row each in `featureData` (Entrez 1364 / 29126 / 925 / 4283) |
| ORR-evaluable (CR/PR vs SD/PD) | **298** | package `binaryResponse`; NE dropped |
| CR/PR | **68** | 25 CR + 43 PR |
| SD/PD | **230** | 63 SD + 167 PD |
| NE after dedup | **49** | excluded from ORR |
| Tissue == bladder | **194** | package `Tissue`; other sites are urothelial mets |
| bladder + ORR-evaluable | **168** | sensitivity, not primary |
| Tissue == lung | **10** | urothelial **met** in lung, not a lung primary |
| lung-primary NSCLC RNA | **0** | not this cohort |

Primary tests use **n=347** (CLDN4 vs CD274) and **n=298** (response). Do not write n=348 for patient-level tests. Do not write n=194 as the primary n.

---

## One-row table

| dataset | analog | n patients | n ORR (CR/PR vs SD/PD) | CLDN4–CD274 ρ (p) | CLDN4–CD274 \| CD8A ρ (p) | CLDN4 vs ORR AUC (p) | CLDN4 median-split Fisher OR (p) | CD274 vs ORR AUC (p) |
|---|---|---:|---:|---|---|---|---|---|
| IMvigor210 mUC atezolizumab | **urothelial, not lung** | 347 | **298** (68 vs 230) | **−0.169 (0.0016)** | −0.087 (0.13) | **0.547 (0.24)** | **1.47 (0.21)** | 0.566 (0.099) |

Full numeric row: `tables/one_row.tsv`. Axis tests: `tables/spearman.tsv`, `tables/orr.tsv`.

---

## CLDN4 vs CD274

Scale: `log2(TPM+1)` from package gene lengths. Spearman; Fisher-z 95% CI. Partial = Pearson of average-rank residuals after `CD8A` (df = n−3; Fisher-z SE = 1/√(n−4)).

| subset | n | ρ | 95% CI | p | partial ρ \| CD8A | partial p |
|---|---:|---:|---|---:|---:|---:|
| all patients | **347** | **−0.169** | −0.269 to −0.064 | 0.0016 | −0.096 | 0.074 |
| ORR-evaluable | **298** | **−0.168** | −0.277 to −0.056 | 0.0036 | −0.087 | 0.13 |
| Tissue == bladder | 194 | −0.115 | −0.252 to +0.026 | 0.11 | −0.030 | 0.68 |
| bladder + ORR | 168 | −0.140 | −0.285 to +0.012 | 0.071 | −0.040 | 0.61 |

DESeq size-factor `log2(count/sf+1)` on the same 347 patients is ρ = −0.198 (p = 2.1×10⁻⁴). Same sign.

CLDN4 is weakly **CD274-low** on the full urothelial matrix. After residualising on CD8A the residual is **null**. Do not call CLDN4 a PD-L1-high or PD-L1-low marker that is independent of T-cell transcript.

Companion ranks (n=347):

| pair | ρ | p |
|---|---:|---:|
| CLDN4 vs CD8A | −0.146 | 0.0066 |
| CLDN4 vs CXCL9 | −0.108 | 0.044 |
| CD274 vs CD8A | +0.678 | 5.2×10⁻⁴⁸ |
| CD274 vs CXCL9 | +0.697 | 7.9×10⁻⁵² |
| CD8A vs CXCL9 (positive control) | +0.833 | 7.5×10⁻⁹¹ |

CD274 tracks CD8A / CXCL9. CLDN4 does not.

---

## CLDN4 vs response (and CD274 beside it)

Author `binaryResponse`: CR/PR vs SD/PD. NE excluded. Median split is prespecified on `log2(TPM+1)`, not optimized. Fisher OR is high vs low. Woolf logit 95% CI. AUC: CR/PR higher. Bootstrap AUC CI, 2,000 resamples, seed `20260817`.

| gene | n | median CR/PR | median SD/PD | AUC (95% CI) | MWU p | Fisher OR (95% CI) | 2×2 (R/n high vs low) | Fisher p |
|---|---:|---:|---:|---|---:|---|---|---:|
| **CLDN4** | **298** | 4.73 | 4.41 | **0.547 (0.47–0.63)** | **0.24** | **1.47 (0.85–2.53)** | 39/149 vs 29/149 | **0.21** |
| CD274 | 298 | 1.49 | 1.19 | 0.566 (0.49–0.64) | 0.099 | 1.49 (0.86–2.58) | 39/148 vs 29/150 | 0.17 |
| CD8A (control) | 298 | 3.07 | 2.77 | 0.585 (0.50–0.66) | **0.033** | 1.28 (0.74–2.20) | 37/148 vs 31/150 | 0.41 |
| CXCL9 (control) | 298 | 4.62 | 3.39 | 0.645 (0.56–0.72) | **2.8×10⁻⁴** | 1.94 (1.11–3.39) | 43/151 vs 25/147 | **0.019** |

CLDN4-high ORR is 26.2% vs 19.5% in CLDN4-low. Point estimate goes **slightly the other way** from a CLDN4-high → worse-ORR claim. The interval includes 1. Continuous logistic OR per 1 SD = 1.15 (0.87–1.52), p = 0.34. Q4 vs Q1 OR = 1.43 (0.68–3.03), n=75 vs 75, p = 0.45.

CD274 vs ORR is also **null** (MWU p = 0.099). CXCL9 and CD8A move in the expected direction on the same 298 rows, so the endpoint is not inert.

### Tissue == bladder only (sensitivity)

| test | n | result |
|---|---:|---|
| CLDN4 vs CD274 | 194 | ρ = −0.115, p = 0.11 |
| CLDN4 vs ORR, median Fisher | **168** (42 vs 126) | OR **0.97** (0.48–1.95), p = **1.0**; AUC 0.494, p = 0.91 |
| CD274 vs ORR, MWU | 168 | p = 0.11 |
| CXCL9 vs ORR, MWU | 168 | p = 6.8×10⁻⁴ |

Restricting to bladder-site RNA does not create a CLDN4 response signal. It **removes** the weak CLDN4–CD274 correlation (CI includes 0).

---

## What this is not

- **Not lung / NSCLC.** Ten `Tissue == lung` rows are urothelial metastases. They stay in the primary 347 / 298. They are not a lung ICI cohort.
- Not a predictive atezolizumab-vs-chemo interaction. Single-arm PD-L1 blockade.
- Not a re-audit of TACSTD2 vs ORR / OS (already null on this package).
- No GEO accession. No invented IMvigor210 GSE.
- ORCESTRA / PredictIO IMvigor210 CLDN4 is a known floor artifact and is **not** used.

---

## Data

| item | value |
|---|---|
| Package | IMvigor210CoreBiologies **1.0.0** |
| URL | `http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz` |
| SHA256 | `cfdd3176d7b34de5b04fb9416bfd2b20fa4b6e238aaad5f20b048a34329ea178` |
| Object | `cds` CountDataSet, 31,286 genes × 348 RNA samples |
| Trials | NCT02108652, NCT02951767 |
| Paper | Mariathasan et al., *Nature* 2018, DOI 10.1038/nature25501 |
| Raw RNA | EGAS00001002556 / EGAD00001003977 (controlled; not downloaded) |
| CLDN4 | Entrez 1364, length 1831 nt, TPM 0.14–136 |
| CD274 | Entrez 29126, length 8095 nt, TPM 0.040–42 |

---

## Files

- `FINDING.md` — this page
- `analyze.py` / `extract_cds.R` — download official tarball, dump four genes + clinical, score
- `tables/one_row.tsv`, `spearman.tsv`, `orr.tsv`, `coverage.tsv`, `tissue_counts.tsv`, `sample_level.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd274.png`, `fig2_cldn4_cd274_by_response.png`, `fig3_orr_forest.png`

```bash
python3 methods/imvigor210_cldn4/analyze.py
```
