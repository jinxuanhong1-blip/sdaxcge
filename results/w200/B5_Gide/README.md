# B5 — Gide 2019 melanoma ICI: CLDN4 vs response

Public pretreatment RNA from **Gide et al., Cancer Cell 2019** (anti-PD-1 ±
anti-CTLA-4, metastatic melanoma). Question: is **CLDN4-high** associated with
objective response?

**Honest headline: no. The locked primary test is null.**

## OR / n / p (primary, prespecified)

| Endpoint | Contrast | Effect | 95% CI | n | p |
|---|---|---|---|---|---|
| **ORR** (CR/PR vs SD/PD) | CLDN4-high vs low (median) | **Fisher OR 0.503** | 0.145–1.750 | **41** (8/21 vs 11/20) | **0.354** |

2×2 ORR (all 41 pretreatment anti-PD-1 patients are RECIST-evaluable):

| CLDN4 | SD/PD | CR/PR | ORR |
|---|---|---|---|
| low | 9 | 11 | 55.0% |
| high | 13 | 8 | 38.1% |

Direction is a modest **lower** ORR in CLDN4-high. The interval includes the
null. Continuous CLDN4 (per SD of `log2(TPM+1)`) is also null: logistic OR
0.719 (0.365–1.414), n=41, p=0.339. Rank AUC if higher CLDN4 predicted
response was 0.366 (bootstrap 95% CI 0.201–0.542), two-sided Mann–Whitney
p=0.140.

Do not read this as evidence that CLDN4-high predicts ICI benefit or
resistance in melanoma.

## What was opened

The authors deposited raw RNA as ENA/SRA **PRJEB23709** (91 runs). There is no
GEO series matrix. This slice uses two independent public processed objects
and requires them to name the same 41 pretreatment anti-PD-1 patients.

| Item | Value | How verified |
|---|---|---|
| Paper | Gide et al., *Cancer Cell* 2019 | DOI `10.1016/j.ccell.2019.01.003`, PMID 30753825 |
| Raw RNA | **PRJEB23709** (91 Illumina RNA-seq runs) | ENA / NCBI BioProject |
| Primary matrix | Liulab `Gide.rda` TPM, gene **CLDN4** once | SHA256 `170a6079505310ddc7cd88a1283bb415331d89f5a144675fcc96d1602e83f9a5` |
| Cross-check object | ORCESTRA **ICB_Gide** v1.0, Zenodo `10.5281/zenodo.7332096` | 41 patients; CLDN4 = `ENSG00000189143.10` |
| ORCESTRA RNA | Kallisto / Gencode v40, log2 TPM | `expr_gene_tpm` assay |
| Overlap | 41/41 pretreatment pembrolizumab or nivolumab biopsies | RECIST strings match exactly |
| Expression concordance | Spearman ρ = 0.910 | Liulab `log2(TPM+1)` vs ORCESTRA log2 TPM |

Liulab is a third-party public TPM matrix, not an author supplement. ORCESTRA
is an independent FASTQ requantification. Raw FASTQ was not re-aligned here.

## Prespecified analysis

1. **Cohort.** Pretreatment (`PRE`) anti-PD-1 monotherapy only (pembrolizumab
   or nivolumab). n=41, one biopsy per patient. This is the open ICB_Gide
   object. Combination ipilimumab + PD-1 (n=32 PRE) is sensitivity, not
   primary.
2. **Exposure.** `log2(Liulab TPM + 1)`. High = ≥ median (0.0887). Median is
   prespecified, not optimized. CLDN4 is low and zero-inflated (13/41 TPM = 0;
   median TPM = 0.063).
3. **ORR.** RECIST CR/PR vs SD/PD. Fisher exact OR is primary; logistic OR per
   SD is reported beside it.
4. **Survival (secondary).** ORCESTRA PFS/OS in months. Unadjusted Cox HR for
   the same median split.

## Sensitivities (do not promote these)

All keep the same direction (CLDN4-high → fewer responses) and remain
compatible with a null or very small effect once the locked primary is
respected.

| Analysis | OR or AUC | n | p | Why it is not primary |
|---|---:|---:|---:|---|
| CR/PR vs PD (SD dropped) | Fisher OR 0.436 | 35 | 0.315 | Loses 6 SD; same conclusion |
| Paper Yes/No (SD split by clinical benefit) | Fisher OR 0.404 | 41 | 0.215 | Mixed SD rule; one PR coded No |
| ORCESTRA R vs NR, Liulab expression | Fisher OR 0.340 | 38 | 0.194 | Drops 3 long-PFS SD |
| ORCESTRA R vs NR, ORCESTRA expression | Fisher OR 0.213 | 38 | 0.050 | Different pipeline + different SD coding; continuous MW still p=0.114 |
| All PRE (PD-1 + combo) median-split OR | Fisher OR 0.480 | 73 | 0.160 | Mixes two regimens |
| All PRE continuous MW | AUC 0.361 | 73 | 0.038 | Not the locked OR test; unadjusted |
| Combo PRE only | Fisher OR 0.429 | 32 | 0.458 | Small and not the open ICB object |
| On-treatment EDT | Fisher OR 0.250 | 18 | 0.347 | Post-baseline; cannot establish prediction |
| Detected vs undetected CLDN4 | Fisher OR 0.404 | 41 | 0.313 | 13 zeros; still null |

The ORCESTRA-expression median split (p=0.050), the pooled-PRE Mann–Whitney
(p=0.038), and pooled-PRE detected-vs-undetected Fisher (OR 0.329, n=73,
p=0.047) are the only nominal p<0.05 cells. All are cutoff- or
cohort-sensitive, were not locked, and do not survive even a two-test view
with the primary. They are not findings.

Secondary survival on the 41-patient primary set:

| Endpoint | Contrast | Effect | 95% CI | n (events) | p |
|---|---|---|---|---|---|
| PFS | CLDN4-high vs low | Cox HR 2.06 | 0.96–4.40 | 41 (29) | 0.062 |
| PFS | per SD | Cox HR 1.08 | 0.78–1.49 | 41 (29) | 0.631 |
| OS | CLDN4-high vs low | Cox HR 1.25 | 0.56–2.80 | 41 (24) | 0.583 |
| OS | per SD | Cox HR 0.99 | 0.67–1.47 | 41 (24) | 0.963 |

Median-split PFS is another near-null. The better-powered continuous model is
centered on 1.

## Limitations (read these)

- Single-arm / no randomized control. Even a real association would be
  prognostic, not predictive of PD-1 benefit versus chemotherapy or targeted
  therapy.
- n=41 is small. An OR of 0.50 with CI 0.15–1.75 rules out only large
  effects.
- CLDN4 RNA is barely expressed in this melanoma cohort. A median split on a
  zero-inflated epithelial tight-junction gene is a weak exposure.
- Liulab TPM and ORCESTRA Kallisto TPM agree in rank (ρ=0.91) but are not
  interchangeable; one pipeline-specific cutoff p-value is not a replication.
- Public metadata do not support a credible adjusted model (purity, biopsy
  site completeness, LDH, prior ipilimumab on the PD-1-only subset).
- Combination-therapy and on-treatment biopsies are reported only as
  sensitivity.
- This is one marker in one cohort. It does not test a lung/TROP2 claim and
  does not validate CLDN4 as an ICI biomarker.

## Reproduce

```bash
python3 -m pip install -r results/w200/B5_Gide/requirements.txt
python3 results/w200/B5_Gide/analyze.py
```

The script downloads Liulab `Gide.rda` into `raw/` when absent, checks the
SHA-256, joins the committed ORCESTRA CLDN4 extract, and regenerates tables
and figures.

Outputs:

- `tables/headline_OR_n_p.tsv`: locked OR / n / p
- `tables/orr_2x2_primary.tsv`: primary 2×2
- `tables/summary.csv`: every binary contrast
- `tables/survival.csv`: PFS/OS Cox
- `tables/sample_level.csv`: public sample-level values
- `figures/cldn4_vs_response.png` / `.svg`
- `figures/orr_bar_median.png` / `.svg`
- `data/provenance.tsv`: source URLs and checksums

## References

- Gide TN, Quek C, Menzies AM, et al. Distinct Immune Cell Populations Define
  Response to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4 Combined
  Therapy. *Cancer Cell*. 2019;35:238–255.e6.
  [doi:10.1016/j.ccell.2019.01.003](https://doi.org/10.1016/j.ccell.2019.01.003)
- ENA/SRA BioProject
  [PRJEB23709](https://www.ebi.ac.uk/ena/browser/view/PRJEB23709)
- ORCESTRA ICB_Gide
  [10.5281/zenodo.7332096](https://doi.org/10.5281/zenodo.7332096)
- Liulab processed object
  [PD1highCD8Tscore `Gide.rda`](https://github.com/Liulab/PD1highCD8Tscore)
