# Claim B5 — 11-cohort ICI meta, CLDN4-high OR = 0.42

**Verdict: not reproduced.**

The claim is that CLDN4-high tumours have an ICI response odds ratio of **0.42** in an **11-cohort** meta-analysis spanning **GBM, NSCLC, melanoma, RCC, and urothelial** cancer. I could not find a published paper that names those 11 cohorts or reports that OR. I did not invent an 11-cohort list to match the claim.

What *is* reconstructable from public processed data is a **different** set: **12 named PD-1/PD-L1 cohorts** with CLDN4 and a binary response label in NSCLC / melanoma / RCC / urothelial. **GBM = 0 cohorts.** The pooled random-effects OR in that set is **1.12 (95% CI 0.79–1.58)**, p = 0.53. The claimed 0.42 lies **outside** that interval, and the point estimate is in the **opposite** direction (CLDN4-high slightly more, not less, likely to respond).

## What was searched for provenance

PubMed and Europe PMC queries for `CLDN4` plus immunotherapy / checkpoint / ICI / PD-1 did not return a multi-cohort ICI response meta-analysis with OR ≈ 0.42. Hits were gastric-cancer CLDN4 clinicopathology, ovarian/RCC immune-risk models that mention CLDN4 as one of many genes, and unrelated antibody work. Bareche et al., *Ann Oncol* 2022 (PredictIO) is the large public ICI expression meta; it does not highlight CLDN4 as a 11-cohort OR=0.42 biomarker.

## Public cohorts that actually contain CLDN4

Source for almost all objects: ORCESTRA / PredictIO ICB TSV extract, Zenodo [10.5281/zenodo.7199344](https://doi.org/10.5281/zenodo.7199344) (Bareche et al. 2022; Mammoliti et al. 2021). IMvigor210 CLDN4 was rescued from the official Genentech counts (see below).

### Entered in the primary meta (k = 12, PD-1/PD-L1, median split, PredictIO R vs NR)

| Study | Indication | n (R/NR) | OR (CLDN4-high vs low) | 95% CI | p | Citation |
|---|---|---|---|---|---|---|
| Jung | NSCLC | 27 (8/19) | 1.85 | 0.34–10.05 | 0.48 | Jung et al., *Nat Commun* 2019, PMID 31383964 |
| Fumet1 (Limagne1) | NSCLC | 69 (17/52) | 3.03 | 0.93–9.83 | 0.065 | ORCESTRA ICB_Fumet1; GEO GSE190266 |
| Fumet2 (Limagne2) | NSCLC | 22 (6/16) | 1.00 | 0.15–6.53 | 1.00 | ORCESTRA ICB_Fumet2; same GEO family (GSM5718771) |
| Gide | melanoma | 38 (19/19) | 1.23 | 0.35–4.41 | 0.75 | Gide et al., *Cancer Cell* 2019, PMID 30753825 |
| Hugo | melanoma | 27 (14/13) | 1.56 | 0.34–7.11 | 0.57 | Hugo et al., *Cell* 2016, PMID 26997480 |
| Liu | melanoma | 112 (52/60) | 1.15 | 0.55–2.43 | 0.70 | Liu et al., *Nat Med* 2019, PMID 31792460 |
| Riaz | melanoma | 30 (9/21) | 2.67 | 0.52–13.66 | 0.24 | Riaz et al., *Cell* 2017, PMID 29033130 |
| Puch | melanoma | 49 (14/35) | 2.14 | 0.59–7.68 | 0.24 | ORCESTRA ICB_Puch (PUCH melanoma) |
| Braun | RCC | 142 (39/103) | 0.81 | 0.39–1.69 | 0.57 | Braun et al., *Nat Med* 2020, PMID 32472114 |
| Shiuan | RCC | 13 (6/7) | 0.20 | 0.02–2.12 | 0.18 | Shiuan et al., *JCI Insight* 2020, PMID 32315290 |
| Mariathasan | urothelial | 235 (68/167) | 0.63 | 0.36–1.12 | 0.12 | Mariathasan et al., *Nature* 2018, PMID 29443960 |
| Snyder | urothelial | 18 (8/10) | 2.50 | 0.37–16.89 | 0.35 | Snyder et al., *PLoS Med* 2017, PMID 28542453 |

Pooled DerSimonian–Laird random effects: **OR = 1.12 (0.79–1.58)**, p = 0.53, I² = 12%, τ² = 0.05, n = 782 labelled patients. Fixed-effect OR is similar. Leave-one-out ORs stay between 0.98 and 1.31; none approach 0.42.

`response` is the PredictIO binary label (CR/PR, or SD without a 6-month PFS event = R; PD, or SD with a 6-month PFS event = NR). CLDN4-high = at or above the **within-cohort median**. Zero cells get a Haldane–Anscombe 0.5 correction (none of the primary 12 needed it).

### Named but not entered in the primary meta

See `excluded_or_unusable.tsv`. Short version:

- **GBM:** no public processed bulk expression + binary ICI response table. Zhao 2019 is SRA + “processed data on request.” GSE121810 series matrix is metadata-only. GSE154795 is immune-cell scRNA-seq. GSE264695 is neoadjuvant pharmacodynamic RNA-seq without a usable public R/NR table here.
- **Miao 2018 RCC:** CLDN4 is at the TPM floor in 32/33 RNA samples. A median split has no low arm.
- **Hwang, Jerby-Arnon, Roh:** panel/signature objects; CLDN4 is not on the panel.
- **Van Allen, Nathanson:** CTLA-4, not PD-1/PD-L1 (Nathanson OR = 1.00 if forced in).
- **Kim 2018 gastric:** CLDN4 present, median-split OR = 0.12 (p = 0.010) — the only nominally significant cohort, and **outside** the claimed indications.
- **Padron pancreas:** out of scope; 3 responders.

I did not pad the list to 11 with gastric, CTLA-4, or unusable GBM.

## IMvigor210 CLDN4 had to be rescued

In the ORCESTRA `ICB_Mariathasan` expression object, `ENSG00000189143` / CLDN4 is a **constant** −9.96578 in all 348 samples. GAPDH, ACTB, KRT8, EPCAM, and PTPRC in the same matrix vary, so this is not a row-shift. It is a gene-specific quantification failure. Using that vector would have silently dropped the largest urothelial ICI RNA cohort.

Replacement: official IMvigor210CoreBiologies counts (`counts(cds)`), log2(count+1), matched to PredictIO IDs as `P{ANONPT_ID}`. CSV dump: [mimifp/tfm_mUC](https://github.com/mimifp/tfm_mUC) from [research-pub.gene.com/IMvigor210CoreBiologies](http://research-pub.gene.com/IMvigor210CoreBiologies/). After replacement: 43 unique values, log2 IQR = 1.28, median 3.32 in responders vs 3.58 in non-responders (OR 0.63, p = 0.12). Documented in `imvigor210_cldn4_rescue.tsv`.

## Sensitivity (still not 0.42)

| Analysis | k | RE OR (95% CI) | p |
|---|---|---|---|
| Primary (above) | 12 | 1.12 (0.79–1.58) | 0.53 |
| Drop Braun (CLDN4 IQR < 0.5) | 11 | 1.22 (0.82–1.82) | 0.33 |
| RECIST CR/PR vs PD (SD dropped) | 9 | 0.97 (0.64–1.45) | 0.86 |
| Tertile high vs low | 12 | 0.93 (0.64–1.35) | 0.70 |
| Add Nathanson CTLA-4 | 13 | 1.08 (0.79–1.48) | 0.64 |
| Add Kim gastric | 14 | 1.06 (0.71–1.58) | 0.77 |
| Epithelial only (drop melanoma) | 7 | 1.01 (0.58–1.74) | 0.98 |

By indication (primary definition): NSCLC 2.11 (0.89–4.99), k=3; melanoma 1.44 (0.87–2.40), k=5; RCC 0.64 (0.23–1.77), k=2; urothelial 0.92 (0.28–3.02), k=2; **GBM k=0**.

Continuous logistic OR per +1 SD CLDN4 is non-significant in every primary cohort (see `primary_cohort_table.tsv`). Mann–Whitney CLDN4 in R vs NR is significant only in Fumet1 (higher CLDN4 in *responders*, p = 0.041) — again the opposite of the claim.

## Why an 11-cohort OR=0.42 should not be treated as established

1. **The 11 are not named** in any source I could find. Forcing k=11 by dropping Mariathasan (the largest UC cohort) or Shiuan (n=13) is post-hoc and still yields OR ≈ 1.2–1.3, not 0.42.
2. **GBM is not in the public processed set.** Including it would require unpublished processed Zhao data or a custom SRA realignment plus a response table that is not in GEO supplements used here.
3. **CLDN4 is a poor pan-cancer ICI feature.** It is an epithelial tight-junction gene. Melanoma cohorts contribute 256/782 primary patients; CLDN4 biology there is not the same as in urothelium or lung. Braun RCC has almost no CLDN4 dynamic range (IQR 0.09). IMvigor210 was all-zero in the harmonized object until rescued.
4. **The reconstructable estimate is null-to-opposite.** Every pre-specified pooling that I can name has a CI that includes 1 and excludes 0.42, except the tiny Shiuan cohort alone (n=13, OR=0.20, CI 0.02–2.12).

## Methods (short)

- Expression: per-study ORCESTRA TSV, CLDN4 = `ENSG00000189143`. IMvigor210 replaced as above.
- One RNA sample per PredictIO `patientid` (inner join of expression columns to metadata).
- Primary cut: within-cohort median. Sensitivities: tertiles, quartiles, RECIST CR/PR vs PD.
- OR = (high R / high NR) / (low R / low NR). Random-effects meta on logOR, DerSimonian–Laird; also inverse-variance fixed effect.
- Code: `code/03_extract_cldn4.py`, `code/03b_imvigor210_official.py`, `code/04_meta_analysis.py`.
- Raw ORCESTRA zips are not committed (see `.gitignore`). Re-download Zenodo 7199344 and the IMvigor210 CSVs to repeat.

## Outputs

- `cohort_inventory.tsv` — every ORCESTRA object inspected
- `excluded_or_unusable.tsv` — why GBM and others are out
- `cldn4_clinical.tsv` — per-sample CLDN4 + labels
- `primary_cohort_table.tsv`, `per_cohort_or.tsv`, `leave_one_out_primary.tsv`
- `meta_summaries.json`, `claim_verdict.json`
- `forest_primary_pd1_inscope_median_predictio.png` — main forest plot
- `imvigor210_cldn4_rescue.tsv`
