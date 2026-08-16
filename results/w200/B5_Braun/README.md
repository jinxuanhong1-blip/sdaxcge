# B5 Braun: open RCC ICI, CLDN4 versus response

## OR / n / p

| Field | Value |
|---|---|
| Cohort | Open RCC ICI: Braun 2020 CheckMate 009/010/025 nivolumab RNA |
| Endpoint | Objective response, CR/PR vs SD/PD |
| Exposure | CLDN4 high vs low (pre-specified median split) |
| **OR** | **0.72 (95% CI 0.35–1.47)** |
| **n** | **172** (86 high / 86 low) |
| **p** | **0.47** (two-sided Fisher exact) |

2×2 counts: high 17/86 responders; low 22/86 responders.

This is a **null** result. CLDN4-high tumors do not have higher odds of objective response. The point estimate is in the opposite direction (OR < 1) and the interval includes 1.

## Continuous check (not the OR/n/p headline)

The same 172 samples, without a cutpoint:

- Responders n = 39, median 28.426 (IQR 27.493–29.030)
- Nonresponders n = 133, median 28.814 (IQR 27.471–29.811)
- Median difference (responder − nonresponder): −0.388 (bootstrap 95% CI −1.007 to 0.210)
- Two-sided Mann–Whitney p = 0.098
- Logistic OR per 1 SD higher CLDN4: 0.75 (95% CI 0.54–1.05), Wald p = 0.095

The continuous test is more powerful than the median split and is still not significant. Cohort-specific Mann–Whitney p-values remain descriptive: CM-009 p = 0.590 (3 responders), CM-010 p = 0.066 (11 responders), CM-025 p = 0.491 (25 responders).

## Cohort and endpoint

The source is Supplementary Tables 1–4 from Braun et al., *Nature Medicine*
2020, DOI [10.1038/s41591-020-0839-y](https://doi.org/10.1038/s41591-020-0839-y).
This is the main **open** RCC ICI RNA/response resource (CheckMate 009, 010, and 025).
The analysis joins the normalized CLDN4 row in `S4A_RNA_Expression` to clinical
records in `S1_Clinical_and_Immune_Data` by `RNA_ID`.

The analysis set includes all RNA-profiled, nivolumab-treated subjects with
evaluable objective response (n = 172). It excludes 130 everolimus-treated
subjects and 9 nivolumab-treated subjects whose response is `NE`. Responders
are source labels `CR`, `PR`, or the pooled CheckMate-025 label `CRPR`;
nonresponders are `SD` or `PD`. Each included subject contributes one sample.

The source matrix contains upper-quartile-normalized, log2-transformed,
ComBat-corrected expression with a positive offset, so its values should not be
interpreted as raw TPM. The odds ratio uses a pre-specified median split of
included CLDN4 values. No response-driven cutpoint was searched.

## Files

- `or_n_p.json` / `or_n_p.tsv`: headline OR, n, and p
- `analysis_samples.csv`: auditable sample-level analysis set
- `results.json`: source manifest, inclusion counts, pooled statistics, and
  descriptive cohort-specific statistics
- `CLDN4_vs_ORR.png`: box-and-jitter plot annotated with OR/n/p
- `run_analysis.py`: complete download, validation, analysis, and plotting code
- `requirements.txt`: Python dependencies

## Reproduce

From the repository root:

```bash
python3 -m pip install -r results/w200/B5_Braun/requirements.txt
python3 results/w200/B5_Braun/run_analysis.py
```

The script downloads the 119 MB public supplementary workbook to `.cache/`,
checks its SHA-256 digest, validates sample joins and uniqueness, and recreates
the CSV, JSON, TSV, and PNG outputs. A local workbook can be supplied with
`--workbook PATH`.

## Limitations

This is a post hoc, unadjusted, single-gene association test in a sequenced
subset of three trials. The median-split OR is the requested headline; it
discards rank information and is less powerful than the continuous tests, which
are also null. The pooled analysis does not adjust for trial-level differences.
Response association in nivolumab-treated patients alone cannot distinguish a
predictive treatment interaction from a general prognostic association; that
would require a prespecified treatment-by-CLDN4 comparison against the
everolimus arm. No multiplicity correction was applied because only the
requested CLDN4 hypothesis was tested.
