# B5 Braun: CLDN4 versus objective response

## Result

There is **no statistically supported association** between pretreatment tumor
CLDN4 RNA expression and objective response to nivolumab in this exploratory
analysis of the public Braun 2020 RCC data.

- Responders (CR/PR): n = 39, median 28.426 (IQR 27.493–29.030)
- Nonresponders (SD/PD): n = 133, median 28.814 (IQR 27.471–29.811)
- Median difference (responder − nonresponder): −0.388
  (bootstrap 95% CI −1.007 to 0.210)
- Two-sided Mann–Whitney U = 2141, p = 0.098
- Common-language effect probability = 0.413 (the probability that a randomly
  selected responder has higher CLDN4 than a randomly selected nonresponder,
  with ties split evenly)

The point estimate trends toward *lower*, not higher, CLDN4 in responders, but
the confidence interval includes zero and the pooled p-value is above 0.05.
Cohort-specific results are heterogeneous: CM-009 p = 0.590 (3 responders),
CM-010 p = 0.066 (11 responders), and CM-025 p = 0.491 (25 responders).
These data do not establish CLDN4 as a response biomarker.

## Cohort and endpoint

The source is Supplementary Tables 1–4 from Braun et al., *Nature Medicine*
2020, DOI [10.1038/s41591-020-0839-y](https://doi.org/10.1038/s41591-020-0839-y).
The analysis joins the normalized CLDN4 row in `S4A_RNA_Expression` to clinical
records in `S1_Clinical_and_Immune_Data` by `RNA_ID`.

The analysis set includes all RNA-profiled, nivolumab-treated subjects from
CheckMate 009, 010, and 025 with evaluable objective response (n = 172). It
excludes 130 everolimus-treated subjects and 9 nivolumab-treated subjects whose
response is `NE`. Responders are source labels `CR`, `PR`, or the pooled
CheckMate-025 label `CRPR`; nonresponders are `SD` or `PD`. Each included
subject contributes one sample.

The source matrix contains upper-quartile-normalized, log2-transformed,
ComBat-corrected expression with a positive offset, so its values should not be
interpreted as raw TPM. No new normalization or response-driven cutpoint was
applied here.

## Files

- `analysis_samples.csv`: auditable sample-level analysis set
- `results.json`: source manifest, inclusion counts, pooled statistics, and
  descriptive cohort-specific statistics
- `CLDN4_vs_ORR.png`: box-and-jitter plot
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
the CSV, JSON, and PNG outputs. A local workbook can be supplied with
`--workbook PATH`.

## Limitations

This is a post hoc, unadjusted, single-gene association test in a sequenced
subset of three trials. The pooled test does not adjust for trial-level
differences, and the cohort-specific responder groups are small. The bootstrap
interval is descriptive, not a replacement for external validation. Response
association in nivolumab-treated patients alone cannot distinguish a
predictive treatment interaction from a general prognostic association; that
would require a prespecified treatment-by-CLDN4 comparison against the
everolimus arm. No multiplicity correction was applied because only the
requested CLDN4 hypothesis was tested.
