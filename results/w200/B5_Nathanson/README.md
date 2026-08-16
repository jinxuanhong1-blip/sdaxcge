# B5 leftover: Nathanson melanoma ICI / CLDN4

## Bottom line

An open matrix exists. Pretreatment tumor **CLDN4 is directionally higher** in
patients labeled as having clinical benefit, but the usable RNA set is only
nine biopsies and the result is **not conclusive**.

- Benefit: n = 4, median = 12.06 FPKM
- No benefit: n = 5, median = 0.259 FPKM
- Exact two-sided Mann–Whitney U = 18, p = 0.0635
- Rank AUC if higher CLDN4 predicts benefit = 0.90
  (bootstrap 95% CI 0.60–1.00)

This is leftover public melanoma ICI RNA after the Riaz and Liu B5 analyses.
No UCLA/Hugo matrix was substituted.

## Cohort and endpoint

- Source: Nathanson et al., *Cancer Immunol Res* 2017, public Cufflinks FPKM
  and cohort table from
  [hammerlab/melanoma-reanalysis](https://github.com/hammerlab/melanoma-reanalysis)
  commit `2ce35b9`.
- Treatment: CTLA-4 blockade (ipilimumab or tremelimumab) in melanoma.
- Public RNA: 24 tumors with sufficient tissue; 9 pretreatment and 15
  post-treatment.
- Gene: canonical **CLDN4** `ENSG00000189143`. A second symbol-labeled ID is
  present and is zero in every sample, so a paper-style median collapse only
  halves the values and does not change ranks or p-values.
- Primary analysis: one pretreatment biopsy per RNA-profiled patient.
- Endpoint: the deposited binary `Benefit` label (clinical benefit as released
  by the authors), not reconstructed RECIST.

Post-treatment biopsies were kept descriptive and were not mixed into the
predictive comparison.

## Results

| Analysis | Benefit / no benefit | Median FPKM | AUC (benefit higher) | 95% bootstrap CI | MW p |
|---|---:|---:|---:|---:|---:|
| Pretreatment (primary) | 4 / 5 | 12.06 / 0.259 | 0.90 | 0.60–1.00 | 0.0635 |
| Post-treatment (exploratory) | 4 / 11 | 0.297 / 0.129 | 0.75 | 0.45–0.95 | 0.177 |
| Pretreatment, symbol-median sensitivity | 4 / 5 | 6.03 / 0.129 | 0.90 | 0.60–1.00 | 0.0635 |

Leave-one-out p-values on the primary set range from 0.036 to 0.143. Dropping
the lowest-benefit sample (`SD1494`) makes the contrast nominally significant;
dropping the high expresser (`LSD0167`, 121 FPKM) or other benefit samples
moves p above 0.10. The direction is therefore interesting and fragile.

In this nine-sample pretreatment RNA subset, vital status is perfectly aligned
with the deposited benefit label (4/4 benefit alive; 5/5 no-benefit dead). An
overall-survival test would not be independent of the benefit contrast and was
not treated as a second primary endpoint.

## Important limitations

1. The primary sample is tiny (4 versus 5). A rank AUC of 0.90 with a
   confidence interval that starts at 0.60 is compatible with both a large
   effect and a chance ranking.
2. One pretreatment benefit biopsy dominates the scale (`LSD0167`).
3. Benefit is the authors' deposited binary endpoint, not a reconstructed
   RECIST call, and the paper notes discordant lesions in the broader cohort.
4. FPKM is suitable for this within-gene comparison but is not a raw-count
   differential-expression analysis.
5. This is a leftover single-marker association check. It does not test
   whether CLDN4 is specifically predictive of CTLA-4 blockade versus
   prognostic in untreated melanoma.

## Files

- `cldn4_vs_benefit.png` / `.svg`: primary and exploratory plots
- `statistics.tsv`: tests, effect sizes, and bootstrap intervals
- `leave_one_out.tsv`: pretreatment influence analysis
- `cldn4_samples.tsv`: public sample-level values
- `summary.json`: honest verdict and provenance
- `download_manifest.json`: source URLs and SHA-256 hashes
- `analysis.py`: complete reproducible analysis

## Reproduce

```bash
python3 -m pip install -r requirements.txt
python3 analysis.py
```

The script downloads the two public inputs when absent, checks checksums and
the 24-sample join, and regenerates all tables and figures.

## Source

Nathanson T, Ahuja A, Rubinsteyn A, et al. Somatic mutations and neoepitope
homology in melanomas treated with CTLA-4 blockade. *Cancer Immunol Res*.
2017;5:84–91.
[doi:10.1158/2326-6066.CIR-16-0019](https://doi.org/10.1158/2326-6066.CIR-16-0019)
