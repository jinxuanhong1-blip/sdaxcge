# B5: leftover public HNSCC ICI cohorts

This is an eligibility audit and marker analysis for **TACSTD2** and **CLDN4**.
It avoids re-running the public lung ICI cohorts already covered elsewhere.
The primary analyzable leftover is Liu et al. 2021
([PMID 34755131](https://pubmed.ncbi.nlm.nih.gov/34755131/),
[GSE179730](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179730)).

## Reproduce

Python 3 is the only requirement:

```sh
python3 results/w200/B5_HNSC_extra/analyze.py
```

The script audits marker presence directly in each deposited matrix. For
GSE179730 it uses pretreatment tumors only and locks outcomes to Supplementary
Table S2. The primary contrast follows the authors' definition: Responder +
Stable (clinical benefit) versus Progressor. A separately labeled sensitivity
contrast compares pathologic Responder against Stable + Progressor.

GEO calls the GSE179730 matrix “log2 CPM,” but each deposited column sums to
approximately one million and values reach hundreds of thousands. It is
therefore demonstrably linear CPM; the analysis uses `log2(CPM+1)`. Inference
uses all possible group-label allocations for an exact, tie-preserving
Mann–Whitney permutation p-value. BH correction covers the two locked primary
marker tests. No threshold or subgroup was optimized.

## Outputs

- `cohort_catalog.tsv`: eligibility and exclusions
- `per_sample_expression.csv`: auditable marker values and outcomes
- `marker_statistics.csv`: primary and sensitivity comparisons
- `summary.json`: machine-readable methods, checks, and honest verdict
- `RESULTS.md`: ≤200-word result
- `data/`: unmodified compressed GEO downloads

Prat/GSE93157 and Foy/GSE159067 cannot answer the marker question: both
targeted panels omit TACSTD2 and CLDN4. GSE179730 measures both but has only 11
pretreatment tumors and highly sparse expression, so it is exploratory, not a
validation cohort.
