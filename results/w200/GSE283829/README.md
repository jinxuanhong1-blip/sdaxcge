# GSE283829 leftover: TACSTD2 / CLDN4 vs ICI response (NSCLC)

**Bottom line: no persuasive association.** In this public leftover ICI
RNA-seq cohort (n = 27; CR 7 / SD 10 / PD 10), neither
TACSTD2 nor CLDN4 tumor expression separates complete responders from
progressors. Point estimates are small, confidence intervals include
chance, and both genes sit in the middle of the genome-wide rank — not
near a response-associated tail.

This is a leftover GEO series: the series matrix has **no expression
table**. Counts were taken from the public supplementary raw-count
matrix. The analysis is open and was run from those files.

## Cohort

- GEO [GSE283829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283829),
  Lindberg et al., *J Thorac Oncol* (2025), PMID 39743139.
- 27 ICI-treated lung-cancer RNA-seq samples (Illumina NovaSeq 6000,
  GPL24676). The parent study profiled PD1–PD-L1 interactions by in situ
  PLA in a larger biopsy series; this deposit is the RNA-seq subset used
  to look for resistance programs in PLA-high non-responders.
- GEO field `disease stage` is **not TNM**. Values are CR / SD / PD and
  are treated as best response. There is **no PR** category in the
  deposit (either none were sequenced or PR was not labeled).
- Histology: AC 16, SqCC 9, other 2. PLA: high 15,
  low 12. Two library batches (batch 2 / batch 3).
- Genes: TACSTD2 = ENSG00000184292, CLDN4 = ENSG00000189143. Expression
  is log2(CPM+1) from the deposited raw counts. No gene-length / TPM
  conversion is possible from this file alone.

## Pre-specified tests

Primary: two-sided exact Mann–Whitney U on log2(CPM+1), **CR vs PD**.
BH q is across the two primary genes only. AUC treats CR as the positive
class (higher expression predicting CR). Sensitivity contrasts and the
PLA split are reported as such and were not folded into the two-gene FDR.

| Feature | Contrast | n | Median pos | Median neg | Δ median | Exact MWU p | BH q (2 genes) | AUC (95% CI) |
|---|---|---|---:|---:|---:|---:|---:|---|
| TACSTD2 | CR vs PD | 7 vs 10 | 7.090 | 7.423 | -0.333 | 0.813 | 0.962 | 0.457 (0.157–0.771) |
| CLDN4 | CR vs PD | 7 vs 10 | 6.248 | 6.659 | -0.412 | 0.962 | 0.962 | 0.486 (0.171–0.786) |
| combined z (exploratory) | CR vs PD | 7 vs 10 | 0.100 | 0.177 | -0.077 | 0.887 | — | 0.471 (0.171–0.786) |
| TACSTD2 | CR vs SD+PD | 7 vs 20 | 7.090 | 7.423 | -0.333 | 0.646 | — | 0.436 (0.171–0.743) |
| CLDN4 | CR vs SD+PD | 7 vs 20 | 6.248 | 6.636 | -0.388 | 0.978 | — | 0.493 (0.207–0.807) |
| TACSTD2 | CR+SD vs PD | 17 vs 10 | 7.202 | 7.423 | -0.221 | 0.941 | — | 0.512 (0.288–0.712) |
| CLDN4 | CR+SD vs PD | 17 vs 10 | 6.574 | 6.659 | -0.086 | 0.824 | — | 0.471 (0.241–0.700) |

Kruskal–Wallis across CR/SD/PD: TACSTD2 p=0.826, CLDN4
p=0.967. Spearman vs ordinal response (CR=2, SD=1, PD=0):
TACSTD2 ρ=-0.037 (p=0.856), CLDN4
ρ=-0.038 (p=0.851).

TACSTD2 vs CLDN4 Spearman ρ=0.542 (p=0.003).
The two genes are correlated, so they are not independent tests of a
tight-junction / TROP2 axis.

Direction on the primary contrast: TACSTD2 is lower in the positive class;
CLDN4 is lower in the positive class. Neither is a large, consistent shift.
Two CR samples (105691_047 and 105691_065; both tumor type `other`,
PLA-high, batch 2) sit near the floor for both genes and pull the CR
*mean* down; the **medians** still overlap PD. All 7 CRs are PLA-high,
which is the authors' PLA–response observation, not a TACSTD2/CLDN4
result.

## PLA split (exploratory; authors' grouping)

The deposit is grouped by PD1–PD-L1 PLA, not by a TACSTD2/CLDN4
hypothesis. PLA-high vs PLA-low:

| Feature | Median high | Median low | Exact MWU p | AUC (high as +) |
|---|---:|---:|---:|---|
| TACSTD2 | 7.090 | 7.530 | 0.256 | 0.367 (0.167–0.583) |
| CLDN4 | 6.090 | 6.670 | 0.277 | 0.372 (0.156–0.606) |

This is not a response test. It is included so a leftover scan does not
confuse the authors' PLA grouping with RECIST.

## Genome-wide calibration (protein-coding, median CPM ≥ 1)

CR vs PD: TACSTD2 rank 12938 / 16187 (p=0.813, Δmedian=-0.333);
CLDN4 rank 14731 / 16187 (p=0.962, Δmedian=-0.412).

CR vs SD+PD: TACSTD2 rank 12610 / 16747 (p=0.646, Δmedian=-0.333);
CLDN4 rank 16057 / 16747 (p=0.978, Δmedian=-0.388).

A leftover claim that these two genes mark ICI response would require
them to stand out. They do not.

## Honest caveats

- **n = 27, 7 CRs.** Only a large effect (AUC ≳ 0.85) would be reliably
  detected. A null here does not prove no association; it also does not
  support one. AUC intervals all include 0.5 on the primary contrast.
- No PR labels, no PFS/OS in GEO, no PD-L1 IHC percent, no treatment
  line, no EGFR/ALK, no purity. Single-arm ICI series: prognostic vs
  predictive cannot be separated.
- Two library batches. No batch-adjusted primary model was fit; that
  would be underpowered and easy to overfit. Library sizes and
  per-sample values are in `sample_data.csv`.
- Bulk diagnostic-biopsy RNA. TACSTD2/CLDN4 are epithelial; composition
  and purity can move both genes without a cell-intrinsic ICI effect.
- Context genes (CD274, PDCD1, EOMES, HAVCR1, JAML, FCRL1) are in
  `context_stats.csv` as a sanity check against the paper, not as a
  new biomarker hunt. On CR vs PD, JAML and FCRL1 are nominally
  lower in CR (exact p=0.033 and p=0.033); that would not
  survive a 6-gene BH correction and is not the leftover TACSTD2/CLDN4
  question.
- No cutpoint search, no multivariable classifier, no in-sample ROC
  optimization.

## Files

- `sample_data.csv` — per-sample GEO metadata, library size, CPM and
  log2(CPM+1) for the primary and context genes, combined z.
- `stats.csv` — Mann–Whitney / AUC for all reported contrasts.
- `kruskal_spearman.csv` — three-level and ordinal tests.
- `gene_correlation.csv` — TACSTD2 vs CLDN4 Spearman.
- `genomewide_rank_summary.csv` — TACSTD2/CLDN4 ranks among
  protein-coding genes with median CPM ≥ 1.
- `context_stats.csv` — context-gene CR vs PD only.
- `boxplots_response.png`, `boxplots_cr_vs_pd.png`,
  `scatter_tacstd2_cldn4.png`.
- `analysis_manifest.txt` — checksums, versions, seed.
- Code: `scripts/gse283829_tacstd2_cldn4_response.py` (downloads if
  needed; seed 20260816).

## Sources

1. GEO [GSE283829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283829).
2. Lindberg et al. *J Thorac Oncol* (2025), PMID 39743139.
