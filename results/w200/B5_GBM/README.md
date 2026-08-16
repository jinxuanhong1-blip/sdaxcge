# B5 analog: GBM/glioma ICI public RNA, CLDN4 versus outcome

## Result

Public GBM/glioma ICI RNA **does not support** a usable CLDN4-high → worse
response biomarker, and it does **not** reproduce a radiographic-response
odds ratio.

There is **no public RECIST/ORR table** joined to these RNA matrices. The
available clinical endpoint is survival.

In the only strand-consistent public bulk tumor cohort (GSE121810, n = 29
recurrent GBM, pembrolizumab):

- CLDN4 is at the detection floor: median 0.48 CPM (IQR 0.30–0.64). Only
  3/29 tumors have CPM ≥ 1. GFAP median is 11,344 CPM in the same matrix.
- Unadjusted continuous Cox OS: HR 2.52 per log2(CPM+1) (95% CI 1.16–5.47),
  p = 0.019, concordance 0.59. Direction matches B5 (higher CLDN4, worse OS).
- That Cox result is leverage from two samples above the floor (Pt15_B,
  8.03 CPM, OS 102 d; Pt33_A, 3.04 CPM, OS 180 d). After removing those two
  samples, HR 1.75 (95% CI 0.29–10.7), p = 0.55. After removing all three
  CPM ≥ 1 samples, HR 1.63, p = 0.67.
- Median-split OS log-rank p = 0.130; Cox HR(high vs low) 1.98
  (95% CI 0.80–4.86), p = 0.138.
- Spearman CLDN4 vs OS: ρ = −0.26, p = 0.181.
- OS ≥ 365 d vs < 365 d (not RECIST): Mann–Whitney p = 0.120.
- Pretreatment adjuvant arm only (n = 15): continuous Cox p = 0.171.
- PFS continuous Cox p = 0.157; median-split log-rank p = 0.045 (exploratory,
  not confirmed by the continuous model).

Two larger public count matrices (GSE264695 n = 20; GSE318645 n = 47) have
OS labels but **fail strandedness QC**. Their CLDN4 values are not
biologically trustworthy. Both are null for OS (GSE318645 continuous Cox
p = 0.636).

Zhao 2019 (the GBM ICI cohort that actually used a radiographic-style
responder definition) deposits raw SRA only. Processed counts and a public
sample–response map were not available, so it was not analyzed.

## What was tested

Claim B5, as stated in this project, is an 11-cohort ICI meta-analysis with
CLDN4-high worse response (user-reported OR = 0.42). This folder is the
**GBM/glioma analog only**. It does not assemble or pool those 11 cohorts.

Pre-specified choices:

1. Survey public human glioma ICI transcriptomes with a response or survival
   label.
2. Use bulk tumor RNA with a joinable public clinical table.
3. Gate expression matrices on strandedness (MT-ND6 / H-strand MT and
   antisense-gene fraction) before interpreting CLDN4.
4. Primary endpoint: OS. Secondary: PFS. Exploratory binary: OS ≥ 365 days,
   labeled as not RECIST.
5. Primary predictor: log2(CLDN4 CPM + 1). Secondary: cohort-internal median
   split.
6. Report detectability before association.

## Cohort and endpoint

Primary source: GSE121810, Cloughesy et al., *Nature Medicine* 2019,
doi:10.1038/s41591-018-0337-7. Counts are the GEO HUGO matrix
`GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.xlsx`.

Survival and covariates come from Supplementary Table sheet `SFig4C` of
de Groot et al., *Nature Communications* 2024,
doi:10.1038/s41467-024-54326-7, which includes the original 29 GSE121810
patients plus the GSE264695 expansion. All 29 + 20 sample IDs join
one-to-one.

Arm A is neoadjuvant pembrolizumab: tumor RNA is after one dose. Arm B is
adjuvant pembrolizumab: tumor RNA is pretreatment. Neoadjuvant patients
lived longer in this public table (OS Mann–Whitney p = 0.003). CLDN4 did
not differ by arm (p = 0.585). Arm-adjusted Cox still gives CLDN4 HR 2.32,
p = 0.035, because the two high-CLDN4 early deaths remain.

## Files

- `analysis_samples.csv`: sample-level CLDN4, arm, OS/PFS, QC flag
- `cohort_inventory.csv`: datasets considered and why they were used or not
- `results.json`: QC, detectability, primary stats, robustness, inventory
- `CLDN4_detectability.png`: GSE121810 log2(CPM+1) for CLDN4 vs control genes
- `CLDN4_vs_OS_GSE121810.png`: Kaplan–Meier by CLDN4 median split
- `CLDN4_vs_OS365.png`: CLDN4 vs 1-year OS (not RECIST)
- `run_analysis.py`: download, checksum, QC, analysis, and plotting
- `requirements.txt`: Python dependencies

## Reproduce

From the repository root:

```bash
python3 -m pip install -r results/w200/B5_GBM/requirements.txt
python3 results/w200/B5_GBM/run_analysis.py
```

The script downloads GEO matrices and the Nature Communications workbook
into `.cache/B5_GBM/`, checks SHA-256 digests, and recreates the CSV, JSON,
and PNG outputs.

## Limitations

This is a post hoc single-gene test in 29 tumors. GBM is a glial, not
epithelial, histology; CLDN4 is barely expressed, so a CLDN4-high split
does not create two biologically distinct groups. Survival under
pembrolizumab is not the same as RECIST response, and it cannot separate a
predictive ICI interaction from a general prognostic association. Neoadjuvant
RNA is on-treatment. The two larger GEO count matrices are reverse-stranded
and should not be used to rescue or refute the signal. Zhao 2019 remains
the missing public radiographic-response RNA cohort. No multiplicity
correction was applied beyond the planned detectability-gated sensitivity
checks.
