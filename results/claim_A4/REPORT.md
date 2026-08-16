# Claim A4 — TISMO Tacstd2 after ICB (honest recount)

## Claim (user)

In TISMO, mouse **Tacstd2** is up in **49 of 64** paired model/treatment
comparisons after ICB, p = **5.8×10⁻⁵**.

## Verdict

**REPRODUCED** for Tacstd2 on TISMO's official ICB gene-module universe.

| Item | User | This recompute |
|---|---|---|
| Comparison universe | 64 models paired | **64** official TISMO ICB groups (All treatments × All tumors) |
| Tacstd2 up / down / tie | 49 / (implied 15) / 0 | **49 / 15 / 0** |
| p = 5.8e-5 | quoted without test name | **Wilcoxon signed-rank two-sided p = 5.84e-05** |
| Binomial two-sided (sign test) | not the quoted p | 2.44e-05 |
| Distinct cell lines inside the 64 | implied “64 models” | **17** cell lines, **22** studies |
| Lung-only | not stated | **2/2** up (LLC only). p_wilcoxon = 0.5 |
| Cldn4, same design | not claimed | native **34/65** up (2 ties); p_wilcoxon = 0.123 |

The quoted p-value is the Wilcoxon signed-rank p on the 64 paired mean
differences, not a binomial/sign-test p (that one is 2.44×10⁻⁵).

## What “64 models” actually is

TISMO's gene module (`/rtismo/gene/downVivoExprn`, gene=`Tacstd2`,
`icbList=["All"]`, `tumorList=["All"]`) returns **64 named comparison
groups**, each with Baseline samples and ICB-treated samples (Responders
and/or Non-responders). That is the universe that yields 49/64.

Those 64 groups are **not** 64 cell lines. They are 17 syngeneic lines
split by study, timepoint, genotype, diet, or combination ICB:

- 17 cell lines, 22 studies
- Cancer types in the 64: mammary (29 groups), melanoma (14), colorectal
  (10), gastric (5), sarcoma (2), HCC (2), **lung (2)**
- 62 unique baseline sample-sets; 2 of them are reused once (64 groups)

Per-group rows: `tables/tacstd2_64_comparisons.tsv`.

## Method (pre-specified for this recount)

1. Download TISMO's official per-sample gene table for Tacstd2 and Cldn4
   (All ICB × All tumors). Snapshot in `data/`.
2. Within each `cell_line` group, take the mean of TISMO's quantified
   expression in treated samples (`Responder` ≠ Baseline) minus the mean
   in Baseline samples.
3. Count up / down / tie. Run:
   - two-sided Wilcoxon signed-rank on the paired deltas (user p)
   - two-sided binomial sign test (mentioned in the claim page; not the
     quoted 5.8e-5)
4. Repeat for Cldn4 (native groups, and aligned to the 64 Tacstd2 stems).
5. Restrict to lung (`Cancer_type == Lung carcinoma` → LLC only).

Values in the TISMO CSV match the public
`TISMO_expressionvivo_profiles.RDS` to floating-point noise on overlapping
sample IDs (Pearson r = 1). A few CSV rows are not in that RDS dump
(non-GEO / extra TISMO samples); they stay in the official 64 because
that is the user-facing table.

## Tacstd2 results

**Official 64:** 49 up, 15 down, 0 tie.
Mean Δ = 0.256; median Δ = 0.116.
Wilcoxon p = 5.84e-05. Binomial two-sided p =
2.44e-05.

Sensitivity (same direction rule, collapsed first):

| Collapse | n | up | Wilcoxon p | Binomial p |
|---|---|---|---|---|
| Official TISMO groups | 64 | 49 | 5.84e-05 | 2.44e-05 |
| Mean Δ per cell line | 17 | 14 | 0.00385 | 0.0127 |
| Mean Δ per study | 22 | 18 | 0.00415 | 0.00434 |
| Non-lung groups | 62 | 47 | 1.35e-04 | 5.78e-05 |
| Lung only (LLC) | 2 | 2 | 0.5 | 0.5 |

The 49/64 count is real on TISMO's grouping. Calling them “64 models” overstates
independence: 29/64 groups are mammary, and several share a study or baseline.

## Lung-only subset

TISMO ICB gene-module lung coverage is **two LLC groups** from GSE155972
(anti-CTLA4 + anti-PD1, parental and Setdb1-KO):

| Group | Tacstd2 Δ | Cldn4 Δ |
|---|---|---|
| LLC_GSE155972_antiCTLA4&antiPD1 | +0.133 | +0.081 |
| LLC_GSE155972_Setdb1_KO_antiCTLA4&antiPD1 | +0.725 | +0.334 |

Both genes go up in both LLC groups. **n = 2. That cannot carry a 49/64,
p = 5.8e-5 claim in lung.** CMT-167 is the only other lung line in TISMO
annotations and is not in this ICB gene-module export.

## Cldn4 (same design, not claimed as 49/64)

TISMO's Cldn4 All×All export has **65** groups, not 64. The extra
stem is `MOC22_RU31562203_antiPD1` (HNSCC; absent from the Tacstd2 export).
Three YTN16 groups have a few more Cldn4 samples than Tacstd2, so TISMO
prints a different `(n=)` suffix.

| Cldn4 universe | n | up | down | tie | Wilcoxon p | Binomial p |
|---|---|---|---|---|---|---|
| Native TISMO groups | 65 | 34 | 29 | 2 | 0.123 | 0.615 |
| Aligned to Tacstd2 64 stems | 64 | 34 | 28 | 2 | 0.0779 | 0.526 |
| Lung only (LLC) | 2 | 2 | 0 | 0 | 0.5 | 0.5 |

Cldn4 is **not** directionally consistent after ICB in this universe
(~half up). Do not recycle the Tacstd2 49/64 sentence for Cldn4.

## What this does **not** show

- Not 64 independent models, and not a lung result.
- Not a paired-mouse longitudinal test. TISMO pairs Baseline vs treated
  *groups* (often isotype / no-treatment vs ICB in the same study).
- Not DESeq2 Wald significance. TISMO's site uses DESeq2 for star labels;
  the 49/64 count is the sign of the mean-expression difference, including
  tiny deltas near zero.
- Not in-vivo vs in-vitro pairing (that intersection is 31 cell lines, not 64).
- Wilcoxon / binomial treat the 64 groups as independent, which they are not.

## Reproduce

```
pip install -r requirements.txt
python scripts/claim_A4_tismo.py
```

Source CSVs are already in `results/claim_A4/data/` (TISMO public API,
2026-08-16 snapshot). The script will use them; it only re-downloads if
a file is missing.

## Files

- `data/tacstd2_tismo_vivo.csv`, `data/cldn4_tismo_vivo.csv` — official exports
- `tables/tacstd2_64_comparisons.tsv` — one row per Tacstd2 group
- `tables/cldn4_comparisons.tsv` — native Cldn4 groups
- `tables/summary.tsv`, `tables/summary.json`
- `figures/tacstd2_waterfall.png`, `figures/cldn4_waterfall.png`

Claim-match flag (64 groups, 49 up, Wilcoxon p ≈ 5.8e-5): **TRUE**
