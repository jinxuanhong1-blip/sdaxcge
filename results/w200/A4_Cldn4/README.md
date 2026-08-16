# W200-A4 · TISMO Cldn4 after ICB (same 64 groups as Tacstd2)

**Task:** A4 analog. User claim A4 is *Tacstd2 up in 49/64 TISMO ICB models, p = 5.8×10⁻⁵*.
This slice asks the same question for **Cldn4**, locked to that 64-group universe.

## Verdict

**The Tacstd2 49/64 sentence is not true for Cldn4.**

On TISMO’s official ICB gene-module export, using the **same 64 comparison
groups** that recover Tacstd2 49/64 (Wilcoxon p = 5.84×10⁻⁵):

| Gene | Up | Down | Tie | Wilcoxon signed-rank p | Sign-test p |
|---|---:|---:|---:|---:|---:|
| Tacstd2 (lock) | 49 | 15 | 0 | **5.84×10⁻⁵** | 2.44×10⁻⁵ |
| **Cldn4 (primary)** | **34** | **28** | **2** | **0.078** | **0.53** |

Cldn4 mean treated − baseline is a weak positive (mean Δ = 0.13 log(TPM);
median Δ = 0.040). That is not 49/64 and not p = 5.8×10⁻⁵. Do not recycle
the Tacstd2 claim for Cldn4.

## What “64 models” actually are

Not 64 cell lines. TISMO’s gene page (All ICB × All tumors, `type=3`) returns
**64 named comparison groups**: `CellLine_GSE_…_treatment`. Each group has
Baseline samples and ICB-treated samples (Responders and/or Non-responders).

- **17 cell lines, 22 studies**
- 29/64 groups are mammary (4T1, T11, KPB25L, p53-2225L, p53-2336R, EMT6, E0771)
- Lung is **2/2 LLC** groups from GSE155972 only
- Contrast used for the claim: **treated vs Baseline**, not R vs NR
- “Up” = sign of the **mean** difference (same rule that yields Tacstd2 49/64)

Tacstd2 lock on this download: **49/15/0, Wilcoxon p = 5.839848×10⁻⁵**. That
is the user-reported p. The binomial/sign test on 49/64 is 2.4×10⁻⁵ and is
**not** the quoted number.

Cldn4’s raw export has **65** groups. The extra stem
`MOC22_RU31562203_antiPD1` is Cldn4-down and is **excluded** so the universe
stays the Tacstd2 64. Including it makes Cldn4 34/29/2, Wilcoxon p = 0.12.

## Results (Cldn4, same 64)

- **34 up / 28 down / 2 tie**
- Wilcoxon signed-rank on the 64 paired mean differences (zeros dropped):
  p = 0.078
- Two-sided sign test excluding ties (34/62): p = 0.53
- Concordance with Tacstd2 direction: **41/64** (Spearman ρ = 0.43, p = 4.2×10⁻⁴)
- Collapse to 17 cell lines (mean of that line’s group Δ): 11 up / 6 down,
  Wilcoxon p = 0.051
- Lung (LLC only): both groups up (Δ = +0.33 Setdb1-KO, +0.081 parental).
  **n = 2 cannot carry 49/64**
- R vs NR: only **2** groups have both classes; both have lower Cldn4 in
  responders (p = 0.5). Not a response test.

Largest Cldn4 increases: `p53-2336R_GSE124821_end_antiCTLA4&antiPD1`
(Δ = +2.53), `YTN16_GSE146027_day7_antiCTLA4` (+2.04),
`YTN16_GSE146027_day14_antiCTLA4` (+1.85).
Largest decreases: `YTN16_GSE146027_day21_antiCTLA4` (−1.72),
`KPB25L_GSE124821_day3_antiCTLA4&antiPD1` (−1.45),
`KPB25L_GSE124821_UV_end_antiCTLA4&antiPD1` (−1.15).

By cancer bucket, no tissue is 49/64-like. Mammary is 14 up / 15 down
(p = 0.64). Melanoma 7/5/2. Colorectal 5/5.

## Sensitivity (still not 49/64)

| Rule | Up / down / tie | Wilcoxon p |
|---|---|---:|
| Primary (mean Δ, all 64) | 34 / 28 / 2 | 0.078 |
| Median Δ | 35 / 21 / 8 | 0.084 |
| Keep only \|Δ\| ≥ 0.1 | 25 / 14 / 0 | 0.048 |
| Keep only \|Δ\| ≥ 0.25 | 19 / 9 / 0 | 0.077 |
| 17 cell-line means | 11 / 6 / 0 | 0.051 |
| All 65 Cldn4 groups | 34 / 29 / 2 | 0.12 |

The \|Δ\| ≥ 0.1 cut is a post-hoc filter. It is barely p < 0.05 and is still
25/39, not 49/64. Primary remains the Tacstd2-matched mean-difference count.

## Honest limits

- These are TISMO group-level Baseline vs treated contrasts, not paired
  pre/post mice.
- Groups are not independent (shared studies, shared baselines, 29 mammary).
- Direction ignores DESeq2 / FDR; tiny near-zero deltas count, as they do
  for Tacstd2 49/64.
- TISMO table download requires `type=3`. Plot-only `getVivoExprn` does not
  return numbers. `downVivoExprn` without `type` returns HTTP 500.
- No FASTQ and no file >2 GB. Official processed export only.

## Data

- TISMO, Zeng et al. *Nucleic Acids Research* 2022 (PMID 34534350)
- `POST https://tismo.pku-genomics.org/rtismo/gene/downVivoExprn`
  `gene={Cldn4|Tacstd2}&icbList=["All"]&tumorList=["All"]&type=3`
- Values are TISMO log(TPM) as served on 2026-08-16

## Files

| File | Role |
|---|---|
| `data/cldn4_tismo_vivo.csv` | Official Cldn4 ICB export (65 groups) |
| `data/tacstd2_tismo_vivo.csv` | Official Tacstd2 ICB export (64 groups; lock) |
| `data/tismo_vivo_meta.json` | TISMO vivo metadata (cancer type map) |
| `tables/cldn4_on_tacstd2_64.tsv` | Per-group Cldn4 Δ on the 64 |
| `tables/tacstd2_64.tsv` | Per-group Tacstd2 Δ (lock) |
| `tables/cldn4_vs_tacstd2_64.tsv` | Paired directions |
| `tables/cldn4_by_cell_line.tsv` | 17-line collapse |
| `tables/cldn4_by_cancer.tsv` | Cancer-type split |
| `tables/cldn4_extra_groups_not_in_64.tsv` | MOC22 only |
| `tables/cldn4_R_vs_NR.tsv` | The 2 groups with both R and NR |
| `tables/summary_stats.tsv` | Primary + sensitivity |
| `figures/cldn4_waterfall_64.png` | Cldn4 Δ waterfall |
| `figures/cldn4_vs_tacstd2_scatter.png` | Concordance scatter |
| `summary.json` / `audit.json` | Verdict + source audit |

## Reproduce

```bash
python3 scripts/w200/A4_Cldn4/download.py
python3 scripts/w200/A4_Cldn4/analyze.py
```
