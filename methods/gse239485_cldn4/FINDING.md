# GSE239485 — Cldn4 and IFN/MHC on the Tacstd2 +2.09 contrast

**Additive** mouse LLC bulk (GEO [GSE239485](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239485); syngeneic Lewis lung carcinoma). Public processed `GSE239485_Processed_data.xlsx` (`DataNorm`, already-log). Unit is the **tumor / mouse**.

**Taken as given (not re-audited):** Tacstd2 **+2.09** on Poly I:C + anti-PD-1 vs vehicle; **TISMO 49/64** Tacstd2 up after ICB. This page scores **Cldn4** and **IFN / MHC-I** on that same contrast.

There is **no aPD-1 monotherapy arm** and **no R vs NR label**. Every treated sample also received Poly I:C.

Reproduce: `python3 methods/gse239485_cldn4/analyze.py` (GEO files cached under `/tmp/gse239485/`).

## One-row table (honest n = 8/8)

Primary contrast: **Poly I:C + anti-PD-1 vs vehicle**. C_ = vehicle (n=8), D_ = doublet (n=8). T_ triplet is extra only.

| dataset | model | contrast | n | Tacstd2 (given) | Cldn4 | IFN mean-z | MHC-I mean-z | Cd8a | R/NR |
|---|---|---|---:|---|---|---|---|---|---|
| GSE239485 | LLC subcutaneous | Poly I:C + aPD-1 vs vehicle | **8 / 8** | **+2.09**, Welch p=6.3×10⁻⁴ | **−0.77**, p=0.24 | **+0.24**, p=0.49 | **+0.81**, p=0.0068 | +0.37, p=0.32 | **not deposited** |

Tacstd2 re-measured on the same 16 columns is +2.093 (Welch p=6.3×10⁻⁴, MWU p=6.2×10⁻⁴) — that locks the given contrast, it is not a new discovery. Cldn4 Cliff δ (treated − vehicle) = −0.34. Full row: `tables/one_row.tsv`.

## How it was scored

- Tacstd2 = ENSMUSG00000051397; Cldn4 = ENSMUSG00000047501 (symbol-confirmed on the matrix).
- IFN = leftover 6-gene list in mouse: **Ifng, Stat1, Cxcl9, Cxcl10, Ido1, H2-Aa** (6/6 present). Score = mean of gene-wise z across the 24 tumors.
- MHC-I = leftover HLA-A/B/C / B2M / TAP1 / TAP2 in mouse: **H2-K1, H2-D1, H2-Q4, B2m, Tap1, Tap2** (6/6). Same mean-z.
- Tests: Welch t and two-sided MWU on the author log matrix (mean difference = log2FC). Spearman on the primary 16 and on all 24.
- Cd8a is the infiltrate proxy (no ESTIMATE list applied here). Partial Spearman is Cldn4 vs IFN/MHC after Cd8a.

## Cldn4 does not track Tacstd2

| feature | mean vehicle | mean doublet | log2FC | Welch p | MWU p | n |
|---|---:|---:|---:|---:|---:|---:|
| Tacstd2 (given) | −3.09 | −1.00 | **+2.09** | **6.3×10⁻⁴** | 6.2×10⁻⁴ | 8/8 |
| Cldn4 | −3.08 | −3.85 | **−0.77** | **0.24** | 0.28 | 8/8 |
| Cd8a | 2.36 | 2.73 | +0.37 | 0.32 | 0.44 | 8/8 |
| Cd274 | 5.15 | 7.10 | **+1.94** | **1.9×10⁻⁶** | 1.6×10⁻⁴ | 8/8 |
| Ptprc | 7.15 | 7.05 | −0.10 | 0.51 | 0.44 | 8/8 |

Cldn4 sits near the same low floor as vehicle Tacstd2 and does not rise. Cd274 (PD-L1) does rise — expected for Poly I:C. Cd8a / Ptprc do not, so the Tacstd2 call is not a global immune-infiltrate shift.

On the primary 16 tumors, Cldn4 vs IFN ρ=+0.09 (p=0.75); vs MHC-I ρ=−0.003 (p=0.99); vs Cd8a ρ=+0.05 (p=0.85). After Cd8a, partial ρ vs IFN = +0.14 (p=0.61) and vs MHC-I = −0.09 (p=0.75). Tacstd2 vs MHC-I ρ=+0.61 (p=0.012); vs Cd274 ρ=+0.76 (p=6.1×10⁻⁴).

## IFN / MHC-I on the same contrast

| score or gene | log2FC or Δz | Welch p | n | note |
|---|---:|---:|---:|---|
| **IFN mean-z (6/6)** | **+0.24** | **0.49** | 8/8 | mixed genes cancel |
| Ifng | +0.18 | 0.75 | 8/8 | |
| Stat1 | +0.43 | 0.14 | 8/8 | |
| Cxcl9 | +0.50 | 0.27 | 8/8 | |
| Cxcl10 | **+1.42** | **0.0016** | 8/8 | chemokine up |
| Ido1 | +0.43 | 0.56 | 8/8 | |
| H2-Aa (MHC-II / HLA-DRA) | **−1.62** | **8.7×10⁻⁷** | 8/8 | pulls IFN mean-z down |
| **MHC-I mean-z (6/6)** | **+0.81** | **0.0068** | 8/8 | |
| H2-K1 | +0.57 | 0.0026 | 8/8 | |
| H2-D1 | +0.37 | 0.0051 | 8/8 | |
| H2-Q4 | +0.95 | 3.3×10⁻⁴ | 8/8 | |
| B2m | +0.61 | 0.0011 | 8/8 | |
| Tap1 | +0.59 | 0.012 | 8/8 | |
| Tap2 | −0.04 | 0.65 | 8/8 | |

IFN as a 6-gene mean is **null**. MHC-I as a 6-gene mean is **up**. That MHC-I rise is a Poly I:C + aPD-1 treatment effect, not a Cldn4-high state: Cldn4 is flat-to-down and uncorrelated with MHC-I on n=16.

## Extra: triplet (not the given contrast)

Poly I:C + anti-PD-1 + C5aR1i vs vehicle, still **n=8/8**, still no monotherapy.

| feature | log2FC / Δz | Welch p |
|---|---:|---:|
| Tacstd2 | +2.58 | 9.7×10⁻⁴ |
| Cldn4 | +0.59 | 0.36 |
| IFN mean-z | +0.84 | 0.034 |
| MHC-I mean-z | +1.87 | 1.2×10⁻⁴ |
| Cd8a | +1.82 | 4.7×10⁻⁴ |
| Ifng | +2.12 | 2.1×10⁻⁴ |

Cldn4 remains non-significant. The triplet moves IFN / Cd8a; the primary doublet does not.

## Two axes (not pooled)

1. **Barrier / immune-low.** Cldn4 does not follow the given Tacstd2 +2.09 rise. Cldn4 vs Cd8a / IFN / MHC-I is null on n=16. This slice does not add a Cldn4-high / immune-low treatment effect.
2. **ADC+ICI candidate (Cldn4-high + IFN/MHC-I).** MHC-I is up on the combo; Cldn4 is not. There is no Cldn4-high / IFN-high or Cldn4-high / MHC-I-high tumor set on the primary 16. This is ICI + TLR3 agonist, not an ADC arm.

## TISMO (given)

**TISMO 49/64** Tacstd2 up after ICB is the external all-cancer prior and is taken as given. This page is extra LLC bulk on GSE239485 only. It does not re-count TISMO.

## What is not here

- No public R vs NR / MPR / PFS field (24 tumors, three treatment arms).
- No aPD-1-only arm — Poly I:C is on every D_ and T_ sample.
- LLC is subcutaneous syngeneic, not an orthotopic GEMM.
- No human data.

## Files

- `analyze.py` — download-once cache + gene lock + Welch/MWU + mean-z + figures
- `tables/one_row.tsv`, `contrasts.tsv`, `per_sample.tsv`, `spearman.tsv`, `partial_spearman.tsv`, `gene_coverage.tsv`, `summary.json`
- `figures/primary_boxes.png` — Tacstd2 / Cldn4 / IFN / MHC-I, n=8/8/8
- `figures/cldn4_vs_ifn_mhc.png` — Cldn4 vs IFN and MHC-I on the primary 16
