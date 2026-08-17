# GSE297630 — Cldn4 and IFN/MHC on the Tacstd2 +2.05 matrix

Additive mouse slice. **Tacstd2 +2.05 is given** and is reproduced on the public processed table. This note only adds **Cldn4** and **IFN/MHC** on that same matrix.

**Series:** Sumii et al., GEO [GSE297630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297630). LLC subcutaneous tumors; Clariom S Mouse. Control (C, no anti-PD-1) vs cells/tissue recovered after anti-PD-1 (P, “tolerant”). Processed file: `GSE297630_processed_data.xlsx`.

## Honest n

| Unit | Count |
|---|---:|
| Samples | **6** |
| Control (C) | **3** |
| Anti-PD-1 tolerant (P) | **3** |
| Smallest two-sided Mann–Whitney U p at 3 vs 3 | **0.1** |

There is no larger public per-sample matrix for this series. Welch *p* at n=3 can look tiny when within-group SD is small (array technical tightness). **MWU cannot reject at α=0.05.** Spearman on n=6 is treatment-driven and one-or-two-point fragile. Do not write “significant IFN opening” from the union score alone.

## Verdict

| Readout | n | log2FC or Δz (P − C) | Welch p | MWU p |
|---|---|---:|---:|---:|
| **Tacstd2 (given)** | 3 vs 3 | **+2.05** | 4.2×10⁻⁵ | **0.1** |
| **Cldn4** | 3 vs 3 | **−0.19** | 0.016 | **0.1** |
| IFN ISG set (54/54 genes) | 3 vs 3 | Δz **+0.13** | 0.21 | **0.1** |
| MHC-I + B2M (7/7) | 3 vs 3 | Δz **+1.72** | 0.0014 | **0.1** |
| MHC-II + CIITA (8/8) | 3 vs 3 | Δz **−0.06** | 0.84 | 1.0 |
| APM / immunoproteasome (12/12) | 3 vs 3 | Δz **+0.04** | 0.62 | 0.7 |
| IFN/MHC union (81/81) | 3 vs 3 | Δz **+0.24** | 0.018 | **0.1** |

**What holds on this matrix**

- Tacstd2 is higher in the P (anti-PD-1 tolerant) Signal columns than in C. Magnitude **+2.05** matches the given number exactly (C 13.10, 13.12, 13.19; P 15.09, 15.18, 15.19).
- Cldn4 moves the **other way**, and only a little (**−0.19**). It does not co-induce with Tacstd2.
- Classical **MHC-I + B2M** is the only coherent IFN/MHC block: all seven genes up (HLA-A/B/C/E/F/G, B2M).
- The broad **IFN ISG** set is a coin flip (**28 up / 26 down**). CORE6 is split (IFI27 +2.07, HLA-A +0.57, MX1 +0.77, OAS2 +0.69; IFIT1 −0.64, ISG15 −0.18).

**What does not hold**

- Cldn4 is **not** +2-like. Calling this a Tacstd2/Cldn4 co-up program is false on these six samples.
- IFN/MHC as one program is **not** uniformly open. The union Δz is MHC-I-driven. Immunoproteasome genes PSMB8/PSMB9 and CIITA go **down**.
- No MWU test here is <0.1. n=3 vs 3 cannot support a confirmatory rank test.

## Contrast (same columns as the given +2.05)

Per-sample Signal columns `C-1/C-2/C-3` and `P-1/P-2/P-3` match the GEO CEL names (`GSM8995899_C-1_…` = control; `GSM8995902_P-1_…` = anti-PD-1). log2FC = mean(P) − mean(C) on those six columns.

The TAC summary columns `C Avg (log2)` / `P Avg (log2)` are **swapped** relative to those Signal/CEL labels (for TACSTD2, “C Avg”=15.18 equals the P-Signal mean). Author linear FC 4.2 is |2^2.05|. This note follows the Signal+CEL labeling that produces the **given +2.05**, not the swapped summary pair.

The deposited Gene Symbol column is **human** (TACSTD2, CLDN4, HLA-A). Mouse H2-* symbols are absent. Scoring uses that column so Cldn4 / IFN/MHC stay on the same public matrix as Tacstd2 +2.05.

## Per-sample (log2 RMA Signal)

| Sample | Group | Tacstd2 | Cldn4 |
|---|---|---:|---:|
| C1 | control | 13.19 | 14.12 |
| C2 | control | 13.12 | 14.21 |
| C3 | control | 13.00 | 14.18 |
| P1 | anti-PD-1 tolerant | 15.18 | 13.91 |
| P2 | anti-PD-1 tolerant | 15.19 | 14.03 |
| P3 | anti-PD-1 tolerant | 15.09 | 14.00 |

## IFN/MHC gene-level (pre-specified; not mined)

CORE6:

| Gene | log2FC (P−C) | Welch p | MWU p |
|---|---:|---:|---:|
| IFI27 | +2.07 | 1.3×10⁻⁵ | 0.1 |
| MX1 | +0.77 | 0.031 | 0.1 |
| OAS2 | +0.69 | 0.068 | 0.1 |
| HLA-A | +0.57 | 0.0043 | 0.1 |
| ISG15 | −0.18 | 0.19 | 0.1 |
| IFIT1 | −0.64 | 0.0039 | 0.1 |

MHC-I + B2M (all up): HLA-B +1.02, HLA-C +0.88, HLA-A +0.57, HLA-E +0.57, B2M +0.51, HLA-G +0.51, HLA-F +0.47.

Selected mixed/down genes that block a blanket “IFN/MHC up” sentence: OASL −2.43, IFIT3 −1.75, IFI44 −1.63, **CIITA −1.30**, PSMB9 −0.74, TAP2 −0.53, PSMB8 −0.40. Ifng +0.26 (p=0.35); Cxcl9 +0.07; Cxcl10 −0.16.

## Correlations (n=6, exploratory)

| Pair | Spearman ρ | p |
|---|---:|---:|
| Tacstd2 vs Cldn4 | −0.77 | 0.072 |
| Tacstd2 vs IFN/MHC z | +1.00 | — |
| Cldn4 vs MHC-I z | −0.83 | 0.042 |
| Cldn4 vs MHC_APM z | −0.94 | 0.0048 |

The Tacstd2–IFN/MHC ρ=+1 is the two treatment clusters lining up, not an independent within-group slope. Do not cite n=6 Spearman as a second cohort.

## Methods

1. Download `GSE297630_processed_data.xlsx` from GEO FTP (Expression sheet, skip 4 header rows).
2. Collapse duplicate Gene Symbol rows by mean. Match case-insensitively (TACSTD2 / CLDN4).
3. Groups from Signal column prefixes C vs P (CEL-confirmed).
4. Gene test: mean log2 difference, Welch two-sample *t*, two-sided Mann–Whitney U.
5. Set score: z-score each gene across the six samples, then mean. Same two tests on the three vs three scores.
6. Gene sets are pre-specified in `gene_sets.py` (CORE6 / IFN_ISG / MHC-I / MHC-II / APM). No Hallmark GMT was required; every listed symbol is present in the table.

```bash
python3 methods/gse297630_cldn4/score_gse297630.py
```

## Read this as additive, not as a new Tacstd2 claim

- Tacstd2 **+2.05** is the given number on this public matrix. Reproduced.
- Cldn4 on the same six samples is **−0.19** (Welch p=0.016, MWU p=0.1). Small, opposite sign.
- IFN/MHC: **MHC-I up, ISGs split, MHC-II/APM mixed.** Honest n is **3 vs 3**.

## Files

- `score_gse297630.py` / `gene_sets.py` — download + score
- `results/targets.tsv` — Tacstd2, Cldn4
- `results/set_scores.tsv` — IFN/MHC set scores
- `results/ifn_mhc_genes.tsv` — per-gene IFN/MHC
- `results/per_sample.tsv` / `correlations.tsv` / `summary.json`
