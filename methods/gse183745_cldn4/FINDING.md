# GSE183745 leftover — Cldn4 vs immune / Cd274

**Additive only.** Public GEO processed counts. Unit is the **kidney / mouse**. This series is **not lung**, **not ICI**, and **not a tumour RNA matrix**. It is leftover because the Fh1 kidney innate-immunity paper (Zecchini / Frezza, *Nature* 2023, PMID 36890229) was not previously scored for **Cldn4 vs Cd274 / immune**.

GEO [GSE183745](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183745). Status: **Public on Feb 17 2023**. Processed file: `GSE183745_RNAseq_rawCounts.csv.gz` (24,421 symbols × 20 samples). Scale used here: **log2(CPM+1)** from those raw counts. No SRA re-alignment.

Reproduce: `python3 methods/gse183745_cldn4/analyze.py` (GEO files cached under `$GSE183745_CLDN4_DATA`, default `/tmp/gse183745_cldn4`).

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Kidneys in the processed count matrix | yes | **20** | 5 Fh1−/− d5, 5 Fh1−/− d10, 5 Fh1+/+ d5, 5 Fh1+/+ d10 |
| Unique GSM | yes | **20** | GSM5570173–GSM5570192 |
| Tissue | yes | 20 | whole kidney after tamoxifen; mixed C57BL/6 129/SvJ |
| Lung / LUAD / LUSC | **no** | **0** | mouse kidney only |
| ICI / R vs NR / MPR / PFS / OS | **no** | **0** | Fh1 induction only (3 × 2 mg IP tamoxifen) |
| Human samples | **no** | **0** | *Mus musculus* |
| Cldn4 finite | yes | **20** | symbol `Cldn4`; mean CPM 67.4; 20/20 nonzero |
| Cd274 finite | yes | **20** | symbol `Cd274`; mean CPM 8.1; 20/20 nonzero |
| Ifng detected | yes | **3** | **3 total counts** across 20 kidneys (1+1+1); 17 zeros |
| **Primary pairwise n (Cldn4 + Cd274)** | yes | **20** | this is the n used below |

Do not write a lung-ICI n. The computable public n is **20 kidneys**. Ifng is effectively **empty**.

## One-row table

| dataset | model | n | Cldn4–Cd274 ρ (p) | Cldn4–Cd8a ρ (p) | Cldn4–IFN ρ (p) | Cldn4–MHC-I ρ (p) | MHC-I \| EPI ρ (p) | Ifng | lung ICI |
|---|---|---:|---|---|---|---|---|---|---|
| GSE183745 | kidney Fh1 inducible KO | **20** | **+0.195 (0.41)** | +0.077 (0.75) | +0.274 (0.24) | +0.626 (0.0032) | **−0.047 (0.85)** | **3/20** | **no** |

Cd274 verdict: **NO_EVIDENCE**. Full row: `tables/one_row.tsv`.

## How it was scored

- Cldn4 / Cd274 / immune genes are **symbols on the author count matrix** (no Ensembl column deposited).
- IFN = leftover 6-gene list in mouse: **Ifng, Stat1, Cxcl9, Cxcl10, Ido1, H2-Aa** (6/6 present). Score = mean of gene-wise z across the 20 kidneys. **Ifng does not contribute** (3 reads total).
- MHC-I = leftover HLA-A/B/C / B2M / TAP1 / TAP2 in mouse: **H2-K1, H2-D1, H2-Q4, B2m, Tap1, Tap2** (6/6).
- Epithelial mean-z = **Epcam, Cdh1, Krt8, Krt18** (4/4). Used only as a residual, not as a claim.
- Tests: Spearman on n=20 (bootstrap 95% CI, 2,000 resamples, seed `20260817`). Partial Spearman residualises ranks on Cd8a or on EPI. Welch / MWU for Fh1−/− vs Fh1+/+ at day 5 (5/5), day 10 (5/5), and pooled (10/10).
- Fh1 counts are **not** a KO lock. The floxed model still maps residual Fh1 exons (mean CPM 262 in both genotypes). Genotype is GEO metadata.

## Cldn4 vs Cd274 (primary leftover)

| pair | n | ρ | 95% CI | p | ρ \| Cd8a | p | ρ \| EPI | p | verdict |
|---|---:|---:|---|---:|---:|---:|---:|---:|---|
| Cldn4 vs **Cd274** | **20** | **+0.195** | −0.230 to +0.572 | 0.41 | +0.110 | 0.65 | −0.308 | 0.19 | **NO_EVIDENCE** |
| Cldn4 Q4 vs Q1 on Cd274 | 5 vs 5 | — | — | MWU 0.55 | — | — | — | — | null (Δ +0.13) |

Cldn4-high is not Cd274-high or Cd274-low. The CI includes both a modest negative and a modest positive effect. Within Fh1+/+ (n=10) ρ=+0.61 (p=0.060); within Fh1−/− (n=10) ρ=−0.15 (p=0.68). Do not pool those two subset ρ values into a claim.

## Cldn4 vs immune (pre-specified leftover)

| pair | n | ρ | 95% CI | p | ρ \| EPI | p | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| Cldn4 vs **Cd8a** | 20 | +0.077 | −0.420 to +0.523 | 0.75 | +0.113 | 0.64 | **NO_EVIDENCE** |
| Cldn4 vs **IFN mean-z (6/6)** | 20 | +0.274 | −0.205 to +0.637 | 0.24 | −0.123 | 0.60 | **NO_EVIDENCE** |
| Cldn4 vs **Ifng** | 20 | −0.370 | — | 0.11 | −0.224 | 0.34 | **EMPTY** (3 counts) |
| Cldn4 vs **MHC-I mean-z (6/6)** | 20 | **+0.626** | +0.242 to +0.858 | 0.0032 | **−0.047** | **0.85** | crude +; **null after epithelium** |
| Cldn4 vs Ptprc | 20 | −0.006 | −0.530 to +0.523 | 0.98 | −0.259 | 0.27 | NO_EVIDENCE |

The crude MHC-I correlation is an epithelial / tight-junction track, not an immune-infiltrate track:

| pair | n | ρ | p |
|---|---:|---:|---:|
| Cldn4 vs epithelial mean-z | 20 | **+0.904** | 4.7×10⁻⁸ |
| Cldn4 vs Tacstd2 | 20 | **+0.923** | 6.6×10⁻⁹ |
| Cldn4 vs Cldn7 | 20 | +0.821 | 9.1×10⁻⁶ |
| Cldn4 vs Cldn3 | 20 | +0.627 | 0.0031 |

Dropping the two high-Cldn4 Fh1−/− day-10 kidneys (K_HOM_10_1, K_HOM_10_3; CPM 313 and 351 vs ~24–88 in the other 18) leaves Cldn4–Cd274 ρ=+0.23 (p=0.35, n=18) and Cldn4–MHC ρ=+0.49 (p=0.041), still killed by the EPI residual on the full 20.

Q4 vs Q1 (5 vs 5) on MHC-I is MWU p=0.016, but the same cut on EPI is MWU p=0.0079 (Cliff δ = +1.0). That is not an immune claim.

## Fh1−/− vs Fh1+/+ (not the leftover claim)

Author contrast is acute Fh1 loss in adult kidney. Cldn4 and Cd274 do **not** move with genotype.

| feature | day 5 (5/5) Δ | Welch p | day 10 (5/5) Δ | Welch p | pooled (10/10) Δ | Welch p |
|---|---:|---:|---:|---:|---:|---:|
| Cldn4 | −0.18 | 0.51 | +0.98 | 0.30 | +0.40 | 0.43 |
| Cd274 | −0.06 | 0.84 | −0.52 | 0.10 | −0.29 | 0.17 |
| Cd8a | **+0.83** | **0.0077** | −0.10 | 0.73 | +0.37 | 0.12 |
| IFN mean-z | −0.66 | 0.13 | −0.64 | 0.041 | **−0.65** | **0.0086** |
| MHC-I mean-z | −0.80 | 0.22 | +0.54 | 0.41 | −0.13 | 0.78 |
| Ptprc | −0.41 | 0.28 | −1.49 | 0.020 | −0.95 | 0.0059 |

Day-5 Fh1−/− kidneys are Cd8a-higher (5/5). That is a genotype/time effect, not a Cldn4-high state. Pooled IFN mean-z is lower in Fh1−/−; Ifng itself is empty on both arms. This page does not re-audit the paper’s cGAS–STING innate-immunity claim.

## Two axes (not pooled)

1. **Barrier / immune-low.** Cldn4 vs Cd8a / Cd274 / IFN is null on n=20. Cldn4 vs MHC-I is crude-positive and **null after epithelium**. This leftover does not add a Cldn4-high / immune-low kidney.
2. **ADC+ICI candidate (Cldn4-high + IFN/MHC-I / Cd274).** Cd274 is null. Ifng is empty. MHC-I tracks epithelium, not a Cldn4-high IFN-high set. This is not an ADC arm and not ICI.

## What is not here

- No public lung sample, no LUAD/LUSC split, no ICI label.
- No human HLRCC tumour RNA in this series (the paper’s human pieces are elsewhere).
- No pathologist IHC for CLDN4 or PD-L1.
- Do not treat n=20 mouse kidneys as a lung-ICI prior.

## Files

- `analyze.py` — download-once cache + CPM + Spearman / partial / Welch + figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/one_row.tsv` — one leftover row
- `tables/gene_coverage.tsv` — symbol presence, mean CPM, nonzero n
- `tables/per_sample.tsv` — 20 kidneys
- `tables/spearman.tsv` — Cldn4 vs Cd274 / immune / MHC / EPI
- `tables/spearman_subsets.tsv` — genotype / day / drop-two sensitivity
- `tables/contrasts.tsv` — Fh1−/− vs Fh1+/+
- `tables/q4_vs_q1.tsv` — Cldn4 Q4 vs Q1 (5 vs 5)
- `tables/summary.json`
- `figures/fig1_group_boxes.png` — Cldn4 / Cd274 / IFN / MHC-I, n=5/5/5/5
- `figures/fig2_scatters.png` — Cldn4 vs Cd274, Cd8a, MHC-I
- `figures/fig3_forest.png` — Spearman forest, n=20
