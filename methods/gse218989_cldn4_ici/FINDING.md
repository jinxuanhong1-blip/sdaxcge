# GSE218989 — CLDN4 / TACSTD2 in public NSCLC ICI bulk TPM

**Additive only.** User thesis (given, not re-audited): CLDN4-high = tight-junction (TJ) barrier; CLDN4-high may mark ADC+ICI candidates; TROP2-high = immune-low.

Public SMC–KAIST PD-1/PD-L1 NSCLC RNA-seq (Kang et al., *Nat Commun* 2024, PMID 38744958; GEO [GSE218989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218989)). Unit is the **patient**. No slide was re-scored.

## One-row table

| dataset | n | genes | CLDN4 / TACSTD2 | GEO R / NR | CLDN4 vs R MWU p (AUC) | TACSTD2 vs R MWU p (AUC) | CLDN4–CD8A ρ (p) | TACSTD2–CD8A ρ (p) | CLDN4 OS HR/SD (p) | histology |
|---|---:|---:|---|---:|---|---|---|---|---|---|
| GSE218989 SMC-KAIST | 355 | 19916 | both on matrix | 168 / 187 | 0.40 (0.47) | 0.67 (0.51) | −0.172 (0.0011) | −0.262 (5.3×10⁻⁷) | 1.19 (0.0081) | **not deposited** |

Full numeric row: `tables/one_row.tsv`.

## Matrix and labels (nothing invented)

`GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz`: **19,916** protein-coding genes × **355** patients (`SMC__Pat*`). **CLDN4** and **TACSTD2** are present. One column per patient.

GEO series-matrix characteristics (all 355):

| field | public? | n | what is there |
|---|---|---:|---|
| treatment | yes | 355 | PD-1/PD-L1 inhibitor |
| treatment outcome | yes | 355 | Responder 168 / Non-responder 187 |
| ethnicity | yes | 355 | Korean |
| histology (LUAD vs LUSC) | **no** | 0 | not a GEO characteristic |
| PFS / OS | **no** | 0 | not on GEO |

Published Supplementary Data 8 (author-correction XLSX, *Nat Commun* 2025; sheet `Clinical_table_v230613`; 497 rows) joins all 355 TPM IDs and adds **OS days + Death** and **PFS days**. It has **no histology column** and **no PFS event column**. Three patients differ from GEO outcome (GEO Non-responder, Supp8 Responder=1): `SMC__Pat50`, `SMC__Pat118`, `SMC__Pat239` — the three IDs named in the author correction. Primary response tests use **GEO**. Survival uses Supp8 OS/Death. **LUAD vs LUSC is not split** because no public histology label exists.

Paper text mentions 497 ICI lung transcriptomes; only 355 are on the public TPM matrix. Tests stay on those 355.

## Immune axes (computable; n=355)

Scores: `log2(TPM+1)` for single genes. Ayers IFN-γ 6-gene and MHC-I (HLA-A/B/C) = mean of gene-wise z-scores (6/6 and 3/3 genes present). ESTIMATE Immune/Stromal = Yoshihara 141-gene lists (140/141 immune, 136/141 stromal present): primary **mean-z**; companion rank ssGSEA (not the R `estimate` purity transform). TJ companion = mean-z of CLDN3/CLDN7/OCLN/TJP1/F11R/MARVELD2/CRB3/PARD3/CGN (9/9). Spearman; BH *q* on the 10 pre-specified gene × {CD8A, IFNG, Ayers6, MHC-I, Immune mean-z} tests.

| pair | n | ρ | p | q (BH, core 10) |
|---|---:|---:|---:|---:|
| CLDN4 vs CD8A | 355 | **−0.172** | 0.00113 | 0.0056 |
| TACSTD2 vs CD8A | 355 | **−0.262** | 5.3×10⁻⁷ | 5.3×10⁻⁶ |
| CLDN4 vs IFNG | 355 | +0.080 | 0.13 | 0.34 |
| TACSTD2 vs IFNG | 355 | −0.015 | 0.78 | 0.98 |
| CLDN4 vs Ayers6 | 355 | +0.029 | 0.58 | 0.97 |
| TACSTD2 vs Ayers6 | 355 | −0.003 | 0.95 | 0.98 |
| CLDN4 vs MHC-I | 355 | +0.021 | 0.70 | 0.98 |
| TACSTD2 vs MHC-I | 355 | −0.001 | 0.98 | 0.98 |
| CLDN4 vs Immune mean-z | 355 | +0.044 | 0.41 | 0.82 |
| TACSTD2 vs Immune mean-z | 355 | **−0.133** | 0.012 | 0.041 |
| CLDN4 vs Immune ssGSEA | 355 | **−0.246** | 2.7×10⁻⁶ | — |
| TACSTD2 vs Immune ssGSEA | 355 | **−0.300** | 7.8×10⁻⁹ | — |
| CLDN4 vs Stromal mean-z | 355 | −0.118 | 0.026 | — |
| TACSTD2 vs Stromal mean-z | 355 | −0.178 | 0.00077 | — |
| CLDN4 vs TJ (no CLDN4) | 355 | **+0.740** | 7.4×10⁻⁶³ | — |
| TACSTD2 vs TJ (no CLDN4) | 355 | +0.494 | 2.9×10⁻²³ | — |
| CLDN4 vs TACSTD2 | 355 | **+0.610** | 1.6×10⁻³⁷ | — |

TROP2-high tracks a CD8-low / ESTIMATE-immune-low axis. CLDN4 tracks the TJ module and is CD8-low; it is near-zero vs Ayers6 / MHC-I / Immune mean-z (ssGSEA immune is negative). That is extra weight for **TROP2-high = immune-low** and **CLDN4-high = TJ barrier**, not a genome-wide IFN/MHC claim for CLDN4.

## ICI response (GEO labels exist)

Two-sided MWU on `log2(TPM+1)` or signature; AUC = P(score_R > score_NR).

| feature | median R | median NR | n | MWU p | AUC |
|---|---:|---:|---:|---:|---:|
| CLDN4 | 1.87 | 1.94 | 168 vs 187 | 0.40 | 0.47 |
| TACSTD2 | 8.28 | 8.20 | 168 vs 187 | 0.67 | 0.51 |
| CD8A | 3.21 | 2.78 | 168 vs 187 | **0.0027** | 0.59 |
| IFNG | 0.86 | 0.51 | 168 vs 187 | **0.00083** | 0.60 |
| Ayers6 | +0.088 | −0.115 | 168 vs 187 | **0.00070** | 0.60 |
| Immune mean-z | +0.130 | −0.021 | 168 vs 187 | 0.016 | 0.57 |
| MHC-I | −0.071 | −0.070 | 168 vs 187 | 0.27 | 0.47 |

CLDN4 and TACSTD2 do not separate GEO responders. CD8A / IFN-γ do, so the endpoint is not dead. Supp8 labels (171 vs 184 after the three flips) leave CLDN4 p=0.31 and TACSTD2 p=0.79.

This does not test ADC+ICI (no ADC arm). It is extra ICI-only context: CLDN4-high is not an ICI-response marker here; it still sits on the TJ / CD8-low side.

## OS / PFS

OS: Supp8 `Overall survival (days)` + `Death` (265 events / 355). Cox HR per 1 SD of the score. Median-split log-rank is also reported.

| feature | OS HR (95% CI) | Cox p | events | median-split log-rank p |
|---|---|---:|---:|---:|
| CLDN4 | **1.19 (1.05–1.36)** | **0.0081** | 265 | 0.54 |
| TACSTD2 | 1.07 (0.94–1.22) | 0.28 | 265 | 0.13 |
| CD8A | 0.84 (0.75–0.94) | 0.0027 | 265 | 0.0053 |
| Ayers6 | 0.84 (0.75–0.94) | 0.0032 | 265 | 0.0011 |

Higher CLDN4, shorter OS on the continuous Cox (median split is NS). TACSTD2 OS is NS. Immune controls go the expected direction.

**PFS event is not deposited.** Spearman vs PFS *time* (ignores censoring): CLDN4 ρ=−0.090 p=0.091; TACSTD2 ρ=−0.035 p=0.51. A reconstructed event (PFS < OS or Death; 349/355 events) is extra only: CLDN4 HR 1.13 p=0.037; TACSTD2 HR 1.06 p=0.24. Do not treat the reconstructed PFS Cox as an author endpoint.

## Histology

LUAD vs LUSC was requested if labels exist. They do not, on GEO or Supplementary Data 8. No inferred histology.

## How to read this against the given thesis

- **CLDN4-high = TJ barrier:** CLDN4 vs the other nine TJ genes ρ=+0.74 (n=355).
- **TROP2-high = immune-low:** TACSTD2 vs CD8A ρ=−0.26; vs Immune mean-z ρ=−0.13; vs Immune ssGSEA ρ=−0.30.
- **CLDN4-high may mark ADC+ICI candidates:** ICI-only cohort; CLDN4 does not predict GEO response (p=0.40) while CD8A/IFN-γ do. CLDN4-high is CD8-low and has worse continuous OS. That is compatible with a combination-candidate hypothesis; it is not an ADC trial.

## Files

- `analyze.py` — download GEO TPM + series matrix + Supp8; patient-level tests
- `SI_geneset.gmt` — Yoshihara ESTIMATE Immune/Stromal lists
- `tables/one_row.tsv`, `stats.tsv`, `per_patient.tsv`, `label_inventory.tsv`, `signature_coverage.tsv`, `summary.json`
- `figures/fig1_response_boxplots.png`
- `figures/fig2_cldn4_immune_scatter.png`
- `figures/fig3_tacstd2_immune_scatter.png`
- `figures/fig4_km_os.png`
- `figures/fig5_correlation_heatmap.png`
- `figures/fig6_spearman_forest.png`
- `figures/fig7_cldn4_vs_tacstd2.png`
