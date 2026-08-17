# GSE218989 — CLDN4 vs CXCL9 / CXCL10 / CXCL13 and GEP-like

**Additive only.** User thesis (given, not re-audited): CLDN4-high = tight-junction (TJ) barrier; CLDN4-high may mark ADC+ICI candidates; TROP2-high = immune-low.

Public SMC–KAIST PD-1/PD-L1 NSCLC RNA-seq (Kang et al., *Nat Commun* 2024, PMID 38744958; GEO [GSE218989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218989)). Unit is the **patient**. No slide was re-scored.

This folder does **not** audit or retract `methods/gse218989_cldn4_ici` or `methods/gse218989_cldn4_gsea`. **CLDN4 vs CD274 / HLA is out of scope** (separate launch). GEP18 includes CD274 as 1/18 genes; that gene is not isolated here.

## Honest n

Paper text mentions 497 ICI lung transcriptomes. Tests stay on the public TPM matrix.

| item | n | rule |
|---|---:|---|
| patients on public TPM | **355** | `GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz` columns (`SMC__Pat*`) |
| protein-coding genes | 19916 | matrix rows after first-duplicate drop |
| paper-text ICI transcriptomes | 497 | **not** the analysis n; 142 lack a public TPM column |
| GEO treatment outcome | 355 | Responder 168 / Non-responder 187 |
| histology (LUAD vs LUSC) | **0** | not a GEO characteristic; not split |
| CLDN4 | 355 | present (4 TPM=0) |
| CXCL9 | 355 | present (8 TPM=0) |
| CXCL10 | 355 | present (15 TPM=0) |
| CXCL13 | 355 | present (25 TPM=0; median TPM 0.55 — low) |
| GEP18 mean-z | 355 | **18/18** Ayers 2017 genes present; unweighted within-cohort z-mean |

Primary Spearman n=355 for every pre-specified pair. No gene was imputed. Zeros are kept (Spearman uses ranks).

## One-row table

| dataset | n | CXCL9 / CXCL10 / CXCL13 | GEP18 | CLDN4–CXCL9 ρ (p) | CLDN4–CXCL10 ρ (p) | CLDN4–CXCL13 ρ (p) | CLDN4–GEP18 ρ (p) | histology |
|---|---:|---|---|---|---|---|---|---|
| GSE218989 SMC-KAIST | 355 | all 3 on matrix | 18/18 | +0.027 (0.61) | +0.091 (0.088) | −0.010 (0.85) | +0.051 (0.34) | **not deposited** |

Full numeric row: `tables/one_row.tsv`.

Honest read: on the public 355-patient TPM, **CLDN4 does not track CXCL9, CXCL10, CXCL13, or GEP-like**. All four primary Spearman tests are NS after BH inside those four pairs. That matches the already-published ICI page (CLDN4 vs Ayers IFN-γ 6-gene ρ=+0.029, p=0.58). It is extra weight for **CLDN4-high = TJ barrier, not an IFN-chemokine / GEP-high state** in this bulk ICI matrix. It is not a CD274/HLA result.

## Locked design (before ρ)

| Item | Choice |
|---|---|
| Matrix | Public author TPM only. No FASTQ / SRA. |
| File | `GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz` |
| Transform | `log2(TPM+1)` for single genes |
| Chemokines | CXCL9, CXCL10, CXCL13 as three pre-specified single genes |
| GEP-like | Ayers *JCI* 2017 public unweighted 18-gene T-cell-inflamed GEP (`gep18_genes.tsv`); mean of gene-wise z-scores |
| Extra (not headline) | mean-z of CXCL9/CXCL10/CXCL13; GEP18 after dropping CXCL9; Ayers IFN-γ 6-gene (already on the ICI page) |
| Test | Spearman; BH *q* inside the four primary CLDN4 pairs (3 chemokines + GEP18) |
| Companion | TACSTD2 vs the same axes (same-run only) |
| Out of scope | Isolated CD274 / HLA-A/B/C; headline ICI response OR; LUAD vs LUSC; ADC+ICI |

CXCL9 is in GEP18. CXCL10 is in Ayers6, **not** in GEP18. CXCL13 is in neither (TLS-associated chemokine).

## Primary Spearman (n=355)

BH *q* is inside the four primary pairs only.

| pair | n | ρ | p | q (BH, 4) |
|---|---:|---:|---:|---:|
| CLDN4 vs CXCL9 | 355 | +0.027 | 0.61 | 0.82 |
| CLDN4 vs CXCL10 | 355 | +0.091 | 0.088 | 0.35 |
| CLDN4 vs CXCL13 | 355 | −0.010 | 0.85 | 0.85 |
| CLDN4 vs GEP18 mean-z | 355 | +0.051 | 0.34 | 0.67 |

Extra / sensitivity (not in BH):

| pair | n | ρ | p |
|---|---:|---:|---:|
| CLDN4 vs CXCL9/10/13 mean-z | 355 | +0.030 | 0.58 |
| CLDN4 vs GEP18 minus CXCL9 | 355 | +0.055 | 0.30 |
| CLDN4 vs Ayers6 mean-z | 355 | +0.029 | 0.58 |
| CLDN4 vs TACSTD2 | 355 | **+0.610** | 1.6×10⁻³⁷ |

Positive-control pairs (axis is alive): CXCL9 vs CXCL10 ρ=+0.843 (p=6.9×10⁻⁹⁷); CXCL9 vs GEP18 ρ=+0.759 (p=1.3×10⁻⁶⁷).

## TACSTD2 companion (same 355; not a headline)

| pair | n | ρ | p |
|---|---:|---:|---:|
| TACSTD2 vs CXCL9 | 355 | −0.006 | 0.91 |
| TACSTD2 vs CXCL10 | 355 | +0.078 | 0.14 |
| TACSTD2 vs CXCL13 | 355 | −0.002 | 0.96 |
| TACSTD2 vs GEP18 mean-z | 355 | −0.071 | 0.18 |

TROP2-high = immune-low on this cohort was the CD8A / ESTIMATE axis (`methods/gse218989_cldn4_ici`: TACSTD2–CD8A ρ=−0.262). It is **not** a CXCL9/10/13 or GEP18 anti-correlation here.

## CLDN4 Q4 vs Q1 (extra; 89 vs 89)

Quartiles of CLDN4 `log2(TPM+1)`: Q1 n=89, Q2 n=89, Q3 n=88, Q4 n=89. Two-sided MWU on the chemokine / GEP scores (Q4 vs Q1 only):

| feature | median Q4 | median Q1 | n | MWU p | AUC |
|---|---:|---:|---:|---:|---:|
| CXCL9 | 4.59 | 4.79 | 89 vs 89 | 0.69 | 0.52 |
| CXCL10 | 5.70 | 5.32 | 89 vs 89 | 0.16 | 0.56 |
| CXCL13 | 0.68 | 0.77 | 89 vs 89 | 0.60 | 0.48 |
| GEP18 mean-z | +0.043 | +0.048 | 89 vs 89 | 0.63 | 0.52 |

Same null as the continuous Spearman. AUC = P(score_Q4 > score_Q1).

## GEO response — already known for CLDN4; chemokine/GEP note only

CLDN4 vs GEO ICI response is already NS in `methods/gse218989_cldn4_ici` (MWU p=0.40, AUC=0.47; 168 R / 187 NR). CD8A / Ayers IFN-γ do separate responders, so the endpoint is not dead.

Same GEO labels, chemokine / GEP as **positive-control note** (not a new CLDN4-response claim):

| feature | median R | median NR | n | MWU p | AUC |
|---|---:|---:|---:|---:|---:|
| CXCL9 | 4.78 | 4.36 | 168 vs 187 | 0.0074 | 0.58 |
| CXCL10 | 5.76 | 5.10 | 168 vs 187 | 0.00088 | 0.60 |
| CXCL13 | 0.69 | 0.48 | 168 vs 187 | 0.013 | 0.58 |
| GEP18 mean-z | +0.124 | −0.115 | 168 vs 187 | 0.00086 | 0.60 |
| CLDN4 | 1.87 | 1.94 | 168 vs 187 | 0.40 | 0.47 |

The chemokine / GEP axis separates GEO responders. CLDN4 does not sit on that axis.

## How to read this against the given thesis

- **CLDN4-high = TJ barrier, not IFN-chemokine-high:** CLDN4 vs CXCL9/10/13 and vs GEP18 are all near zero (n=355). The ICI page already has CLDN4 vs the other nine TJ genes ρ=+0.74.
- **TROP2-high = immune-low:** not supported on CXCL9/10/13 or GEP18 here; the CD8A / ESTIMATE inverse remains on the ICI page.
- **CLDN4-high may mark ADC+ICI candidates:** ICI-only cohort. This page does not re-test response. Chemokine / GEP track response; CLDN4 does not track those chemokines.

## What this is not

- Not a re-run or retraction of `methods/gse218989_cldn4_ici` or `methods/gse218989_cldn4_gsea`.
- Not CLDN4 vs CD274 / HLA (separate launch).
- Not a 497-patient analysis. Public TPM n=355.
- Not a LUAD vs LUSC split — histology is not deposited.
- Not an ADC+ICI trial. ICI-only cohort.
- Not a tumour-cell-intrinsic chemokine call. Bulk ICI RNA.
- Not a claim that CXCL13 is well measured: median TPM 0.55, 25 zeros.

## Files

- `analyze.py` — download GEO TPM + series matrix; patient-level Spearman
- `gep18_genes.tsv` — Ayers 2017 18-gene list
- `tables/one_row.tsv`, `stats.tsv`, `per_patient.tsv`, `label_inventory.tsv`, `signature_coverage.tsv`, `gene_presence.tsv`, `quartile_n.tsv`, `quartile_mwu.tsv`, `response_mwu_note.tsv`, `summary.json`
- `figures/fig1_cldn4_chemokine_scatter.png`
- `figures/fig2_spearman_forest.png`
