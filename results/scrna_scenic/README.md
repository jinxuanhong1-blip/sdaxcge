# results/scrna_scenic — GSE207422 AUCell ρ / n

**Overall verdict: `PARTIAL`.**

A10 ELF3–CLDN4 **bulk RNA** is taken as given. This folder tests a different claim: whether **ELF3 / GRHL1 / KLF4 regulon AUCell** tracks **TACSTD2 / CLDN4** in public lung epithelium.

No ChIP peaks were created. Public priors do **not** list TACSTD2 or CLDN4 as targets of these TFs. pySCENIC 0.12.1 imports; cisTarget was not run; GRNBoost2 failed (dask `diagnostics_port` incompatibility).

## Primary table (sample means, malignant-like, min 20 cells)

| split | regulon | TACSTD2 ρ (p) | CLDN4 ρ (p) | n |
| --- | --- | --- | --- | ---: |
| all | ELF3_prior | 0.084 (0.80) | 0.007 (0.98) | 12 |
| all | ELF3_pearson | **0.825** (9.5e-4) | **0.937** (7.0e-6) | 12 |
| all | GRHL1_prior | −0.480 (0.11) | −0.393 (0.21) | 12 |
| all | GRHL1_pearson | **0.853** (4.2e-4) | **0.972** (1.3e-7) | 12 |
| all | KLF4_prior | **0.601** (0.039) | **0.671** (0.017) | 12 |
| all | KLF4_pearson | **0.671** (0.017) | **0.839** (6.4e-4) | 12 |
| all | ELF3+GRHL1+KLF4 prior union | **0.734** (0.0065) | **0.720** (0.0082) | 12 |
| all | ELF3+GRHL1+KLF4 pearson union | **0.755** (0.0045) | **0.895** (8.4e-5) | 12 |
| all | ELF3_RNA (not a regulon) | **0.832** (7.9e-4) | **0.797** (0.0019) | 12 |
| LUAD | ELF3_prior | −0.257 (0.62) | −0.143 (0.79) | 6 |
| LUAD | ELF3_pearson | 0.600 (0.21) | 0.771 (0.072) | 6 |
| LUAD | GRHL1_pearson | 0.714 (0.11) | **0.829** (0.042) | 6 |
| LUAD | KLF4_prior | 0.143 (0.79) | 0.429 (0.40) | 6 |
| LUAD | KLF4_pearson | 0.029 (0.96) | 0.314 (0.54) | 6 |
| LUAD | ELF3+GRHL1+KLF4 prior union | 0.429 (0.40) | 0.657 (0.16) | 6 |
| LUAD | ELF3_RNA (not a regulon) | **0.886** (0.019) | **0.943** (0.0048) | 6 |

Bold = pre-specified SUPPORT (ρ ≥ 0.30 and p < 0.05). GRHL1_prior is 1 gene (`IGF2`); treat as unspecified.

## What this supports / does not

- **Does not support** “the published ELF3 target list is active where TACSTD2/CLDN4 are high.” ELF3_prior AUCell is **null** at the sample level (all n=12 ρ≈0; LUAD n=6 slightly negative). ELF3 prior ∩ Pearson top-50 = **0 genes**.
- **Does support**, in the **pooled Adeno+Squamous** split only, that **Pearson co-expression neighborhoods** of ELF3/GRHL1/KLF4 (TACSTD2/CLDN4 held out) track both genes. That is the same epithelial program, not binding.
- **KLF4 public prior** (112 genes, includes EPCAM/KRT19/CDH1, not TACSTD2/CLDN4) meets SUPPORT in the pooled split and fails in LUAD-only.
- **LUAD-only n=6** is fragile. ELF3 **RNA** still tracks both genes (A10, given). Regulon AUCell does not meet SUPPORT except GRHL1_pearson vs CLDN4.
- **Cell-level** n=9782 (all) / 2360 (LUAD) makes every weak ρ “significant.” Example: ELF3_prior vs CLDN4 cell ρ=0.11, p=1e-28, sample ρ=0.007, n=12. Use the sample column.

## Cohort

GSE207422 Hu et al. *Genome Med* 2023 (PMID 36869384). 92,330 cells; 13,043 epithelial; 9,782 malignant-like (marker proxy; CopyKAT labels not on GEO). 12 samples with ≥20 malignant-like cells (6 LUAD, 6 LUSC). LUSC sample BD_immune07 contributes 4,483 of those cells — the pooled split is not a balanced LUAD/LUSC mix.

## Files

| file | content |
| --- | --- |
| `summary.json` / `verdict.json` | machine-readable call |
| `correlations.tsv` | every unit × split × compartment × regulon × gene |
| `regulons.tsv` | gene sets (held-out TACSTD2/CLDN4) |
| `pearson_tf_targets.tsv` | full Pearson TF–gene table |
| `sample_means_*.tsv` | per-sample means used for primary ρ |
| `fig_sample_rho_forest.png` | primary forest |
| `fig_elf3_prior_vs_cldn4.png` | ELF3_prior vs CLDN4 (null) |
| `fig_cell_vs_sample_rho.png` | why cell-level n misleads |
