# CLDN4 virtual knockout on concordant-4 malignant cells

Additive. The locked patient-level result stands: malignant CLDN4 percent-positive versus T/NK fraction, DerSimonian–Laird ρ = −0.531 on 65 units. This note does not re-estimate that correlation. It asks whether removing CLDN4’s outgoing edges from a malignant-cell gene regulatory network moves the same genes as the observed CLDN4-low versus CLDN4-high contrast.

Cohorts are only GSE123902, GSE131907, GSE205335, and GSE189357.

## What was run

`scTenifoldKnk` 1.1 / `scTenifoldNet` 1.4 on R 4.3.3. Each cohort is a separate raw-count matrix of QC malignant cells (`n_genes` ≥ 200, `n_UMI` ≥ 500, mitochondrial fraction < 20%), at most 100 cells per locked unit, and only units with at least 30 malignant cells. The package then builds 10 principal-component-regression networks (500 cells drawn with replacement, 3 components, absolute-edge quantile 0.9), denoises them by CP decomposition (rank 3), zeros the CLDN4 outgoing row, and ranks genes by manifold distance. `distance`, `Z`, `FC`, and the chi-square `p.value` describe that distance. They are not expression log fold-changes. The signed network quantity is the denoised weight `WT[CLDN4, target]`.

The gene universe in each cohort is every gene detected in at least 5% of the GRN subsample, restricted to the 400 most variable non-family genes plus every IFN, MHC-I/APM, and TJ gene that passes the same detection filter, plus CLDN4. Family definitions are the locked ones: Hallmark IFN-α ∪ IFN-γ; the custom MHC-I/APM list; KEGG tight junction ∪ GOBP tight junction organization plus CDH1, VIM, and ZEB1, with CLDN4 removed. Twelve genes sit in both IFN and MHC (including B2M and the classical HLA-A/B/C genes), so those two set tests share members.

| cohort | GRN cells | genes in the knockout | FDR < 0.05 including CLDN4 | FDR < 0.05 excluding CLDN4 | CP norm explained |
|---|---:|---:|---:|---:|---:|
| GSE123902 | 1222 | 745 | 2 | 1 | 42.2% |
| GSE131907 | 2062 | 750 | 4 | 3 | 34.1% |
| GSE205335 | 2100 | 750 | 3 | 2 | 41% |
| GSE189357 | 900 | 723 | 3 | 2 | 44.3% |

CLDN4 itself is the largest shift in every cohort because its outgoing row was set to zero. It is left out of the enrichment tests and out of `knk_diffregulation_*.tsv`. Cell counts above are the GRN subsample, not the test n.

The observed contrast uses the locked units and the locked within-cohort CLDN4 percent-positive quartiles. Malignant counts, T/NK counts, and CLDN4 detection were recomputed from the public matrices and match the locked unit table (maximum absolute difference 0 on the integer counts). Pseudobulk is the sum of all malignant cells in the unit, TMM-CPM as in the locked script, then `log2(CPM+1) ~ cohort + Q4` on units with at least 30 malignant cells. That leaves 64 units; the OLS uses Q1 = 18 and Q4 = 16. Positive logFC is higher in CLDN4-high.

## Observed contrast

| family | genes present | logFC high−low | p | locked logFC | locked p |
|---|---:|---:|---:|---:|---:|
| IFN | 222 | −0.600 | 0.057 | −0.609 | 0.052 |
| MHC-I/APM | 21 | −0.793 | 0.054 | −0.803 | 0.051 |
| TJ (CLDN4 out) | 211 | +0.018 | 0.94 | +0.011 | 0.96 |

IFN and MHC-I/APM are higher in CLDN4-low, with the same gene counts and nearly the same coefficients as the locked family OLS. The TJ family score stays flat. The expected sign for TJ/barrier was positive; the point estimate is positive and the test is null, which is the locked result again.

## Knockout versus that contrast

Primary set test: Stouffer meta-Z across cohorts for genes present in at least two knockout universes (794 genes), then a two-sided Wilcoxon of family versus non-family genes. Positive rank-biserial would mean the family is more shifted than the other genes in the universe.

| family | n in meta | median meta-Z | background median | rank-biserial | Wilcoxon p | BH q (3 tests) |
|---|---:|---:|---:|---:|---:|---:|
| IFN | 195 | −1.19 | 0.55 | −0.695 | 5.1×10⁻⁴⁴ | 1.5×10⁻⁴³ |
| MHC-I/APM | 21 | −0.045 | 0.55 | −0.225 | 0.082 | 0.082 |
| TJ | 164 | −0.89 | 0.55 | −0.590 | 9.7×10⁻²⁹ | 1.5×10⁻²⁸ |

The same direction is present in every cohort (`knk_family_wilcoxon_by_cohort.tsv`). IFN and TJ genes are less shifted than the rest of the universe. MHC is in the same direction and is not significant at 0.05.

That comparison is tilted by construction. The background is the top of the variance list, while family genes were added whenever they cleared 5% detection. In each cohort only about 23–27 of ~170–200 IFN genes, and 11–17 of the TJ genes, are at least as variable as the least variable background gene. Restricting each family to that variance floor (`knk_family_wilcoxon_variance_matched.tsv`) brings the rank-biserials to about −0.16 to +0.39, and they do not agree in sign across cohorts. So the full-set deficit in Z is mostly the low-variance tail. The variable part of IFN, MHC, and TJ is not the knockout-shifted set either.

Genes with package FDR < 0.05, excluding CLDN4:

| cohort | gene | Z | observed logFC high−low | family |
|---|---|---:|---:|---|
| GSE123902 | CD74 | 3.54 | −1.20 | IFN |
| GSE131907 | MT-ND4 | 4.07 | +0.68 | — |
| GSE131907 | MT-CYB | 3.56 | +0.81 | — |
| GSE131907 | MT-CO1 | 3.55 | +0.62 | — |
| GSE205335 | MT-ND4 | 3.85 | +0.68 | — |
| GSE205335 | MT-ND1 | 3.74 | +0.95 | — |
| GSE189357 | SFTPB | 4.34 | +2.62 | — |
| GSE189357 | SFTPC | 3.72 | +1.93 | — |

The calls that survive FDR are mitochondrial genes and, in GSE189357, surfactant genes. Those genes are higher in CLDN4-high. CD74, an IFN-set gene, is the one FDR call in GSE123902 and is higher in CLDN4-low. Across cohorts the top of the meta-Z list is the same mitochondrial and surfactant group (MT-ND4, MT-ND1, MT-ND2, MT-CYB, SFTPB, SFTA2), again with positive observed logFC. B2M and HLA-B (in both IFN and MHC) and VIM (in the TJ set) also land in the top 15, with negative observed logFC and negative CLDN4 edges.

Overlap of the top decile of meta-Z with the observed contrast, on the 794-gene meta set:

| observed slice | overlap | odds ratio | Fisher p |
|---|---:|---:|---:|
| top decile of \|logFC\| | 19/80 | 3.33 | 1.2×10⁻⁴ |
| top decile up in CLDN4-high | 19/80 | 3.33 | 1.2×10⁻⁴ |
| top decile up in CLDN4-low | 11/80 | 1.49 | 0.24 |

The absolute-logFC decile and the up-in-high decile are different gene sets. Both share 19 genes with the knockout decile, so the odds ratios match.

Spearman of meta-Z with \|observed logFC\| is 0.093 (p = 0.008). The knockout magnitude only weakly tracks how large the expression difference is. The top of the shift list lines up with genes that are higher when CLDN4 is high, not with the IFN-up-in-low side.

Signed edges, separate from the distance: Spearman of the mean CLDN4→gene weight with observed logFC is 0.457 (p = 2.6×10⁻⁴²). Sign agreement is 492/781 (0.630, binomial p = 3.6×10⁻¹³). Inside families the edge is near zero on average (IFN −0.003, MHC −0.003, TJ +0.004). IFN sign agreement is 104/189 (p = 0.19) and MHC is 12/21 (p = 0.66). TJ sign agreement is 107/160 (0.669, p = 2.4×10⁻⁵), which includes both the epithelial genes that rise in CLDN4-high and genes such as VIM that fall. The TJ family score can stay flat while individual edges still point the same way as the gene-level logFC.

## How to read this

On these public malignant cells, a CLDN4 virtual knockout shifts mitochondrial and surfactant genes that are already higher in CLDN4-high cells. It does not concentrate the shift in the IFN or MHC-I/APM sets, and the TJ family score remains the flat result from the locked OLS. The outgoing weights do track the observed high-versus-low direction, which is the co-expression structure the network was built from. This is not an experimental knockdown, and it is not a substitute for the locked 65-unit T/NK correlation.

A ribosomal control was run on the same GSE131907 count matrix and the same 750-gene universe (`scTenifoldKnk` 1.1, same network settings, CP norm explained 34.1%, 3497 s). RPL13A, RPLP0, and RPS18 are outside the three families but did not enter that universe. The control is RPL17, the alphabetically first background gene matching `RPL` or `RPS`. It was not chosen from knockout scores.

On GSE131907 the family-versus-background Wilcoxon, with the knocked-out gene removed, is:

| knockout | family | n | median Z | background median | rank-biserial | p |
|---|---|---:|---:|---:|---:|---:|
| RPL17 | IFN | 189 | −0.583 | 0.403 | −0.527 | 6.0×10⁻²⁵ |
| CLDN4 | IFN | 189 | −0.720 | 0.508 | −0.570 | 6.8×10⁻²⁹ |
| RPL17 | MHC-I/APM | 21 | +0.054 | 0.403 | −0.222 | 0.087 |
| CLDN4 | MHC-I/APM | 21 | +0.278 | 0.508 | −0.108 | 0.40 |
| RPL17 | TJ | 156 | −0.566 | 0.403 | −0.430 | 3.9×10⁻¹⁵ |
| CLDN4 | TJ | 156 | −0.513 | 0.508 | −0.478 | 2.2×10⁻¹⁸ |

RPL17 has the same negative family rank-biserials. Four genes have FDR < 0.05, including RPL17; the other three are MT-ND4, MT-CYB, and MT-CO1, the same three genes called for the CLDN4 knockout in this cohort. Spearman of Z between the two knockouts, on the 748 genes that are neither gene, is 0.889 (p = 6.4×10⁻²⁵⁶). The mitochondrial shift and the full-set family deficit on this matrix are shared with a ribosomal knockout. They are not a CLDN4-specific move of the IFN set.

## Files

- `results/tables/family_q4q1_ols.tsv` — this run’s family OLS
- `results/tables/family_ko_vs_observed.tsv` — family score, knockout rank test, edge sign test
- `results/tables/knk_meta_genes.tsv` — gene-level meta-Z, mean edge, observed logFC
- `results/tables/knk_fdr05_genes.tsv` — FDR < 0.05 calls excluding CLDN4
- `results/tables/knk_diffregulation_*.tsv`, `knk_edges_*.tsv` — per-cohort scTenifoldKnk output
- `results/tables/knk_diffregulation_GSE131907_RPL17.tsv` — ribosomal control on the GSE131907 matrix
- `results/tables/comparison_tests.tsv` — edge, decile, and RPL17-versus-CLDN4 tests
- `results/tables/unit_sanity_summary.tsv` — match to the locked 65 units
- `results/figures/family_ko_vs_observed.png` — family logFC and knockout rank-biserial
