# Concordant-4 malignant pseudobulk: Q4 vs Q1 IFN / MHC / chemokine sweep

ADDITIVE. **CLDN4-only.** Cohorts are the locked concordant four:
GSE123902 + GSE131907 + GSE205335 + GSE189357. Quartiles are the locked
within-cohort malignant CLDN4 %pos labels. The unit is the patient,
donor, or sample. P4001 is already out of the count matrix. Expression
n is **64**. Q1 n=**18**, Q4 n=**16**. GSE189357 Q4 has 2 units, so that
cohort stays in the pooled model and is skipped for within-cohort ranks.
CLDN4 is removed from every set.

The question is whether another pre-specified ranking or gene set makes
interferon downregulation in CLDN4-high malignant pseudobulks larger than
the locked Hallmark IFN-gamma result (NES about −3.85), and what remains
after GSE205335 is removed.

## Reproduction

Zero-filled Q1/Q4 counts, one TMM, OLS `log2(CPM+1) ~ cohort + Q4`.
Genes ranked: **21604**. CLDN4 logFC **1.673** (t 2.67, p 0.012).

Hallmark IFN-gamma on that ranking: ES **−0.687**, NES **−3.87**.
The prior fgsea point was ES −0.687 and NES −3.85. The enrichment score
matches. The NES differs in the second decimal because this null is a
gene-set permutation (20,000), not fgsea's multilevel null. The
permutation p sits on the floor, **5.0×10⁻⁵** (1/20,001). It is not a
multilevel p. Patient bootstrap of the same rank, TMM held fixed, 95%
(195 of 200 draws): **−4.00 to −2.58**.

## The grid does not beat that NES

Eligible interferon tests on the full four: **310** (sets of 15–500 genes
after intersection with the ranked universe). Ranks were OLS t, OLS logFC,
signed −log10 p, limma-style moderated t, cohort-residualized
signal-to-noise, and cohort-residualized Wilcoxon z, each on the
zero-filled matrix and on the shared-gene within-cohort TMM matrix, plus
equal-cohort Stouffer z, mean logFC, inverse-variance logFC, and the
within-cohort maximum logFC.

The most negative interferon NES is the reproduced row:

| | |
|---|---|
| set | Hallmark IFN-gamma |
| rank | OLS t, zero-filled Q4 matrix |
| size | 198 |
| NES | **−3.87** |
| bootstrap 95% | **−4.00 to −2.58** |

The next rows are the same set on other ranks: Wilcoxon z −3.85,
signal-to-noise −3.82, moderated t −3.80. Shared-gene OLS t is −3.70.

Removing a frozen lymphocyte/myeloid identity list (GZMA, CD69, CD86,
LCP2, CSF2RB, and the other genes named in `analyze.py`) leaves Hallmark
IFN-gamma at NES **−3.68** (size 184, OLS t). Lineage genes are 0.11 of
the unfiltered leading edge. They are not the whole enrichment. The
leading edge still starts CCL5, GZMA, CSF2RB, GBP4, CD69, TAP1.

A 44-gene epithelial ISG/APM list with no chemokine ligands and no lineage
genes (STAT, IRF, MX, OAS, IFIT, GBP, HLA-A/B/C/E/F, B2M, TAP, PSMB8/9/10)
is NES **−3.49** on the same OLS t.

## logFC is larger on the tight ISG list

Signature logFC is the cohort-adjusted OLS beta of the mean
`log2(TMM-CPM+1)` of the set. It is the pre-specified program-level
logFC. Leading-edge means were not used to pick it.

| set | size | signature logFC | 95% | median gene logFC |
|---|---:|---:|---|---:|
| ISG/APM + CXCL9/10/11/13, CCL5 | 49 | **−0.971** | −1.455 to −0.480 | −0.917 |
| ISG/APM only | 44 | −0.914 | −1.378 to −0.411 | −0.820 |
| Hallmark IFN-gamma | 198 | −0.548 | −0.788 to −0.311 | −0.540 |

The most negative signature is the 49-gene ISG + CXCR3-chemokine list,
p 0.004. Mean and median gene logFC pick the same set.

The most negative non-sparse canonical gene (nonzero in at least 8 of 34
Q1/Q4 units; CCL5 is nonzero in all 34) is **CCL5, logFC −2.451**
(p 4.7×10⁻⁶). Then GBP2 −2.343, GBP5 −2.300, GBP1 −2.274, GBP4 −2.146,
CXCL10 −1.674. Two canonical genes are up: USP18 +0.244 and IRF7 +0.232.

## GSE205335 removed

The matrix is re-filtered and TMM is re-fit. Q4 vs Q1 is then 10 vs 13.

Hallmark IFN-gamma, OLS t: NES **−2.40** (prior leave-one-out −2.37).
The signature logFC falls to −0.248 (95% −0.551 to 0.043).

The most negative interferon NES on this three-cohort grid is the 44-gene
ISG + CXCR3 list, equal-cohort Stouffer z, NES **−2.95** (size 44).
The patient bootstrap that matches that rank — resample inside each
cohort's own Q1/Q4 TMM; GSE123902 n=7 and GSE131907 n=11; GSE189357 still
skipped — has median −2.48 and 95% **−3.27 to +0.57** (B=200). The
interval includes a positive NES. The point estimate is not a stable
three-cohort result.

The same set's signature logFC without GSE205335 is −0.458 (95% −0.982 to
0.078). CCL5 remains the most negative canonical gene, at **−1.712**.

## MHC and chemokine

Locked 21-gene MHC-I panel, same OLS t as the IFN reproduction: NES
**−2.71** (prior fgsea −2.75). Leading edge is TAP2, TAP1, PSMB9, HLA-E,
HLA-B, B2M, HLA-F, HLA-C. Lineage fraction 0.

The most negative set tagged MHC in the grid is KEGG antigen processing,
Stouffer z, NES **−3.21** (size 59). Its leading edge starts CD8A, KLRD1,
PSME2, TAP2, KLRC1, IFNG and also includes CD8B and CD4. That is not a
clean malignant MHC-I result. Without GSE205335 the same configuration
is NES −2.48.

The most negative chemokine set is Reactome "chemokine receptors bind
chemokines", signal-to-noise, NES **−3.26** (size 44). The leading edge
mixes ligands (CCL5, CCL3, CCL4, CXCL10, CXCL9, CXCL11) with receptors
(CXCR3, CXCR6, CCR7, CCR5). Without GSE205335: NES −2.64. The frozen
ligand-only list on OLS t is NES −2.90 (size 22), led by CCL5, XCL2,
CCL3, CCL4, CXCL10, CXCL9.

## What the maximum is

On the full concordant four, |NES| for interferon down is maximized by
the locked contrast: Hallmark IFN-gamma, cohort-adjusted t, NES −3.87,
bootstrap −4.00 to −2.58. |logFC| is maximized by the 49-gene ISG +
CXCR3 signature (−0.97) and, at gene level, by CCL5 (−2.45). Neither
tighter sets nor other ranks pushed the NES past that Hallmark row.

Dropping GSE205335 cuts the Hallmark NES to −2.40. A Stouffer NES of
−2.95 on the tight ISG list does not survive its own patient bootstrap.
FDR is Benjamini–Hochberg inside one rank, not across the 310 tests.
Normalization factors are not re-estimated inside a bootstrap.

Tables: `results/tables/nes_sweep.tsv`, `signature_logfc.tsv`,
`gene_logfc.tsv`, `headline.tsv`, `nes_bootstrap_matched.tsv`.
