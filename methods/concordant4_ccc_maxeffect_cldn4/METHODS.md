# Methods — barrier ligand communication, maximized absolute Δ

ADDITIVE. Concordant four only: GSE123902, GSE131907, GSE205335, GSE189357.
CLDN4 only. No dual-high. No TACSTD2 gate.
Do not add GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526.
The patient (or locked sample) is the unit.

This does not replace PR 616. That run reported a CellChat family sum of +0.0037
with `population.size = TRUE`. Population size multiplies each edge by the sender
proportion times the receiver proportion, so the absolute probability is small
even when the ligand gap is not.

## Question

Which standard communication score gives the largest absolute high−low delta for
F11R, NECTIN2, CDH1, and LGALS9 from malignant CLDN4-high cells to T/NK, with
the same sign in every cohort?

IFN/recruit chemokines (CXCL9, CXCL10, CXCL11, CCL5, CXCL16) are a separate family.
HLA-A/B/C–CD8A is a third family. It is antigen presentation, not recruitment.

## Cells

- GSE123902 and GSE189357: epithelium marker-malignant
  (EPCAM, KRT8, KRT18, or KRT19 count > 0, and PTPRC = 0).
  T/NK = CD3D, CD3E, CD8A, NKG7, GNLY, or KLRD1, and not malignant.
  GSE123902 uses the primary tumor when a donor also has a metastasis. Normal is out.
- GSE131907: author subtype in {Malignant cells, tS1, tS2, tS3}.
  T/NK = author cell type in {T lymphocytes, NK cells}.
- GSE205335: author `lineage.sub` malignant. T/NK = author `lineage.total`.
  Normal tissue is out.

CLDN4 rank is log1p(count / full UMI library size × 10,000) inside malignant cells.
Primary split is Q4 vs Q1 by ordinal rank (ties broken by position), matching the
CellChat run: low = rank ≤ floor(n/4), high = rank > ceil(0.75 n).
A unit enters when malignant cells ≥ 40, each arm ≥ 10, and T/NK ≥ 20.
Decile (outer 10%) is the same rule with arm floor 10, so n_mal ≥ 100.

Library size is the full transcriptome, not the ligand submatrix.
PVRL2 and JAM1 are renamed to NECTIN2 and F11R when the official symbol is absent.
Expression is not capped. There is no cell downsample.

## Scores

Every edge is zero unless the receptor complex is in ≥10% of T/NK cells and the
ligand is in ≥10% of at least one malignant arm. A complex receptor (ITGAL_ITGB2,
ITGAE_ITGB7) uses the geometric mean of subunit means and the minimum subunit proportion.

Ligand Δ is the mean of that ligand's detected edges. Undetected ligands contribute 0.
Family sum adds the four barrier ligands once each. It does not add the same ligand
once per receptor.

| method | definition |
|---|---|
| cellchat_hill | 10% trimmed log1p mean. P = [L/(0.5+L)] [R/(0.5+R)]. No population-size weight. |
| cellchat_hill_pop | Same Hill probability times (n_sender/N) (n_TNK/N), N = n_high + n_low + n_TNK. |
| cellchat_prod | 10% trimmed log1p mean. P = (L R)/(0.5 + L R). |
| cpdb_means | CellPhoneDB `lr_means`: arithmetic mean of log1p, (L+R)/2. |
| conn_prod | Connectome product of arithmetic log1p means. |
| natmi_spec | NATMI edge specificity across high, low, and T/NK. |
| sca_lrscore | SingleCellSignalR LRscore. μ is the mean of the log1p block. |
| liana_logfc | (mean log1p high − mean log1p low) / ln(2). Receptor term cancels. |
| cp10k_log2fc | log2((mean CP10k high + 1) / (mean CP10k low + 1)). |
| cp10k_log2fc_trim | Same with a 10% trimmed mean. |
| cp10k_delta | mean CP10k high − mean CP10k low. |
| cp10k_delta_trim | 10% trimmed-mean CP10k difference. |
| expr_prop_pp | 100 × (fraction ligand-positive high − fraction ligand-positive low). |
| conn_z | Connectome-style z-score product of the three group means. |

The seven-pair sum (F11R–LFA1, NECTIN2–TIGIT, CDH1–integrin, CDH1–KLRG1,
LGALS9–HAVCR2, LGALS9–CD44, LGALS9–CD45) is reported only for the Hill scores,
so it can be set next to the PR 616 family sum. CDH1 and LGALS9 are counted once
in the primary family sum and more than once in that seven-pair sum.

## Winner rule

Fixed before ranking:

1. Barrier family sum mean > 0.
2. Wilcoxon signed-rank p < 0.05 across patients (two-sided).
3. All four cohorts have n ≥ 3 and a positive mean.
4. Barrier per-ligand mean is larger than the chemokine per-ligand mean on the same score.
5. Untrimmed CP10k scores are dropped when the trimmed sibling is not positive and at least half as large.
6. Prefer scores where each of the four ligands is positive in all four cohorts, and the chemokine family is not itself positive in all four cohorts.
7. Among the passing scores, take the largest barrier family sum. The unit is reported with the number. Sums on different scales are not fold-changes of each other.

Sign-flip null: 10,000 flips, seed 3979, one-sided for high > low.
Cohort heterogeneity: DerSimonian–Laird I² on the four cohort means.

## Not done

Dual-high, extra cohorts, spatial distance, cell-pooled p-values, a claim that a CP10k difference is a CellChat probability.
