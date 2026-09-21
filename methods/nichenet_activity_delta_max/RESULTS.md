# RESULTS — NicheNet activity Δ, concordant-4

CLDN4-high versus CLDN4-low **malignant senders**, T/NK **receiver axis**.
Cohorts are only GSE123902, GSE131907, GSE205335, and GSE189357.
The locked CLDN4 percent-positive versus T/NK fraction correlation is not re-estimated.

Two contrasts were scored because they are different quantities.

- **Outgoing Δ** asks whether CLDN4-high malignant cells express the ligand more, in units that also have the NicheNet receptor on T/NK.
- **Target activity** asks whether that ligand's NicheNet v2 targets sit on the CLDN4-high or the CLDN4-low side of the T/NK axis.

The classification activity Δ (AUROC on a CLDN4-high gene list minus AUROC on a CLDN4-low gene list) is reported as well. On this axis it is large for the barrier ligands and almost as large for a typical ligand.

## Outgoing Δ is specific for the barrier ligands

Eligible units for the quartile split: **59**. Sender Δ is mean `log1p(CP10k)` in the CLDN4 top quartile minus the bottom quartile. The receptor fraction is the share of those units in which the best NicheNet v2 receptor is detected in T/NK. LR-gated Δ multiplies the two. Ranks are among **286** potential ligands (1 = highest outgoing from CLDN4-high).

| ligand | receptor | sender Δ | receptor fraction | LR-gated Δ | sender Wilcoxon p | rank |
|---|---|---:|---:|---:|---:|---:|
| F11R | ITGB2 | +0.350 | 0.983 | +0.344 | 2.52×10⁻¹¹ | 6 |
| CDH1 | ITGAE | +0.382 | 0.898 | +0.343 | 3.25×10⁻¹¹ | 7 |
| NECTIN2 | CD96 | +0.253 | 0.949 | +0.240 | 2.65×10⁻¹¹ | 18 |
| LGALS9 | CD44 | +0.052 | 1.000 | +0.052 | 0.00192 | 88 |
| **family mean** |  | **+0.259** |  | **+0.245** |  |  |

A draw of four potential ligands beats this LR-gated mean in 8 of 4,000 permutations (one-sided p = 0.0022). The same test on sender Δ alone, without the receptor weight, is p = 0.013. LGALS9 is the weak member on expression. F11R, CDH1, and NECTIN2 are not the top three ligands in the matrix: LGALS3BP, CD55, and CEACAM6 also have large sender deltas. The pre-specified set is enriched. It is not a unique surface ranking.

## IFN / recruit in CLDN4-low splits in two

**Target activity, no expression filter.** Top 50 prior targets inside the 6,398-gene T/NK universe. Target AUROC above 0.5 means those targets sit toward CLDN4-high. The universe mean meta-z is −0.097, so centered z is the mean of the 50 targets minus that background.

| ligand | target AUROC | centered z | rank of AUROC toward CLDN4-low, of 1,226 ligands | potential ligand | sender Δ |
|---|---:|---:|---:|---|---:|
| IFNG | 0.346 | −0.080 | 6 | no | −0.032 |
| CXCL9 | 0.579 | +0.052 | 345 | no | −0.000 |
| CXCL10 | 0.611 | +0.070 | 650 | no | −0.005 |
| CCL5 | 0.555 | +0.025 | 230 | yes | −0.121 |

IFNG is the IFN ligand whose targets sit in CLDN4-low T/NK cells. The six lowest target AUROCs in the whole prior are IFNB1 (0.281), IFNA1 (0.289), IFNL3 (0.317), IFNL1 (0.322), IL27 (0.331), and IFNG (0.346). That tail is the interferon prior, not a random low rank. CXCL9 and CXCL10 do not join it.

IFNG, CXCL9, and CXCL10 are not potential ligands. They are detected in at least 10% of malignant cells in fewer than 10% of eligible units. Their sender deltas are near zero. The target shift is not an outgoing malignant signal.

**CCL5 is the expressed recruit ligand.** Sender Δ = −0.121 (Wilcoxon p = 1.16×10⁻⁵). Among 286 potential ligands it is 6th from the bottom on the LR-gated score (CXCR3 detected in 69% of units; LR-gated Δ = −0.084; one-sided permutation p = 0.018). Its own top-50 targets are not the CLDN4-low program (target AUROC 0.555). The ligands below CCL5 on the LR-gated list are mostly leukocyte genes (LGALS1, TYROBP, CD48, CD2, PTPRC), not a second chemokine.

## Classification activity Δ is large and not specific

Same-sign genes in both author cohorts. Activity Δ = AUROC(regulatory potential, top-N CLDN4-high genes) − AUROC(regulatory potential, top-N CLDN4-low genes).

Primary N = 200. Barrier mean Δ AUROC = **+0.135** (CDH1 +0.163, F11R +0.135, LGALS9 +0.125, NECTIN2 +0.116). Permutation versus other potential ligands, p = 0.175. AUPR-corrected Δ at this N has mean **+0.019**, p = 0.42. With 200 positives in 6,398 genes the average precision stays small. That compression is why a corrected AUPR on this axis looks null even when the AUROC difference does not.

Grid maximum for the barrier mean, N = 50:

| ligand | AUROC high | AUROC low | Δ AUROC | Δ AUPR-corrected |
|---|---:|---:|---:|---:|
| CDH1 | 0.665 | 0.350 | +0.315 | +0.011 |
| F11R | 0.665 | 0.377 | +0.287 | +0.011 |
| LGALS9 | 0.637 | 0.351 | +0.286 | +0.008 |
| NECTIN2 | 0.639 | 0.386 | +0.253 | +0.009 |
| **mean** |  |  | **+0.285** |  |

At that same N the median potential ligand has Δ AUROC **+0.273**. Excess over the median is +0.013. Permutation p = 0.135 (`specificity_by_n.tsv`). N = 100, 200, and 400 give the same pattern (excess over the median +0.013, +0.009, and +0.007; permutation p = 0.124, 0.171, and 0.155). The family-table test at the primary N = 200, drawn separately, is p = 0.175. The CLDN4-high gene list is easier for the prior to rank than the CLDN4-low list, for barrier ligands and for most other ligands. The grid was not used to drop a cohort or to add a discordant accession.

IFN-family mean Δ AUROC is also positive at every N (+0.235 at N = 50, +0.086 at N = 200). It does not flip to the CLDN4-low list. The CLDN4-low direction for IFN is the IFNG target-set result above, not this geneset difference.

## Joint score

Realized activity Δ = sender Δ × (target AUROC − 0.5) × 2, K = 50.

Barrier mean **+0.067**, permutation p = 0.0090. CDH1 is +0.116 (rank 5 of 286). F11R is +0.079 (rank 15). NECTIN2 is +0.061 (rank 25). LGALS9 is +0.012 (rank 109). CEACAM6, CD55, CTSD, and LGALS3BP rank above CDH1. The joint score moves the pre-specified set because sender Δ and a mild target shift point the same way. It does not make F11R or NECTIN2 the top ligand.

Centered target z, which does not use sender expression, is positive for all four barrier ligands at K = 50 (family mean +0.074) and the permutation p is 0.13. The grid maximum of barrier centered z minus IFN centered z is at K = 25 (contrast +0.075). CXCL9, CXCL10, and CCL5 have positive centered z, so they pull the IFN family mean up. IFNG alone is the negative member (centered z −0.080).

## Patient-level T/NK genes, not contact

On the same author-T/NK axis, GZMB meta-z = −0.37 (p = 0.030), PRF1 = −0.34 (p = 0.045), NKG7 = −0.41 (p = 0.016), IFNG = −0.22 (p = 0.35, signs disagree across the two cohorts). These are correlations of a T/NK pseudobulk with malignant CLDN4 percent across patients. They are not a measurement of effector cells that touch a CLDN4-high tumor cell. The spatial result stays the CosMx one: exclusion, not a seal on GZMB/PRF1/NKG7/IFNG at the interface.

ITGB2, the F11R receptor used in the LR gate, has meta-z −0.46 (p = 0.0076) inside T/NK. The gate above is a detection rate across units, not a claim that ITGB2 rises in CLDN4-high patients.

## What this does not claim

- It does not re-estimate the locked concordant-4 ρ.
- It does not treat a NicheNet target score as spatial exclusion.
- It does not score CXCL9, CXCL10, or IFNG as malignant outgoing ligands. They fail the expression filter. IFNG's target rank is reported separately and labeled as prior-target alignment.
- It does not drop a cohort to enlarge the barrier Δ. Every N in {50, 100, 200, 400} is in `specificity_by_n.tsv`.
- Cell-pooled p-values are not reported. The ligand–target matrix was not fit on these four cohorts.
