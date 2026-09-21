# Methods — CellChat v2 full network and roles

ADDITIVE. Thesis already correct: CLDN4-high malignant cells are the barrier senders.
CLDN4 only. No dual-high. Concordant four only
(GSE123902 + GSE131907 + GSE205335 + GSE189357).
Do not add GSE148071 / GSE127465 / CD45-only.
Patient / locked sample is the unit. This does not replace the locked 14-pair table.

## Engine

R + Seurat + CellChat (`jinworks/CellChat`).
`computeCommunProb` uses `type = "truncatedMean"`, `trim = 0.1`, `population.size = TRUE`.
Python is used only to stream the GSE131907 genes-by-cells text into a slim matrix.
If Seurat or CellChat cannot install, the analysis stops.

## Matrices

Same public files as the locked concordant-four CellChat run.

- GSE123902 dense UMI CSV. Primary preferred over metastasis. Normal skipped.
- GSE131907 raw UMI text + author annotation. Locked samples with malignant cells only.
- GSE205335 UMI RDS (double-gzip peeled) + CellIdentity + GSM map. Normal tissue excluded.
- GSE189357 10x MTX.

## Malignant / T/NK and the CLDN4 split

- GSE123902 and GSE189357: epithelium marker-malignant
  (EPCAM or KRT8 or KRT18 or KRT19 count > 0, and PTPRC == 0).
  T/NK = CD3D or CD3E or CD8A or NKG7 or GNLY or KLRD1, and not malignant.
- GSE131907: author subtype in {Malignant cells, tS1, tS2, tS3};
  T/NK = author cell type in {T lymphocytes, NK cells}.
- GSE205335: author `lineage.sub` malignant; T/NK = author `lineage.total`.

TACSTD2 is never a gate.

CLDN4 rank inside malignant cells is log1p(count / full UMI library size × 10,000).
Primary split is Q4 vs Q1 (`rank(..., ties.method = "first")`).
A unit enters the test only if malignant cells ≥ 40, each arm ≥ 10, and T/NK ≥ 20.
Each identity is then capped at 200 cells so the permutation loop stays on the same cell budget.
The cap is a patient-specific random sample. Population-size scaling uses the capped counts.

CellChat sees Seurat LogNormalize data (scale 10,000) computed with that same full library size,
not the column sums of a ligand-only submatrix.

## Full network

`CellChatDB.human`, restricted to Secreted Signaling, ECM-Receptor, and Cell-Cell Contact.
Any pre-specified pair that carries a different annotation label is kept.

Per unit:

1. `identifyOverExpressedGenes` (Wilcoxon, `thresh.p = 0.05`) and `identifyOverExpressedInteractions`.
2. `computeCommunProb` on that overexpressed set (`nboot = 20`).
3. Pathway probability = sum of L-R probabilities in the pathway (the continuous network, not zeroed by the p-value gate).
4. `netAnalysis_computeCentrality` on each pathway network:
   out-degree (sender), in-degree (receiver), flow betweenness (mediator), information centrality (influencer).
5. A global control sums every overexpressed interaction from the sender to T/NK, and the sender out-degree on that whole network.

With three nodes, mediator and influencer are secondary. The sender role is out-degree.

## Extracted outgoing

Pre-specified, not a screen.

Barrier / inhibitory (expect high > low):

- JAM1–ITGAL/ITGB2 (F11R), NECTIN2–TIGIT, CDH1–ITGAE/ITGB7, CDH1–KLRG1,
  LGALS9–HAVCR2, LGALS9–CD44, LGALS9–PTPRC/CD45.

IFN / recruit (expect low > high; the KD-like arm):

- CXCL9–CXCR3, CXCL10–CXCR3, CCL5–CCR5, CCL5–CCR1, HLA-A/B/C–CD8A.

These pairs are scored even when they are not in the overexpressed set, so the thesis ligands are not dropped by the discovery filter.
`nboot = 100` is the CellChat label-shuffle p-value for each directed edge (high→T/NK and low→T/NK).

Family outgoing `delta_sum` = sum of high→T/NK probabilities minus sum of low→T/NK probabilities.
Undetected pairs contribute 0.
`mean_detected_delta` averages only pairs with a positive probability on at least one arm.
Family sender role = out-degree of CLDN4-high minus out-degree of CLDN4-low on the subnetwork built from that family's pairs.

## Permutation tests

1. **Edge permutation.** CellChat `nboot = 100` shuffles cell labels and asks whether each extracted edge is larger than that null.
2. **Within-unit sender permutation.** Among malignant cells only, shuffle the Q4/Q1 labels, keep T/NK fixed, and recompute the family `delta_sum` and the family out-degree with CellChat (`B = 49`). One-sided p-values follow the thesis direction. Stouffer combines units. This is the direct high-vs-low test.
3. **Across-unit sign-flip.** Flip the sign of each patient delta (`B = 10000`, seed 3979). One-sided p is the thesis test. Wilcoxon signed-rank vs 0 is reported beside it.

The pathway table's BH q-values are exploratory. The thesis tests are the two families.

## Not done

Dual-high, GSE148071, GSE127465, CD45-only, 7-pool, cell-pooled tests, a discovery claim from the pathway ranking.
