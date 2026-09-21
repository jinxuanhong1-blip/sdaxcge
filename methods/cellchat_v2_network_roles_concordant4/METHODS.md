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
PVRL2, JAM1, and CD45 in a matrix are renamed to the CellChat symbols NECTIN2, F11R, and PTPRC.
Complex receptors (ITGAL_ITGB2, ITGAE_ITGB7) stay complexes; only their subunits are read from the matrix.

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

## Max-effect sweep

`scripts/run_max_effect_sweep.R` recomputes the same extracted outgoing probabilities under a grid. It does not replace the primary Q4/Q1 table above.

The probability engine is the CellChat Hill model (`Kh = 0.5`, `n = 1`) on full-library log-normalized expression. The first scored Q4/Q1 T/NK unit is checked against `computeCommunProb`. Group fractions for `population.size = TRUE` are applied as CellChat does: probability × (n_sender / N) × (n_receiver / N), with N the three groups after the cap of 200.

Axes:

- CLDN4 split inside malignant cells: Q4 vs Q1; top vs bottom 10%, 20%, 25%, 30% (`rank(..., ties.method = "first")` on log1p CP10k); detected vs undetected (raw count > 0 vs 0).
- Receiver: T/NK, or CD8 only. GSE205335 CD8 is `lineage.sub == "CD8+ T cells"`. GSE131907 CD8 is author subtype in {CD8 low T, Cytotoxic CD8+ T, Exhausted CD8+ T, Naive CD8+ T}. Marker cohorts require CD8A > 0, CD3D or CD3E > 0, and not malignant.
- Average: `truncatedMean` trim 0.1 and 0.05, `triMean`, `thresholdedMean` trim 0.1 and 0.05, `median`.
- `population.size` TRUE and FALSE.
- Aggregation: sum of the seven barrier edges; mean of detected barrier edges; sum of every interaction in the JAM, NECTIN, CDH1, and GALECTIN pathways; mean of those pathways when the pathway sum is nonzero.
- Gene universe: the 14 pre-specified pairs in one object (the primary table's max-normalization), or the barrier-pathway interactions in one object.
- Unit filters: all units that clear the floors; drop a zero family Δ; at least three detected pairs; each malignant arm at least 30 or at least 50 cells before the cap; receiver at least 50; arm at least 30 and at least one detected pair.

A unit is scored only when each sender arm has at least 10 cells and the receiver has at least 20. The reported maximum has to be positive in all four cohorts and have a one-sided sign-flip p ≤ 0.05 (B = 10000, seed 3979). That p is the p of the chosen spec. The script also records how many specs were scored.

Log2 fold is log2(sum_high / sum_low) on units where both sums are positive. Log-odds is the mean logit difference on edges that are positive on both arms, with probabilities clipped to [1e-6, 1 − 1e-6]. Those scales are summarized separately from the probability Δ.

IFN/recruit stays the seven pre-specified IFN edges on the winning split, receiver, mean, and population-size setting. Those edges are not added into the barrier score.

If a winning barrier Δ is compared with total sender strength, both pieces are taken from one overexpressed full network (`thresh.p = 0.05`) on that same split, so they share one max-normalization. The matched residual is barrier_high − (barrier_low / total_low) × total_high. A sign-flip is run on that residual.

## Not done

Dual-high, GSE148071, GSE127465, CD45-only, 7-pool, cell-pooled tests, a discovery claim from the pathway ranking.
