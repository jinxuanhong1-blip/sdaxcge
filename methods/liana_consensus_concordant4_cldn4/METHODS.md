# Methods — LIANA+ / CellPhoneDB / Connectome consensus

ADDITIVE. CLDN4 only. No dual-high. Not a spatial test.
Concordant four: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Do not add GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526.
The unit is the locked patient or sample.

## Question

Which pre-specified ligand–receptor edges from malignant CLDN4-high vs CLDN4-low
cells toward T/NK and toward myeloid agree across LIANA+, CellPhoneDB, and Connectome,
and which of those edges support a barrier/exclusion program versus a recruitment program?

## Engine

Python `liana` 1.10 `rank_aggregate` (Dimitrov et al., Nat Commun 2022), one call per
unit per receiver. That call runs CellPhoneDB, Connectome, log2FC, NATMI, and
SingleCellSignalR, then RobustRankAggregate.

Reported scores, all taken from that call:

| name | column | contrast |
|---|---|---|
| CellPhoneDB | `lr_means` | high − low |
| Connectome | `expr_prod` | high − low |
| LIANA+ | `magnitude_rank` (ρ) | ρ(low) − ρ(high) |

Positive Δ means the edge is stronger from CLDN4-high. ρ is smaller for stronger
consensus magnitude, so the LIANA+ Δ is flipped to the same sign.

LIANA+ ρ already includes the CellPhoneDB and Connectome magnitudes. The three
columns are corroboration, not three independent experiments. NATMI’s magnitude
is also an expression product, so the aggregate is weighted toward that product.

Settings: `expr_prop=0.1`, `min_cells=5`, `n_perms=1000`, `seed=1337`,
`aggregate_method='rra'`, `return_all_lrs=True`. Each group is capped at 400 cells
after the cell-number gate, with one Generator(1337) walked in cohort then patient order.
Expression fed to LIANA is log1p counts-per-10k using the **full-matrix** library size.
The earlier CellChat run ranked CLDN4 on a ligand-gene-subset library; the quartile
rule is the same (ordinal rank, low = rank ≤ floor(n/4), high = rank > ceil(3n/4)),
the library size is not.

GSE205335 has no MTX on GEO. `export_gse205335.R` peels the double-gzipped RDS
(dgCMatrix, 33714 × 96505) in R and writes a gene-subset Matrix Market.

## Labels

Same malignant and T/NK rules as the CellChat concordant-4 run. Myeloid is new.

- GSE123902 and GSE189357: malignant = (EPCAM or KRT8 or KRT18 or KRT19) > 0 and PTPRC = 0.
  T/NK = (CD3D or CD3E or CD8A or NKG7 or GNLY or KLRD1) > 0 and not malignant.
  Myeloid = (LYZ or CD68 or CD14 or FCGR3A or CSF1R or AIF1) > 0 and not malignant or T/NK.
  GSE123902 uses the donor tumour file (PRIMARY over METASTASIS). Normal libraries are out.
- GSE131907: author subtype in {Malignant cells, tS1, tS2, tS3};
  T/NK = Cell_type in {T lymphocytes, NK cells};
  myeloid = Cell_type == Myeloid cells. Locked samples with malignant cells only.
- GSE205335: author `lineage.sub` == Malignant cells; `lineage.total` == T/NK cells or Myeloid cells.
  Tissue starting with "Normal" is out. Locked patients only.

The script stops if n_mal or n_tnk disagrees with the CellChat inventory on any locked unit.
TACSTD2 is never a gate. Aliases applied only when the official symbol is absent:
PVRL2→NECTIN2, JAM1→F11R, IL8→CXCL8, PDL1/PDCD1LG1→CD274, PD1→PDCD1, CD155→PVR, CD112→NECTIN2.

## Gates

Primary split is within-unit malignant Q4 vs Q1. A unit is scored for a receiver when
n_mal ≥ 40, both arms ≥ 10, and that receiver ≥ 20. T/NK and myeloid are separate runs
so the rank universe is not mixed across receivers.

## Edge panel

Not a discovery screen. Edges are listed in `data/edges.tsv`.

1. Barrier / exclusion, expect CLDN4-high > low. Junction, checkpoint, don’t-eat-me, and
   suppressive ligands: F11R–F11R, F11R–ITGAL_ITGB2, NECTIN2/PVR–TIGIT/CD96, CDH1–CDH1,
   CDH1–ITGAE_ITGB7, CDH1–KLRG1, LGALS9–HAVCR2/CD44/PTPRC, CD274–PDCD1, PDCD1LG2–PDCD1,
   HLA-E–KLRC1_KLRD1, CD47–SIRPA, CD24–SIGLEC10, TGFB1–TGFBR1_TGFBR2, MIF–CD74.
2. Effector recruitment, expect CLDN4-low > high (KD-like arm): CXCL9/10/11–CXCR3,
   CXCL16–CXCR6, CCL5–CCR5, CCL5–CCR1, CX3CL1–CX3CR1.
3. Myeloid recruitment, no directional thesis: CCL2–CCR2, CCL7–CCR2, CXCL8–CXCR1,
   CXCL8–CXCR2, CXCL1–CXCR2, CXCL2–CXCR2, CCL3–CCR1, CSF1–CSF1R.
4. Scored but not used as barrier or recruitment support: HLA-A/B/C–CD8A (continuity
   with the CellChat MHC arm) and CXCL12–CXCR4 (ambiguous retention vs recruitment).

A missing arm is 0 for `lr_means` and `expr_prod` when the other arm is finite.
If both arms are missing, Δ is NA. LIANA+ Δ is NA unless both ρ values are finite.
Complexes use LIANA’s underscore syntax; a missing subunit drops the edge.

## Tests and BH

Three BH families, not one pooled family:

1. **Edges.** For each method and receiver, Wilcoxon signed-rank of Δ vs 0 across units
   with a finite Δ. Tested when n ≥ 8. BH across the pre-specified panel (all classes,
   including HLA–CD8 and CXCL12–CXCR4) within that method × receiver.
2. **Families.** Within a unit, the family score is the mean of finite edge Δs in that
   class. Wilcoxon across units. BH across the 18 tests (3 classes × 2 receivers × 3 methods).
   Myeloid recruitment has no expected sign. Barrier expects high>low. Effector recruitment
   expects low>high.
3. **CellPhoneDB permutations.** Within each unit × sender (high or low) × receiver, BH
   across edges with a finite permutation p. The summary is the fraction of units with
   q<0.05. That p tests specificity against a label shuffle, not the high-vs-low contrast.

Consensus call (a hit): at least two of the three methods have BH q<0.05 with the same
sign, and none of the three is BH-significant in the opposite direction.

Support:

- Barrier hit with high>low → supports barrier/exclusion from CLDN4-high.
- Effector-recruitment hit with low>high → supports recruitment from CLDN4-low.
- Effector-recruitment hit with high>low → recruitment from CLDN4-high (opposite the KD-like arm).
- Myeloid-recruitment hit → supports myeloid recruitment from the arm that is higher.

Ranking uses a Stouffer combination of the signed Wilcoxon z values. That z orders edges.
It is not reported as a calibrated meta-analytic p-value, because the methods share cells
and LIANA+ already aggregates the other two magnitudes.

Zeros are dropped inside the Wilcoxon, matching `wilcox.test` / `zero_method='wilcox'`.

## Not done

Dual-high, non-concordant cohorts, cell-pooled tests as the honest n, spatial neighbours,
Visium same-spot correlation, and an unconstrained OmniPath screen.
