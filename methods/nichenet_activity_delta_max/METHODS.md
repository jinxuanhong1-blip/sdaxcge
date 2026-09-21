# Methods

Concordant four only: GSE123902, GSE131907, GSE205335, GSE189357. GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, and GSE207422 are not included. TACSTD2 is not a gate.

## Already computed inputs

`inputs/ligand_sender_delta.tsv` is the within-unit malignant split. Malignant cells are ranked by CLDN4 `log1p(CP10k)`. High is the top quartile, low is the bottom quartile. A unit is eligible when it has at least 40 malignant cells and 20 T/NK cells, mean CLDN4 is higher in the high arm, and at least 25% of malignant cells have a CLDN4 UMI. Eligible n = 59. Sender Δ is mean log expression in the high arm minus the low arm. The per-ligand Wilcoxon signed-rank p in that table is the patient-level test and is not recomputed here.

`inputs/tnk_gene_vs_cldn4.tsv` is the receiver axis. It uses author-annotated T/NK only (GSE131907 and GSE205335). For each gene, Spearman correlation of the T/NK pseudobulk with malignant CLDN4 percent positive, Fisher-z combined across the two cohorts. Marker-gated T/NK in GSE123902 and GSE189357 is not used to choose the axis.

`inputs/background_genes.tsv` and `inputs/potential_ligands.tsv` use the detection rule from that extract. Background genes are detected in at least 10% of T/NK cells in at least 10% of eligible units. Potential ligands are NicheNet v2 ligands detected in at least 10% of malignant cells in at least 10% of eligible units, with at least one v2 receptor detected in T/NK at the same threshold (286 ligands).

This run does not re-estimate the locked concordant-4 correlation of malignant CLDN4 percent positive with the T/NK fraction.

## Prior

Human NicheNet v2 `ligand_target_matrix_nsga2r_final.rds` and `lr_network_human_21122021.rds`, Zenodo record 7074291. The matrix is 33,354 genes by 1,226 ligands. Regulatory potential ranges from 0 to 0.337.

The scored universe is the intersection of matrix rows, background genes, and receiver-axis genes that are not on the epithelial-leak list and have a finite meta-z. That is 6,398 genes. Mean meta-z on that universe is −0.097 (SD 0.168). The axis is shifted, so a raw mean z near zero is above the transcriptome, not null.

Epithelial, secretory, and ciliated markers removed from the universe: EPCAM, TACSTD2, CDH1, NKX2-1, FOXA2, GATA6, ELF3, SOX2, HNF1B, KRT5/7/8/17/18/19, CLDN3/4/7, MUC1/4/5B/16, CEACAM5/6, WFDC2, SCGB1A1, SCGB3A1/2, BPIFA1/B1, SFTPA1/A2, SFTPB/C/D, SFTA2/3, AGER, AQP1/4, FOXJ1, TMC5, ARMC3, EMP2, RAB25, GRB7, NAPSA, AGR2.

## Activity definitions

All of these are reported. The grid is not used to pick a silent winner.

1. **Classification activity Δ** (NicheNet direction). Genes with the same sign in both author cohorts are ranked by meta-z. The high set is the top N positive genes. The low set is the top N negative genes. N is 50, 100, 200, or 400. Primary N is 200. For each ligand, AUROC and average precision use regulatory potential as the score and set membership as the label. Average precision is `sklearn.metrics.average_precision_score` (not the ROCR trapezoid). AUPR-corrected is average precision minus the positive-class prevalence. Activity Δ is the high-set value minus the low-set value. Positive means the ligand's prior targets match the CLDN4-high T/NK list better than the CLDN4-low list.

2. **Target-set shift.** For K in {25, 50, 100, 200}, primary K = 50, take the K universe genes with the highest regulatory potential. `target_auroc` is the AUROC of meta-z with those K genes as the positive class. Values above 0.5 mean the targets sit toward CLDN4-high. Centered z is the mean meta-z of those targets minus the universe mean (−0.097).

3. **LR-gated outgoing Δ.** Sender Δ times the fraction of eligible units in which the ligand's best NicheNet v2 receptor is detected in T/NK. The receptor does not differ by arm. The product rescales the sender Δ by whether a receiver exists. It is not a second AUPR.

4. **Realized activity Δ.** Sender Δ × (target AUROC − 0.5) × 2. Positive when the ligand is higher in CLDN4-high cells and its targets also sit toward CLDN4-high. The product of two negatives (higher in CLDN4-low, targets toward CLDN4-low) is also positive, so the sign of the sender Δ has to be read beside it.

## Families and the null

Core barrier / inhibitory, expect higher outgoing from CLDN4-high: F11R, NECTIN2, CDH1, LGALS9.

Core IFN / recruit, expect higher in CLDN4-low: CXCL9, CXCL10, CCL5, IFNG.

The ligand-level null is a draw of the same number of ligands from the 286 potential ligands (4,000 draws, seed 5533). One-sided. CXCL9, CXCL10, and IFNG are not potential ligands, so the IFN permutation on sender-gated metrics scores CCL5 only. Target-set AUROC is still reported for all four, because that score does not require malignant expression.

The N and K grid is written in full. The largest barrier-minus-IFN contrast on that grid is labeled as a grid maximum. A permutation p at the N that was chosen by looking at the grid is descriptive.

## What is a unit

Cells are not the independent unit. Cohort means are not pooled into one cell-level p-value. The ligand–target matrix is the published prior, not a model fit on these four cohorts.
