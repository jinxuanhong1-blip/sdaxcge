# Methods — multi-sample NicheNet (CLDN4-only)

ADDITIVE. **CLDN4 only.** TACSTD2 is extracted as a companion column and never
defines senders or patient groups. Patient is the unit. Cells are counts.

## Software / prior

R is not available in this environment. Packages `nichenetr` and
`multinichenetr` were **not** run.

Ligand activity is the published **NicheNet-v2** human ligand–target matrix
`ligand_target_matrix_nsga2r_final.rds` and LR network
`lr_network_human_21122021` (Zenodo [10.5281/zenodo.7074291](https://doi.org/10.5281/zenodo.7074291);
Browaeys et al., *Nat Methods* 2020). RDS → parquet via Python `rdata`.

Prioritization follows the MultiNicheNet vignette (Bonte / Browaeys): equal-weight
mean of min-max scaled (1) sender ligand Δ, (2) ligand activity, (3) receptor
expression, (4) fraction of samples/patients expressing the ligand. This is a
documented reimplementation, not a call to `multinichenetr`.

## Cohorts

| Cohort | Public input | Malignant | T/NK | Patient rule |
| --- | --- | --- | --- | --- |
| GSE131907 | Kim et al. raw UMI + author annotation + series matrix | tumor-origin ∩ {Malignant cells, tS1, tS2, tS3} | T lymphocytes + NK, tumor-origin | pool tumor-origin samples per `patient_id` |
| GSE205335 | Hu/Ahn/Lee UMI RDS + CellIdentity + SOFT | tumor sample ∩ `lineage.sub == Malignant cells` | `lineage.total == T/NK cells` on tumor samples | pool tumor samples; drop Normal-* |
| GSE207422 | Hu et al. UMI + DRMref labels | `celltype == Malignant cells` | CD4+T + CD8+T + NK | post-treatment only |

Skipped: EGA raw FASTQ, GSE131907 2.86 GB log2TPM, GSE253013 9.3 GB RDS.

## Normalization and CLDN4 split

Gene score = `log1p(1e4 * UMI / cell total)`. **High** = CLDN4 ≥ the
**cohort-wide malignant median**; **low** = below. Same rule as the CLDN4
LIANA folders. No dual-high gate.

## Two n's (both honest)

1. **Program n** — patient has ≥20 malignant and ≥20 T/NK. Used for Spearman
   of malignant CLDN4 %pos vs same-patient T/NK fraction / cytotoxicity / IFN /
   exhaustion, and for Q4 vs Q1 (tails only; thin if n<8 or a tail <3).
2. **Paired n** — patient has ≥10 CLDN4-high, ≥10 CLDN4-low, and ≥20 T/NK.
   Used for within-patient ligand Δ (Wilcoxon signed-rank).

GSE131907 is **patient-pooled**. That differs from the sample-level T/NK extract
in PR #320; unique-patient n is the header n here.

## Ligand activity

Potential ligand: expressed in ≥10% of CLDN4-high malignant cells in ≥25% of
paired patients, ≥1 receptor expressed in the T/NK panel, and a column in the
v2 matrix.

- Gene-set activity: Pearson / AUROC / AUPR of prior target scores vs a priori
  IFN, cytotoxicity, or exhaustion membership on T/NK-expressed background genes.
- Empirical activity (primary for ranking): Pearson of prior columns vs the
  patient-level T/NK gene contrast (rank-biserial of CLDN4-high vs low **patients**,
  median split on malignant CLDN4 %pos). Same-patient T/NK only.

Background is panel ∩ prior ∩ T/NK-expressed genes — **not** the full transcriptome.

## Meta

Spearman ρ is Fisher-z pooled across cohorts with weights n−3. I² is
DerSimonian–Laird on the z scale. p-values are descriptive.

## Reuse

Add a cohort by writing an extractor in `scripts/02_extract_panels.py` that
emits `{COHORT}_panel.parquet` (malignant + T/NK rows, UMI + `ncount`) and
`{COHORT}_patient_n.tsv`, then append the name to `COHORTS` in `03_analyze.py`.
