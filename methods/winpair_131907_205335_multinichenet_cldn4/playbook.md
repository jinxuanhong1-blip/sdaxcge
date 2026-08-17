# Playbook — MultiNicheNet on winning pair GSE131907+GSE205335 (CLDN4 only)

> Public processed files only. No EGA FASTQ. No invented statistics.
> English first; 中文 at the end of `FINDING.md`.

## Why this folder

The winning pair is **GSE131907 + GSE205335**. This is **additive MultiNicheNet / multi-sample NicheNet**, not a re-run of GSE207422 NicheNet (#334 NS) and not a merge that adds 207422.

Claim to rank (not to assume):

> Ligands from **CLDN4-high vs CLDN4-low malignant** cells, received by **same-patient T/NK**, can be ranked with a public NicheNet-v2 prior when **patient is the unit**.

## Stop rules

- Do not download or analyze GSE207422 here.
- Do not AND TACSTD2 into the sender call.
- Cells are not replicates.
- The prior is unsigned. Direction of IFN / cytotoxicity is patient-level only.
- GSE131907 has no MPR. Do not swap RECIST for MPR on GSE205335.

## Estimands

| ID | Estimand | Unit |
| --- | --- | --- |
| E1 | T/NK IFN and cytotoxicity vs malignant CLDN4 | patient (combined + per dataset) |
| E2 | Same vs T/NK fraction | patient |
| E3 | Ligand activity (Pearson / AUROC of v2 prior vs gene-set membership) | ligand rank |
| E4 | Paired CLDN4-high vs low ligand DE + MultiNicheNet-like prioritization | ligand × patient |

## Run

See `README.md`. Outputs must include `results/ligand_activity_table.tsv` and `FINDING.md`.
