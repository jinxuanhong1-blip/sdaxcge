# NicheNet-style ligand activity on GSE207422

**ADDITIVE scRNA.** Ligands from **TACSTD2-high / CLDN4-high malignant cells** scored against **T/NK** gene sets (a priori cytotoxicity and exhaustion; empirical NMPR vs MPR). Combinatorial **MPR vs NMPR** receivers. Public files only.

| | |
| --- | --- |
| Cohort | [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) (Hu et al., *Genome Med* 2023; PMID 36869384) |
| Prior | NicheNet-v2 ligand–target matrix + LR network ([Zenodo 7074291](https://doi.org/10.5281/zenodo.7074291)) |
| Labels | DRMref public cell types (Liu et al., *NAR* 2024), not Hu CopyKAT |
| Unit | **Patient** (n=4 MPR including pCR, n=8 NMPR). Cells are counts, not replicates. |
| Not used | GSE253013 (9.3 GB RDS; no MPR/NMPR) |

- Playbook: [`playbook.md`](playbook.md)
- What was found: [`FINDING.md`](FINDING.md)
- Methods: [`METHODS.md`](METHODS.md)
- Results: [`results/`](results/)
- Scripts: [`scripts/`](scripts/)
