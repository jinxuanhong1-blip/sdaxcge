# NicheNet-style ligand activity — CLDN4-only senders on GSE207422

**ADDITIVE scRNA.** Ligands from **CLDN4-high malignant cells** (CLDN4 only; TACSTD2 is not used to call senders) scored against **T/NK IFN** and **cytotoxicity** gene sets. Combinatorial **MPR vs NMPR** receivers. Public files only. Python NicheNet-v2 prior (no R).

| | |
| --- | --- |
| Cohort | [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) (Hu et al., *Genome Med* 2023; PMID 36869384) |
| Prior | NicheNet-v2 ligand–target matrix + LR network ([Zenodo 7074291](https://doi.org/10.5281/zenodo.7074291)) |
| Labels | DRMref public cell types (Liu et al., *NAR* 2024), not Hu CopyKAT |
| Sender | Malignant **CLDN4** `log1p(CP10k)` ≥ malignant median |
| Unit | **Patient** (n=4 MPR including pCR, n=8 NMPR). Cells are counts, not replicates. |
| Not used | GSE253013 (9.3 GB RDS; no MPR/NMPR); TACSTD2 in the sender call |

This folder does **not** replace `methods/scrna_nichenet/` (TACSTD2∩CLDN4 dual-high).

- Playbook: [`playbook.md`](playbook.md)
- What was found: [`FINDING.md`](FINDING.md)
- Methods: [`METHODS.md`](METHODS.md)
- Results: [`results/`](results/)
- Scripts: [`scripts/`](scripts/)
