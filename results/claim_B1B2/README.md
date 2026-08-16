# Claims B1+B2: public-data audit

## Bottom line

| Claim | Verdict | Public-data result |
|---|---|---|
| **B1. CLDN4 is the top surface gene coexpressed with TACSTD2 in TCGA pan-cancer.** | **No.** | Under a prespecified SURFY surfaceome, **PVRL4 (NECTIN4)** is #1. CLDN4 is #14 by pooled Spearman correlation (ρ=0.615, n=11,060), #29 in primary tumors (ρ=0.571, n=9,701), and #5 after a cancer-type fixed-effect sensitivity analysis (ρ=0.398). |
| **B2. CCLE NSCLC CLDN4–TROP2 protein correlation is ρ=0.69 at n=118.** | **No as stated.** | In the CCLE TMT mass-spectrometry dataset, there are only 64 NSCLC models and 36 have both proteins measured: ρ=0.731 (p=4.24×10^-7). The claimed n=118 cannot be a complete-pair result from this 375-model protein study. DepMap RNA gives ρ=0.660 (n=137; 20Q2) and ρ=0.662 (n=141; 24Q2), not a protein measurement. |

## B1 methods

- Expression: UCSC Xena `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`, batch-effects-normalized TCGA PanCan mRNA, `log2(norm_value+1)`.
- Surface genes: the representative human surfaceome from SURFY (5% false-positive-rate threshold plus GPI-anchored proteins), 2,886 proteins / 2,799 unique gene symbols. Of these, 2,619 symbols occur in the Xena matrix.
- Ranking: all eligible surface genes except the query gene TACSTD2 were ranked without a significance or expression prefilter. The primary result is two-sided Spearman correlation over every Xena sample. Sensitivities use primary tumors only and global-rank residualization by TCGA cancer type.
- “Surface gene” is not intrinsic to an RNA matrix. The result is therefore conditional on the stated, independently published surfaceome rather than a hand-selected target list.
- Xena reports some legacy symbols (for example PVRL4, now NECTIN4); symbols were not silently remapped because that could alter the eligible set.

The complete, unsuppressed ranking is in `B1_surface_gene_rankings.csv`. P-values shown as zero are floating-point underflow at this sample size, not literal probabilities of zero.

## B2 methods and modality check

The wording “CCLE protein” was tested against the public measured-protein resources:

1. **Nusinow et al. CCLE proteomics (TMT-MS):** 375 unique models in the source study. NSCLC was defined by the matching DepMap 20Q2 `lineage_subtype`. There are 64 NSCLC models in the protein matrix; CLDN4 is missing in 28, leaving 36 complete pairs. Both proteins are present, and the complete-pair Spearman result is ρ=0.731.
2. **Sanger ProCan-DepMap (DIA-MS):** TROP2/TACSTD2 has 948 records, but CLDN4 has zero records, so the pair cannot be calculated.
3. **CCLE RPPA500:** neither TACSTD2 nor CLDN4 is among the 447 antibodies.

For diagnosis of the likely modality mix-up, the same genes were tested in DepMap RNA-seq:

- 20Q2: NSCLC `lineage_subtype`, n=137, Spearman ρ=0.660.
- 24Q2: OncoTree primary disease `Non-Small Cell Lung Cancer`, n=141, Spearman ρ=0.662.

Thus, a strong association is reproducible at both RNA and protein levels, but **ρ=0.69, n=118, protein** is not. The release, model list, modality, and missing-value rule would need to be supplied to rescue that exact numerical claim.

## Reproduce

```bash
python3 -m pip install -r results/claim_B1B2/requirements.txt
python3 results/claim_B1B2/analyze.py
```

Raw public files are cached under the ignored `data/` directory. The run makes API calls to UCSC Xena, Cell Model Passports, and TCPA and downloads archived DepMap 20Q2/24Q2 expression matrices plus the Gygi CCLE protein matrix.

## Outputs

- `summary.json`: machine-readable verdict statistics and dataset availability
- `B1_surface_gene_rankings.csv`: complete honest surfaceome ranking
- `B1_top25_surface_genes.csv`: compact view of the leading genes
- `B2_CCLE_NSCLC_protein_pairs.csv`: 36 complete measured-protein pairs
- `B2_DepMap_NSCLC_RNA_pairs.csv`: release-specific RNA sensitivity data
- `B2_DepMap_NSCLC_RNA_scatter.png`: RNA scatter plots, clearly labeled as RNA
- `analyze.py`, `requirements.txt`: executable provenance

## Sources

- UCSC Xena PanCan Atlas hub: <https://pancanatlas.xenahubs.net>
- Bausch-Fluck et al., *The in silico human surfaceome* (PNAS 2018): <https://doi.org/10.1073/pnas.1808790115>
- DepMap 20Q2: <https://doi.org/10.6084/m9.figshare.12280541.v4>
- DepMap 24Q2: <https://doi.org/10.25452/figshare.plus.25880521.v1>
- Nusinow et al., *Quantitative Proteomics of the Cancer Cell Line Encyclopedia* (Cell 2020): <https://doi.org/10.1016/j.cell.2019.12.023>
- Gonçalves et al., *Pan-cancer proteomic map of 949 human cell lines* (Cancer Cell 2022): <https://doi.org/10.1016/j.ccell.2022.06.010>
- RPPA500: <https://doi.org/10.1038/s43018-024-00817-x>
