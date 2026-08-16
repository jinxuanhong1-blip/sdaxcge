# Claims B1+B2: public-data audit

## Bottom line

| Claim | Verdict | Public-data result |
|---|---|---|
| **B1. CLDN4 is the top surface gene coexpressed with TACSTD2 in TCGA pan-cancer.** | **No.** | Under a prespecified SURFY surfaceome, **PVRL4 (NECTIN4)** is #1. CLDN4 is #14 by pooled Spearman correlation (ρ=0.615, n=11,060), #29 in primary tumors (ρ=0.571, n=9,701), and #5 after a cancer-type fixed-effect sensitivity analysis (ρ=0.398). |
| **B2. Gygi CCLE mass spec: CLDN4–TROP2 ρ=0.69 at n=45.** | **Yes.** | Nusinow/Gygi TMT-MS, Gygi `Tissue of Origin == Lung`, complete pairs: **n=45, Spearman ρ=0.693, Pearson r=0.687**. The earlier n=118 wording does not describe this protein matrix. n=45 is all Gygi lung, not NSCLC-only: 36 NSCLC + 9 SCLC. NSCLC-only complete pairs are n=36, ρ=0.731. |

## B1 methods

- Expression: UCSC Xena `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`, batch-effects-normalized TCGA PanCan mRNA, `log2(norm_value+1)`.
- Surface genes: the representative human surfaceome from SURFY (5% false-positive-rate threshold plus GPI-anchored proteins), 2,886 proteins / 2,799 unique gene symbols. Of these, 2,619 symbols occur in the Xena matrix.
- Ranking: all eligible surface genes except the query gene TACSTD2 were ranked without a significance or expression prefilter. The primary result is two-sided Spearman correlation over every Xena sample. Sensitivities use primary tumors only and global-rank residualization by TCGA cancer type.
- “Surface gene” is not intrinsic to an RNA matrix. The result is therefore conditional on the stated, independently published surfaceome rather than a hand-selected target list.
- Xena reports some legacy symbols (for example PVRL4, now NECTIN4); symbols were not silently remapped because that could alter the eligible set.

The complete, unsuppressed ranking is in `B1_surface_gene_rankings.csv`. P-values shown as zero are floating-point underflow at this sample size, not literal probabilities of zero.

## B2 methods

The corrected claim is **Gygi/Nusinow CCLE TMT mass spectrometry, n=45, ρ=0.69**, not RPPA and not n=118.

1. **Primary pairing (recovers the claim):** Gygi `Table_S1` `Tissue of Origin == Lung`. TACSTD2 is measured in 77 unique lung models; CLDN4 is measured in 45. Complete-pair Spearman ρ=0.693 (p=1.31×10^-7). Pearson r=0.687. The 45 models are 36 NSCLC and 9 SCLC by DepMap 20Q2 `lineage_subtype`.
2. **NSCLC-only sensitivity:** the same matrix restricted to DepMap 20Q2 `lineage_subtype == NSCLC` has 64 lung-NSCLC models but only 36 complete pairs (CLDN4 missingness). Spearman ρ=0.731.
3. **Other public protein matrices cannot produce this pair at n=45:** Sanger ProCan-DepMap DIA-MS has TROP2 (948 records) but zero CLDN4 records. CCLE RPPA500 has neither antibody.

DepMap RNA is reported only as a modality check, not as the B2 claim:

- 20Q2 NSCLC `lineage_subtype`, n=137, Spearman ρ=0.660.
- 24Q2 OncoTree primary disease `Non-Small Cell Lung Cancer`, n=141, Spearman ρ=0.662.

So the association is real at protein and RNA. The exact quoted number **ρ=0.69, n=45** is the Gygi lung complete-pair result. Calling that set “NSCLC” is slightly loose because it includes 9 SCLC lines. n=118 is not a complete-pair count in this protein study.

## Reproduce

```bash
python3 -m pip install -r results/claim_B1B2/requirements.txt
python3 results/claim_B1B2/analyze.py
```

Raw public files are cached under the ignored `data/` directory. The run queries UCSC Xena, Cell Model Passports, and TCPA, and downloads archived DepMap 20Q2/24Q2 expression matrices plus the Gygi CCLE protein matrix.

## Outputs

- `summary.json`: machine-readable verdict statistics and dataset availability
- `B1_surface_gene_rankings.csv`: complete honest surfaceome ranking
- `B1_top25_surface_genes.csv`: compact view of the leading genes
- `B2_Gygi_lung_protein_pairs.csv`: 45 complete Gygi lung protein pairs with lineage labels
- `B2_Gygi_protein_sensitivity.csv`: lung vs NSCLC-only protein correlations
- `B2_Gygi_lung_protein_scatter.png`: primary B2 scatter, NSCLC vs SCLC marked
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
