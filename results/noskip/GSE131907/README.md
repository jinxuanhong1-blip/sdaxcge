# GSE131907 no-skip — epithelial TACSTD2 vs immune

**Dataset:** Kim et al., *Nat Commun* 2020 (PMID 32385277); GEO [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907).
**Cells used:** 208,506 / 58 samples / 44 patients. Nothing subsampled.

## Verdict

TACSTD2 is epithelial-restricted versus immune (cell-level median log2TPM 1.75 vs 0.00, %pos 74.5 vs 3.1, MWU p=<1e-300, n_epi=36467, n_imm=164525). Sample-level epithelial TACSTD2 vs immune fraction is a null (tumor-site ρ=0.06, p=0.719, n=36; tLung ρ=-0.28, p=0.401, n=11). No ICI / MPR labels. Size was not a skip reason: both processed text matrices were used.

## What was used (size was not a skip reason)

- Author cell annotation (`GSE131907_Lung_Cancer_cell_annotation.txt.gz`, 208,506 cells).
- Processed raw UMI matrix (`GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`, 409 MB gzip; 29634 genes).
- Author normalized log2(TPM+1) matrix (`GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz`, 3.07 GB gzip; 29634 genes). This is the file a prior A3 pass skipped for size.
- GEO series matrix for patient / stage / origin.
- Selected genes streamed from both full matrices (TACSTD2, CLDN4, lineage markers). Other gene rows were scanned for inventory, not loaded.
- RDS copies were **not** needed: they hold the same two matrices. This environment has no R; the text supplements were used instead.
- EGA FASTQ (`EGAD00001005054`) was **not** used: it is controlled-access raw data, not a processed GEO supplement.

## Definitions

- Epithelial = author `Cell_type == Epithelial cells` (36,467 cells). In tLung these are mostly tS1/tS2/tS3; author `Malignant cells` are labeled in metastases / PE / tL-B.
- Immune = T lymphocytes + NK cells + B lymphocytes + Myeloid cells + MAST cells (164,525 cells).
- T/NK is reported as a secondary split to match the earlier A3 write-up.
- Primary expression metric = author log2(TPM+1). Sensitivity = raw UMI / log1p(UMI).
- Eligible sample for sample-level tests: ≥20 epithelial and ≥20 immune (or T/NK) cells.
- Detection = value > 0 in the stated matrix.

## Honest limits

- Treatment-naive LUAD atlas. No ICI, RECIST, or MPR/NMPR labels in GEO.
- Author annotations are used as-is. Doublets / ambient RNA are not re-called.
- Sample-level tLung n is small (11 tumors). A ρ≈−0.45 claim is underpowered here.
- Cell-level p-values are tiny because n is ~200k; effect size and detection rate are the meaningful numbers.

## Cell-level results

- all_cells: epithelial vs immune TACSTD2 log2TPM: median 1.754 vs 0.000; %pos 74.5 vs 3.1; Cliff's δ=0.722; MWU p=<1e-300; n=36,467/164,525
- all_cells: epithelial vs T/NK TACSTD2 log2TPM: median 1.754 vs 0.000; %pos 74.5 vs 1.6; Cliff's δ=0.729; MWU p=<1e-300; n=36,467/91,227
- all_cells: epithelial vs myeloid TACSTD2 log2TPM: median 1.754 vs 0.000; %pos 74.5 vs 7.3; Cliff's δ=0.700; MWU p=<1e-300; n=36,467/42,245
- all_cells: epithelial vs immune TACSTD2 raw UMI: median 2.000 vs 0.000; %pos 74.5 vs 3.1; Cliff's δ=0.730; MWU p=<1e-300; n=36,467/164,525
- tumor_sites: epithelial vs immune TACSTD2 log2TPM: median 1.756 vs 0.000; %pos 75.0 vs 4.6; Cliff's δ=0.714; MWU p=<1e-300; n=32,764/91,053
- tumor_sites: epithelial vs T/NK TACSTD2 log2TPM: median 1.756 vs 0.000; %pos 75.0 vs 2.8; Cliff's δ=0.721; MWU p=<1e-300; n=32,764/48,012
- tumor_sites: epithelial vs myeloid TACSTD2 log2TPM: median 1.756 vs 0.000; %pos 75.0 vs 9.7; Cliff's δ=0.687; MWU p=<1e-300; n=32,764/24,285
- tumor_sites: epithelial vs immune TACSTD2 raw UMI: median 2.000 vs 0.000; %pos 75.0 vs 4.6; Cliff's δ=0.727; MWU p=<1e-300; n=32,764/91,053
- tLung: epithelial vs immune TACSTD2 log2TPM: median 2.503 vs 0.000; %pos 86.0 vs 4.9; Cliff's δ=0.833; MWU p=<1e-300; n=7,270/35,506
- tLung: epithelial vs T/NK TACSTD2 log2TPM: median 2.503 vs 0.000; %pos 86.0 vs 3.2; Cliff's δ=0.835; MWU p=<1e-300; n=7,270/19,591
- tLung: epithelial vs myeloid TACSTD2 log2TPM: median 2.503 vs 0.000; %pos 86.0 vs 10.1; Cliff's δ=0.818; MWU p=<1e-300; n=7,270/8,794
- tLung: epithelial vs immune TACSTD2 raw UMI: median 5.000 vs 0.000; %pos 86.0 vs 4.9; Cliff's δ=0.843; MWU p=<1e-300; n=7,270/35,506
- all_cells: epithelial vs immune CLDN4 log2TPM: median 2.490 vs 0.000; %pos 83.6 vs 2.7; Cliff's δ=0.822; MWU p=<1e-300; n=36,467/164,525

- Lineage sanity: epithelial EPCAM mean log2TPM=2.095 vs immune 0.043; epithelial PTPRC=0.031 vs immune 1.501.
- UMI vs log2TPM concordance: cell-level TACSTD2 Spearman ρ=0.996; sample epithelial means ρ=0.834.

## Sample-level results

- tumor_sites: epi TACSTD2 mean log2TPM vs immune fraction: ρ=0.062, p=0.719, n=36
- tumor_sites: epi TACSTD2 mean log2TPM vs T/NK fraction: ρ=0.102, p=0.555, n=36
- tumor_sites: epi TACSTD2 %pos vs immune fraction: ρ=-0.012, p=0.945, n=36
- tumor_sites eligible-T/NK: epi TACSTD2 mean log2TPM vs T/NK fraction: ρ=0.102, p=0.555, n=36
- tLung: epi TACSTD2 mean log2TPM vs immune fraction: ρ=-0.282, p=0.401, n=11
- tLung: epi TACSTD2 mean log2TPM vs T/NK fraction: ρ=0.055, p=0.873, n=11
- nLung: epi TACSTD2 mean log2TPM vs immune fraction: ρ=-0.445, p=0.17, n=11
- mets/PE/tL-B: malignant TACSTD2 mean log2TPM vs immune fraction: ρ=-0.021, p=0.929, n=21
- mets/PE/tL-B: malignant TACSTD2 mean log2TPM vs T/NK fraction: ρ=0.013, p=0.955, n=21
- tumor_sites: epi TACSTD2 mean log1p UMI vs immune fraction: ρ=0.142, p=0.407, n=36
- tumor_sites: epi TACSTD2 mean log1p UMI vs T/NK fraction: ρ=0.168, p=0.326, n=36
- paired tumor samples: epi TACSTD2 > immune TACSTD2 (Wilcoxon signed-rank, log2TPM): n=36, p=1.02e-10; median %pos epithelial=76.1 vs immune=5.1 (T/NK=2.9).

## Files

| File | Role |
|---|---|
| `summary.json` / `audit.json` | Verdict and methods audit |
| `sample_metadata.tsv` | GEO patient / stage / origin |
| `cell_type_composition.tsv` / `sample_composition.tsv` | Author labels, full atlas |
| `expr_by_celltype_origin.tsv` / `expr_by_subtype_origin.tsv` | TACSTD2/CLDN4 by label |
| `per_sample_tacstd2.tsv` | Sample-level epithelial and immune TACSTD2 |
| `cell_level_contrasts.tsv` | MWU epithelial vs immune / T/NK / myeloid |
| `sample_level_associations.tsv` | Spearman tests |
| `paired_compartment_test.json` / `sanity_checks.json` / `umi_vs_log2tpm_concordance.json` | Extra tests |
| `fig_tacstd2_by_celltype.png` | Mean and %pos by author cell type |
| `fig_tacstd2_epithelial_vs_immune.png` | Primary compartment boxplot |
| `fig_epi_tacstd2_vs_immune_fraction.png` | Sample-level correlation panels |
| `fig_tacstd2_epithelial_subtypes.png` | tS / malignant / normal epi subtypes |

## Reproduce

```bash
bash analysis/01_download.sh
python3 analysis/02_extract.py
python3 analysis/03_analyze.py
```

