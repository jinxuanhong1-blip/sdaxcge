# FINDING — LINCS/L1000 connectivity of CLDN4-loss profiles to NHEJ loss and STING/IFN

Additive public evidence. Numbers are written from `tables/class_summary.tsv`.

## Result

On T47D and MCF7, CLDN4-knockout connectivity to LINCS NHEJ-loss signatures and to IFN ligand signatures is centered at zero. Primary median cell-line Spearman values sit between -0.003 and 0.003. BH q versus 0, inside that primary family, is at least 0.647. The same rows’ cell-line median weighted connectivity scores are 0. Matched-background tests are in the table below.

Shared genes per scored signature are the intersection of the query rank with the characteristic-direction table (8018 genes on the GSE207704 queries). Lung lines in this signature set are A549, HCC515, and LCLC103H. That subset has two lines at most, so the lung Wilcoxon is not computed.

## What was queried

SigCom LINCS (characteristic-direction Level 5, metadata retrieved with this run) has **0** signatures with `pert_name = CLDN4`. iLINCS `SignatureMeta` treatment CLDN4 returned **0**. GSE92742 Phase I `pert_info` has **0** CLDN4 perturbagen. CLDN4 in the Phase I gene table: pr_gene_id=1364 pr_is_lm=0 pr_is_bing=1. TACSTD2, for context in the same catalog, has **21** signatures.

The connectivity query is therefore each public CLDN4-loss transcriptome, not a LINCS CLDN4 hairpin. Positive Spearman means the CLDN4-loss profile points the same way as the reference signature.

| Query | Contrast | Species | n | CLDN4 / Cldn4 log2FC | Genes in rank |
|---|---|---|---|---:|---:|
| GSE207704_T47D | T47D CLDN4-/- vs WT | human | 2 vs 2 | -1.039 | 12635 |
| GSE207704_MCF7 | MCF7 CLDN4-/- vs WT | human | 2 vs 2 | -0.750 | 12635 |
| GSE50927_naive_lung | Cldn4 KO vs WT, naive whole lung | mouse | 1 vs 1 | -6.061 | 21588 |
| GSE22493_SKOV3 | SKOV-3 CLDN4 siRNA vs CLDN4-overexpressing control | human | 3 vs 3 | -1.225 | 18682 |

GSE207704 deposits one FPKM per group (2 vs 2 already collapsed). GSE22493 control arm is CLDN4 overexpression. GSE50927 is naive whole-lung Cldn4 knockout, n=1 vs 1; mouse symbols were uppercased to overlap human LINCS symbols.

## Reference signatures

Core NHEJ loss: PRKDC, LIG4, XRCC4, XRCC5, XRCC6, NHEJ1, DCLRE1C. CRISPR knockout and shRNA are separate classes because Ku80 (XRCC5) and XRCC4 are shRNA in this build.

STING/IFN: ligand signatures for IFNA, IFNA1, IFNA2, IFNB1, and IFNG. STING agonist names diABZI, ADU-S100, cGAMP, 2'3'-cGAMP, vadimezan, and DMXAA have **0 / 0 / 0 / 0 / 0 / 0** signatures. TMEM173 and MB21D1 CRISPR knockouts are reported as STING/cGAS loss, a different class from IFN ligand treatment.

Background, seed 42: up to 8 other signatures of the same `pert_type` and cell line (CRISPR, shRNA, or ligand). The paired test uses cell-line median Spearman minus the matched background median.

## Cell-line summary

Unit is the cell line. Wilcoxon is two-sided on cell-line medians and needs at least 5 lines. BH q is inside the primary family only: T47D and MCF7, classes NHEJ CRISPR, NHEJ shRNA, and IFN ligand, scope all cell lines. Other rows are sensitivity and keep raw p.

| Query | Class | Cell lines | Signatures | Median Spearman | IQR | P vs 0 | q vs 0 | Median minus background | P vs background |
|---|---|---:|---:|---:|---|---:|---:|---:|---:|
| T47D CLDN4 KO | NHEJ CRISPR KO | 19 | 112 | 0.001 | -0.005–0.012 | 0.332 | 0.647 | 0.003 | 0.679 |
| T47D CLDN4 KO | NHEJ shRNA | 10 | 245 | -0.001 | -0.002–0.001 | 0.695 | 0.695 | 0.001 | 0.770 |
| T47D CLDN4 KO | IFN ligand | 14 | 201 | -0.003 | -0.009–0.003 | 0.296 | 0.647 | 0.007 | 0.153 |
| T47D CLDN4 KO | STING/cGAS KO | 10 | 40 | 0.004 | -0.014–0.020 | 0.770 | — | -0.014 | 0.770 |
| MCF7 CLDN4 KO | NHEJ CRISPR KO | 19 | 112 | 0.003 | -0.010–0.018 | 0.541 | 0.649 | 0.008 | 0.145 |
| MCF7 CLDN4 KO | NHEJ shRNA | 10 | 245 | -0.001 | -0.003–0.001 | 0.432 | 0.647 | -0.002 | 0.695 |
| MCF7 CLDN4 KO | IFN ligand | 14 | 201 | 0.002 | -0.001–0.009 | 0.194 | 0.647 | 0.006 | 0.135 |
| MCF7 CLDN4 KO | STING/cGAS KO | 10 | 40 | -0.005 | -0.014–0.001 | 0.232 | — | -0.009 | 0.695 |
| SKOV-3 siCLDN4 | NHEJ CRISPR KO | 19 | 112 | 0.001 | -0.003–0.012 | 0.241 | — | 0.002 | 0.293 |
| SKOV-3 siCLDN4 | NHEJ shRNA | 10 | 245 | -0.001 | -0.002–0.001 | 0.625 | — | -0.003 | 0.064 |
| SKOV-3 siCLDN4 | IFN ligand | 14 | 201 | 0.002 | 0.000–0.005 | 0.011 | — | 0.004 | 0.058 |
| SKOV-3 siCLDN4 | STING/cGAS KO | 10 | 40 | -0.000 | -0.010–0.004 | 0.625 | — | -0.001 | 0.695 |
| Mouse lung Cldn4 KO | NHEJ CRISPR KO | 19 | 112 | 0.001 | -0.008–0.010 | 0.922 | — | 0.005 | 0.922 |
| Mouse lung Cldn4 KO | NHEJ shRNA | 10 | 245 | 0.000 | -0.001–0.001 | 0.922 | — | -0.001 | 0.322 |
| Mouse lung Cldn4 KO | IFN ligand | 14 | 201 | 0.008 | -0.001–0.012 | 0.035 | — | 0.012 | 0.119 |
| Mouse lung Cldn4 KO | STING/cGAS KO | 10 | 40 | -0.013 | -0.028–0.010 | 0.232 | — | -0.016 | 0.432 |

### Lung lines only

Phase I `primary_site = lung`, hyphen-stripped ids, and LCLC103H (large-cell lung carcinoma in the 2021 CRISPR panel; absent from the 2017 cell table).

| Query | Class | Cell lines | Signatures | Median Spearman | IQR | P vs 0 | q vs 0 | Median minus background | P vs background |
|---|---|---:|---:|---:|---|---:|---:|---:|---:|
| T47D CLDN4 KO | NHEJ CRISPR KO | 2 | 12 | 0.003 | -0.005–0.010 | NA | — | 0.009 | NA |
| T47D CLDN4 KO | NHEJ shRNA | 2 | 54 | -0.002 | -0.003–-0.000 | NA | — | -0.007 | NA |
| T47D CLDN4 KO | IFN ligand | 2 | 15 | -0.003 | -0.006–-0.001 | NA | — | 0.009 | NA |
| T47D CLDN4 KO | STING/cGAS KO | 1 | 4 | 0.001 | 0.001–0.001 | NA | — | 0.005 | NA |
| MCF7 CLDN4 KO | NHEJ CRISPR KO | 2 | 12 | -0.005 | -0.009–-0.001 | NA | — | -0.012 | NA |
| MCF7 CLDN4 KO | NHEJ shRNA | 2 | 54 | 0.003 | 0.002–0.003 | NA | — | 0.000 | NA |
| MCF7 CLDN4 KO | IFN ligand | 2 | 15 | -0.005 | -0.011–-0.000 | NA | — | -0.003 | NA |
| MCF7 CLDN4 KO | STING/cGAS KO | 1 | 4 | 0.012 | 0.012–0.012 | NA | — | 0.015 | NA |

Secondary scores (median of cell-line medians) are `median_cell_cosine` and `median_cell_wtcs` in `tables/class_summary.tsv`. WTCS uses the top and bottom 100 CLDN4-loss genes, weighted by |log2FC|, and is 0 when the up and down enrichments are not opposite.

## NHEJ gene medians

### T47D

| Gene | Class | n sig | n cell lines | median Spearman | median WTCS |
|---|---|---:|---:|---:|---:|
| DCLRE1C | NHEJ CRISPR KO | 20 | 10 | 0.003 | 0.000 |
| LIG4 | NHEJ CRISPR KO | 20 | 10 | 0.017 | 0.000 |
| NHEJ1 | NHEJ CRISPR KO | 20 | 10 | -0.001 | 0.000 |
| PRKDC | NHEJ CRISPR KO | 20 | 10 | -0.012 | 0.000 |
| PRKDC | NHEJ shRNA | 45 | 10 | -0.004 | 0.000 |
| XRCC4 | NHEJ shRNA | 56 | 8 | -0.003 | 0.000 |
| XRCC5 | NHEJ shRNA | 88 | 9 | 0.002 | 0.000 |
| XRCC6 | NHEJ CRISPR KO | 32 | 19 | 0.002 | 0.000 |
| XRCC6 | NHEJ shRNA | 56 | 8 | -9.45e-04 | 0.000 |

### MCF7

| Gene | Class | n sig | n cell lines | median Spearman | median WTCS |
|---|---|---:|---:|---:|---:|
| DCLRE1C | NHEJ CRISPR KO | 20 | 10 | -0.006 | 0.000 |
| LIG4 | NHEJ CRISPR KO | 20 | 10 | 0.017 | 0.000 |
| NHEJ1 | NHEJ CRISPR KO | 20 | 10 | -0.016 | 0.000 |
| PRKDC | NHEJ CRISPR KO | 20 | 10 | -0.013 | 0.000 |
| PRKDC | NHEJ shRNA | 45 | 10 | -3.90e-04 | 0.000 |
| XRCC4 | NHEJ shRNA | 56 | 8 | 0.001 | 0.000 |
| XRCC5 | NHEJ shRNA | 88 | 9 | -0.002 | 0.000 |
| XRCC6 | NHEJ CRISPR KO | 32 | 19 | 0.004 | 0.000 |
| XRCC6 | NHEJ shRNA | 56 | 8 | 0.002 | 0.000 |

## Reading rule

The headline is the cell-line median Spearman and its comparison with the matched background. A small Spearman can produce a low p when many cell lines share a sign. Signature-level p-values are not used: each signature has thousands of genes, so a tiny correlation is automatically “significant” at gene resolution.

Signature files dropped at parse: 0.

## Reproduce

```bash
python3 -m pip install -r methods/lincs_cldn4_connectivity/requirements.txt
python3 scripts/lincs_cldn4_connectivity/analyze.py
```

Raw GEO and LINCS signature files are cached under `methods/lincs_cldn4_connectivity/data/` and are gitignored. Tables and figures are the committed result.
