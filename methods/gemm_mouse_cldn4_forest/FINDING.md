# FINDING — Mouse-level forest of epithelial Cldn4 %pos vs immune / IFN

Public autochthonous GEMM lung scRNA only. Private 8 KL matrices were not used. Cell-line transplants (LKR13, LLC, tail-vein KP lines) and bulk RNA were not put in the forest.

The primary specification was locked before the pooled p-value was read. Cutoffs were not searched to maximize a thesis-aligned (negative) association. The human concordant-4 result (malignant CLDN4+ % vs T/NK, negative) is the scientific question, not a target the mouse p-value had to hit.

## Primary result

Exposure is the percent of epithelial cells with Cldn4 detected. Immune endpoint is the T/NK fraction in unsorted tumor digests. IFN endpoint is an IFN-only score inside epithelial cells (MHC genes are not in that score). One row is one mouse. Studies enter the pool only with at least 4 mice. The pool is REML on Fisher z with a Hartung-Knapp interval. A pooled p-value is reported only when at least 3 studies are estimable.

### Cldn4 %pos vs T/NK fraction

Pooled REML-HK ρ=+0.477, 95% CI -0.604 to +0.940, p=0.302, I²=72%, k=5, mice=34.

| Study | Mice | Spearman ρ | p | Genotype-partial ρ |
|---|---:|---:|---:|---|
| GSE165641 | 2 | not estimable (n<4 or no variance) | | |
| GSE179501 | 4 | +0.200 | 0.800 | NA |
| GSE180963 | 2 | not estimable (n<4 or no variance) | | |
| GSE201247 | 4 | +0.800 | 0.200 | NA |
| GSE264739 | 6 | +0.943 | 0.005 | +0.750 (p=0.086) |
| GSE266323 | 4 | +0.400 | 0.600 | NA |
| GSE295824 | 16 | -0.471 | 0.066 | -0.030 (p=0.911) |

### Cldn4 %pos vs epithelial IFN

Pooled REML-HK ρ=+0.092, 95% CI -0.682 to +0.769, p=0.796, I²=65%, k=5, mice=56.

| Study | Mice | Spearman ρ | p | Genotype-partial ρ |
|---|---:|---:|---:|---|
| GSE149813 | 2 | not estimable (n<4 or no variance) | | |
| GSE154977 | 2 | not estimable (n<4 or no variance) | | |
| GSE154989 | 24 | -0.442 | 0.031 | -0.222 (p=0.297) |
| GSE165641 | 2 | not estimable (n<4 or no variance) | | |
| GSE179502 | 6 | -0.600 | 0.208 | +0.750 (p=0.086) |
| GSE180963 | 2 | not estimable (n<4 or no variance) | | |
| GSE264739 | 6 | +0.771 | 0.072 | +0.000 (p=1.000) |
| GSE295824 | 16 | +0.312 | 0.240 | +0.131 (p=0.629) |
| GSE319598 | 4 | +0.800 | 0.200 | NA |

## What the primary numbers say

The pooled T/NK association is positive (higher epithelial Cldn4 %pos, higher T/NK fraction) and the interval covers zero. The pooled IFN association is near zero. That is the opposite direction from the human concordant-4 patient-level result, and it is not significant.

The only within-study immune Spearman below 0.05 is GSE264739 (KP vs KPP, ρ positive). Its genotype-partial correlation stays positive and is no longer below 0.05. The only within-study IFN Spearman below 0.05 is GSE154989. Its genotype-partial correlation (K vs KP) is weaker and not below 0.05, which repeats the earlier public note that the unadjusted plate-seq inverse track is largely genotype composition.

GSE295824 is the largest unsorted study (16 Sox2-GEMM mice). Unadjusted Cldn4 %pos vs T/NK is negative and not below 0.05; the genotype-partial correlation is about zero. Dropping studies whose Cldn4 %pos range is under 5 percentage points leaves only two immune studies, so that spec has no pooled p-value by the pre-specified k≥3 rule.

## Sensitivity grid

Every spec below was named in the script before looking at which p was smallest. The smallest Hartung-Knapp p in the pre-listed grid is 0.283 for `include_WT_GSE201247` / immune (ρ=+0.476, k=5). That spec is not the primary. It was not promoted because it won a search.

See `tables/sensitivity_grid.tsv` and `figures/sensitivity_pooled.png`.

## What was scored

- GSE149813: 2 rows in the assembled table
- GSE154977: 4 rows in the assembled table
- GSE154989: 30 rows in the assembled table
- GSE165641: 2 rows in the assembled table
- GSE179501: 4 rows in the assembled table
- GSE179502: 6 rows in the assembled table
- GSE180963: 2 rows in the assembled table
- GSE188436: 4 rows in the assembled table
- GSE201247: 6 rows in the assembled table
- GSE264739: 6 rows in the assembled table
- GSE266323: 4 rows in the assembled table
- GSE281744: 2 rows in the assembled table
- GSE281964: 4 rows in the assembled table
- GSE295824: 16 rows in the assembled table
- GSE317576: 6 rows in the assembled table
- GSE319598: 6 rows in the assembled table
- GSE338088: 5 rows in the assembled table

Newly scored matrices use one marker rule: epithelial = Epcam>0 or Sftpc>0 or (Krt8>0 and Ptprc==0); T/NK = Cd3d, Cd3e, Nkg7, or Ncr1 > 0; IFN = mean log1p(CP10k) of Stat1, Stat2, Irf1, Irf7, Irf9, Isg15, Ifit1, Ifit2, Ifit3, Mx1, Oasl2, Rsad2, Ifih1, Ddx58, Ifnb1. QC is 200–8000 genes, at least 500 UMIs, mitochondrial fraction under 25%.

Already published mouse tables (GSE154989, GSE179502, GSE179501, GSE201247, GSE266323, GSE154977, GSE165641, GSE180963) are reused at the mouse unit. Their IFN columns are used only when they are IFN-only. The three-study composite IFN/MHC score stays in the sensitivity grid.

## Kept out of the mouse forest

- Pools: GSE281744 (2 mice pooled per library), GSE319598 stromal pools, GSE133604 treatment-arm pools.
- Libraries with no mouse id: GSE188436, GSE281964 (pre/post), GSE317576 (condition labels).
- Epithelial sorts cannot support a T/NK fraction: GSE154989, GSE179502, GSE154977, GSE149813, GSE319598 eGFP+ tumor cells. They stay in the IFN forest when a mouse id and an IFN-only score exist.
- GSE338088 snRNA was scored. The epithelial marker rule returned fewer than 20 cells in every mouse and no Cldn4-positive epithelial cells, so neither endpoint is estimable.
- Not opened: GSE297023 (9.9 Gb Parse), GSE322632 (1.9 Gb Seurat, CD45-enriched), GSE277777 per-mouse tumor h5ads (the 252 Mb `luadAdata.h5ad` is the GSE154989 plate object, not a second cohort; combined file is 7.5 Gb and several GSMs are tail-vein transplants or HTO hashes), GSE176185 (scaled matrix, percent-positive undefined), GSE203447 (rds, sort fractions, n<4), private 8 KL.
- Not autochthonous GEMM lung scRNA: LLC, LKR13, bulk GSE6135.

GSE295824 (Sox2 GEMM, 16 mice named in the archive) is in the primary forest for both endpoints.
