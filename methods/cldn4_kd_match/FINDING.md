# CLDN4 KD match: lung Chronos effects and public knockout omics

DepMap lung lines do not depend on CLDN4. Expression Public on those same lines is baseline RNA, not a transcriptome measured after CLDN4 knockout. No public CLDN4 knockout has paired RNA-seq and proteomics.

Cohort rule, same as PR #331: `OncotreeLineage == Lung` and `ModelType == Cell Line`. Chronos is scaled so the median common essential sits near −1 and the median nonessential sits near 0. A line is counted as dependent when the gene effect is below −0.5. CLDN4 (Entrez 1364) is absent from `CRISPRInferredCommonEssentials.csv`.

## Chronos in lung lines

260 lung cell-line models are in DepMap Public 24Q4 `Model.csv`. 126 have a CLDN4 score in `CRISPRGeneEffect.csv`.

| cohort | n | median | IQR | min | max | n < −0.5 | n < −1 |
|---|---:|---:|---|---:|---:|---:|---:|
| lung cell lines | 126 | +0.047 | −0.037 to +0.109 | −0.294 | +0.394 | 0 | 0 |
| NSCLC | 98 | +0.044 | −0.045 to +0.101 | −0.278 | +0.394 | 0 | 0 |
| LUAD | 53 | +0.055 | +0.013 to +0.111 | −0.278 | +0.250 | 0 | 0 |
| LUSC | 21 | −0.037 | −0.131 to +0.044 | −0.241 | +0.093 | 0 | 0 |
| other NSCLC | 24 | +0.078 | −0.044 to +0.148 | −0.255 | +0.394 | 0 | 0 |
| SCLC / NET | 26 | +0.059 | −0.009 to +0.120 | −0.294 | +0.381 | 0 | 0 |

84 of 126 scores are positive. The lowest score is NCI-H82 (−0.294). The highest is LCLC-97TM1 (+0.394). LUSC is the lowest subgroup median and is still above −0.5.

The 123 lines that also have Expression Public have median Chronos **+0.050**. That is the same overlap n and the same median already stated in PR #331. This page does not recompute the IFN / MHC-I / CD274 correlations from that PR.

Per-line values: `tables/lung_lines_cldn4.tsv`. Figure: `figures/cldn4_chronos_lung.png` (dashed line at −0.5).

## Expression Public

The portal name Expression Public is the protein-coding log2(TPM+1) matrix. For 24Q4 that file is `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (DepMap forum: the 24Q2 name maps to this filename; 24Q4 keeps the same filename). It is baseline RNA-seq of the models that were screened. It is not RNA collected after CLDN4 CRISPR.

The no-captcha catalog (`https://depmap.org/portal/api/no-captcha/download/files`, queried 2026-09-21) lists later releases, and their bulk URLs are empty:

| release | date | CRISPRGeneEffect URL | Expression file | Expression URL |
|---|---|---|---|---|
| DepMap Public 26Q1 | 2026-04-01 | no | `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | no |
| DepMap Public 25Q3 | 2025-09-25 | no | same new filename | no |
| DepMap Public 25Q2 | 2025-06-27 | no | old and new filenames | no |
| DepMap Public 24Q4 | 2024-12-16 | yes | `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | yes |

24Q4 is the newest DepMap Public release that serves both files. The catalog has no matrix whose name is expression after CRISPR. The only filename containing “GeCKO” is `gecko_ceres_gene_effects.csv`, a dependency score.

On the 24Q4 file, 214 lung cell lines have CLDN4 expression (median log2(TPM+1) = 6.46). 197 of 214 are at least 1. Among the 123 lines with both assays, Spearman correlation of Chronos versus Expression Public is **+0.159** (p = 0.079, bootstrap 95% CI −0.031 to +0.335, BH q = 0.119 across the three pre-specified cohorts). NSCLC n = 95, ρ = +0.186, CI −0.020 to +0.384. SCLC/NET n = 26, ρ = +0.153, p = 0.46. 115 of 123 lines are both expressed (log2(TPM+1) ≥ 1) and above −0.5.

NCI-H1688 (ACH-002170), the SCLC line in the 2025 CLDN4-knockout paper, has Expression Public 6.02 and no Chronos score in 24Q4 or in the 26Q1 parameter file.

## 26Q1 Chronos parameters

Portal `CRISPRGeneEffect.csv` for 26Q1 has an empty bulk URL. Figshare 31660582 (“Chronos parameters (Public 26Q1)”) does include `gene_effect.csv` and does not include an expression matrix. The first 12 inferred common essentials that are present in both matrices have mean-of-medians −1.22 (24Q4 `CRISPRGeneEffect`) and −1.04 (26Q1 `gene_effect`), so the parameter file is on the same scaled gene-effect axis.

Joined to 24Q4 lung labels, 26Q1 n = 126, median **+0.184**, minimum −0.218, and still **0/126 below −0.5**. Spearman agreement with 24Q4 on these 126 lines is ρ = +0.809 (p = 2.0×10⁻³⁰, CI +0.721 to +0.872). No line crosses −0.5 between releases. Sixteen 26Q1 model IDs are absent from the 24Q4 model table, so they are not assigned to lung here (`tables/chronos_26Q1_models_absent_from_24Q4_metadata.tsv`).

## Paired RNA-seq and proteomics of a CLDN4 knockout

Queries on 2026-09-21 (GEO, SRA, PubMed, Europe PMC, PRIDE) are in `tables/search_log.json`. PRIDE keywords `CLDN4`, `claudin-4`, and `Cldn4` each returned 0 projects. The 24Q4 Figshare file set (73 files) contains 0 proteomics matrices.

Public CLDN4 loss-of-function transcriptomes already scored in earlier PRs, none with a proteome:

| accession | system | RNA | proteomics |
|---|---|---|---|
| GSE207704 | MCF7 and T47D CRISPR knockout | yes | no |
| GSE50927 | mouse lung germline knockout | yes | no (BAL protein leak is not mass spec) |
| GSE22493 | SKOV-3 siRNA microarray | array | no |

PMID 41016339 reports RNA-seq after CLDN4 knockout in NCI-H1688. Europe PMC `hasData` is N, there is no PMCID, and GEO/SRA queries returned no matching accession. The abstract does not describe proteomics. The only GEO series returned by `H1688 AND (CLDN4 OR SAA1 OR knockout)` is GSE174462, which is FOXM1 knockdown.

PubMed hits for knockout plus proteome (PMID 38345099, 37889067, 19318328) and the other lung RNA hits (PMID 36840413 ELF3 knockdown, PMID 33006362 CRAD knockdown) are different perturbations. They are listed in `tables/ko_omics_inventory.tsv`.

Non-GEO RNA deposits were already searched in PR #592 (none open beyond the three GSE series). Sci Rep 2025 CLDN4 CRISPRi (PMID 41214101) already has no open raw RNA-seq or proteomics (PR #586). Those matrices were not re-downloaded.

## PR inventory

| PR | What it already covers |
|---|---|
| #331 | CLDN4 Chronos vs IFN/MHC-I/CD274 on the n=123 overlap. Median Chronos +0.050, 0/123 < −0.5. |
| #304 | CLDN4 RNA vs IFN/MHC-I/CD274, lung n=214. |
| #388, #575 | Baseline Gygi protein, not a knockout proteome. #575 also recorded blank 26Q1 expression URLs. |
| #120, #106, #59, #144 | Earlier DepMap lung CRISPR / IFN cuts, mostly TACSTD2. |
| #43, #98, #104, #131, #155, #287, #310, #512, #517, #580, #585 | GSE207704, GSE50927, GSE22493 transcriptomes. |
| #592 | No open CLDN4-knockdown RNA-seq outside those GEO series. H1688 RNA-seq has no accession. |
| #586 | Sci Rep 2025 CLDN4 CRISPRi: no open raw RNA-seq or proteomics. |
| #558 | Spatial proteomics hunt for CLDN4 protein + CD8, empty. Not a knockout. |

This wave adds the lung-line Chronos distribution (n=126, including lines without RNA), the Expression Public catalog query, the Chronos-versus-expression test, and the 26Q1 parameter-file check.

## 中文

DepMap 肺细胞系里 CLDN4 不是依赖性基因。24Q4 Chronos：126 系，中位数 +0.047，0/126 低于 −0.5。Expression Public 是同一批细胞的基线 RNA（log2(TPM+1)），不是 CLDN4 敲除之后的转录组；26Q1 目录里该文件和 CRISPRGeneEffect 的批量链接都是空的，能直接下到两者的最新发行版是 24Q4。两者重叠 n=123，Spearman +0.159，置信区间含 0。公开检索没有 CLDN4 敲除的配对 RNA-seq + 蛋白组。H1688（PMID 41016339）的 RNA-seq 没有入库号，该系有表达、没有 Chronos。
