# Methods — scRNA CytoTRACE-like stemness / cycling (public only)

Additive module. Lives under `methods/scrna_cytotrace/`. No unpublished objects.

## What was (and was not) run

**CytoTRACE2 was not run.** The CytoTRACE2 R package (Kang et al. *Nat Commun* 2024) and its pretrained model are not installed here, and this box has no R/Bioconductor stack. The prompt allowed “CytoTRACE2 **or** a simple stemness/cycling score.” This module uses the latter, documented as such.

Stemness proxy (Gulati et al. *Science* 2020 CytoTRACE idea, not the R package):

1. Per cell, count genes with UMI > 0 (`n_genes`) and total UMI.
   GSE207422: counted while streaming the public UMI TSV.
   GSE241934: author Seurat `nFeature_RNA` / `nCount_RNA` (same quantities).
2. Among **malignant / epithelial** cells of that dataset (GSE241934: within IIT and REAL separately), OLS: `n_genes ~ log1p(total_UMI)`.
3. Residual = observed − fitted. Rank-scale residuals to `[0, 1]`. **Higher = more genes than expected for depth = more stem-like.**

This is a **CytoTRACE-like** gene-count potency score. It is not `cytotrace2()` output.

Cycling: Tirosh et al. *Science* 2016 S-phase and G2/M lists (Seurat `cc.genes`). Score = mean `log1p(CP10k)` of genes present in the matrix. `cycle_score = (S + G2M) / 2`.

- **Cycling** = top quartile of `cycle_score` among malignant cells (per dataset / cohort).
- **Non-cycling** = bottom quartile.

Keratin / differentiation: mean `log1p(CP10k)` of simple keratins (`KRT7/8/18/19`), basal/squamous keratins (`KRT5/6A/6B/14/17`), and alveolar/club markers (`SFTPA1/2`, `SFTPB/C/D`, `NAPSA`, `AGER`, `SCGB1A1`, `SCGB3A2`).

## Data (public GEO supplementary only)

| Accession | What | Malignant definition | Response |
|---|---|---|---|
| **GSE207422** (Hu et al. *Genome Med* 2023, PMID 36869384) | 24,292 genes × 92,330 barcodes UMI TSV | Marker epithelial (argmax lineage on log1p CP10k) **and** normal-lung score ≤ epithelial 75th percentile. Author CopyKAT IDs are **not** on GEO. | GEO `Pathologic Response`. pCR counted as MPR. **Primary = 12 post-treatment resections.** Pre-biopsy not pooled into the primary test. |
| **GSE241934** (Zhang/Zhong et al. *Cell Rep Med* 2024, PMID 38897205; NEOTIDE + real-world) | Public MTX + author `major.cell.type` | `major.cell.type == Epi` (author). IIT = 11 EGFR-mutant trial tumors; REAL = 34 WT LUAD/ASC neoadjuvant IO+chemo. | Author `Pathological Response`. `non-MPR` → NMPR; pCR → MPR. |

Raw FASTQ / GSA-Human was not used. No fabricated accessions.

QC: drop cells with total UMI < 200. Samples with fewer than **20** malignant/epithelial cells are excluded from the primary sample-level tests; a **min-10** sensitivity is stored.

## Tests (pre-specified)

**Unit = sample / patient.** Cell-level Spearman and MWU are computed and labeled **exploratory** (pseudoreplication).

1. Spearman ρ of sample-mean malignant `TACSTD2` / `CLDN4` (`log1p` CP10k) vs CytoTRACE-like, cycling score, keratin, alveolar/diff, and fraction cycling.
2. **Split by MPR:** two-sided Mann–Whitney NMPR vs MPR on those sample means. Also Spearman **within** MPR and within NMPR.
3. **Combinatorial cycling vs non-cycling:** paired Wilcoxon of sample-mean TACSTD2 (and CLDN4) in cycling vs non-cycling malignant cells. Same for CytoTRACE-like.
4. **Extra — is TACSTD2-high more differentiated / keratin?** Paired Wilcoxon of keratin (and CytoTRACE-like, alveolar) in within-sample TACSTD2 top vs bottom quartile.

Two-sided p-values. Honest n. No multiple-testing theater; every test is in `results/stats.tsv`.

## Re-run

```bash
python3 methods/scrna_cytotrace/scripts/download.py
python3 methods/scrna_cytotrace/scripts/extract_gse207422.py
python3 methods/scrna_cytotrace/scripts/extract_gse241934.py
python3 methods/scrna_cytotrace/scripts/analyze.py
```

Matrices stay in `/tmp/scrna_cytotrace/` and are not committed.
