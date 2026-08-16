# Hunt GSE131907 — epithelial TACSTD2/CLDN4 vs CD8/NK/TLS; TJ/keratin GSEA

Kim et al., *Nat Commun* 2020 (PMID 32385277); GEO [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907).

208,506 cells / 58 samples. Author `normalized_log2TPM` streamed. TACSTD2 and CLDN4 used for association and GSEA grouping are **epithelial-only**.

## Verdict (honest)

| Claim | Result |
|---|---|
| TACSTD2 / CLDN4 are epithelial-restricted vs immune | **Supported.** Tumor-site TACSTD2 %pos 75.0 epithelial vs 2.8 T / 3.2 NK / 2.6 B / 9.7 myeloid. CLDN4 %pos 84.9 vs 2.9 / 3.0 / 3.0 / 8.0. Mann–Whitney p ≈ 0 because n is tens of thousands of cells; the evidence is the %pos gap, not the p-value. |
| Epithelial TACSTD2 anti-correlates with CD8 / NK / TLS | **Not supported.** Tumor-site n=36: CD8 ρ=0.11 p=0.54; NK ρ=−0.01 p=0.94; B ρ=−0.03 p=0.87; GC-B ρ=0.18 p=0.29; 12-chemokine ρ=0.16 p=0.35. All BH-FDR > 0.75. |
| Epithelial CLDN4 anti-correlates with CD8 / NK / TLS | **Nominal only; not supported after BH-FDR.** Tumor-site CLDN4 vs NK ρ=−0.40 p=0.015; vs B ρ=−0.37 p=0.025; vs T/NK ρ=−0.38 p=0.022; vs CD8 ρ=−0.32 p=0.055; vs 12-chemokine ρ=−0.27 p=0.11. Minimum BH-FDR among these = 0.15. Do not treat as a confirmed cold-neighborhood result. |
| TJ program up in TROP2-high epithelium | **Supported (sample-level prerank).** CUSTOM_TJ_CORE NES=2.21 FDR=0; still NES=2.17 FDR=0 after removing CLDN4 and TACSTD2 from the set. KEGG TJ NES=1.74 FDR=0.0075; GOBP TJ organization NES=1.79 FDR=0.0080. |
| Keratin / barrier up in TROP2-high epithelium | **Partially supported.** GOBP_KERATINIZATION NES=1.48 FDR=0.031; CUSTOM_BARRIER_KERATIN NES=1.50 FDR=0.027 (unchanged without CLDN4/TACSTD2). GOBP_CORNIFICATION is **not** up (NES=−1.06 FDR=0.33). Cell-level keratin is weaker (FDR 0.11–0.16). |
| Spatial neighborhood or histological TLS | **Cannot be tested.** Dissociated 10x. No coordinates. TLS here is B / GC-B / 12-chemokine proxy only. |
| ICI / MPR | **Cannot be tested.** Treatment-naive atlas. No response labels. |

Primary tLung n=11 is underpowered. One tLung TACSTD2–NK ρ=−0.69 p=0.019 does **not** survive BH-FDR (FDR=0.15) and is one of 24 tests.

## What was measured

- Epithelial = author `Cell_type == Epithelial cells`. TROP2 quartiles use tumor epithelium (`Malignant cells` + `tS1/tS2/tS3`) in tumor-site samples: Q4 ≥ 2.67 log2(TPM+1) n=7,784; Q1 ≤ 0.40 n=7,784.
- CD8 = author subtypes Cytotoxic / Exhausted / Naive / CD8-low CD8+ T. NK = `Cell_type == NK cells`.
- TLS proxies: B-cell fraction, GC-B (DZ+LZ) fraction, and the 12-chemokine mean (CCL2/3/4/5/8/18/19/21, CXCL9/10/11/13) averaged over **all cells** in the sample.
- Eligible tumor-site sample: ≥20 epithelial cells (36 samples: 11 tLung + mets / PE / tL-B).
- Primary GSEA: prerank by Spearman ρ of each gene's epithelial sample-mean vs epithelial TACSTD2, 1,000 permutations. This is the honest GSEA (n=36 samples).
- Secondary GSEA: Welch t of TROP2 Q4 vs Q1 tumor epithelial cells. Cells from the same sample are not independent; anti-conservative.
- 24 neighborhood Spearmans were BH-adjusted together. None pass FDR < 0.05.
- Genes missing from the log2TPM matrix: none.

## Honest limits

- Dissociated droplet scRNA-seq: no spatial neighborhood and no histologically scored TLS.
- Immune fractions are author-label proportions within each dissociated sample.
- The 12-chemokine TLS score is a bulk-style proxy, not a follicle.
- Cell-level GSEA inflates n and is secondary.
- Tumor-site n=36 pools primary lung, lymph-node, brain met, and pleural effusion.
- CLDN4–immune inverses are exploratory and fail BH-FDR.
- EMT is not a clean hit: Hallmark EMT FDR=0.098; CUSTOM_EMT_CORE NES=0.78 FDR=0.87.

## Neighborhood associations (sample-level Spearman + BH-FDR)

- tumor_sites: epi TACSTD2 mean vs CD8 fraction: ρ=0.106, p=0.538, FDR=0.89, n=36
- tumor_sites: epi TACSTD2 mean vs NK fraction: ρ=−0.012, p=0.943, FDR=0.96, n=36
- tumor_sites: epi TACSTD2 mean vs B fraction: ρ=−0.027, p=0.874, FDR=0.95, n=36
- tumor_sites: epi TACSTD2 mean vs GC-B fraction: ρ=0.183, p=0.286, FDR=0.76, n=36
- tumor_sites: epi TACSTD2 mean vs T/NK fraction: ρ=0.102, p=0.555, FDR=0.89, n=36
- tumor_sites: epi TACSTD2 mean vs 12-chemokine TLS score: ρ=0.161, p=0.347, FDR=0.76, n=36
- tumor_sites: epi CLDN4 mean vs CD8 fraction: ρ=−0.322, p=0.055, FDR=0.27, n=36
- tumor_sites: epi CLDN4 mean vs NK fraction: ρ=−0.402, p=0.015, FDR=0.15, n=36
- tumor_sites: epi CLDN4 mean vs B fraction: ρ=−0.374, p=0.025, FDR=0.15, n=36
- tumor_sites: epi CLDN4 mean vs GC-B fraction: ρ=−0.126, p=0.464, FDR=0.86, n=36
- tumor_sites: epi CLDN4 mean vs T/NK fraction: ρ=−0.380, p=0.022, FDR=0.15, n=36
- tumor_sites: epi CLDN4 mean vs 12-chemokine TLS score: ρ=−0.274, p=0.106, FDR=0.36, n=36
- tLung: epi TACSTD2 mean vs CD8 fraction: ρ=0.327, p=0.326, FDR=0.76, n=11
- tLung: epi TACSTD2 mean vs NK fraction: ρ=−0.691, p=0.019, FDR=0.15, n=11
- tLung: epi TACSTD2 mean vs B fraction: ρ=0.127, p=0.709, FDR=0.90, n=11
- tLung: epi TACSTD2 mean vs GC-B fraction: ρ=−0.112, p=0.744, FDR=0.90, n=11
- tLung: epi TACSTD2 mean vs T/NK fraction: ρ=0.055, p=0.873, FDR=0.95, n=11
- tLung: epi TACSTD2 mean vs 12-chemokine TLS score: ρ=0.109, p=0.750, FDR=0.90, n=11
- tLung: epi CLDN4 mean vs CD8 fraction: ρ=0.245, p=0.467, FDR=0.86, n=11
- tLung: epi CLDN4 mean vs NK fraction: ρ=−0.555, p=0.077, FDR=0.31, n=11
- tLung: epi CLDN4 mean vs B fraction: ρ=0.018, p=0.958, FDR=0.96, n=11
- tLung: epi CLDN4 mean vs GC-B fraction: ρ=−0.154, p=0.652, FDR=0.90, n=11
- tLung: epi CLDN4 mean vs T/NK fraction: ρ=−0.318, p=0.340, FDR=0.76, n=11
- tLung: epi CLDN4 mean vs 12-chemokine TLS score: ρ=0.109, p=0.750, FDR=0.90, n=11

## Sample-level GSEA (primary)

- CUSTOM_TJ_CORE: NES=2.21, FDR=0
- CUSTOM_TJ_CORE without CLDN4/TACSTD2: NES=2.17, FDR=0
- GOBP_TIGHT_JUNCTION_ORGANIZATION: NES=1.79, FDR=0.0080
- HALLMARK_APICAL_JUNCTION: NES=1.79, FDR=0.0060
- KEGG_TIGHT_JUNCTION: NES=1.74, FDR=0.0075
- CUSTOM_BARRIER_KERATIN: NES=1.50, FDR=0.027
- GOBP_KERATINIZATION: NES=1.48, FDR=0.031
- HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION: NES=1.30, FDR=0.098
- GOBP_CORNIFICATION: NES=−1.06, FDR=0.33
- CUSTOM_EMT_CORE: NES=0.78, FDR=0.87

## Cell-level GSEA (secondary, anti-conservative)

- CUSTOM_TJ_CORE: NES=1.81, FDR=0.0021
- KEGG_TIGHT_JUNCTION: NES=1.54, FDR=0.0046
- GOBP_KERATINIZATION: NES=1.22, FDR=0.11
- CUSTOM_BARRIER_KERATIN: NES=1.16, FDR=0.16
- HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION: NES=1.21, FDR=0.11

## Files

| File | Role |
|---|---|
| `summary.json` | Machine-readable verdict |
| `epithelial_vs_immune.tsv` | TACSTD2/CLDN4 epithelial vs T/NK/B/myeloid |
| `per_sample.tsv` | Sample-level epithelial scores and immune fractions |
| `neighborhood_associations.tsv` | Spearman tests + BH-FDR |
| `gene_rank_metrics.tsv.gz` | Per-gene sample Spearman and cell t |
| `gsea_sample_spearman_TACSTD2.tsv` | Primary prerank GSEA |
| `gsea_cell_t_TROP2_Q4vsQ1.tsv` | Secondary cell-level GSEA |
| `gsea_all.tsv` / `gsea_key_sets_sample.tsv` | Combined / key-set slices |
| `fig_epithelial_vs_immune.png` | Restriction boxplots |
| `fig_neighborhood_cooccurrence.png` | TACSTD2 vs CD8/NK/B/TLS (null) |
| `fig_neighborhood_cldn4.png` | CLDN4 vs CD8/NK/B/TLS (nominal) |
| `fig_gsea_sample_TACSTD2.png` | Primary NES bars |
| `fig_gsea_cell_TROP2Q4Q1.png` | Secondary NES bars |
| `per_cell_selected_genes.csv.gz` | Streamed marker genes |
| `sample_metadata.tsv` / `cell_type_composition.tsv` | GEO + author labels |

## Reproduce

```bash
python3 scripts/hunt_gse131907/download.py
python3 scripts/hunt_gse131907/analyze.py
```
