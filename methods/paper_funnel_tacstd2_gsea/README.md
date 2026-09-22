# PAPER FUNNEL — Tacstd2-high vs low malignant enrichment GSEA

Public-only re-run. No private 8KL matrices. No LLC-as-KL.

## Question

In GSE137244 epithelial (GEMM tumor cell lines) and public GEMM bulk tumors, do Tacstd2-high samples enrich tight-junction / cell-adhesion / Claudin-family programs, and do those themes rank **#1–3** among positive-NES pathways?

## Cohorts

| tag | accession | compartment | n |
|---|---|---|---|
| GSE137244_epithelium_cell_lines | GSE137244 | GEMM nodule–derived epithelial cell lines (KP+KL) | 10 |
| GSE164758_GEMM_KL_KP_primary | GSE164758 | in-vivo primary NSCLC, untreated KL+KP | 17 |
| GSE137396_GEMM_KL_KP_nodules | GSE137396 | in-vivo lung tumor nodules, KL+KP | 10 |

## Method

1. Median-split samples on Tacstd2 (anchor gene removed from the ranked list).
2. Welch *t* (high − low) preranked GSEA (`gseapy`, 2000 gene-set permutations, seed 49901).
3. Ranking universe = MSigDB Hallmark Mm + Enrichr KEGG_2019_Mouse + WikiPathways_2019_Mouse + GO TJ/adhesion-focused terms + curated Claudin family.
4. Force-report full NES / NOM p / FDR tables. Highlight ranks of tight-junction, cell-adhesion, and Claudin-family terms among **positive NES** sets.

## Reproduce

```bash
pip install -r requirements.txt
# matrices already under cache/ (or re-download inside analyze.py loaders' expected paths)
python3 analyze.py
```

Outputs: `tables/`, `figures/`, `FINDING.md`.
