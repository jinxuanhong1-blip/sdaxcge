# Claim A8/A9 — TROP2-high keratin / tight-junction GSEA

Honest overlap of an a priori tight-junction panel
(**CLDN1, CLDN4, CLDN7, F11R, PARD3**) with **TACSTD2 (TROP2)**
across TCGA-LUAD, TCGA-LUSC, GSE207422, and GSE131907.

Outputs live in [`results/claim_A8A9/`](results/claim_A8A9/).
Start with [`results/claim_A8A9/REPORT.md`](results/claim_A8A9/REPORT.md)
and [`results/claim_A8A9/intersection_summary.json`](results/claim_A8A9/intersection_summary.json).

## What is computed

Same rules in every dataset:

1. Rank every gene by Spearman correlation with TACSTD2.
2. Pre-ranked GSEA (gseapy, 1000 permutations) on a local GMT of
   keratin / keratinization and tight-junction gene sets
   (Enrichr GO BP/CC 2021, KEGG 2021 Human, Reactome 2022;
   see `data/genesets/claim_A8A9_genesets.provenance.json`).
3. TROP2-high vs TROP2-low = top vs bottom tertile of TACSTD2.
4. A panel gene is `up_in_trop2_high` only if rho > 0 **and**
   log2FC(high vs low) > 0 **and** Welch BH-FDR < 0.05.

Intersection is the subset of those five genes that pass in the
named datasets. Empty intersections are reported as empty — they
are negatives, not missing analyses.

## Datasets

| Dataset | Assay | Role |
|---|---|---|
| TCGA-LUAD | STAR log2(TPM+1), primary tumour, one sample/patient | primary bulk |
| GSE207422 | author log2TPM, n=24 pre-treatment NSCLC biopsies | primary bulk |
| TCGA-LUSC | same as LUAD | squamous histology check |
| GSE131907 | LUAD scRNA-seq, epithelial cells | supporting (tLung pseudobulk) + exploratory cell-level |

Raw matrices are not in git (`data/raw/` is gitignored). They are
public: UCSC Xena GDC hub `TCGA-{LUAD,LUSC}.star_tpm.tsv.gz` and
GEO supplementary files for GSE207422 / GSE131907.

## Reproduce

```bash
python3 -m pip install pandas numpy scipy matplotlib openpyxl gseapy
# place the public matrices in data/raw/ (see scripts for filenames)
PYTHONPATH=scripts python3 scripts/00_build_genesets.py
PYTHONPATH=scripts python3 scripts/01_tcga.py LUAD
PYTHONPATH=scripts python3 scripts/01_tcga.py LUSC
PYTHONPATH=scripts python3 scripts/02_gse207422.py
PYTHONPATH=scripts python3 scripts/03_gse131907.py
PYTHONPATH=scripts python3 scripts/04_intersect_and_report.py
```
