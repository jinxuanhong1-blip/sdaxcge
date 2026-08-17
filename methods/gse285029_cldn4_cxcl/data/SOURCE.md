# GSE285029 processed source

Public author WTS matrix only. No FASTQ / SRA. No RECIST / histology / PD-L1 IHC / TMB / purity on GEO.

| File | URL |
|---|---|
| `GSE285029_WTS_expr_count_235_032820.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285029/suppl/GSE285029_WTS_expr_count_235_032820.txt.gz |

Series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029

Paper: Koh et al., *J Immunother Cancer* 2025; PMID 40050048.

Filename says “count” and 235 columns (gene id + **234** samples, Case1–Case234). Values are author-processed (not raw integer counts; not TPM). Primary transform in this folder: `log2(pmax(x,0)+1)`.

CLDN4 only. TACSTD2 is not an anchor here.

IFN-compact Spearman ρ vs CLDN4 is already known from the sibling GSE285029 score folder and is **not** re-claimed here. This folder is chemokine extra: CXCL9, CXCL10, CXCL13, and GEP-like (Ayers 2017 GEP18, unweighted z-mean).
