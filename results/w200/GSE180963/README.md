# GSE180963: Tacstd2 / Cldn4 in KL vs K lung GEMM scRNA (w200)

**Status.** Public processed 10x matrices **are available** on GEO and were
downloaded. Analysis script is in `scripts/w200_gse180963/`. This writeup is
updated after the run with real numbers.

**Question.** In the public KrasG12D/+ (K) vs KrasG12D/+;Lkb1fl/fl (KL) lung
GEMM scRNA series `GSE180963`, what is the expression of `Tacstd2` (TROP2) and
`Cldn4` in tumor epithelium vs immune cells, and is there any honest
genotype-level difference given **n = 2**?

## Honest constraints (read first)

1. **n = 2 samples, 1 mouse per genotype.** GSM5481386 = K, GSM5481387 = KL.
   The two samples were mixed into one 10x library and demultiplexed by label.
   Genotype is fully confounded with mouse and library. Any K-vs-KL p-value is
   cell-level **pseudoreplication**, not a biological-replicate test. We still
   extract the genes because that is what this series can support.
2. **Processed matrices exist.** GEO supplementary `GSE180963_RAW.tar` (88 MB)
   contains author-filtered count matrices (`K/{matrix.mtx,genes.tsv,barcodes.tsv}`
   and the KL counterpart). No raw/unfiltered matrices. We do **not** reprocess
   SRA FASTQs.
3. **Author QC already applied** (Seurat 3.1.5): 500–6000 features, mito < 20%,
   genes in ≥3 cells. 20,304 genes × 6,696 K cells + 7,564 KL cells.
4. **`Cldn4` is scored as a target, not used to call epithelium** (avoids
   circularity). Epithelial calling uses `Epcam` / keratins / `Cldn18` / `Cdh1`.
5. No private 8-KL-mouse cohort. Public GEO only.

## Dataset

| GSM | label | genotype | cells (author-filtered) |
|---|---|---|---|
| GSM5481386 | K | KrasG12D/+ | 6,696 |
| GSM5481387 | KL | KrasG12D/+;Lkb1fl/fl | 7,564 |

Series title: *Single cell RNA sequencing of tumor sections from GEMM harboring
KrasG12D/+ or KrasG12D/+Lkb1fl/fl (KL) mutation.* Summary explicitly frames
LKB1 loss as producing an immune-desert / “cold” TME. Platform GPL24247
(Illumina, mm10), 10x Chromium 3′ v3.1, FVB strain. Contact: Xue Bai, Southern
Medical University. Public 2022-07-16.

## Methods (to be executed)

`scripts/w200_gse180963/01_download.py` → GEO FTP supplementary tar.
`scripts/w200_gse180963/02_analyze.py` → load MTX → re-apply author QC →
CP10k + log1p → HVG/PCA/neighbors/UMAP/Leiden → marker lineage annotation →
`Tacstd2` and `Cldn4` %positive and mean log-norm by compartment and genotype →
cluster-free epithelial rule as a robustness check → Tacstd2–Cldn4
co-expression in epithelium. **No batch integration** (would erase the
genotype=sample contrast). Seed 0.

## Results

*Filled after `02_analyze.py` completes. Do not treat this draft as a result.*

## Reproduce

```bash
python3 scripts/w200_gse180963/01_download.py
python3 scripts/w200_gse180963/02_analyze.py
```

Requires scanpy, anndata, leidenalg, python-igraph, pandas, scipy, numpy,
matplotlib.
