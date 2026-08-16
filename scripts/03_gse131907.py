#!/usr/bin/env python3
"""GSE131907 LUAD scRNA-seq: TROP2 (TACSTD2) vs keratin / tight-junction programs.

Two complementary views, both restricted to epithelial cells:

  1. Sample-level pseudobulk of tLung (primary tumour) epithelial cells.
     Comparable to the bulk datasets, but n is small (~11 samples).
  2. Cell-level ranking among tLung epithelial cells (exploratory; cells are
     not independent, so p-values / FDRs are anti-conservative).

A sensitivity pseudobulk that adds other tumour-site epithelial cells
(tL/B, mLN, mBrain, PE) is also written.

Raw UMI matrix is streamed so the 208k-cell file is never fully loaded.
"""
import gzip
import os
import numpy as np
import pandas as pd
import lib_analysis as L

RAW = "data/raw"
OUT = "results/claim_A8A9/gse131907"
os.makedirs(OUT, exist_ok=True)

TUMOR_ORIGINS_PRIMARY = {"tLung"}
TUMOR_ORIGINS_ALL = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
MIN_CELLS_PSEUDOBULK = 50


def load_annotation():
    ann = pd.read_csv(f"{RAW}/GSE131907_cell_annotation.txt.gz", sep="\t")
    # Index column is the matrix barcode (Barcode_Sample)
    if "Index" in ann.columns:
        ann = ann.set_index("Index")
    return ann


def stream_umi(keep_barcodes):
    """Return genes x cells UMI DataFrame for the requested barcodes."""
    keep = set(keep_barcodes)
    path = f"{RAW}/GSE131907_raw_UMI.txt.gz"
    print(f"Streaming UMI matrix for {len(keep)} cells ...")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        idx = [i for i, b in enumerate(barcodes) if b in keep]
        names = [barcodes[i] for i in idx]
        print(f"  matched {len(idx)} / {len(keep)} requested barcodes")
        genes, rows = [], []
        for n, line in enumerate(fh, 1):
            parts = line.rstrip("\n").split("\t")
            gene = parts[0]
            vals = np.fromiter((parts[i + 1] for i in idx), dtype=np.float32, count=len(idx))
            genes.append(gene)
            rows.append(vals)
            if n % 5000 == 0:
                print(f"  ... {n} genes")
    mat = pd.DataFrame(np.vstack(rows), index=genes, columns=names)
    print("  UMI matrix:", mat.shape)
    return mat


def log2_cpm(umi, scale=1e4):
    lib = umi.sum(axis=0).replace(0, np.nan)
    return np.log2(umi.div(lib, axis=1) * scale + 1.0)


def pseudobulk_log2cpm(umi, sample_of_cell, min_cells=MIN_CELLS_PSEUDOBULK):
    counts = pd.Series(sample_of_cell).value_counts()
    keep_samples = [s for s, n in counts.items() if n >= min_cells]
    cols = [c for c in umi.columns if sample_of_cell[c] in keep_samples]
    umi = umi[cols]
    grouped = umi.T.groupby([sample_of_cell[c] for c in umi.columns]).sum().T
    print(f"  pseudobulk samples (min {min_cells} epi cells): {grouped.shape[1]}")
    return log2_cpm(grouped, scale=1e6)


def run_one(label, expr, trop2, outdir, min_frac):
    os.makedirs(outdir, exist_ok=True)
    stats = L.per_gene_stats(expr, trop2, min_frac_expressed=min_frac)
    stats.to_csv(f"{outdir}/per_gene_stats.csv")
    panel = stats.reindex(L.PANEL)
    panel.to_csv(f"{outdir}/panel_stats.csv")
    print(f"--- {label} ---")
    print("  genes", stats.shape[0],
          "n_high/n_low/n_total",
          int(stats["n_high"].iloc[0]), int(stats["n_low"].iloc[0]), int(stats["n_total"].iloc[0]))
    print(panel[["log2FC_high_vs_low", "spearman_rho", "welch_fdr", "up_in_trop2_high"]])
    gsea = L.run_prerank_gsea(stats, outdir)
    print(gsea[["Term", "NES", "NOM p-val", "FDR q-val"]].to_string(index=False))
    return stats, panel, gsea


def main():
    ann = load_annotation()
    epi = ann[ann["Cell_type"] == "Epithelial cells"].copy()
    print("epithelial cells:", epi.shape[0])
    print(epi["Sample_Origin"].value_counts().to_string())

    tlung = epi[epi["Sample_Origin"].isin(TUMOR_ORIGINS_PRIMARY)]
    tumor = epi[epi["Sample_Origin"].isin(TUMOR_ORIGINS_ALL)]
    print("tLung epi:", tlung.shape[0], "all-tumor-site epi:", tumor.shape[0])

    # Stream once for the union of barcodes we need
    umi = stream_umi(list(tumor.index))
    umi = umi.groupby(level=0).sum()

    # ---- tLung cell-level (exploratory) ----
    tlung_cells = [c for c in umi.columns if c in tlung.index]
    umi_t = umi[tlung_cells]
    expr_t = log2_cpm(umi_t)
    trop2_t = expr_t.loc[L.TROP2]
    trop2_t.to_csv(f"{OUT}/tacstd2_tlung_cells.csv")
    run_one("tLung epithelial, cell-level (exploratory)",
            expr_t, trop2_t, f"{OUT}/tlung_cell", min_frac=0.01)

    # ---- tLung sample pseudobulk ----
    sample_of = tlung["Sample"].to_dict()
    pb = pseudobulk_log2cpm(umi_t, sample_of)
    if L.TROP2 not in pb.index:
        raise SystemExit("TACSTD2 missing from GSE131907")
    pb.loc[L.TROP2].to_csv(f"{OUT}/tacstd2_tlung_pseudobulk.csv")
    n_cells = tlung.loc[tlung_cells, "Sample"].value_counts()
    n_cells.to_csv(f"{OUT}/tlung_epithelial_cells_per_sample.csv")
    run_one("tLung epithelial, sample pseudobulk",
            pb, pb.loc[L.TROP2], f"{OUT}/tlung_pseudobulk", min_frac=0.0)

    # ---- all tumour-site sample pseudobulk (sensitivity) ----
    tumor_cells = [c for c in umi.columns if c in tumor.index]
    umi_u = umi[tumor_cells]
    sample_of_u = tumor["Sample"].to_dict()
    pb_u = pseudobulk_log2cpm(umi_u, sample_of_u)
    pb_u.loc[L.TROP2].to_csv(f"{OUT}/tacstd2_alltumor_pseudobulk.csv")
    run_one("all tumour-site epithelial, sample pseudobulk (sensitivity)",
            pb_u, pb_u.loc[L.TROP2], f"{OUT}/alltumor_pseudobulk", min_frac=0.0)
    print("GSE131907 done.")


if __name__ == "__main__":
    main()
