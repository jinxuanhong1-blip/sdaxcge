#!/usr/bin/env python3
"""Analyze TACSTD2 / CLDN4 expression in the 8-LUAD scRNA-seq matrix (Zenodo 10731914).

Matrix layout: genes in rows, cells in columns (raw UMI counts).
Steps:
  1. Stream the gzipped matrix, pull rows for genes of interest + reference markers.
  2. Join to per-cell annotation (Cell_type / Cell_subtype).
  3. Report mean expression and % of cells expressing (count>0) per cell type,
     using CP10K normalisation for the mean (per-cell library size from nCount_RNA).
  4. Save tidy CSV tables and a bar/dot summary figure.
"""
import gzip
import csv
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
DL = ROOT / "results" / "fable_zenodo" / "downloads" / "10731914"
OUT = ROOT / "results" / "fable_zenodo"
MAT = DL / "cell_matrix_sce_8LUADs.csv.gz"
ANN = DL / "cell_annotation_sce_8LUADs.csv"

GENES_OF_INTEREST = ["TACSTD2", "CLDN4"]
REFERENCE = ["EPCAM", "PTPRC", "CD3D", "SFTPC", "KRT19", "KRT18"]
WANT = set(GENES_OF_INTEREST + REFERENCE)


def load_annotation():
    ann = pd.read_csv(ANN)
    ann = ann.rename(columns={ann.columns[0]: "barcode"})
    ann["barcode"] = ann["barcode"].astype(str)
    return ann


def stream_gene_rows(path, wanted):
    """Return dict gene -> np.array(counts) and the ordered list of cell barcodes."""
    found = {}
    with gzip.open(path, "rt") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        barcodes = header[1:]
        for row in reader:
            g = row[0]
            if g in wanted:
                found[g] = np.array(row[1:], dtype=np.float64)
                if len(found) == len(wanted):
                    break
    return found, barcodes


def main():
    ann = load_annotation()
    gene_counts, barcodes = stream_gene_rows(MAT, WANT)
    print("Genes found:", sorted(gene_counts))
    missing = WANT - set(gene_counts)
    if missing:
        print("Genes NOT in matrix:", sorted(missing))

    counts = pd.DataFrame({g: gene_counts[g] for g in gene_counts}, index=barcodes)
    counts.index.name = "barcode"

    # Align annotation to matrix column order
    ann_idx = ann.set_index("barcode")
    common = [b for b in barcodes if b in ann_idx.index]
    print(f"Cells in matrix: {len(barcodes)}; with annotation: {len(common)}")
    counts = counts.loc[common]
    meta = ann_idx.loc[common]

    lib = meta["nCount_RNA"].astype(float).replace(0, np.nan)
    cp10k = counts.div(lib, axis=0) * 1e4  # normalised expression

    # Choose the most informative cell-type column
    ct_col = "Cell_type" if "Cell_type" in meta.columns else meta.columns[0]
    # 'Cell_type' has many 'nan' strings; prefer 'Cell_type.unimodel' if richer
    for cand in ["Cell_type.unimodel", "Cell_type"]:
        if cand in meta.columns and meta[cand].astype(str).ne("nan").mean() > 0.5:
            ct_col = cand
            break
    print("Using cell-type column:", ct_col)
    meta = meta.copy()
    meta["_ct"] = meta[ct_col].astype(str)

    rows = []
    for gene in GENES_OF_INTEREST + [g for g in REFERENCE if g in counts.columns]:
        for ct, idx in meta.groupby("_ct").groups.items():
            sub_raw = counts.loc[idx, gene]
            sub_norm = cp10k.loc[idx, gene]
            n = len(idx)
            pct = float((sub_raw > 0).mean() * 100)
            rows.append({
                "gene": gene,
                "cell_type": ct,
                "n_cells": int(n),
                "pct_expressing": round(pct, 2),
                "mean_cp10k": round(float(np.nanmean(sub_norm)), 4),
                "mean_raw": round(float(sub_raw.mean()), 4),
            })
    tab = pd.DataFrame(rows).sort_values(["gene", "mean_cp10k"], ascending=[True, False])
    tab.to_csv(OUT / "luad_expression_by_celltype.csv", index=False)
    print("\n=== TACSTD2 / CLDN4 by cell type (top) ===")
    for g in GENES_OF_INTEREST:
        print(f"\n-- {g} --")
        print(tab[tab.gene == g].head(12).to_string(index=False))

    # Global summary for the two genes
    glob = []
    for gene in GENES_OF_INTEREST + [g for g in REFERENCE if g in counts.columns]:
        glob.append({
            "gene": gene,
            "n_cells": int(counts.shape[0]),
            "pct_expressing": round(float((counts[gene] > 0).mean() * 100), 2),
            "mean_cp10k": round(float(np.nanmean(cp10k[gene])), 4),
            "max_raw": float(counts[gene].max()),
        })
    pd.DataFrame(glob).to_csv(OUT / "luad_expression_global.csv", index=False)

    # Figure: dot plot (size=%expressing, color=mean cp10k) for the 2 genes across cell types
    fig_genes = GENES_OF_INTEREST
    cts = (tab[tab.gene == "TACSTD2"].sort_values("mean_cp10k", ascending=False)
           ["cell_type"].tolist())
    fig, ax = plt.subplots(figsize=(7, max(3, 0.4 * len(cts))))
    for xi, gene in enumerate(fig_genes):
        sub = tab[tab.gene == gene].set_index("cell_type").reindex(cts)
        sizes = sub["pct_expressing"].fillna(0) * 6
        colors = sub["mean_cp10k"].fillna(0)
        sc = ax.scatter([xi] * len(cts), range(len(cts)), s=sizes,
                        c=colors, cmap="Reds", edgecolors="grey", linewidths=0.5)
    ax.set_xticks(range(len(fig_genes)))
    ax.set_xticklabels(fig_genes)
    ax.set_yticks(range(len(cts)))
    ax.set_yticklabels(cts, fontsize=8)
    ax.invert_yaxis()
    ax.set_title("TACSTD2 / CLDN4 in 8-LUAD scRNA-seq (Zenodo 10731914)\nsize=% expressing, color=mean CP10K")
    fig.colorbar(sc, ax=ax, label="mean CP10K", shrink=0.6)
    fig.tight_layout()
    fig.savefig(OUT / "luad_tacstd2_cldn4_dotplot.png", dpi=150)
    print(f"\nWrote figure {OUT / 'luad_tacstd2_cldn4_dotplot.png'}")
    print(f"Wrote {OUT / 'luad_expression_by_celltype.csv'} and luad_expression_global.csv")


if __name__ == "__main__":
    main()
