#!/usr/bin/env python3
"""Extract a target + reference gene panel from the GSE285888 scRNA-seq matrix.

The full matrix (GSM8712033_matrix.txt.gz) is a dense genes x cells text table
(~28.9k genes x 222,144 cells). Loading it fully is memory-heavy, so we stream
it line-by-line and keep only the rows for a small gene panel, then join the
per-cell response / irAE annotation from the metadata file.

Outputs (intermediate, written to $GEO_DIR, NOT committed):
    gse285888_panel_cells.parquet   cells x (panel genes + annotation)
"""
import gzip
import os
import sys
import pandas as pd
import numpy as np

GEO_DIR = os.environ.get("GEO_DIR", "/tmp/geo")
MATRIX = os.path.join(GEO_DIR, "GSM8712033_matrix.txt.gz")
META = os.path.join(GEO_DIR, "GSM8712033_metadata.csv.gz")
OUT = os.path.join(GEO_DIR, "gse285888_panel_cells.parquet")

# Targets of interest + reference markers to contextualise detectability.
TARGETS = ["TACSTD2", "CLDN4"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18"]          # epithelial / tumour context
IMMUNE = ["PTPRC", "CD3D", "CD8A", "MS4A1", "CD14", "NKG7"]  # PBMC lineages
EFFECTOR = ["GZMB", "PRF1", "IL1B", "CXCL8"]     # markers highlighted by the study
PANEL = TARGETS + EPITHELIAL + IMMUNE + EFFECTOR


def main():
    if not os.path.exists(MATRIX):
        sys.exit(f"Missing {MATRIX}; run 01_download.sh first.")

    panel_set = set(PANEL)
    with gzip.open(MATRIX, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        barcodes = [b.strip('"') for b in header]
        n_cells = len(barcodes)
        print(f"cells (header columns): {n_cells}", flush=True)

        collected = {}
        for i, line in enumerate(fh, start=1):
            tab = line.find("\t")
            if tab == -1:
                continue
            gene = line[:tab].strip('"')
            if gene not in panel_set:
                continue
            vals = line.rstrip("\n").split("\t")[1:]
            if len(vals) != n_cells:
                sys.exit(f"Field mismatch for {gene}: {len(vals)} vs {n_cells}")
            collected[gene] = np.asarray(vals, dtype=np.float32)
            print(f"  extracted {gene} ({len(collected)}/{len(PANEL)})", flush=True)
            if len(collected) == len(PANEL):
                # keep scanning is unnecessary once all found
                pass
    missing = panel_set - set(collected)
    if missing:
        print(f"WARNING: genes absent from matrix feature list: {sorted(missing)}",
              flush=True)

    expr = pd.DataFrame(collected, index=barcodes)
    expr.index.name = "cell"

    meta = pd.read_csv(META, index_col=0)
    meta.index.name = "cell"
    keep_meta = ["orig.ident", "nCount_RNA", "nFeature_RNA", "SubType",
                 "SampleID", "irAE"]
    meta = meta[keep_meta]

    df = meta.join(expr, how="inner")
    print(f"joined cells: {df.shape[0]} (meta {meta.shape[0]}, expr {expr.shape[0]})",
          flush=True)

    df.to_parquet(OUT)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
