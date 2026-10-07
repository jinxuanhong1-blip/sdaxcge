#!/usr/bin/env python3
"""Assemble streamed E-MTAB-16433 panel rows (SG/Trodelvy vs vehicle PDOX tumours).

Pseudobulk unit is the mouse, so each arm has 4 independent biological replicates.
"""
import glob
import os

import numpy as np
import pandas as pd

STREAM = "/tmp/ae/stream433"
META = "/tmp/ae/treatment_HD4246_SG_metadata.tsv"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)

barcodes = open(f"{STREAM}/header.tsv").readline().rstrip("\n").split("\t")[1:]
n = len(barcodes)

rows, genes = {}, []
for f in sorted(glob.glob(f"{STREAM}/genes_*.tsv")):
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        g, vals = p[0], p[1:]
        assert len(vals) == n, (g, len(vals), n)
        if g in rows:
            continue
        v = np.asarray(vals, dtype=np.float64)
        assert np.all(v == np.round(v)), f"{g} non-integer"
        rows[g] = v.astype(np.int32)
        genes.append(g)
genes.sort()
panel = pd.DataFrame({g: rows[g] for g in genes}, index=barcodes)

libsize = np.zeros(n, dtype=np.int64)
for f in glob.glob(f"{STREAM}/colsum_*.tsv"):
    a = np.loadtxt(f, dtype=np.int64, ndmin=2)
    libsize[a[:, 0] - 1] += a[:, 1]
assert (libsize > 0).all()

meta = pd.read_csv(META, sep="\t").set_index("barcode").reindex(barcodes)
assert meta["treatment"].notna().all(), "barcode/metadata mismatch"

cells = meta[["sample_id", "treatment", "mouse_ID", "batch", "leiden"]].copy()
cells["total_counts"] = libsize
cells = pd.concat([cells, panel], axis=1)
cells.index.name = "barcode"
cells.to_csv(f"{OUT}/emtab16433_panel_per_cell.csv.gz", compression="gzip")

grp = cells.groupby("mouse_ID")
pb = grp[genes].sum()
info = pd.DataFrame({
    "n_cells": grp.size(), "total_counts": grp["total_counts"].sum(),
    "treatment": grp["treatment"].first(), "batch": grp["batch"].first(),
})
cpm = pb.div(info["total_counts"], axis=0) * 1e6
out = pd.concat([info, cpm.add_suffix("__cpm"), pb.add_suffix("__count")], axis=1)
out.index.name = "mouse_ID"
out.to_csv(f"{OUT}/emtab16433_pseudobulk_cpm.csv")

print(f"cells={n} genes={len(genes)} mice={len(out)}")
print(info.sort_values(["treatment", "batch"]).to_string())
print("\nmedian library size:", int(np.median(libsize)))
