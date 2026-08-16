#!/usr/bin/env python3
"""Leftover CosMx GSE276083 (mycobacterial lung). Public files lack x/y; FOV is the spatial unit."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import B_GENES, BROAD_EPI, EPI, T_GENES, dump_json, present, spearman  # noqa: E402

OUT = ROOT / "results" / "leftover_spatial"
WANTED = list(dict.fromkeys(EPI + BROAD_EPI + T_GENES + B_GENES + ["PTPRC", "CD68"]))


def load_wanted_rows(path):
    rows = {}
    header = None
    with gzip.open(path, "rt") as f:
        header = [c.strip().strip('"') for c in f.readline().rstrip("\n").split(",")]
        for line in f:
            gene = line.split(",", 1)[0].strip().strip('"').upper()
            if gene in {w.upper() for w in WANTED}:
                vals = line.rstrip("\n").split(",")
                rows[gene] = np.array([float(x) if x not in {"", "NA"} else np.nan for x in vals[1:]], dtype=float)
    cells = header[1:]
    return cells, rows


def main():
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    cells, rows = load_wanted_rows(ROOT / "data" / "GSE276083" / "GSE276083_CosMx_raw_counts.csv.gz")
    ann = pd.read_csv(ROOT / "data" / "GSE276083" / "GSE276083_CosMx_annotations.csv.gz")
    # join on cell name
    key = ann.columns[0]
    ann[key] = ann[key].astype(str)
    cell_s = pd.Series(cells, dtype=str)
    # map annotation fov
    fov_map = ann.set_index(key)["fov"] if "fov" in ann.columns else None
    ncount = ann.set_index(key)["nCount_Nanostring"] if "nCount_Nanostring" in ann.columns else None

    expr = pd.DataFrame(rows)
    expr.index = cell_s.to_numpy()
    # QC: require annotation + counts
    keep = expr.index.isin(ann[key])
    expr = expr.loc[keep]
    if ncount is not None:
        nc = ncount.reindex(expr.index)
        expr = expr.loc[nc.fillna(0) >= 20]
    fov = fov_map.reindex(expr.index) if fov_map is not None else pd.Series("all", index=expr.index)

    # same-cell
    def score(df, genes):
        cols = [g for g in genes if g in df.columns]
        if not cols:
            return pd.Series(np.nan, index=df.index)
        return df[cols].mean(axis=1)

    cl4 = expr["CLDN4"] if "CLDN4" in expr.columns else pd.Series(np.nan, index=expr.index)
    tac = expr["TACSTD2"] if "TACSTD2" in expr.columns else pd.Series(np.nan, index=expr.index)
    tsc = score(expr, T_GENES)
    bsc = score(expr, B_GENES)
    tb = pd.concat([tsc, bsc], axis=1).mean(axis=1)
    same = {
        "n_cells": int(len(expr)),
        "frac_CLDN4_pos": float((cl4 > 0).mean()) if cl4.notna().any() else None,
        "frac_TACSTD2_pos": float((tac > 0).mean()) if tac.notna().any() else None,
        "same_CLDN4_TB": spearman(cl4, tb),
        "same_TACSTD2_TB": spearman(tac, tb),
        "same_CLDN4_T": spearman(cl4, tsc),
        "same_TACSTD2_T": spearman(tac, tsc),
        "panel_genes_used": sorted(expr.columns.tolist()),
    }

    # FOV-level means (the only public spatial unit; x/y not deposited)
    tmp = pd.DataFrame({"fov": fov.to_numpy(), "CLDN4": cl4.to_numpy(), "TACSTD2": tac.to_numpy(), "T": tsc.to_numpy(), "B": bsc.to_numpy(), "TB": tb.to_numpy()})
    fov_df = tmp.groupby("fov").mean(numeric_only=True)
    fov_n = tmp.groupby("fov").size()
    fov_df = fov_df.loc[fov_n >= 20]
    fov_stats = {
        "n_fov": int(len(fov_df)),
        "median_cells_per_fov": float(fov_n.loc[fov_df.index].median()) if len(fov_df) else None,
        "CLDN4_vs_TB": spearman(fov_df["CLDN4"], fov_df["TB"]),
        "TACSTD2_vs_TB": spearman(fov_df["TACSTD2"], fov_df["TB"]),
        "CLDN4_vs_T": spearman(fov_df["CLDN4"], fov_df["T"]),
        "TACSTD2_vs_T": spearman(fov_df["TACSTD2"], fov_df["T"]),
        "note": "GSE276083_CosMx_spatial_coordinates.csv.gz contains only Width/Height; no CenterX/Y. Neighborhood is FOV-level.",
    }
    fov_df.assign(n_cells=fov_n.reindex(fov_df.index)).to_csv(OUT / "tables" / "gse276083_fov_means.csv")
    payload = {"series": "GSE276083", "platform": "CosMx SMI", "same_cell": same, "fov": fov_stats}
    dump_json(OUT / "tables" / "cosmx_summary.json", payload)
    print(payload)


if __name__ == "__main__":
    main()
