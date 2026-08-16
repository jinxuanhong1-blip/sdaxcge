#!/usr/bin/env python3
"""Stream TACSTD2/CLDN4 from GSE241934 MTX and write per-sample tables.

Public GEO processed files only. Author cell types and Histology are used as-is.
Does not invent LUAD/LUSC labels. ASC is kept as ASC.
"""

from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


GENES = ("TACSTD2", "CLDN4")


def load_feature_rows(path: Path) -> dict[str, int]:
    rows = {}
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split("\t")
            name = parts[1] if len(parts) > 1 else parts[0]
            if name in GENES:
                rows[name] = i
    missing = [g for g in GENES if g not in rows]
    if missing:
        raise SystemExit(f"Missing genes in {path}: {missing}")
    return rows


def stream_mtx(path: Path, keep_rows: dict[str, int], n_cells: int) -> dict[str, np.ndarray]:
    want = {idx: gene for gene, idx in keep_rows.items()}
    out = {gene: np.zeros(n_cells, dtype=np.float32) for gene in keep_rows}
    with gzip.open(path, "rt") as fh:
        header = fh.readline()
        if not header.startswith("%%MatrixMarket"):
            raise SystemExit(f"Not MTX: {path}")
        dims = fh.readline()
        while dims.startswith("%"):
            dims = fh.readline()
        n_genes, n_cols, _nnz = map(int, dims.split())
        if n_cols != n_cells:
            raise SystemExit(f"{path}: MTX cells {n_cols} != barcodes {n_cells}")
        for line in fh:
            r_s, c_s, v_s = line.split()
            r = int(r_s)
            if r in want:
                out[want[r]][int(c_s) - 1] = float(v_s)
    return out


def per_sample(meta: pd.DataFrame, counts: dict[str, np.ndarray], cohort: str) -> pd.DataFrame:
    meta = meta.reset_index(drop=True)
    ncount = meta["nCount_RNA"].to_numpy(dtype=float)
    ncount = np.where(ncount > 0, ncount, np.nan)
    expr = {}
    for gene, umi in counts.items():
        cp10k = 1e4 * umi / ncount
        expr[f"{gene}_log1p_cp10k"] = np.log1p(cp10k)
        expr[f"{gene}_pos"] = umi > 0
    lineage = meta["major.cell.type"].astype(str)
    is_epi = lineage.eq("Epi")
    is_tnk = lineage.isin(["T", "NK"])
    rows = []
    for sample, idx in meta.groupby("sampleID", sort=True).groups.items():
        idx = np.asarray(list(idx))
        sub = meta.iloc[idx]
        epi = is_epi.iloc[idx].to_numpy()
        tnk = is_tnk.iloc[idx].to_numpy()
        rec = {
            "cohort": cohort,
            "sample": sample,
            "histology_raw": sub["Histology"].iloc[0],
            "path_response": sub["Pathological Response"].iloc[0],
            "egfr": sub["EGFR"].iloc[0],
            "n_cells": int(len(idx)),
            "n_epi": int(epi.sum()),
            "n_tnk": int(tnk.sum()),
            "frac_tnk": float(tnk.mean()) if len(idx) else np.nan,
        }
        for gene in GENES:
            vals = expr[f"{gene}_log1p_cp10k"][idx]
            pos = expr[f"{gene}_pos"][idx]
            rec[f"epi_{gene}_mean_log1p_cp10k"] = float(np.nanmean(vals[epi])) if epi.any() else np.nan
            rec[f"epi_{gene}_pct_pos"] = float(pos[epi].mean() * 100) if epi.any() else np.nan
            rec[f"tnk_{gene}_mean_log1p_cp10k"] = float(np.nanmean(vals[tnk])) if tnk.any() else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def run_one(meta_path: Path, feat_path: Path, mtx_path: Path, barcodes_path: Path, cohort: str) -> pd.DataFrame:
    barcodes = [ln.strip() for ln in gzip.open(barcodes_path, "rt")]
    meta = pd.read_csv(meta_path, sep="\t")
    if list(meta["cellID"]) != barcodes:
        raise SystemExit(f"{cohort}: cellID order does not match barcodes")
    keep = load_feature_rows(feat_path)
    counts = stream_mtx(mtx_path, keep, n_cells=len(barcodes))
    return per_sample(meta, counts, cohort)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", default="/tmp/gse241934")
    p.add_argument("--out", required=True)
    args = p.parse_args()
    data = Path(args.data_dir)
    iit = run_one(
        data / "IIT_Meta.txt.gz",
        data / "IIT_features.tsv.gz",
        data / "IIT_Matrix.mtx.gz",
        data / "IIT_barcodes.tsv.gz",
        "GSE241934_IIT",
    )
    rwc = run_one(
        data / "Real_Meta.txt.gz",
        data / "IIT_features.tsv.gz",  # same 27693 genes
        data / "Real_Matrix.mtx.gz",
        data / "RWC_barcodes.tsv.gz",
        "GSE241934_RWC",
    )
    out = pd.concat([iit, rwc], ignore_index=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out} n={len(out)}")


if __name__ == "__main__":
    main()
