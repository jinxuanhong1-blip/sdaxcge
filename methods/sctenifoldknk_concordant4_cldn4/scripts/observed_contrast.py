#!/usr/bin/env python3
"""Sanity-check malignant counts against the locked 65 units, then refit the family OLS.

This does not replace the locked patient-level CLDN4 vs T/NK result.
Positive logFC means higher in CLDN4-high (Q4) than CLDN4-low (Q1).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from contrastlib import (
    LOCKED_UNITS,
    ROOT,
    align_zero_fill,
    load_pseudobulk,
    load_units,
    log2_tmm_cpm,
    ols_q4,
)

SETS = json.loads((ROOT / "data" / "gene_sets.json").read_text())
FAMILY = {
    "IFN": SETS["IFN"],
    "MHC-I/APM": SETS["MHC_I_APM"],
    "TJ": SETS["TJ"],
}


def sanity(units: pd.DataFrame) -> pd.DataFrame:
    locked = pd.read_csv(LOCKED_UNITS, sep="\t")
    locked["unit_id"] = locked["unit_id"].astype(str)
    units = units.copy()
    units["unit_id"] = units["unit_id"].astype(str)
    m = units.merge(locked, on=["dataset", "unit_id"], suffixes=("_new", "_locked"))
    if len(m) != len(locked) or len(m) != len(units):
        raise SystemExit(f"unit id mismatch: new {len(units)} locked {len(locked)} joined {len(m)}")
    rows = []
    for col in ["n_cells", "n_malignant", "n_tnk", "mal_CLDN4_pct", "mal_CLDN4_mean"]:
        a = m[f"{col}_new"].to_numpy(dtype=float)
        b = m[f"{col}_locked"].to_numpy(dtype=float)
        rows.append(
            {
                "field": col,
                "max_abs_diff": float(np.nanmax(np.abs(a - b))),
                "n_mismatch_gt_tol": int(np.sum(np.abs(a - b) > (0.05 if "CLDN4_pct" in col else 1e-3 if "mean" in col else 0))),
            }
        )
    out = pd.DataFrame(rows)
    detail = m[
        [
            "dataset",
            "unit_id",
            "n_malignant_new",
            "n_malignant_locked",
            "mal_CLDN4_pct_new",
            "mal_CLDN4_pct_locked",
            "mal_CLDN4_mean_new",
            "mal_CLDN4_mean_locked",
            "cldn4_quartile",
        ]
    ].copy()
    return out, detail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, default=Path("/tmp/c4_work"))
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "tables")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    units = load_units(args.work)
    summary, detail = sanity(units)
    summary.to_csv(args.out / "unit_sanity_summary.tsv", sep="\t", index=False)
    detail.to_csv(args.out / "unit_sanity_detail.tsv", sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)
    bad = summary.loc[summary["n_mismatch_gt_tol"] > 0]
    if len(bad):
        raise SystemExit("sanity failed; refusing to treat this extract as the locked contrast")

    frames = load_pseudobulk(args.work)
    genes, counts, meta, _panel = align_zero_fill(frames)
    locked = pd.read_csv(LOCKED_UNITS, sep="\t")
    locked["unit_id"] = locked["unit_id"].astype(str)
    meta["unit_id"] = meta["unit_id"].astype(str)
    meta = meta.merge(
        locked[["dataset", "unit_id", "n_malignant", "cldn4_quartile"]],
        on=["dataset", "unit_id"],
        how="left",
    )
    if meta["cldn4_quartile"].isna().any():
        raise SystemExit("missing locked quartile")
    keep = meta["n_malignant"].to_numpy() >= 30
    counts = counts[:, keep]
    meta = meta.loc[keep].reset_index(drop=True)
    logcpm, sf = log2_tmm_cpm(counts)
    qmask = meta["cldn4_quartile"].isin(["Q1", "Q4"]).to_numpy()
    y_ds = meta.loc[qmask, "dataset"].to_numpy()
    y_q4 = (meta.loc[qmask, "cldn4_quartile"] == "Q4").to_numpy()
    rows = []
    for name, members in FAMILY.items():
        present = [g for g in members if g in set(genes)]
        ix = [genes.index(g) for g in present]
        score = logcpm[ix].mean(axis=0)
        logfc, p, df, rank = ols_q4(score[qmask], y_ds, y_q4)
        rows.append(
            {
                "family": name,
                "n_q1": int((meta.loc[qmask, "cldn4_quartile"] == "Q1").sum()),
                "n_q4": int((meta.loc[qmask, "cldn4_quartile"] == "Q4").sum()),
                "n_genes_in_set": len(members),
                "n_genes_present": len(present),
                "logFC_high_minus_low": logfc,
                "p": p,
                "df": df,
                "rank": rank,
            }
        )
        print(name, "logFC", logfc, "p", p, "genes", len(present), flush=True)
    pd.DataFrame(rows).to_csv(args.out / "family_q4q1_ols.tsv", sep="\t", index=False)
    np.savez_compressed(
        args.work / "logcpm_de_units.npz",
        genes=np.asarray(genes),
        logcpm=logcpm,
        dataset=meta["dataset"].to_numpy(),
        unit_id=meta["unit_id"].to_numpy(),
        quartile=meta["cldn4_quartile"].to_numpy(),
        n_malignant=meta["n_malignant"].to_numpy(),
        tmm=sf,
    )
    print("wrote logcpm", logcpm.shape, flush=True)


if __name__ == "__main__":
    main()
