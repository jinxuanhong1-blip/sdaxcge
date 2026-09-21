#!/usr/bin/env python3
"""Search core/margin cutoffs, LUSC-6-out, and cytotoxic definitions.

The 50 µm / 150 µm CD8+NK analysis is unchanged. This grid asks which
alternate specifications, if any, make both the absolute count gap and the
high/low fold ratio stronger in the core than at the margin.

Nothing here replaces the locked unstratified ratios 0.36 (50 µm) and 0.52
(100 µm). Selection across the grid is post hoc.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cosmx_core_margin_exclusion as base

MIN_ARM = 40
RADII = (25.0, 50.0, 100.0)
MARGIN_CUTS = (20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0)
CORE_CUTS = (60.0, 80.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0, 400.0)
QUANTILES = ((0.20, 0.80), (0.25, 0.75), (0.33, 0.67), (0.40, 0.60), (0.50, 0.50))

LYMPHOID = (
    "NK",
    "T CD4 memory",
    "T CD4 naive",
    "T CD8 memory",
    "T CD8 naive",
    "Treg",
)
IMMUNE = LYMPHOID + (
    "B-cell",
    "mDC",
    "macrophage",
    "mast",
    "monocyte",
    "neutrophil",
    "pDC",
    "plasmablast",
)


def extract_genes(path: str, names: list[str]) -> dict[str, np.ndarray]:
    import h5py

    with h5py.File(path, "r") as f:
        genes = base._decode(f["var"]["_index"][:])
        want = {nm: genes.index(nm) for nm in names}
        n = int(f["obs"]["_index"].shape[0])
        indptr = f["X"]["indptr"][:]
        indices = f["X"]["indices"]
        data = f["X"]["data"]
        out = {nm: np.zeros(n, dtype=np.float32) for nm in names}
        step = 40000
        for start in range(0, n, step):
            end = min(n, start + step)
            a = int(indptr[start])
            b = int(indptr[end])
            if b <= a:
                continue
            cols = np.asarray(indices[a:b])
            vals = np.asarray(data[a:b])
            for nm, gi in want.items():
                hit = np.flatnonzero(cols == gi)
                if hit.size == 0:
                    continue
                rows = np.searchsorted(indptr[start : end + 1], a + hit, side="right") - 1
                out[nm][start + rows] = vals[hit]
    return out


def prepare(bundle, sample_code, rich_ids, cyto_global):
    obs = bundle["obs"]
    cats = bundle["cats"]
    idx = np.flatnonzero(obs["sample"] == sample_code)
    ct = obs["cell_type"][idx]
    niche = obs["niche"][idx]
    xy = bundle["xy"][idx]
    tumor_ids = base.tumor_codes(cats["cell_type"])
    is_tumor = np.isin(ct, tumor_ids)
    is_rich = np.isin(niche, rich_ids)
    is_cyto = cyto_global[idx]
    if is_tumor.sum() < 100 or is_rich.sum() < 20 or is_cyto.sum() < 10:
        return None
    tum_xy = xy[is_tumor]
    dist, _ = cKDTree(xy[is_rich]).query(tum_xy, k=1, workers=4)
    dist = np.asarray(dist, dtype=np.float64)
    dist[is_rich[is_tumor]] = 0.0
    tree_cyto = cKDTree(xy[is_cyto])
    expr = bundle["cldn4"][idx][is_tumor]
    high = expr > float(np.median(expr))
    counts = {}
    for radius in RADII:
        counts[radius] = np.asarray(
            tree_cyto.query_ball_point(tum_xy, r=radius, return_length=True, workers=4),
            dtype=np.float64,
        )
    patient = cats["patient"][int(obs["patient"][idx][0])]
    return {
        "dist": dist,
        "high": high,
        "counts": counts,
        "patient": patient,
        "n_cyto": int(is_cyto.sum()),
        "n_tumor": int(is_tumor.sum()),
    }


def arm(count, high, mask):
    h = mask & high
    lo = mask & ~high
    nh, nl = int(h.sum()), int(lo.sum())
    if nh < MIN_ARM or nl < MIN_ARM:
        return None
    mh = float(count[h].mean())
    ml = float(count[lo].mean())
    ratio = mh / ml if ml >= 0.05 else np.nan
    return mh, ml, mh - ml, ratio


def wilcox(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 5 or np.allclose(x, 0):
        return np.nan
    return float(stats.wilcoxon(x, alternative="two-sided", zero_method="wilcox").pvalue)


def eval_sections(preps, radius, core_m, margin_m, drop_lusc):
    rows = []
    for name, prep in preps.items():
        if drop_lusc and name == "LUSC-6":
            continue
        got_c = arm(prep["counts"][radius], prep["high"], core_m(prep["dist"]))
        got_m = arm(prep["counts"][radius], prep["high"], margin_m(prep["dist"]))
        if got_c is None or got_m is None:
            continue
        mh_c, ml_c, d_c, r_c = got_c
        mh_m, ml_m, d_m, r_m = got_m
        rows.append({
            "sample": name,
            "patient": prep["patient"],
            "d_core": d_c,
            "d_margin": d_m,
            "r_core": r_c,
            "r_margin": r_m,
            "inter_d": d_c - d_m,
            "inter_r": (r_c - r_m) if np.isfinite(r_c) and np.isfinite(r_m) else np.nan,
            "mh_core": mh_c,
            "ml_core": ml_c,
            "mh_margin": mh_m,
            "ml_margin": ml_m,
            "both": bool(
                d_c < 0 and d_c < d_m and np.isfinite(r_c) and np.isfinite(r_m) and r_c < 1 and r_c < r_m
            ),
        })
    return rows


def summarize(rows, meta):
    if len(rows) < 4:
        return None
    df = pd.DataFrame(rows)
    inter_d = df["inter_d"].to_numpy()
    inter_r = df["inter_r"].to_numpy()
    pat_both = 0
    n_pat = 0
    for _, g in df.groupby("patient"):
        n_pat += 1
        d_c = float(g["d_core"].mean())
        d_m = float(g["d_margin"].mean())
        r_c = float(np.nanmean(g["r_core"]))
        r_m = float(np.nanmean(g["r_margin"]))
        if d_c < 0 and d_c < d_m and np.isfinite(r_c) and np.isfinite(r_m) and r_c < 1 and r_c < r_m:
            pat_both += 1
    n = len(df)
    n_both = int(df["both"].sum())
    mean_inter_d = float(np.mean(inter_d))
    mean_inter_r = float(np.nanmean(inter_r)) if np.isfinite(inter_r).any() else np.nan
    aligned = bool(
        mean_inter_d < 0
        and np.isfinite(mean_inter_r)
        and mean_inter_r < 0
        and float(df["d_core"].mean()) < 0
        and float(np.nanmean(df["r_core"])) < 1
    )
    return {
        **meta,
        "n_sections": n,
        "n_both": n_both,
        "frac_both": n_both / n,
        "n_pat": n_pat,
        "n_pat_both": pat_both,
        "mean_d_core": float(df["d_core"].mean()),
        "mean_d_margin": float(df["d_margin"].mean()),
        "mean_inter_d": mean_inter_d,
        "mean_r_core": float(np.nanmean(df["r_core"])),
        "mean_r_margin": float(np.nanmean(df["r_margin"])),
        "mean_inter_r": mean_inter_r,
        "n_core_excl": int((df["d_core"] < 0).sum()),
        "n_abs_stronger": int((df["inter_d"] < 0).sum()),
        "n_fold_stronger": int((df["inter_r"] < 0).sum()),
        "wilcox_d": wilcox(inter_d),
        "wilcox_r": wilcox(inter_r),
        "mean_count_core_low": float(df["ml_core"].mean()),
        "mean_count_margin_low": float(df["ml_margin"].mean()),
        "aligned": aligned,
    }


def main():
    print("load", flush=True)
    bundle = base.load_obs_and_cldn4(base.H5)
    genes = extract_genes(base.H5, ["GZMB", "PRF1", "GNLY", "NKG7"])
    cats = bundle["cats"]["cell_type"]
    ct = bundle["obs"]["cell_type"]

    def ids(names):
        return base.codes_for(cats, names)

    immune = np.isin(ct, ids(IMMUNE))
    defs = {
        "cd8_nk": np.isin(ct, ids(("NK", "T CD8 memory", "T CD8 naive"))),
        "cd8": np.isin(ct, ids(("T CD8 memory", "T CD8 naive"))),
        "nk": np.isin(ct, ids(("NK",))),
        "lymphoid_tnk": np.isin(ct, ids(LYMPHOID)),
        "effector_rna": immune & ((genes["GZMB"] > 0) | (genes["PRF1"] > 0) | (genes["GNLY"] > 0)),
    }
    rich_ids = base.codes_for(bundle["cats"]["niche"], base.IMMUNE_RICH)
    sample_cats = bundle["cats"]["sample"]
    out_rows = []
    stored = {}
    for dname, mask in defs.items():
        print(f"definition {dname} n={int(mask.sum())}", flush=True)
        preps = {}
        for sname in base.SAMPLE_ORDER:
            prep = prepare(bundle, sample_cats.index(sname), rich_ids, mask)
            if prep is None:
                print(f"  skip {sname}", flush=True)
                continue
            preps[sname] = prep
            print(f"  {sname} tumor={prep['n_tumor']} cyto_cells={prep['n_cyto']}", flush=True)
        stored[dname] = preps
        for radius in RADII:
            for drop in (False, True):
                for mcut in MARGIN_CUTS:
                    for ccut in CORE_CUTS:
                        if ccut <= mcut:
                            continue
                        meta = {
                            "cyto": dname,
                            "radius_um": radius,
                            "cut": "absolute",
                            "margin_um": mcut,
                            "core_um": ccut,
                            "drop_lusc": drop,
                        }
                        rows = eval_sections(
                            preps, radius,
                            lambda d, c=ccut: d >= c,
                            lambda d, m=mcut: d <= m,
                            drop,
                        )
                        rec = summarize(rows, meta)
                        if rec:
                            out_rows.append(rec)
                for q_lo, q_hi in QUANTILES:
                    meta = {
                        "cyto": dname,
                        "radius_um": radius,
                        "cut": f"quantile_{q_lo:.2f}_{q_hi:.2f}",
                        "margin_um": np.nan,
                        "core_um": np.nan,
                        "drop_lusc": drop,
                    }

                    def core_q(dist, hi=q_hi):
                        return dist >= np.quantile(dist, hi)

                    def margin_q(dist, lo=q_lo):
                        return dist <= np.quantile(dist, lo)

                    rows = eval_sections(preps, radius, core_q, margin_q, drop)
                    rec = summarize(rows, meta)
                    if rec:
                        out_rows.append(rec)
    df = pd.DataFrame(out_rows)
    os.makedirs(base.TAB, exist_ok=True)
    path = os.path.join(base.TAB, "cutoff_grid.csv")
    df.to_csv(path, index=False)
    print("wrote", path, "rows", len(df), flush=True)
    aligned = df[df.aligned].copy()
    print("aligned specs", len(aligned), "of", len(df), flush=True)
    if aligned.empty:
        # closest: fold-aligned, smallest (least positive / most negative) absolute interaction
        fold = df[(df.mean_inter_r < 0) & (df.mean_d_core < 0)].copy()
        fold = fold.sort_values(["mean_inter_d", "frac_both"])
        print("NO spec with both mean gaps core-stronger. Closest absolute interactions:")
        print(fold.head(15).to_string(index=False))
        return
    aligned["score"] = (
        aligned["frac_both"]
        + 0.25 * (aligned["n_pat_both"] / aligned["n_pat"])
        + 0.05 * ((aligned["wilcox_d"] < 0.05) & (aligned["wilcox_r"] < 0.05)).astype(float)
        - 0.02 * aligned["drop_lusc"].astype(float)
    )
    aligned = aligned.sort_values(
        ["score", "n_both", "n_pat_both"], ascending=False
    )
    cols = [
        "cyto", "radius_um", "cut", "margin_um", "core_um", "drop_lusc",
        "n_sections", "n_both", "n_pat_both", "n_pat",
        "mean_d_core", "mean_d_margin", "mean_inter_d",
        "mean_r_core", "mean_r_margin", "mean_inter_r",
        "n_abs_stronger", "n_fold_stronger", "wilcox_d", "wilcox_r",
        "mean_count_core_low", "mean_count_margin_low", "score",
    ]
    print(aligned[cols].head(25).to_string(index=False))
    # keep preps only implicitly via rerun of the winner in a later step
    aligned.head(40)[cols].to_csv(os.path.join(base.TAB, "cutoff_grid_top.csv"), index=False)


if __name__ == "__main__":
    main()
