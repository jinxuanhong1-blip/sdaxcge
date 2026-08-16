"""CLAIM B6 neighborhood test on GSE287472 CosMx LUAD (1 patient, ~64k cells).

This is the only public CosMx LUAD matrix we found with CLDN4 + TACSTD2 +
T/B markers and cell coordinates. n=1 patient: treat as a single-sample
exploratory test, not a cohort.

Pre-specified:
  * QC: nCount_RNA ≥ 20, Area.um2 in [20, 2000].
  * Pixel scale from Area.um2 / Area (median sqrt).
  * Radii 25 / 50 / 100 µm (self excluded, ≥5 neighbors).
  * Index = log1p(CP10K CLDN4); outcome = mean T / B / TB score in the ball.
  * Also report CLDN4-high (top tertile) vs rest neighbor-TB contrast
    (rank-biserial).
  * Sensitivity: restrict index cells to PanCK-high (Mean.PanCK ≥ median).
"""
import os
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

import common

warnings.simplefilter("ignore")

DATA = "data/gse287472"
OUT = "results/claim_B6"
os.makedirs(f"{OUT}/tables", exist_ok=True)

RADII_UM = (25, 50, 100)


def load():
    expr = pd.read_csv(f"{DATA}/exprMat.csv.gz")
    meta = pd.read_csv(f"{DATA}/metadata.csv.gz")
    # CosMx exprMat: first cols are fov, cell_ID
    id_cols = [c for c in expr.columns if c.lower() in ("fov", "cell_id", "cell_id")]
    # standard: fov, cell_ID
    fov_col = [c for c in expr.columns if c.lower() == "fov"][0]
    cid_col = [c for c in expr.columns if c.lower() in ("cell_id", "cell_id")][0]
    genes = [c for c in expr.columns if c not in (fov_col, cid_col)]
    expr["_key"] = expr[fov_col].astype(str) + "_" + expr[cid_col].astype(str)
    meta["_key"] = meta["fov"].astype(str) + "_" + meta["cell_ID"].astype(str)
    meta = meta.drop_duplicates("_key").set_index("_key")
    expr = expr.drop_duplicates("_key").set_index("_key")
    shared = expr.index.intersection(meta.index)
    expr = expr.loc[shared]
    meta = meta.loc[shared]
    return expr[genes], meta


def log_cp10k(counts):
    tot = counts.sum(axis=1).replace(0, np.nan)
    return np.log1p(counts.div(tot, axis=0) * 1e4)


def score_from_df(logx, genes):
    used = [g for g in genes if g in logx.columns]
    if not used:
        return pd.Series(np.nan, index=logx.index), used
    # simple z-mean (scanpy-like without background)
    z = (logx[used] - logx[used].mean()) / logx[used].std(ddof=0).replace(0, np.nan)
    return z.mean(axis=1), used


def main():
    expr, meta = load()
    print(f"raw cells={len(meta)} genes={expr.shape[1]}")
    print("CLDN4 in panel:", "CLDN4" in expr.columns, "TACSTD2:", "TACSTD2" in expr.columns)

    qc = (
        (meta["nCount_RNA"] >= 20)
        & (meta["Area.um2"] >= 20)
        & (meta["Area.um2"] <= 2000)
        & np.isfinite(meta["CenterX_global_px"])
        & np.isfinite(meta["CenterY_global_px"])
    )
    expr, meta = expr.loc[qc], meta.loc[qc]
    print(f"after QC: {len(meta)} cells")

    # µm per pixel from Area.um2 / Area (px)
    scale = np.sqrt(meta["Area.um2"] / meta["Area"]).median()
    print(f"inferred µm/px = {scale:.4f}")

    xy = np.column_stack([
        meta["CenterX_global_px"].values * scale,
        meta["CenterY_global_px"].values * scale,
    ])

    logx = log_cp10k(expr)
    cldn4 = logx["CLDN4"] if "CLDN4" in logx.columns else pd.Series(np.nan, index=logx.index)
    t_sc, t_used = score_from_df(logx, common.T_MARKERS)
    b_sc, b_used = score_from_df(logx, common.B_MARKERS)
    tb_sc, tb_used = score_from_df(logx, common.T_MARKERS + common.B_MARKERS)
    print("T used", t_used)
    print("B used", b_used)

    tree = cKDTree(xy)
    rows = []
    # also PanCK-high subset
    panck_hi = meta["Mean.PanCK"] >= meta["Mean.PanCK"].median()

    for radius in RADII_UM:
        # query_ball is ok at 64k for 3 radii
        lists = tree.query_ball_point(xy, r=radius)
        t_nb = np.full(len(xy), np.nan)
        b_nb = np.full(len(xy), np.nan)
        tb_nb = np.full(len(xy), np.nan)
        n_nb = np.zeros(len(xy), int)
        t_vals = t_sc.values
        b_vals = b_sc.values
        tb_vals = tb_sc.values
        for i, nbrs in enumerate(lists):
            nbrs = [j for j in nbrs if j != i]
            n_nb[i] = len(nbrs)
            if len(nbrs) < 5:
                continue
            t_nb[i] = np.nanmean(t_vals[nbrs])
            b_nb[i] = np.nanmean(b_vals[nbrs])
            tb_nb[i] = np.nanmean(tb_vals[nbrs])

        def pack(tag, mask):
            idx = cldn4.values[mask]
            rec = dict(radius_um=radius, scope=tag, n_index=int(mask.sum()),
                       median_neighbors=float(np.median(n_nb[mask])))
            for lab, nb in [("T", t_nb), ("B", b_nb), ("TB", tb_nb)]:
                rho, p, n = common.spearman(idx, nb[mask])
                rec[f"rho_CLDN4_vs_{lab}nb"] = rho
                rec[f"p_CLDN4_vs_{lab}nb"] = p
                rec[f"n_{lab}"] = n
            # tertile contrast
            if np.isfinite(idx).sum() >= 50:
                thr = np.nanquantile(idx, 2 / 3)
                hi = tb_nb[mask][idx >= thr]
                lo = tb_nb[mask][idx < thr]
                rb, p = common.rank_biserial(hi, lo)
                rec["rb_TBnb_in_CLDN4high"] = rb
                rec["p_rb"] = p
            print(f"r={radius:3d}µm [{tag:16s}] n={rec['n_index']:6d}  "
                  f"ρ(CLDN4,TBnb)={rec['rho_CLDN4_vs_TBnb']:+.3f} p={rec['p_CLDN4_vs_TBnb']:.1e}  "
                  f"ρ(T)={rec['rho_CLDN4_vs_Tnb']:+.3f} ρ(B)={rec['rho_CLDN4_vs_Bnb']:+.3f}  "
                  f"rb={rec.get('rb_TBnb_in_CLDN4high', np.nan):+.3f}")
            return rec

        rows.append(pack("all_QC_cells", np.ones(len(xy), bool)))
        rows.append(pack("PanCK_high", panck_hi.values))

    pd.DataFrame(rows).to_csv(f"{OUT}/tables/gse287472_cosmx_neighborhood.csv", index=False)

    # same-cell (not neighborhood) as a reference
    same = []
    for tag, mask in [("all_QC_cells", np.ones(len(xy), bool)),
                      ("PanCK_high", panck_hi.values)]:
        rec = dict(scope=tag, n=int(mask.sum()))
        for lab, scv in [("T", t_sc.values), ("B", b_sc.values), ("TB", tb_sc.values)]:
            rho, p, n = common.spearman(cldn4.values[mask], scv[mask])
            rec[f"rho_CLDN4_vs_{lab}"] = rho
            rec[f"p_CLDN4_vs_{lab}"] = p
        same.append(rec)
        print(f"same-cell [{tag}] ρ(CLDN4,TB)={rec['rho_CLDN4_vs_TB']:+.3f}")
    pd.DataFrame(same).to_csv(f"{OUT}/tables/gse287472_cosmx_samecell.csv", index=False)
    print("wrote CosMx tables")


if __name__ == "__main__":
    main()
