#!/usr/bin/env python3
"""One real Palantir run on the concordant-4 epithelial object.

Early cell is the same IFN-high, CLDN4-low root used for the PAGA headline.
Waypoints are 400 (package default is 1200) so the run finishes on this machine.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import stats

import palantir

from sweep_trajectory import AT2, BARRIER, DATA, TAB, aucell, hvg_idx, load_ifn, present, zmean

def main() -> None:
    blob = np.load(DATA, allow_pickle=True)
    X = blob["X"].astype(np.float32)
    genes = [str(g) for g in blob["genes"]]
    dataset = blob["dataset"].astype(str)
    unit = blob["unit"].astype(str)
    pool = blob["pool"].astype(str)
    uid = np.array([f"{d}:{u}" for d, u in zip(dataset, unit)])
    tumor = pool == "tumor"
    cldn4 = X[:, present(["CLDN4"], genes)[0]]
    ifn_cols = present(load_ifn(), genes)
    b_cols = present(list(BARRIER), genes)
    ifn_z = zmean(X, ifn_cols)
    cand = np.flatnonzero(cldn4 <= np.median(cldn4[tumor]))
    root = int(cand[np.argmax(ifn_z[cand])])
    cols = hvg_idx(X[tumor], 1000, [])
    ad = AnnData(X[:, cols].copy())
    ad.obs_names = [str(i) for i in range(ad.n_obs)]
    sc.pp.scale(ad, max_value=10)
    sc.tl.pca(ad, n_comps=20, svd_solver="arpack")
    import harmonypy
    ho = harmonypy.run_harmony(ad.obsm["X_pca"], pd.DataFrame({"dataset": dataset}), ["dataset"], max_iter_harmony=8, verbose=False)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] == 20:
        Z = Z.T
    ad.obsm["X_pca"] = Z
    print("diffusion", flush=True)
    palantir.utils.run_diffusion_maps(ad, n_components=10, knn=30)
    palantir.utils.determine_multiscale_space(ad)
    print("palantir", flush=True)
    palantir.core.run_palantir(ad, early_cell=str(root), knn=30, num_waypoints=400, n_jobs=4, seed=4)
    pt = ad.obs["palantir_pseudotime"].to_numpy().astype(float)
    barrier = aucell(X, b_cols)
    ifn = aucell(X, ifn_cols)
    rows = []
    for u in pd.unique(uid[tumor]):
        m = tumor & (uid == u) & np.isfinite(pt)
        if m.sum() < 12:
            continue
        rows.append((dataset[m][0], cldn4[m].mean(), pt[m].mean(), barrier[m].mean(), ifn[m].mean()))
    df = pd.DataFrame(rows, columns=["dataset", "CLDN4", "PT", "barrier", "IFN"])
    out = {}
    for col in ("CLDN4", "barrier", "IFN"):
        rho, p = stats.spearmanr(df["PT"], df[col])
        out[col] = {"n": int(len(df)), "rho": float(rho), "p": float(p)}
    out["root"] = uid[root]
    out["n_datasets"] = int(df.dataset.nunique())
    out["waypoints"] = 400
    print(json.dumps(out, indent=2))
    (TAB / "palantir_once.json").write_text(json.dumps(out, indent=2) + "\n")
    df.to_csv(TAB / "palantir_unit_means.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
