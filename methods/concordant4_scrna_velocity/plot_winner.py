#!/usr/bin/env python3
"""Replay the selected full-cohort clocks and write figures."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from sweep_trajectory import (
    AT2, BARRIER, DATA, FIG, TAB, add_module, aucell, hvg_idx, load_ifn,
    present, spearman, ucell, zmean, cluster_pt,
)

FIG.mkdir(parents=True, exist_ok=True)


def build(X, dataset, tumor, n_top, n_pcs, n_nb, harmony, drop_cols):
    cols = hvg_idx(X[tumor], n_top, drop_cols)
    ad = AnnData(X[:, cols].copy())
    sc.pp.scale(ad, max_value=10)
    n_comps = min(n_pcs, cols.size - 1, ad.n_obs - 1)
    sc.tl.pca(ad, n_comps=n_comps, svd_solver="arpack")
    if harmony:
        import harmonypy
        meta = pd.DataFrame({"dataset": dataset})
        ho = harmonypy.run_harmony(ad.obsm["X_pca"], meta, ["dataset"], max_iter_harmony=8, verbose=False)
        Z = np.asarray(ho.Z_corr)
        if Z.shape[0] == ad.obsm["X_pca"].shape[1]:
            Z = Z.T
        ad.obsm["X_pca"] = Z
    sc.pp.neighbors(ad, n_neighbors=n_nb, n_pcs=n_comps, random_state=4)
    sc.tl.leiden(ad, resolution=0.5, key_added="leiden", flavor="igraph", random_state=4)
    sc.tl.paga(ad, groups="leiden")
    sc.tl.diffmap(ad)
    sc.tl.umap(ad, random_state=4)
    conn = np.asarray(ad.uns["paga"]["connectivities"].todense())
    clusters = ad.obs["leiden"].astype(int).to_numpy()
    return ad, conn, clusters


def unit_means(uid, tumor, dataset, cldn4, pt, barrier, ifn, tj, apm):
    rows = []
    for u in pd.unique(uid[tumor]):
        m = tumor & (uid == u) & np.isfinite(pt)
        if m.sum() < 12:
            continue
        rows.append({
            "unit": u,
            "dataset": dataset[m][0],
            "CLDN4": float(cldn4[m].mean()),
            "PT": float(pt[m].mean()),
            "barrier": float(barrier[m].mean()),
            "IFN": float(ifn[m].mean()),
            "TJ": float(tj[m].mean()),
            "APM": float(apm[m].mean()),
        })
    return pd.DataFrame(rows)


def scatter(df, path, title):
    colors = {
        "GSE123902": "#1b9e77",
        "GSE131907": "#d95f02",
        "GSE205335": "#7570b3",
        "GSE189357": "#e7298a",
    }
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.4))
    for ax, col, ylab in zip(axes, ("CLDN4", "barrier", "IFN"), ("CLDN4 (log norm)", "barrier AUCell", "IFN AUCell")):
        for ds, sub in df.groupby("dataset"):
            ax.scatter(sub["PT"], sub[col], s=28, c=colors.get(ds, "gray"), label=ds, alpha=0.9)
        n, rho, p = spearman(df["PT"].to_numpy(), df[col].to_numpy())
        ax.set_title(f"ρ={rho:.3f}  p={p:.2e}  n={n}")
        ax.set_xlabel("PAGA pseudotime")
        ax.set_ylabel(ylab)
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close()


def main():
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
    t_cols = present(["OCLN", "TJP1", "TJP2", "TJP3", "CLDN3", "CLDN7", "CDH1", "F11R", "MARVELD2", "CGN", "CRB3"], genes)
    a_cols = present(["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "NLRC5"], genes)
    auc = {"barrier": aucell(X, b_cols), "IFN": aucell(X, ifn_cols), "TJ": aucell(X, t_cols), "APM": aucell(X, a_cols)}
    uc = {"barrier": ucell(X, b_cols), "IFN": ucell(X, ifn_cols)}

    # Primary clock: PAGA, IFN-high root, the selected grid point.
    at2 = zmean(X, present(list(AT2), genes))
    ifn_z = zmean(X, ifn_cols)
    med_c = np.median(cldn4[tumor])
    cand = np.flatnonzero(cldn4 <= med_c)
    root_ifn = int(cand[np.argmax(ifn_z[cand])])
    root_at2 = int(np.flatnonzero(pool == "root")[np.argmax(at2[pool == "root"])])

    ad, conn, clusters = build(X, dataset, tumor, 1000, 20, 30, True, [])
    pt = cluster_pt(conn, clusters, root_ifn, "paga", cldn4)
    df = unit_means(uid, tumor, dataset, cldn4, pt, auc["barrier"], auc["IFN"], auc["TJ"], auc["APM"])
    df.to_csv(TAB / "winner_paga_ifnroot_unit_means.tsv", sep="\t", index=False)
    scatter(df, FIG / "fig_paga_ifnroot_patient.png",
            "Trajectory, not splicing velocity. PAGA from an IFN-high, CLDN4-low root. n=65 units.")

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
    xy = ad.obsm["X_umap"]
    axes[0].scatter(xy[:, 0], xy[:, 1], c=cldn4, s=3, cmap="viridis", linewidths=0)
    axes[0].set_title("CLDN4")
    axes[1].scatter(xy[:, 0], xy[:, 1], c=pt, s=3, cmap="plasma", linewidths=0)
    axes[1].scatter([xy[root_ifn, 0]], [xy[root_ifn, 1]], c="red", s=28, label="IFN-high root")
    axes[1].legend(fontsize=7, frameon=False)
    axes[1].set_title("PAGA pseudotime")
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Expression graph only. Not RNA velocity.", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "fig_umap_cldn4_pt.png", dpi=160)
    fig.savefig(FIG / "fig_umap_cldn4_pt.pdf")
    plt.close()

    # External AT2 root, DPT
    ad2, conn2, clusters2 = build(X, dataset, tumor, 2000, 40, 15, False, [])
    ad2.uns["iroot"] = root_at2
    sc.tl.dpt(ad2, n_dcs=min(10, ad2.obsm["X_diffmap"].shape[1]))
    pt2 = ad2.obs["dpt_pseudotime"].to_numpy().astype(float)
    df2 = unit_means(uid, tumor, dataset, cldn4, pt2, uc["barrier"], uc["IFN"], uc["barrier"], uc["IFN"])
    df2 = df2.rename(columns={"barrier": "barrier", "IFN": "IFN"})
    # unit_means used auc names; we passed ucell. relabel columns already barrier/IFN.
    df2.to_csv(TAB / "winner_dpt_at2root_unit_means.tsv", sep="\t", index=False)
    # custom scatter labels
    colors = {"GSE123902": "#1b9e77", "GSE131907": "#d95f02", "GSE205335": "#7570b3", "GSE189357": "#e7298a"}
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.4))
    for ax, col, ylab in zip(axes, ("CLDN4", "barrier", "IFN"), ("CLDN4 (log norm)", "barrier UCell", "IFN UCell")):
        for ds, sub in df2.groupby("dataset"):
            ax.scatter(sub["PT"], sub[col], s=28, c=colors[ds], label=ds, alpha=0.9)
        n, rho, p = spearman(df2["PT"].to_numpy(), df2[col].to_numpy())
        ax.set_title(f"ρ={rho:.3f}  p={p:.2e}  n={n}")
        ax.set_xlabel("DPT from nLung AT2")
        ax.set_ylabel(ylab)
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("External root: uninvolved-lung AT2. Trajectory, not splicing velocity.", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "fig_dpt_at2root_patient.png", dpi=160)
    fig.savefig(FIG / "fig_dpt_at2root_patient.pdf")
    plt.close()
    print("figures written", root_ifn, uid[root_ifn], root_at2, uid[root_at2])
    print("paga n", len(df), "dpt n", len(df2))


if __name__ == "__main__":
    main()
