#!/usr/bin/env python3
"""Palantir + PAGA + Slingshot on concordant-4 malignant cells.

Cell-state path: AT2-like malignant → barrier-like malignant.
CLDN4 is a gene tracked along that path. It does not define the root or the
barrier terminal.

Not an ICI clock. Pseudotime is not treatment time. RECIST/MPR are not tested.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.sparse.csgraph import shortest_path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import AT2, BARRIER, CONTROL, FOCAL, IFN  # noqa: E402

SEED = 1
N_HVG = 2000
N_PCS = 30
N_NEIGH = 30
LEIDEN_RES = 0.6
MIN_CLUSTER = 40
MIN_UNIT_CELLS = 10


def _score(adata, name: str, genes: tuple[str, ...]) -> list[str]:
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    if len(present) < 2:
        raise SystemExit(f"{name} has <2 genes present: missing {missing}")
    sc.tl.score_genes(adata, present, score_name=f"score_{name}", use_raw=False)
    adata.uns[f"genes_{name}"] = present
    adata.uns[f"genes_{name}_missing"] = missing
    return missing


def _preprocess(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=N_HVG, batch_key="dataset", flavor="seurat"
    )
    sc.tl.pca(adata, n_comps=N_PCS, use_highly_variable=True, svd_solver="arpack", random_state=SEED)
    import harmonypy

    ho = harmonypy.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        "dataset",
        max_iter_harmony=30,
        random_state=SEED,
    )
    Z = np.asarray(ho.Z_corr)
    if Z.shape == (adata.n_obs, N_PCS):
        adata.obsm["X_pca_harmony"] = Z
    elif Z.shape == (N_PCS, adata.n_obs):
        adata.obsm["X_pca_harmony"] = Z.T
    else:
        raise SystemExit(f"unexpected harmony shape {Z.shape}")
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=N_NEIGH, random_state=SEED)
    sc.tl.umap(adata, min_dist=0.3, random_state=SEED)
    sc.tl.leiden(
        adata,
        resolution=LEIDEN_RES,
        key_added="leiden",
        random_state=SEED,
        flavor="igraph",
        n_iterations=2,
    )
    sc.tl.paga(adata, groups="leiden")


def _cluster_means(adata) -> pd.DataFrame:
    rows = []
    for lab, idx in adata.obs.groupby("leiden", observed=True).groups.items():
        sub = adata.obs.loc[idx]
        rows.append(
            {
                "leiden": str(lab),
                "n": int(len(sub)),
                "AT2": float(sub["score_AT2"].mean()),
                "barrier": float(sub["score_barrier"].mean()),
                "IFN": float(sub["score_IFN"].mean()),
                "CLDN4": float(sub["expr_CLDN4"].mean()),
                "SFTPB": float(sub["expr_SFTPB"].mean()) if "expr_SFTPB" in sub else np.nan,
            }
        )
    return pd.DataFrame(rows).set_index("leiden")


def _nearest(obs: pd.DataFrame, col: str, target: float) -> str:
    d = (obs[col] - target).abs()
    return str(d.idxmin())


def _choose_poles(adata) -> dict:
    means = _cluster_means(adata)
    eligible = means.loc[means["n"] >= MIN_CLUSTER].copy()
    if eligible.empty:
        raise SystemExit("no Leiden cluster passed the size floor")
    cldn4_max = str(eligible["CLDN4"].idxmax())
    at2_order = eligible.sort_values("AT2", ascending=False).index.tolist()
    root = next(c for c in at2_order if c != cldn4_max)
    barrier = str(eligible.drop(index=root)["barrier"].idxmax())
    rest = eligible.drop(index=[root, barrier], errors="ignore")
    if rest.empty:
        raise SystemExit("need a second terminal cluster besides root and barrier")
    alt = str(rest["IFN"].idxmax())
    root_cells = adata.obs.loc[adata.obs["leiden"].astype(str) == root]
    med_cl = float(root_cells["expr_CLDN4"].median())
    low = root_cells.loc[root_cells["expr_CLDN4"] <= med_cl]
    if low.empty:
        low = root_cells
    early = _nearest(low, "score_AT2", float(low["score_AT2"].median()))
    bar_cells = adata.obs.loc[adata.obs["leiden"].astype(str) == barrier]
    q80 = float(bar_cells["score_barrier"].quantile(0.80))
    top = bar_cells.loc[bar_cells["score_barrier"] >= q80]
    barrier_cell = _nearest(top if len(top) else bar_cells, "score_barrier", q80)
    alt_cells = adata.obs.loc[adata.obs["leiden"].astype(str) == alt]
    alt_cell = _nearest(alt_cells, "score_IFN", float(alt_cells["score_IFN"].median()))
    info = {
        "root_cluster": root,
        "barrier_cluster": barrier,
        "alt_cluster": alt,
        "cldn4_max_cluster_excluded_from_root": cldn4_max,
        "early_cell": early,
        "barrier_cell": barrier_cell,
        "alt_cell": alt_cell,
        "early_CLDN4": float(adata.obs.loc[early, "expr_CLDN4"]),
        "early_AT2": float(adata.obs.loc[early, "score_AT2"]),
        "early_CLDN4_tertile_rule": "at or below root-cluster median CLDN4",
        "barrier_cell_rule": "nearest 80th percentile of barrier score; CLDN4 not used",
        "alt_rule": "highest mean IFN among clusters that are neither root nor barrier",
        "cluster_means": means.reset_index().to_dict(orient="records"),
    }
    return info


def _paga_path(adata, root: str, barrier: str) -> dict:
    cats = [str(c) for c in adata.obs["leiden"].cat.categories]
    conn = adata.uns["paga"]["connectivities"]
    arr = conn.toarray() if hasattr(conn, "toarray") else np.asarray(conn)
    # distance = inverse connectivity; zero edges are absent
    dist = np.full_like(arr, np.inf, dtype=float)
    pos = arr > 0
    dist[pos] = 1.0 / arr[pos]
    np.fill_diagonal(dist, 0.0)
    i = cats.index(root)
    j = cats.index(barrier)
    sp = shortest_path(dist, directed=False, indices=i)
    # reconstruct path
    if not np.isfinite(sp[j]):
        return {"connected": False, "path": [], "weight": 0.0}
    path = [j]
    cur = j
    # greedy predecessor
    while cur != i:
        nbrs = np.where(np.isfinite(dist[cur]) & (dist[cur] < np.inf) & (np.arange(len(cats)) != cur))[0]
        if len(nbrs) == 0:
            break
        nxt = nbrs[np.argmin(sp[nbrs] + 1e-9 * dist[cur, nbrs])]
        # only step if it decreases distance-to-root
        better = [n for n in nbrs if sp[n] < sp[cur] - 1e-8]
        if not better:
            break
        nxt = min(better, key=lambda n: sp[n])
        path.append(int(nxt))
        cur = int(nxt)
        if len(path) > len(cats):
            break
    path = path[::-1]
    nodes = [cats[k] for k in path]
    return {
        "connected": nodes[0] == root and nodes[-1] == barrier,
        "path": nodes,
        "n_edges": max(0, len(nodes) - 1),
    }


def _run_palantir(adata, poles: dict) -> dict:
    import palantir

    palantir.utils.run_diffusion_maps(
        adata,
        n_components=10,
        knn=N_NEIGH,
        pca_key="X_pca_harmony",
        seed=SEED,
    )
    palantir.utils.determine_multiscale_space(adata)
    n_wp = int(min(400, max(150, adata.n_obs // 20)))
    # Palantir 1.4 accepts {cell_id: label} (docstring order is inverted).
    terminals = {
        poles["barrier_cell"]: "barrier",
        poles["alt_cell"]: "alt",
    }
    pr = palantir.core.run_palantir(
        adata,
        early_cell=poles["early_cell"],
        terminal_states=terminals,
        num_waypoints=n_wp,
        knn=N_NEIGH,
        n_jobs=2,
        seed=SEED,
        use_early_cell_as_start=True,
    )
    branch = None
    if pr is not None and getattr(pr, "branch_probs", None) is not None:
        branch = pr.branch_probs
    if branch is None and "palantir_fate_probabilities" in adata.obsm:
        raw = adata.obsm["palantir_fate_probabilities"]
        branch = raw if isinstance(raw, pd.DataFrame) else pd.DataFrame(raw, index=adata.obs_names)
    if not isinstance(branch, pd.DataFrame):
        raise RuntimeError(f"Palantir fate table type {type(branch)}")
    branch = branch.reindex(adata.obs_names)
    # Map columns onto barrier / alt. Columns may be labels or cell ids.
    colmap = {}
    for col in branch.columns:
        s = str(col)
        if s == "barrier" or s == poles["barrier_cell"]:
            colmap[col] = "barrier"
        elif s == "alt" or s == poles["alt_cell"]:
            colmap[col] = "alt"
    if "barrier" not in colmap.values():
        # last resort: column whose terminal cell sits in the barrier cluster
        raise RuntimeError(f"barrier fate column missing from {list(branch.columns)}")
    branch = branch.rename(columns=colmap)
    adata.obs["palantir_pseudotime"] = pd.to_numeric(
        adata.obs["palantir_pseudotime"] if "palantir_pseudotime" in adata.obs
        else (pr.pseudotime if pr is not None else np.nan),
        errors="coerce",
    )
    if "palantir_pseudotime" not in adata.obs or adata.obs["palantir_pseudotime"].isna().all():
        pt = getattr(pr, "pseudotime", None)
        adata.obs["palantir_pseudotime"] = pd.to_numeric(pd.Series(pt, index=adata.obs_names), errors="coerce")
    adata.obs["fate_barrier"] = pd.to_numeric(branch["barrier"], errors="coerce")
    if "alt" in branch.columns:
        adata.obs["fate_alt"] = pd.to_numeric(branch["alt"], errors="coerce")
    else:
        adata.obs["fate_alt"] = 1.0 - adata.obs["fate_barrier"]
    n_fin = int(np.isfinite(adata.obs["palantir_pseudotime"]).sum())
    if n_fin < adata.n_obs * 0.5:
        raise RuntimeError(f"Palantir pseudotime finite for only {n_fin} cells")
    return {
        "version": getattr(palantir, "__version__", "unknown"),
        "n_waypoints": n_wp,
        "n_finite_pt": n_fin,
        "terminals": {"barrier": poles["barrier_cell"], "alt": poles["alt_cell"]},
        "engine": "palantir",
    }


def _load_slingshot_module():
    """Load Street et al. 2018 port without pyslingshot/__init__.py (needs ggplot2_py)."""
    import importlib.util
    import site
    import sys
    import sysconfig

    candidates = [
        Path(sysconfig.get_paths()["purelib"]) / "pyslingshot" / "core.py",
        Path(site.getusersitepackages()) / "pyslingshot" / "core.py",
    ]
    for sp in site.getsitepackages():
        candidates.append(Path(sp) / "pyslingshot" / "core.py")
    path = next((c for c in candidates if c.is_file()), None)
    if path is None:
        raise RuntimeError(f"pyslingshot/core.py not found in {candidates}")
    spec = importlib.util.spec_from_file_location("pyslingshot_core_street", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load slingshot core at {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _run_slingshot(adata, poles: dict) -> dict:
    mod = _load_slingshot_module()
    X = np.asarray(adata.obsm["X_pca_harmony"][:, :10], dtype=float)
    labels = adata.obs["leiden"].astype(str).to_numpy()
    # getLineages is the MST. Curves are then fit on lineage cells only:
    # the port's getCurves smooths every cell and does not apply its weights.
    sr = mod.getLineages(X, labels, start_cluster=poles["root_cluster"])
    n = X.shape[0]
    lineages_raw = sr.lineages
    L = len(lineages_raw)
    weights = np.zeros((n, L), dtype=float)
    pseudotime = np.full((n, L), np.nan)
    for li, lineage in enumerate(lineages_raw):
        in_lin = np.isin(sr.cluster_labels, list(lineage))
        weights[in_lin, li] = 1.0
        if int(in_lin.sum()) < 15:
            continue
        init = np.stack([sr.cluster_centers[c] for c in lineage])
        _curve, lam, _ord = mod._principal_curve(
            X[in_lin], init, max_iter=4, smoother_span=0.4
        )
        pseudotime[in_lin, li] = lam
    sr.weights = weights
    sr.pseudotime = pseudotime
    lineages = [[str(x) for x in lin] for lin in sr.lineages]
    pt = np.asarray(sr.pseudotime, dtype=float)
    wt = np.asarray(sr.weights, dtype=float)
    if pt.ndim == 1:
        pt = pt[:, None]
    if wt.ndim == 1:
        wt = wt[:, None]
    barrier = poles["barrier_cluster"]
    hits = [i for i, lin in enumerate(lineages) if barrier in lin]
    if not hits:
        # lineage whose leaf has the highest barrier score
        means = (
            adata.obs.groupby(adata.obs["leiden"].astype(str), observed=True)["score_barrier"].mean()
        )
        def leaf_score(lin):
            return float(means.get(lin[-1], -np.inf))
        best = int(np.argmax([leaf_score(lin) for lin in lineages]))
        chosen_rule = "no lineage contains the barrier cluster; used highest-barrier leaf"
    else:
        # prefer the lineage where barrier is a leaf, else the shortest path that includes it
        leaves = [i for i in hits if lineages[i][-1] == barrier]
        if leaves:
            best = leaves[0]
            chosen_rule = "lineage leaf is the barrier cluster"
        else:
            best = min(hits, key=lambda i: lineages[i].index(barrier))
            chosen_rule = "lineage passes through the barrier cluster"
    adata.obs["slingshot_pt_barrier"] = pt[:, best]
    adata.obs["slingshot_w_barrier"] = wt[:, best]
    # shared pseudotime = weight-averaged finite times
    with np.errstate(invalid="ignore"):
        wsum = np.nansum(np.where(np.isfinite(pt), wt, 0.0), axis=1)
        ssum = np.nansum(np.where(np.isfinite(pt), wt * pt, 0.0), axis=1)
        shared = np.divide(ssum, wsum, out=np.full(len(ssum), np.nan), where=wsum > 0)
    adata.obs["slingshot_pt_shared"] = shared
    # progress toward barrier: pseudotime clipped at the barrier-cluster median
    on_prefix = adata.obs["leiden"].astype(str).isin(lineages[best][: lineages[best].index(barrier) + 1])
    adata.obs["on_at2_barrier_prefix"] = on_prefix.to_numpy()
    return {
        "engine": "pyslingshot.core getLineages (Street et al. 2018 MST) + lineage-restricted principal curves",
        "n_lineages": len(lineages),
        "lineages": lineages,
        "barrier_lineage_index": int(best),
        "barrier_lineage": lineages[best],
        "chosen_rule": chosen_rule,
        "barrier_is_leaf": lineages[best][-1] == barrier,
        "start_cluster": poles["root_cluster"],
    }


def _spearman(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


def _dl_meta(pairs: list[tuple[float, int]]) -> dict:
    zs, vs = [], []
    for rho, n in pairs:
        if n < 5 or not np.isfinite(rho) or abs(rho) >= 1:
            if n >= 5 and np.isfinite(rho):
                rho = float(np.clip(rho, -0.999999, 0.999999))
            else:
                continue
        zs.append(np.arctanh(rho))
        vs.append(1.0 / (n - 3))
    if len(zs) < 2:
        return {"k": len(zs), "rho": np.nan, "p": np.nan, "I2": np.nan, "tau2": np.nan}
    zs = np.asarray(zs)
    vs = np.asarray(vs)
    w = 1.0 / vs
    zbar = np.sum(w * zs) / np.sum(w)
    Q = float(np.sum(w * (zs - zbar) ** 2))
    k = len(zs)
    df = k - 1
    C = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    wstar = 1.0 / (vs + tau2)
    zmeta = float(np.sum(wstar * zs) / np.sum(wstar))
    se = float(np.sqrt(1.0 / np.sum(wstar)))
    p = float(2 * stats.norm.sf(abs(zmeta / se))) if se > 0 else np.nan
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    return {"k": k, "rho": float(np.tanh(zmeta)), "p": p, "I2": float(I2), "tau2": float(tau2), "Q": Q}


def _patient_table(adata, locked: pd.DataFrame) -> pd.DataFrame:
    g = adata.obs.groupby(["dataset", "unit_id"], observed=True)
    rows = []
    for (ds, unit), sub in g:
        rows.append(
            {
                "dataset": ds,
                "unit_id": str(unit),
                "n_cells": int(len(sub)),
                "mean_fate_barrier": float(sub["fate_barrier"].mean()),
                "mean_palantir_pt": float(sub["palantir_pseudotime"].mean()),
                "mean_slingshot_pt": float(sub["slingshot_pt_barrier"].mean()),
                "mean_slingshot_w": float(sub["slingshot_w_barrier"].mean()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier": float(sub["score_barrier"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "frac_on_prefix": float(sub["on_at2_barrier_prefix"].mean()),
            }
        )
    pt = pd.DataFrame(rows)
    key = locked.copy()
    key["unit_id"] = key["unit_id"].astype(str)
    key["dataset"] = key["dataset"].astype(str)
    keep = key[["dataset", "unit_id", "frac_tnk", "n_malignant", "n_tnk", "mal_CLDN4_pct"]].copy()
    out = pt.merge(keep, on=["dataset", "unit_id"], how="left", validate="one_to_one")
    out["eligible"] = (out["n_cells"] >= MIN_UNIT_CELLS) & out["frac_tnk"].notna()
    return out


def _bin_curve(df: pd.DataFrame, pt_col: str, n_bins: int = 8) -> pd.DataFrame:
    x = df[pt_col].to_numpy(dtype=float)
    m = np.isfinite(x)
    rows = []
    if m.sum() < n_bins * 5:
        return pd.DataFrame(rows)
    qs = np.quantile(x[m], np.linspace(0, 1, n_bins + 1))
    qs[0] -= 1e-9
    qs[-1] += 1e-9
    for i in range(n_bins):
        sel = m & (x > qs[i]) & (x <= qs[i + 1])
        if sel.sum() == 0:
            continue
        sub = df.loc[sel]
        rows.append(
            {
                "bin": i,
                "n_cells": int(sel.sum()),
                "pt_mean": float(sub[pt_col].mean()),
                "CLDN4": float(sub["expr_CLDN4"].mean()),
                "AT2": float(sub["score_AT2"].mean()),
                "barrier": float(sub["score_barrier"].mean()),
                "IFN": float(sub["score_IFN"].mean()),
            }
        )
    return pd.DataFrame(rows)


def _fmt(x, nd=3):
    if x is None or not np.isfinite(x):
        return "NA"
    if abs(x) >= 0.001 or x == 0:
        return f"{x:.{nd}f}"
    return f"{x:.2e}"


def _fmt_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3f}"


def _tests(patients: pd.DataFrame) -> pd.DataFrame:
    elig = patients.loc[patients["eligible"]].copy()
    specs = [
        ("fate_barrier vs T/NK", "mean_fate_barrier", "frac_tnk"),
        ("slingshot barrier PT vs T/NK", "mean_slingshot_pt", "frac_tnk"),
        ("slingshot barrier weight vs T/NK", "mean_slingshot_w", "frac_tnk"),
        ("CLDN4 vs Palantir PT", "mean_CLDN4", "mean_palantir_pt"),
        ("CLDN4 vs Slingshot barrier PT", "mean_CLDN4", "mean_slingshot_pt"),
        ("barrier score vs Palantir PT", "mean_barrier", "mean_palantir_pt"),
        ("AT2 score vs Palantir PT", "mean_AT2", "mean_palantir_pt"),
        ("SFTPC-program AT2 vs fate_barrier", "mean_AT2", "mean_fate_barrier"),
        ("barrier score vs fate_barrier", "mean_barrier", "mean_fate_barrier"),
        ("CLDN4 vs fate_barrier", "mean_CLDN4", "mean_fate_barrier"),
    ]
    rows = []
    for name, a, b in specs:
        sp = _spearman(elig[a], elig[b])
        per = []
        for ds, sub in elig.groupby("dataset"):
            s2 = _spearman(sub[a], sub[b])
            per.append((ds, s2))
        meta = _dl_meta([(s2["rho"], s2["n"]) for _, s2 in per])
        row = {
            "contrast": name,
            "n": sp["n"],
            "rho_pooled": sp["rho"],
            "p_pooled": sp["p"],
            "k_cohorts": meta["k"],
            "rho_DL": meta["rho"],
            "p_DL": meta["p"],
            "I2": meta["I2"],
        }
        for ds, s2 in per:
            row[f"n_{ds}"] = s2["n"]
            row[f"rho_{ds}"] = s2["rho"]
            row[f"p_{ds}"] = s2["p"]
        rows.append(row)
    return pd.DataFrame(rows)


def _style():
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _scatter_umap(ax, adata, color, title, cmap="viridis"):
    xy = adata.obsm["X_umap"]
    vals = adata.obs[color]
    if pd.api.types.is_numeric_dtype(vals):
        sca = ax.scatter(xy[:, 0], xy[:, 1], c=vals, s=2, cmap=cmap, linewidths=0, rasterized=True)
        plt.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
    else:
        cats = pd.Categorical(vals)
        cmap_d = plt.get_cmap("tab20")
        for i, cat in enumerate(cats.categories):
            m = cats == cat
            ax.scatter(xy[m, 0], xy[m, 1], s=2, color=cmap_d(i % 20), linewidths=0, rasterized=True, label=str(cat))
        ax.legend(markerscale=4, fontsize=6, frameon=False, loc="best")
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])


def _figures(adata, patients, bins_pal, bins_sl, paga_info, outdir: Path):
    _style()
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(11.2, 7.2))
    _scatter_umap(axes[0, 0], adata, "dataset", "Dataset")
    _scatter_umap(axes[0, 1], adata, "score_AT2", "AT2 score", "YlGnBu")
    _scatter_umap(axes[0, 2], adata, "score_barrier", "Barrier score (no CLDN4)", "YlOrBr")
    _scatter_umap(axes[1, 0], adata, "expr_CLDN4", "CLDN4", "magma")
    _scatter_umap(axes[1, 1], adata, "palantir_pseudotime", "Palantir pseudotime", "viridis")
    _scatter_umap(axes[1, 2], adata, "fate_barrier", "Palantir P(barrier)", "cividis")
    fig.suptitle("Concordant-4 malignant cells — not an ICI clock", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "fig_umap_overview.png", dpi=200)
    fig.savefig(figdir / "fig_umap_overview.pdf")
    plt.close(fig)

    # PAGA
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    try:
        sc.pl.paga(
            adata,
            color="expr_CLDN4",
            ax=ax,
            show=False,
            title="PAGA colored by CLDN4",
            frameon=False,
        )
    except Exception:
        ax.text(0.1, 0.5, "PAGA draw failed; connectivities are in the table")
    fig.tight_layout()
    fig.savefig(figdir / "fig_paga_cldn4.png", dpi=200)
    fig.savefig(figdir / "fig_paga_cldn4.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), sharey=True)
    for ax, bins, title in (
        (axes[0], bins_pal, "Along Palantir pseudotime"),
        (axes[1], bins_sl, "Along Slingshot barrier lineage"),
    ):
        if bins.empty:
            ax.set_title(title + " (empty)")
            continue
        ax.plot(bins["bin"], bins["CLDN4"], marker="o", label="CLDN4")
        ax.plot(bins["bin"], bins["barrier"], marker="o", label="Barrier (no CLDN4)")
        ax.plot(bins["bin"], bins["AT2"], marker="o", label="AT2")
        ax.set_xlabel("Bin (0 = early)")
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=7)
    axes[0].set_ylabel("Mean")
    fig.suptitle("AT2 → barrier path. CLDN4 is the readout, not the terminal definition.")
    fig.tight_layout()
    fig.savefig(figdir / "fig_along_path.png", dpi=200)
    fig.savefig(figdir / "fig_along_path.pdf")
    plt.close(fig)

    elig = patients.loc[patients["eligible"]]
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    datasets = list(elig["dataset"].unique())
    cmap = plt.get_cmap("tab10")
    for i, ds in enumerate(datasets):
        sub = elig.loc[elig["dataset"] == ds]
        ax.scatter(sub["mean_fate_barrier"], sub["frac_tnk"], s=28, color=cmap(i), label=f"{ds} (n={len(sub)})")
    sp = _spearman(elig["mean_fate_barrier"], elig["frac_tnk"])
    ax.set_xlabel("Patient mean Palantir P(barrier terminal)")
    ax.set_ylabel("Locked T/NK fraction")
    ax.set_title(f"Barrier fate vs T/NK  ρ={_fmt(sp['rho'])}  p={_fmt_p(sp['p'])}  n={sp['n']}")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(figdir / "fig_fate_vs_tnk.png", dpi=200)
    fig.savefig(figdir / "fig_fate_vs_tnk.pdf")
    plt.close(fig)

    # honest n
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    counts = (
        patients.groupby("dataset")
        .apply(lambda d: pd.Series({"units": d["unit_id"].nunique(), "eligible": int(d["eligible"].sum()), "cells": int(d["n_cells"].sum())}), include_groups=False)
        .reset_index()
    )
    x = np.arange(len(counts))
    ax.bar(x - 0.15, counts["units"], width=0.3, label="Units in object")
    ax.bar(x + 0.15, counts["eligible"], width=0.3, label=f"Eligible (≥{MIN_UNIT_CELLS} cells + locked T/NK)")
    ax.set_xticks(x)
    ax.set_xticklabels(counts["dataset"], rotation=15)
    ax.set_ylabel("Patients / samples")
    ax.set_title("Honest n is units, not cells")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(figdir / "fig_honest_n.png", dpi=200)
    fig.savefig(figdir / "fig_honest_n.pdf")
    plt.close(fig)
    return str(paga_info.get("path"))


def _write_finding(path: Path, summary: dict, tests: pd.DataFrame, bins_pal: pd.DataFrame, bins_sl: pd.DataFrame):
    tlookup = tests.set_index("contrast")

    def row(name):
        r = tlookup.loc[name]
        return (
            f"n={int(r['n'])}, pooled ρ={_fmt(r['rho_pooled'])} (p={_fmt_p(r['p_pooled'])}); "
            f"DL ρ={_fmt(r['rho_DL'])} (p={_fmt_p(r['p_DL'])}, I²={_fmt(r['I2'], 2)}, k={int(r['k_cohorts'])})"
        )

    primary = tlookup.loc["fate_barrier vs T/NK"]
    pal_bins = bins_pal.to_dict(orient="records") if not bins_pal.empty else []
    sl_bins = bins_sl.to_dict(orient="records") if not bins_sl.empty else []

    def ends(bins):
        if len(bins) < 2:
            return "NA"
        a, b = bins[0], bins[-1]
        return (
            f"CLDN4 {_fmt(a['CLDN4'])} → {_fmt(b['CLDN4'])}; "
            f"barrier {_fmt(a['barrier'])} → {_fmt(b['barrier'])}; "
            f"AT2 {_fmt(a['AT2'])} → {_fmt(b['AT2'])}"
        )

    lines = []
    lines.append("# Finding — concordant-4 malignant Palantir + PAGA + Slingshot")
    lines.append("")
    lines.append("**This is not an ICI clock.** Pseudotime orders malignant cell state from an AT2-like pole to a barrier-like pole. It is not immunotherapy exposure, not time-on-treatment, and not MPR/RECIST. GSE205335 is an ICI cohort in the locked four, but RECIST was not a variable in any test here.")
    lines.append("")
    lines.append("ADDITIVE. **CLDN4 only.** Locked concordant-4 is GSE123902 + GSE131907 + GSE205335 + GSE189357. The locked malignant CLDN4 %pos vs T/NK result (n=65, ρ=−0.531, p=1.6×10⁻⁵, I²=0%) is not re-derived and is not re-audited. GSE148071, GSE127465, GSE154826, GSE200563, and E-MTAB-13526 are not added.")
    lines.append("")
    lines.append("Object: **malignant cells only**, same gates as the locked table. GSE131907 uses `Cell_subtype == Malignant cells` (tS1/tS2/tS3 and nLung AT2 are outside that gate and were not added). GSE205335 uses author `Malignant cells` on non-normal tissue. GSE123902 and GSE189357 use marker malignant `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`. Normal GSE123902 libraries were dropped.")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(
        f"Primary test is the patient/sample mean **Palantir probability of the barrier terminal** versus the **locked T/NK fraction**. {row('fate_barrier vs T/NK')}."
    )
    lines.append("")
    lines.append(
        f"CLDN4 along Palantir pseudotime (patient means): {row('CLDN4 vs Palantir PT')}. "
        f"CLDN4 along the Slingshot barrier lineage: {row('CLDN4 vs Slingshot barrier PT')}. "
        f"AT2 score vs Palantir PT (direction control): {row('AT2 score vs Palantir PT')}."
    )
    lines.append("")
    lines.append(f"Cell-bin ends, Palantir: {ends(pal_bins)}. Slingshot barrier lineage: {ends(sl_bins)}.")
    lines.append("")
    lines.append("## What was run")
    lines.append("")
    lines.append("- **PAGA** (`scanpy.tl.paga` on Leiden 0.6 after Harmony `batch=dataset`).")
    lines.append(
        f"- **Palantir** {summary['palantir']['version']}: diffusion maps on Harmony PCs, "
        f"{summary['palantir']['n_waypoints']} waypoints, early cell in the AT2-like cluster, "
        f"two terminals (barrier program; IFN-high alternative). Fate is an absorption probability, not a treatment clock."
    )
    lines.append(
        f"- **Slingshot** via `{summary['slingshot']['engine']}`. "
        f"Start cluster = AT2-like Leiden {summary['poles']['root_cluster']}. "
        f"Barrier lineage rule: {summary['slingshot']['chosen_rule']}. "
        f"Lineage: {' → '.join(summary['slingshot']['barrier_lineage'])}."
    )
    lines.append("")
    lines.append("The barrier terminal is the cell nearest the 80th percentile of the barrier/keratin score inside the highest-barrier Leiden cluster. **CLDN4 is not in that score and is not used to pick the cell.** The root cell is inside the highest-AT2 Leiden cluster after excluding the max-CLDN4 cluster, and is at or below that cluster's median CLDN4.")
    lines.append("")
    lines.append("## Honest n")
    lines.append("")
    lines.append(
        f"- Cells after QC and cap ≤{summary['cap']}/unit: **{summary['n_cells']}**."
    )
    lines.append(
        f"- Units in the object: **{summary['n_units']}**. Eligible for Spearman (n_cells≥{MIN_UNIT_CELLS} and locked T/NK present): **{summary['n_eligible']}**."
    )
    lines.append(f"- By dataset (cells / units / eligible): {summary['by_dataset']}.")
    lines.append(f"- PAGA path root → barrier: {summary['paga_path']}. Connected={summary['paga_connected']}.")
    lines.append(
        f"- Early cell `{summary['poles']['early_cell']}` CLDN4={_fmt(summary['poles']['early_CLDN4'])}, AT2={_fmt(summary['poles']['early_AT2'])}."
    )
    lines.append(f"- Slingshot lineages: {summary['slingshot']['n_lineages']}. Barrier cluster is a lineage leaf: {summary['slingshot']['barrier_is_leaf']}.")
    lines.append("")
    lines.append("## Primary and pre-specified tests")
    lines.append("")
    lines.append("Unit = patient (GSE123902 donor, GSE189357, GSE205335) or sample (GSE131907), matching the locked table. Pooled Spearman treats units equally. DL is Fisher-z DerSimonian–Laird across the four cohorts. Cell-level p-values are not the claim.")
    lines.append("")
    lines.append("| Contrast | n | pooled ρ | pooled p | DL ρ | DL p | I² |")
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for _, r in tests.iterrows():
        lines.append(
            f"| {r['contrast']} | {int(r['n'])} | {_fmt(r['rho_pooled'])} | {_fmt_p(r['p_pooled'])} | "
            f"{_fmt(r['rho_DL'])} | {_fmt_p(r['p_DL'])} | {_fmt(r['I2'], 2)} |"
        )
    lines.append("")
    lines.append("### Within cohort — barrier fate vs locked T/NK")
    lines.append("")
    lines.append("| Dataset | n | ρ | p |")
    lines.append("| --- | ---: | ---: | ---: |")
    for ds in ("GSE123902", "GSE131907", "GSE205335", "GSE189357"):
        ncol, rcol, pcol = f"n_{ds}", f"rho_{ds}", f"p_{ds}"
        if ncol not in primary:
            continue
        lines.append(f"| {ds} | {int(primary[ncol])} | {_fmt(primary[rcol])} | {_fmt_p(primary[pcol])} |")
    lines.append("")
    lines.append("## What this does not claim")
    lines.append("")
    lines.append("- This is not an ICI clock, not a response clock, and not a test of RECIST or MPR.")
    lines.append("- The locked CLDN4 %pos vs T/NK ρ=−0.531 was not recomputed. Mean CLDN4 on this capped object is a different score.")
    lines.append("- Palantir terminals and Slingshot lineages are phenotypic orderings, not demonstrated developmental lineages.")
    lines.append("- Do not write “AT2 differentiates into a barrier tumor because PAGA is connected.”")
    lines.append("- Malignant labels are author calls (GSE131907, GSE205335) or the locked marker gate. They are not CNV re-calls.")
    lines.append("- No TACSTD2∩CLDN4 dual-high gate. No GSE148071.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/concordant4_malignant_traj_cldn4/scripts/download.py --out /tmp/c4raw")
    lines.append("python3 methods/concordant4_malignant_traj_cldn4/scripts/extract_malignant.py --data /tmp/c4raw --out /tmp/c4_malignant.h5ad")
    lines.append("python3 methods/concordant4_malignant_traj_cldn4/scripts/analyze.py --input /tmp/c4_malignant.h5ad")
    lines.append("```")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    here = Path(__file__).resolve().parents[1]
    p.add_argument("--input", type=Path, default=Path("/tmp/c4_malignant.h5ad"))
    p.add_argument("--locked", type=Path, default=here / "data" / "locked_patient_units.tsv")
    p.add_argument("--outdir", type=Path, default=here / "results")
    p.add_argument("--finding", type=Path, default=here / "FINDING.md")
    args = p.parse_args()
    sc.settings.verbosity = 1
    np.random.seed(SEED)
    adata = sc.read_h5ad(args.input)
    print(f"loaded {adata.n_obs} cells {adata.n_vars} genes", flush=True)
    # expression of focal genes on log-normalized scale is filled after normalize
    _preprocess(adata)
    missing = {
        "AT2": _score(adata, "AT2", AT2),
        "barrier": _score(adata, "barrier", BARRIER),
        "IFN": _score(adata, "IFN", IFN),
    }
    for g in ("CLDN4", "SFTPB", "SFTPC"):
        if g not in adata.var_names:
            if g == "CLDN4":
                raise SystemExit("missing focal gene CLDN4")
            print(f"gene {g} absent from the shared matrix; not used as an expression track", flush=True)
            continue
        vec = adata[:, g].X
        if hasattr(vec, "toarray"):
            vec = np.asarray(vec.toarray()).ravel()
        else:
            vec = np.asarray(vec).ravel()
        adata.obs[f"expr_{g}"] = vec
    poles = _choose_poles(adata)
    print(json.dumps({k: poles[k] for k in poles if k != "cluster_means"}, indent=2), flush=True)
    paga = _paga_path(adata, poles["root_cluster"], poles["barrier_cluster"])
    print("PAGA", paga, flush=True)
    pal = _run_palantir(adata, poles)
    print("Palantir", {k: pal[k] for k in pal if k != "terminals"}, flush=True)
    sling = _run_slingshot(adata, poles)
    print("Slingshot", {k: sling[k] for k in sling if k != "lineages"}, flush=True)
    locked = pd.read_csv(args.locked, sep="\t")
    patients = _patient_table(adata, locked)
    tests = _tests(patients)
    # curves: all cells for palantir; prefix cells for slingshot barrier path
    bins_pal = _bin_curve(adata.obs, "palantir_pseudotime")
    prefix = adata.obs.loc[adata.obs["on_at2_barrier_prefix"]]
    bins_sl = _bin_curve(prefix, "slingshot_pt_barrier")
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "tables").mkdir(parents=True, exist_ok=True)
    patients.to_csv(args.outdir / "tables" / "patient_fate_tnk.tsv", sep="\t", index=False)
    cell_cols = [
        "dataset",
        "unit_id",
        "leiden",
        "expr_CLDN4",
        "expr_SFTPB",
        "score_AT2",
        "score_barrier",
        "score_IFN",
        "palantir_pseudotime",
        "fate_barrier",
        "fate_alt",
        "slingshot_pt_barrier",
        "slingshot_w_barrier",
        "on_at2_barrier_prefix",
    ]
    adata.obs[cell_cols].to_csv(
        args.outdir / "tables" / "cell_trajectory.tsv.gz", sep="\t", index=True
    )
    tests.to_csv(args.outdir / "tables" / "tests.tsv", sep="\t", index=False)
    bins_pal.to_csv(args.outdir / "tables" / "bins_palantir.tsv", sep="\t", index=False)
    bins_sl.to_csv(args.outdir / "tables" / "bins_slingshot.tsv", sep="\t", index=False)
    pd.DataFrame(poles["cluster_means"]).to_csv(args.outdir / "tables" / "leiden_means.tsv", sep="\t", index=False)
    pd.DataFrame({"leiden": paga.get("path", [])}).to_csv(
        args.outdir / "tables" / "paga_path.tsv", sep="\t", index=False
    )
    # paga connectivities
    cats = [str(c) for c in adata.obs["leiden"].cat.categories]
    conn = adata.uns["paga"]["connectivities"]
    arr = conn.toarray() if hasattr(conn, "toarray") else np.asarray(conn)
    conn_rows = []
    for i, a in enumerate(cats):
        for j, b in enumerate(cats):
            if j <= i:
                continue
            if arr[i, j] > 0:
                conn_rows.append({"a": a, "b": b, "connectivity": float(arr[i, j])})
    pd.DataFrame(conn_rows).to_csv(args.outdir / "tables" / "paga_connectivities.tsv", sep="\t", index=False)
    # lineage table
    pd.DataFrame(
        {
            "lineage": [f"Lineage{i+1}" for i in range(len(sling["lineages"]))],
            "clusters": [" → ".join(lin) for lin in sling["lineages"]],
            "is_barrier_lineage": [i == sling["barrier_lineage_index"] for i in range(len(sling["lineages"]))],
        }
    ).to_csv(args.outdir / "tables" / "slingshot_lineages.tsv", sep="\t", index=False)
    by = {}
    for ds, sub in patients.groupby("dataset"):
        by[ds] = {
            "cells": int(sub["n_cells"].sum()),
            "units": int(sub["unit_id"].nunique()),
            "eligible": int(sub["eligible"].sum()),
        }
    summary = {
        "not_an_ici_clock": True,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_units": int((patients["dataset"] + ":" + patients["unit_id"]).nunique()),
        "n_eligible": int(patients["eligible"].sum()),
        "cap": 180,
        "by_dataset": by,
        "missing_genes": missing,
        "poles": {k: poles[k] for k in poles if k != "cluster_means"},
        "paga_path": paga.get("path"),
        "paga_connected": paga.get("connected"),
        "palantir": pal,
        "slingshot": {k: sling[k] for k in sling if k != "lineages"} | {
            "barrier_lineage": sling["barrier_lineage"],
            "n_lineages": sling["n_lineages"],
            "chosen_rule": sling["chosen_rule"],
            "barrier_is_leaf": sling["barrier_is_leaf"],
            "engine": sling["engine"],
        },
        "tests": tests.to_dict(orient="records"),
    }
    # slingshot lineages are large but useful
    summary["slingshot_lineages"] = sling["lineages"]
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    _figures(adata, patients, bins_pal, bins_sl, paga, args.outdir)
    _write_finding(args.finding, summary, tests, bins_pal, bins_sl)
    print(tests[["contrast", "n", "rho_pooled", "p_pooled", "rho_DL", "p_DL"]].to_string(index=False), flush=True)
    print(f"wrote {args.finding}", flush=True)


if __name__ == "__main__":
    main()
