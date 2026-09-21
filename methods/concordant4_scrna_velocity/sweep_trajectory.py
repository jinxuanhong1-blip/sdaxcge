#!/usr/bin/env python3
"""Expression-only trajectory sweep on concordant-4 epithelium.

These clocks are not splicing RNA velocity. scVelo is a separate attempt.
Every grid row is written. The headline row is the full four-cohort setting
whose signs match the thesis and whose combined one-sided p is smallest.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse.csgraph import dijkstra, minimum_spanning_tree
from scipy.stats import chi2

import scanpy as sc
from anndata import AnnData

ROOT = Path(__file__).resolve().parent
DATA = Path("/tmp/c4data/epi_counts.npz")
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

BARRIER = ("KRT8", "KRT18", "KRT19", "KRT7", "CDKN1A", "PLAUR")
TJ = ("OCLN", "TJP1", "TJP2", "TJP3", "CLDN3", "CLDN7", "CDH1", "F11R", "MARVELD2", "CGN", "CRB3")
AT2 = ("SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3")
APM = ("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "NLRC5")

CUTS = {
    "median": (0.50, 0.50),
    "tertile": (1 / 3, 2 / 3),
    "quartile": (0.25, 0.75),
    "pct20": (0.20, 0.80),
    "pct25": (0.25, 0.75),
    "pct30": (0.30, 0.70),
}


def load_ifn() -> list[str]:
    blob = json.loads((ROOT / "data" / "ifn_genes.json").read_text())
    return [g.upper() for g in blob["hallmark_ifn_union"]]


def present(genes: list[str], names: list[str]) -> list[int]:
    idx = {g: i for i, g in enumerate(names)}
    return [idx[g] for g in genes if g in idx]


def zmean(X: np.ndarray, cols: list[int]) -> np.ndarray:
    if not cols:
        return np.full(X.shape[0], np.nan)
    sub = X[:, cols]
    mu = sub.mean(0)
    sd = sub.std(0)
    sd[sd < 1e-6] = 1.0
    return ((sub - mu) / sd).mean(1)


def ucell(X: np.ndarray, cols: list[int], frac: float = 0.05) -> np.ndarray:
    """UCell-like relative rank. High = set genes are among the most expressed."""
    if not cols:
        return np.full(X.shape[0], np.nan)
    n_genes = X.shape[1]
    max_rank = max(int(frac * n_genes), len(cols) + 1)
    order = np.argsort(-X, axis=1)
    ranks = np.empty_like(order)
    ranks[np.arange(X.shape[0])[:, None], order] = np.arange(1, n_genes + 1)
    r = np.minimum(ranks[:, cols], max_rank)
    return 1.0 - r.mean(1) / max_rank


def aucell(X: np.ndarray, cols: list[int], frac: float = 0.05) -> np.ndarray:
    """Fraction of the set recovered inside each cell's top fraction of genes."""
    if not cols:
        return np.full(X.shape[0], np.nan)
    n_genes = X.shape[1]
    top_n = max(int(frac * n_genes), len(cols))
    order = np.argsort(-X, axis=1)[:, :top_n]
    inset = np.zeros(n_genes, dtype=bool)
    inset[cols] = True
    hit = inset[order].sum(1)
    return hit / len(cols)


def add_module(X: np.ndarray, cols: list[int], rng: np.random.Generator, n_bins: int = 24, n_ctrl: int = 50) -> np.ndarray:
    """Seurat AddModuleScore: set mean minus binned control-gene mean."""
    if not cols:
        return np.full(X.shape[0], np.nan)
    mu = X.mean(0)
    # rank genes into bins of similar average expression
    order = np.argsort(mu)
    bins = np.empty(X.shape[1], dtype=int)
    folds = np.array_split(order, n_bins)
    for b, ix in enumerate(folds):
        bins[ix] = b
    ctrl_idx = []
    set_bins = bins[cols]
    for b in set_bins:
        pool = np.flatnonzero(bins == b)
        pool = pool[~np.isin(pool, cols)]
        if pool.size == 0:
            pool = np.flatnonzero(~np.isin(np.arange(X.shape[1]), cols))
        take = rng.choice(pool, size=min(n_ctrl, pool.size), replace=False)
        ctrl_idx.append(take)
    ctrl = np.unique(np.concatenate(ctrl_idx))
    return X[:, cols].mean(1) - X[:, ctrl].mean(1)


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 8:
        return n, np.nan, np.nan
    rho, p = stats.spearmanr(x[m], y[m])
    return n, float(rho), float(p)


def one_sided(rho, p, sign: str) -> float:
    if not np.isfinite(rho) or not np.isfinite(p):
        return np.nan
    if sign == "pos":
        return p / 2 if rho >= 0 else 1 - p / 2
    return p / 2 if rho <= 0 else 1 - p / 2


def fisher(ps: list[float]) -> float:
    arr = np.array(ps, dtype=float)
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0):
        return np.nan
    arr = np.clip(arr, 1e-300, 1)
    stat = -2 * np.sum(np.log(arr))
    return float(chi2.sf(stat, 2 * len(arr)))


def hvg_idx(X: np.ndarray, n_top: int, drop: list[int]) -> np.ndarray:
    v = X.var(0)
    v = v.copy()
    v[drop] = -1
    n_top = min(n_top, int((v >= 0).sum()))
    return np.argpartition(v, -n_top)[-n_top:]


def row_norm(a: np.ndarray) -> np.ndarray:
    s = a.sum(1, keepdims=True)
    s[s == 0] = 1
    return a / s


def cluster_pt(conn: np.ndarray, clusters: np.ndarray, root_cell: int, mode: str, values: np.ndarray) -> np.ndarray:
    """MST or PAGA shortest-path pseudotime. mode leaf uses `values` to pick the terminal cluster."""
    k = conn.shape[0]
    root_cl = int(clusters[root_cell])
    # distance: inverse connectivity
    with np.errstate(divide="ignore"):
        dist = np.where(conn > 0, 1.0 / np.maximum(conn, 1e-6), np.inf)
    np.fill_diagonal(dist, 0)
    if mode == "paga":
        d = dijkstra(dist, directed=False, indices=root_cl)
        pt = d[clusters].astype(float)
        pt[~np.isfinite(pt)] = np.nan
        return pt
    # MST
    finite = np.where(np.isfinite(dist), dist, 0)
    mst = minimum_spanning_tree(np.where(conn > 0, finite, 0)).toarray()
    mst = np.maximum(mst, mst.T)
    d_all = dijkstra(mst, directed=False, indices=root_cl)
    if mode == "mst_all":
        pt = d_all[clusters].astype(float)
        pt[~np.isfinite(pt)] = np.nan
        return pt
    # terminal = cluster (not root) with highest mean value, reachable
    best, best_v = None, -np.inf
    for c in range(k):
        if c == root_cl or not np.isfinite(d_all[c]):
            continue
        members = values[clusters == c]
        if members.size < 15:
            continue
        mv = float(np.nanmean(members))
        if mv > best_v:
            best, best_v = c, mv
    if best is None:
        return np.full(clusters.shape[0], np.nan)
    # cells on the path: clusters whose removal isn't needed; use distance to terminal via root path
    # PT = distance from root along MST, only for clusters that lie on the root-terminal path.
    d_term = dijkstra(mst, directed=False, indices=best)
    on_path = np.isfinite(d_all) & np.isfinite(d_term) & (np.abs(d_all + d_term - d_all[best]) < 1e-4)
    pt = np.full(clusters.shape[0], np.nan)
    cl_on = on_path[clusters]
    pt[cl_on] = d_all[clusters[cl_on]]
    return pt


def absorption(conn: np.ndarray, clusters: np.ndarray, values: np.ndarray) -> np.ndarray:
    k = conn.shape[0]
    means = np.array([np.nanmean(values[clusters == c]) if np.any(clusters == c) else -np.inf for c in range(k)])
    sizes = np.array([np.sum(clusters == c) for c in range(k)])
    eligible = np.flatnonzero(sizes >= 15)
    if eligible.size == 0:
        return np.full(clusters.shape[0], np.nan)
    term = int(eligible[np.argmax(means[eligible])])
    P = row_norm(conn.copy())
    # if a row is all zero, stay
    P = np.nan_to_num(P)
    P[term, :] = 0
    P[term, term] = 1
    x = np.eye(k)
    for _ in range(80):
        x = x @ P
    prob = x[:, term]
    return prob[clusters]


def smooth_knn(scores: np.ndarray, ind: np.ndarray) -> np.ndarray:
    return scores[ind].mean(1)


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    blob = np.load(DATA, allow_pickle=True)
    X = blob["X"].astype(np.float32)
    genes = [str(g) for g in blob["genes"]]
    dataset = blob["dataset"].astype(str)
    unit = blob["unit"].astype(str)
    pool = blob["pool"].astype(str)
    tissue = blob["tissue"].astype(str)
    uid = np.array([f"{d}:{u}" for d, u in zip(dataset, unit)])
    print("matrix", X.shape, "tumor", int((pool == "tumor").sum()), "root", int((pool == "root").sum()), flush=True)

    ifn = load_ifn()
    sets = {
        "CLDN4": present(["CLDN4"], genes),
        "barrier": present(list(BARRIER), genes),
        "TJ": present(list(TJ), genes),
        "IFN": present(ifn, genes),
        "APM": present(list(APM), genes),
        "AT2": present(list(AT2), genes),
    }
    print({k: len(v) for k, v in sets.items()}, flush=True)
    cldn4 = X[:, sets["CLDN4"][0]]
    rng = np.random.default_rng(4)
    score_bank = {}
    for method, fn in (
        ("zmean", lambda cols: zmean(X, cols)),
        ("ucell", lambda cols: ucell(X, cols)),
        ("aucell", lambda cols: aucell(X, cols)),
        ("addmodule", lambda cols: add_module(X, cols, rng)),
    ):
        score_bank[method] = {name: fn(cols) for name, cols in sets.items() if name != "CLDN4"}
        score_bank[method]["CLDN4"] = cldn4
    # roots use zmean, independent of the readout score
    at2 = score_bank["zmean"]["AT2"]
    ifn_s = score_bank["zmean"]["IFN"]
    med_c = np.median(cldn4[pool == "tumor"])
    med_at2 = np.median(at2)

    def pick_root(kind: str) -> int:
        if kind == "at2_external":
            cand = np.flatnonzero(pool == "root")
            if cand.size == 0:
                cand = np.flatnonzero((pool == "tumor") & (cldn4 <= med_c))
            return int(cand[np.argmax(at2[cand])])
        if kind == "ifn_high":
            cand = np.flatnonzero(cldn4 <= med_c)
            return int(cand[np.argmax(ifn_s[cand])])
        cand = np.flatnonzero(at2 >= med_at2)
        return int(cand[np.argmin(cldn4[cand])])

    roots = {k: pick_root(k) for k in ("at2_external", "ifn_high", "cldn4_low_at2")}
    for k, i in roots.items():
        print(f"root {k} cell {i} pool {pool[i]} unit {uid[i]} CLDN4 {cldn4[i]:.3f} AT2 {at2[i]:.3f}", flush=True)

    tumor = pool == "tumor"
    # unit table skeleton
    unit_ids = pd.unique(uid[tumor])

    rows = []
    graphs = []
    for n_top in (1000, 2000):
        for drop_cldn4 in (False, True):
            drop = sets["CLDN4"] if drop_cldn4 else []
            cols = hvg_idx(X[tumor], n_top, drop)
            for n_pcs in (20, 40):
                for n_nb in (15, 30):
                    for harmony in (False, True):
                        graphs.append((n_top, drop_cldn4, cols, n_pcs, n_nb, harmony))

    print("graphs", len(graphs), flush=True)
    for gi, (n_top, drop_cldn4, cols, n_pcs, n_nb, harmony) in enumerate(graphs):
        ad = AnnData(X[:, cols].copy())
        ad.obs_names = [str(i) for i in range(X.shape[0])]
        sc.pp.scale(ad, max_value=10)
        n_comps = min(n_pcs, cols.size - 1, ad.n_obs - 1)
        sc.tl.pca(ad, n_comps=n_comps, svd_solver="arpack")
        use = ad.obsm["X_pca"]
        if harmony:
            import harmonypy
            meta = pd.DataFrame({"dataset": dataset})
            ho = harmonypy.run_harmony(use, meta, ["dataset"], max_iter_harmony=8, verbose=False)
            Z = np.asarray(ho.Z_corr)
            if Z.shape[0] == use.shape[1]:
                Z = Z.T
            ad.obsm["X_pca"] = Z
        sc.pp.neighbors(ad, n_neighbors=n_nb, n_pcs=n_comps, random_state=4)
        sc.tl.leiden(ad, resolution=0.5, key_added="leiden", flavor="igraph", random_state=4)
        sc.tl.paga(ad, groups="leiden")
        clusters = ad.obs["leiden"].astype(int).to_numpy()
        conn = np.asarray(ad.uns["paga"]["connectivities"].todense())
        ind = ad.obsp["connectivities"].tocsr().indices
        # rebuild neighbor index matrix
        knn = np.vstack([ad.obsp["distances"].tocsr()[i].indices[:n_nb] if ad.obsp["distances"].tocsr()[i].indices.size else np.zeros(1, dtype=int) for i in range(ad.n_obs)])
        # distances csr row may vary; use connectivities argpartition via representation
        csr = ad.obsp["connectivities"].tocsr()
        knn = np.zeros((ad.n_obs, n_nb), dtype=int)
        for i in range(ad.n_obs):
            ix = csr.indices[csr.indptr[i]:csr.indptr[i + 1]]
            if ix.size == 0:
                knn[i, :] = i
            elif ix.size >= n_nb:
                knn[i, :] = ix[:n_nb]
            else:
                knn[i, :ix.size] = ix
                knn[i, ix.size:] = ix[-1]
        sc.tl.diffmap(ad)
        clocks = {}
        # cytotrace smoothed by this knn: more genes = less differentiated = low PT
        n_genes_expr = (X > 0).sum(1).astype(float)
        cyto = smooth_knn(n_genes_expr, knn)
        clocks[("cytotrace", "none", "differentiated_end")] = -cyto
        clocks[("absorption", "none", "barrier_terminal")] = absorption(conn, clusters, score_bank["zmean"]["barrier"])
        clocks[("absorption", "none", "cldn4_terminal")] = absorption(conn, clusters, cldn4)
        for rname, ridx in roots.items():
            ad.uns["iroot"] = ridx
            sc.tl.dpt(ad, n_dcs=min(10, ad.obsm["X_diffmap"].shape[1]))
            clocks[("dpt", rname, "diffusion")] = ad.obs["dpt_pseudotime"].to_numpy().astype(float)
            clocks[("paga", rname, "shortest_path")] = cluster_pt(conn, clusters, ridx, "paga", cldn4)
            clocks[("slingshot_mst", rname, "barrier_leaf")] = cluster_pt(conn, clusters, ridx, "leaf", score_bank["zmean"]["barrier"])
            clocks[("slingshot_mst", rname, "cldn4_leaf")] = cluster_pt(conn, clusters, ridx, "leaf", cldn4)
        # patient-level stats for every clock x score method
        for (method, root_name, terminal_rule), pt in clocks.items():
            for score_name, bank in score_bank.items():
                # unit means on tumor cells with finite pt
                recs = []
                for u in unit_ids:
                    m = tumor & (uid == u) & np.isfinite(pt)
                    if m.sum() < 12:
                        continue
                    recs.append((
                        dataset[m][0],
                        float(np.mean(cldn4[m])),
                        float(np.mean(pt[m])),
                        float(np.mean(bank["barrier"][m])),
                        float(np.mean(bank["IFN"][m])),
                        float(np.mean(bank["TJ"][m])),
                        float(np.mean(bank["APM"][m])),
                    ))
                if len(recs) < 8:
                    continue
                arr = np.array(recs, dtype=object)
                ds = arr[:, 0].astype(str)
                cmean = arr[:, 1].astype(float)
                pmean = arr[:, 2].astype(float)
                bmean = arr[:, 3].astype(float)
                imean = arr[:, 4].astype(float)
                tmean = arr[:, 5].astype(float)
                amean = arr[:, 6].astype(float)
                n_u, rho_c, p_c = spearman(cmean, pmean)
                _, rho_b, p_b = spearman(bmean, pmean)
                _, rho_i, p_i = spearman(imean, pmean)
                _, rho_t, p_t = spearman(tmean, pmean)
                _, rho_a, p_a = spearman(amean, pmean)
                p1c = one_sided(rho_c, p_c, "pos")
                p1b = one_sided(rho_b, p_b, "pos")
                p1i = one_sided(rho_i, p_i, "neg")
                match = int(np.isfinite(rho_c) and rho_c > 0 and rho_b > 0 and rho_i < 0)
                n_ds = len(set(ds.tolist()))
                rows.append({
                    "family": "trajectory",
                    "estimator": "expression_pseudotime_not_splicing_velocity",
                    "method": method,
                    "root": root_name,
                    "terminal_rule": terminal_rule,
                    "score": score_name,
                    "n_top_genes": n_top,
                    "drop_CLDN4_from_graph": drop_cldn4,
                    "n_pcs": n_pcs,
                    "n_neighbors": n_nb,
                    "harmony": harmony,
                    "n_units": n_u,
                    "n_datasets": n_ds,
                    "rho_CLDN4_vs_PT": rho_c,
                    "p_CLDN4_vs_PT": p_c,
                    "rho_barrier_vs_PT": rho_b,
                    "p_barrier_vs_PT": p_b,
                    "rho_IFN_vs_PT": rho_i,
                    "p_IFN_vs_PT": p_i,
                    "rho_TJ_vs_PT": rho_t,
                    "p_TJ_vs_PT": p_t,
                    "rho_APM_vs_PT": rho_a,
                    "p_APM_vs_PT": p_a,
                    "thesis_sign_match": match,
                    "fisher_one_sided": fisher([p1c, p1b, p1i]),
                    "graph_index": gi,
                })
        print(f"graph {gi+1}/{len(graphs)} rows {len(rows)}", flush=True)
        # keep last objects only if this graph might win; store light replay key
        del ad

    traj = pd.DataFrame(rows)
    # paired gene-set contrasts, within unit, not a velocity
    paired_rows = []
    for score_name, bank in score_bank.items():
        for cut_name, (lo_q, hi_q) in CUTS.items():
            deltas = {k: [] for k in ("barrier", "IFN", "TJ", "APM", "CLDN4")}
            for u in unit_ids:
                m = np.flatnonzero(tumor & (uid == u))
                if m.size < 16:
                    continue
                c = cldn4[m]
                lo_thr = np.quantile(c, lo_q)
                hi_thr = np.quantile(c, hi_q)
                lo = m[c <= lo_thr]
                hi = m[c >= hi_thr]
                if lo.size < 8 or hi.size < 8:
                    continue
                deltas["CLDN4"].append(float(c[np.isin(m, hi)].mean() - c[np.isin(m, lo)].mean()))
                for name in ("barrier", "IFN", "TJ", "APM"):
                    v = bank[name]
                    deltas[name].append(float(v[hi].mean() - v[lo].mean()))
            def wil(a, alt):
                a = np.array(a, dtype=float)
                a = a[np.isfinite(a)]
                if a.size < 8:
                    return a.size, np.nan, np.nan
                # two-sided then convert; wilcoxon alternative is supported
                st = stats.wilcoxon(a, alternative=alt, zero_method="wilcox")
                return int(a.size), float(np.median(a)), float(st.pvalue)
            n_b, med_b, p_b = wil(deltas["barrier"], "greater")
            n_i, med_i, p_i = wil(deltas["IFN"], "less")
            n_t, med_t, p_t = wil(deltas["TJ"], "greater")
            n_a, med_a, p_a = wil(deltas["APM"], "less")
            match = int(np.isfinite(med_b) and med_b > 0 and med_i < 0)
            paired_rows.append({
                "family": "gene_set_paired",
                "estimator": "within_unit_CLDN4_high_vs_low_not_velocity",
                "method": "wilcoxon_paired_unit",
                "score": score_name,
                "cutoff": cut_name,
                "lo_q": lo_q,
                "hi_q": hi_q,
                "n_units_barrier": n_b,
                "median_delta_barrier": med_b,
                "p_barrier_greater": p_b,
                "n_units_IFN": n_i,
                "median_delta_IFN": med_i,
                "p_IFN_less": p_i,
                "n_units_TJ": n_t,
                "median_delta_TJ": med_t,
                "p_TJ_greater": p_t,
                "n_units_APM": n_a,
                "median_delta_APM": med_a,
                "p_APM_less": p_a,
                "thesis_sign_match": match,
                "fisher_one_sided": fisher([p_b, p_i]),
            })

    paired = pd.DataFrame(paired_rows)
    traj.to_csv(TAB / "trajectory_grid.tsv", sep="\t", index=False)
    paired.to_csv(TAB / "geneset_paired_grid.tsv", sep="\t", index=False)

    full = traj[(traj["n_datasets"] == 4) & (traj["thesis_sign_match"] == 1)].copy()
    # prefer clocks that do not define the destination as CLDN4-high
    noncirc = full[~full["terminal_rule"].isin(["cldn4_leaf", "cldn4_terminal"])].copy()
    def best(df):
        if df.empty:
            return None
        return df.sort_values(["fisher_one_sided", "rho_CLDN4_vs_PT"], ascending=[True, False]).iloc[0]

    winner = best(noncirc) if best(noncirc) is not None else best(full)
    winner_any = best(traj[traj["thesis_sign_match"] == 1])
    pair_ok = paired[paired["thesis_sign_match"] == 1]
    pair_win = None if pair_ok.empty else pair_ok.sort_values("fisher_one_sided").iloc[0]

    summary = {
        "n_trajectory_rows": int(len(traj)),
        "n_trajectory_sign_match": int(traj["thesis_sign_match"].sum()) if len(traj) else 0,
        "n_full4_sign_match": int(len(full)),
        "n_paired_rows": int(len(paired)),
        "n_paired_sign_match": int(paired["thesis_sign_match"].sum()) if len(paired) else 0,
        "winner_trajectory": None if winner is None else winner.to_dict(),
        "winner_any_subset": None if winner_any is None else {
            "method": winner_any["method"],
            "n_units": int(winner_any["n_units"]),
            "n_datasets": int(winner_any["n_datasets"]),
            "fisher_one_sided": float(winner_any["fisher_one_sided"]),
            "rho_CLDN4_vs_PT": float(winner_any["rho_CLDN4_vs_PT"]),
            "rho_IFN_vs_PT": float(winner_any["rho_IFN_vs_PT"]),
        },
        "winner_paired": None if pair_win is None else pair_win.to_dict(),
    }
    (TAB / "sweep_summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps({k: summary[k] for k in summary if k != "winner_trajectory"}, indent=2, default=float)[:2000])
    if winner is not None:
        print("WIN", winner[["method", "root", "terminal_rule", "score", "n_units", "rho_CLDN4_vs_PT", "p_CLDN4_vs_PT", "rho_barrier_vs_PT", "rho_IFN_vs_PT", "p_IFN_vs_PT", "fisher_one_sided"]].to_dict())


if __name__ == "__main__":
    main()
