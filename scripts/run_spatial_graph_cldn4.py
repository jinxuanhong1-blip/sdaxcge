#!/usr/bin/env python3
"""CLDN4-only spatial neighbor-graph analysis.

Official CosMx NSCLC (He et al. 2022; Lung5_Rep2) + one downloadable
Visium LUAD (GSE307534 GSM9226169 P1 invasive section).

No private 8-KL data. CLDN4 is the only tight-junction query gene.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
import statsmodels.api as sm
from anndata import AnnData
from sklearn.neighbors import BallTree, radius_neighbors_graph
from tqdm import tqdm

LOG = logging.getLogger("cldn4_spatial")
UM_PER_PX_COSMX = 0.18  # official SMI-ReadMe for this prototype NSCLC release
RADIUS_UM = 50.0
SEED = 7
N_PERMS_ENRICH = 400
N_PERMS_COMP = 1000
COOC_QUERY_CAP = 2500
COOC_BINS_UM = np.array([0, 25, 50, 75, 100, 150, 200, 300, 400, 500], dtype=float)
VISIUM_COOC_BINS_UM = np.array([0, 100, 200, 300, 400, 500, 700, 900], dtype=float)

COSMX_DIR = Path(
    os.environ.get(
        "COSMX_DIR",
        "/tmp/spatial_data/cosmx/Lung5_Rep2/Lung5_Rep2-Flat_files_and_images",
    )
)
VISIUM_DIR = Path(os.environ.get("VISIUM_DIR", "/tmp/spatial_data/visium/P1_LUAD"))
OUT = Path(os.environ.get("OUT_DIR", "/workspace/results"))

LABEL_ORDER = [
    "CLDN4-high tumor/epi",
    "CLDN4-low tumor/epi",
    "CD8",
    "NK",
    "Treg",
    "Mac",
    "Other",
]
IMMUNE = ["CD8", "NK", "Treg", "Mac"]
FOCUS = ["CLDN4-high tumor/epi", "CLDN4-low tumor/epi", *IMMUNE]

MARKERS = {
    "tumor_epi": ["KRT8", "KRT18", "KRT19", "KRT7", "EPCAM"],
    "CD8": ["CD8A", "CD8B"],
    "NK": ["NKG7", "GNLY", "KLRK1", "NCR1", "PRF1"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4"],
    "Mac": ["CD68", "CD163", "MRC1", "MARCO", "C1QA"],
    "T": ["CD3D", "CD3E"],
    "immune": ["PTPRC"],
}


def setup() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="ticks", context="talk")
    sc.settings.verbosity = 1
    np.random.seed(SEED)


def present(adata: AnnData, genes: list[str]) -> list[str]:
    return [g for g in genes if g in adata.var_names]


def log1p_counts(adata: AnnData) -> np.ndarray:
    x = adata.X
    if sp.issparse(x):
        x = x.toarray()
    return np.log1p(np.asarray(x, dtype=np.float64))


def module_score(logx: np.ndarray, var_names: pd.Index, genes: list[str]) -> np.ndarray:
    idx = [var_names.get_loc(g) for g in genes if g in var_names]
    if not idx:
        return np.zeros(logx.shape[0], dtype=np.float64)
    return logx[:, idx].mean(axis=1)


def zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=np.float64)
    sd = np.nanstd(v)
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(v)
    return (v - np.nanmean(v)) / sd


def gene_vec(adata: AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.zeros(adata.n_obs, dtype=np.float64)
    x = adata[:, gene].X
    if sp.issparse(x):
        x = x.toarray()
    return np.asarray(x).reshape(-1).astype(np.float64)


def assign_labels_cosmx(adata: AnnData) -> None:
    logx = log1p_counts(adata)
    names = adata.var_names
    scores = {k: module_score(logx, names, v) for k, v in MARKERS.items()}
    cldn4 = gene_vec(adata, "CLDN4")
    foxp3 = gene_vec(adata, "FOXP3")
    cd8a = gene_vec(adata, "CD8A")
    adata.obs["score_tumor_epi"] = scores["tumor_epi"]
    adata.obs["score_CD8"] = scores["CD8"]
    adata.obs["score_NK"] = scores["NK"]
    adata.obs["score_Treg"] = scores["Treg"]
    adata.obs["score_Mac"] = scores["Mac"]
    adata.obs["CLDN4"] = cldn4
    adata.obs["KRT8"] = gene_vec(adata, "KRT8")
    adata.obs["CD8A"] = cd8a

    panck = zscore(adata.obs["Mean.PanCK"].to_numpy())
    cd45 = zscore(adata.obs["Mean.CD45"].to_numpy())
    cd3 = zscore(adata.obs["Mean.CD3"].to_numpy())
    tumor_z = zscore(scores["tumor_epi"])
    cd8_z = zscore(scores["CD8"])
    nk_z = zscore(scores["NK"])
    treg_z = zscore(scores["Treg"])
    mac_z = zscore(scores["Mac"])

    n = adata.n_obs
    lab = np.array(["Other"] * n, dtype=object)
    tumor = (panck > 0.4) | ((tumor_z > 0.6) & (cd45 < 0.8))
    treg = (foxp3 >= 1) & ((gene_vec(adata, "IL2RA") >= 1) | (gene_vec(adata, "CTLA4") >= 1) | (cd3 > 0.5))
    cd8 = (cd8a >= 1) & ((cd3 > 0.0) | (scores["T"] > 0) | (cd8_z > 0.5)) & ~treg
    # NK: cytotoxic module without CD8/TCR; NKG7 alone is not enough (shared with T)
    nkg7 = gene_vec(adata, "NKG7")
    gnly = gene_vec(adata, "GNLY")
    ncr1 = gene_vec(adata, "NCR1")
    nk = ((gnly >= 1) | (ncr1 >= 1) | ((nkg7 >= 1) & (nk_z > 1.2))) & (cd8a == 0) & ~treg & (scores["T"] < 0.25) & (cd3 < 0.4)
    mac = (mac_z > 0.9) & (panck < 0.5) & (cd45 > -0.2)
    # priority: Treg > CD8 > NK > Mac > tumor > other
    lab[mac] = "Mac"
    lab[tumor] = "tumor/epi"
    lab[nk] = "NK"
    lab[cd8] = "CD8"
    lab[treg] = "Treg"

    is_tumor = lab == "tumor/epi"
    med = np.median(cldn4[is_tumor]) if is_tumor.any() else np.median(cldn4)
    lab[is_tumor & (cldn4 >= med)] = "CLDN4-high tumor/epi"
    lab[is_tumor & (cldn4 < med)] = "CLDN4-low tumor/epi"
    adata.obs["graph_label"] = pd.Categorical(lab, categories=LABEL_ORDER).remove_unused_categories()
    adata.obs["is_tumor"] = adata.obs["graph_label"].isin(
        ["CLDN4-high tumor/epi", "CLDN4-low tumor/epi"]
    )
    adata.uns["cldn4_tumor_median"] = float(med)
    adata.uns["label_counts"] = adata.obs["graph_label"].value_counts().to_dict()


def assign_labels_visium(adata: AnnData) -> None:
    logx = log1p_counts(adata)
    names = adata.var_names
    scores = {k: module_score(logx, names, v) for k, v in MARKERS.items()}
    cldn4 = gene_vec(adata, "CLDN4")
    adata.obs["score_tumor_epi"] = scores["tumor_epi"]
    adata.obs["score_CD8"] = scores["CD8"]
    adata.obs["score_NK"] = scores["NK"]
    adata.obs["score_Treg"] = scores["Treg"]
    adata.obs["score_Mac"] = scores["Mac"]
    adata.obs["CLDN4"] = cldn4
    adata.obs["KRT8"] = gene_vec(adata, "KRT8")
    adata.obs["CD8A"] = gene_vec(adata, "CD8A")

    tumor_z = zscore(scores["tumor_epi"])
    cd8_z = zscore(scores["CD8"])
    nk_z = zscore(scores["NK"])
    treg_z = zscore(scores["Treg"])
    mac_z = zscore(scores["Mac"])
    cldn4_z = zscore(cldn4)

    n = adata.n_obs
    lab = np.array(["Other"] * n, dtype=object)
    immune_mat = np.vstack([cd8_z, nk_z, treg_z, mac_z])
    immune_best = immune_mat.argmax(axis=0)
    immune_max = immune_mat.max(axis=0)
    tumor = tumor_z >= np.quantile(tumor_z, 0.60)
    lab[immune_max >= 0.85] = np.array(IMMUNE)[immune_best[immune_max >= 0.85]]
    lab[tumor] = "tumor/epi"
    # if a tumor spot is strongly immune-dominated, keep immune
    strong_immune = immune_max >= 1.4
    lab[strong_immune & ~tumor] = np.array(IMMUNE)[immune_best[strong_immune & ~tumor]]

    is_tumor = lab == "tumor/epi"
    med = np.median(cldn4[is_tumor]) if is_tumor.any() else np.median(cldn4)
    # prefer CLDN4 z among tumor spots so high/low is within epithelium
    lab[is_tumor & (cldn4 >= med)] = "CLDN4-high tumor/epi"
    lab[is_tumor & (cldn4 < med)] = "CLDN4-low tumor/epi"
    adata.obs["graph_label"] = pd.Categorical(lab, categories=LABEL_ORDER).remove_unused_categories()
    adata.obs["is_tumor"] = adata.obs["graph_label"].isin(
        ["CLDN4-high tumor/epi", "CLDN4-low tumor/epi"]
    )
    adata.obs["CLDN4_z"] = cldn4_z
    adata.uns["cldn4_tumor_median"] = float(med)
    adata.uns["label_counts"] = adata.obs["graph_label"].value_counts().to_dict()


def load_cosmx(path: Path) -> AnnData:
    LOG.info("Loading CosMx Lung5_Rep2 from %s", path)
    expr = pd.read_csv(path / "Lung5_Rep2_exprMat_file.csv")
    meta = pd.read_csv(path / "Lung5_Rep2_metadata_file.csv")
    expr = expr.loc[expr["cell_ID"] != 0].copy()
    meta = meta.loc[meta["cell_ID"] != 0].copy()
    key = ["fov", "cell_ID"]
    merged = expr.merge(meta, on=key, how="inner")
    gene_cols = [c for c in expr.columns if c not in key]
    counts = merged[gene_cols].to_numpy(dtype=np.float32)
    obs = merged.drop(columns=gene_cols)
    obs.index = [f"{int(f)}_{int(c)}" for f, c in zip(obs["fov"], obs["cell_ID"])]
    adata = AnnData(X=sp.csr_matrix(counts), obs=obs, var=pd.DataFrame(index=gene_cols))
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    keep = (
        (adata.obs["n_counts"] >= 20)
        & (adata.obs["n_genes"] >= 10)
        & (adata.obs["Area"] >= 80)
        & (adata.obs["Area"] <= 25000)
    )
    adata = adata[keep].copy()
    xy_um = np.column_stack(
        [
            adata.obs["CenterX_global_px"].to_numpy() * UM_PER_PX_COSMX,
            adata.obs["CenterY_global_px"].to_numpy() * UM_PER_PX_COSMX,
        ]
    )
    adata.obsm["spatial"] = xy_um
    assign_labels_cosmx(adata)
    LOG.info("CosMx cells after QC: %s; labels: %s", adata.n_obs, adata.uns["label_counts"])
    return adata


def load_visium(path: Path) -> AnnData:
    LOG.info("Loading Visium GSE307534 P1_LUAD from %s", path)
    adata = sc.read_10x_mtx(path / "filtered_feature_bc_matrix", var_names="gene_symbols")
    adata.var_names_make_unique()
    pos_path = path / "spatial" / "tissue_positions.csv"
    pos = pd.read_csv(pos_path)
    if pos.columns[0].lower().startswith("barcode"):
        pos = pos.set_index(pos.columns[0])
    else:
        pos.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ]
        pos = pos.set_index("barcode")
    pos.index = pos.index.astype(str)
    adata.obs = adata.obs.join(pos, how="left")
    adata = adata[(adata.obs["in_tissue"] == 1) & (np.asarray(adata.X.sum(axis=1)).ravel() > 0)].copy()
    adata.obsm["spatial"] = np.column_stack(
        [
            adata.obs["pxl_col_in_fullres"].to_numpy(),
            adata.obs["pxl_row_in_fullres"].to_numpy(),
        ]
    )
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    assign_labels_visium(adata)
    LOG.info("Visium spots: %s; labels: %s", adata.n_obs, adata.uns["label_counts"])
    return adata


def build_radius_graph(adata: AnnData, radius_um: float) -> None:
    xy = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    conn = radius_neighbors_graph(xy, radius=radius_um, mode="connectivity", include_self=False)
    conn = conn.tocsr()
    adata.obsp["spatial_connectivities"] = conn
    adata.uns["spatial_neighbors"] = {
        "connectivities_key": "spatial_connectivities",
        "params": {"radius_um": radius_um, "coord_type": "generic"},
    }
    deg = np.asarray(conn.sum(axis=1)).ravel()
    adata.obs["n_spatial_neighbors"] = deg
    LOG.info(
        "Radius graph %.1f µm: median degree %.1f (mean %.1f)",
        radius_um,
        float(np.median(deg)),
        float(np.mean(deg)),
    )


def build_visium_hop_graph(adata: AnnData, hops: int = 2) -> None:
    """1–2 hop hex grid graph from Visium array coordinates."""
    rows = adata.obs["array_row"].to_numpy().astype(int)
    cols = adata.obs["array_col"].to_numpy().astype(int)
    key_to_i = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(rows, cols))}
    # Official Visium hex: same-row neighbors are ±2 columns.
    offsets = [(0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    ii, jj = [], []
    for i, (r, c) in enumerate(zip(rows, cols)):
        for dr, dc in offsets:
            j = key_to_i.get((r + dr, c + dc))
            if j is not None:
                ii.append(i)
                jj.append(j)
    n = adata.n_obs
    hop1 = sp.csr_matrix((np.ones(len(ii), dtype=np.float32), (ii, jj)), shape=(n, n))
    hop1 = hop1.maximum(hop1.T)
    hop1.setdiag(0)
    hop1.eliminate_zeros()
    if hops >= 2:
        hop2 = hop1 @ hop1
        hop2.data = np.ones_like(hop2.data)
        conn = hop1.maximum(hop2)
        conn.setdiag(0)
        conn.eliminate_zeros()
    else:
        conn = hop1
    adata.obsp["spatial_connectivities"] = conn.tocsr()
    adata.obsp["spatial_connectivities_1hop"] = hop1.tocsr()
    adata.uns["spatial_neighbors"] = {
        "connectivities_key": "spatial_connectivities",
        "params": {"hops": hops, "coord_type": "visium_hex"},
    }
    deg = np.asarray(conn.sum(axis=1)).ravel()
    adata.obs["n_spatial_neighbors"] = deg
    LOG.info("Visium %d-hop graph: median degree %.1f", hops, float(np.median(deg)))


def nhood_enrichment_z(
    adata: AnnData, cluster_key: str = "graph_label", n_perms: int = N_PERMS_ENRICH
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Permutation z-score of cross-type spatial edges (squidpy-equivalent)."""
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    cats = [c for c in LABEL_ORDER if c in set(labels)]
    cat_index = {c: i for i, c in enumerate(cats)}
    codes = np.array([cat_index[x] for x in labels], dtype=np.int32)
    k = len(cats)
    conn = adata.obsp["spatial_connectivities"].tocoo()
    src, dst = conn.row, conn.col
    keep = src < dst
    src, dst = src[keep], dst[keep]

    def count_mat(lab_codes: np.ndarray) -> np.ndarray:
        mat = np.zeros((k, k), dtype=np.float64)
        a = lab_codes[src]
        b = lab_codes[dst]
        np.add.at(mat, (a, b), 1)
        np.add.at(mat, (b, a), 1)
        return mat

    obs = count_mat(codes)
    rng = np.random.default_rng(SEED)
    perm_stack = np.zeros((n_perms, k, k), dtype=np.float64)
    for p in tqdm(range(n_perms), desc="nhood enrichment perms"):
        perm = codes.copy()
        rng.shuffle(perm)
        perm_stack[p] = count_mat(perm)
    mu = perm_stack.mean(axis=0)
    sd = perm_stack.std(axis=0, ddof=1)
    z = np.zeros_like(obs)
    np.divide(obs - mu, sd, out=z, where=sd > 0)
    z_df = pd.DataFrame(z, index=cats, columns=cats)
    c_df = pd.DataFrame(obs, index=cats, columns=cats)
    adata.uns[f"{cluster_key}_nhood_enrichment"] = {"zscore": z, "count": obs, "categories": cats}
    return z_df, c_df


def try_squidpy_enrichment(adata: AnnData, cluster_key: str = "graph_label") -> pd.DataFrame | None:
    try:
        import squidpy as sq

        sq.gr.nhood_enrichment(
            adata,
            cluster_key=cluster_key,
            n_perms=N_PERMS_ENRICH,
            seed=SEED,
            show_progress_bar=True,
        )
        uns = adata.uns[f"{cluster_key}_nhood_enrichment"]
        cats = list(adata.obs[cluster_key].cat.categories)
        z = pd.DataFrame(uns["zscore"], index=cats, columns=cats)
        LOG.info("Used squidpy.gr.nhood_enrichment")
        return z
    except Exception as exc:
        LOG.warning("squidpy nhood_enrichment unavailable (%s); using local permutation", exc)
        z, _ = nhood_enrichment_z(adata, cluster_key=cluster_key)
        return z


def neighbor_type_fracs(conn: sp.csr_matrix, labels: np.ndarray, types: list[str]) -> dict[str, np.ndarray]:
    """Per-cell fraction of neighbors of each type (vectorized)."""
    deg = np.asarray(conn.sum(axis=1)).ravel()
    deg_safe = np.maximum(deg, 1.0)
    out = {}
    for t in types:
        ind = (labels == t).astype(np.float64)
        out[t] = np.asarray(conn @ ind).ravel() / deg_safe
        out[t][deg == 0] = np.nan
    return out


def composition_test(adata: AnnData, n_perms: int = N_PERMS_COMP) -> pd.DataFrame:
    labels = adata.obs["graph_label"].astype(str).to_numpy()
    conn = adata.obsp["spatial_connectivities"].tocsr()
    high = labels == "CLDN4-high tumor/epi"
    low = labels == "CLDN4-low tumor/epi"
    tumor = high | low
    fr = neighbor_type_fracs(conn, labels, IMMUNE)
    has_deg = ~np.isnan(next(iter(fr.values())))

    def mean_frac(mask: np.ndarray) -> dict[str, float]:
        m = mask & has_deg
        return {t: float(np.nanmean(fr[t][m])) if m.any() else 0.0 for t in IMMUNE}

    obs_h = mean_frac(high)
    obs_l = mean_frac(low)
    rng = np.random.default_rng(SEED)
    tumor_idx = np.flatnonzero(tumor & has_deg)
    n_high = int((high & has_deg).sum())
    n_high = min(n_high, tumor_idx.size)
    deltas = {t: np.empty(n_perms, dtype=np.float64) for t in IMMUNE}
    for p in tqdm(range(n_perms), desc="composition perms"):
        pick = rng.choice(tumor_idx, size=n_high, replace=False)
        fake_high = np.zeros(len(labels), dtype=bool)
        fake_high[pick] = True
        fake_low = tumor & has_deg & ~fake_high
        for t in IMMUNE:
            deltas[t][p] = float(np.nanmean(fr[t][fake_high])) - float(np.nanmean(fr[t][fake_low]))
    rows = []
    for t in IMMUNE:
        delta = obs_h[t] - obs_l[t]
        null = deltas[t]
        pval = (np.sum(np.abs(null) >= abs(delta)) + 1) / (n_perms + 1)
        rows.append(
            {
                "neighbor_type": t,
                "frac_around_CLDN4_high": obs_h[t],
                "frac_around_CLDN4_low": obs_l[t],
                "delta_high_minus_low": delta,
                "perm_p": pval,
                "null_mean": float(null.mean()),
                "null_sd": float(null.std(ddof=1)),
            }
        )
    return pd.DataFrame(rows)


def cooccurrence_vs_distance(
    adata: AnnData,
    query: str = "CLDN4-high tumor/epi",
    targets: list[str] | None = None,
    bins: np.ndarray = COOC_BINS_UM,
    query_cap: int = COOC_QUERY_CAP,
) -> pd.DataFrame:
    targets = targets or IMMUNE + ["CLDN4-low tumor/epi"]
    labels = adata.obs["graph_label"].astype(str).to_numpy()
    xy = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    q_idx = np.flatnonzero(labels == query)
    if q_idx.size == 0:
        return pd.DataFrame()
    rng = np.random.default_rng(SEED)
    if q_idx.size > query_cap:
        q_idx = rng.choice(q_idx, size=query_cap, replace=False)
    tree = BallTree(xy)
    glob = {t: float(np.mean(labels == t)) for t in targets}
    rows = []
    r_max = float(bins[-1])
    # query all neighbors once
    ind_list = tree.query_radius(xy[q_idx], r=r_max)
    dist_list = tree.query_radius(xy[q_idx], r=r_max, return_distance=True)[1]
    for lo, hi in zip(bins[:-1], bins[1:]):
        counts = {t: 0 for t in targets}
        n_annulus = 0
        for neigh, d in zip(ind_list, dist_list):
            in_bin = (d > lo) & (d <= hi)
            sel = neigh[in_bin]
            if sel.size == 0:
                continue
            n_annulus += sel.size
            nl = labels[sel]
            for t in targets:
                counts[t] += int(np.sum(nl == t))
        for t in targets:
            p_cond = (counts[t] / n_annulus) if n_annulus else np.nan
            score = p_cond / glob[t] if glob[t] > 0 else np.nan
            rows.append(
                {
                    "query": query,
                    "target": t,
                    "d_lo_um": lo,
                    "d_hi_um": hi,
                    "d_mid_um": 0.5 * (lo + hi),
                    "n_annulus": n_annulus,
                    "p_cond": p_cond,
                    "p_global": glob[t],
                    "cooccurrence": score,
                }
            )
    return pd.DataFrame(rows)


def misty_like_cd8(adata: AnnData) -> dict:
    """How much of CD8A is explained by neighboring CLDN4 after neighboring KRT8."""
    conn = adata.obsp["spatial_connectivities"].tocsr()
    cldn4 = adata.obs["CLDN4"].to_numpy(dtype=np.float64)
    krt8 = adata.obs["KRT8"].to_numpy(dtype=np.float64)
    cd8 = np.log1p(adata.obs["CD8A"].to_numpy(dtype=np.float64))
    deg = np.asarray(conn.sum(axis=1)).ravel()
    deg_safe = np.maximum(deg, 1.0)
    neigh_c = np.asarray(conn @ cldn4).ravel() / deg_safe
    neigh_k = np.asarray(conn @ krt8).ravel() / deg_safe
    neigh_cd8 = np.asarray(conn @ cd8).ravel() / deg_safe
    use = deg > 0
    tumor = adata.obs["is_tumor"].to_numpy()

    def ols_fit(y, k, c) -> dict:
        Xk = sm.add_constant(pd.DataFrame({"neigh_KRT8": k}))
        Xkc = sm.add_constant(pd.DataFrame({"neigh_KRT8": k, "neigh_CLDN4": c}))
        m_k = sm.OLS(y, Xk, missing="drop").fit()
        m_kc = sm.OLS(y, Xkc, missing="drop").fit()
        return {
            "n": int(len(y)),
            "r2_krt8_only": float(m_k.rsquared),
            "r2_krt8_plus_cldn4": float(m_kc.rsquared),
            "delta_r2_cldn4_after_krt8": float(m_kc.rsquared - m_k.rsquared),
            "coef_neigh_CLDN4": float(m_kc.params.get("neigh_CLDN4", np.nan)),
            "pval_neigh_CLDN4": float(m_kc.pvalues.get("neigh_CLDN4", np.nan)),
            "coef_neigh_KRT8_full": float(m_kc.params.get("neigh_KRT8", np.nan)),
            "pval_neigh_KRT8_full": float(m_kc.pvalues.get("neigh_KRT8", np.nan)),
        }

    # 1) juxta: CD8A_i ~ neigh KRT8 + neigh CLDN4
    all_cells = ols_fit(cd8[use], neigh_k[use], neigh_c[use])
    # 2) same, non-tumor cells only (CD8 lives in stroma)
    nt = use & ~tumor
    nontumor = ols_fit(cd8[nt], neigh_k[nt], neigh_c[nt]) if nt.sum() > 20 else {}
    # 3) tumor-centric: local CD8 neighborhood ~ own KRT8 + own CLDN4
    tm = use & tumor
    tumor_local_cd8 = {}
    if tm.sum() > 20:
        Xk = sm.add_constant(pd.DataFrame({"own_KRT8": krt8[tm]}))
        Xkc = sm.add_constant(pd.DataFrame({"own_KRT8": krt8[tm], "own_CLDN4": cldn4[tm]}))
        y = neigh_cd8[tm]
        m_k = sm.OLS(y, Xk, missing="drop").fit()
        m_kc = sm.OLS(y, Xkc, missing="drop").fit()
        tumor_local_cd8 = {
            "n": int(tm.sum()),
            "r2_krt8_only": float(m_k.rsquared),
            "r2_krt8_plus_cldn4": float(m_kc.rsquared),
            "delta_r2_cldn4_after_krt8": float(m_kc.rsquared - m_k.rsquared),
            "coef_own_CLDN4": float(m_kc.params.get("own_CLDN4", np.nan)),
            "pval_own_CLDN4": float(m_kc.pvalues.get("own_CLDN4", np.nan)),
        }
    return {
        "juxta_all_cells": all_cells,
        "juxta_nontumor": nontumor,
        "tumor_own_cldn4_explains_neigh_cd8": tumor_local_cd8,
        # convenience aliases used by RESULTS.md (prefer non-tumor juxta)
        **{f"primary_{k}": v for k, v in (nontumor or all_cells).items()},
    }


def plot_enrichment(z: pd.DataFrame, title: str, out: Path) -> None:
    show = [c for c in FOCUS if c in z.index]
    mat = z.loc[show, show]
    fig, ax = plt.subplots(figsize=(8.2, 6.8))
    vmax = np.nanmax(np.abs(mat.to_numpy()))
    vmax = 8 if not np.isfinite(vmax) or vmax < 8 else min(vmax, 20)
    sns.heatmap(
        mat,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        annot=True,
        fmt=".1f",
        linewidths=0.4,
        cbar_kws={"label": "Neighborhood enrichment z"},
    )
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_cooccurrence(df: pd.DataFrame, title: str, out: Path) -> None:
    if df.empty:
        return
    df = df.dropna(subset=["cooccurrence"]).copy()
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    for t, sub in df.groupby("target", sort=False):
        ax.plot(sub["d_mid_um"], sub["cooccurrence"], marker="o", label=t, lw=2)
    ax.axhline(1.0, color="0.4", ls="--", lw=1, label="random (score = 1)")
    ax.set_xlabel("Distance (µm)")
    ax.set_ylabel("Co-occurrence score  P(target | query, d) / P(target)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=11)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_composition(comp: pd.DataFrame, title: str, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 5.0))
    x = np.arange(len(comp))
    w = 0.38
    axes[0].bar(x - w / 2, comp["frac_around_CLDN4_high"], w, label="CLDN4-high tumor/epi", color="#b2182b")
    axes[0].bar(x + w / 2, comp["frac_around_CLDN4_low"], w, label="CLDN4-low tumor/epi", color="#2166ac")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(comp["neighbor_type"])
    axes[0].set_ylabel("Mean neighbor fraction")
    axes[0].set_title("Per-type composition")
    axes[0].legend(frameon=False, fontsize=9)
    pal = {"CD8": "#2166ac", "NK": "#67a9cf", "Treg": "#762a83", "Mac": "#4d9221"}
    bottoms = np.zeros(2)
    for _, row in comp.iterrows():
        vals = [row["frac_around_CLDN4_high"], row["frac_around_CLDN4_low"]]
        axes[1].bar([0, 1], vals, bottom=bottoms, color=pal.get(row["neighbor_type"], "0.5"), label=row["neighbor_type"])
        bottoms = bottoms + np.asarray(vals)
    axes[1].set_xticks([0, 1])
    axes[1].set_xticklabels(["CLDN4-high", "CLDN4-low"])
    axes[1].set_ylabel("Stacked mean neighbor fraction")
    axes[1].set_title("Stacked immune neighborhood")
    axes[1].legend(frameon=False, fontsize=9)
    colors = ["#b2182b" if p < 0.05 else "0.55" for p in comp["perm_p"]]
    axes[2].bar(x, comp["delta_high_minus_low"], color=colors)
    axes[2].axhline(0, color="0.3", lw=1)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(comp["neighbor_type"])
    axes[2].set_ylabel("Δ fraction (high − low)")
    axes[2].set_title("Permutation-tested difference")
    ymax = max(abs(float(comp["delta_high_minus_low"].min())), abs(float(comp["delta_high_minus_low"].max())), 1e-6)
    for i, p in enumerate(comp["perm_p"]):
        axes[2].text(i, 0.08 * ymax, f"p={p:.3f}", ha="center", fontsize=8)
    sns.despine(fig)
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_spatial_labels(adata: AnnData, title: str, out: Path, s: float = 2.0) -> None:
    xy = np.asarray(adata.obsm["spatial"])
    colors = {
        "CLDN4-high tumor/epi": "#b2182b",
        "CLDN4-low tumor/epi": "#f4a582",
        "CD8": "#2166ac",
        "NK": "#67a9cf",
        "Treg": "#762a83",
        "Mac": "#4d9221",
        "Other": "#d9d9d9",
    }
    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    for lab, col in colors.items():
        m = adata.obs["graph_label"].astype(str) == lab
        if not m.any():
            continue
        ax.scatter(xy[m, 0], xy[m, 1], s=s, c=col, label=lab, linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(title)
    ax.legend(markerscale=4, frameon=False, fontsize=9, loc="upper right")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_focus_z(z_cos: pd.DataFrame, z_vis: pd.DataFrame, out: Path) -> None:
    rows = []
    for name, z in [("CosMx Lung5_Rep2", z_cos), ("Visium GSE307534 P1_LUAD", z_vis)]:
        if "CLDN4-high tumor/epi" not in z.index:
            continue
        for t in IMMUNE:
            if t in z.columns:
                rows.append({"dataset": name, "pair": f"CLDN4-high vs {t}", "z": float(z.loc["CLDN4-high tumor/epi", t])})
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    sns.barplot(data=df, x="pair", y="z", hue="dataset", ax=ax)
    ax.axhline(0, color="0.3", lw=1)
    ax.set_ylabel("Neighborhood enrichment z")
    ax.set_xlabel("")
    ax.set_title("CLDN4-high tumor/epi vs immune neighborhoods")
    ax.legend(frameon=False)
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fmt_p(p: float) -> str:
    try:
        p = float(p)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.4f}"


def write_results(
    cos: dict,
    vis: dict,
    path: Path,
) -> None:
    zc, zv = cos["z"], vis["z"]
    cc, cv = cos["comp"], vis["comp"]
    mc, mv = cos["misty"], vis["misty"]

    def zf(z: pd.DataFrame, a: str, b: str) -> float:
        if a not in z.index or b not in z.columns:
            return float("nan")
        return float(z.loc[a, b])

    def zline(z: pd.DataFrame, a: str, b: str) -> str:
        v = zf(z, a, b)
        return "NA" if not np.isfinite(v) else f"{v:.2f}"

    def cline(comp: pd.DataFrame, t: str, col: str) -> str:
        row = comp.loc[comp["neighbor_type"] == t]
        if row.empty:
            return "NA"
        return f"{float(row.iloc[0][col]):.4f}"

    def misty_row(m: dict) -> dict:
        if "primary_n" in m:
            return {
                "n": m["primary_n"],
                "r2k": m["primary_r2_krt8_only"],
                "r2c": m["primary_r2_krt8_plus_cldn4"],
                "dr2": m["primary_delta_r2_cldn4_after_krt8"],
                "coef": m["primary_coef_neigh_CLDN4"],
                "p": m["primary_pval_neigh_CLDN4"],
            }
        j = m.get("juxta_nontumor") or m.get("juxta_all_cells") or m
        return {
            "n": j.get("n") or j.get("n_cells"),
            "r2k": j.get("r2_krt8_only", float("nan")),
            "r2c": j.get("r2_krt8_plus_cldn4", float("nan")),
            "dr2": j.get("delta_r2_cldn4_after_krt8", float("nan")),
            "coef": j.get("coef_neigh_CLDN4", float("nan")),
            "p": j.get("pval_neigh_CLDN4", float("nan")),
        }

    mcr, mvr = misty_row(mc), misty_row(mv)
    z_cd8 = zf(zc, "CLDN4-high tumor/epi", "CD8")
    z_nk = zf(zc, "CLDN4-high tumor/epi", "NK")
    z_mac = zf(zc, "CLDN4-high tumor/epi", "Mac")
    lead_bits = []
    lead_bits.append(
        f"On CosMx Lung5_Rep2 (50 µm), CLDN4-high tumor/epi vs CD8 neighborhood z = {z_cd8:.2f}"
        + (" (depleted vs random)" if z_cd8 < 0 else " (enriched vs random)")
        + f", vs NK z = {z_nk:.2f}, vs Mac z = {z_mac:.2f}."
    )
    lead_bits.append(
        f"Visium P1_LUAD 1–2 hop CLDN4-high vs CD8 z = {zf(zv, 'CLDN4-high tumor/epi', 'CD8'):.2f}."
    )
    if np.isfinite(mcr["dr2"]):
        lead_bits.append(
            f"After neighboring KRT8, neighboring CLDN4 adds ΔR² = {mcr['dr2']:.4f} to non-tumor CD8A on CosMx"
            f" (coef {mcr['coef']:.4f}, p={fmt_p(mcr['p'])})."
        )

    text = f"""# RESULTS: CLDN4-only spatial neighbor-graph analysis

{' '.join(lead_bits)}

## Datasets (public only)

| Assay | Accession / source | Sample | Resolution | Graph |
| --- | --- | --- | --- | --- |
| CosMx SMI NSCLC FFPE | Official NanoString / He et al. 2022 *Nat Biotechnol* public release | **Lung5_Rep2** | single cell | radius **50 µm** (0.18 µm/px from the official SMI-ReadMe for this prototype) |
| Visium CytAssist LUAD | **GSE307534** / GSM9226169 | **P1_LUAD invasive** section | spot (~55 µm) | Visium hex **1–2 hop** |

No private 8-KL material was used. TACSTD2/TROP2 was not used as a query. Cell labels are marker- and IF-informed (CosMx PanCK/CD45/CD3 plus RNA modules); they are not a published InSituType atlas re-import.

CosMx QC: drop `cell_ID==0`, `n_counts≥20`, `n_genes≥10`, area 80–25000 px. Tumor/epi cells were split by **median CLDN4 counts within tumor/epi**. Visium spots were library-size normalized to 10k and log1p-transformed; tumor spots = top 40% epithelial-score spots, then median-split on CLDN4.

## Cell / spot counts

**CosMx Lung5_Rep2** (n={cos['n']} cells)

{cos['counts_md']}

**Visium GSE307534 P1_LUAD** (n={vis['n']} in-tissue spots)

{vis['counts_md']}

## Neighborhood enrichment z-scores

Permutation test of spatial edges vs label shuffles (n={N_PERMS_ENRICH}). Positive z = more neighbors than expected; negative = avoidance.

### CosMx 50 µm (Lung5_Rep2)

CLDN4-high tumor/epi vs CD8 z = **{zline(zc, 'CLDN4-high tumor/epi', 'CD8')}**; vs NK **{zline(zc, 'CLDN4-high tumor/epi', 'NK')}**; vs Treg **{zline(zc, 'CLDN4-high tumor/epi', 'Treg')}**; vs Mac **{zline(zc, 'CLDN4-high tumor/epi', 'Mac')}**.

CLDN4-low tumor/epi vs CD8 z = **{zline(zc, 'CLDN4-low tumor/epi', 'CD8')}**.

![CosMx enrichment](results/figures/cosmx_nhood_enrichment.png)

### Visium 1–2 hop (P1_LUAD)

CLDN4-high tumor/epi vs CD8 z = **{zline(zv, 'CLDN4-high tumor/epi', 'CD8')}**; vs NK **{zline(zv, 'CLDN4-high tumor/epi', 'NK')}**; vs Treg **{zline(zv, 'CLDN4-high tumor/epi', 'Treg')}**; vs Mac **{zline(zv, 'CLDN4-high tumor/epi', 'Mac')}**.

![Visium enrichment](results/figures/visium_nhood_enrichment.png)

![Cross-platform focus](results/figures/cldn4_high_vs_immune_z.png)

## Co-occurrence vs distance

Score = P(target in annulus d | query) / P(target). Score > 1 = attraction at that scale.

![CosMx co-occurrence](results/figures/cosmx_cooccurrence_vs_distance.png)

![Visium co-occurrence](results/figures/visium_cooccurrence_vs_distance.png)

Short-range CosMx scores (0–50 µm) and Visium first occupied bin are tabulated in `*_cooccurrence.csv`. Visium graph is hop-based; the curve uses Euclidean µm from hex pitch (~100 µm center-to-center). Empty short bins are expected on Visium because spots are not 25 µm apart.

## Neighborhood composition (CLDN4-high vs low tumor) + permutation p

Mean fraction of CD8 / NK / Treg / Mac among spatial neighbors of CLDN4-high vs CLDN4-low **tumor/epi** cells. Labels of high vs low were shuffled among tumor cells ({N_PERMS_COMP} permutations).

### CosMx

| Neighbor | High | Low | Δ (high−low) | perm p |
| --- | ---: | ---: | ---: | ---: |
| CD8 | {cline(cc,'CD8','frac_around_CLDN4_high')} | {cline(cc,'CD8','frac_around_CLDN4_low')} | {cline(cc,'CD8','delta_high_minus_low')} | {cline(cc,'CD8','perm_p')} |
| NK | {cline(cc,'NK','frac_around_CLDN4_high')} | {cline(cc,'NK','frac_around_CLDN4_low')} | {cline(cc,'NK','delta_high_minus_low')} | {cline(cc,'NK','perm_p')} |
| Treg | {cline(cc,'Treg','frac_around_CLDN4_high')} | {cline(cc,'Treg','frac_around_CLDN4_low')} | {cline(cc,'Treg','delta_high_minus_low')} | {cline(cc,'Treg','perm_p')} |
| Mac | {cline(cc,'Mac','frac_around_CLDN4_high')} | {cline(cc,'Mac','frac_around_CLDN4_low')} | {cline(cc,'Mac','delta_high_minus_low')} | {cline(cc,'Mac','perm_p')} |

![CosMx composition](results/figures/cosmx_neighborhood_composition.png)

### Visium

| Neighbor | High | Low | Δ (high−low) | perm p |
| --- | ---: | ---: | ---: | ---: |
| CD8 | {cline(cv,'CD8','frac_around_CLDN4_high')} | {cline(cv,'CD8','frac_around_CLDN4_low')} | {cline(cv,'CD8','delta_high_minus_low')} | {cline(cv,'CD8','perm_p')} |
| NK | {cline(cv,'NK','frac_around_CLDN4_high')} | {cline(cv,'NK','frac_around_CLDN4_low')} | {cline(cv,'NK','delta_high_minus_low')} | {cline(cv,'NK','perm_p')} |
| Treg | {cline(cv,'Treg','frac_around_CLDN4_high')} | {cline(cv,'Treg','frac_around_CLDN4_low')} | {cline(cv,'Treg','delta_high_minus_low')} | {cline(cv,'Treg','perm_p')} |
| Mac | {cline(cv,'Mac','frac_around_CLDN4_high')} | {cline(cv,'Mac','frac_around_CLDN4_low')} | {cline(cv,'Mac','delta_high_minus_low')} | {cline(cv,'Mac','perm_p')} |

![Visium composition](results/figures/visium_neighborhood_composition.png)

## Spatial variance partitioning (MISTy-style)

OLS of log1p(CD8A) on **neighbor-mean KRT8**, then neighbor-mean KRT8 + **neighbor-mean CLDN4**. Incremental R² is the unique spatial contribution of neighboring CLDN4 after epithelial density (KRT8).

| Dataset | n | R² (KRT8 neigh) | R² (+ CLDN4 neigh) | ΔR² CLDN4 after KRT8 | CLDN4 coef | CLDN4 p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CosMx 50 µm (non-tumor) | {mcr['n']} | {mcr['r2k']:.4f} | {mcr['r2c']:.4f} | {mcr['dr2']:.4f} | {mcr['coef']:.4f} | {fmt_p(mcr['p'])} |
| Visium 1–2 hop (non-tumor) | {mvr['n']} | {mvr['r2k']:.4f} | {mvr['r2c']:.4f} | {mvr['dr2']:.4f} | {mvr['coef']:.4f} | {fmt_p(mvr['p'])} |

Primary model is juxta-view OLS on **non-tumor** cells: log1p(CD8A) ~ neighbor-mean KRT8 + neighbor-mean CLDN4. Sign of the CLDN4 coefficient is the direction of the residual association after KRT8. Tumor-centric and all-cell models are in `results/tables/misty_cd8_partition.json`. ΔR² can be significant and still tiny; report the coefficient with the increment.

## Spatial maps

![CosMx map](results/figures/cosmx_spatial_labels.png)

![Visium map](results/figures/visium_spatial_labels.png)

## Methods notes

- CosMx graph: `sklearn.neighbors.radius_neighbors_graph` on global µm coordinates, radius 50 µm, no self-loops.
- Visium graph: official hex 6-neighbors (`±2` columns on the same row; `±1` row/`±1` col), then 2-hop closure (`A ∨ A²`).
- Enrichment engine: `squidpy.gr.nhood_enrichment` when it completes; otherwise the same permutation z-score on undirected edges (implemented in this repo).
- Co-occurrence: BallTree annuli around a cap of {COOC_QUERY_CAP} CLDN4-high query cells.
- Composition p-values are two-sided label-shuffle tests among tumor/epi only (keeps tumor geography, breaks CLDN4 high/low).
- Squidpy {cos.get('squidpy_version', 'installed')} was available in this run.

## Interpretation (CLDN4-only)

Numbers above are the result. CosMx is the single-cell test: CLDN4-high epithelium is **self-clustered**, **CD8-poor vs random** (z={z_cd8:.2f}), and **macrophage-poor vs random** (z={z_mac:.2f}). NK is not treated as equivalent to CD8. Visium P1_LUAD labels are signature-dominant mixed spots; use hop-graph z-scores and composition p-values, not cell-pure calls. Composition tests ask a different question than enrichment (high vs low tumor, not vs random).

## Files

- `results/tables/cosmx_nhood_zscore.csv`, `visium_nhood_zscore.csv`
- `results/tables/cosmx_composition.csv`, `visium_composition.csv`
- `results/tables/cosmx_cooccurrence.csv`, `visium_cooccurrence.csv`
- `results/tables/misty_cd8_partition.json`
- `results/figures/*.png`
"""
    path.write_text(text)


def counts_markdown(adata: AnnData) -> str:
    vc = adata.obs["graph_label"].value_counts().reindex(LABEL_ORDER).fillna(0).astype(int)
    lines = ["| Label | n |", "| --- | ---: |"]
    for k, v in vc.items():
        lines.append(f"| {k} | {v} |")
    return "\n".join(lines)


def run_one(name: str, adata: AnnData, graph: str) -> dict:
    LOG.info("=== %s graph=%s ===", name, graph)
    if graph == "radius50":
        build_radius_graph(adata, RADIUS_UM)
        # keep µm for co-occurrence
        xy_um = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    else:
        build_visium_hop_graph(adata, hops=2)
        # estimate µm from hex pitch ~100 µm
        xy_um = np.column_stack(
            [
                adata.obs["array_col"].to_numpy() * 100.0,
                adata.obs["array_row"].to_numpy() * np.sqrt(3) / 2.0 * 100.0,
            ]
        )
        adata.obsm["spatial_um_grid"] = xy_um
        # co-occurrence uses spatial; temporarily swap to µm-grid
        adata.obsm["spatial_px"] = adata.obsm["spatial"]
        adata.obsm["spatial"] = xy_um

    z = try_squidpy_enrichment(adata)
    if z is None:
        z, _ = nhood_enrichment_z(adata)
    z.to_csv(OUT / "tables" / f"{name}_nhood_zscore.csv")
    comp = composition_test(adata)
    comp.to_csv(OUT / "tables" / f"{name}_composition.csv", index=False)
    bins = VISIUM_COOC_BINS_UM if graph == "hop2" else COOC_BINS_UM
    cooc = cooccurrence_vs_distance(adata, bins=bins)
    cooc.to_csv(OUT / "tables" / f"{name}_cooccurrence.csv", index=False)
    misty = misty_like_cd8(adata)
    plot_enrichment(
        z,
        f"{name}: neighborhood enrichment z",
        OUT / "figures" / f"{name}_nhood_enrichment.png",
    )
    plot_cooccurrence(
        cooc,
        f"{name}: CLDN4-high tumor/epi co-occurrence vs distance",
        OUT / "figures" / f"{name}_cooccurrence_vs_distance.png",
    )
    plot_composition(
        comp,
        f"{name}: CLDN4-high vs low tumor neighborhoods",
        OUT / "figures" / f"{name}_neighborhood_composition.png",
    )
    plot_spatial_labels(
        adata,
        f"{name}: graph labels",
        OUT / "figures" / f"{name}_spatial_labels.png",
        s=1.4 if name == "cosmx" else 10,
    )
    return {
        "z": z,
        "comp": comp,
        "cooc": cooc,
        "misty": misty,
        "n": int(adata.n_obs),
        "counts_md": counts_markdown(adata),
        "label_counts": adata.uns["label_counts"],
    }


def main() -> None:
    setup()
    try:
        import squidpy as sq

        sq_ver = getattr(sq, "__version__", "1.x")
    except Exception:
        sq_ver = "unavailable"

    adata_c = load_cosmx(COSMX_DIR)
    adata_v = load_visium(VISIUM_DIR)
    cos = run_one("cosmx", adata_c, "radius50")
    vis = run_one("visium", adata_v, "hop2")
    cos["squidpy_version"] = sq_ver
    vis["squidpy_version"] = sq_ver
    plot_focus_z(cos["z"], vis["z"], OUT / "figures" / "cldn4_high_vs_immune_z.png")
    (OUT / "tables" / "misty_cd8_partition.json").write_text(
        json.dumps({"cosmx": cos["misty"], "visium": vis["misty"]}, indent=2)
    )
    write_results(cos, vis, Path("/workspace/RESULTS.md"))
    LOG.info("Wrote RESULTS.md and figures")


if __name__ == "__main__":
    main()
