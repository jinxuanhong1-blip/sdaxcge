#!/usr/bin/env python3
"""Spatial ligand-receptor CCC: CLDN4-high vs CLDN4-low tumor niches vs CD8.

Official CosMx SMI NSCLC 960-plex (He et al. 2022; NanoString public S3;
same 8 slides as figshare 25976224). CLDN4-only. No private 8-KL. No TACSTD2 gate.

Three installable tools, same cell definitions:

- squidpy.gr.ligrec: CellPhoneDB permutation on cells within 40 um of the
  other compartment. The high-vs-low contrast with a shared CD8 receiver is a
  ligand-abundance contrast (the receptor mean cancels). It is not distance-weighted.
- SpatialDM global bivariate Moran I / z, single-cell RBF, with ligand kept on
  the sender class and receptor kept on the receiver class, and weights kept
  only on sender-to-receiver edges.
- COMMOT collective optimal transport (dis_thr=50 um) on proximity-conditioned
  cells. The reported score is transport into the opposite class, per sender.

Primary inference is the sign of (CLDN4-high - CLDN4-low) across the 8 slides,
not a cell-level permutation p-value.

F11R and NECTIN2 are not on this 960-plex and are not tested.
"""

from __future__ import annotations

import argparse
import json
import time
import zlib
from pathlib import Path

import numpy as np

# commot 0.0.3 evaluates np.Inf at import. NumPy 2 removed that alias.
if not hasattr(np, "Inf"):
    np.Inf = np.inf

import anndata as ad
import commot as ct
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import spatialdm as sdm
import squidpy as sq
from scipy import sparse, stats
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cosmx"
OUT = ROOT / "results" / "cosmx_cldn4_spatial_ccc"

SAMPLES = [
    "Lung5_Rep1",
    "Lung5_Rep2",
    "Lung5_Rep3",
    "Lung6",
    "Lung9_Rep1",
    "Lung9_Rep2",
    "Lung12",
    "Lung13",
]
PATIENT = {
    "Lung5_Rep1": "Lung5",
    "Lung5_Rep2": "Lung5",
    "Lung5_Rep3": "Lung5",
    "Lung6": "Lung6",
    "Lung9_Rep1": "Lung9",
    "Lung9_Rep2": "Lung9",
    "Lung12": "Lung12",
    "Lung13": "Lung13",
}

PX_TO_UM = 0.18
MIN_COUNTS = 20
MIN_GENES = 5
PROX_UM = 40.0
COMMOT_DIS_UM = 50.0
SDM_L_UM = 30.0
SDM_CUTOFF = 0.15
SDM_N_NEIGHBORS = 80
SDM_N_NEAREST = 10
MIN_SENDER = 8
MIN_RECEIVER = 5
GEOM_CAP_SENDER = 4000
GEOM_CAP_RECEIVER = 1500
COMMOT_CAP = {"tumor_hi": 140, "tumor_lo": 140, "cd8": 120}
COMMOT_NITER = 1000
SEED = 20260921
CODE_STAMP = "v1-prox40-sdm30-commot50-niter1000"
# TIGIT itself is on the 960-plex; its CellChat ligands (PVR, NECTIN2) and CD226 are not.
ABSENT_GENES = ["F11R", "NECTIN2", "PVR", "SIRPA", "CD226"]

# ligand, receptor, pathway, annotation, family, direction
# direction tumor_to_cd8: ligand on tumor, receptor on CD8
# direction cd8_to_tumor: ligand on CD8, receptor on tumor
PAIRS = [
    ("CDH1", "CDH1", "CDH", "Cell-Cell Contact", "barrier", "tumor_to_cd8"),
    ("CDH1", "ITGA2_ITGB1", "CDH1", "Cell-Cell Contact", "barrier", "tumor_to_cd8"),
    ("ICAM1", "ITGAL_ITGB2", "ICAM", "Cell-Cell Contact", "barrier", "tumor_to_cd8"),
    ("ICAM1", "ITGAL", "ICAM", "Cell-Cell Contact", "barrier", "tumor_to_cd8"),
    ("CD274", "PDCD1", "PD-L1", "Cell-Cell Contact", "barrier", "tumor_to_cd8"),
    ("PDCD1LG2", "PDCD1", "PDL2", "Cell-Cell Contact", "barrier", "tumor_to_cd8"),
    ("LGALS9", "HAVCR2", "GALECTIN", "Secreted Signaling", "barrier", "tumor_to_cd8"),
    ("LGALS9", "PTPRC", "GALECTIN", "Secreted Signaling", "barrier", "tumor_to_cd8"),
    ("LGALS9", "CD44", "GALECTIN", "Secreted Signaling", "barrier", "tumor_to_cd8"),
    ("MIF", "CD74_CD44", "MIF", "Secreted Signaling", "barrier", "tumor_to_cd8"),
    ("MIF", "CD74_CXCR4", "MIF", "Secreted Signaling", "barrier", "tumor_to_cd8"),
    ("TGFB1", "TGFBR1_TGFBR2", "TGFb", "Secreted Signaling", "barrier", "tumor_to_cd8"),
    ("CXCL12", "CXCR4", "CXCL", "Secreted Signaling", "chemokine", "tumor_to_cd8"),
    ("CXCL16", "CXCR6", "CXCL", "Secreted Signaling", "chemokine", "tumor_to_cd8"),
    ("CXCL9", "CXCR3", "CXCL", "Secreted Signaling", "chemokine", "tumor_to_cd8"),
    ("CXCL10", "CXCR3", "CXCL", "Secreted Signaling", "chemokine", "tumor_to_cd8"),
    ("CCL5", "CCR5", "CCL", "Secreted Signaling", "chemokine", "tumor_to_cd8"),
    ("HLA-A", "CD8A", "MHC-I", "Cell-Cell Contact", "mhc_control", "tumor_to_cd8"),
    ("HLA-B", "CD8A", "MHC-I", "Cell-Cell Contact", "mhc_control", "tumor_to_cd8"),
    ("HLA-C", "CD8A", "MHC-I", "Cell-Cell Contact", "mhc_control", "tumor_to_cd8"),
    ("HLA-A", "CD8B", "MHC-I", "Cell-Cell Contact", "mhc_control", "tumor_to_cd8"),
    ("HLA-B", "CD8B", "MHC-I", "Cell-Cell Contact", "mhc_control", "tumor_to_cd8"),
    ("HLA-C", "CD8B", "MHC-I", "Cell-Cell Contact", "mhc_control", "tumor_to_cd8"),
    ("APP", "CD74", "APP", "Cell-Cell Contact", "other", "tumor_to_cd8"),
    ("FASLG", "FAS", "FASLG", "Secreted Signaling", "other", "tumor_to_cd8"),
    ("SPP1", "CD44", "SPP1", "Secreted Signaling", "other", "tumor_to_cd8"),
    ("IFNG", "IFNGR1_IFNGR2", "IFN-II", "Secreted Signaling", "effector_reverse", "cd8_to_tumor"),
]

ANN_RANK = {"Cell-Cell Contact": 0, "ECM-Receptor": 1, "Secreted Signaling": 2}

EPI_MARKERS = [
    "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT17", "KRT15",
    "KRT14", "KRT13", "KRT16", "KRT6A", "KRT5", "KRT7", "CEACAM6",
    "MUC1", "ELF3", "SFTPC", "SFTPB", "NAPSA", "S100A14",
]
IMM_MARKERS = [
    "PTPRC", "CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "CD4",
    "IL7R", "CCL5", "GZMA", "GZMB", "GZMK", "GNLY", "KLRB1", "KLRK1",
    "CD68", "CD163", "C1QA", "C1QB", "C1QC", "LYZ", "FCGR3A", "ITGAX",
    "ITGAM", "CD14", "CD79A", "MS4A1", "CD19", "CD37", "CD52", "CD53",
]
STR_MARKERS = [
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "DCN",
    "LUM", "FN1", "BGN", "PDGFRB", "ACTA2", "MYH11", "PECAM1", "VWF",
    "CDH5", "CLEC14A", "FLT1", "KDR", "ESAM", "RAMP2",
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def pair_frame() -> pd.DataFrame:
    rows = []
    for lig, rec, path, ann, family, direction in PAIRS:
        rows.append(
            {
                "pair_id": f"{lig}|{rec}",
                "ligand": lig,
                "receptor": rec,
                "pathway": path,
                "annotation": ann,
                "family": family,
                "direction": direction,
            }
        )
    return pd.DataFrame(rows)


def subunits(name: str) -> list[str]:
    return [g for g in str(name).split("_") if g]


def atomic_name(name: str) -> str:
    """Squidpy splits complexes on '_'. Keep the per-cell minimum under a '+' name."""
    parts = subunits(name)
    return parts[0] if len(parts) == 1 else "+".join(parts)


def find_file(sample: str, kind: str, data_dir: Path) -> Path:
    hits = list(data_dir.rglob(f"{sample}_{kind}_file.csv"))
    if not hits:
        raise FileNotFoundError(f"{sample} {kind} not under {data_dir}")
    return hits[0]


def read_panel(expr_path: Path) -> list[str]:
    header = pd.read_csv(expr_path, nrows=0).columns.tolist()
    return [c for c in header if c not in ("fov", "cell_ID") and not str(c).lower().startswith("neg")]


def confirm_pairs_in_cellchat(pairs: pd.DataFrame) -> pd.DataFrame:
    db = ct.pp.ligand_receptor_database(database="CellChat", species="human", signaling_type=None)
    db = db.copy()
    db.columns = ["ligand", "receptor", "pathway", "annotation"][: db.shape[1]]
    keys = set(zip(db["ligand"].astype(str), db["receptor"].astype(str)))
    out = pairs.copy()
    out["in_cellchat_human"] = [
        (r.ligand, r.receptor) in keys for r in out.itertuples(index=False)
    ]
    return out


def gene_universe(pairs: pd.DataFrame) -> list[str]:
    genes = set(EPI_MARKERS + IMM_MARKERS + STR_MARKERS + ["CLDN4", "CD8A", "CD8B", "CD3D", "CD3E", "CD3G"])
    for r in pairs.itertuples(index=False):
        genes.update(subunits(r.ligand))
        genes.update(subunits(r.receptor))
    return sorted(genes)


def rng_for(sample: str, fov: int) -> np.random.Generator:
    salt = zlib.adler32(f"{sample}:{fov}:{CODE_STAMP}".encode()) % 1_000_003
    return np.random.default_rng(SEED + salt)


def load_sample(sample: str, data_dir: Path, keep_genes: list[str], panel: list[str]) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    expr_path = find_file(sample, "exprMat", data_dir)
    meta_path = find_file(sample, "metadata", data_dir)
    use_genes = [g for g in keep_genes if g in set(panel)]
    pieces = []
    ncount_parts = []
    ngene_parts = []
    log(f"{sample}: reading expression")
    for chunk in pd.read_csv(expr_path, chunksize=20000):
        chunk = chunk.loc[chunk["cell_ID"] != 0]
        if chunk.empty:
            continue
        mat = chunk.loc[:, panel].to_numpy(dtype=np.float32, copy=False)
        ncount_parts.append(mat.sum(axis=1))
        ngene_parts.append((mat > 0).sum(axis=1))
        pieces.append(chunk.loc[:, ["fov", "cell_ID", *use_genes]])
        del mat
    expr = pd.concat(pieces, ignore_index=True)
    expr["nCount"] = np.concatenate(ncount_parts)
    expr["nGene"] = np.concatenate(ngene_parts)
    del pieces, ncount_parts, ngene_parts

    meta = pd.read_csv(meta_path)
    meta = meta.loc[meta["cell_ID"] != 0, ["fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"]].copy()
    expr["fov"] = expr["fov"].astype(int)
    expr["cell_ID"] = expr["cell_ID"].astype(int)
    meta["fov"] = meta["fov"].astype(int)
    meta["cell_ID"] = meta["cell_ID"].astype(int)
    df = expr.merge(meta, on=["fov", "cell_ID"], how="inner")
    df = df.loc[(df["nCount"] >= MIN_COUNTS) & (df["nGene"] >= MIN_GENES)].reset_index(drop=True)
    raw = df.loc[:, use_genes].to_numpy(dtype=np.float32, copy=True)
    tot = df["nCount"].to_numpy(dtype=np.float64)
    med = float(np.median(tot)) if len(tot) else 1.0
    scale = med / np.maximum(tot, 1.0)
    X = np.log1p(raw * scale[:, None]).astype(np.float32)

    gmap = {g: i for i, g in enumerate(use_genes)}

    def score(markers: list[str]) -> np.ndarray:
        idx = [gmap[g] for g in markers if g in gmap]
        if not idx:
            return np.zeros(len(df), dtype=np.float32)
        return X[:, idx].mean(axis=1)

    epi, imm, stro = score(EPI_MARKERS), score(IMM_MARKERS), score(STR_MARKERS)
    S = np.vstack([epi, imm, stro])
    lab = np.array(["epithelial", "immune", "stromal"])[S.argmax(axis=0)]
    lab[S.max(axis=0) <= 0] = "unassigned"

    def raw_col(name: str) -> np.ndarray:
        if name not in gmap:
            return np.zeros(len(df), dtype=np.float32)
        return raw[:, gmap[name]]

    out = pd.DataFrame(
        {
            "sample": sample,
            "patient": PATIENT[sample],
            "fov": df["fov"].to_numpy(),
            "cell_ID": df["cell_ID"].to_numpy(),
            "x_um": df["CenterX_global_px"].to_numpy(dtype=float) * PX_TO_UM,
            "y_um": df["CenterY_global_px"].to_numpy(dtype=float) * PX_TO_UM,
            "nCount": df["nCount"].to_numpy(dtype=float),
            "nGene": df["nGene"].to_numpy(dtype=float),
            "compartment": lab,
        }
    )
    out["is_tumor"] = out["compartment"].eq("epithelial")
    cd8_raw = raw_col("CD8A") + raw_col("CD8B")
    cd3_raw = raw_col("CD3D") + raw_col("CD3E") + raw_col("CD3G")
    out["is_cd8"] = (cd8_raw > 0) & (cd3_raw > 0) & (~out["is_tumor"].to_numpy())
    out["cldn4"] = X[:, gmap["CLDN4"]] if "CLDN4" in gmap else 0.0
    out["cldn4_tertile"] = 0
    log(f"{sample}: QC cells={len(out)} tumor={int(out.is_tumor.sum())} cd8={int(out.is_cd8.sum())}")
    return out, X, use_genes


def tertile_split(values: np.ndarray) -> np.ndarray:
    lab = np.zeros(len(values), dtype=int)
    if len(values) < 6:
        return lab
    ranks = pd.Series(values).rank(method="first").to_numpy()
    q1, q2 = np.quantile(ranks, [1 / 3, 2 / 3])
    lab[ranks <= q1] = -1
    lab[ranks >= q2] = 1
    return lab


def add_classes_and_distances(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["nn_cd8_um"] = np.nan
    df["nn_tumor_um"] = np.nan
    df["nn_tumor_hi_um"] = np.nan
    df["nn_tumor_lo_um"] = np.nan
    for fov, idx in df.groupby("fov").groups.items():
        pos = df.index.get_indexer(idx)
        tumor = df.loc[idx, "is_tumor"].to_numpy()
        if tumor.sum() >= 6:
            labs = tertile_split(df.loc[idx, "cldn4"].to_numpy()[tumor])
            tpos = pos[tumor]
            df.loc[df.index[tpos], "cldn4_tertile"] = labs
        xy = df.loc[idx, ["x_um", "y_um"]].to_numpy(float)
        cd8 = df.loc[idx, "is_cd8"].to_numpy()
        hi = tumor & (df.loc[idx, "cldn4_tertile"].to_numpy() == 1)
        lo = tumor & (df.loc[idx, "cldn4_tertile"].to_numpy() == -1)
        if cd8.any() and tumor.any():
            d_t, _ = cKDTree(xy[cd8]).query(xy[tumor], k=1)
            df.loc[df.index[pos[tumor]], "nn_cd8_um"] = d_t
            d_c, _ = cKDTree(xy[tumor]).query(xy[cd8], k=1)
            df.loc[df.index[pos[cd8]], "nn_tumor_um"] = d_c
        if cd8.any() and hi.any():
            d_c, _ = cKDTree(xy[hi]).query(xy[cd8], k=1)
            df.loc[df.index[pos[cd8]], "nn_tumor_hi_um"] = d_c
        if cd8.any() and lo.any():
            d_c, _ = cKDTree(xy[lo]).query(xy[cd8], k=1)
            df.loc[df.index[pos[cd8]], "nn_tumor_lo_um"] = d_c
    return df


def per_cell_min(X: np.ndarray, gmap: dict[str, int], name: str) -> np.ndarray | None:
    genes = subunits(name)
    if any(g not in gmap for g in genes):
        return None
    cols = np.column_stack([X[:, gmap[g]] for g in genes])
    return cols.min(axis=1).astype(np.float64)


def cap_indices(mask_s: np.ndarray, mask_r: np.ndarray, cap_s: int, cap_r: int, rng: np.random.Generator) -> np.ndarray:
    s = np.flatnonzero(mask_s)
    r = np.flatnonzero(mask_r)
    if len(s) > cap_s:
        s = rng.choice(s, size=cap_s, replace=False)
    if len(r) > cap_r:
        r = rng.choice(r, size=cap_r, replace=False)
    if len(s) == 0 or len(r) == 0:
        return np.array([], dtype=int)
    return np.sort(np.concatenate([s, r]))


def _mask_weight(W, receiver: np.ndarray, sender: np.ndarray):
    W = W.tolil(copy=True)
    drop_r = np.flatnonzero(~receiver)
    drop_s = np.flatnonzero(~sender)
    if len(drop_r):
        W[drop_r, :] = 0
    if len(drop_s):
        W[:, drop_s] = 0
    W = W.tocsr()
    W.eliminate_zeros()
    total = float(W.sum())
    if total <= 0:
        return None
    return (W * (W.shape[0] / total)).tocsr()


def spatialdm_directional(
    xy: np.ndarray,
    sender: np.ndarray,
    receiver: np.ndarray,
    lig: np.ndarray,
    rec: np.ndarray,
    annotations: list[str],
    pair_ids: list[str],
) -> dict[str, tuple[float, float]]:
    """Return pair_id -> (Moran I, z). I=0 when the arm exists but no edge or no expression."""
    nan = {pid: (np.nan, np.nan) for pid in pair_ids}
    n = xy.shape[0]
    if int(sender.sum()) < MIN_SENDER or int(receiver.sum()) < MIN_RECEIVER or n < 12:
        return nan
    order = np.argsort([ANN_RANK.get(a, 2) for a in annotations], kind="mergesort")
    ids = [pair_ids[i] for i in order]
    anns = [annotations[i] for i in order]
    L = np.asarray(lig[:, order], dtype=np.float64).copy()
    R = np.asarray(rec[:, order], dtype=np.float64).copy()
    L[~sender, :] = 0.0
    R[~receiver, :] = 0.0
    keep = []
    absent = []
    for j in range(L.shape[1]):
        if L[:, j].sum() > 0 and R[:, j].sum() > 0 and L[:, j].std() > 1e-8 and R[:, j].std() > 1e-8:
            keep.append(j)
        else:
            absent.append(j)
    out: dict[str, tuple[float, float]] = {ids[j]: (0.0, np.nan) for j in absent}
    if not keep:
        return out
    k = min(SDM_N_NEIGHBORS, n - 1)
    k0 = min(SDM_N_NEAREST, n - 1)
    if k < 2 or k0 < 1:
        return {pid: (0.0, np.nan) for pid in pair_ids}
    L = L[:, keep]
    R = R[:, keep]
    anns_k = [anns[j] for j in keep]
    ids_k = [ids[j] for j in keep]
    # Rebuild so non-secreted pairs stay in front after the expression filter.
    reorder = np.argsort([ANN_RANK.get(a, 2) for a in anns_k], kind="mergesort")
    L, R = L[:, reorder], R[:, reorder]
    anns_k = [anns_k[j] for j in reorder]
    ids_k = [ids_k[j] for j in reorder]
    var_names = [f"Lg{j}" for j in range(len(ids_k))] + [f"Rg{j}" for j in range(len(ids_k))]
    X = np.hstack([L, R])
    adata = ad.AnnData(X=X)
    adata.var_names = var_names
    adata.obsm["spatial"] = np.asarray(xy, dtype=np.float64)
    adata.uns["mean"] = "algebra"
    adata.uns["ligand"] = pd.DataFrame({"Ligand0": [f"Lg{j}" for j in range(len(ids_k))]}, index=ids_k)
    adata.uns["receptor"] = pd.DataFrame({"Receptor0": [f"Rg{j}" for j in range(len(ids_k))]}, index=ids_k)
    adata.uns["geneInter"] = pd.DataFrame({"annotation": anns_k}, index=ids_k)
    sdm.weight_matrix(
        adata,
        l=SDM_L_UM,
        cutoff=SDM_CUTOFF,
        n_neighbors=k,
        n_nearest_neighbors=k0,
        single_cell=True,
    )
    Ww = _mask_weight(adata.obsp["weight"], receiver, sender)
    Wk = _mask_weight(adata.obsp["nearest_neighbors"], receiver, sender)
    if Ww is None and Wk is None:
        for pid in ids_k:
            out[pid] = (0.0, np.nan)
        return out
    n_obs = adata.n_obs
    adata.obsp["weight"] = Ww if Ww is not None else sparse.csr_matrix((n_obs, n_obs))
    adata.obsp["nearest_neighbors"] = Wk if Wk is not None else sparse.csr_matrix((n_obs, n_obs))
    try:
        sdm.spatialdm_global(adata, n_perm=0, method="z-score", nproc=1)
    except Exception as exc:  # noqa: BLE001 — record the tool failure, do not invent a score
        log(f"SpatialDM failed: {exc}")
        for pid in ids_k:
            out[pid] = (np.nan, np.nan)
        return out
    got_ids = list(adata.uns["global_res"].index)
    I = np.asarray(adata.uns["global_I"], dtype=float).reshape(-1)
    z = np.asarray(adata.uns["global_stat"]["z"]["z"], dtype=float).reshape(-1)
    for pid, Ii, zi in zip(got_ids, I, z):
        if not np.isfinite(Ii):
            Ii = np.nan
        if not np.isfinite(zi):
            zi = np.nan
        out[str(pid)] = (float(Ii), float(zi))
    return out


def expr_block(X: np.ndarray, gmap: dict[str, int], pairs: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, list[str]]:
    ligs = []
    recs = []
    ids = []
    for r in pairs.itertuples(index=False):
        L = per_cell_min(X, gmap, r.ligand)
        R = per_cell_min(X, gmap, r.receptor)
        if L is None or R is None:
            continue
        ligs.append(L)
        recs.append(R)
        ids.append(r.pair_id)
    if not ids:
        return np.zeros((len(X), 0)), np.zeros((len(X), 0)), []
    return np.column_stack(ligs), np.column_stack(recs), ids


def score_fov_spatial(
    xy: np.ndarray,
    labels: np.ndarray,
    lig: np.ndarray,
    rec: np.ndarray,
    pairs_used: pd.DataFrame,
    rng: np.random.Generator,
) -> dict[str, dict[str, float]]:
    """labels in {tumor_hi, tumor_lo, cd8, other}. lig/rec columns match pairs_used order."""
    pair_ids = pairs_used["pair_id"].tolist()
    anns = pairs_used["annotation"].tolist()
    directions = pairs_used["direction"].tolist()
    out = {
        pid: {
            "I_hi_geom": np.nan,
            "z_hi_geom": np.nan,
            "I_lo_geom": np.nan,
            "z_lo_geom": np.nan,
            "I_hi_prox": np.nan,
            "z_hi_prox": np.nan,
            "I_lo_prox": np.nan,
            "z_lo_prox": np.nan,
        }
        for pid in pair_ids
    }
    is_hi = labels == "tumor_hi"
    is_lo = labels == "tumor_lo"
    is_cd8 = labels == "cd8"
    if lig.shape[1] == 0:
        return out

    # Split forward and reverse so the sender/receiver masks match the pair.
    for direction, sender_hi, sender_lo, recv_name in (
        ("tumor_to_cd8", is_hi, is_lo, "cd8"),
        ("cd8_to_tumor", is_cd8, is_cd8, "tumor"),
    ):
        sel = [i for i, d in enumerate(directions) if d == direction]
        if not sel:
            continue
        sub_ids = [pair_ids[i] for i in sel]
        sub_ann = [anns[i] for i in sel]
        sub_lig = lig[:, sel]
        sub_rec = rec[:, sel]
        if direction == "tumor_to_cd8":
            arms = (("hi", sender_hi, is_cd8), ("lo", sender_lo, is_cd8))
        else:
            arms = (("hi", is_cd8, is_hi), ("lo", is_cd8, is_lo))
        for arm, sender_mask, receiver_mask in arms:
            idx = cap_indices(sender_mask, receiver_mask, GEOM_CAP_SENDER, GEOM_CAP_RECEIVER, rng)
            if len(idx) == 0:
                continue
            scored = spatialdm_directional(
                xy[idx],
                sender_mask[idx],
                receiver_mask[idx],
                sub_lig[idx],
                sub_rec[idx],
                sub_ann,
                sub_ids,
            )
            for pid, (Ii, zi) in scored.items():
                out[pid][f"I_{arm}_geom"] = Ii
                out[pid][f"z_{arm}_geom"] = zi
            # Proximity-conditioned: same masks, but only cells already within PROX_UM
            # are passed in `labels` for the prox call. The caller uses a full-FOV label
            # vector here; proximity is applied by the caller via a second invocation
            # only when mode requests it. Geometry is this function's first pass.
            _ = recv_name
    return out


def score_modes_for_fov(
    df_fov: pd.DataFrame,
    lig: np.ndarray,
    rec: np.ndarray,
    pairs_used: pd.DataFrame,
    rng: np.random.Generator,
) -> dict[str, dict[str, float]]:
    labels = np.full(len(df_fov), "other", dtype=object)
    tert = df_fov["cldn4_tertile"].to_numpy()
    tumor = df_fov["is_tumor"].to_numpy()
    cd8 = df_fov["is_cd8"].to_numpy()
    labels[tumor & (tert == 1)] = "tumor_hi"
    labels[tumor & (tert == -1)] = "tumor_lo"
    labels[cd8] = "cd8"
    xy = df_fov[["x_um", "y_um"]].to_numpy(float)
    pair_ids = pairs_used["pair_id"].tolist()
    base = {
        pid: {
            "I_hi_geom": np.nan,
            "z_hi_geom": np.nan,
            "I_lo_geom": np.nan,
            "z_lo_geom": np.nan,
            "I_hi_prox": np.nan,
            "z_hi_prox": np.nan,
            "I_lo_prox": np.nan,
            "z_lo_prox": np.nan,
            "expr_hi": np.nan,
            "expr_lo": np.nan,
        }
        for pid in pair_ids
    }
    directions = pairs_used["direction"].tolist()
    anns = pairs_used["annotation"].tolist()

    def run_mode(mode: str, cell_idx: np.ndarray, sender_hi, sender_lo, receiver_for_forward_is_cd8=True):
        if len(cell_idx) < 12:
            return
        # sender/receiver masks are in the FULL fov index space; slice after.
        for direction in ("tumor_to_cd8", "cd8_to_tumor"):
            sel = [i for i, d in enumerate(directions) if d == direction]
            if not sel:
                continue
            sub_ids = [pair_ids[i] for i in sel]
            sub_ann = [anns[i] for i in sel]
            sub_lig = lig[:, sel]
            sub_rec = rec[:, sel]
            if direction == "tumor_to_cd8":
                arms = (("hi", sender_hi, cd8), ("lo", sender_lo, cd8))
            else:
                arms = (("hi", cd8, sender_hi), ("lo", cd8, sender_lo))
            for arm, smask, rmask in arms:
                both = smask | rmask
                idx = np.intersect1d(cell_idx, np.flatnonzero(both), assume_unique=False)
                if len(idx) < 12:
                    continue
                # Cap inside the selected idx by rebuilding masks on the subset via cap_indices
                # in subset coordinates.
                local_s = smask[idx]
                local_r = rmask[idx]
                take_local = cap_indices(local_s, local_r, GEOM_CAP_SENDER, GEOM_CAP_RECEIVER, rng)
                if len(take_local) < 12:
                    continue
                scored = spatialdm_directional(
                    xy[idx][take_local],
                    local_s[take_local],
                    local_r[take_local],
                    sub_lig[idx][take_local],
                    sub_rec[idx][take_local],
                    sub_ann,
                    sub_ids,
                )
                key_I = f"I_{arm}_{mode}"
                key_z = f"z_{arm}_{mode}"
                for pid, (Ii, zi) in scored.items():
                    base[pid][key_I] = Ii
                    base[pid][key_z] = zi

    hi = labels == "tumor_hi"
    lo = labels == "tumor_lo"
    geom_idx = np.flatnonzero(hi | lo | cd8)
    run_mode("geom", geom_idx, hi, lo)
    prox_tumor = tumor & np.isfinite(df_fov["nn_cd8_um"].to_numpy()) & (df_fov["nn_cd8_um"].to_numpy() <= PROX_UM)
    prox_cd8 = cd8 & np.isfinite(df_fov["nn_tumor_um"].to_numpy()) & (df_fov["nn_tumor_um"].to_numpy() <= PROX_UM)
    # Class-specific proximity for the arms: tumor within 40 um of any CD8,
    # CD8 within 40 um of that tumor class.
    hi_prox = hi & prox_tumor
    lo_prox = lo & prox_tumor
    cd8_hi = cd8 & np.isfinite(df_fov["nn_tumor_hi_um"].to_numpy()) & (df_fov["nn_tumor_hi_um"].to_numpy() <= PROX_UM)
    cd8_lo = cd8 & np.isfinite(df_fov["nn_tumor_lo_um"].to_numpy()) & (df_fov["nn_tumor_lo_um"].to_numpy() <= PROX_UM)
    # For prox mode the receiver set differs by arm, so call spatialdm per arm
    # rather than the shared run_mode receiver. Reuse run_mode with arm-specific
    # CD8 masks by passing cd8_hi as the receiver stand-in only for the hi arm.
    # Easiest correct path: two idx sets.
    run_mode_prox_split(base, xy, lig, rec, pair_ids, directions, anns, hi_prox, lo_prox, cd8_hi, cd8_lo, rng)

    # Expression delta among proximal senders. Forward: tumor ligand. Reverse: CD8 ligand.
    for j, pid in enumerate(pair_ids):
        if directions[j] == "tumor_to_cd8":
            if hi_prox.sum() >= MIN_SENDER:
                base[pid]["expr_hi"] = float(lig[hi_prox, j].mean())
            if lo_prox.sum() >= MIN_SENDER:
                base[pid]["expr_lo"] = float(lig[lo_prox, j].mean())
        else:
            if cd8_hi.sum() >= MIN_RECEIVER:
                base[pid]["expr_hi"] = float(lig[cd8_hi, j].mean())
            if cd8_lo.sum() >= MIN_RECEIVER:
                base[pid]["expr_lo"] = float(lig[cd8_lo, j].mean())
    return base


def run_mode_prox_split(base, xy, lig, rec, pair_ids, directions, anns, hi_prox, lo_prox, cd8_hi, cd8_lo, rng):
    for direction in ("tumor_to_cd8", "cd8_to_tumor"):
        sel = [i for i, d in enumerate(directions) if d == direction]
        if not sel:
            continue
        sub_ids = [pair_ids[i] for i in sel]
        sub_ann = [anns[i] for i in sel]
        sub_lig = lig[:, sel]
        sub_rec = rec[:, sel]
        if direction == "tumor_to_cd8":
            arms = (("hi", hi_prox, cd8_hi), ("lo", lo_prox, cd8_lo))
        else:
            arms = (("hi", cd8_hi, hi_prox), ("lo", cd8_lo, lo_prox))
        for arm, smask, rmask in arms:
            idx = np.flatnonzero(smask | rmask)
            if len(idx) < 12:
                continue
            local_s = smask[idx]
            local_r = rmask[idx]
            take = cap_indices(local_s, local_r, GEOM_CAP_SENDER, GEOM_CAP_RECEIVER, rng)
            if len(take) < 12:
                continue
            scored = spatialdm_directional(
                xy[idx][take],
                local_s[take],
                local_r[take],
                sub_lig[idx][take],
                sub_rec[idx][take],
                sub_ann,
                sub_ids,
            )
            for pid, (Ii, zi) in scored.items():
                base[pid][f"I_{arm}_prox"] = Ii
                base[pid][f"z_{arm}_prox"] = zi


def commot_fov(
    df_fov: pd.DataFrame,
    X: np.ndarray,
    gmap: dict[str, int],
    pairs_used: pd.DataFrame,
    rng: np.random.Generator,
) -> dict[str, dict[str, float]]:
    pair_ids = pairs_used["pair_id"].tolist()
    out = {
        pid: {
            "per_sender_hi": np.nan,
            "per_sender_lo": np.nan,
            "per_supply_hi": np.nan,
            "per_supply_lo": np.nan,
            "frac_pos_hi": np.nan,
            "frac_pos_lo": np.nan,
        }
        for pid in pair_ids
    }
    tumor = df_fov["is_tumor"].to_numpy()
    cd8 = df_fov["is_cd8"].to_numpy()
    tert = df_fov["cldn4_tertile"].to_numpy()
    nn_cd8 = df_fov["nn_cd8_um"].to_numpy()
    nn_hi = df_fov["nn_tumor_hi_um"].to_numpy()
    nn_lo = df_fov["nn_tumor_lo_um"].to_numpy()
    hi = tumor & (tert == 1) & np.isfinite(nn_cd8) & (nn_cd8 <= PROX_UM)
    lo = tumor & (tert == -1) & np.isfinite(nn_cd8) & (nn_cd8 <= PROX_UM)
    cd = cd8 & (
        (np.isfinite(nn_hi) & (nn_hi <= PROX_UM))
        | (np.isfinite(nn_lo) & (nn_lo <= PROX_UM))
    )
    labels = np.full(len(df_fov), "other", dtype=object)
    labels[hi] = "tumor_hi"
    labels[lo] = "tumor_lo"
    labels[cd] = "cd8"
    chosen = []
    for lab, cap in COMMOT_CAP.items():
        w = np.flatnonzero(labels == lab)
        if len(w) > cap:
            w = rng.choice(w, size=cap, replace=False)
        chosen.append(w)
    if min(len(c) for c in chosen) < MIN_SENDER and len(chosen[2]) < MIN_RECEIVER:
        return out
    if len(chosen[0]) < MIN_SENDER or len(chosen[1]) < MIN_SENDER or len(chosen[2]) < MIN_RECEIVER:
        return out
    idx = np.sort(np.concatenate(chosen))
    sub_lab = labels[idx]
    # Genes actually required
    needed = []
    for r in pairs_used.itertuples(index=False):
        needed.extend(subunits(r.ligand))
        needed.extend(subunits(r.receptor))
    needed = [g for g in dict.fromkeys(needed) if g in gmap]
    if len(needed) < 2:
        return out
    mat = np.asarray(X[idx][:, [gmap[g] for g in needed]], dtype=np.float64)
    adata = ad.AnnData(X=sparse.csr_matrix(mat))
    adata.var_names = needed
    adata.obsm["spatial"] = df_fov.iloc[idx][["x_um", "y_um"]].to_numpy(float)
    df_lr = pairs_used.loc[:, ["ligand", "receptor", "pathway"]].reset_index(drop=True)
    try:
        ct.tl.spatial_communication(
            adata,
            database_name="cosmx",
            df_ligrec=df_lr,
            pathway_sum=False,
            heteromeric=True,
            heteromeric_rule="min",
            dis_thr=COMMOT_DIS_UM,
            cot_eps_mu=1e-1,
            cot_eps_nu=1e-1,
            cot_nitermax=COMMOT_NITER,
        )
    except Exception as exc:  # noqa: BLE001
        log(f"COMMOT failed: {exc}")
        return out
    ix = {lab: np.flatnonzero(sub_lab == lab) for lab in ("tumor_hi", "tumor_lo", "cd8")}
    for r in pairs_used.itertuples(index=False):
        key = f"commot-cosmx-{r.ligand}-{r.receptor}"
        if key not in adata.obsp:
            continue
        S = adata.obsp[key].tocsr()
        L = per_cell_min(mat, {g: i for i, g in enumerate(needed)}, r.ligand)
        if L is None:
            continue
        if r.direction == "tumor_to_cd8":
            blocks = (("hi", ix["tumor_hi"]), ("lo", ix["tumor_lo"]))
            recv = ix["cd8"]
        else:
            blocks = (("hi", ix["cd8"]), ("lo", ix["cd8"]))
            # receiver differs: tumor class. Handle inside loop.
            recv = None
        for arm, senders in blocks:
            if r.direction == "tumor_to_cd8":
                receivers = recv
            else:
                receivers = ix["tumor_hi"] if arm == "hi" else ix["tumor_lo"]
            if len(senders) == 0 or len(receivers) == 0:
                continue
            block = S[senders][:, receivers]
            mass = float(block.sum())
            out[r.pair_id][f"per_sender_{arm}"] = mass / float(len(senders))
            supply = float(L[senders].sum())
            out[r.pair_id][f"per_supply_{arm}"] = mass / supply if supply > 0 else np.nan
            row_sum = np.asarray(block.sum(axis=1)).ravel()
            out[r.pair_id][f"frac_pos_{arm}"] = float((row_sum > 0).mean()) if len(row_sum) else np.nan
    return out


def squidpy_slide(
    df: pd.DataFrame,
    X: np.ndarray,
    gmap: dict[str, int],
    pairs_used: pd.DataFrame,
) -> list[dict]:
    tumor = df["is_tumor"].to_numpy()
    cd8 = df["is_cd8"].to_numpy()
    tert = df["cldn4_tertile"].to_numpy()
    nn_cd8 = df["nn_cd8_um"].to_numpy()
    nn_tumor = df["nn_tumor_um"].to_numpy()
    prox_tumor = tumor & np.isfinite(nn_cd8) & (nn_cd8 <= PROX_UM) & np.isin(tert, [-1, 1])
    prox_cd8 = cd8 & np.isfinite(nn_tumor) & (nn_tumor <= PROX_UM)
    labels = np.full(len(df), "other", dtype=object)
    labels[prox_tumor & (tert == 1)] = "tumor_hi"
    labels[prox_tumor & (tert == -1)] = "tumor_lo"
    labels[prox_cd8] = "cd8"
    keep = labels != "other"
    rows = []
    n_hi = int((labels == "tumor_hi").sum())
    n_lo = int((labels == "tumor_lo").sum())
    n_cd = int((labels == "cd8").sum())
    if n_hi < 20 or n_lo < 20 or n_cd < 15:
        log(f"{df['sample'].iloc[0]}: squidpy skipped (interface hi/lo/cd8={n_hi}/{n_lo}/{n_cd})")
        return rows
    colmap = {}
    for r in pairs_used.itertuples(index=False):
        for uname in (r.ligand, r.receptor):
            an = atomic_name(uname)
            if an not in colmap:
                vec = per_cell_min(X, gmap, uname)
                if vec is None:
                    continue
                colmap[an] = vec
    if len(colmap) < 2:
        return rows
    names = list(colmap)
    mat = np.column_stack([colmap[n] for n in names])[keep]
    obs = pd.DataFrame({"niche": pd.Categorical(labels[keep], categories=["tumor_hi", "tumor_lo", "cd8"])})
    adata = ad.AnnData(X=mat, obs=obs)
    adata.var_names = names
    inter_rows = []
    for r in pairs_used.itertuples(index=False):
        inter_rows.append({"source": atomic_name(r.ligand), "target": atomic_name(r.receptor), "pair_id": r.pair_id})
    inter = pd.DataFrame(inter_rows)
    log(f"{df['sample'].iloc[0]}: squidpy ligrec n={adata.n_obs} pairs={len(inter)}")
    res = sq.gr.ligrec(
        adata,
        cluster_key="niche",
        interactions=inter.loc[:, ["source", "target"]],
        complex_policy="min",
        n_perms=1000,
        threshold=0.0,
        seed=SEED,
        copy=True,
        use_raw=False,
        n_jobs=1,
        numba_parallel=False,
        corr_method="fdr_bh",
        corr_axis="interactions",
        clusters=[("tumor_hi", "cd8"), ("tumor_lo", "cd8"), ("cd8", "tumor_hi"), ("cd8", "tumor_lo")],
    )
    means, pvalues = res["means"], res["pvalues"]

    def grab(frame, src, tgt, c1, c2):
        try:
            v = frame.loc[(src, tgt), (c1, c2)]
        except KeyError:
            return np.nan
        arr = np.asarray(v, dtype=float).reshape(-1)
        return float(arr[0]) if len(arr) else np.nan

    for r in inter.itertuples(index=False):
        meta = pairs_used.loc[pairs_used["pair_id"] == r.pair_id].iloc[0]
        if meta.direction == "tumor_to_cd8":
            m_hi = grab(means, r.source, r.target, "tumor_hi", "cd8")
            m_lo = grab(means, r.source, r.target, "tumor_lo", "cd8")
            p_hi = grab(pvalues, r.source, r.target, "tumor_hi", "cd8")
            p_lo = grab(pvalues, r.source, r.target, "tumor_lo", "cd8")
        else:
            m_hi = grab(means, r.source, r.target, "cd8", "tumor_hi")
            m_lo = grab(means, r.source, r.target, "cd8", "tumor_lo")
            p_hi = grab(pvalues, r.source, r.target, "cd8", "tumor_hi")
            p_lo = grab(pvalues, r.source, r.target, "cd8", "tumor_lo")
        rows.append(
            {
                "sample": df["sample"].iloc[0],
                "patient": df["patient"].iloc[0],
                "pair_id": r.pair_id,
                "n_tumor_hi_interface": n_hi,
                "n_tumor_lo_interface": n_lo,
                "n_cd8_interface": n_cd,
                "mean_hi": m_hi,
                "mean_lo": m_lo,
                "delta_mean": m_hi - m_lo if np.isfinite(m_hi) and np.isfinite(m_lo) else np.nan,
                "p_hi": p_hi,
                "p_lo": p_lo,
            }
        )
    return rows


def fov_records(
    df: pd.DataFrame,
    X: np.ndarray,
    gmap: dict[str, int],
    pairs_used: pd.DataFrame,
    max_fovs: int | None,
    skip_commot: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    lig, rec, ids = expr_block(X, gmap, pairs_used)
    id_to_col = {p: i for i, p in enumerate(ids)}
    used = pairs_used[pairs_used["pair_id"].isin(ids)].reset_index(drop=True)
    # Reorder lig/rec to used order
    order = [id_to_col[p] for p in used["pair_id"]]
    lig, rec = lig[:, order], rec[:, order]
    sample = df["sample"].iloc[0]
    fovs = list(df.groupby("fov").groups)
    if max_fovs is not None:
        # Prefer FOVs with enough of both arms so the pilot is informative.
        ranked = []
        for fov, idx in df.groupby("fov").groups.items():
            g = df.loc[idx]
            n_hi = int(((g["is_tumor"]) & (g["cldn4_tertile"] == 1)).sum())
            n_lo = int(((g["is_tumor"]) & (g["cldn4_tertile"] == -1)).sum())
            n_cd = int(g["is_cd8"].sum())
            ranked.append((min(n_hi, n_lo, n_cd), fov))
        ranked.sort(reverse=True)
        fovs = [f for _, f in ranked[:max_fovs]]
    rows = []
    geos = []
    for i, fov in enumerate(fovs, start=1):
        g = df.loc[df["fov"] == fov]
        n_hi = int(((g["is_tumor"]) & (g["cldn4_tertile"] == 1)).sum())
        n_lo = int(((g["is_tumor"]) & (g["cldn4_tertile"] == -1)).sum())
        n_cd = int(g["is_cd8"].sum())
        if n_hi < 10 or n_lo < 10 or n_cd < MIN_RECEIVER:
            continue
        pos = g.index.to_numpy()
        rng = rng_for(sample, int(fov))
        t0 = time.time()
        spatial = score_modes_for_fov(g, lig[pos], rec[pos], used, rng)
        comm = (
            {pid: {"per_sender_hi": np.nan, "per_sender_lo": np.nan, "per_supply_hi": np.nan, "per_supply_lo": np.nan} for pid in used["pair_id"]}
            if skip_commot
            else commot_fov(g, X[pos], gmap, used, rng)
        )
        nn_hi = g.loc[g["is_tumor"] & (g["cldn4_tertile"] == 1), "nn_cd8_um"]
        nn_lo = g.loc[g["is_tumor"] & (g["cldn4_tertile"] == -1), "nn_cd8_um"]
        geos.append(
            {
                "sample": sample,
                "patient": PATIENT[sample],
                "fov": int(fov),
                "n_tumor_hi": n_hi,
                "n_tumor_lo": n_lo,
                "n_cd8": n_cd,
                "mean_nn_cd8_um_hi": float(nn_hi.mean()) if len(nn_hi) else np.nan,
                "mean_nn_cd8_um_lo": float(nn_lo.mean()) if len(nn_lo) else np.nan,
            }
        )
        for r in used.itertuples(index=False):
            sp = spatial[r.pair_id]
            cm = comm[r.pair_id]
            rows.append(
                {
                    "sample": sample,
                    "patient": PATIENT[sample],
                    "fov": int(fov),
                    "pair_id": r.pair_id,
                    "ligand": r.ligand,
                    "receptor": r.receptor,
                    "family": r.family,
                    "direction": r.direction,
                    "annotation": r.annotation,
                    "n_tumor_hi": n_hi,
                    "n_tumor_lo": n_lo,
                    "n_cd8": n_cd,
                    **{k: sp[k] for k in sp},
                    **{f"commot_{k}": cm[k] for k in cm},
                }
            )
        log(f"{sample} fov {fov} ({i}/{len(fovs)}) {time.time()-t0:.1f}s cells={len(g)}")
    return pd.DataFrame(rows), pd.DataFrame(geos)


def delta(a, b):
    if np.isfinite(a) and np.isfinite(b):
        return float(a) - float(b)
    return np.nan


def add_deltas(fov: pd.DataFrame) -> pd.DataFrame:
    if fov.empty:
        return fov
    fov = fov.copy()
    fov["d_sdm_prox"] = [delta(a, b) for a, b in zip(fov["I_hi_prox"], fov["I_lo_prox"])]
    fov["d_sdm_geom"] = [delta(a, b) for a, b in zip(fov["I_hi_geom"], fov["I_lo_geom"])]
    fov["d_commot"] = [delta(a, b) for a, b in zip(fov["commot_per_sender_hi"], fov["commot_per_sender_lo"])]
    fov["d_commot_supply"] = [delta(a, b) for a, b in zip(fov["commot_per_supply_hi"], fov["commot_per_supply_lo"])]
    fov["d_expr"] = [delta(a, b) for a, b in zip(fov["expr_hi"], fov["expr_lo"])]
    return fov


def wilcoxon_vec(x: np.ndarray) -> dict:
    v = np.asarray(x, dtype=float)
    v = v[np.isfinite(v)]
    n_neg = int((v < 0).sum())
    n_pos = int((v > 0).sum())
    n_zero = int((v == 0).sum())
    out = {"n": int(len(v)), "n_neg": n_neg, "n_pos": n_pos, "n_zero": n_zero, "stat": np.nan, "p": np.nan}
    if len(v) < 5 or np.allclose(v, 0):
        return out
    try:
        res = stats.wilcoxon(v, zero_method="wilcox", alternative="two-sided", method="auto")
        out["stat"] = float(res.statistic)
        out["p"] = float(res.pvalue)
    except ValueError:
        return out
    return out


def summarize(fov: pd.DataFrame, squid: pd.DataFrame, geo: pd.DataFrame, pairs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if fov.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    value_cols = ["d_sdm_prox", "d_sdm_geom", "d_commot", "d_commot_supply", "d_expr"]
    slide = fov.groupby(["sample", "patient", "pair_id"], as_index=False)[value_cols].mean()
    nf = fov.groupby(["sample", "pair_id"])["fov"].nunique().rename("n_fov").reset_index()
    slide = slide.merge(nf, on=["sample", "pair_id"], how="left")
    if not squid.empty:
        slide = slide.merge(
            squid.loc[:, ["sample", "pair_id", "delta_mean", "mean_hi", "mean_lo", "p_hi", "p_lo", "n_tumor_hi_interface", "n_cd8_interface"]],
            on=["sample", "pair_id"],
            how="left",
        )
        slide = slide.rename(columns={"delta_mean": "d_squidpy"})
    else:
        slide["d_squidpy"] = np.nan
    endpoints = ["d_sdm_prox", "d_sdm_geom", "d_commot", "d_commot_supply", "d_expr", "d_squidpy"]
    rows = []
    for pid, g in slide.groupby("pair_id"):
        meta = pairs.loc[pairs["pair_id"] == pid].iloc[0]
        rec = {
            "pair_id": pid,
            "ligand": meta.ligand,
            "receptor": meta.receptor,
            "pathway": meta.pathway,
            "annotation": meta.annotation,
            "family": meta.family,
            "direction": meta.direction,
            "in_cellchat_human": bool(meta.in_cellchat_human),
        }
        for ep in endpoints:
            w = wilcoxon_vec(g[ep].to_numpy())
            rec[f"{ep}_median"] = float(np.nanmedian(g[ep].to_numpy())) if np.isfinite(g[ep]).any() else np.nan
            rec[f"{ep}_mean"] = float(np.nanmean(g[ep].to_numpy())) if np.isfinite(g[ep]).any() else np.nan
            rec[f"{ep}_n_slides"] = w["n"]
            rec[f"{ep}_n_neg"] = w["n_neg"]
            rec[f"{ep}_n_pos"] = w["n_pos"]
            rec[f"{ep}_p"] = w["p"]
            # Patient collapse: mean of slide deltas, then sign count. n=5, report signs.
            pg = g.groupby("patient")[ep].mean()
            pv = pg.to_numpy(dtype=float)
            pv = pv[np.isfinite(pv)]
            rec[f"{ep}_n_patients"] = int(len(pv))
            rec[f"{ep}_n_patients_neg"] = int((pv < 0).sum())
            rec[f"{ep}_n_patients_pos"] = int((pv > 0).sum())
        rows.append(rec)
    summary = pd.DataFrame(rows)
    # BH within family for the two spatial endpoints and squidpy
    for ep in ("d_sdm_prox", "d_sdm_geom", "d_commot", "d_squidpy"):
        summary[f"{ep}_q"] = np.nan
        for fam, idx in summary.groupby("family").groups.items():
            p = summary.loc[idx, f"{ep}_p"].to_numpy(dtype=float)
            ok = np.isfinite(p)
            if ok.sum() == 0:
                continue
            q = np.full(len(p), np.nan)
            _, qv, _, _ = multipletests(p[ok], method="fdr_bh")
            q[ok] = qv
            summary.loc[idx, f"{ep}_q"] = q
    # Geometry sanity: FOV distance delta vs SpatialDM geometry delta, MHC control and barrier mean
    sanity_rows = []
    if not geo.empty:
        geo = geo.copy()
        geo["d_nn"] = geo["mean_nn_cd8_um_hi"] - geo["mean_nn_cd8_um_lo"]
        for fam in ("mhc_control", "barrier"):
            sub = fov[fov["family"] == fam]
            if sub.empty:
                continue
            md = sub.groupby(["sample", "fov"])["d_sdm_geom"].mean().rename("d_sdm_geom").reset_index()
            m = geo.merge(md, on=["sample", "fov"], how="inner")
            ok = m["d_nn"].notna() & m["d_sdm_geom"].notna()
            if ok.sum() >= 10:
                rho, p = stats.spearmanr(m.loc[ok, "d_nn"], m.loc[ok, "d_sdm_geom"])
            else:
                rho, p = np.nan, np.nan
            sanity_rows.append({"family": fam, "spearman_dnn_vs_dI_geom": rho, "p": p, "n_fov": int(ok.sum())})
    sanity = pd.DataFrame(sanity_rows)
    return slide, summary, sanity


def _fmt(x, nd=3):
    if x is None or not np.isfinite(x):
        return "NA"
    ax = abs(float(x))
    if ax != 0 and ax < 1e-3:
        return f"{float(x):.2e}"
    return f"{float(x):.{nd}f}"


def _fmt_p(x):
    if x is None or not np.isfinite(x):
        return "NA"
    if x < 1e-3:
        return f"{float(x):.2e}"
    return f"{float(x):.3f}"


def _pair_row(summary: pd.DataFrame, pair_id: str) -> pd.Series:
    hit = summary.loc[summary["pair_id"] == pair_id]
    if hit.empty:
        raise KeyError(pair_id)
    return hit.iloc[0]


def _slide_phrase(row: pd.Series, ep: str) -> str:
    n = int(row[f"{ep}_n_slides"])
    n_pos = int(row[f"{ep}_n_pos"])
    p = _fmt_p(row[f"{ep}_p"])
    q = _fmt_p(row[f"{ep}_q"]) if f"{ep}_q" in row.index and np.isfinite(row[f"{ep}_q"]) else None
    pt_pos = int(row[f"{ep}_n_patients_pos"])
    pt_n = int(row[f"{ep}_n_patients"])
    qbit = f", BH q={q}" if q else ""
    return (
        f"{n_pos}/{n} slides Δ>0 (median {_fmt_signed(row[f'{ep}_median'])}, "
        f"Wilcoxon p={p}{qbit}; {pt_pos}/{pt_n} tissues)"
    )


def _fmt_signed(x, nd=3):
    if x is None or not np.isfinite(x):
        return "NA"
    ax = abs(float(x))
    body = f"{ax:.2e}" if ax != 0 and ax < 1e-3 else f"{ax:.{nd}f}"
    if float(x) > 0:
        return "+" + body
    if float(x) < 0:
        return "−" + body
    return "0"


def write_results_md(path: Path, summary: pd.DataFrame, sanity: pd.DataFrame, inventory: dict) -> None:
    def line_for(ep: str, family: str) -> str:
        sub = summary[summary["family"] == family]
        if sub.empty:
            return "no pairs"
        bits = []
        for r in sub.itertuples(index=False):
            qv = getattr(r, f"{ep}_q", np.nan)
            qbit = f", BH q={_fmt_p(qv)}" if np.isfinite(qv) else ""
            bits.append(
                f"{r.pair_id} {int(getattr(r, f'{ep}_n_pos'))}/{int(getattr(r, f'{ep}_n_slides'))} slides Δ>0, "
                f"median Δ={_fmt_signed(getattr(r, f'{ep}_median'))}, Wilcoxon p={_fmt_p(getattr(r, f'{ep}_p'))}{qbit}, "
                f"tissues {int(getattr(r, f'{ep}_n_patients_pos'))}/{int(getattr(r, f'{ep}_n_patients'))} Δ>0"
            )
        return "; ".join(bits)

    def count_consistent(ep: str, family: str, sign: str) -> tuple[int, int]:
        sub = summary[summary["family"] == family]
        n = len(sub)
        if sign == "neg":
            k = int((sub[f"{ep}_n_neg"] >= 6).sum()) if n else 0
        else:
            k = int((sub[f"{ep}_n_pos"] >= 6).sum()) if n else 0
        return k, n

    k_neg, n_b = count_consistent("d_sdm_prox", "barrier", "neg")
    k_pos, _ = count_consistent("d_sdm_prox", "barrier", "pos")
    c_neg, _ = count_consistent("d_commot", "barrier", "neg")
    c_pos, _ = count_consistent("d_commot", "barrier", "pos")
    lines = []
    lines.append("# CosMx NSCLC — spatial ligand–receptor CCC, CLDN4-high vs CLDN4-low niches vs CD8")
    lines.append("")
    lines.append("CLDN4-only. Official CosMx SMI 960-plex NSCLC (He et al. 2022; NanoString public S3; the same 8 slides / 5 tissues as figshare 25976224). No private 8-KL. No TACSTD2 gate. No ICI labels.")
    lines.append("")
    lines.append("This file does **not** replace the locked contact odds ratio or the Ripley g(r) result. It asks a different question: among CellChat ligand–receptor pairs that are actually on the 960-plex, is spatial coupling to CD8 different for CLDN4-high tumor than for CLDN4-low tumor?")
    lines.append("")
    lines.append("## Definitions")
    lines.append("")
    lines.append("| Item | Choice |")
    lines.append("|---|---|")
    lines.append("| QC | Drop `cell_ID==0`; ≥20 counts and ≥5 genes on the 960-plex |")
    lines.append("| Tumor | RNA epithelial marker-score argmax |")
    lines.append("| CD8 | `CD8A\\|CD8B > 0` and `CD3D\\|CD3E\\|CD3G > 0` and not epithelial |")
    lines.append("| CLDN4 split | Per-FOV rank tertiles among tumor cells (Q3 vs Q1) |")
    lines.append("| Pixel size | 0.18 µm |")
    lines.append("| Proximity gate | Centroid distance ≤ 40 µm, within FOV |")
    lines.append("| Usable FOV | ≥10 CLDN4-high tumor, ≥10 CLDN4-low tumor, ≥5 CD8 |")
    lines.append("| Expression | log1p, median-library normalized |")
    lines.append("| Tools | squidpy 1.6.5 `gr.ligrec`; SpatialDM 0.3.1 global z-score; COMMOT 0.0.3 collective OT |")
    lines.append("")
    absent = ", ".join(inventory.get("absent_genes") or [])
    lines.append(f"**{absent} are not on this panel.** No score is reported for them. TIGIT is on the panel, but PVR, NECTIN2, and CD226 are not, so the TIGIT axis is not scored.")
    lines.append("")
    lines.append("HLA-A/B/C → CD8A/CD8B is a **proximity control**, not an independent checkpoint result: CD8A/CD8B are part of the CD8 gate, so the receptor is nearly constant on CD8 cells and the spatial score mostly tracks whether tumor cells sit near CD8 cells.")
    lines.append("")
    lines.append("### What each tool is allowed to mean")
    lines.append("")
    lines.append("- **SpatialDM** (primary spatial test). Single-cell RBF, `l=30` µm, cutoff 0.15, up to 80 neighbors for secreted pairs and 10 nearest neighbors for contact pairs. Ligand expression is zeroed outside the sender class and receptor expression outside the receiver class. Weights are kept only on sender→receiver edges. Two cell sets: **geometry** (all cells of the two classes in the FOV) and **proximity** (only cells within 40 µm of the partner class). The slide endpoint is the mean over FOVs of I(high) − I(low).")
    lines.append("- **COMMOT**. Collective optimal transport on proximity-conditioned cells (cap 140/140/120), `dis_thr=50` µm, heteromeric subunits combined by per-cell minimum, `cot_nitermax=1000`. Endpoint: transport mass into the receiver class divided by the number of senders. A supply-normalized version (mass / ligand sum) is secondary.")
    lines.append("- **Squidpy `gr.ligrec`**. CellPhoneDB permutation (1,000 label shuffles, seed 20260921) on the 40 µm interface, clusters `tumor_hi`, `tumor_lo`, `cd8`. Heteromeric genes are pre-collapsed with a per-cell minimum so the tool does not replace a complex by a single globally low subunit. The score is the average of the sender-cluster ligand mean and the receiver-cluster receptor mean. With one shared CD8 receiver, high-vs-low **cancels the receptor** and equals half the ligand-mean difference, so pairs that share a ligand share one squidpy delta. It is a niche-restricted abundance test, not a distance-weighted transport score. The tool p-value is a cell-level permutation and is not the primary claim.")
    lines.append("")
    lines.append("Primary inference is a Wilcoxon signed-rank test on the **8 slide-level deltas** (high − low). Negative means CLDN4-high scores lower than CLDN4-low. Patient signs collapse Lung5 and Lung9 technical replicates. BH q-values are within family.")
    lines.append("")
    lines.append("## Inventory")
    lines.append("")
    lines.append(f"- QC cells: {inventory.get('n_qc', 'NA')}")
    lines.append(f"- Usable FOVs / slides: {inventory.get('n_fov', 'NA')} / {inventory.get('n_slides', 'NA')}")
    lines.append(f"- Pairs on panel and in the CellChat human table shipped with COMMOT: {inventory.get('n_pairs_in_db', 'NA')} / {inventory.get('n_pairs', 'NA')}")
    lines.append("")
    lines.append("## Result")
    lines.append("")
    lines.append(f"Barrier pairs with SpatialDM proximity ΔI lower on at least 6 slides: **{k_neg}/{n_b}**. Higher on at least 6 slides: **{k_pos}/{n_b}**.")
    lines.append(f"Barrier pairs with COMMOT per-sender transport lower on at least 6 slides: **{c_neg}/{n_b}**. Higher on at least 6 slides: **{c_pos}/{n_b}**.")
    lines.append("")
    lines.append("### Reading")
    lines.append("")
    lines.append("CLDN4-high tumor is not a ligand–receptor desert next to CD8, and the galectin-9–TIM-3 and PD-L1–PD-1 pairs are not the pairs that move. Where a COMMOT per-sender increase survives BH within the barrier family, dividing that mass by ligand supply removes it. The increase tracks extra ligand on the CLDN4-high tumor cells that already have a CD8 neighbor, which is the same kind of result as the earlier contact-level CDH1 abundance contrast, not a higher delivery fraction and not evidence that CD8 cells next to CLDN4-high tumor are transcriptionally shut down.")
    lines.append("")
    if len(summary):
        cdh = _pair_row(summary, "CDH1|CDH1")
        cdh_int = _pair_row(summary, "CDH1|ITGA2_ITGB1")
        icam = _pair_row(summary, "ICAM1|ITGAL_ITGB2")
        icam1 = _pair_row(summary, "ICAM1|ITGAL")
        mif = _pair_row(summary, "MIF|CD74_CD44")
        lg = _pair_row(summary, "LGALS9|HAVCR2")
        pdl1 = _pair_row(summary, "CD274|PDCD1")
        pdl2 = _pair_row(summary, "PDCD1LG2|PDCD1")
        tgfb = _pair_row(summary, "TGFB1|TGFBR1_TGFBR2")
        ifng = _pair_row(summary, "IFNG|IFNGR1_IFNGR2")
        cx9 = _pair_row(summary, "CXCL9|CXCR3")
        cx16 = _pair_row(summary, "CXCL16|CXCR6")
        hla = _pair_row(summary, "HLA-A|CD8A")
        lines.append(f"**CDH1 abundance, not homophilic transport.** On tumor cells within 40 µm of CD8, CLDN4-high minus CLDN4-low CDH1 is {_slide_phrase(cdh, 'd_expr')}. Squidpy, which reduces to that ligand contrast, agrees ({_slide_phrase(cdh, 'd_squidpy')}). COMMOT CDH1–CDH1 transport per sender does not ({_slide_phrase(cdh, 'd_commot')}). COMMOT CDH1–ITGA2/ITGB1 per sender does ({_slide_phrase(cdh_int, 'd_commot')}), and the supply-normalized version does not ({_slide_phrase(cdh_int, 'd_commot_supply')}).")
        lines.append("")
        lines.append(f"**ICAM1–LFA-1 and MIF–CD74/CD44 are the COMMOT hits inside the barrier family.** ICAM1–ITGAL/ITGB2 per sender {_slide_phrase(icam, 'd_commot')}; ICAM1–ITGAL {_slide_phrase(icam1, 'd_commot')}; MIF–CD74/CD44 {_slide_phrase(mif, 'd_commot')}. BH q-values are within the 12 barrier pairs. SpatialDM proximity is in the same direction for ICAM1 (ITGAL/ITGB2 {_slide_phrase(icam, 'd_sdm_prox')}; ITGAL {_slide_phrase(icam1, 'd_sdm_prox')}) but those q-values stay above 0.05. MIF SpatialDM proximity does not ({_slide_phrase(mif, 'd_sdm_prox')}). MIF ligand abundance is higher ({_slide_phrase(mif, 'd_expr')}). ICAM1 abundance is only a trend ({_slide_phrase(icam, 'd_expr')}). Supply-normalized COMMOT is not significant for these pairs (ICAM1–ITGAL/ITGB2 {_slide_phrase(icam, 'd_commot_supply')}; MIF–CD74/CD44 {_slide_phrase(mif, 'd_commot_supply')}).")
        lines.append("")
        lines.append(f"**Checkpoint pairs do not mark CLDN4-high contacts.** LGALS9–HAVCR2 SpatialDM proximity {_slide_phrase(lg, 'd_sdm_prox')}; COMMOT {_slide_phrase(lg, 'd_commot')}; LGALS9 abundance {_slide_phrase(lg, 'd_expr')}. CD274–PDCD1 SpatialDM {_slide_phrase(pdl1, 'd_sdm_prox')}; CD274 abundance {_slide_phrase(pdl1, 'd_expr')}. PDCD1LG2–PDCD1 is the barrier pair that leans lower (SpatialDM {_slide_phrase(pdl2, 'd_sdm_prox')}; abundance {_slide_phrase(pdl2, 'd_expr')}) and is not significant. TGFB1 abundance {_slide_phrase(tgfb, 'd_expr')}; COMMOT {_slide_phrase(tgfb, 'd_commot')}.")
        lines.append("")
        lines.append(f"**CD8 IFNG next to CLDN4-high tumor is flat.** IFNG on CD8 within 40 µm of CLDN4-high tumor minus CD8 within 40 µm of CLDN4-low tumor: {_slide_phrase(ifng, 'd_expr')}. SpatialDM {_slide_phrase(ifng, 'd_sdm_prox')}. COMMOT {_slide_phrase(ifng, 'd_commot')}.")
        lines.append("")
        lines.append(f"**HLA–CD8 is a proximity control, and the two spatial tools disagree.** SpatialDM proximity HLA-A–CD8A is higher around CLDN4-high tumor ({_slide_phrase(hla, 'd_sdm_prox')}). COMMOT transport per sender is lower ({_slide_phrase(hla, 'd_commot')}), as is the supply-normalized mass ({_slide_phrase(hla, 'd_commot_supply')}). CD8A/CD8B are part of the CD8 definition, so this is not a checkpoint result. Lower COMMOT mass is what a distance-limited transporter does when fewer CD8 cells sit in range. Higher Moran I is a correlation among the cells that were kept, not a delivered-mass estimate.")
        lines.append("")
        lines.append(f"**Chemokines are not one sign.** Tumor CXCL9 is lower ({_slide_phrase(cx9, 'd_expr')}) without a significant CXCL9–CXCR3 COMMOT or SpatialDM contrast ({_slide_phrase(cx9, 'd_commot')}; {_slide_phrase(cx9, 'd_sdm_prox')}). Tumor CXCL16 is higher ({_slide_phrase(cx16, 'd_expr')}); CXCL16–CXCR6 COMMOT per sender {_slide_phrase(cx16, 'd_commot')} and does not survive BH inside the chemokine family (q above). Supply normalization removes it ({_slide_phrase(cx16, 'd_commot_supply')}).")
        lines.append("")
    lines.append("A high-side increase that remains after supply normalization would be the spatial-CCC version of more delivery, not just more ligand. That is not what the barrier pairs do. A high-side decrease on the geometry run, especially for HLA–CD8, would have restated exclusion as missing coupling. SpatialDM does not show that decrease. The proximity run asks about coupling given a neighbor within 40 µm.")
    lines.append("")
    lines.append("![Slides with a lower CLDN4-high score](results/cosmx_cldn4_spatial_ccc/figures/barrier_slides_delta_negative.png)")
    lines.append("")
    lines.append("![SpatialDM proximity slide deltas](results/cosmx_cldn4_spatial_ccc/figures/spatialdm_prox_delta_heatmap.png)")
    lines.append("")
    lines.append("![COMMOT slide deltas](results/cosmx_cldn4_spatial_ccc/figures/commot_delta_heatmap.png)")
    lines.append("")
    lines.append("")
    lines.append("### Barrier family — SpatialDM proximity (I high − I low)")
    lines.append("")
    lines.append(line_for("d_sdm_prox", "barrier"))
    lines.append("")
    lines.append("### Barrier family — SpatialDM geometry")
    lines.append("")
    lines.append(line_for("d_sdm_geom", "barrier"))
    lines.append("")
    lines.append("### Barrier family — COMMOT transport per sender")
    lines.append("")
    lines.append(line_for("d_commot", "barrier"))
    lines.append("")
    lines.append("### Barrier family — squidpy interface mean (ligand-driven if the CD8 receiver is shared)")
    lines.append("")
    lines.append(line_for("d_squidpy", "barrier"))
    lines.append("")
    lines.append("### Barrier family — sender expression on the 40 µm interface (not a CCC score)")
    lines.append("")
    lines.append(line_for("d_expr", "barrier"))
    lines.append("")
    lines.append("### MHC–CD8 proximity control")
    lines.append("")
    lines.append("Geometry: " + line_for("d_sdm_geom", "mhc_control"))
    lines.append("")
    lines.append("Proximity: " + line_for("d_sdm_prox", "mhc_control"))
    lines.append("")
    lines.append("### Chemokine, other, and CD8→tumor IFNG")
    lines.append("")
    lines.append("Chemokine, SpatialDM proximity: " + line_for("d_sdm_prox", "chemokine"))
    lines.append("")
    lines.append("Chemokine, COMMOT: " + line_for("d_commot", "chemokine"))
    lines.append("")
    lines.append("Other, SpatialDM proximity: " + line_for("d_sdm_prox", "other"))
    lines.append("")
    lines.append("CD8→tumor IFNG, SpatialDM proximity: " + line_for("d_sdm_prox", "effector_reverse"))
    lines.append("")
    lines.append("CD8→tumor IFNG, COMMOT: " + line_for("d_commot", "effector_reverse"))
    lines.append("")
    lines.append("CD8→tumor IFNG, sender expression (IFNG on CD8 near CLDN4-high tumor minus CD8 near CLDN4-low tumor): " + line_for("d_expr", "effector_reverse"))
    lines.append("")
    if sanity is not None and len(sanity):
        lines.append("### Distance sanity check")
        lines.append("")
        lines.append("Spearman correlation, across FOVs, of (mean nearest-CD8 distance, CLDN4-high minus CLDN4-low) versus the mean SpatialDM geometry ΔI. A negative correlation means FOVs where CLDN4-high tumor sits farther from CD8 also show less geometry-inclusive spatial coupling.")
        lines.append("")
        for r in sanity.itertuples(index=False):
            lines.append(f"- {r.family}: ρ={_fmt(r.spearman_dnn_vs_dI_geom)} , p={_fmt_p(r.p)} , n FOV={int(r.n_fov)}")
        lines.append("")
    lines.append("## What is not claimed")
    lines.append("")
    lines.append("- No CellChat `computeCommunProb` and no causal barrier.")
    lines.append("- No F11R result and no NECTIN2–TIGIT result. F11R, NECTIN2, PVR, and CD226 are absent. TIGIT is on the panel without those partners.")
    lines.append("- Squidpy high-vs-low is not a distance-weighted communication probability.")
    lines.append("- HLA–CD8 scores are not evidence of a CD8-receptor checkpoint. CD8A/CD8B define the receiver class.")
    lines.append("- Slide Wilcoxon n=8. Patient signs are 5 tissues, with Lung5 and Lung9 replicates averaged.")
    lines.append("- COMMOT was capped (140 high / 140 low / 120 CD8 proximal cells per FOV) and stopped at 1,000 OT iterations.")
    lines.append("- No private 8-KL and no ICI response labels.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 -m venv .venv && . .venv/bin/activate")
    lines.append("pip install -r scripts/requirements-cosmx-spatial-ccc.txt")
    lines.append("bash scripts/download_cosmx_nsclc.sh")
    lines.append("python3 scripts/cosmx_cldn4_spatial_ccc.py")
    lines.append("```")
    lines.append("")
    lines.append("Tables: `results/cosmx_cldn4_spatial_ccc/`. Machine summary: `stats.json`.")
    lines.append("")
    path.write_text("\n".join(lines))


def plot_heatmap(slide: pd.DataFrame, summary: pd.DataFrame, value: str, title: str, path: Path) -> None:
    if slide.empty or value not in slide.columns:
        return
    order = summary.sort_values(["family", "pair_id"])["pair_id"].tolist()
    samples = [s for s in SAMPLES if s in set(slide["sample"])]
    mat = np.full((len(order), len(samples)), np.nan)
    for i, pid in enumerate(order):
        for j, s in enumerate(samples):
            hit = slide[(slide["pair_id"] == pid) & (slide["sample"] == s)]
            if len(hit):
                mat[i, j] = hit.iloc[0][value]
    # Scale each pair by its own max abs so the color is the sign pattern.
    scale = np.nanmax(np.abs(mat), axis=1)
    scale[~np.isfinite(scale) | (scale == 0)] = 1.0
    show = mat / scale[:, None]
    fig_h = max(6.0, 0.28 * len(order) + 1.5)
    fig, ax = plt.subplots(figsize=(8.5, fig_h))
    im = ax.imshow(show, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("slide Δ / max|Δ| within pair\n(red: CLDN4-high higher)")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_sign_bars(summary: pd.DataFrame, path: Path) -> None:
    sub = summary[summary["family"].isin(["barrier", "effector_reverse"])].copy()
    if sub.empty:
        return
    endpoints = [
        ("d_sdm_prox", "SpatialDM prox"),
        ("d_commot", "COMMOT"),
        ("d_squidpy", "Squidpy"),
    ]
    labels = sub["pair_id"].tolist()
    x = np.arange(len(labels))
    width = 0.25
    fig, ax = plt.subplots(figsize=(12, 4.8))
    for i, (ep, name) in enumerate(endpoints):
        vals = sub[f"{ep}_n_neg"].to_numpy(dtype=float)
        ax.bar(x + (i - 1) * width, vals, width=width, label=name)
    ax.axhline(4, color="0.4", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Slides with Δ < 0 (of 8)")
    ax.set_ylim(0, 8.5)
    ax.set_title("CLDN4-high minus CLDN4-low: slides with a lower score")
    ax.legend(frameon=False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def self_test() -> None:
    rng = np.random.default_rng(0)
    # Close tumor ligand next to CD8 receptor should score a higher Moran I
    # than the same cells placed far apart.
    n_t, n_c = 30, 20
    xy_close = np.zeros((n_t + n_c, 2))
    xy_close[:n_t, 0] = np.linspace(0, 40, n_t)
    xy_close[n_t:, 0] = np.linspace(1, 41, n_c)
    xy_far = xy_close.copy()
    xy_far[n_t:, 0] += 400
    lig = np.zeros((n_t + n_c, 1))
    rec = np.zeros((n_t + n_c, 1))
    lig[:n_t, 0] = 1.0
    rec[n_t:, 0] = 1.0
    sender = np.arange(n_t + n_c) < n_t
    receiver = ~sender
    close = spatialdm_directional(xy_close, sender, receiver, lig, rec, ["Secreted Signaling"], ["P"])
    far = spatialdm_directional(xy_far, sender, receiver, lig, rec, ["Secreted Signaling"], ["P"])
    log(f"self-test SpatialDM close I={close['P'][0]:.3f} far I={far['P'][0]:.3f}")
    if not (close["P"][0] > far["P"][0]):
        raise SystemExit("SpatialDM self-test failed: close coupling was not higher than far coupling")
    # COMMOT toy from the package docstring, plus a distance contrast.
    X = np.array([[1.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    adata = ad.AnnData(X=sparse.csr_matrix(X), var=pd.DataFrame(index=["L", "R", "Z"]))
    adata.obsm["spatial"] = np.array([[0.0, 0], [5, 0], [10, 0], [80, 0]])
    df = pd.DataFrame([["L", "R", "P"]])
    ct.tl.spatial_communication(
        adata,
        database_name="toy",
        df_ligrec=df,
        heteromeric=True,
        dis_thr=30,
        cot_eps_mu=0.1,
        cot_eps_nu=0.1,
        cot_nitermax=400,
    )
    mass_near = float(adata.obsp["commot-toy-L-R"][0, 2])
    mass_far = float(adata.obsp["commot-toy-L-R"][0, 3])
    log(f"self-test COMMOT near={mass_near:.4f} far={mass_far:.4f}")
    if not (mass_near > mass_far):
        raise SystemExit("COMMOT self-test failed: transport ignored the distance threshold")
    log("self-test passed")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=DATA)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--samples", nargs="*", default=SAMPLES)
    ap.add_argument("--max-fovs", type=int, default=None)
    ap.add_argument("--skip-commot", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return

    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    (out / "checkpoints").mkdir(exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)
    pairs = confirm_pairs_in_cellchat(pair_frame())
    pairs.to_csv(out / "pairs_registered.tsv", sep="\t", index=False)
    missing_db = pairs.loc[~pairs["in_cellchat_human"], "pair_id"].tolist()
    if missing_db:
        log(f"pairs missing from bundled CellChat table: {missing_db}")

    # Panel from the first available sample.
    first = next(s for s in args.samples if list(args.data.rglob(f"{s}_exprMat_file.csv")))
    panel = read_panel(find_file(first, "exprMat", args.data))
    panel_set = set(panel)
    want = gene_universe(pairs)
    status_rows = []
    for g in sorted(set(want) | set(ABSENT_GENES)):
        status_rows.append({"gene": g, "on_960_panel": g in panel_set, "role": "queried"})
    pd.DataFrame(status_rows).to_csv(out / "panel_gene_status.tsv", sep="\t", index=False)
    absent = [g for g in ABSENT_GENES if g not in panel_set]
    log(f"confirmed absent: {absent}")

    fov_parts = []
    geo_parts = []
    squid_parts = []
    n_qc = 0
    for sample in args.samples:
        stamp_path = out / "checkpoints" / f"{sample}.{CODE_STAMP}.fov.tsv.gz"
        geo_path = out / "checkpoints" / f"{sample}.{CODE_STAMP}.geo.tsv"
        sq_path = out / "checkpoints" / f"{sample}.{CODE_STAMP}.squidpy.tsv"
        nqc_path = out / "checkpoints" / f"{sample}.nqc"
        if stamp_path.exists() and geo_path.exists() and sq_path.exists() and not args.force and args.max_fovs is None:
            log(f"{sample}: resume checkpoint")
            fov_parts.append(pd.read_csv(stamp_path, sep="\t"))
            geo_parts.append(pd.read_csv(geo_path, sep="\t"))
            squid_parts.append(pd.read_csv(sq_path, sep="\t"))
            if nqc_path.exists():
                n_qc += int(nqc_path.read_text().strip())
            continue
        df, X, genes = load_sample(sample, args.data, want, panel)
        n_qc += len(df)
        nqc_path.write_text(str(len(df)))
        gmap = {g: i for i, g in enumerate(genes)}
        df = add_classes_and_distances(df)
        # Positional alignment: load_sample returns reset index; distances keep it.
        fov, geo = fov_records(df, X, gmap, pairs, args.max_fovs, args.skip_commot)
        squid_rows = squidpy_slide(df, X, gmap, pairs)
        squid = pd.DataFrame(squid_rows)
        if len(fov):
            fov.to_csv(stamp_path, sep="\t", index=False)
        if len(geo):
            geo.to_csv(geo_path, sep="\t", index=False)
        squid.to_csv(sq_path, sep="\t", index=False)
        fov_parts.append(fov)
        geo_parts.append(geo)
        squid_parts.append(squid)
        del df, X

    fov = add_deltas(pd.concat([p for p in fov_parts if p is not None and len(p)], ignore_index=True)) if any(len(p) for p in fov_parts) else pd.DataFrame()
    geo = pd.concat([p for p in geo_parts if p is not None and len(p)], ignore_index=True) if any(len(p) for p in geo_parts) else pd.DataFrame()
    squid = pd.concat([p for p in squid_parts if p is not None and len(p)], ignore_index=True) if any(len(p) for p in squid_parts) else pd.DataFrame()
    if len(fov):
        fov.to_csv(out / "fov_scores.tsv.gz", sep="\t", index=False)
    if len(geo):
        geo.to_csv(out / "fov_geometry.tsv", sep="\t", index=False)
    if len(squid):
        squid.to_csv(out / "squidpy_slide.tsv", sep="\t", index=False)
    slide, summary, sanity = summarize(fov, squid, geo, pairs)
    if len(slide):
        slide.to_csv(out / "slide_deltas.tsv", sep="\t", index=False)
    if len(summary):
        summary.to_csv(out / "pair_summary.tsv", sep="\t", index=False)
    if len(sanity):
        sanity.to_csv(out / "distance_sanity.tsv", sep="\t", index=False)
    inventory = {
        "n_qc": int(n_qc) if n_qc else None,
        "n_fov": int(fov.groupby(["sample", "fov"]).ngroups) if len(fov) else 0,
        "n_slides": int(fov["sample"].nunique()) if len(fov) else 0,
        "n_pairs": int(len(pairs)),
        "n_pairs_in_db": int(pairs["in_cellchat_human"].sum()),
        "code_stamp": CODE_STAMP,
        "prox_um": PROX_UM,
        "commot_dis_um": COMMOT_DIS_UM,
        "sdm_l_um": SDM_L_UM,
        "commot_niter": COMMOT_NITER,
        "absent_genes": absent,
        "versions": {
            "commot": ct.__version__ if hasattr(ct, "__version__") else "0.0.3",
            "spatialdm": sdm.__version__ if hasattr(sdm, "__version__") else "0.3.1",
            "squidpy": sq.__version__,
        },
    }
    (out / "stats.json").write_text(json.dumps(inventory, indent=2))
    if len(summary):
        write_results_md(out / "RESULTS.md", summary, sanity, inventory)
        if out.resolve() == OUT.resolve():
            write_results_md(ROOT / "RESULTS.md", summary, sanity, inventory)
    if len(slide) and len(summary):
        plot_heatmap(
            slide,
            summary,
            "d_sdm_prox",
            "SpatialDM proximity: CLDN4-high − CLDN4-low Moran I",
            out / "figures" / "spatialdm_prox_delta_heatmap.png",
        )
        plot_heatmap(
            slide,
            summary,
            "d_commot",
            "COMMOT: CLDN4-high − CLDN4-low transport per sender",
            out / "figures" / "commot_delta_heatmap.png",
        )
        plot_heatmap(
            slide,
            summary,
            "d_squidpy",
            "Squidpy ligrec: CLDN4-high − CLDN4-low interface mean",
            out / "figures" / "squidpy_delta_heatmap.png",
        )
        plot_sign_bars(summary, out / "figures" / "barrier_slides_delta_negative.png")
    log(f"wrote {out}")


if __name__ == "__main__":
    main()
