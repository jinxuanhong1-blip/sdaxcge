#!/usr/bin/env python3
"""CLDN4-only public Visium analysis (10x NSCLC demos + open GEO LUAD).

Pre-specified per section (not tuned after seeing results):
  1. Spearman of log1p-CPM CLDN4 vs CD8A on QC spots.
  2. Epithelial-like spots: epi score (mean of available KRT8/KRT18/KRT19/EPCAM)
     >= section median.
  3. Among epithelial-like spots, CLDN4 quartiles (Q4 vs Q1):
       - Euclidean distance (µm) to nearest CD8A-high spot
         (CD8A-high = CD8A >= 75th percentile of all QC spots; fallback CD8A>0).
       - Hex ring-1 mean neighbor CD8A.
  4. KRT8 residual: OLS residual of CLDN4 on KRT8 (log1p-CPM), then Spearman
     vs CD8A and the same Q4 vs Q1 distance / neighbor tests on the residual.
  5. Spatial maps.

No TACSTD2 dual-high filter. No private 8-KL. Numbers only; no claim language.
"""
from __future__ import annotations

import json
import os
import tarfile
import warnings
from collections import deque
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread
from scipy.sparse import csc_matrix, csr_matrix, issparse
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent.parent  # /workspace
DATA = Path(os.environ.get("VISIUM_DATA", "/workspace/data"))
OUT = ROOT / "results"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

EPI_GENES = ["KRT8", "KRT18", "KRT19", "EPCAM"]
SPOT_UM = 55.0
NN_CENTER_UM = 100.0  # Visium adjacent-spot center-to-center


def decode(x):
    if isinstance(x, bytes):
        return x.decode("utf-8")
    return str(x)


def make_unique(names):
    seen = {}
    out = []
    for n in names:
        if n in seen:
            seen[n] += 1
            out.append(f"{n}.{seen[n]}")
        else:
            seen[n] = 0
            out.append(n)
    return out


def read_10x_h5(path):
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        data = m["data"][:]
        indices = m["indices"][:]
        indptr = m["indptr"][:]
        shape = tuple(int(x) for x in m["shape"][:])
        barcodes = [decode(x) for x in m["barcodes"][:]]
        genes = make_unique([decode(x) for x in m["features"]["name"][:]])
        X = csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
    return X, barcodes, genes


def read_mtx_dir(mtx, features, barcodes):
    X = mmread(mtx).tocsr().T
    feat = pd.read_csv(features, sep="\t", header=None)
    genes = make_unique(feat.iloc[:, 1 if feat.shape[1] > 1 else 0].astype(str).tolist())
    bc = pd.read_csv(barcodes, header=None)[0].astype(str).tolist()
    return X, bc, genes


def read_positions(path):
    # gzip or plain; header or not
    df = pd.read_csv(path, header=None)
    first = str(df.iloc[0, 0]).lower()
    if first in ("barcode", "barcodes"):
        df = pd.read_csv(path)
        mapping = {}
        for c in df.columns:
            cl = c.lower()
            if cl in ("barcode", "barcodes"):
                mapping[c] = "barcode"
            elif "in_tissue" in cl or cl == "in-tissue":
                mapping[c] = "in_tissue"
            elif "array_row" in cl:
                mapping[c] = "array_row"
            elif "array_col" in cl:
                mapping[c] = "array_col"
            elif "pxl_row" in cl:
                mapping[c] = "pxl_row"
            elif "pxl_col" in cl:
                mapping[c] = "pxl_col"
        df = df.rename(columns=mapping)
    else:
        df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"][: df.shape[1]]
    df["barcode"] = df["barcode"].astype(str)
    for c in ("in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("barcode")


def _load_scalefactors(spatial_dir):
    if not spatial_dir:
        return None
    d = Path(spatial_dir)
    cands = [
        d / "scalefactors_json.json",
        d / "scalefactors_json.json.gz",
    ]
    cands += list(d.glob("*scalefactors_json.json*"))
    for p in cands:
        if not p.exists():
            continue
        raw = p.read_bytes()
        if p.name.endswith(".gz") or raw[:2] == b"\x1f\x8b":
            import gzip

            raw = gzip.decompress(raw)
        return json.loads(raw.decode("utf-8"))
    return None


def microns_per_pixel(spatial_dir, pxl):
    sf = _load_scalefactors(spatial_dir)
    if sf:
        d = float(sf.get("spot_diameter_fullres") or 0)
        if d > 0:
            return SPOT_UM / d
    # calibrate so median nearest-neighbor = 100 µm
    tree = cKDTree(pxl)
    nn = tree.query(pxl, k=2)[0][:, 1]
    med = float(np.median(nn[nn > 0])) if np.any(nn > 0) else np.nan
    if not np.isfinite(med) or med <= 0:
        return 1.0
    return NN_CENTER_UM / med


def gene_col(X, genes, name):
    try:
        i = genes.index(name)
    except ValueError:
        return None
    v = X[:, i]
    return np.asarray(v.todense() if issparse(v) else v).ravel().astype(float)


def logcpm(X):
    if issparse(X):
        lib = np.asarray(X.sum(axis=1)).ravel()
        lib[lib == 0] = 1.0
        # avoid densifying whole matrix: we only need a few genes later
        return lib
    lib = X.sum(axis=1)
    lib[lib == 0] = 1.0
    return lib


def log1p_cpm_gene(counts, lib):
    return np.log1p(counts / lib * 1e4)


def spearman(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < 20 or np.std(a[m]) == 0 or np.std(b[m]) == 0:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(a[m], b[m])
    return float(rho), float(p), n


def ols_residual(y, x):
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    m = np.isfinite(y) & np.isfinite(x)
    r = np.full_like(y, np.nan, dtype=float)
    if m.sum() < 10 or np.std(x[m]) == 0:
        return r
    A = np.column_stack([np.ones(m.sum()), x[m]])
    coef, *_ = np.linalg.lstsq(A, y[m], rcond=None)
    r[m] = y[m] - A @ coef
    return r


def mw(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 5 or len(b) < 5:
        return np.nan, np.nan, np.nan, np.nan, len(a), len(b)
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(np.median(a)), float(np.median(b)), float(u), float(p), len(a), len(b)


def hex_adj(rows, cols):
    present = set(zip(rows.astype(int), cols.astype(int)))
    adj = {}
    for r, c in present:
        cand = [
            (r, c - 2),
            (r, c + 2),
            (r - 1, c - 1),
            (r - 1, c + 1),
            (r + 1, c - 1),
            (r + 1, c + 1),
        ]
        adj[(r, c)] = [p for p in cand if p in present]
    return adj


def neighbor_mean(keys, values, adj, min_n=1):
    key2i = {k: i for i, k in enumerate(keys)}
    out = np.full(len(values), np.nan)
    for i, src in enumerate(keys):
        nbrs = adj.get(src, [])
        if len(nbrs) < min_n:
            continue
        idx = [key2i[n] for n in nbrs if n in key2i]
        if not idx:
            continue
        out[i] = float(np.nanmean(values[idx]))
    return out


def quartiles(x, mask):
    """Equal-sized bottom/top 25% among mask, ties broken by stable sort."""
    idx = np.where(mask & np.isfinite(x))[0]
    if len(idx) < 20:
        return None
    xv = x[idx]
    order = np.argsort(xv, kind="mergesort")
    n = len(idx)
    k = max(8, n // 4)
    q1_idx = idx[order[:k]]
    q4_idx = idx[order[-k:]]
    q1m = np.zeros(len(x), dtype=bool)
    q4m = np.zeros(len(x), dtype=bool)
    q1m[q1_idx] = True
    q4m[q4_idx] = True
    return q1m, q4m, float(xv[order[k - 1]]), float(xv[order[-k]])


def analyze_section(name, dataset, histology, source_url, X, barcodes, genes, pos, spatial_dir=None):
    rec = dict(
        section=name,
        dataset=dataset,
        histology=histology,
        source=source_url,
        n_raw=int(X.shape[0]),
        has_CLDN4=int("CLDN4" in genes),
        has_CD8A=int("CD8A" in genes),
        has_KRT8=int("KRT8" in genes),
        status="ok",
    )
    if "CLDN4" not in genes or "CD8A" not in genes:
        rec["status"] = "missing_gene"
        return rec, None

    shared = [b for b in barcodes if b in pos.index]
    if not shared:
        rec["status"] = "no_coord_overlap"
        return rec, None
    b2i = {b: i for i, b in enumerate(barcodes)}
    idx = np.array([b2i[b] for b in shared])
    X = X[idx]
    barcodes = shared
    pos = pos.loc[barcodes]
    if "in_tissue" in pos.columns and (pos["in_tissue"] == 1).any():
        keep = (pos["in_tissue"] == 1).to_numpy()
        X = X[keep]
        pos = pos.loc[keep]
        barcodes = list(pos.index)

    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    n_umi = np.asarray(X.sum(axis=1)).ravel()
    qc = (n_genes >= 200) & (n_umi >= 100)
    if qc.sum() < 80:
        rec["status"] = f"qc_fail_n={int(qc.sum())}"
        rec["n_qc"] = int(qc.sum())
        return rec, None
    X = X[qc]
    pos = pos.loc[qc]
    n_genes = n_genes[qc]
    n_umi = n_umi[qc]
    lib = n_umi.astype(float)
    lib[lib == 0] = 1.0

    def gexp(g):
        c = gene_col(X, genes, g)
        if c is None:
            return np.full(X.shape[0], np.nan)
        return log1p_cpm_gene(c, lib)

    cldn4 = gexp("CLDN4")
    cd8a = gexp("CD8A")
    krt8 = gexp("KRT8")
    epi_mat = []
    epi_used = []
    for g in EPI_GENES:
        v = gexp(g)
        if np.isfinite(v).any() and np.nanstd(v) > 0:
            epi_mat.append(v)
            epi_used.append(g)
    if not epi_mat:
        rec["status"] = "no_epi_genes"
        return rec, None
    epi = np.nanmean(np.vstack(epi_mat), axis=0)
    epi_like = epi >= np.nanmedian(epi)

    pxl = pos[["pxl_col", "pxl_row"]].to_numpy(float)
    mpp = microns_per_pixel(spatial_dir, pxl)
    xy_um = pxl * mpp

    rec.update(
        n_qc=int(X.shape[0]),
        median_genes=float(np.median(n_genes)),
        median_umi=float(np.median(n_umi)),
        n_epi=int(epi_like.sum()),
        epi_genes=",".join(epi_used),
        um_per_px=float(mpp),
        pct_CLDN4_pos=float((cldn4 > 0).mean() * 100),
        pct_CD8A_pos=float((cd8a > 0).mean() * 100),
        pct_KRT8_pos=float((krt8 > 0).mean() * 100) if np.isfinite(krt8).any() else np.nan,
    )

    rho, p, n = spearman(cldn4, cd8a)
    rec.update(rho_CLDN4_CD8A=rho, p_CLDN4_CD8A=p, n_CLDN4_CD8A=n)
    rho, p, n = spearman(cldn4[epi_like], cd8a[epi_like])
    rec.update(rho_CLDN4_CD8A_epi=rho, p_CLDN4_CD8A_epi=p, n_CLDN4_CD8A_epi=n)

    # CD8A-high: if the 75th percentile is 0 (sparse CD8A), use CD8A>0.
    # Otherwise the top quartile. Never label every spot CD8A-high.
    q75 = float(np.nanpercentile(cd8a, 75))
    if q75 > 0:
        cd8_hi = cd8a >= q75
        rec["CD8A_high_rule"] = "q75"
    else:
        cd8_hi = cd8a > 0
        rec["CD8A_high_rule"] = "gt0"
    rec["n_CD8A_high"] = int(cd8_hi.sum())

    dist = np.full(len(cldn4), np.nan)
    if cd8_hi.sum() >= 10:
        xy_hi = xy_um[cd8_hi]
        tree = cKDTree(xy_hi)
        # k=2 so CD8A-high index spots can skip self (distance 0)
        k = 2 if cd8_hi.sum() >= 2 else 1
        d, _ = tree.query(xy_um, k=k)
        d = np.asarray(d, float)
        if k == 1:
            dist = d
        else:
            dist = d[:, 0]
            self_hit = cd8_hi & (d[:, 0] <= 1e-8)
            dist[self_hit] = d[self_hit, 1]
        rec["med_nn_CD8A_high_um"] = float(np.nanmedian(dist))
    else:
        rec["med_nn_CD8A_high_um"] = np.nan

    keys = list(zip(pos["array_row"].astype(int).to_numpy(), pos["array_col"].astype(int).to_numpy()))
    adj = hex_adj(pos["array_row"].to_numpy(), pos["array_col"].to_numpy())
    nb_cd8 = neighbor_mean(keys, cd8a, adj, min_n=1)
    rho, p, n = spearman(cldn4[epi_like], nb_cd8[epi_like])
    rec.update(rho_CLDN4_nbCD8A_epi=rho, p_CLDN4_nbCD8A_epi=p, n_CLDN4_nbCD8A_epi=n)

    q = quartiles(cldn4, epi_like)
    if q is None:
        rec["status"] = "epi_quartile_fail"
        return rec, None
    q1m, q4m, q1v, q4v = q
    rec.update(n_Q1=int(q1m.sum()), n_Q4=int(q4m.sum()), CLDN4_Q1=q1v, CLDN4_Q3=q4v)

    med_q4, med_q1, u, p, n4, n1 = mw(dist[q4m], dist[q1m])
    rec.update(
        med_dist_Q4_um=med_q4,
        med_dist_Q1_um=med_q1,
        delta_dist_Q4_minus_Q1_um=(med_q4 - med_q1) if np.isfinite(med_q4) and np.isfinite(med_q1) else np.nan,
        mw_u_dist=u,
        p_dist_Q4_vs_Q1=p,
    )
    med_q4, med_q1, u, p, n4, n1 = mw(nb_cd8[q4m], nb_cd8[q1m])
    rec.update(
        med_nbCD8A_Q4=med_q4,
        med_nbCD8A_Q1=med_q1,
        delta_nbCD8A_Q4_minus_Q1=(med_q4 - med_q1) if np.isfinite(med_q4) and np.isfinite(med_q1) else np.nan,
        mw_u_nbCD8A=u,
        p_nbCD8A_Q4_vs_Q1=p,
    )

    # KRT8 residual
    resid = ols_residual(cldn4, krt8)
    rho, p, n = spearman(resid, cd8a)
    rec.update(rho_KRT8resid_CD8A=rho, p_KRT8resid_CD8A=p, n_KRT8resid_CD8A=n)
    rho, p, n = spearman(resid[epi_like], cd8a[epi_like])
    rec.update(rho_KRT8resid_CD8A_epi=rho, p_KRT8resid_CD8A_epi=p, n_KRT8resid_CD8A_epi=n)
    rho, p, n = spearman(resid[epi_like], nb_cd8[epi_like])
    rec.update(rho_KRT8resid_nbCD8A_epi=rho, p_KRT8resid_nbCD8A_epi=p)

    rq = quartiles(resid, epi_like)
    if rq is not None:
        r1, r4, _, _ = rq
        rec["n_resid_Q1"] = int(r1.sum())
        rec["n_resid_Q4"] = int(r4.sum())
        med_q4, med_q1, u, p, *_ = mw(dist[r4], dist[r1])
        rec.update(
            med_dist_residQ4_um=med_q4,
            med_dist_residQ1_um=med_q1,
            delta_dist_residQ4_minus_Q1_um=(med_q4 - med_q1)
            if np.isfinite(med_q4) and np.isfinite(med_q1)
            else np.nan,
            p_dist_residQ4_vs_Q1=p,
        )
        med_q4, med_q1, u, p, *_ = mw(nb_cd8[r4], nb_cd8[r1])
        rec.update(
            med_nbCD8A_residQ4=med_q4,
            med_nbCD8A_residQ1=med_q1,
            delta_nbCD8A_residQ4_minus_Q1=(med_q4 - med_q1)
            if np.isfinite(med_q4) and np.isfinite(med_q1)
            else np.nan,
            p_nbCD8A_residQ4_vs_Q1=p,
        )

    maps = dict(
        xy=xy_um,
        cldn4=cldn4,
        cd8a=cd8a,
        krt8=krt8,
        resid=resid,
        epi_like=epi_like,
        q1=q1m,
        q4=q4m,
        cd8_hi=cd8_hi,
        dist=dist,
        nb_cd8=nb_cd8,
    )
    return rec, maps


def plot_maps(name, maps):
    xy = maps["xy"]
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 8.2))
    items = [
        (axes[0, 0], maps["cldn4"], "CLDN4 (log1p CPM)", "viridis"),
        (axes[0, 1], maps["cd8a"], "CD8A (log1p CPM)", "inferno"),
        (axes[0, 2], maps["krt8"], "KRT8 (log1p CPM)", "cividis"),
        (axes[1, 0], maps["resid"], "CLDN4 | KRT8 residual", "coolwarm"),
        (axes[1, 1], maps["nb_cd8"], "ring-1 neighbor CD8A", "inferno"),
        (axes[1, 2], None, "epi Q4 vs Q1 CLDN4", None),
    ]
    for ax, val, title, cmap in items[:-1]:
        sc = ax.scatter(xy[:, 0], xy[:, 1], c=val, s=4, cmap=cmap, linewidths=0)
        ax.set_title(title, fontsize=10)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    ax = axes[1, 2]
    ax.scatter(xy[:, 0], xy[:, 1], c="#dddddd", s=3, linewidths=0, label="_other")
    ax.scatter(xy[maps["q1"], 0], xy[maps["q1"], 1], c="#3b6ea5", s=8, linewidths=0, label="epi Q1 CLDN4")
    ax.scatter(xy[maps["q4"], 0], xy[maps["q4"], 1], c="#c0392b", s=8, linewidths=0, label="epi Q4 CLDN4")
    ax.scatter(
        xy[maps["cd8_hi"], 0],
        xy[maps["cd8_hi"], 1],
        facecolors="none",
        edgecolors="#222222",
        s=10,
        linewidths=0.3,
        label="CD8A-high",
    )
    ax.set_title("epi Q4 vs Q1 CLDN4", fontsize=10)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(loc="best", fontsize=7, frameon=False, markerscale=1.4)
    fig.suptitle(name, fontsize=12)
    fig.tight_layout()
    out = FIG / f"map_{name}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def load_10x_section(h5, spatial_dir):
    X, bc, genes = read_10x_h5(h5)
    pos_path = Path(spatial_dir) / "tissue_positions.csv"
    if not pos_path.exists():
        pos_path = Path(spatial_dir) / "tissue_positions_list.csv"
    pos = read_positions(pos_path)
    return X, bc, genes, pos, spatial_dir


def load_gse189_section(raw_dir, gsm, sid):
    raw_dir = Path(raw_dir)
    prefix = f"{gsm}_{sid}"
    X, bc, genes = read_mtx_dir(
        raw_dir / f"{prefix}_matrix.mtx.gz",
        raw_dir / f"{prefix}_features.tsv.gz",
        raw_dir / f"{prefix}_barcodes.tsv.gz",
    )
    pos = read_positions(raw_dir / f"{prefix}_tissue_positions_list.csv.gz")
    return X, bc, genes, pos, None


def iter_geo_sections_from_tar(tar_path, extract_dir):
    """Yield (sample_id, kind, files) after extracting matrix-like members only."""
    extract_dir = Path(extract_dir)
    extract_dir.mkdir(parents=True, exist_ok=True)
    keep_suffixes = (
        ".h5",
        ".h5ad",
        "matrix.mtx.gz",
        "matrix.mtx",
        "features.tsv.gz",
        "features.tsv",
        "barcodes.tsv.gz",
        "barcodes.tsv",
        "tissue_positions_list.csv.gz",
        "tissue_positions.csv.gz",
        "tissue_positions_list.csv",
        "tissue_positions.csv",
        "scalefactors_json.json.gz",
        "scalefactors_json.json",
        "spatial.tar.gz",
        "spatial.tar",
    )
    with tarfile.open(tar_path) as tf:
        members = [m for m in tf.getmembers() if m.isfile()]
        names = [m.name for m in members]
        wanted = []
        for m in members:
            bn = os.path.basename(m.name).lower()
            if any(bn.endswith(s) or bn.endswith(s.replace(".gz", "")) for s in keep_suffixes):
                if bn.endswith((".tif", ".tiff", ".png", ".jpg", ".ndpi", ".btf", ".cloupe")):
                    continue
                wanted.append(m)
        # always extract small text/matrix; skip huge tifs already filtered
        for m in wanted:
            dest = extract_dir / os.path.basename(m.name)
            if dest.exists() and dest.stat().st_size > 0:
                continue
            src = tf.extractfile(m)
            if src is None:
                continue
            dest.write_bytes(src.read())
    return extract_dir


def discover_geo_sections(extract_dir):
    """Group extracted files into sections."""
    extract_dir = Path(extract_dir)
    files = list(extract_dir.iterdir())
    by_stem = {}
    for p in files:
        name = p.name
        # GSM123_sample_...
        if name.endswith("_filtered_feature_bc_matrix.h5"):
            sid = name.replace("_filtered_feature_bc_matrix.h5", "")
            by_stem.setdefault(sid, {})["h5"] = p
        elif name.endswith(".h5") and "matrix" in name.lower():
            sid = name.replace(".h5", "")
            by_stem.setdefault(sid, {})["h5"] = p
        elif name.endswith("_matrix.mtx.gz") or name.endswith("_matrix.mtx"):
            sid = name.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
            by_stem.setdefault(sid, {})["mtx"] = p
        elif "_features.tsv" in name:
            sid = name.split("_features.tsv")[0]
            by_stem.setdefault(sid, {})["features"] = p
        elif "_barcodes.tsv" in name:
            sid = name.split("_barcodes.tsv")[0]
            by_stem.setdefault(sid, {})["barcodes"] = p
        elif "tissue_positions" in name:
            sid = name.split("_tissue_positions")[0]
            by_stem.setdefault(sid, {})["pos"] = p
        elif name.endswith("_spatial.tar.gz") or name.endswith("-spatial.tar") or name.endswith("_spatial.tar"):
            sid = name.replace("_spatial.tar.gz", "").replace("-spatial.tar", "").replace("_spatial.tar", "")
            by_stem.setdefault(sid, {})["spatial_tar"] = p
        elif name.endswith("scalefactors_json.json") or name.endswith("scalefactors_json.json.gz"):
            sid = name.split("_scalefactors")[0].split("-scalefactors")[0]
            by_stem.setdefault(sid, {})["sf"] = p
    return by_stem


def load_discovered(sid, files, work):
    work = Path(work)
    spatial_dir = None
    if "spatial_tar" in files:
        sdir = work / f"{sid}_spatial"
        sdir.mkdir(exist_ok=True)
        with tarfile.open(files["spatial_tar"]) as tf:
            tf.extractall(sdir)
        # find folder containing positions
        cands = list(sdir.rglob("tissue_positions*.csv")) + list(sdir.rglob("tissue_positions*.csv.gz"))
        if cands:
            spatial_dir = str(cands[0].parent)
    if "h5" in files:
        X, bc, genes = read_10x_h5(files["h5"])
        pos = None
        if spatial_dir:
            for fn in ("tissue_positions.csv", "tissue_positions_list.csv"):
                p = Path(spatial_dir) / fn
                if p.exists():
                    pos = read_positions(p)
                    break
        if pos is None and "pos" in files:
            pos = read_positions(files["pos"])
        if pos is None:
            raise FileNotFoundError(f"no positions for {sid}")
        return X, bc, genes, pos, spatial_dir
    if "mtx" in files and "features" in files and "barcodes" in files and "pos" in files:
        X, bc, genes = read_mtx_dir(files["mtx"], files["features"], files["barcodes"])
        pos = read_positions(files["pos"])
        if spatial_dir is None and "sf" in files:
            sdir = work / f"{sid}_spatial"
            sdir.mkdir(exist_ok=True)
            dest = sdir / files["sf"].name
            if not dest.exists():
                dest.write_bytes(Path(files["sf"]).read_bytes())
            spatial_dir = str(sdir)
        return X, bc, genes, pos, spatial_dir
    raise FileNotFoundError(f"incomplete files for {sid}: {sorted(files)}")


GSE189_META = {
    "TD1": ("IAC", "GSM5702473"),
    "TD2": ("IAC", "GSM5702474"),
    "TD3": ("MIA", "GSM5702475"),
    "TD5": ("AIS", "GSM5702476"),
    "TD6": ("MIA", "GSM5702477"),
    "TD8": ("AIS", "GSM5702478"),
}


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_rho(r):
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def fmt_num(x, nd=1):
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{nd}f}"


def plot_summary(ok: pd.DataFrame):
    if ok.empty:
        return
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.2))
    specs = [
        (axes[0], "rho_CLDN4_CD8A", "Spearman CLDN4 vs CD8A"),
        (axes[1], "delta_dist_Q4_minus_Q1_um", "Δ dist Q4−Q1 (µm)"),
        (axes[2], "delta_nbCD8A_Q4_minus_Q1", "Δ neighbor CD8A Q4−Q1"),
    ]
    colors = {"10x": "#1f4e79", "GSE189487": "#c0392b", "GSE273378": "#1e8449", "GSE300676": "#8e44ad"}
    for ax, col, title in specs:
        y = np.arange(len(ok))
        v = ok[col].astype(float).to_numpy()
        c = []
        for s in ok["dataset"]:
            key = "10x" if str(s).startswith("10x") else str(s)
            c.append(colors.get(key, "#555555"))
        ax.axvline(0, color="#888888", lw=0.8)
        ax.scatter(v, y, c=c, s=18, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(ok["section"].tolist(), fontsize=5)
        ax.set_title(title, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("value")
    fig.tight_layout()
    fig.savefig(FIG / "summary_per_section.png", dpi=140)
    plt.close(fig)


def write_results(rows):
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "per_section.tsv", sep="\t", index=False)
    ok = df[df["status"] == "ok"].copy()
    plot_summary(ok)
    lines = []
    lines.append("# RESULTS: public Visium CLDN4 vs CD8A (additive)")
    lines.append("")
    lines.append("**CLDN4-only. Public data only. No private 8-KL. No claim language.**")
    lines.append("")
    lines.append(f"- Sections attempted: **{len(df)}**")
    lines.append(f"- Sections with numbers (QC pass + CLDN4/CD8A present): **{len(ok)}**")
    if len(ok):
        vc = ok["dataset"].map(lambda s: "10x" if str(s).startswith("10x") else str(s)).value_counts()
        bits = ", ".join(f"{k} n={int(v)}" for k, v in vc.items())
        lines.append(f"- By source: {bits}")
    lines.append("")
    lines.append("## Methods (pre-specified)")
    lines.append("")
    lines.append("- Counts: Space Ranger filtered matrix (in-tissue barcodes).")
    lines.append("- QC: ≥200 genes and ≥100 UMI per spot.")
    lines.append("- Normalization: log1p(CPM to 10,000).")
    lines.append("- Epithelial-like: mean of available `KRT8/KRT18/KRT19/EPCAM` ≥ section median.")
    lines.append("- CD8A-high: CD8A ≥ 75th percentile of QC spots when that cutoff is >0; otherwise CD8A > 0. Distance excludes self.")
    lines.append("- Distance: Euclidean µm from pixel coordinates (55 µm / `spot_diameter_fullres`, or nearest-neighbor calibrated to 100 µm).")
    lines.append("- Neighbor CD8A: mean CD8A of Visium hex ring-1 neighbors.")
    lines.append("- KRT8 residual: OLS residual of log1p-CPM CLDN4 on KRT8.")
    lines.append("- Q4 vs Q1: Mann–Whitney U, two-sided, among epithelial-like spots.")
    lines.append("- Unit: **section** (one capture area).")
    lines.append("")
    lines.append("## Inventory")
    lines.append("")
    lines.append("| section | dataset | histology | n_qc | n_epi | n_Q1 | n_Q4 | n_CD8A-high | status |")
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---|")
    for _, r in df.iterrows():
        lines.append(
            f"| {r['section']} | {r['dataset']} | {r.get('histology','')} | "
            f"{r.get('n_qc','')} | {r.get('n_epi','')} | {r.get('n_Q1','')} | {r.get('n_Q4','')} | "
            f"{r.get('n_CD8A_high','')} | {r['status']} |"
        )
    lines.append("")
    if len(ok):
        lines.append("## Per-section numbers")
        lines.append("")
        lines.append(
            "| section | n_qc | ρ CLDN4–CD8A | p | ρ (epi) | p_epi | "
            "med dist Q4 / Q1 µm | Δdist | p_dist | "
            "med nbCD8A Q4 / Q1 | Δnb | p_nb | "
            "ρ KRT8-resid–CD8A | p_resid |"
        )
        lines.append("|---|---:|---:|---:|---:|---:|---|---:|---:|---|---:|---:|---:|---:|")
        for _, r in ok.iterrows():
            lines.append(
                f"| {r['section']} | {int(r['n_qc'])} | "
                f"{fmt_rho(r.get('rho_CLDN4_CD8A'))} | {fmt_p(r.get('p_CLDN4_CD8A'))} | "
                f"{fmt_rho(r.get('rho_CLDN4_CD8A_epi'))} | {fmt_p(r.get('p_CLDN4_CD8A_epi'))} | "
                f"{fmt_num(r.get('med_dist_Q4_um'))} / {fmt_num(r.get('med_dist_Q1_um'))} | "
                f"{fmt_num(r.get('delta_dist_Q4_minus_Q1_um'))} | {fmt_p(r.get('p_dist_Q4_vs_Q1'))} | "
                f"{fmt_num(r.get('med_nbCD8A_Q4'),3)} / {fmt_num(r.get('med_nbCD8A_Q1'),3)} | "
                f"{fmt_num(r.get('delta_nbCD8A_Q4_minus_Q1'),3)} | {fmt_p(r.get('p_nbCD8A_Q4_vs_Q1'))} | "
                f"{fmt_rho(r.get('rho_KRT8resid_CD8A'))} | {fmt_p(r.get('p_KRT8resid_CD8A'))} |"
            )
        lines.append("")
        lines.append("## Cross-section summary (honest n = number of QC-pass sections)")
        lines.append("")
        n = len(ok)
        def med_iqr(col):
            v = ok[col].astype(float).to_numpy()
            v = v[np.isfinite(v)]
            if len(v) == 0:
                return "NA"
            return f"{np.median(v):+.3f} (IQR {np.percentile(v,25):+.3f} to {np.percentile(v,75):+.3f}; n={len(v)})"

        def sign_test(col):
            v = ok[col].astype(float).to_numpy()
            v = v[np.isfinite(v)]
            if len(v) == 0:
                return "NA"
            pos = int((v > 0).sum())
            neg = int((v < 0).sum())
            # two-sided binomial vs 0.5 on non-zero
            k = min(pos, neg)
            p = stats.binomtest(pos, pos + neg, 0.5, alternative="two-sided").pvalue if (pos + neg) else np.nan
            return f"{pos} positive / {neg} negative / {n - pos - neg} zero; two-sided sign-test p={fmt_p(p)}"

        lines.append(f"- n sections = **{n}**")
        lines.append(f"- Median Spearman CLDN4 vs CD8A (all QC spots): {med_iqr('rho_CLDN4_CD8A')}")
        lines.append(f"- Signs of ρ CLDN4–CD8A: {sign_test('rho_CLDN4_CD8A')}")
        lines.append(f"- Median Spearman CLDN4 vs CD8A (epithelial-like): {med_iqr('rho_CLDN4_CD8A_epi')}")
        lines.append(f"- Median Δ distance Q4−Q1 (µm; positive = Q4 farther from CD8A-high): {med_iqr('delta_dist_Q4_minus_Q1_um')}")
        lines.append(f"- Signs of Δ distance: {sign_test('delta_dist_Q4_minus_Q1_um')}")
        lines.append(f"- Median Δ neighbor CD8A Q4−Q1 (negative = Q4 has lower neighbor CD8A): {med_iqr('delta_nbCD8A_Q4_minus_Q1')}")
        lines.append(f"- Signs of Δ neighbor CD8A: {sign_test('delta_nbCD8A_Q4_minus_Q1')}")
        lines.append(f"- Median Spearman KRT8-residual CLDN4 vs CD8A: {med_iqr('rho_KRT8resid_CD8A')}")
        lines.append(f"- Signs of KRT8-residual ρ: {sign_test('rho_KRT8resid_CD8A')}")
        lines.append("")
        lines.append("Maps: `methods/visium_10x_nsclc_cldn4/results/figures/map_<section>.png`. Summary: `summary_per_section.png`.")
        lines.append("")
    lines.append("## Skipped / not used")
    lines.append("")
    lines.append("- Private 8-KL: not used.")
    lines.append("- JGAS / HUM0394 / EGA controlled: skipped.")
    lines.append("- Visium HD 10x LUAD (fixed frozen) binned_outputs.tar.gz is 8.9 GB; not downloaded here (separate HD hunt).")
    lines.append("- GSE307534 (56 Visium LUAD/precursor, 9.4 GB) and GSE277206: left to dedicated hunts.")
    lines.append("- GSE303162: mouse LLC Visium, not human.")
    lines.append("- GSE301973 / GSE288758: Visium HD (EGFR NSCLC), not standard Visium.")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (ROOT / "RESULTS.md").write_text(text)
    (REPO / "RESULTS.md").write_text(text)
    (OUT / "summary.json").write_text(json.dumps({
        "n_attempted": int(len(df)),
        "n_ok": int(len(ok)),
        "sections": ok["section"].tolist() if len(ok) else [],
    }, indent=2))
    return df


def run():
    rows = []

    def add(name, dataset, histo, url, loader):
        print(f"== {name} ==", flush=True)
        try:
            X, bc, genes, pos, sdir = loader()
        except Exception as e:
            rows.append(dict(section=name, dataset=dataset, histology=histo, source=url, status=f"load_fail:{e}"))
            print("  LOAD FAIL", e, flush=True)
            return
        rec, maps = analyze_section(name, dataset, histo, url, X, bc, genes, pos, sdir)
        rows.append(rec)
        print(" ", rec.get("status"), "n_qc", rec.get("n_qc"), "rho", rec.get("rho_CLDN4_CD8A"), flush=True)
        if maps is not None:
            plot_maps(name, maps)

    # Official 10x Visium (CytAssist FFPE) lung cancer demos
    add(
        "10x_LUSC_FFPE",
        "10x_CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma",
        "LUSC",
        "https://www.10xgenomics.com/datasets/human-lung-cancer-ffpe-2-standard",
        lambda: load_10x_section(
            DATA / "10x_visium/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5",
            DATA / "10x_visium/LUSC_FFPE/spatial",
        ),
    )
    add(
        "10x_NEC_11mm_FFPE",
        "10x_CytAssist_11mm_FFPE_Human_Lung_Cancer",
        "lung_neuroendocrine",
        "https://www.10xgenomics.com/datasets/human-lung-cancer-11-mm-capture-area-ffpe-2-standard",
        lambda: load_10x_section(
            DATA / "10x_visium/CytAssist_11mm_FFPE_Human_Lung_Cancer_filtered_feature_bc_matrix.h5",
            DATA / "10x_visium/NEC_11mm/spatial",
        ),
    )

    # GSE189487 LUAD Visium
    raw189 = DATA / "geo/GSE189487/raw"
    if raw189.exists():
        for sid, (stage, gsm) in GSE189_META.items():
            add(
                f"GSE189487_{sid}_{stage}",
                "GSE189487",
                f"LUAD_{stage}",
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189487",
                lambda gsm=gsm, sid=sid: load_gse189_section(raw189, gsm, sid),
            )

    # Additional GEO tars if present
    for acc, histo_default, url in [
        ("GSE273378", "LUAD", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273378"),
        ("GSE300676", "LUAD", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300676"),
    ]:
        tar = DATA / f"geo/{acc}/{acc}_RAW.tar"
        if not tar.exists() or tar.stat().st_size < 1_000_000:
            print(f"skip {acc}: tar missing/incomplete", flush=True)
            continue
        ex = DATA / f"geo/{acc}/extracted"
        print(f"== extracting {acc} ==", flush=True)
        try:
            iter_geo_sections_from_tar(tar, ex)
            discovered = discover_geo_sections(ex)
        except Exception as e:
            print(" extract fail", acc, e, flush=True)
            continue
        print(f"  discovered {len(discovered)} stems", flush=True)
        for sid, files in sorted(discovered.items()):
            add(
                f"{acc}_{sid}",
                acc,
                histo_default,
                url,
                lambda sid=sid, files=files, ex=ex: load_discovered(sid, files, ex),
            )

    write_results(rows)
    print("wrote RESULTS.md n=", sum(1 for r in rows if r.get("status") == "ok"), flush=True)


if __name__ == "__main__":
    run()
