#!/usr/bin/env python3
"""Spatial-lag CLDN4 exclusion readout on public Visium LUAD/NSCLC.

Same-spot Spearman and nearest-CD8A µm are the wrong punch: a 55 µm Visium
spot mixes epithelium and T cells, so within-spot ρ and nearest-spot distance
collapse toward 0. Primary metrics use *neighbors*, not self.

Pre-specified (not tuned after seeing results)
----------------------------------------------
1. Spatial LAG: among epithelial-like index spots, Spearman of index CLDN4
   vs mean CD8A of hex ring-1 neighbors (~100 µm center-to-center;
   Euclidean 80–150 µm annulus as a paired sensitivity). Expect negative
   if CLDN4-high epithelium excludes nearby CD8.
2. Among epithelial-like spots, CLDN4 Q4 vs Q1: *mean* neighbor CD8A
   (not self). Section-level Δ = mean(Q4) − mean(Q1).
3. Morisita–Horn overlap on a 200 µm grid between CLDN4-high epi (Q4) and
   CD8A-high spots (low = segregation).
4. KRT8 residual *on the lag*: OLS residual of index CLDN4 ~ KRT8, then
   Spearman vs neighbor CD8A (not same-spot CD8A).
5. Drop sections with almost no CD8A or no CLDN4 contrast.
6. Optional: restrict index spots to the tumor–stroma interface
   (epi-like with ≥1 non-epi ring-1 neighbor).

Unit = section. Cross-section test = Wilcoxon signed-rank on section
lag-ρ and on section Δ. CLDN4-only. No private KL. No fabrication.
"""
from __future__ import annotations

import json
import os
import tarfile
from collections import Counter
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread
from scipy.sparse import csc_matrix, issparse
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent.parent
DATA = Path(os.environ.get("VISIUM_DATA", "/workspace/data"))
OUT = ROOT / "results"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

EPI_GENES = ["KRT8", "KRT18", "KRT19", "EPCAM"]
SPOT_UM = 55.0
NN_CENTER_UM = 100.0
ANNULUS_LO = 80.0
ANNULUS_HI = 150.0
GRID_UM = 200.0

# Contrast filters (pre-specified). Sections failing these have no exclusion test.
MIN_QC = 80
MIN_CD8A_POS = 50
MIN_CD8A_PCT = 5.0
MIN_CLDN4_POS_EPI_PCT = 10.0
MIN_EPI_WITH_NBR = 80
MIN_Q = 20


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
    path = Path(path)
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
    cands = [d / "scalefactors_json.json", d / "scalefactors_json.json.gz"]
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


def spearman_ci(rho, n, alpha=0.05):
    if n < 4 or not np.isfinite(rho) or abs(rho) >= 1:
        return np.nan, np.nan
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3)
    zc = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - zc * se)), float(np.tanh(z + zc * se))


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


def annulus_mean(xy, values, rmin=ANNULUS_LO, rmax=ANNULUS_HI):
    tree = cKDTree(xy)
    out = np.full(len(values), np.nan)
    for i, nbrs in enumerate(tree.query_ball_point(xy, rmax)):
        if len(nbrs) <= 1:
            continue
        dxy = xy[nbrs] - xy[i]
        d = np.sqrt((dxy ** 2).sum(axis=1))
        keep = [j for j, dist in zip(nbrs, d) if rmin < dist <= rmax]
        if keep:
            out[i] = float(np.nanmean(values[keep]))
    return out


def quartiles(x, mask):
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


def mean_delta_ci(a, b):
    """Δ = mean(a) − mean(b) with Welch SE / 1.96 CI."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 5 or len(b) < 5:
        return np.nan, np.nan, np.nan, len(a), len(b)
    da, db = float(np.mean(a)), float(np.mean(b))
    delta = da - db
    va = float(np.var(a, ddof=1)) if len(a) > 1 else 0.0
    vb = float(np.var(b, ddof=1)) if len(b) > 1 else 0.0
    se = np.sqrt(va / len(a) + vb / len(b))
    if not np.isfinite(se):
        return delta, np.nan, np.nan, len(a), len(b)
    return delta, delta - 1.96 * se, delta + 1.96 * se, len(a), len(b)


def morisita_horn(xy, mask_a, mask_b, bin_um=GRID_UM):
    xy = np.asarray(xy, float)
    if mask_a.sum() < 5 or mask_b.sum() < 5:
        return np.nan
    origin = xy.min(axis=0)
    ix = np.floor((xy[:, 0] - origin[0]) / bin_um).astype(int)
    iy = np.floor((xy[:, 1] - origin[1]) / bin_um).astype(int)
    keys = list(zip(ix, iy))
    ca = Counter(k for k, m in zip(keys, mask_a) if m)
    cb = Counter(k for k, m in zip(keys, mask_b) if m)
    bins = set(ca) | set(cb)
    if not bins:
        return np.nan
    a = np.array([ca[b] for b in bins], float)
    b = np.array([cb[b] for b in bins], float)
    A, B = a.sum(), b.sum()
    if A <= 0 or B <= 0:
        return np.nan
    lam_a = np.sum(a ** 2) / (A ** 2)
    lam_b = np.sum(b ** 2) / (B ** 2)
    denom = (lam_a + lam_b) * A * B
    if denom <= 0:
        return np.nan
    return float(2.0 * np.sum(a * b) / denom)


def interface_mask(keys, epi_like, adj):
    key2i = {k: i for i, k in enumerate(keys)}
    out = np.zeros(len(keys), dtype=bool)
    for i, src in enumerate(keys):
        if not epi_like[i]:
            continue
        for nk in adj.get(src, []):
            j = key2i.get(nk)
            if j is not None and not epi_like[j]:
                out[i] = True
                break
    return out


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
        contrast_ok=0,
    )
    if "CLDN4" not in genes or "CD8A" not in genes:
        rec["status"] = "missing_gene"
        return rec

    shared = [b for b in barcodes if b in pos.index]
    if not shared:
        rec["status"] = "no_coord_overlap"
        return rec
    b2i = {b: i for i, b in enumerate(barcodes)}
    idx = np.array([b2i[b] for b in shared])
    X = X[idx]
    pos = pos.loc[shared]
    if "in_tissue" in pos.columns and (pos["in_tissue"] == 1).any():
        keep = (pos["in_tissue"] == 1).to_numpy()
        X = X[keep]
        pos = pos.loc[keep]

    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    n_umi = np.asarray(X.sum(axis=1)).ravel()
    qc = (n_genes >= 200) & (n_umi >= 100)
    rec["n_qc"] = int(qc.sum())
    if qc.sum() < MIN_QC:
        rec["status"] = f"qc_fail_n={int(qc.sum())}"
        return rec
    X = X[qc]
    pos = pos.loc[qc]
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
    epi_mat, epi_used = [], []
    for g in EPI_GENES:
        v = gexp(g)
        if np.isfinite(v).any() and np.nanstd(v) > 0:
            epi_mat.append(v)
            epi_used.append(g)
    if not epi_mat:
        rec["status"] = "no_epi_genes"
        return rec
    epi = np.nanmean(np.vstack(epi_mat), axis=0)
    epi_like = epi >= np.nanmedian(epi)

    pxl = pos[["pxl_col", "pxl_row"]].to_numpy(float)
    mpp = microns_per_pixel(spatial_dir, pxl)
    xy_um = pxl * mpp

    rec.update(
        median_umi=float(np.median(n_umi)),
        n_epi=int(epi_like.sum()),
        epi_genes=",".join(epi_used),
        um_per_px=float(mpp),
        n_CLDN4_pos=int((cldn4 > 0).sum()),
        n_CD8A_pos=int((cd8a > 0).sum()),
        pct_CLDN4_pos=float((cldn4 > 0).mean() * 100),
        pct_CD8A_pos=float((cd8a > 0).mean() * 100),
        pct_CLDN4_pos_epi=float((cldn4[epi_like] > 0).mean() * 100) if epi_like.any() else np.nan,
        iqr_CLDN4_epi=float(np.subtract(*np.nanpercentile(cldn4[epi_like], [75, 25])))
        if epi_like.sum() >= 8
        else np.nan,
    )

    q75 = float(np.nanpercentile(cd8a, 75))
    if q75 > 0:
        cd8_hi = cd8a >= q75
        rec["CD8A_high_rule"] = "q75"
    else:
        cd8_hi = cd8a > 0
        rec["CD8A_high_rule"] = "gt0"
    rec["n_CD8A_high"] = int(cd8_hi.sum())

    keys = list(zip(pos["array_row"].astype(int).to_numpy(), pos["array_col"].astype(int).to_numpy()))
    adj = hex_adj(pos["array_row"].to_numpy(), pos["array_col"].to_numpy())
    nb_cd8 = neighbor_mean(keys, cd8a, adj, min_n=1)
    nb_cd8_ann = annulus_mean(xy_um, cd8a)

    epi_nbr = epi_like & np.isfinite(nb_cd8)
    rec["n_epi_with_nbr"] = int(epi_nbr.sum())
    rec["sd_nbCD8A_epi"] = float(np.nanstd(nb_cd8[epi_nbr])) if epi_nbr.sum() else np.nan

    # Contrast filter
    reasons = []
    if rec["n_CD8A_pos"] < MIN_CD8A_POS or rec["pct_CD8A_pos"] < MIN_CD8A_PCT:
        reasons.append("no_CD8A_contrast")
    if (
        not np.isfinite(rec["pct_CLDN4_pos_epi"])
        or rec["pct_CLDN4_pos_epi"] < MIN_CLDN4_POS_EPI_PCT
        or not np.isfinite(rec["iqr_CLDN4_epi"])
        or rec["iqr_CLDN4_epi"] <= 0
    ):
        reasons.append("no_CLDN4_contrast")
    if rec["n_epi_with_nbr"] < MIN_EPI_WITH_NBR:
        reasons.append("too_few_epi_neighbors")
    if not np.isfinite(rec["sd_nbCD8A_epi"]) or rec["sd_nbCD8A_epi"] <= 0:
        reasons.append("neighbor_CD8A_flat")
    if reasons:
        rec["status"] = "no_contrast:" + ",".join(reasons)
        rec["contrast_ok"] = 0
        return rec
    rec["contrast_ok"] = 1

    # Primary: spatial lag (ring-1)
    rho, p, n = spearman(cldn4[epi_nbr], nb_cd8[epi_nbr])
    lo, hi = spearman_ci(rho, n)
    rec.update(
        lag_rho_ring1=rho,
        lag_p_ring1=p,
        lag_n_ring1=n,
        lag_rho_ring1_lo=lo,
        lag_rho_ring1_hi=hi,
    )

    epi_ann = epi_like & np.isfinite(nb_cd8_ann)
    rho, p, n = spearman(cldn4[epi_ann], nb_cd8_ann[epi_ann])
    lo, hi = spearman_ci(rho, n)
    rec.update(
        lag_rho_annulus=rho,
        lag_p_annulus=p,
        lag_n_annulus=n,
        lag_rho_annulus_lo=lo,
        lag_rho_annulus_hi=hi,
    )

    q = quartiles(cldn4, epi_nbr)
    if q is None:
        rec["status"] = "epi_quartile_fail"
        rec["contrast_ok"] = 0
        return rec
    q1m, q4m, q1v, q4v = q
    rec.update(n_Q1=int(q1m.sum()), n_Q4=int(q4m.sum()), CLDN4_Q1=q1v, CLDN4_Q4=q4v)
    if rec["n_Q1"] < MIN_Q or rec["n_Q4"] < MIN_Q:
        rec["status"] = "quartile_too_small"
        rec["contrast_ok"] = 0
        return rec

    dlt, lo, hi, n4, n1 = mean_delta_ci(nb_cd8[q4m], nb_cd8[q1m])
    rec.update(
        mean_nbCD8A_Q4=float(np.nanmean(nb_cd8[q4m])),
        mean_nbCD8A_Q1=float(np.nanmean(nb_cd8[q1m])),
        delta_nbCD8A_Q4_minus_Q1=dlt,
        delta_nbCD8A_lo=lo,
        delta_nbCD8A_hi=hi,
        n_delta_Q4=n4,
        n_delta_Q1=n1,
    )
    # Mann–Whitney as a within-section companion (not the punch)
    a = nb_cd8[q4m]
    b = nb_cd8[q1m]
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) >= 5 and len(b) >= 5:
        rec["p_nbCD8A_Q4_vs_Q1"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    else:
        rec["p_nbCD8A_Q4_vs_Q1"] = np.nan

    # Mixing / Morisita–Horn (CLDN4-high epi vs CD8A-high)
    rec["morisita_horn_Q4epi_CD8hi"] = morisita_horn(xy_um, q4m, cd8_hi, GRID_UM)
    key2i = {k: i for i, k in enumerate(keys)}
    hits = []
    for i, src in enumerate(keys):
        if not q4m[i]:
            continue
        nbrs = adj.get(src, [])
        hits.append(any(cd8_hi[key2i[nk]] for nk in nbrs if nk in key2i))
    rec["frac_Q4_with_CD8hi_nbr"] = float(np.mean(hits)) if hits else np.nan

    # KRT8 residual on the LAG
    resid = ols_residual(cldn4, krt8)
    rho, p, n = spearman(resid[epi_nbr], nb_cd8[epi_nbr])
    lo, hi = spearman_ci(rho, n)
    rec.update(
        lag_rho_KRT8resid_ring1=rho,
        lag_p_KRT8resid_ring1=p,
        lag_n_KRT8resid_ring1=n,
        lag_rho_KRT8resid_lo=lo,
        lag_rho_KRT8resid_hi=hi,
    )
    rq = quartiles(resid, epi_nbr)
    if rq is not None:
        r1, r4, _, _ = rq
        dlt, lo, hi, *_ = mean_delta_ci(nb_cd8[r4], nb_cd8[r1])
        rec.update(
            delta_nbCD8A_residQ4_minus_Q1=dlt,
            delta_nbCD8A_resid_lo=lo,
            delta_nbCD8A_resid_hi=hi,
        )

    # Optional interface
    iface = interface_mask(keys, epi_like, adj) & np.isfinite(nb_cd8)
    rec["n_interface_epi"] = int(iface.sum())
    if iface.sum() >= MIN_EPI_WITH_NBR and np.nanstd(nb_cd8[iface]) > 0:
        rho, p, n = spearman(cldn4[iface], nb_cd8[iface])
        lo, hi = spearman_ci(rho, n)
        rec.update(
            lag_rho_interface=rho,
            lag_p_interface=p,
            lag_n_interface=n,
            lag_rho_interface_lo=lo,
            lag_rho_interface_hi=hi,
        )
        iq = quartiles(cldn4, iface)
        if iq is not None:
            i1, i4, _, _ = iq
            dlt, lo, hi, *_ = mean_delta_ci(nb_cd8[i4], nb_cd8[i1])
            rec.update(
                delta_nbCD8A_interface_Q4_minus_Q1=dlt,
                delta_nbCD8A_interface_lo=lo,
                delta_nbCD8A_interface_hi=hi,
                n_interface_Q1=int(i1.sum()),
                n_interface_Q4=int(i4.sum()),
            )
    else:
        rec["lag_rho_interface"] = np.nan
        rec["delta_nbCD8A_interface_Q4_minus_Q1"] = np.nan

    return rec


def load_10x_section(h5, spatial_dir):
    X, bc, genes = read_10x_h5(h5)
    pos_path = Path(spatial_dir) / "tissue_positions.csv"
    if not pos_path.exists():
        pos_path = Path(spatial_dir) / "tissue_positions_list.csv"
    pos = read_positions(pos_path)
    return X, bc, genes, pos, spatial_dir


def load_mtx_section(raw_dir, prefix, sf_name=None):
    raw_dir = Path(raw_dir)
    X, bc, genes = read_mtx_dir(
        raw_dir / f"{prefix}_matrix.mtx.gz",
        raw_dir / f"{prefix}_features.tsv.gz",
        raw_dir / f"{prefix}_barcodes.tsv.gz",
    )
    pos = read_positions(raw_dir / f"{prefix}_tissue_positions_list.csv.gz")
    sdir = None
    if sf_name and (raw_dir / sf_name).exists():
        sdir_p = raw_dir / f"{prefix}_spatial"
        sdir_p.mkdir(exist_ok=True)
        dest = sdir_p / "scalefactors_json.json.gz"
        if not dest.exists():
            dest.write_bytes((raw_dir / sf_name).read_bytes())
        sdir = str(sdir_p)
    return X, bc, genes, pos, sdir


def load_gse307_section(out_dir):
    out_dir = Path(out_dir)
    mtx = out_dir / "filtered_feature_bc_matrix"
    X, bc, genes = read_mtx_dir(mtx / "matrix.mtx.gz", mtx / "features.tsv.gz", mtx / "barcodes.tsv.gz")
    sdir = out_dir / "spatial"
    pos_path = sdir / "tissue_positions.csv"
    if not pos_path.exists():
        pos_path = sdir / "tissue_positions_list.csv"
    pos = read_positions(pos_path)
    return X, bc, genes, pos, str(sdir)


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


def fmt_num(x, nd=3):
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{nd}f}"


def wilcoxon_vs_zero(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) < 6:
        return dict(n=len(v), n_neg=int((v < 0).sum()), n_pos=int((v > 0).sum()), n_zero=int((v == 0).sum()), p=np.nan, median=float(np.median(v)) if len(v) else np.nan)
    try:
        p = float(stats.wilcoxon(v, alternative="two-sided", zero_method="pratt").pvalue)
    except ValueError:
        p = np.nan
    return dict(
        n=len(v),
        n_neg=int((v < 0).sum()),
        n_pos=int((v > 0).sum()),
        n_zero=int((v == 0).sum()),
        p=p,
        median=float(np.median(v)),
        iqr_lo=float(np.percentile(v, 25)),
        iqr_hi=float(np.percentile(v, 75)),
    )


def forest(ax, labels, est, lo, hi, colors, title, xlabel):
    y = np.arange(len(labels))
    ax.axvline(0, color="#888888", lw=0.8)
    for i, (e, a, b, c) in enumerate(zip(est, lo, hi, colors)):
        if np.isfinite(a) and np.isfinite(b):
            ax.plot([a, b], [i, i], color=c, lw=1.1)
        if np.isfinite(e):
            ax.scatter([e], [i], color=c, s=16, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.invert_yaxis()


def plot_forests(ok: pd.DataFrame):
    if ok.empty:
        return
    colors_map = {
        "10x": "#1f4e79",
        "GSE189487": "#c0392b",
        "GSE273378": "#1e8449",
        "GSE300676": "#8e44ad",
        "GSE307534": "#d35400",
    }
    ds = ok["dataset"].map(lambda s: "10x" if str(s).startswith("10x") else str(s))
    cols = [colors_map.get(s, "#555555") for s in ds]
    labels = ok["section"].tolist()

    fig, axes = plt.subplots(1, 2, figsize=(11.5, max(4.2, 0.18 * len(ok) + 1.4)))
    forest(
        axes[0],
        labels,
        ok["lag_rho_ring1"].astype(float).to_numpy(),
        ok["lag_rho_ring1_lo"].astype(float).to_numpy(),
        ok["lag_rho_ring1_hi"].astype(float).to_numpy(),
        cols,
        "Lag-ρ  (index CLDN4 vs ring-1 CD8A)",
        "Spearman ρ (epi-like index spots)",
    )
    forest(
        axes[1],
        labels,
        ok["delta_nbCD8A_Q4_minus_Q1"].astype(float).to_numpy(),
        ok["delta_nbCD8A_lo"].astype(float).to_numpy(),
        ok["delta_nbCD8A_hi"].astype(float).to_numpy(),
        cols,
        "Q4−Q1 mean neighbor CD8A",
        "Δ mean ring-1 CD8A (log1p CPM)",
    )
    fig.tight_layout()
    fig.savefig(FIG / "forest_lag_and_delta.png", dpi=150)
    plt.close(fig)

    # compact strip by dataset
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.0))
    specs = [
        (axes[0], "lag_rho_ring1", "Lag-ρ ring-1"),
        (axes[1], "delta_nbCD8A_Q4_minus_Q1", "Δ neighbor CD8A Q4−Q1"),
        (axes[2], "morisita_horn_Q4epi_CD8hi", "Morisita–Horn (Q4 epi ∩ CD8-high)"),
    ]
    for ax, col, title in specs:
        ax.axvline(0 if col != "morisita_horn_Q4epi_CD8hi" else 0.5, color="#888", lw=0.8)
        y = np.arange(len(ok))
        ax.scatter(ok[col].astype(float), y, c=cols, s=16)
        ax.set_yticks([])
        ax.set_title(title, fontsize=9)
        ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIG / "section_strips.png", dpi=140)
    plt.close(fig)


def write_results(rows):
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "per_section.tsv", sep="\t", index=False)
    ok = df[df["contrast_ok"] == 1].copy()
    dropped = df[df["contrast_ok"] != 1].copy()
    plot_forests(ok)

    lines = []
    lines.append("# RESULTS: Visium spatial lag (CLDN4 vs neighbor CD8A)")
    lines.append("")
    lines.append("**CLDN4-only. Public Visium LUAD/NSCLC. No private KL. No claim-failed. No fabrication.**")
    lines.append("")
    lines.append("Primary readout is **spatial lag**, not same-spot Spearman and not nearest-spot µm.")
    lines.append("A 55 µm Visium spot mixes tumor and T cells; within-spot ρ and nearest-spot distance")
    lines.append("are the wrong exclusion test and are not reported as the punch.")
    lines.append("")
    lines.append(f"- Sections attempted: **{len(df)}**")
    lines.append(f"- Sections kept (QC + CLDN4/CD8A contrast): **{len(ok)}**")
    lines.append(f"- Sections dropped (no contrast / load / QC): **{len(dropped)}**")
    if len(ok):
        vc = ok["dataset"].map(lambda s: "10x" if str(s).startswith("10x") else str(s)).value_counts()
        bits = ", ".join(f"{k} n={int(v)}" for k, v in vc.items())
        lines.append(f"- Kept by source: {bits}")
    lines.append("")
    lines.append("## Methods (pre-specified)")
    lines.append("")
    lines.append("- Counts: Space Ranger filtered matrix; in-tissue barcodes.")
    lines.append("- QC: ≥200 genes and ≥100 UMI.")
    lines.append("- Normalization: log1p(CPM to 10,000).")
    lines.append("- Epithelial-like: mean of available `KRT8/KRT18/KRT19/EPCAM` ≥ section median.")
    lines.append("- Neighbors: Visium hex ring-1 (self excluded). Sensitivity: Euclidean 80–150 µm annulus.")
    lines.append("- Lag-ρ: Spearman(index CLDN4, mean neighbor CD8A) on epithelial-like spots with ≥1 neighbor.")
    lines.append("- Q4 vs Q1: among those spots, equal-sized bottom/top 25% of CLDN4; Δ = mean(neighbor CD8A | Q4) − mean(neighbor CD8A | Q1).")
    lines.append("- Morisita–Horn: 200 µm grid counts of CLDN4-Q4 epi vs CD8A-high spots (q75, or CD8A>0 if q75=0).")
    lines.append("- KRT8 residual: OLS residual of index CLDN4 on KRT8, then Spearman vs *neighbor* CD8A.")
    lines.append("- Contrast drop: CD8A+ spots <50 or <5% of QC, **or** CLDN4+ among epi <10% or IQR=0, **or** <80 epi spots with a neighbor, **or** neighbor-CD8A SD=0.")
    lines.append("- Optional interface: epi-like spots with ≥1 non-epi ring-1 neighbor.")
    lines.append("- Cross-section: Wilcoxon signed-rank vs 0 on section lag-ρ and on section Δ.")
    lines.append("- Unit: **section** (one capture area).")
    lines.append("")

    if len(ok):
        w_lag = wilcoxon_vs_zero(ok["lag_rho_ring1"])
        w_dlt = wilcoxon_vs_zero(ok["delta_nbCD8A_Q4_minus_Q1"])
        w_ann = wilcoxon_vs_zero(ok["lag_rho_annulus"])
        w_res = wilcoxon_vs_zero(ok["lag_rho_KRT8resid_ring1"])
        w_mh = wilcoxon_vs_zero(ok["morisita_horn_Q4epi_CD8hi"] - 0.5)
        w_iface = wilcoxon_vs_zero(ok["lag_rho_interface"]) if "lag_rho_interface" in ok else None
        w_idlt = (
            wilcoxon_vs_zero(ok["delta_nbCD8A_interface_Q4_minus_Q1"])
            if "delta_nbCD8A_interface_Q4_minus_Q1" in ok
            else None
        )

        lines.append("## Primary: section-level spatial lag")
        lines.append("")
        lines.append(f"- n sections = **{w_lag['n']}**")
        lines.append(
            f"- Lag-ρ (ring-1): median {fmt_rho(w_lag['median'])} "
            f"(IQR {fmt_rho(w_lag.get('iqr_lo'))} to {fmt_rho(w_lag.get('iqr_hi'))})"
        )
        lines.append(
            f"- Signs: **{w_lag['n_neg']} negative / {w_lag['n_pos']} positive** / {w_lag['n_zero']} zero"
        )
        lines.append(f"- Wilcoxon signed-rank on section lag-ρ vs 0: p={fmt_p(w_lag['p'])}")
        lines.append("")
        lines.append(
            f"- Q4−Q1 Δ mean neighbor CD8A: median {fmt_num(w_dlt['median'])} "
            f"(IQR {fmt_num(w_dlt.get('iqr_lo'))} to {fmt_num(w_dlt.get('iqr_hi'))})"
        )
        lines.append(
            f"- Signs of Δ (negative = Q4 has lower neighbor CD8A): "
            f"**{w_dlt['n_neg']} negative / {w_dlt['n_pos']} positive** / {w_dlt['n_zero']} zero"
        )
        lines.append(f"- Wilcoxon signed-rank on section Δ vs 0: p={fmt_p(w_dlt['p'])}")
        lines.append("")
        lines.append("## Mixing (Morisita–Horn)")
        lines.append("")
        mh = ok["morisita_horn_Q4epi_CD8hi"].astype(float)
        mh = mh[np.isfinite(mh)]
        lines.append(
            f"- Morisita–Horn (200 µm grid; 0=segregated, 1=mixed): "
            f"median {fmt_num(float(np.median(mh)) if len(mh) else np.nan)} "
            f"(IQR {fmt_num(float(np.percentile(mh,25)) if len(mh) else np.nan)} to "
            f"{fmt_num(float(np.percentile(mh,75)) if len(mh) else np.nan)}; n={len(mh)})"
        )
        lines.append(
            f"- Wilcoxon of (Morisita−0.5) vs 0 (negative = more segregated than 0.5): p={fmt_p(w_mh['p'])} "
            f"({w_mh['n_neg']} below 0.5 / {w_mh['n_pos']} above)"
        )
        lines.append("")
        lines.append("## KRT8 residual on the lag (not same-spot)")
        lines.append("")
        lines.append(
            f"- Lag-ρ of (CLDN4 | KRT8 residual) vs ring-1 CD8A: median {fmt_rho(w_res['median'])} "
            f"(IQR {fmt_rho(w_res.get('iqr_lo'))} to {fmt_rho(w_res.get('iqr_hi'))})"
        )
        lines.append(
            f"- Signs: **{w_res['n_neg']} negative / {w_res['n_pos']} positive**; Wilcoxon p={fmt_p(w_res['p'])}"
        )
        lines.append("")
        lines.append("## Sensitivity")
        lines.append("")
        lines.append(
            f"- 80–150 µm annulus lag-ρ: median {fmt_rho(w_ann['median'])}; "
            f"{w_ann['n_neg']} neg / {w_ann['n_pos']} pos; Wilcoxon p={fmt_p(w_ann['p'])}"
        )
        if w_iface and w_iface["n"] >= 6:
            lines.append(
                f"- Interface-only lag-ρ: n={w_iface['n']}, median {fmt_rho(w_iface['median'])}; "
                f"{w_iface['n_neg']} neg / {w_iface['n_pos']} pos; Wilcoxon p={fmt_p(w_iface['p'])}"
            )
        if w_idlt and w_idlt["n"] >= 6:
            lines.append(
                f"- Interface Q4−Q1 Δ neighbor CD8A: median {fmt_num(w_idlt['median'])}; "
                f"{w_idlt['n_neg']} neg / {w_idlt['n_pos']} pos; Wilcoxon p={fmt_p(w_idlt['p'])}"
            )
        lines.append("")
        lines.append("Forest: `methods/visium_spatial_lag_cldn4/results/figures/forest_lag_and_delta.png`.")
        lines.append("")
        lines.append("## Per-section lag table")
        lines.append("")
        lines.append(
            "| section | dataset | histo | n_qc | n_epi+nbr | "
            "lag-ρ ring-1 (95% CI) | p | "
            "mean nbCD8A Q4 / Q1 | Δ | "
            "Morisita | resid lag-ρ |"
        )
        lines.append("|---|---|---|---:|---:|---|---:|---|---:|---:|---:|")
        for _, r in ok.iterrows():
            lines.append(
                f"| {r['section']} | {r['dataset']} | {r.get('histology','')} | "
                f"{int(r['n_qc'])} | {int(r['n_epi_with_nbr'])} | "
                f"{fmt_rho(r.get('lag_rho_ring1'))} "
                f"({fmt_rho(r.get('lag_rho_ring1_lo'))} to {fmt_rho(r.get('lag_rho_ring1_hi'))}) | "
                f"{fmt_p(r.get('lag_p_ring1'))} | "
                f"{fmt_num(r.get('mean_nbCD8A_Q4'))} / {fmt_num(r.get('mean_nbCD8A_Q1'))} | "
                f"{fmt_num(r.get('delta_nbCD8A_Q4_minus_Q1'))} | "
                f"{fmt_num(r.get('morisita_horn_Q4epi_CD8hi'))} | "
                f"{fmt_rho(r.get('lag_rho_KRT8resid_ring1'))} |"
            )
        lines.append("")

    lines.append("## Dropped / not used")
    lines.append("")
    if len(dropped):
        for _, r in dropped.iterrows():
            lines.append(f"- `{r['section']}` ({r.get('dataset','')}): {r.get('status','')}")
    lines.append("- Private 8-KL: not used.")
    lines.append("- Same-spot Spearman and nearest-spot µm: computed nowhere as a result; they are the discarded readout.")
    lines.append("- GSE300676 sample IDs are `CRC_mPAP*` as deposited; the GEO series is LUAD micropapillary Visium.")
    lines.append("- GSE307534 normal lung sections were not requested as LUAD/NSCLC and were not downloaded.")
    lines.append("")

    text = "\n".join(lines) + "\n"
    (ROOT / "RESULTS.md").write_text(text)
    (REPO / "RESULTS.md").write_text(text)
    summary = {
        "n_attempted": int(len(df)),
        "n_kept": int(len(ok)),
        "n_dropped": int(len(dropped)),
        "sections_kept": ok["section"].tolist() if len(ok) else [],
    }
    if len(ok):
        summary["lag_ring1"] = wilcoxon_vs_zero(ok["lag_rho_ring1"])
        summary["delta_nbCD8A"] = wilcoxon_vs_zero(ok["delta_nbCD8A_Q4_minus_Q1"])
        summary["lag_krt8resid"] = wilcoxon_vs_zero(ok["lag_rho_KRT8resid_ring1"])
        mh = ok["morisita_horn_Q4epi_CD8hi"].astype(float).to_numpy()
        mh = mh[np.isfinite(mh)]
        summary["morisita_median"] = float(np.median(mh)) if len(mh) else None
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    return df


def run():
    rows = []

    def add(name, dataset, histo, url, loader):
        print(f"== {name} ==", flush=True)
        try:
            X, bc, genes, pos, sdir = loader()
        except Exception as e:
            rows.append(
                dict(
                    section=name,
                    dataset=dataset,
                    histology=histo,
                    source=url,
                    status=f"load_fail:{e}",
                    contrast_ok=0,
                )
            )
            print("  LOAD FAIL", e, flush=True)
            return
        rec = analyze_section(name, dataset, histo, url, X, bc, genes, pos, sdir)
        rows.append(rec)
        print(
            " ",
            rec.get("status"),
            "lag",
            rec.get("lag_rho_ring1"),
            "dlt",
            rec.get("delta_nbCD8A_Q4_minus_Q1"),
            flush=True,
        )

    # 10x demos
    h5_lusc = DATA / "10x_visium/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5"
    s_lusc = DATA / "10x_visium/LUSC_FFPE/spatial"
    if h5_lusc.exists() and s_lusc.exists():
        add(
            "10x_LUSC_FFPE",
            "10x_CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma",
            "LUSC",
            "https://www.10xgenomics.com/datasets/human-lung-cancer-ffpe-2-standard",
            lambda: load_10x_section(h5_lusc, s_lusc),
        )
    h5_nec = DATA / "10x_visium/CytAssist_11mm_FFPE_Human_Lung_Cancer_filtered_feature_bc_matrix.h5"
    s_nec = DATA / "10x_visium/NEC_11mm/spatial"
    if h5_nec.exists() and s_nec.exists():
        add(
            "10x_NEC_11mm_FFPE",
            "10x_CytAssist_11mm_FFPE_Human_Lung_Cancer",
            "lung_neuroendocrine",
            "https://www.10xgenomics.com/datasets/human-lung-cancer-11-mm-capture-area-ffpe-2-standard",
            lambda: load_10x_section(h5_nec, s_nec),
        )

    # GSE189487
    raw189 = DATA / "geo/GSE189487"
    for gsm, sid, stage in [
        ("GSM5702473", "TD1", "IAC"),
        ("GSM5702474", "TD2", "IAC"),
        ("GSM5702475", "TD3", "MIA"),
        ("GSM5702476", "TD5", "AIS"),
        ("GSM5702477", "TD6", "MIA"),
        ("GSM5702478", "TD8", "AIS"),
    ]:
        prefix = f"{gsm}_{sid}"
        if (raw189 / f"{prefix}_matrix.mtx.gz").exists():
            add(
                f"GSE189487_{sid}_{stage}",
                "GSE189487",
                f"LUAD_{stage}",
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189487",
                lambda prefix=prefix: load_mtx_section(raw189, prefix),
            )

    # GSE273378
    raw273 = DATA / "geo/GSE273378"
    for gsm, sid in [
        ("GSM8427428", "LM_SD_1216_1"),
        ("GSM8427429", "LM_SD_16"),
        ("GSM8427430", "LM_SD_11"),
        ("GSM8427431", "LM_SD_2"),
        ("GSM8427432", "LM_SD_3"),
        ("GSM8427433", "LM_SD_4"),
        ("GSM8427434", "LM_SD_5"),
        ("GSM8427435", "LM_SD_6"),
        ("GSM8427436", "LM_SD_7"),
        ("GSM8427437", "LM_SD_1216_8"),
        ("GSM8427438", "LM_SD_9"),
        ("GSM8427439", "LM_SD_10"),
        ("GSM8427440", "LM_SD_1216_12"),
        ("GSM8427441", "LM_SD_13"),
        ("GSM8427442", "LM_SD_1216_14"),
        ("GSM8427443", "LM_SD_15"),
    ]:
        prefix = f"{gsm}_{sid}"
        if (raw273 / f"{prefix}_matrix.mtx.gz").exists():
            add(
                f"GSE273378_{sid}",
                "GSE273378",
                "LUAD_stageI",
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273378",
                lambda prefix=prefix, gsm=gsm, sid=sid: load_mtx_section(
                    raw273, prefix, sf_name=f"{gsm}_{sid}_scalefactors_json.json.gz"
                ),
            )

    # GSE300676
    raw300 = DATA / "geo/GSE300676"
    for gsm, sid in [
        ("GSM9066288", "CRC_mPAP3_A"),
        ("GSM9066289", "CRC_mPAP3_B"),
        ("GSM9066290", "CRC_mPAP4_A"),
        ("GSM9066291", "CRC_mPAP4_B"),
        ("GSM9066292", "CRC_mPAP1_A"),
        ("GSM9066293", "CRC_mPAP1_B"),
        ("GSM9066294", "CRC_mPAP2_A"),
        ("GSM9066295", "CRC_mPAP2_B"),
    ]:
        h5 = raw300 / f"{gsm}_{sid}_filtered_feature_bc_matrix.h5"
        sdir = raw300 / f"{sid}_spatial"
        if h5.exists() and sdir.exists():
            add(
                f"GSE300676_{sid}",
                "GSE300676",
                "LUAD_micropapillary",
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300676",
                lambda h5=h5, sdir=sdir: load_10x_section(h5, sdir),
            )

    # GSE307534 (invasive + precursor if present)
    root307 = DATA / "geo/GSE307534"
    man = ROOT / "scripts" / "sample_manifest_gse307534.tsv"
    if root307.exists() and man.exists():
        man_df = pd.read_csv(man, sep="\t")
        for _, row in man_df.iterrows():
            if row["class"] == "normal":
                continue
            d = root307 / f"{row['gsm']}_{row['label']}"
            if not (d / "filtered_feature_bc_matrix" / "matrix.mtx.gz").exists():
                continue
            add(
                f"GSE307534_{row['label']}",
                "GSE307534",
                f"{row['histology']}_{row['class']}",
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307534",
                lambda d=d: load_gse307_section(d),
            )

    if not rows:
        print("no sections found; download first", flush=True)
        return
    write_results(rows)
    print("wrote RESULTS.md n_kept=", sum(1 for r in rows if r.get("contrast_ok") == 1), flush=True)


if __name__ == "__main__":
    run()
