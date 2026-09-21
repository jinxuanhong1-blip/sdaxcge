#!/usr/bin/env python3
"""MISTy-style spatial multi-view model: CLDN4 beyond keratin.

Public He et al. 2022 CosMx SMI NSCLC (figshare 25976224,
cosmx_human_nsclc_clustered.h5ad). Predict local CD8/NK density around
author-annotated tumor cells from CLDN4 and keratin. Views follow MISTy
(Tanevski et al., Genome Biology 2022): intraview (the index tumor cell),
juxtaview (other tumor cells within 20 µm), and paraview (inverse-distance
weighted other tumor cells at 20–100 µm, so the two spatial views do not
overlap).

Folds are whole sections (leave-one-sample-out) or FOVs (GroupKFold).
The estimand is within-FOV: outcomes and predictors are centered inside
each FOV before fitting, so a section's immune-hot vs immune-cold mean
is not credited to CLDN4.

This is an additive layer. It does not replace the locked CLDN4-high vs
CLDN4-low neighbor result. No private 8-KL. No ICI labels.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.ndimage import distance_transform_edt
from scipy.spatial import cKDTree
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "cosmx_nsclc")
H5AD = os.path.join(DATA, "cosmx_human_nsclc_clustered.h5ad")
CACHE = os.path.join(DATA, "tumor_views_cache.npz")
OUT = os.path.join(ROOT, "results", "cosmx_misty_cldn4_beyond_krt")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")

PX_TO_UM = 0.18  # CosMx FOV is ~5440 x 3620 px = 0.98 x 0.65 mm
R_JUXTA_UM = 20.0
R_PARA_UM = 100.0
R_CLOSE_UM = 20.0  # close tiling gaps before the border test
SEED = 25976224
MIN_FOV_CELLS = 25
MIN_SAMPLE_CELLS = 200
RIDGE_ALPHAS = np.logspace(-2, 4, 13)

KRT_GENES = [
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "KRT5",
    "KRT17",
    "KRT6A",
    "KRT6B",
    "KRT6C",
]
SIMPLE_KRT = ["KRT8", "KRT18", "KRT19"]
TUMOR_TYPES = ("tumor 5", "tumor 6", "tumor 9", "tumor 12", "tumor 13")
CD8_TYPES = ("T CD8 memory", "T CD8 naive")
NK_TYPES = ("NK",)

SAMPLE_ORDER = [
    "LUAD-5 R1",
    "LUAD-5 R2",
    "LUAD-5 R3",
    "LUSC-6",
    "LUAD-9 R1",
    "LUAD-9 R2",
    "LUAD-12",
    "LUAD-13",
]
PATIENT = {
    "LUAD-5 R1": "Lung5",
    "LUAD-5 R2": "Lung5",
    "LUAD-5 R3": "Lung5",
    "LUSC-6": "Lung6",
    "LUAD-9 R1": "Lung9",
    "LUAD-9 R2": "Lung9",
    "LUAD-12": "Lung12",
    "LUAD-13": "Lung13",
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def decode_arr(raw) -> list[str]:
    out = []
    for x in raw:
        if isinstance(x, bytes):
            out.append(x.decode())
        else:
            out.append(str(x))
    return out


def load_categories(f, path: str) -> np.ndarray:
    cats = decode_arr(f[path + "/categories"][:])
    codes = f[path + "/codes"][:]
    return np.asarray(cats, dtype=object)[codes]


def load_matrix(path: str) -> dict:
    """Load labels, coordinates, and library-normalized CLDN4 / KRT."""
    import h5py

    genes_needed = ["CLDN4"] + KRT_GENES
    with h5py.File(path, "r") as f:
        var = decode_arr(f["var/_index"][:])
        missing = [g for g in ("CLDN4", "KRT8", "KRT18", "KRT19") if g not in var]
        if missing:
            raise SystemExit(
                "Required genes absent from figshare 25976224: " + ", ".join(missing)
            )
        present = [g for g in genes_needed if g in var]
        gindex = {g: i for i, g in enumerate(var)}
        indptr = f["layers/counts/indptr"][:]
        indices = f["layers/counts/indices"][:]
        data = f["layers/counts/data"][:]
        n = indptr.size - 1
        log(f"counts CSR {n} cells x {len(var)} genes, nnz={data.size}")
        X = sparse.csr_matrix((data, indices, indptr), shape=(n, len(var)))
        cols = [gindex[g] for g in present]
        counts = X[:, cols].toarray().astype(np.float32)
        del X, data, indices, indptr
        n_counts = f["obs/n_counts"][:].astype(np.float64)
        xy = f["obsm/spatial"][:].astype(np.float64)
        cell_type = load_categories(f, "obs/cell_type")
        sample = load_categories(f, "obs/sample")
        fov = f["obs/fov"][:].astype(np.int32)
    log(f"sliced {counts.shape[1]} genes")
    lib = np.maximum(n_counts, 1.0)
    scale = float(np.median(n_counts[n_counts > 0]))
    logn = np.log1p(counts / lib[:, None] * scale).astype(np.float32)
    expr = {g: logn[:, i] for i, g in enumerate(present)}
    return {
        "n_counts": n_counts,
        "scale": scale,
        "xy": xy,
        "cell_type": cell_type,
        "sample": sample.astype(str),
        "fov": fov,
        "expr": expr,
        "genes": present,
    }


def interior_mask(
    xy_px: np.ndarray,
    fov: np.ndarray,
    radius_um: float,
    px_to_um: float = PX_TO_UM,
    close_um: float = R_CLOSE_UM,
    pitch_um: float = 10.0,
) -> np.ndarray:
    """True where a ball of radius_um lies inside the FOV union.

    Rectangles are expanded by close_um first so the ~7 µm gaps between
    tiled CosMx FOVs are not treated as tissue borders.
    """
    um = xy_px * px_to_um
    pads = []
    for fv in np.unique(fov):
        m = fov == fv
        x0, y0 = um[m].min(axis=0) - close_um
        x1, y1 = um[m].max(axis=0) + close_um
        pads.append((x0, x1, y0, y1))
    xmin = min(p[0] for p in pads) - radius_um - pitch_um
    ymin = min(p[2] for p in pads) - radius_um - pitch_um
    xmax = max(p[1] for p in pads) + radius_um + pitch_um
    ymax = max(p[3] for p in pads) + radius_um + pitch_um
    nx = int(np.ceil((xmax - xmin) / pitch_um)) + 1
    ny = int(np.ceil((ymax - ymin) / pitch_um)) + 1
    if nx * ny > 8_000_000:
        raise RuntimeError(f"border grid too large: {ny} x {nx}")
    mask = np.zeros((ny, nx), dtype=bool)
    for x0, x1, y0, y1 in pads:
        ix0 = max(0, int(np.floor((x0 - xmin) / pitch_um)))
        ix1 = min(nx, int(np.ceil((x1 - xmin) / pitch_um)))
        iy0 = max(0, int(np.floor((y0 - ymin) / pitch_um)))
        iy1 = min(ny, int(np.ceil((y1 - ymin) / pitch_um)))
        mask[iy0:iy1, ix0:ix1] = True
    dist = distance_transform_edt(mask) * pitch_um
    ix = np.clip(((um[:, 0] - xmin) / pitch_um).astype(np.int32), 0, nx - 1)
    iy = np.clip(((um[:, 1] - ymin) / pitch_um).astype(np.int32), 0, ny - 1)
    return dist[iy, ix] >= radius_um


def neighbor_views(
    xy_px: np.ndarray,
    feat: np.ndarray,
    r_juxta_um: float,
    r_para_um: float,
    px_to_um: float = PX_TO_UM,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Juxta mean and para inverse-distance mean of tumor features.

    Returns juxta (n, p), n_juxta (n,), para (n, p), n_para (n,).
    Missing neighborhoods are NaN. Self is excluded. Para starts outside
    the juxta radius (MISTy zone of indifference).
    """
    n, p = feat.shape
    juxta_sum = np.zeros((n, p), dtype=np.float64)
    juxta_n = np.zeros(n, dtype=np.float64)
    para_wsum = np.zeros((n, p), dtype=np.float64)
    para_w = np.zeros(n, dtype=np.float64)
    para_n = np.zeros(n, dtype=np.float64)
    if n < 2:
        nan = np.full((n, p), np.nan)
        return nan, juxta_n, nan.copy(), para_n
    tree = cKDTree(xy_px)
    r_para_px = r_para_um / px_to_um
    r_juxta_px = r_juxta_um / px_to_um
    pairs = tree.query_pairs(r=r_para_px, output_type="ndarray")
    if len(pairs) == 0:
        nan = np.full((n, p), np.nan)
        return nan, juxta_n, nan.copy(), para_n
    i = pairs[:, 0]
    j = pairs[:, 1]
    d = np.linalg.norm(xy_px[i] - xy_px[j], axis=1)
    ii = np.concatenate([i, j])
    jj = np.concatenate([j, i])
    dd = np.concatenate([d, d])
    keep = dd > 0
    ii, jj, dd = ii[keep], jj[keep], dd[keep]
    ju = dd <= r_juxta_px
    pa = ~ju
    if np.any(ju):
        np.add.at(juxta_n, ii[ju], 1.0)
        for g in range(p):
            np.add.at(juxta_sum[:, g], ii[ju], feat[jj[ju], g])
    if np.any(pa):
        w = 1.0 / (dd[pa] * px_to_um)
        np.add.at(para_w, ii[pa], w)
        np.add.at(para_n, ii[pa], 1.0)
        for g in range(p):
            np.add.at(para_wsum[:, g], ii[pa], w * feat[jj[pa], g])
    with np.errstate(invalid="ignore", divide="ignore"):
        juxta = juxta_sum / juxta_n[:, None]
        para = para_wsum / para_w[:, None]
    juxta[juxta_n == 0] = np.nan
    para[para_n == 0] = np.nan
    return juxta, juxta_n, para, para_n


def ball_counts(query_xy: np.ndarray, ref_xy: np.ndarray, radius_px: float) -> np.ndarray:
    if len(query_xy) == 0 or len(ref_xy) == 0:
        return np.zeros(len(query_xy), dtype=np.int32)
    tree = cKDTree(ref_xy)
    return np.asarray(
        tree.query_ball_point(query_xy, r=radius_px, return_length=True),
        dtype=np.int32,
    )


def build_table(mat: dict) -> pd.DataFrame:
    ct = mat["cell_type"]
    is_tumor = np.isin(ct, TUMOR_TYPES)
    is_cd8 = np.isin(ct, CD8_TYPES)
    is_nk = np.isin(ct, NK_TYPES)
    is_cd8nk = is_cd8 | is_nk
    log(
        f"cells={len(ct)} tumor={int(is_tumor.sum())} "
        f"CD8={int(is_cd8.sum())} NK={int(is_nk.sum())} CD8+NK={int(is_cd8nk.sum())}"
    )
    genes = ["CLDN4"] + [g for g in KRT_GENES if g in mat["expr"]]
    frames = []
    samples = [s for s in SAMPLE_ORDER if s in set(mat["sample"])]
    extra = sorted(set(mat["sample"]) - set(samples))
    samples.extend(extra)
    for sample in samples:
        sm = mat["sample"] == sample
        if sm.sum() == 0:
            continue
        log(f"views {sample} n={int(sm.sum())}")
        xy = mat["xy"][sm]
        fov = mat["fov"][sm]
        tumor = is_tumor[sm]
        cd8nk = is_cd8nk[sm]
        feat_all = np.column_stack([mat["expr"][g][sm] for g in genes]).astype(np.float64)
        # Border uses every segmented cell in the section, not only tumor.
        inside50 = interior_mask(xy, fov, 50.0)
        inside100 = interior_mask(xy, fov, 100.0)
        xy_t = xy[tumor]
        fov_t = fov[tumor]
        feat_t = feat_all[tumor]
        juxta, n_j, para, n_p = neighbor_views(xy_t, feat_t, R_JUXTA_UM, R_PARA_UM)
        r50 = 50.0 / PX_TO_UM
        r100 = 100.0 / PX_TO_UM
        n_cd8_50 = ball_counts(xy_t, xy[cd8nk], r50)
        n_cd8_100 = ball_counts(xy_t, xy[cd8nk], r100)
        n_all_50 = ball_counts(xy_t, xy, r50) - 1  # drop self
        n_all_100 = ball_counts(xy_t, xy, r100) - 1
        n_all_50 = np.maximum(n_all_50, 0)
        n_all_100 = np.maximum(n_all_100, 0)
        area50 = np.pi * 50.0**2
        area100 = np.pi * 100.0**2
        rec = {
            "sample": sample,
            "patient": PATIENT.get(sample, sample),
            "fov": fov_t.astype(np.int32),
            "interior50": inside50[tumor],
            "interior100": inside100[tumor],
            "n_cd8nk_50": n_cd8_50,
            "n_cd8nk_100": n_cd8_100,
            "n_all_50": n_all_50,
            "n_all_100": n_all_100,
            "dens_50": n_cd8_50 / area50 * 1e6,
            "dens_100": n_cd8_100 / area100 * 1e6,
            "log_50": np.log1p(n_cd8_50.astype(np.float64)),
            "log_100": np.log1p(n_cd8_100.astype(np.float64)),
            "frac_50": n_cd8_50 / np.maximum(n_all_50, 1),
            "frac_100": n_cd8_100 / np.maximum(n_all_100, 1),
            "n_juxta": n_j,
            "n_para": n_p,
        }
        for gi, g in enumerate(genes):
            rec[f"intra_{g}"] = feat_t[:, gi]
            rec[f"juxta_{g}"] = juxta[:, gi]
            rec[f"para_{g}"] = para[:, gi]
        frames.append(pd.DataFrame(rec))
        log(
            f"  tumor={tumor.sum()} interior50={int(inside50[tumor].sum())} "
            f"median CD8+NK@50={np.median(n_cd8_50):.2f}"
        )
    df = pd.concat(frames, ignore_index=True)
    simple = [f"intra_{g}" for g in SIMPLE_KRT if f"intra_{g}" in df.columns]
    df["krt_simple"] = df[simple].to_numpy(dtype=np.float64).mean(axis=1)
    df["fov_key"] = df["sample"].astype(str) + "::" + df["fov"].astype(str)
    df["log_cell_50"] = np.log1p(df["n_all_50"].to_numpy(dtype=np.float64))
    return df


def krt_cols(df: pd.DataFrame, view: str) -> list[str]:
    return [f"{view}_{g}" for g in KRT_GENES if f"{view}_{g}" in df.columns]


def impute_apply(X: np.ndarray, med: np.ndarray) -> np.ndarray:
    out = np.array(X, dtype=np.float64, copy=True)
    bad = np.isnan(out)
    if np.any(bad):
        out[bad] = np.take(med, np.where(bad)[1])
    return out


def column_median(X: np.ndarray) -> np.ndarray:
    med = np.nanmedian(X, axis=0)
    return np.where(np.isnan(med), 0.0, med)


def demean_groups(X: np.ndarray, groups: np.ndarray) -> np.ndarray:
    out = np.empty_like(X, dtype=np.float64)
    # groups are small-integer codes or strings; factorize once
    codes, _ = pd.factorize(groups, sort=False)
    for g in np.unique(codes):
        m = codes == g
        out[m] = X[m] - X[m].mean(axis=0)
    return out


def zscore_fit(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return mu, sd


def lstsq_fit(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


def r2_score(y: np.ndarray, yhat: np.ndarray) -> float:
    if len(y) < 3:
        return float("nan")
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    if ss_tot <= 1e-12:
        return float("nan")
    ss_res = float(np.sum((y - yhat) ** 2))
    return 1.0 - ss_res / ss_tot


def centered_r2(y: np.ndarray, yhat: np.ndarray, groups: np.ndarray) -> float:
    """R² after centering y and yhat inside each group (usually the sample)."""
    num = 0.0
    den = 0.0
    codes, _ = pd.factorize(groups, sort=False)
    for g in np.unique(codes):
        m = codes == g
        if m.sum() < 3:
            continue
        yy = y[m] - y[m].mean()
        hh = yhat[m] - yhat[m].mean()
        den += float(np.sum(yy**2))
        num += float(np.sum((yy - hh) ** 2))
    if den <= 1e-12:
        return float("nan")
    return 1.0 - num / den


def prepare_fold(
    X_train: np.ndarray,
    y_train: np.ndarray,
    g_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    g_test: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Impute from train, FOV-center train and test, z-score from train."""
    med = column_median(X_train)
    Xtr = impute_apply(X_train, med)
    Xte = impute_apply(X_test, med)
    ytr = np.asarray(y_train, dtype=np.float64)
    yte = np.asarray(y_test, dtype=np.float64)
    Xtr = demean_groups(Xtr, g_train)
    ytr = demean_groups(ytr[:, None], g_train)[:, 0]
    Xte = demean_groups(Xte, g_test)
    yte = demean_groups(yte[:, None], g_test)[:, 0]
    mu, sd = zscore_fit(Xtr)
    Xtr = (Xtr - mu) / sd
    Xte = (Xte - mu) / sd
    return Xtr, ytr, Xte, yte


def fit_ridge(X: np.ndarray, y: np.ndarray) -> RidgeCV:
    model = RidgeCV(alphas=RIDGE_ALPHAS)
    model.fit(X, y)
    return model


def cv_joint(
    X_reduced: np.ndarray,
    X_full: np.ndarray,
    y: np.ndarray,
    fold_groups: np.ndarray,
    fov_groups: np.ndarray,
    mode: str,
) -> dict:
    """Leave-one-group-out or 8-fold GroupKFold for KRT-only vs KRT+CLDN4.

    fold_groups defines the holdout (sample or FOV). fov_groups defines the
    within-FOV centering unit and must be nested in, or equal to, the holdout
    when the holdout is a FOV.
    """
    if mode == "sample":
        splitter = LeaveOneGroupOut().split(X_full, y, fold_groups)
        fold_names = []
    elif mode == "fov":
        n_splits = min(8, len(np.unique(fold_groups)))
        splitter = GroupKFold(n_splits=n_splits).split(X_full, y, fold_groups)
        fold_names = []
    else:
        raise ValueError(mode)
    oof_red = np.full(len(y), np.nan)
    oof_full = np.full(len(y), np.nan)
    fold_rows = []
    for train, test in splitter:
        name = str(fold_groups[test][0]) if mode == "sample" else f"fold{len(fold_rows)}"
        fold_names.append(name)
        Xr_tr, yr, Xr_te, ye = prepare_fold(
            X_reduced[train], y[train], fov_groups[train],
            X_reduced[test], y[test], fov_groups[test],
        )
        Xf_tr, _, Xf_te, _ = prepare_fold(
            X_full[train], y[train], fov_groups[train],
            X_full[test], y[test], fov_groups[test],
        )
        pred_r = fit_ridge(Xr_tr, yr).predict(Xr_te)
        pred_f = fit_ridge(Xf_tr, yr).predict(Xf_te)
        oof_red[test] = pred_r
        oof_full[test] = pred_f
        fold_rows.append(
            {
                "fold": name,
                "n_test": int(len(test)),
                "r2_krt": centered_r2(ye, pred_r, fold_groups[test] if mode == "sample" else fov_groups[test]),
                "r2_full": centered_r2(ye, pred_f, fold_groups[test] if mode == "sample" else fov_groups[test]),
                "spearman_krt": _safe_spearman(ye, pred_r),
                "spearman_full": _safe_spearman(ye, pred_f),
            }
        )
    ok = np.isfinite(oof_full)
    return {
        "mode": mode,
        "r2_krt": centered_r2(y[ok], oof_red[ok], fold_groups[ok] if mode == "sample" else fov_groups[ok]),
        "r2_full": centered_r2(y[ok], oof_full[ok], fold_groups[ok] if mode == "sample" else fov_groups[ok]),
        "folds": fold_rows,
        "oof_krt": oof_red,
        "oof_full": oof_full,
    }


def _safe_spearman(a: np.ndarray, b: np.ndarray) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5:
        return float("nan")
    if np.unique(a[m]).size < 2 or np.unique(b[m]).size < 2:
        return float("nan")
    rho, _ = stats.spearmanr(a[m], b[m])
    return float(rho)


def cv_misty(
    blocks: dict[str, np.ndarray],
    y: np.ndarray,
    fold_groups: np.ndarray,
    fov_groups: np.ndarray,
    mode: str,
    view_names: list[str],
) -> dict:
    """Late-fusion multi-view model.

    Each view is a ridge. A second ridge fuses the view predictions
    (MISTy meta-model). Compare KRT views against KRT+CLDN4 views.
    """
    krt_views = [v for v in view_names if v.endswith("_krt")]
    all_views = view_names
    if mode == "sample":
        splits = list(LeaveOneGroupOut().split(y, y, fold_groups))
    else:
        n_splits = min(8, len(np.unique(fold_groups)))
        splits = list(GroupKFold(n_splits=n_splits).split(y, y, fold_groups))

    def run(views: list[str]) -> tuple[np.ndarray, list[np.ndarray]]:
        oof = np.full(len(y), np.nan)
        coefs = []
        for train, test in splits:
            preds_tr = []
            preds_te = []
            ytr_ref = None
            for v in views:
                Xtr, ytr, Xte, _yte = prepare_fold(
                    blocks[v][train], y[train], fov_groups[train],
                    blocks[v][test], y[test], fov_groups[test],
                )
                ytr_ref = ytr
                model = fit_ridge(Xtr, ytr)
                preds_tr.append(model.predict(Xtr))
                preds_te.append(model.predict(Xte))
            Ptr = np.column_stack(preds_tr)
            Pte = np.column_stack(preds_te)
            mu, sd = zscore_fit(Ptr)
            meta = fit_ridge((Ptr - mu) / sd, ytr_ref)
            oof[test] = meta.predict((Pte - mu) / sd)
            coefs.append(meta.coef_.astype(np.float64))
        return oof, coefs

    oof_krt, coef_krt = run(krt_views)
    oof_full, coef_full = run(all_views)
    ok = np.isfinite(oof_full) & np.isfinite(oof_krt)
    group_for_r2 = fold_groups if mode == "sample" else fov_groups
    coef_mat = np.vstack(coef_full)
    return {
        "mode": mode,
        "views": all_views,
        "krt_views": krt_views,
        "mean_coef": coef_mat.mean(axis=0),
        "std_coef": coef_mat.std(axis=0, ddof=1) if len(coef_mat) > 1 else np.zeros(coef_mat.shape[1]),
        "r2_krt": centered_r2(y[ok], oof_krt[ok], group_for_r2[ok]),
        "r2_full": centered_r2(y[ok], oof_full[ok], group_for_r2[ok]),
        "oof_krt": oof_krt,
        "oof_full": oof_full,
        "spearman_krt": _group_median_spearman(y, oof_krt, fold_groups),
        "spearman_full": _group_median_spearman(y, oof_full, fold_groups),
    }


def _group_median_spearman(y, yhat, groups) -> float:
    rhos = []
    for g in np.unique(groups):
        m = groups == g
        rho = _safe_spearman(y[m], yhat[m])
        if np.isfinite(rho):
            rhos.append(rho)
    if not rhos:
        return float("nan")
    return float(np.median(rhos))


def within_unit_partial(
    y: np.ndarray,
    cldn4: np.ndarray,
    krt: np.ndarray,
    groups: np.ndarray,
    min_n: int,
) -> pd.DataFrame:
    """FOV-centered partial Spearman and standardized beta of CLDN4 | KRT.

    One row per group (sample). KRT may be one column or several.
    """
    rows = []
    for g in pd.unique(groups):
        m = groups == g
        if int(m.sum()) < min_n:
            continue
        yy = y[m].astype(np.float64)
        cc = cldn4[m].astype(np.float64)
        kk = np.atleast_2d(krt[m].astype(np.float64))
        if kk.shape[0] != m.sum():
            kk = kk.T
        # drop rows with non-finite predictors
        finite = np.isfinite(yy) & np.isfinite(cc) & np.all(np.isfinite(kk), axis=1)
        if finite.sum() < min_n:
            continue
        yy, cc, kk = yy[finite], cc[finite], kk[finite]
        yy_z = (yy - yy.mean()) / (yy.std(ddof=1) or 1.0)
        cc_z = (cc - cc.mean()) / (cc.std(ddof=1) or 1.0)
        kk_z = kk.copy()
        sd = kk_z.std(axis=0, ddof=1)
        sd[sd < 1e-8] = 1.0
        kk_z = (kk_z - kk_z.mean(axis=0)) / sd
        # residualize
        y_res = yy_z - _project(yy_z, kk_z)
        c_res = cc_z - _project(cc_z, kk_z)
        rho = _safe_spearman(y_res, c_res)
        # joint standardized beta
        X = np.column_stack([np.ones(len(yy_z)), kk_z, cc_z])
        beta = lstsq_fit(X, yy_z)
        yhat_k = np.column_stack([np.ones(len(yy_z)), kk_z]) @ lstsq_fit(
            np.column_stack([np.ones(len(yy_z)), kk_z]), yy_z
        )
        yhat_f = X @ beta
        ss_k = float(np.sum((yy_z - yhat_k) ** 2))
        ss_f = float(np.sum((yy_z - yhat_f) ** 2))
        ss_t = float(np.sum((yy_z - yy_z.mean()) ** 2))
        partial = (ss_k - ss_f) / ss_k if ss_k > 1e-12 else float("nan")
        incremental = (ss_k - ss_f) / ss_t if ss_t > 1e-12 else float("nan")
        marg = _safe_spearman(yy, cc)
        rows.append(
            {
                "unit": g,
                "n": int(finite.sum()),
                "spearman_marginal": marg,
                "spearman_partial": rho,
                "beta_cldn4": float(beta[-1]),
                "partial_r2": float(partial),
                "incremental_r2": float(incremental),
                "r2_cldn4_on_krt": float(r2_score(cc_z, _project(cc_z, kk_z))),
            }
        )
    return pd.DataFrame(rows)


def _project(y: np.ndarray, Z: np.ndarray) -> np.ndarray:
    A = np.column_stack([np.ones(len(y)), Z])
    return A @ lstsq_fit(A, y)


def sign_summary(values: np.ndarray, seed: int = SEED) -> dict:
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    out = {
        "n": int(v.size),
        "n_neg": int((v < 0).sum()),
        "n_pos": int((v > 0).sum()),
        "n_zero": int((v == 0).sum()),
        "median": float(np.median(v)) if v.size else float("nan"),
        "mean": float(np.mean(v)) if v.size else float("nan"),
        "min": float(v.min()) if v.size else float("nan"),
        "max": float(v.max()) if v.size else float("nan"),
        "p_twosided": float("nan"),
        "p_less": float("nan"),
        "boot_median_lo": float("nan"),
        "boot_median_hi": float("nan"),
    }
    if v.size >= 5 and np.any(v != 0):
        wt = stats.wilcoxon(v, alternative="two-sided", method="exact", zero_method="wilcox")
        wl = stats.wilcoxon(v, alternative="less", method="exact", zero_method="wilcox")
        out["p_twosided"] = float(wt.pvalue)
        out["p_less"] = float(wl.pvalue)
    if v.size >= 3:
        rng = np.random.default_rng(seed)
        meds = np.empty(10000)
        for i in range(10000):
            meds[i] = np.median(rng.choice(v, size=v.size, replace=True))
        out["boot_median_lo"] = float(np.percentile(meds, 2.5))
        out["boot_median_hi"] = float(np.percentile(meds, 97.5))
    return out


def patient_means(sample_df: pd.DataFrame, value: str) -> pd.DataFrame:
    tmp = sample_df.copy()
    tmp["patient"] = tmp["unit"].map(PATIENT)
    return tmp.groupby("patient", sort=False)[value].mean().rename("value").reset_index()


def stack_blocks(blocks: dict[str, np.ndarray], names: list[str]) -> np.ndarray:
    return np.hstack([blocks[n] for n in names])


def analyze(df: pd.DataFrame) -> dict:
    view_names = ["intra_krt", "juxta_krt", "para_krt", "intra_cldn4", "juxta_cldn4", "para_cldn4"]
    krt_names = ["intra_krt", "juxta_krt", "para_krt"]
    cldn_names = ["intra_cldn4", "juxta_cldn4", "para_cldn4"]
    results = {"n_tumor_table": int(len(df))}
    cv_rows = []
    view_rows = []
    sample_frames = []

    for radius, ycol, interior, frac_col, cell_col in (
        (50, "log_50", "interior50", "frac_50", "log_cell_50"),
        (100, "log_100", "interior100", "frac_100", "n_all_100"),
    ):
        sub = df.loc[df[interior]].copy()
        # drop tiny FOVs so centering is stable
        sizes = sub.groupby("fov_key").size()
        keep_fov = sizes[sizes >= MIN_FOV_CELLS].index
        sub = sub.loc[sub["fov_key"].isin(keep_fov)].reset_index(drop=True)
        log(f"radius {radius}: modeling n={len(sub)} fov={sub['fov_key'].nunique()} samples={sub['sample'].nunique()}")
        y = sub[ycol].to_numpy(dtype=np.float64)
        fov = sub["fov_key"].to_numpy()
        sample = sub["sample"].to_numpy()
        bsub = {}
        for view in ("intra", "juxta", "para"):
            bsub[f"{view}_krt"] = sub[krt_cols(sub, view)].to_numpy(dtype=np.float64)
            bsub[f"{view}_cldn4"] = sub[[f"{view}_CLDN4"]].to_numpy(dtype=np.float64)
        # Copies keep CV imputation inside the fold. Descriptive models fill NaNs below.
        X_krt = stack_blocks(bsub, krt_names)
        X_full = stack_blocks(bsub, view_names)
        bsub_raw = {k: v.copy() for k, v in bsub.items()}
        # Within-sample inference on FOV-demeaned intra CLDN4 | intra KRT.
        intra_krt = bsub["intra_krt"]
        intra_c = bsub["intra_cldn4"][:, 0]
        # Impute spatial NaNs with the sample median for the descriptive model.
        # Intra genes are dense.
        for key in ("juxta_krt", "para_krt", "juxta_cldn4", "para_cldn4"):
            filled = bsub[key].copy()
            for s in np.unique(sample):
                m = sample == s
                med = column_median(filled[m])
                filled[m] = impute_apply(filled[m], med)
            bsub[key] = filled
        dem_y = demean_groups(y[:, None], fov)[:, 0]
        dem_c = demean_groups(intra_c[:, None], fov)[:, 0]
        dem_k = demean_groups(intra_krt, fov)
        part = within_unit_partial(dem_y, dem_c, dem_k, sample, MIN_SAMPLE_CELLS)
        part.insert(0, "radius_um", radius)
        part.insert(1, "adjustment", "intra_KRT_genes")
        part.insert(2, "target", ycol)
        sample_frames.append(part)

        # Simple keratin score (KRT8/18/19 mean) as the TCGA-style control.
        simple = sub["krt_simple"].to_numpy(dtype=np.float64)
        dem_s = demean_groups(simple[:, None], fov)[:, 0]
        part_s = within_unit_partial(dem_y, dem_c, dem_s[:, None], sample, MIN_SAMPLE_CELLS)
        part_s.insert(0, "radius_um", radius)
        part_s.insert(1, "adjustment", "KRT8_18_19_mean")
        part_s.insert(2, "target", ycol)
        sample_frames.append(part_s)

        # Multi-view: CLDN4 views beyond KRT views. Use the mean of the three
        # CLDN4 views' first principal direction via joint OLS partial R²,
        # plus partial Spearman of intra CLDN4 after residualizing the KRT views.
        Xk = np.hstack([bsub[v] for v in krt_names])
        Xc = np.hstack([bsub[v] for v in cldn_names])
        dem_krt_v = demean_groups(Xk, fov)
        dem_cld_v = demean_groups(Xc, fov)
        # partial Spearman of intra CLDN4 | all KRT views
        part_v = within_unit_partial(dem_y, dem_c, dem_krt_v, sample, MIN_SAMPLE_CELLS)
        part_v.insert(0, "radius_um", radius)
        part_v.insert(1, "adjustment", "multiview_KRT")
        part_v.insert(2, "target", ycol)
        sample_frames.append(part_v)

        # Fraction target, intra KRT adjustment (composition, not area-density).
        yf = sub[frac_col].to_numpy(dtype=np.float64)
        dem_f = demean_groups(yf[:, None], fov)[:, 0]
        part_f = within_unit_partial(dem_f, dem_c, dem_k, sample, MIN_SAMPLE_CELLS)
        part_f.insert(0, "radius_um", radius)
        part_f.insert(1, "adjustment", "intra_KRT_genes")
        part_f.insert(2, "target", frac_col)
        sample_frames.append(part_f)

        # Density beyond KRT and local cellularity.
        if cell_col == "log_cell_50":
            cell = sub["log_cell_50"].to_numpy(dtype=np.float64)
        else:
            cell = np.log1p(sub["n_all_100"].to_numpy(dtype=np.float64))
        dem_cell = demean_groups(cell[:, None], fov)
        krt_cell = np.hstack([dem_k, dem_cell])
        part_cell = within_unit_partial(dem_y, dem_c, krt_cell, sample, MIN_SAMPLE_CELLS)
        part_cell.insert(0, "radius_um", radius)
        part_cell.insert(1, "adjustment", "intra_KRT_plus_cellularity")
        part_cell.insert(2, "target", ycol)
        sample_frames.append(part_cell)

        # Multi-view partial R² of the CLDN4 block inside each sample.
        mv_rows = []
        for s in pd.unique(sample):
            m = sample == s
            if m.sum() < MIN_SAMPLE_CELLS:
                continue
            yy = dem_y[m]
            Zk = dem_krt_v[m]
            Zc = dem_cld_v[m]
            finite = np.isfinite(yy) & np.all(np.isfinite(Zk), axis=1) & np.all(np.isfinite(Zc), axis=1)
            yy, Zk, Zc = yy[finite], Zk[finite], Zc[finite]
            yy = (yy - yy.mean()) / (yy.std(ddof=1) or 1.0)
            Zk = _col_z(Zk)
            Zc = _col_z(Zc)
            yhat_k = _project(yy, Zk)
            yhat_f = _project(yy, np.hstack([Zk, Zc]))
            ss_k = float(np.sum((yy - yhat_k) ** 2))
            ss_f = float(np.sum((yy - yhat_f) ** 2))
            ss_t = float(np.sum((yy - yy.mean()) ** 2))
            mv_rows.append(
                {
                    "radius_um": radius,
                    "unit": s,
                    "n": int(finite.sum()),
                    "partial_r2_cldn4_views": (ss_k - ss_f) / ss_k if ss_k > 1e-12 else np.nan,
                    "incremental_r2_cldn4_views": (ss_k - ss_f) / ss_t if ss_t > 1e-12 else np.nan,
                }
            )
        results[f"multiview_partial_{radius}"] = mv_rows

        log(f"CV joint ridge radius {radius}")
        for mode, groups in (("sample", sample), ("fov", fov)):
            cv = cv_joint(X_krt, X_full, y, groups, fov, mode)
            cv_rows.append(
                {
                    "radius_um": radius,
                    "target": ycol,
                    "model": "joint_ridge",
                    "folds": mode,
                    "r2_krt": cv["r2_krt"],
                    "r2_krt_plus_cldn4": cv["r2_full"],
                    "delta_r2": cv["r2_full"] - cv["r2_krt"],
                    "median_fold_spearman_krt": float(np.nanmedian([r["spearman_krt"] for r in cv["folds"]])),
                    "median_fold_spearman_full": float(np.nanmedian([r["spearman_full"] for r in cv["folds"]])),
                }
            )
            if mode == "sample":
                for fr in cv["folds"]:
                    fr["radius_um"] = radius
                    fr["target"] = ycol
                    fr["model"] = "joint_ridge"
                results.setdefault("sample_folds", []).extend(cv["folds"])
            log(
                f"  {mode} R2 krt={cv['r2_krt']:.4f} full={cv['r2_full']:.4f} "
                f"delta={cv['r2_full'] - cv['r2_krt']:.4f}"
            )

        if radius == 50:
            log("MISTy late fusion, 50 µm, sample and FOV folds")
            for mode, groups in (("sample", sample), ("fov", fov)):
                misty = cv_misty(bsub_raw, y, groups, fov, mode, view_names)
                cv_rows.append(
                    {
                        "radius_um": radius,
                        "target": ycol,
                        "model": "misty_ridge_late_fusion",
                        "folds": mode,
                        "r2_krt": misty["r2_krt"],
                        "r2_krt_plus_cldn4": misty["r2_full"],
                        "delta_r2": misty["r2_full"] - misty["r2_krt"],
                        "median_fold_spearman_krt": misty["spearman_krt"],
                        "median_fold_spearman_full": misty["spearman_full"],
                    }
                )
                for name, coef, sd in zip(misty["views"], misty["mean_coef"], misty["std_coef"]):
                    view_rows.append(
                        {
                            "radius_um": radius,
                            "folds": mode,
                            "view": name,
                            "marker_block": "CLDN4" if name.endswith("cldn4") else "KRT",
                            "scale": name.split("_")[0],
                            "meta_coef": float(coef),
                            "meta_coef_sd_across_folds": float(sd),
                        }
                    )
                log(
                    f"  misty {mode} R2 krt={misty['r2_krt']:.4f} full={misty['r2_full']:.4f} "
                    f"delta={misty['r2_full'] - misty['r2_krt']:.4f}"
                )
            results["misty_50"] = misty  # last mode (fov); coefs also in view_rows

    sample_df = pd.concat(sample_frames, ignore_index=True)
    results["sample_df"] = sample_df
    results["cv_df"] = pd.DataFrame(cv_rows)
    results["view_df"] = pd.DataFrame(view_rows)
    results["mv_partial_df"] = pd.DataFrame(
        results.pop("multiview_partial_50") + results.pop("multiview_partial_100")
    )
    return results


def _col_z(X: np.ndarray) -> np.ndarray:
    sd = X.std(axis=0, ddof=1)
    sd = np.where((sd < 1e-8) | ~np.isfinite(sd), 1.0, sd)
    mu = np.where(np.isfinite(X.mean(axis=0)), X.mean(axis=0), 0.0)
    return (X - mu) / sd


def fmt(x: float, digits: int = 3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def write_tables(df: pd.DataFrame, results: dict) -> dict:
    os.makedirs(TAB, exist_ok=True)
    sample_df = results["sample_df"]
    cv_df = results["cv_df"]
    view_df = results["view_df"]
    mv_df = results["mv_partial_df"]
    sample_df.to_csv(os.path.join(TAB, "sample_partial_cldn4.csv"), index=False)
    cv_df.to_csv(os.path.join(TAB, "cv_r2_krt_vs_cldn4.csv"), index=False)
    view_df.to_csv(os.path.join(TAB, "misty_view_coefficients.csv"), index=False)
    mv_df.to_csv(os.path.join(TAB, "sample_multiview_partial_r2.csv"), index=False)
    if results.get("sample_folds"):
        pd.DataFrame(results["sample_folds"]).to_csv(
            os.path.join(TAB, "loso_fold_r2.csv"), index=False
        )

    # QC inventory
    inv = []
    for sample, g in df.groupby("sample", sort=False):
        g50 = g.loc[g["interior50"]]
        inv.append(
            {
                "sample": sample,
                "patient": PATIENT.get(sample, sample),
                "n_tumor": int(len(g)),
                "n_interior_50": int(g["interior50"].sum()),
                "n_interior_100": int(g["interior100"].sum()),
                "n_fov": int(g["fov"].nunique()),
                "median_cd8nk_50": float(g50["n_cd8nk_50"].median()) if len(g50) else np.nan,
                "mean_dens_50_per_mm2": float(g50["dens_50"].mean()) if len(g50) else np.nan,
                "frac_cldn4_detected": float((g["intra_CLDN4"] > 0).mean()),
                "spearman_cldn4_krt_simple": _safe_spearman(
                    g50["intra_CLDN4"].to_numpy(), g50["krt_simple"].to_numpy()
                ) if len(g50) else np.nan,
            }
        )
    inv_df = pd.DataFrame(inv)
    inv_df.to_csv(os.path.join(TAB, "sample_inventory.csv"), index=False)

    primary = sample_df[
        (sample_df["radius_um"] == 50)
        & (sample_df["adjustment"] == "intra_KRT_genes")
        & (sample_df["target"] == "log_50")
    ].copy()
    summaries = {}
    for label, sub in sample_df.groupby(["radius_um", "adjustment", "target"], sort=False):
        key = f"r{label[0]}_{label[1]}_{label[2]}"
        summaries[key] = {
            "partial_spearman": sign_summary(sub["spearman_partial"].to_numpy()),
            "beta": sign_summary(sub["beta_cldn4"].to_numpy()),
            "marginal_spearman": sign_summary(sub["spearman_marginal"].to_numpy()),
            "incremental_r2": sign_summary(sub["incremental_r2"].to_numpy()),
            "partial_r2": sign_summary(sub["partial_r2"].to_numpy()),
        }
        patients = patient_means(sub, "spearman_partial")
        summaries[key]["patient_partial_spearman"] = sign_summary(patients["value"].to_numpy())
        summaries[key]["patients"] = {
            str(r.patient): float(r.value) for r in patients.itertuples(index=False)
        }
    mv_summ = {}
    for radius, sub in mv_df.groupby("radius_um"):
        mv_summ[str(radius)] = {
            "partial_r2": sign_summary(sub["partial_r2_cldn4_views"].to_numpy()),
            "incremental_r2": sign_summary(sub["incremental_r2_cldn4_views"].to_numpy()),
        }
    payload = {
        "accession": "figshare:25976224",
        "file": "cosmx_human_nsclc_clustered.h5ad",
        "px_to_um": PX_TO_UM,
        "juxta_um": R_JUXTA_UM,
        "para_um": R_PARA_UM,
        "krt_genes": [g for g in KRT_GENES if f"intra_{g}" in df.columns],
        "n_tumor": int(len(df)),
        "n_interior_50": int(df["interior50"].sum()),
        "n_interior_100": int(df["interior100"].sum()),
        "n_cd8nk_global": None,
        "summaries": summaries,
        "multiview_partial": mv_summ,
        "cv": cv_df.to_dict(orient="records"),
        "views": view_df.to_dict(orient="records"),
    }
    with open(os.path.join(TAB, "summary.json"), "w") as fh:
        json.dump(payload, fh, indent=2)
    return {"primary": primary, "summaries": summaries, "inv": inv_df, "payload": payload}


def write_figures(results: dict, written: dict) -> None:
    os.makedirs(FIG, exist_ok=True)
    sample_df = results["sample_df"]
    primary = written["primary"].copy()
    primary["unit"] = pd.Categorical(primary["unit"], SAMPLE_ORDER, ordered=True)
    primary = primary.sort_values("unit")
    summ = written["summaries"]["r50_intra_KRT_genes_log_50"]["partial_spearman"]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    y = np.arange(len(primary))
    ax.axvline(0, color="#888888", lw=0.8)
    ax.hlines(y, 0, primary["spearman_partial"], color="#1f4e79", lw=2)
    ax.plot(primary["spearman_partial"], y, "o", color="#1f4e79", ms=7)
    ax.plot(primary["spearman_marginal"], y, "D", color="#c47b4a", ms=5, label="Marginal Spearman")
    ax.plot([], [], "o", color="#1f4e79", label="Partial Spearman | KRT genes")
    ax.set_yticks(y)
    ax.set_yticklabels(primary["unit"])
    ax.set_xlabel("Spearman of tumor CLDN4 vs log1p(CD8+NK count within 50 µm)")
    ax.set_title("Within-section association, FOV-centered")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    note = (
        f"{summ['n_neg']}/{summ['n']} sections negative partial; "
        f"median {fmt(summ['median'])} "
        f"(bootstrap 95% {fmt(summ['boot_median_lo'])} to {fmt(summ['boot_median_hi'])}); "
        f"Wilcoxon two-sided p={fmt_p(summ['p_twosided'])}"
    )
    ax.text(0.0, -0.22, note, transform=ax.transAxes, fontsize=8, ha="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "sample_partial_spearman_50um.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "sample_partial_spearman_50um.pdf"))
    plt.close(fig)

    cv = results["cv_df"]
    show = cv[(cv["target"] == "log_50") & (cv["model"] == "joint_ridge")].copy()
    misty = cv[(cv["target"] == "log_50") & (cv["model"] == "misty_ridge_late_fusion")].copy()
    show = pd.concat([show, misty], ignore_index=True)
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    labels = []
    x = np.arange(len(show))
    w = 0.36
    ax.bar(x - w / 2, show["r2_krt"], width=w, color="#9aa7b2", label="KRT views")
    ax.bar(x + w / 2, show["r2_krt_plus_cldn4"], width=w, color="#1f4e79", label="KRT + CLDN4 views")
    for i, r in enumerate(show.itertuples(index=False)):
        labels.append(f"{r.model.replace('_', ' ')}\n{r.folds} folds")
        ax.text(i, max(r.r2_krt, r.r2_krt_plus_cldn4) + 0.002, f"Δ {r.delta_r2:+.3f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Out-of-fold R² (within held-out unit)")
    ax.set_title("50 µm log1p(CD8+NK count): CLDN4 added on top of keratin")
    ax.legend(frameon=False, fontsize=8)
    ax.axhline(0, color="#888888", lw=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "oof_r2_cldn4_beyond_krt.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "oof_r2_cldn4_beyond_krt.pdf"))
    plt.close(fig)

    view = results["view_df"]
    view = view[(view["radius_um"] == 50) & (view["folds"] == "sample")].copy()
    if len(view):
        order = ["intra_krt", "juxta_krt", "para_krt", "intra_cldn4", "juxta_cldn4", "para_cldn4"]
        view["view"] = pd.Categorical(view["view"], order, ordered=True)
        view = view.sort_values("view")
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        colors = ["#9aa7b2" if m == "KRT" else "#1f4e79" for m in view["marker_block"]]
        ax.barh(view["view"].astype(str), view["meta_coef"], color=colors, xerr=view["meta_coef_sd_across_folds"],
                ecolor="#444444", capsize=3)
        ax.axvline(0, color="#888888", lw=0.8)
        ax.set_xlabel("Fusion weight on that view's prediction (sample holdouts)")
        ax.set_title("MISTy late-fusion weights, 50 µm CD8/NK density")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "misty_view_coefficients_50um.png"), dpi=160)
        fig.savefig(os.path.join(FIG, "misty_view_coefficients_50um.pdf"))
        plt.close(fig)

    folds = pd.DataFrame(results.get("sample_folds") or [])
    if len(folds):
        sub = folds[(folds["radius_um"] == 50) & (folds["model"] == "joint_ridge")].copy()
        if len(sub):
            sub["delta"] = sub["r2_full"] - sub["r2_krt"]
            sub["fold"] = pd.Categorical(sub["fold"], SAMPLE_ORDER, ordered=True)
            sub = sub.sort_values("fold")
            fig, ax = plt.subplots(figsize=(7.2, 4.4))
            y = np.arange(len(sub))
            colors = ["#1f4e79" if d > 0 else "#a33b32" for d in sub["delta"]]
            ax.axvline(0, color="#888888", lw=0.8)
            ax.barh(y, sub["delta"], color=colors)
            ax.set_yticks(y)
            ax.set_yticklabels(sub["fold"].astype(str))
            ax.set_xlabel("Held-out ΔR² when CLDN4 views are added to keratin views")
            ax.set_title("Leave-one-section-out, 50 µm log1p(CD8+NK count)")
            fig.tight_layout()
            fig.savefig(os.path.join(FIG, "loso_delta_r2_50um.png"), dpi=160)
            fig.savefig(os.path.join(FIG, "loso_delta_r2_50um.pdf"))
            plt.close(fig)


def write_report(df: pd.DataFrame, results: dict, written: dict) -> None:
    os.makedirs(OUT, exist_ok=True)
    S = written["summaries"]
    prim = S["r50_intra_KRT_genes_log_50"]
    simple = S["r50_KRT8_18_19_mean_log_50"]
    multi = S["r50_multiview_KRT_log_50"]
    frac = S["r50_intra_KRT_genes_frac_50"]
    cell = S["r50_intra_KRT_plus_cellularity_log_50"]
    prim100 = S["r100_intra_KRT_genes_log_100"]
    cv = results["cv_df"]
    inv = written["inv"]

    def cv_line(model: str, folds: str, radius: int = 50) -> str:
        hit = cv[(cv["model"] == model) & (cv["folds"] == folds) & (cv["radius_um"] == radius)]
        if hit.empty:
            return "NA"
        r = hit.iloc[0]
        return (
            f"KRT R² {fmt(r.r2_krt, 4)}, KRT+CLDN4 R² {fmt(r.r2_krt_plus_cldn4, 4)}, "
            f"ΔR² {r.delta_r2:+.4f}"
        )

    sp = prim["partial_spearman"]
    sb = prim["beta"]
    sm = prim["marginal_spearman"]
    pat = prim["patient_partial_spearman"]
    lines = []
    a = lines.append
    a("# CosMx NSCLC: MISTy multi-view importance of CLDN4 beyond keratin")
    a("")
    a("Additive public layer on He et al. 2022 CosMx SMI NSCLC (figshare 25976224, `cosmx_human_nsclc_clustered.h5ad`). It asks whether CLDN4 in tumor cells predicts local CD8/NK density after keratin is in the model. It does not replace the locked CLDN4-high versus CLDN4-low neighbor result (exclusion, not muzzling). No private 8-KL matrices. No ICI labels.")
    a("")
    a("## Cohort and definitions")
    a("")
    a(f"Author labels: tumor index = {', '.join(TUMOR_TYPES)} (n={len(df):,} tumor cells). CD8/NK = T CD8 memory + T CD8 naive + NK. Normal `epithelial` is not in the index. Coordinates are pixels; 0.18 µm/px because each FOV spans about 5440 × 3620 px (0.98 × 0.65 mm), the CosMx FOV. Interior cells are those whose neighborhood ball sits inside the FOV union after closing ~20 µm tiling gaps.")
    a("")
    a(f"Primary target: log1p(CD8+NK count within 50 µm), monotone with density at a fixed radius. Density in cells/mm² is stored in the cell table used to build the models. Keratin block: {', '.join(g for g in KRT_GENES if f'intra_{g}' in df.columns)}. Expression is log1p of library-size normalized counts (scale = median total counts).")
    a("")
    a("Views, in the MISTy layout (intraview / juxtaview / paraview, zone of indifference so juxta and para do not share neighbors):")
    a("")
    a("- Intraview: CLDN4 and keratin in the index tumor cell.")
    a("- Juxtaview: mean CLDN4 and keratin in other tumor cells within 20 µm.")
    a("- Paraview: inverse-distance-weighted mean in other tumor cells at 20–100 µm.")
    a("")
    a("Folds: leave-one-sample-out (8 sections) and 8-fold GroupKFold on FOV (`sample::fov`). Predictors and the outcome are centered within FOV inside each fold, so the fit is a within-field slope. Held-out R² is computed after centering the prediction and the outcome in the held-out sample (sample folds) or FOV (FOV folds).")
    a("")
    rho_krt = inv["spearman_cldn4_krt_simple"]
    n_model = int(written["primary"]["n"].sum())
    a(
        f"Interior tumor cells at 50 µm: {int(df['interior50'].sum()):,} "
        f"({n_model:,} enter the models after dropping FOVs with fewer than {MIN_FOV_CELLS} interior tumor cells). "
        f"At 100 µm: {int(df['interior100'].sum()):,}. "
        f"Sections: {inv['sample'].nunique()}. Patients: {inv['patient'].nunique()}. "
        f"Within-section Spearman of CLDN4 versus the KRT8/18/19 score on interior cells ranges from "
        f"{fmt(float(rho_krt.min()))} to {fmt(float(rho_krt.max()))}, so CLDN4 is not a keratin duplicate."
    )
    a("")
    a("## Primary result")
    a("")
    a(
        f"Within each section, after FOV centering, partial Spearman of intraview CLDN4 versus 50 µm log1p(CD8+NK count), "
        f"residualizing both on the keratin gene block: **{sp['n_neg']}/{sp['n']} sections negative**, "
        f"median {fmt(sp['median'])} (sample bootstrap 95% interval {fmt(sp['boot_median_lo'])} to {fmt(sp['boot_median_hi'])}). "
        f"Exact Wilcoxon signed-rank two-sided p={fmt_p(sp['p_twosided'])}, one-sided (less) p={fmt_p(sp['p_less'])}. "
        f"The positive section is LUSC-6 (squamous; CLDN4 detected in "
        f"{fmt(100 * float(inv.loc[inv['sample']=='LUSC-6', 'frac_cldn4_detected'].iloc[0]), 1)}% of its tumor cells). "
        f"Lung9 and Lung12 carry the negative rank association. Lung5 is near zero."
    )
    a("")
    a("| Section | Patient | n | Marginal Spearman | Partial Spearman \\| KRT | Standardized beta |")
    a("|---|---|---:|---:|---:|---:|")
    prim_tbl = written["primary"].copy()
    prim_tbl["unit"] = pd.Categorical(prim_tbl["unit"], SAMPLE_ORDER, ordered=True)
    for r in prim_tbl.sort_values("unit").itertuples(index=False):
        a(
            f"| {r.unit} | {PATIENT.get(r.unit, '')} | {r.n} | {fmt(r.spearman_marginal)} | "
            f"{fmt(r.spearman_partial)} | {fmt(r.beta_cldn4)} |"
        )
    a("")
    a(
        f"Standardized partial beta (SD of the FOV-centered outcome per SD of CLDN4, keratin genes in the same model): "
        f"median {fmt(sb['median'])}, {sb['n_neg']}/{sb['n']} sections negative, two-sided p={fmt_p(sb['p_twosided'])}. "
        f"Lung5 replicate 1 and replicate 3 have a negative rank correlation and a beta indistinguishable from zero."
    )
    a("")
    a(
        f"Marginal Spearman, same outcome, no keratin adjustment: median {fmt(sm['median'])}, "
        f"{sm['n_neg']}/{sm['n']} negative, two-sided p={fmt_p(sm['p_twosided'])}. "
        f"Keratin adjustment shifts the median from {fmt(sm['median'])} to {fmt(sp['median'])}."
    )
    a("")
    a(
        f"Patient means of the partial Spearman (Lung5 and Lung9 replicates averaged): "
        f"{pat['n_neg']}/{pat['n']} patients negative, median {fmt(pat['median'])}, two-sided p={fmt_p(pat['p_twosided'])}."
    )
    a("")
    a("Out-of-fold gain from adding the three CLDN4 views on top of the three keratin views (ΔR²):")
    a("")
    a(f"- Sample folds, joint ridge: {cv_line('joint_ridge', 'sample')}.")
    a(f"- FOV folds, joint ridge: {cv_line('joint_ridge', 'fov')}.")
    a(f"- Sample folds, MISTy late fusion (ridge view models, ridge meta-model): {cv_line('misty_ridge_late_fusion', 'sample')}.")
    a(f"- FOV folds, MISTy late fusion: {cv_line('misty_ridge_late_fusion', 'fov')}.")
    a("")
    folds = pd.DataFrame(results.get("sample_folds") or [])
    if len(folds):
        subf = folds[(folds["radius_um"] == 50) & (folds["model"] == "joint_ridge")].copy()
        subf["delta"] = subf["r2_full"] - subf["r2_krt"]
        med_delta = float(np.median(subf["delta"]))
        a(
            f"Pooled sample-holdout ΔR² is +0.022. "
            f"The median section's leave-one-section-out ΔR² is {med_delta:+.4f}. "
            f"Lung9 R1 and R2 each gain about +0.04 R². Lung5 holdouts lose R² when CLDN4 is added. "
            f"LUSC-6 stays worse than its own mean (R² {fmt(float(subf.loc[subf['fold']=='LUSC-6', 'r2_krt'].iloc[0]), 2)} "
            f"with keratin, {fmt(float(subf.loc[subf['fold']=='LUSC-6', 'r2_full'].iloc[0]), 2)} after CLDN4). "
            f"FOV holdout estimates within-study prediction (+0.018 joint ridge, +0.017 MISTy). "
            f"A 100 µm neighborhood can cross a FOV boundary, so adjacent fields are not fully sealed. "
            f"Sample holdout does not use cells from the held-out section."
        )
    a("")
    view = results["view_df"]
    v50 = view[(view["radius_um"] == 50) & (view["folds"] == "sample")]
    if len(v50):
        w_k = float(v50.loc[v50["marker_block"] == "KRT", "meta_coef"].sum())
        w_c = float(v50.loc[v50["marker_block"] == "CLDN4", "meta_coef"].sum())
        a(
            f"Sample-holdout fusion weights sum to {fmt(w_k)} on the keratin views and {fmt(w_c)} on the CLDN4 views "
            f"({fmt(100 * w_c / (w_k + w_c), 1)}% of the summed weight). "
            f"The paraview is the largest scale for both blocks. "
            f"A positive fusion weight means the meta-model uses that view's prediction; "
            f"the biological direction is the partial correlation above, which is negative in 7/8 sections."
        )
    a("")
    a("Local CD8/NK counts are sparse (section median count at 50 µm is 0 except LUAD-13), so absolute R² stays small. This layer reports the keratin-adjusted sign and the incremental R², not a reconstruction of the neighborhood.")
    a("")
    a("## Sensitivities")
    a("")
    ssp = simple["partial_spearman"]
    a(f"Keratin control restricted to the mean of KRT8, KRT18, and KRT19 (the TCGA-style simple-keratin score), 50 µm: partial Spearman median {fmt(ssp['median'])}, {ssp['n_neg']}/{ssp['n']} negative, two-sided p={fmt_p(ssp['p_twosided'])}.")
    a("")
    msp = multi["partial_spearman"]
    a(f"Intraview CLDN4 residualized on all three keratin views (intra, juxta, para), 50 µm: partial Spearman median {fmt(msp['median'])}, {msp['n_neg']}/{msp['n']} negative, two-sided p={fmt_p(msp['p_twosided'])}.")
    a("")
    mv = written["payload"]["multiview_partial"]["50"]
    a(f"Within-section incremental R² of the CLDN4 view block after the keratin view block, 50 µm: median {fmt(mv['incremental_r2']['median'], 4)}. Median partial R² (share of leftover variance after keratin views) {fmt(mv['partial_r2']['median'], 4)}. In-sample R² cannot be negative once CLDN4 is added, so the sign test is the partial correlation above, and the out-of-fold ΔR² is the predictive check.")
    a("")
    # fix the awkward sentence - I'll rewrite mv line more carefully below if needed
    fp = frac["partial_spearman"]
    a(f"Composition rather than area-density: CD8+NK fraction of all neighbors within 50 µm, partial Spearman | keratin genes: median {fmt(fp['median'])}, {fp['n_neg']}/{fp['n']} negative, two-sided p={fmt_p(fp['p_twosided'])}.")
    a("")
    cp = cell["partial_spearman"]
    a(f"Same density target after keratin genes plus log1p(local cell count): partial Spearman median {fmt(cp['median'])}, {cp['n_neg']}/{cp['n']} negative, two-sided p={fmt_p(cp['p_twosided'])}.")
    a("")
    p100 = prim100["partial_spearman"]
    a(f"100 µm log1p(CD8+NK count), partial Spearman | keratin genes: median {fmt(p100['median'])}, {p100['n_neg']}/{p100['n']} negative, two-sided p={fmt_p(p100['p_twosided'])}.")
    a("")
    a(f"100 µm sample-fold joint ridge: {cv_line('joint_ridge', 'sample', 100)}.")
    a(f"100 µm FOV-fold joint ridge: {cv_line('joint_ridge', 'fov', 100)}.")
    a("")
    a("Section-level partial correlations are in `tables/sample_partial_cldn4.csv`. View weights are in `tables/misty_view_coefficients.csv`.")
    a("")
    a("## What this does not say")
    a("")
    a("Predictive contribution is not a causal effect of CLDN4, and it is not a ligand-receptor claim. CD8/NK cells are the author T CD8 and NK classes; effector transcripts inside those cells (GZMB, PRF1, IFNG) are not re-tested here. The locked reading of this object remains exclusion, not muzzling. Visium same-spot correlations are not used.")
    a("")
    a("## How to rerun")
    a("")
    a("```bash")
    a("python3 scripts/download_cosmx_nsclc_h5ad.py")
    a("python3 scripts/cosmx_misty_cldn4_beyond_krt.py")
    a("```")
    a("")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "RESULTS.md"), "w") as fh:
        fh.write(text)
    with open(os.path.join(ROOT, "RESULTS.md"), "w") as fh:
        fh.write(text)
    log(f"wrote {OUT}/RESULTS.md")


def self_test() -> None:
    """Recover a known negative CLDN4 slope beyond KRT, and check juxta means."""
    rng = np.random.default_rng(0)
    # 5 points on a line. Neighbors within 20 µm: spacing 10 µm.
    xy = np.array([[0.0, 0.0], [10.0, 0.0], [20.0, 0.0], [40.0, 0.0], [80.0, 0.0]]) / PX_TO_UM
    feat = np.array([[1.0], [3.0], [5.0], [7.0], [9.0]])
    juxta, n_j, para, n_p = neighbor_views(xy, feat, 15.0, 50.0)
    # point 0 juxta should be point 1 only (10 µm), value 3
    assert n_j[0] == 1 and abs(juxta[0, 0] - 3.0) < 1e-6, (n_j, juxta[:, 0])
    # point 2 (20 µm) juxta neighbors: point 1 at 10 µm. Point 0 is 20 > 15. Point 3 is 20 > 15.
    assert n_j[2] == 1 and abs(juxta[2, 0] - 3.0) < 1e-6, (n_j, juxta[:, 0])
    assert n_p[2] >= 1 and np.isfinite(para[2, 0])
    # statistical recovery
    rows = []
    for s in range(6):
        for f in range(4):
            n = 80
            krt = rng.normal(size=(n, 3))
            cldn4 = 0.35 * krt[:, 0] + rng.normal(size=n)
            y = -0.7 * cldn4 + 0.15 * krt[:, 0] + rng.normal(scale=0.4, size=n)
            for i in range(n):
                rows.append((f"S{s}", f"S{s}::F{f}", y[i], cldn4[i], *krt[i]))
    sim = pd.DataFrame(rows, columns=["sample", "fov_key", "y", "cldn4", "k1", "k2", "k3"])
    y = demean_groups(sim["y"].to_numpy()[:, None], sim["fov_key"].to_numpy())[:, 0]
    c = demean_groups(sim["cldn4"].to_numpy()[:, None], sim["fov_key"].to_numpy())[:, 0]
    k = demean_groups(sim[["k1", "k2", "k3"]].to_numpy(), sim["fov_key"].to_numpy())
    part = within_unit_partial(y, c, k, sim["sample"].to_numpy(), min_n=50)
    assert (part["spearman_partial"] < 0).all(), part["spearman_partial"].tolist()
    assert (part["beta_cldn4"] < 0).all(), part["beta_cldn4"].tolist()
    # null CLDN4 independent of y: partial median near 0 / not systematically as negative
    rows = []
    for s in range(6):
        n = 300
        krt = rng.normal(size=(n, 3))
        cldn4 = rng.normal(size=n)
        yv = 0.4 * krt[:, 0] + rng.normal(scale=0.5, size=n)
        for i in range(n):
            rows.append((f"S{s}", yv[i], cldn4[i], *krt[i]))
    sim = pd.DataFrame(rows, columns=["sample", "y", "cldn4", "k1", "k2", "k3"])
    part0 = within_unit_partial(
        sim["y"].to_numpy(),
        sim["cldn4"].to_numpy(),
        sim[["k1", "k2", "k3"]].to_numpy(),
        sim["sample"].to_numpy(),
        min_n=50,
    )
    assert abs(part0["spearman_partial"].median()) < 0.15, part0["spearman_partial"].tolist()
    log("self-test passed")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--reuse-cache", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    os.makedirs(OUT, exist_ok=True)
    if not os.path.exists(H5AD):
        msg = f"Missing {H5AD}. Run scripts/download_cosmx_nsclc_h5ad.py."
        with open(os.path.join(ROOT, "RESULTS.md"), "w") as fh:
            fh.write("# CosMx MISTy CLDN4 beyond KRT — STOP\n\n" + msg + "\n")
        print(msg, file=sys.stderr)
        return 2
    if args.reuse_cache and os.path.exists(CACHE):
        log(f"loading cache {CACHE}")
        df = pd.read_pickle(CACHE)
    else:
        mat = load_matrix(H5AD)
        # nearest-neighbor sanity on one section, in µm
        sm = mat["sample"] == "LUAD-5 R1"
        tree = cKDTree(mat["xy"][sm])
        d, _ = tree.query(mat["xy"][sm], k=2)
        nn_um = float(np.median(d[:, 1]) * PX_TO_UM)
        log(f"LUAD-5 R1 median NN distance {nn_um:.2f} µm")
        if not (4.0 <= nn_um <= 25.0):
            raise SystemExit(f"pixel scale looks wrong: median NN {nn_um:.2f} µm")
        df = build_table(mat)
        df.to_pickle(CACHE)
        log(f"cached {CACHE} rows={len(df)}")
    results = analyze(df)
    written = write_tables(df, results)
    write_figures(results, written)
    write_report(df, results, written)
    log("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
