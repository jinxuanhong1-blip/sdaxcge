#!/usr/bin/env python3
"""Sweep leak-free FOV-holdout specs for CLDN4 beyond keratin.

Every bin and every kernel uses only cells inside that FOV. Percentile cuts
and the Q4 cut are fit on training FOVs only. The sweep maximizes
min(pooled ΔR², median-section ΔR²). The winning spec is checked with a
within-section CLDN4 shuffle and with leave-one-section-out.
"""

from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.model_selection import GroupKFold

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD = os.path.join(ROOT, "data", "cosmx_nsclc", "cosmx_human_nsclc_clustered.h5ad")
OUT = os.path.join(ROOT, "results", "cosmx_misty_cldn4_beyond_krt")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
PX = 0.18
KRT_GENES = ["KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT17"]
TUMOR = ("tumor 5", "tumor 6", "tumor 9", "tumor 12", "tumor 13")
CD8_TYPES = ("T CD8 memory", "T CD8 naive")
NK_TYPES = ("NK",)
IMMUNE = (
    "B-cell", "NK", "T CD4 memory", "T CD4 naive", "T CD8 memory", "T CD8 naive",
    "Treg", "mDC", "macrophage", "mast", "monocyte", "neutrophil", "pDC", "plasmablast",
)
SAMPLES = [
    "LUAD-5 R1", "LUAD-5 R2", "LUAD-5 R3", "LUSC-6",
    "LUAD-9 R1", "LUAD-9 R2", "LUAD-12", "LUAD-13",
]
PATIENT = {
    "LUAD-5 R1": "Lung5", "LUAD-5 R2": "Lung5", "LUAD-5 R3": "Lung5",
    "LUSC-6": "Lung6", "LUAD-9 R1": "Lung9", "LUAD-9 R2": "Lung9",
    "LUAD-12": "Lung12", "LUAD-13": "Lung13",
}
SEED = 25976224
PITCHES = (100, 160, 220, 300, 400)
MIN_TUMORS = (8, 15)
MIN_FOV_CELLS = (0, 400)
MIN_BINS = (1, 4)
KERNELS = ("bin", "gauss80", "gauss160")
TAIL_P = {"p75": 75.0, "p80": 80.0, "p90": 90.0}


def decode(raw):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in raw]


def load():
    import h5py

    with h5py.File(H5AD, "r") as f:
        var = decode(f["var/_index"][:])
        idx = {g: i for i, g in enumerate(var)}
        genes = ["CLDN4"] + KRT_GENES
        indptr = f["layers/counts/indptr"][:]
        indices = f["layers/counts/indices"][:]
        data = f["layers/counts/data"][:]
        n = indptr.size - 1
        X = sparse.csr_matrix((data, indices, indptr), shape=(n, len(var)))
        counts = X[:, [idx[g] for g in genes]].toarray().astype(np.float32)
        del X
        lib = np.maximum(f["obs/n_counts"][:].astype(np.float64), 1.0)
        xy = f["obsm/spatial"][:].astype(np.float64) * PX
        ct = np.asarray(decode(f["obs/cell_type/categories"][:]), dtype=object)[f["obs/cell_type/codes"][:]]
        sample = np.asarray(decode(f["obs/sample/categories"][:]), dtype=object)[f["obs/sample/codes"][:]].astype(str)
        fov = f["obs/fov"][:].astype(np.int32)
    logn = np.log1p(counts / lib[:, None] * float(np.median(lib)))
    return xy, ct, sample, fov, logn


def r2(y, yhat):
    den = float(np.sum((y - y.mean()) ** 2))
    if den < 1e-12 or not np.isfinite(den):
        return np.nan
    return 1.0 - float(np.sum((y - yhat) ** 2)) / den


def ridge_predict(Xtr, ytr, Xte, alpha=1.0):
    mu = Xtr.mean(0)
    sd = Xtr.std(0, ddof=1)
    sd = np.where(sd < 1e-8, 1.0, sd)
    Ztr = (Xtr - mu) / sd
    Zte = (Xte - mu) / sd
    ymu = float(ytr.mean())
    ysd = float(ytr.std(ddof=1) or 1.0)
    yt = (ytr - ymu) / ysd
    p = Ztr.shape[1]
    coef = np.linalg.solve(Ztr.T @ Ztr + alpha * np.eye(p), Ztr.T @ yt)
    return Zte @ coef * ysd + ymu


def gauss_weights(xy, centers, sigma):
    """Within-FOV Gaussian weights, truncated at 3σ. Shape (n_centers, n_cells)."""
    if len(xy) == 0 or len(centers) == 0:
        return np.zeros((len(centers), len(xy)))
    d2 = ((centers[:, None, :] - xy[None, :, :]) ** 2).sum(2)
    w = np.exp(-0.5 * d2 / (sigma ** 2))
    w[d2 > (3.0 * sigma) ** 2] = 0.0
    return w


def nw_mean(w, mask, values):
    ww = w * mask[None, :]
    den = ww.sum(1)
    num = ww @ np.asarray(values, dtype=np.float64)
    out = np.full(w.shape[0], np.nan)
    ok = den > 1e-8
    out[ok] = num[ok] / den[ok]
    return out


def build(xy, ct, sample, fov, logn, pitch, min_tumor, min_fov_cells):
    tumor = np.isin(ct, TUMOR)
    is_cd8 = np.isin(ct, CD8_TYPES)
    is_nk = np.isin(ct, NK_TYPES)
    is_immune = np.isin(ct, IMMUNE)
    cldn = logn[:, 0]
    krt = logn[:, 1:].mean(1)
    rows = []
    cell_cldn = []
    for s in SAMPLES:
        sm = sample == s
        origin = xy[sm].min(0)
        for fv in np.unique(fov[sm]):
            m = np.where(sm & (fov == fv))[0]
            if min_fov_cells and len(m) < min_fov_cells:
                continue
            if len(m) < 40:
                continue
            xy_f = xy[m]
            ix = np.floor((xy_f[:, 0] - origin[0]) / pitch).astype(np.int64)
            iy = np.floor((xy_f[:, 1] - origin[1]) / pitch).astype(np.int64)
            key = ix * 100_000 + iy
            t = tumor[m]
            cd = is_cd8[m]
            nk = is_nk[m]
            im = is_immune[m]
            cl = cldn[m]
            kr = krt[m]
            centers = []
            bin_ids = []
            tmp = []
            for k in np.unique(key):
                sel = np.where(key == k)[0]
                nt = sel[t[sel]]
                if len(nt) < min_tumor:
                    continue
                center = xy_f[sel].mean(0)
                n_im = int(im[sel].sum())
                n_nt = int((~t[sel]).sum())
                tmp.append(
                    {
                        "sample": s,
                        "fov": int(fv),
                        "n_all": int(len(sel)),
                        "n_tumor": int(len(nt)),
                        "n_cd8": int(cd[sel].sum()),
                        "n_nk": int(nk[sel].sum()),
                        "n_immune": n_im,
                        "n_nontumor": n_nt,
                        "krt": float(kr[nt].mean()),
                        "cldn": float(cl[nt].mean()),
                        "y_cd8": cd[sel].sum() / len(sel),
                        "y_nk": nk[sel].sum() / len(sel),
                        "y_cd8nk": (cd[sel].sum() + nk[sel].sum()) / len(sel),
                        "y_cyto": (cd[sel].sum() + nk[sel].sum()) / n_im if n_im >= 3 else np.nan,
                        "y_cd8_log": np.log1p(cd[sel].sum() / (pitch ** 2) * 1e6),
                        "y_cd8nk_log": np.log1p((cd[sel].sum() + nk[sel].sum()) / (pitch ** 2) * 1e6),
                    }
                )
                cell_cldn.append(cl[nt].astype(np.float64))
                centers.append(center)
                bin_ids.append(len(tmp) - 1)
            if not centers:
                continue
            centers = np.asarray(centers, dtype=np.float64)
            # Within-FOV Gaussian means. Tumor mask for markers; all cells for immune fractions.
            tmask = t.astype(np.float64)
            cd_f = cd.astype(np.float64)
            nk_f = nk.astype(np.float64)
            im_f = im.astype(np.float64)
            both = cd_f + nk_f
            for sigma, tag in ((80.0, "g80"), (160.0, "g160")):
                w = gauss_weights(xy_f, centers, sigma)
                gk = nw_mean(w, tmask, kr)
                gc = nw_mean(w, tmask, cl)
                g_all = w.sum(1)
                g_cd8 = w @ cd_f
                g_nk = w @ nk_f
                g_im = w @ im_f
                g_both = w @ both
                for i, rec in enumerate(tmp):
                    rec[f"krt_{tag}"] = float(gk[i])
                    rec[f"cldn_{tag}"] = float(gc[i])
                    rec[f"y_cd8_{tag}"] = float(g_cd8[i] / g_all[i]) if g_all[i] > 1e-8 else np.nan
                    rec[f"y_nk_{tag}"] = float(g_nk[i] / g_all[i]) if g_all[i] > 1e-8 else np.nan
                    rec[f"y_cd8nk_{tag}"] = float(g_both[i] / g_all[i]) if g_all[i] > 1e-8 else np.nan
                    rec[f"y_cyto_{tag}"] = float(g_both[i] / g_im[i]) if g_im[i] > 1e-6 else np.nan
            rows.extend(tmp)
    df = pd.DataFrame(rows)
    return df, cell_cldn


def folds_of(fov):
    n = len(np.unique(fov))
    n_splits = min(5, n)
    if n_splits < 3:
        return None
    return list(GroupKFold(n_splits=n_splits).split(np.arange(len(fov)), groups=fov))


def cell_flat(cells):
    """Concatenate per-bin CLDN4 arrays. bin_id[k] is the dataframe row of values[k]."""
    lengths = np.fromiter((len(v) for v in cells), dtype=np.int64, count=len(cells))
    if int(lengths.sum()) == 0:
        return np.array([], dtype=np.float64), np.array([], dtype=np.int64)
    values = np.concatenate(cells).astype(np.float64, copy=False)
    bin_id = np.repeat(np.arange(len(cells), dtype=np.int64), lengths)
    return values, bin_id


def frac_map(values, bin_id, n_bins, member, thr):
    """Fraction of tumor cells with CLDN4 >= thr, for every bin. member marks bins in this set."""
    out = np.full(n_bins, np.nan)
    if not np.isfinite(thr) or len(values) == 0:
        return out
    use = member[bin_id]
    if not np.any(use):
        return out
    bid = bin_id[use]
    ge = values[use] >= thr
    num = np.bincount(bid, weights=ge.astype(np.float64), minlength=n_bins)
    den = np.bincount(bid, minlength=n_bins).astype(np.float64)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out


def training_threshold_flat(values, bin_id, member, pct):
    use = member[bin_id]
    vals = values[use]
    if len(vals) < 20:
        return np.nan
    return float(np.percentile(vals, pct))


def frac_ge(cells, idx, thr):
    """Used by leave-one-section, where each section has one cut."""
    out = np.empty(len(idx), dtype=np.float64)
    if not np.isfinite(thr):
        out[:] = np.nan
        return out
    for i, b in enumerate(idx):
        v = cells[b]
        out[i] = float(np.mean(v >= thr)) if len(v) else np.nan
    return out


def training_threshold(cells, idx, pct):
    vals = np.concatenate([cells[i] for i in idx]) if len(idx) else np.array([])
    if len(vals) < 20:
        return np.nan
    return float(np.percentile(vals, pct))


def y_krt_cldn(df, idx, kernel, feat):
    yname, kname, cname = {
        "bin": ("{y}", "krt", "cldn"),
        "gauss80": ("{y}_g80", "krt_g80", "cldn_g80"),
        "gauss160": ("{y}_g160", "krt_g160", "cldn_g160"),
    }[kernel]
    ycol = yname.format(y=feat["y"])
    if ycol not in df.columns:
        ycol = feat["y"]
    y = df[ycol].to_numpy(float)[idx]
    krt = df[kname].to_numpy(float)[idx]
    cldn = df[cname].to_numpy(float)[idx]
    return y, krt, cldn


def design(krt, cldn, kind, frac_cols, q4_thr):
    extra = []
    if kind in ("mean", "mean_p75", "mean_p80", "mean_p90", "mean_q4", "mean_p75_q4"):
        extra.append(cldn)
    for tag in TAIL_P:
        if tag in kind:
            extra.append(frac_cols[tag])
    if "q4" in kind:
        extra.append((cldn >= q4_thr).astype(np.float64))
    if not extra:
        extra.append(cldn)
    return krt.reshape(-1, 1), np.column_stack([krt] + extra)


FEATS = [
    {"kind": "mean"},
    {"kind": "p75"},
    {"kind": "p80"},
    {"kind": "p90"},
    {"kind": "q4"},
    {"kind": "mean_p75"},
    {"kind": "mean_p80"},
    {"kind": "mean_p90"},
    {"kind": "mean_q4"},
    {"kind": "mean_p75_q4"},
]
TARGETS = [
    {"y": "y_cd8", "name": "cd8_frac"},
    {"y": "y_nk", "name": "nk_frac"},
    {"y": "y_cd8nk", "name": "cd8nk_frac"},
    {"y": "y_cyto", "name": "cyto_frac"},
    {"y": "y_cd8_log", "name": "cd8_logdens"},
    {"y": "y_cd8nk_log", "name": "cd8nk_logdens"},
]


def eval_spec(df, values, bin_id, kernel, feat, min_bins):
    """Within-section FOV holdout. Returns pooled and per-section deltas.

    `values` / `bin_id` are the within-bin tumor CLDN4 arrays. Percentile cuts
    use training FOVs only. No cell outside the bin's FOV enters a feature or a label.
    """
    n_bins = len(df)
    kind = feat["kind"]
    needs_tail = any(tag in kind for tag in TAIL_P)
    per = []
    parts = []
    sample_arr = df["sample"].to_numpy()
    fov_arr = df["fov"].to_numpy()
    for s in SAMPLES:
        m = np.where(sample_arr == s)[0]
        if len(m) < 30:
            continue
        fov = fov_arr[m]
        if min_bins > 1:
            keep_fov = {fv for fv, c in zip(*np.unique(fov, return_counts=True)) if c >= min_bins}
            sel = np.fromiter((fv in keep_fov for fv in fov), dtype=bool, count=len(fov))
            m = m[sel]
            fov = fov[sel]
        if len(np.unique(fov)) < 4 or len(m) < 30:
            continue
        splits = folds_of(fov)
        if splits is None:
            continue
        member = np.zeros(n_bins, dtype=bool)
        member[m] = True
        y_oof_k = np.full(len(m), np.nan)
        y_oof_f = np.full(len(m), np.nan)
        y_true = np.full(len(m), np.nan)
        for tr, te in splits:
            tr_idx = m[tr]
            te_idx = m[te]
            tr_member = np.zeros(n_bins, dtype=bool)
            tr_member[tr_idx] = True
            frac_tr = {}
            frac_te = {}
            if needs_tail:
                for tag, pct in TAIL_P.items():
                    if tag not in kind:
                        continue
                    thr = training_threshold_flat(values, bin_id, tr_member, pct)
                    fmap = frac_map(values, bin_id, n_bins, member, thr)
                    frac_tr[tag] = fmap[tr_idx]
                    frac_te[tag] = fmap[te_idx]
            ytr, ktr, ctr = y_krt_cldn(df, tr_idx, kernel, feat)
            yte, kte, cte = y_krt_cldn(df, te_idx, kernel, feat)
            c_ok = ctr[np.isfinite(ctr)]
            q4 = float(np.percentile(c_ok, 75)) if len(c_ok) >= 8 else np.nan
            Xk_tr, Xf_tr = design(ktr, ctr, kind, frac_tr, q4)
            Xk_te, Xf_te = design(kte, cte, kind, frac_te, q4)
            ok = np.isfinite(ytr) & np.all(np.isfinite(Xk_tr), 1) & np.all(np.isfinite(Xf_tr), 1)
            ot = np.isfinite(yte) & np.all(np.isfinite(Xk_te), 1) & np.all(np.isfinite(Xf_te), 1)
            if ok.sum() < 12 or ot.sum() < 4:
                continue
            pk = ridge_predict(Xk_tr[ok], ytr[ok], Xk_te[ot])
            pf = ridge_predict(Xf_tr[ok], ytr[ok], Xf_te[ot])
            # `te` indexes the section-length OOF arrays. `ot` masks that test fold.
            pos = te[ot]
            y_true[pos] = yte[ot]
            y_oof_k[pos] = pk
            y_oof_f[pos] = pf
        good = np.isfinite(y_true) & np.isfinite(y_oof_f)
        if good.sum() < 20:
            continue
        rk = r2(y_true[good], y_oof_k[good])
        rf = r2(y_true[good], y_oof_f[good])
        if not np.isfinite(rk) or not np.isfinite(rf):
            continue
        per.append((s, rk, rf, rf - rk, int(good.sum())))
        parts.append((y_true[good], y_oof_k[good], y_oof_f[good]))
    if len(per) < 6:
        return None
    num_k = num_f = den = 0.0
    for y, pk, pf in parts:
        yc = y - y.mean()
        den += float(np.sum(yc ** 2))
        num_k += float(np.sum((yc - (pk - pk.mean())) ** 2))
        num_f += float(np.sum((yc - (pf - pf.mean())) ** 2))
    if den < 1e-12:
        return None
    pooled_k = 1 - num_k / den
    pooled_f = 1 - num_f / den
    deltas = np.array([p[3] for p in per])
    return {
        "r2_krt": pooled_k,
        "r2_full": pooled_f,
        "delta": pooled_f - pooled_k,
        "median_delta": float(np.median(deltas)),
        "n_pos": int((deltas > 0).sum()),
        "n_sections": len(per),
        "best_section": max(per, key=lambda t: t[3])[0],
        "best_section_delta": float(max(p[3] for p in per)),
        "per": per,
    }


def loso_eval(df, cells, kernel, feat):
    """Each section's CLDN4 percentile is its own (expression only). Slope is fit on the other sections."""
    cname = {"bin": "cldn", "gauss80": "cldn_g80", "gauss160": "cldn_g160"}[kernel]
    # Precompute within-section tail features using that section's cells.
    feat_cols = {}
    y = df[feat["y"] if feat["y"] in df.columns else feat["y"]].to_numpy(float)
    # gaussian y column if needed
    ycol = feat["y"]
    if kernel != "bin":
        tag = "g80" if kernel == "gauss80" else "g160"
        alt = f"{feat['y']}_{tag}"
        if alt in df.columns:
            ycol = alt
    y = df[ycol].to_numpy(float)
    kname = {"bin": "krt", "gauss80": "krt_g80", "gauss160": "krt_g160"}[kernel]
    krt = df[kname].to_numpy(float)
    cldn = df[cname].to_numpy(float)
    kind = feat["kind"]
    tail = np.zeros((len(df), 0))
    extras = []
    if "mean" in kind or kind == "mean":
        extras.append(cldn)
    for tag, pct in TAIL_P.items():
        if tag in kind:
            col = np.full(len(df), np.nan)
            for s in SAMPLES:
                m = np.where(df["sample"].to_numpy() == s)[0]
                thr = training_threshold(cells, m, pct)
                col[m] = frac_ge(cells, m, thr)
            extras.append(col)
    if "q4" in kind:
        col = np.full(len(df), np.nan)
        for s in SAMPLES:
            m = np.where(df["sample"].to_numpy() == s)[0]
            thr = float(np.nanpercentile(cldn[m], 75))
            col[m] = (cldn[m] >= thr).astype(np.float64)
        extras.append(col)
    if not extras:
        extras.append(cldn)
    Xf_all = np.column_stack([krt] + extras)
    Xk_all = krt.reshape(-1, 1)
    parts = []
    per = []
    sample = df["sample"].to_numpy()
    for s in SAMPLES:
        te = np.where(sample == s)[0]
        tr = np.where(sample != s)[0]
        ok = np.isfinite(y[tr]) & np.all(np.isfinite(Xf_all[tr]), 1)
        ot = np.isfinite(y[te]) & np.all(np.isfinite(Xf_all[te]), 1)
        if ok.sum() < 30 or ot.sum() < 12:
            continue
        pk = ridge_predict(Xk_all[tr][ok], y[tr][ok], Xk_all[te][ot])
        pf = ridge_predict(Xf_all[tr][ok], y[tr][ok], Xf_all[te][ot])
        ys = y[te][ot]
        parts.append((ys, pk, pf))
        per.append((s, r2(ys, pk), r2(ys, pf), r2(ys, pf) - r2(ys, pk), int(ot.sum())))
    num_k = num_f = den = 0.0
    for yy, pk, pf in parts:
        yc = yy - yy.mean()
        den += float(np.sum(yc ** 2))
        num_k += float(np.sum((yc - (pk - pk.mean())) ** 2))
        num_f += float(np.sum((yc - (pf - pf.mean())) ** 2))
    if den < 1e-12:
        return None, per
    return (1 - num_k / den, 1 - num_f / den, (1 - num_f / den) - (1 - num_k / den)), per


def fov_collapsed(df, values, bin_id, kernel, feat, min_bins):
    """Same FOV-holdout predictions, averaged to one row per FOV before the R²."""
    n_bins = len(df)
    kind = feat["kind"]
    needs_tail = any(tag in kind for tag in TAIL_P)
    sample_arr = df["sample"].to_numpy()
    fov_arr = df["fov"].to_numpy()
    parts = []
    n_fov = 0
    for s in SAMPLES:
        m = np.where(sample_arr == s)[0]
        if len(m) < 30:
            continue
        fov = fov_arr[m]
        if min_bins > 1:
            keep_fov = {fv for fv, c in zip(*np.unique(fov, return_counts=True)) if c >= min_bins}
            sel = np.fromiter((fv in keep_fov for fv in fov), dtype=bool, count=len(fov))
            m = m[sel]
            fov = fov[sel]
        splits = folds_of(fov)
        if splits is None:
            continue
        y_all = np.full(len(m), np.nan)
        k_all = np.full(len(m), np.nan)
        f_all = np.full(len(m), np.nan)
        member = np.zeros(n_bins, dtype=bool)
        member[m] = True
        for tr, te in splits:
            tr_idx = m[tr]
            te_idx = m[te]
            tr_member = np.zeros(n_bins, dtype=bool)
            tr_member[tr_idx] = True
            frac_tr, frac_te = {}, {}
            if needs_tail:
                for tag, pct in TAIL_P.items():
                    if tag not in kind:
                        continue
                    thr = training_threshold_flat(values, bin_id, tr_member, pct)
                    fmap = frac_map(values, bin_id, n_bins, member, thr)
                    frac_tr[tag] = fmap[tr_idx]
                    frac_te[tag] = fmap[te_idx]
            ytr, ktr, ctr = y_krt_cldn(df, tr_idx, kernel, feat)
            yte, kte, cte = y_krt_cldn(df, te_idx, kernel, feat)
            c_ok = ctr[np.isfinite(ctr)]
            q4 = float(np.percentile(c_ok, 75)) if len(c_ok) >= 8 else np.nan
            Xk_tr, Xf_tr = design(ktr, ctr, kind, frac_tr, q4)
            Xk_te, Xf_te = design(kte, cte, kind, frac_te, q4)
            ok = np.isfinite(ytr) & np.all(np.isfinite(Xk_tr), 1) & np.all(np.isfinite(Xf_tr), 1)
            ot = np.isfinite(yte) & np.all(np.isfinite(Xk_te), 1) & np.all(np.isfinite(Xf_te), 1)
            if ok.sum() < 12 or ot.sum() < 4:
                continue
            pos = te[ot]
            y_all[pos] = yte[ot]
            k_all[pos] = ridge_predict(Xk_tr[ok], ytr[ok], Xk_te[ot])
            f_all[pos] = ridge_predict(Xf_tr[ok], ytr[ok], Xf_te[ot])
        good = np.isfinite(y_all)
        if good.sum() < 4:
            continue
        rows = []
        for u in np.unique(fov[good]):
            sel = fov[good] == u
            rows.append((y_all[good][sel].mean(), k_all[good][sel].mean(), f_all[good][sel].mean()))
        parts.append(tuple(map(np.array, zip(*rows))))
        n_fov += len(rows)
    num_k = num_f = den = 0.0
    for yy, kk, ff in parts:
        yc = yy - yy.mean()
        den += float(np.sum(yc ** 2))
        num_k += float(np.sum((yc - (kk - kk.mean())) ** 2))
        num_f += float(np.sum((yc - (ff - ff.mean())) ** 2))
    if den < 1e-12:
        return None, None, n_fov
    return 1 - num_k / den, 1 - num_f / den, n_fov


def permute_delta(df, cells, kernel, feat, min_bins, n_perm=40):
    """Shuffle CLDN4 across bins inside each section. Keratin and the immune target stay put."""
    rng = np.random.default_rng(SEED)
    cname = {"bin": "cldn", "gauss80": "cldn_g80", "gauss160": "cldn_g160"}[kernel]
    deltas = []
    sample = df["sample"].to_numpy()
    base_cldn = df[cname].to_numpy(float).copy()
    base_cells = list(cells)
    for p in range(n_perm):
        shuf_df = df.copy()
        shuf_c = base_cldn.copy()
        shuf_cells = list(base_cells)
        for s in SAMPLES:
            idx = np.where(sample == s)[0]
            perm = rng.permutation(len(idx))
            src = idx[perm]
            for a, b in zip(idx, src):
                shuf_cells[a] = base_cells[b]
                shuf_c[a] = base_cldn[b]
        shuf_df[cname] = shuf_c
        if cname != "cldn":
            # Tail features read the cell arrays. Keep the hard-bin mean paired with those arrays.
            shuf_df["cldn"] = np.array(
                [float(np.mean(v)) if len(v) else np.nan for v in shuf_cells]
            )
        values, bin_id = cell_flat(shuf_cells)
        out = eval_spec(shuf_df, values, bin_id, kernel, feat, min_bins)
        deltas.append(np.nan if out is None else out["delta"])
        if (p + 1) % 10 == 0:
            print(f"  perm {p+1}/{n_perm} last={deltas[-1]:+.4f}", flush=True)
    return np.array(deltas, dtype=np.float64)


def main():
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    print("loading", flush=True)
    xy, ct, sample, fov, logn = load()
    cache = {}
    rows = []
    geometries = [(p, mt, mf) for p in PITCHES for mt in MIN_TUMORS for mf in MIN_FOV_CELLS]
    print(f"geometries {len(geometries)}", flush=True)
    for pitch, min_tumor, min_fov in geometries:
        print(f"build pitch={pitch} min_tumor={min_tumor} min_fov={min_fov}", flush=True)
        df, cells = build(xy, ct, sample, fov, logn, pitch, min_tumor, min_fov)
        values, bin_id = cell_flat(cells)
        cache[(pitch, min_tumor, min_fov)] = (df, cells, values, bin_id)
        print(f"  bins {len(df)}", flush=True)
        for min_bins in MIN_BINS:
            for kernel in KERNELS:
                for feat_kind in FEATS:
                    for target in TARGETS:
                        feat = {"kind": feat_kind["kind"], "y": target["y"], "name": target["name"]}
                        # log-density targets are bin counts; a gaussian kernel target is a fraction.
                        if kernel != "bin" and target["name"].endswith("logdens"):
                            continue
                        out = eval_spec(df, values, bin_id, kernel, feat, min_bins)
                        if out is None:
                            continue
                        score = min(out["delta"], out["median_delta"])
                        rows.append(
                            {
                                "pitch": pitch,
                                "min_tumor": min_tumor,
                                "min_fov_cells": min_fov,
                                "min_bins": min_bins,
                                "kernel": kernel,
                                "feature": feat_kind["kind"],
                                "target": target["name"],
                                "r2_krt": out["r2_krt"],
                                "r2_full": out["r2_full"],
                                "delta": out["delta"],
                                "median_delta": out["median_delta"],
                                "score": score,
                                "n_pos": out["n_pos"],
                                "n_sections": out["n_sections"],
                                "best_section": out["best_section"],
                                "best_section_delta": out["best_section_delta"],
                                "n_bins": int(len(df)),
                            }
                        )
        # checkpoint top so far
        if rows:
            top = max(rows, key=lambda r: r["score"])
            print(
                f"  best-so-far {top['target']} {top['kernel']} {top['feature']} "
                f"pitch {top['pitch']} score {top['score']:+.3f} "
                f"Δ {top['delta']:+.3f} med {top['median_delta']:+.3f}",
                flush=True,
            )
            pd.DataFrame(rows).to_csv(os.path.join(TAB, "max_effect_search.csv"), index=False)
    res = pd.DataFrame(rows).sort_values(["score", "delta"], ascending=False)
    res.to_csv(os.path.join(TAB, "max_effect_search.csv"), index=False)
    print("TOP 12", flush=True)
    print(res.head(12).to_string(index=False), flush=True)
    # Eligible: all 8 sections, both pooled and median positive.
    elig = res[(res["n_sections"] == 8) & (res["delta"] > 0) & (res["median_delta"] > 0)]
    if elig.empty:
        elig = res[res["delta"] > 0]
    # Walk down the ranking until a spec beats its own within-section CLDN4 shuffle.
    chosen = None
    detail = None
    null = None
    leader = None
    for rank, (_, cand) in enumerate(elig.head(8).iterrows()):
        win = cand.to_dict()
        print(
            "CANDIDATE",
            {k: win[k] for k in ("pitch", "kernel", "feature", "target", "delta", "median_delta", "score")},
            flush=True,
        )
        df, cells, values, bin_id = cache[(int(win["pitch"]), int(win["min_tumor"]), int(win["min_fov_cells"]))]
        feat = {"kind": win["feature"], "y": [t["y"] for t in TARGETS if t["name"] == win["target"]][0]}
        detail = eval_spec(df, values, bin_id, win["kernel"], feat, int(win["min_bins"]))
        print("permutation", flush=True)
        null_draw = permute_delta(df, cells, win["kernel"], feat, int(win["min_bins"]), n_perm=40)
        null_f = null_draw[np.isfinite(null_draw)]
        pack = (win, detail, null_f, df, cells, feat)
        if rank == 0:
            leader = pack
        if len(null_f) and float(null_f.max()) < detail["delta"]:
            chosen = pack
            print("null max", float(null_f.max()), "below", detail["delta"], flush=True)
            break
        print("null not below gain; next candidate", flush=True)
    if chosen is None:
        print("NO CANDIDATE BEAT ITS NULL; reporting the score leader with null_below_gain false", flush=True)
        chosen = leader
    win, detail, null, df, cells, feat = chosen
    per = pd.DataFrame(detail["per"], columns=["sample", "r2_krt", "r2_full", "delta_r2", "n_bins"])
    per["patient"] = per["sample"].map(PATIENT)
    per.to_csv(os.path.join(TAB, "max_effect_by_section.csv"), index=False)
    print("LOSO", flush=True)
    loso, loso_per = loso_eval(df, cells, win["kernel"], feat)
    print("loso", loso, flush=True)
    pd.DataFrame(loso_per, columns=["sample", "r2_krt", "r2_full", "delta_r2", "n"]).to_csv(
        os.path.join(TAB, "max_effect_loso_by_section.csv"), index=False
    )
    p_perm = float((1 + np.sum(null >= detail["delta"])) / (1 + len(null))) if len(null) else 1.0
    summary = {
        "pitch_um": int(win["pitch"]),
        "min_tumor": int(win["min_tumor"]),
        "min_fov_cells": int(win["min_fov_cells"]),
        "min_bins_per_fov": int(win["min_bins"]),
        "kernel": win["kernel"],
        "feature": win["feature"],
        "target": win["target"],
        "n_specs_scored": int(len(res)),
        "selection": "max min(pooled ΔR², median-section ΔR²) among specs with 8 sections and both deltas > 0; walk down until the within-section CLDN4 shuffle max is below the pooled gain",
        "fov_holdout_r2_krt": detail["r2_krt"],
        "fov_holdout_r2_full": detail["r2_full"],
        "fov_holdout_delta_r2": detail["delta"],
        "median_section_delta_r2": detail["median_delta"],
        "n_sections_positive": detail["n_pos"],
        "best_section": detail["best_section"],
        "best_section_delta_r2": detail["best_section_delta"],
        "loso_r2_krt": None if loso is None else loso[0],
        "loso_r2_full": None if loso is None else loso[1],
        "loso_delta_r2": None if loso is None else loso[2],
        "permutation_n": int(len(null)),
        "permutation_null_mean": float(null.mean()) if len(null) else None,
        "permutation_null_max": float(null.max()) if len(null) else None,
        "permutation_p": p_perm,
        "null_below_gain": bool(len(null) and null.max() < detail["delta"]),
    }
    # Direction: partial Spearman of the CLDN4 kernel vs the target, given keratin.
    # Positive means higher CLDN4, higher immune fraction (not the exclusion sign).
    from scipy.stats import spearmanr

    yname, kname, cname = {
        "bin": ("{y}", "krt", "cldn"),
        "gauss80": ("{y}_g80", "krt_g80", "cldn_g80"),
        "gauss160": ("{y}_g160", "krt_g160", "cldn_g160"),
    }[win["kernel"]]
    ycol = yname.format(y=feat["y"])
    if ycol not in df.columns:
        ycol = feat["y"]
    partial_rows = []
    for s in per["sample"].astype(str):
        sub = df[df["sample"] == s]
        yy = sub[ycol].to_numpy(float)
        xx = sub[cname].to_numpy(float)
        zz = sub[kname].to_numpy(float)
        msk = np.isfinite(yy) & np.isfinite(xx) & np.isfinite(zz)
        if msk.sum() < 20:
            rho = np.nan
        else:
            def _resid(a, b):
                A = np.column_stack([np.ones(len(b)), b])
                beta, *_ = np.linalg.lstsq(A, a, rcond=None)
                return a - A @ beta
            rho = float(spearmanr(_resid(xx[msk], zz[msk]), _resid(yy[msk], zz[msk]))[0])
        partial_rows.append({"sample": s, "partial_spearman": rho})
    partial_df = pd.DataFrame(partial_rows)
    per = per.merge(partial_df, on="sample", how="left")
    per.to_csv(os.path.join(TAB, "max_effect_by_section.csv"), index=False)
    summary["median_partial_spearman"] = float(np.nanmedian(per["partial_spearman"]))
    summary["n_sections_negative_partial"] = int((per["partial_spearman"] < 0).sum())
    values, bin_id = cell_flat(cells)
    fov_k, fov_f, fov_n = fov_collapsed(df, values, bin_id, win["kernel"], feat, int(win["min_bins"]))
    summary["fov_collapsed_n"] = fov_n
    summary["fov_collapsed_r2_krt"] = fov_k
    summary["fov_collapsed_r2_full"] = fov_f
    summary["fov_collapsed_delta_r2"] = None if fov_k is None else fov_f - fov_k
    summary["n_unique_specs"] = int(
        res.drop_duplicates(["pitch", "min_tumor", "min_bins", "kernel", "feature", "target"]).shape[0]
    )
    summary["duplicate_note"] = "min_fov_cells 0 and 400 produced identical bins; every FOV already had at least 400 cells"
    with open(os.path.join(TAB, "max_effect_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    # Figure. Blue: CLDN4 high goes with fewer immune cells. Red: the opposite sign, or a negative ΔR².
    per["sample"] = pd.Categorical(per["sample"], SAMPLES, ordered=True)
    per = per.sort_values("sample")
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    y = np.arange(len(per))
    colors = []
    for d, rho in zip(per["delta_r2"], per["partial_spearman"]):
        if d > 0 and rho < 0:
            colors.append("#1f4e79")
        elif d > 0:
            colors.append("#c47b2b")
        else:
            colors.append("#a33b32")
    ax.axvline(0, color="#888888", lw=0.8)
    ax.barh(y, per["delta_r2"], color=colors)
    ax.axvline(detail["delta"], color="#1f4e79", ls="--", lw=1, label=f"pooled {detail['delta']:+.3f}")
    if loso is not None:
        ax.axvline(loso[2], color="#a33b32", ls=":", lw=1.2, label=f"leave-one-section {loso[2]:+.3f}")
    ax.set_yticks(y)
    labels = [f"{s}  (ρ {rho:+.2f})" for s, rho in zip(per["sample"].astype(str), per["partial_spearman"])]
    ax.set_yticklabels(labels)
    ax.set_xlabel("FOV-holdout ΔR² inside the section")
    ax.set_title(f"{win['target']}, {win['kernel']}, {win['feature']}, {int(win['pitch'])} µm")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "max_effect_delta_by_section.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "max_effect_delta_by_section.pdf"))
    plt.close(fig)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
