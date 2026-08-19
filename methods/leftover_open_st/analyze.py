#!/usr/bin/env python3
"""CLDN4-only spatial exclusion vs CD8A on leftover public lung ST.

For each runnable section:
  - same-spot Spearman(CLDN4, CD8A)
  - partial Spearman residualizing KRT8
  - nearest CD8-high distance from CLDN4-high epithelial spots
  - neighbor CD8-high counts in a local radius

Skip if CLDN4 is absent. Do not write claim-failed pages.
"""

from __future__ import annotations

import json
import tarfile
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import io, stats
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "leftover_open_st"
OUT = ROOT / "results" / "leftover_open_st"
TABLES = OUT / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

MIN_GENES = 100
MIN_N_SPEARMAN = 30
CLDN4_Q = 0.75
CD8_Q = 0.75
EPI_Q = 0.50
NEI_RADIUS_BINS = 2


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if pd.isna(o):
        return None
    raise TypeError(type(o))


def dump_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, default=_json_default) + "\n")


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < MIN_N_SPEARMAN:
        return {"n": n, "rho": None, "p": None}
    rho, p = stats.spearmanr(x[m], y[m])
    if not np.isfinite(rho):
        return {"n": n, "rho": None, "p": None}
    return {"n": n, "rho": float(rho), "p": float(p)}


def partial_spearman(x, y, z):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 40:
        return {"n": n, "rho": None, "p": None}
    rx = stats.rankdata(x[m])
    ry = stats.rankdata(y[m])
    rz = stats.rankdata(z[m])
    rz = (rz - rz.mean()) / (rz.std() + 1e-12)
    bx = np.polyfit(rz, rx, 1)
    by = np.polyfit(rz, ry, 1)
    ex = rx - (bx[0] * rz + bx[1])
    ey = ry - (by[0] * rz + by[1])
    rho, p = stats.spearmanr(ex, ey)
    if not np.isfinite(rho):
        return {"n": n, "rho": None, "p": None}
    return {"n": n, "rho": float(rho), "p": float(p)}


def mw(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "n_a": int(a.size),
        "n_b": int(b.size),
        "median_a": float(np.median(a)) if a.size else None,
        "median_b": float(np.median(b)) if b.size else None,
        "p": None,
    }
    if a.size < 8 or b.size < 8:
        return out
    try:
        out["p"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except ValueError:
        out["p"] = None
    return out


def wilcoxon_signed(values):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    out = {
        "n": int(v.size),
        "median": float(np.median(v)) if v.size else None,
        "n_neg": int((v < 0).sum()) if v.size else 0,
        "n_pos": int((v > 0).sum()) if v.size else 0,
        "p": None,
    }
    if v.size < 3:
        return out
    try:
        out["p"] = float(stats.wilcoxon(v, alternative="two-sided").pvalue)
    except ValueError:
        out["p"] = None
    return out


def log1p_cp10k_from_counts(counts, lib=None):
    counts = np.asarray(counts, dtype=float)
    if lib is None:
        lib = counts  # unused placeholder
    lib = np.asarray(lib, dtype=float)
    lib = np.clip(lib, 1.0, None)
    return np.log1p(counts * (1e4 / lib))


def exclusion_bundle(name, cl4, cd8, krt8, x, y, lib, extra=None, radius=None):
    """CLDN4 vs CD8A spatial exclusion on one section."""
    extra = extra or {}
    n = int(cl4.size)
    if radius is None:
        # ~2 bins in coordinate units (array or pixels)
        if n >= 8:
            step = np.median(np.diff(np.unique(np.round(x, 6)))) if np.unique(x).size > 3 else 1.0
            if not np.isfinite(step) or step <= 0:
                step = 1.0
            radius = float(NEI_RADIUS_BINS * step)
        else:
            radius = 2.0

    cl4_log = log1p_cp10k_from_counts(cl4, lib)
    cd8_log = log1p_cp10k_from_counts(cd8, lib)
    krt_log = log1p_cp10k_from_counts(krt8, lib)

    if np.nanmedian(krt8) <= 0:
        epi = krt8 > 0
    else:
        epi_cut = np.nanquantile(krt_log, EPI_Q)
        epi = krt_log >= epi_cut
    if int(epi.sum()) < 40:
        epi = np.ones(n, dtype=bool)

    # Primary correlations on epithelial spots (avoids empty-spot ties).
    same = spearman(cl4_log[epi], cd8_log[epi])
    krt_resid = partial_spearman(cl4_log[epi], cd8_log[epi], krt_log[epi])
    same_all = spearman(cl4_log, cd8_log)

    cl4_epi = cl4_log[epi]
    hi_cut = np.nanquantile(cl4_epi, CLDN4_Q)
    if hi_cut <= 0:
        cl4_hi = epi & (cl4 > 0)
        cl4_lo = epi & (cl4 == 0)
    else:
        lo_cut = np.nanquantile(cl4_epi, 1.0 - CLDN4_Q)
        cl4_hi = epi & (cl4_log >= hi_cut)
        cl4_lo = epi & (cl4_log <= lo_cut)

    cd8_pos = cd8_log[cd8_log > 0]
    if cd8_pos.size >= 20:
        cd8_cut = np.nanquantile(cd8_pos, CD8_Q)
    else:
        cd8_cut = 0.0
    cd8_hi = cd8_log > max(cd8_cut, 0.0)
    if int(cd8_hi.sum()) < 8:
        cd8_hi = cd8 > 0

    xy = np.column_stack([np.asarray(x, float), np.asarray(y, float)])
    nearest_hi = np.full(n, np.nan)
    nearest_lo = np.full(n, np.nan)
    n_nei_hi = np.full(n, np.nan)
    n_nei_lo = np.full(n, np.nan)
    if int(cd8_hi.sum()) >= 1:
        tree = cKDTree(xy[cd8_hi])
        if int(cl4_hi.sum()):
            d, _ = tree.query(xy[cl4_hi], k=1, workers=-1)
            nearest_hi[cl4_hi] = d
            counts = tree.query_ball_point(xy[cl4_hi], r=radius, return_length=True)
            n_nei_hi[cl4_hi] = np.asarray(counts, float)
        if int(cl4_lo.sum()):
            d, _ = tree.query(xy[cl4_lo], k=1, workers=-1)
            nearest_lo[cl4_lo] = d
            counts = tree.query_ball_point(xy[cl4_lo], r=radius, return_length=True)
            n_nei_lo[cl4_lo] = np.asarray(counts, float)

    dist_mw = mw(nearest_hi, nearest_lo)
    nei_mw = mw(n_nei_hi, n_nei_lo)

    out = {
        "section": name,
        "n_spots": n,
        "n_epi": int(epi.sum()),
        "n_cldn4_hi_epi": int(cl4_hi.sum()),
        "n_cldn4_lo_epi": int(cl4_lo.sum()),
        "n_cd8_hi": int(cd8_hi.sum()),
        "frac_cldn4_pos": float((cl4 > 0).mean()),
        "frac_cd8a_pos": float((cd8 > 0).mean()),
        "same_spot_CLDN4_CD8A": same,
        "same_spot_CLDN4_CD8A_all": same_all,
        "krt8_residual_CLDN4_CD8A": krt_resid,
        "nearest_cd8_from_cldn4hi": {
            "median": float(np.nanmedian(nearest_hi)) if np.isfinite(nearest_hi).any() else None,
            "n": int(np.isfinite(nearest_hi).sum()),
        },
        "nearest_cd8_from_cldn4lo": {
            "median": float(np.nanmedian(nearest_lo)) if np.isfinite(nearest_lo).any() else None,
            "n": int(np.isfinite(nearest_lo).sum()),
        },
        "nearest_cd8_hi_vs_lo_mw": dist_mw,
        "neighbor_cd8_count_hi": {
            "median": float(np.nanmedian(n_nei_hi)) if np.isfinite(n_nei_hi).any() else None,
            "mean": float(np.nanmean(n_nei_hi)) if np.isfinite(n_nei_hi).any() else None,
            "n": int(np.isfinite(n_nei_hi).sum()),
        },
        "neighbor_cd8_count_lo": {
            "median": float(np.nanmedian(n_nei_lo)) if np.isfinite(n_nei_lo).any() else None,
            "mean": float(np.nanmean(n_nei_lo)) if np.isfinite(n_nei_lo).any() else None,
            "n": int(np.isfinite(n_nei_lo).sum()),
        },
        "neighbor_cd8_hi_vs_lo_mw": nei_mw,
        "radius": float(radius),
        **extra,
    }
    return out


def gene_index_from_names(names):
    idx = {}
    for i, g in enumerate(names):
        g = g.decode() if isinstance(g, (bytes, np.bytes_)) else str(g)
        idx.setdefault(g.upper(), i)
    return idx


def load_10x_h5(path):
    with h5py.File(path, "r") as f:
        mat = f["matrix"]
        data = mat["data"][:]
        indices = mat["indices"][:]
        indptr = mat["indptr"][:]
        shape = tuple(int(x) for x in mat["shape"][:])
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in mat["barcodes"][:]]
        genes = [x.decode() if isinstance(x, bytes) else str(x) for x in mat["features"]["name"][:]]
    mtx = csr_matrix((data, indices, indptr), shape=(shape[1], shape[0]))
    return pd.Series(barcodes, dtype=str), genes, mtx


def load_mtx_dir(d):
    d = Path(d)
    bc = next(d.rglob("barcodes.tsv*"))
    feat = next(d.rglob("features.tsv*"))
    mtxp = next(d.rglob("matrix.mtx*"))
    barcodes = pd.read_csv(bc, header=None, compression="infer")[0].astype(str)
    feat_df = pd.read_csv(feat, sep="\t", header=None, compression="infer")
    genes = feat_df[1].astype(str).tolist() if feat_df.shape[1] > 1 else feat_df[0].astype(str).tolist()
    mtx = io.mmread(mtxp).tocsr()
    if mtx.shape[0] == len(barcodes) and mtx.shape[1] == len(genes):
        pass
    elif mtx.shape[1] == len(barcodes) and mtx.shape[0] == len(genes):
        mtx = mtx.T.tocsr()
    else:
        raise ValueError(f"shape mismatch {mtx.shape} {len(barcodes)} {len(genes)}")
    return barcodes, genes, mtx


def read_positions(path):
    path = Path(path)
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path)
        if "barcode" not in {c.lower() for c in df.columns}:
            df = pd.read_csv(path, header=None)
            if str(df.iloc[0, 0]).lower() in {"barcode", "barcodes"}:
                df = df.iloc[1:].reset_index(drop=True)
            df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"][: df.shape[1]]
    cols = {c.lower(): c for c in df.columns}
    rename = {}
    mapping = {
        "barcode": ["barcode"],
        "in_tissue": ["in_tissue"],
        "array_row": ["array_row"],
        "array_col": ["array_col"],
        "pxl_row": ["pxl_row_in_fullres", "pxl_row"],
        "pxl_col": ["pxl_col_in_fullres", "pxl_col"],
    }
    for want, opts in mapping.items():
        for o in opts:
            if o in cols:
                rename[cols[o]] = want
                break
    df = df.rename(columns=rename)
    df["barcode"] = df["barcode"].astype(str)
    if "in_tissue" in df.columns:
        df["in_tissue"] = pd.to_numeric(df["in_tissue"], errors="coerce").fillna(0).astype(int)
    else:
        df["in_tissue"] = 1
    for c in ("array_row", "array_col", "pxl_row", "pxl_col"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("barcode")


def analyze_visium_like(name, barcodes, genes, mtx, pos):
    gidx = gene_index_from_names(genes)
    if "CLDN4" not in gidx:
        return {"section": name, "skip": "CLDN4_absent", "n_genes": len(genes)}
    if "CD8A" not in gidx:
        return {"section": name, "skip": "CD8A_absent", "genes_CLDN4": True}

    n_genes = np.asarray((mtx > 0).sum(axis=1)).ravel()
    pos = pos.reindex(barcodes)
    in_t = pos["in_tissue"].fillna(0).astype(int).to_numpy() if "in_tissue" in pos.columns else np.ones(len(barcodes), int)
    keep = (in_t == 1) & (n_genes >= MIN_GENES)
    if int(keep.sum()) < 50:
        return {"section": name, "skip": "too_few_spots", "n_keep": int(keep.sum())}

    mtx_k = mtx[keep]
    bc_k = np.asarray(barcodes)[keep]
    pos_k = pos.loc[bc_k]
    lib = np.asarray(mtx_k.sum(axis=1)).ravel()
    cl4 = np.asarray(mtx_k[:, gidx["CLDN4"]].todense()).ravel()
    cd8 = np.asarray(mtx_k[:, gidx["CD8A"]].todense()).ravel()
    krt = np.asarray(mtx_k[:, gidx["KRT8"]].todense()).ravel() if "KRT8" in gidx else np.zeros_like(cl4)
    if "array_row" in pos_k.columns and pos_k["array_row"].notna().all():
        x = pos_k["array_col"].to_numpy()
        y = pos_k["array_row"].to_numpy()
    else:
        x = pos_k.get("pxl_col", pos_k.iloc[:, 0]).to_numpy()
        y = pos_k.get("pxl_row", pos_k.iloc[:, 1]).to_numpy()
    extra = {"genes_CLDN4": True, "genes_CD8A": True, "genes_KRT8": "KRT8" in gidx, "platform": "visium_like"}
    return exclusion_bundle(name, cl4, cd8, krt, x, y, lib, extra=extra)


def bin_coo(row, col, data, nrows, ncols, scale):
    br = (row // scale).astype(np.int64)
    bc = (col // scale).astype(np.int64)
    nr = int(np.ceil(nrows / scale))
    nc = int(np.ceil(ncols / scale))
    out = np.zeros(nr * nc, dtype=np.int64)
    np.add.at(out, br * nc + bc, data)
    return out.reshape(nr, nc)


def analyze_visium_hd_feature_slice(path, name="10x_VisiumHD_Human_Lung_Cancer", scale=8):
    """16 µm bins from official 10x feature_slice (scale=8 on 2 µm grid)."""
    path = Path(path)
    with h5py.File(path, "r") as f:
        names = [x.decode() if isinstance(x, bytes) else str(x) for x in f["features"]["name"][:]]
        gidx = gene_index_from_names(names)
        if "CLDN4" not in gidx:
            return {"section": name, "skip": "CLDN4_absent"}
        meta = json.loads(f.attrs["metadata_json"])
        nrows = int(meta["nrows"])
        ncols = int(meta["ncols"])
        # tissue mask at 2 µm
        mask_g = f["masks"]["filtered"]
        mask = np.zeros((nrows, ncols), dtype=bool)
        mask[mask_g["row"][:], mask_g["col"][:]] = True
        mask_b = mask.reshape(
            int(np.ceil(nrows / scale)), scale, int(np.ceil(ncols / scale)), scale
        ).any(axis=(1, 3)) if (nrows % scale == 0 and ncols % scale == 0) else None
        if mask_b is None:
            # conservative: bin mask with OR via max
            tmp = np.zeros((int(np.ceil(nrows / scale)), int(np.ceil(ncols / scale))), dtype=np.int8)
            np.maximum.at(tmp, (mask_g["row"][:] // scale, mask_g["col"][:] // scale), 1)
            mask_b = tmp.astype(bool)

        def gene_grid(symbol):
            i = gidx[symbol]
            key = str(i)
            if key not in f["feature_slices"]:
                return np.zeros_like(mask_b, dtype=np.int64)
            g = f["feature_slices"][key]
            return bin_coo(g["row"][:], g["col"][:], g["data"][:], nrows, ncols, scale)

        umis = f["umis"]["total"]
        lib_g = bin_coo(umis["row"][:], umis["col"][:], umis["data"][:], nrows, ncols, scale)
        cl4 = gene_grid("CLDN4")
        cd8 = gene_grid("CD8A") if "CD8A" in gidx else np.zeros_like(cl4)
        krt = gene_grid("KRT8") if "KRT8" in gidx else np.zeros_like(cl4)

    keep = mask_b & (lib_g >= 50)
    if int(keep.sum()) < 50:
        return {"section": name, "skip": "too_few_bins", "n_keep": int(keep.sum())}
    rr, cc = np.where(keep)
    extra = {
        "platform": "visium_hd_16um",
        "bin_um": 2 * scale,
        "genes_CLDN4": True,
        "genes_CD8A": "CD8A" in gidx,
        "genes_KRT8": "KRT8" in gidx,
        "source": str(path.name),
    }
    if "CD8A" not in gidx:
        extra["skip"] = "CD8A_absent"
        return extra
    return exclusion_bundle(
        name,
        cl4[keep].astype(float),
        cd8[keep].astype(float),
        krt[keep].astype(float),
        cc.astype(float),
        rr.astype(float),
        lib_g[keep].astype(float),
        extra=extra,
        radius=float(NEI_RADIUS_BINS),
    )


def _radius_from_xy(xy, k=3.0):
    xy = np.asarray(xy, float)
    if xy.shape[0] < 20:
        return 2.0
    tree = cKDTree(xy)
    d, _ = tree.query(xy, k=2, workers=-1)
    nn = d[:, 1]
    nn = nn[np.isfinite(nn) & (nn > 0)]
    if nn.size == 0:
        return 2.0
    return float(k * np.median(nn))


def analyze_h5ad(path, name=None, max_obs=250000, sample_col=None):
    import anndata as ad

    path = Path(path)
    name = name or path.stem
    adata = ad.read_h5ad(path, backed=None)
    genes = [str(x) for x in adata.var_names]
    gidx = gene_index_from_names(genes)
    if "CLDN4" not in gidx:
        for col in ("gene_name", "gene_symbols", "symbol", "name"):
            if col in adata.var.columns:
                genes = [str(x) for x in adata.var[col]]
                gidx = gene_index_from_names(genes)
                break
    if "CLDN4" not in gidx:
        return [{"section": name, "skip": "CLDN4_absent", "n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)}]
    if "CD8A" not in gidx:
        return [{"section": name, "skip": "CD8A_absent", "genes_CLDN4": True}]

    xy = None
    for key in ("spatial", "X_spatial", "spatial_loc"):
        if key in adata.obsm:
            xy = np.asarray(adata.obsm[key])[:, :2]
            break
    if xy is None:
        cols = {c.lower(): c for c in adata.obs.columns}
        for xs, ys in (("x", "y"), ("coord_x", "coord_y"), ("centerx", "centery"), ("pxl_col_in_fullres", "pxl_row_in_fullres")):
            if xs in cols and ys in cols:
                xy = np.column_stack(
                    [
                        pd.to_numeric(adata.obs[cols[xs]], errors="coerce"),
                        pd.to_numeric(adata.obs[cols[ys]], errors="coerce"),
                    ]
                )
                break
    if xy is None:
        return [{"section": name, "skip": "no_coordinates", "genes_CLDN4": True}]

    X = adata.X
    if hasattr(X, "tocsc"):
        X = X.tocsc()

    def colvec(gi):
        if hasattr(X, "todense"):
            return np.asarray(X[:, gi].todense()).ravel()
        return np.asarray(X[:, gi]).ravel()

    cl4 = colvec(gidx["CLDN4"])
    cd8 = colvec(gidx["CD8A"])
    krt = colvec(gidx["KRT8"]) if "KRT8" in gidx else np.zeros_like(cl4)
    lib = np.asarray(X.sum(axis=1)).ravel() if hasattr(X, "sum") else np.ones_like(cl4)

    groups = {"all": np.ones(adata.n_obs, dtype=bool)}
    if sample_col and sample_col in adata.obs.columns:
        groups = {str(s): (adata.obs[sample_col].astype(str) == str(s)).to_numpy() for s in adata.obs[sample_col].unique()}

    rows = []
    for gname, gmask in groups.items():
        sec = name if gname == "all" else f"{name}_{gname}"
        keep = gmask & np.isfinite(xy).all(axis=1) & (lib >= 20)
        extra_note = {}
        if int(keep.sum()) > max_obs:
            rng = np.random.default_rng(0)
            idx = np.where(keep)[0]
            keep_i = rng.choice(idx, size=max_obs, replace=False)
            keep = np.zeros_like(keep)
            keep[keep_i] = True
            extra_note["subsampled_to"] = max_obs
        if int((cl4[keep] > 0).sum()) < 20:
            rows.append({"section": sec, "skip": "CLDN4_too_sparse", "n_keep": int(keep.sum()), "n_cldn4_pos": int((cl4[keep] > 0).sum())})
            continue
        xy_k = xy[keep]
        radius = _radius_from_xy(xy_k, k=3.0)
        extra = {
            "platform": "h5ad",
            "n_obs_raw": int(gmask.sum()),
            "genes_CLDN4": True,
            "genes_CD8A": True,
            "genes_KRT8": "KRT8" in gidx,
            **extra_note,
        }
        rows.append(
            exclusion_bundle(
                sec,
                cl4[keep].astype(float),
                cd8[keep].astype(float),
                krt[keep].astype(float),
                xy_k[:, 0],
                xy_k[:, 1],
                lib[keep].astype(float),
                extra=extra,
                radius=radius,
            )
        )
    return rows


def flatten_rows(rows):
    flat = []
    for r in rows:
        if not r or r.get("skip"):
            continue
        same = r.get("same_spot_CLDN4_CD8A") or {}
        resid = r.get("krt8_residual_CLDN4_CD8A") or {}
        dmw = r.get("nearest_cd8_hi_vs_lo_mw") or {}
        nmw = r.get("neighbor_cd8_hi_vs_lo_mw") or {}
        flat.append(
            {
                "section": r.get("section"),
                "n_spots": r.get("n_spots"),
                "n_cldn4_hi_epi": r.get("n_cldn4_hi_epi"),
                "n_cd8_hi": r.get("n_cd8_hi"),
                "frac_cldn4_pos": r.get("frac_cldn4_pos"),
                "frac_cd8a_pos": r.get("frac_cd8a_pos"),
                "same_rho": same.get("rho"),
                "same_p": same.get("p"),
                "krt8_resid_rho": resid.get("rho"),
                "krt8_resid_p": resid.get("p"),
                "nearest_cd8_hi_median": (r.get("nearest_cd8_from_cldn4hi") or {}).get("median"),
                "nearest_cd8_lo_median": (r.get("nearest_cd8_from_cldn4lo") or {}).get("median"),
                "nearest_cd8_mw_p": dmw.get("p"),
                "nei_cd8_hi_median": (r.get("neighbor_cd8_count_hi") or {}).get("median"),
                "nei_cd8_lo_median": (r.get("neighbor_cd8_count_lo") or {}).get("median"),
                "nei_cd8_mw_p": nmw.get("p"),
            }
        )
    return pd.DataFrame(flat)


def main_hd():
    p = DATA / "Visium_HD_Human_Lung_Cancer_feature_slice.h5"
    if not p.exists():
        return None
    print("analyze", p)
    row = analyze_visium_hd_feature_slice(p)
    if float(row.get("frac_cldn4_pos") or 0) < 0.001:
        row["skip"] = "CLDN4_too_sparse"
    dump_json(TABLES / "visiumhd_10x_lung_cancer.json", row)
    return row


def main_stereo():
    d = DATA / "gse328481"
    rows = []
    for p in sorted(d.glob("*.h5ad")):
        print("analyze stereo", p.name, flush=True)
        rows.extend(analyze_h5ad(p, name=p.stem.replace(".h5ad", "")))
    dump_json(TABLES / "gse328481_stereo_xcr.json", rows)
    df = flatten_rows(rows)
    if len(df):
        df.to_csv(TABLES / "gse328481_per_section.csv", index=False)
    return rows


def main_gse301973():
    rows = []
    for p in sorted((DATA / "gse301973").rglob("*.h5ad")):
        print("analyze visiumhd leftover", p, flush=True)
        rows.extend(analyze_h5ad(p, name=p.parent.name, sample_col="sample"))
    dump_json(TABLES / "gse301973_visiumhd.json", rows)
    df = flatten_rows(rows)
    if len(df):
        df.to_csv(TABLES / "gse301973_per_section.csv", index=False)
    return rows


def analyze_visium_folder(root, name):
    root = Path(root)
    h5s = list(root.rglob("filtered_feature_bc_matrix.h5"))
    mtxs = list(root.rglob("matrix.mtx*"))
    poss = list(root.rglob("tissue_positions*"))
    if not poss:
        return [{"section": name, "skip": "no_coordinates"}]
    # pair by parent
    rows = []
    if h5s:
        for h5 in h5s:
            pos = None
            for cand in h5.parent.rglob("tissue_positions*"):
                pos = cand
                break
            if pos is None:
                for cand in h5.parent.parent.rglob("tissue_positions*"):
                    pos = cand
                    break
            if pos is None:
                continue
            barcodes, genes, mtx = load_10x_h5(h5)
            pos_df = read_positions(pos)
            sec = f"{name}_{h5.parent.parent.name}"
            rows.append(analyze_visium_like(sec, barcodes, genes, mtx, pos_df))
    elif mtxs:
        for mtxp in mtxs:
            d = mtxp.parent
            try:
                barcodes, genes, mtx = load_mtx_dir(d)
            except Exception as e:
                rows.append({"section": str(d), "skip": f"mtx_fail:{e}"})
                continue
            pos = next(d.rglob("tissue_positions*"), None) or next(d.parent.rglob("tissue_positions*"), None)
            if pos is None:
                rows.append({"section": str(d), "skip": "no_coordinates"})
                continue
            rows.append(analyze_visium_like(f"{name}_{d.name}", barcodes, genes, mtx, read_positions(pos)))
    else:
        return [{"section": name, "skip": "no_counts"}]
    return rows


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="all")
    args = ap.parse_args()
    if args.only in {"all", "hd"}:
        print(json.dumps(main_hd(), indent=2, default=_json_default))
    if args.only in {"all", "stereo"}:
        print(json.dumps(main_stereo(), indent=2, default=_json_default)[:4000])
    if args.only in {"all", "gse301973"}:
        print(json.dumps(main_gse301973(), indent=2, default=_json_default)[:4000])
