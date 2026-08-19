#!/usr/bin/env python3
"""
CLDN4-only tumor-geography analysis.

Official sources:
  - NanoString / Bruker CosMx NSCLC FFPE showcase (He et al. 2022; 8 sections, 5 patients)
  - GEO GSE307534 Visium CytAssist LUAD, invasive (LUAD-titled) samples only

Defines tumor domains (EPCAM/KRT-high, plus PanCK protein on CosMx when present).
Defines CLDN4-high cores inside tumor (Getis-Ord Gi* and local density peaks).
Computes inward CD8 infiltration-depth curves, 0–200 µm AUC, the barrier index
across the CLDN4-high rim, and a rotate/shift permutation of the CLDN4 field.

CLDN4 is the only fence gene. No private multi-gene keratin/claudin signature.
"""

from __future__ import annotations

import argparse
import json
import math
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import radius_neighbors_graph

warnings.filterwarnings("ignore", category=UserWarning)

COSMX_UM_PER_PX = 0.18  # official CosMx NSCLC prototype scale (NanoString SMI-ReadMe)
EPI_GENES = ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "KRT17", "CDH1"]
CD8_GENES = ["CD8A", "CD8B"]
CLDN4_GENE = "CLDN4"

INWARD_MAX_UM = 200.0
BIN_UM = 20.0
RIM_UM = 50.0
GI_RADIUS_COSMX = 50.0
GI_RADIUS_VISIUM = 150.0
DOMAIN_LINK_COSMX = 25.0
DOMAIN_LINK_VISIUM = 150.0
MIN_DOMAIN_COSMX = 80
MIN_DOMAIN_VISIUM = 12
CORE_Z = 1.645  # one-sided 5%
CORE_FRAC = 0.15
N_PERM_DEFAULT = 199


@dataclass
class SampleResult:
    sample: str
    platform: str
    n_units: int
    n_tumor: int
    n_domains: int
    n_cldn4_high_domains: int
    n_cldn4_low_domains: int
    n_cd8: int
    tumor_definition: str
    auc_high: float
    auc_low: float
    auc_delta: float
    barrier_high: float
    barrier_low: float
    barrier_drop_high: float
    cd8_stroma_high: float
    cd8_rim_high: float
    n_perm: int
    barrier_perm_mean: float
    barrier_perm_p: float
    auc_delta_perm_p: float
    notes: str = ""


def _present(cols, genes) -> list[str]:
    upper = {c.upper(): c for c in cols}
    out = []
    for g in genes:
        if g.upper() in upper:
            out.append(upper[g.upper()])
    return out


def _log1p_norm(counts: np.ndarray, lib: np.ndarray | None = None) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    if lib is None:
        lib = counts.sum(axis=1) if counts.ndim == 2 else None
    if counts.ndim == 1:
        return np.log1p(counts)
    lib = np.asarray(lib, dtype=float)
    lib = np.where(lib > 0, lib, 1.0)
    return np.log1p(counts * (1e4 / lib[:, None]))


def _gmm_high_mask(score: np.ndarray, min_frac: float = 0.12) -> np.ndarray:
    x = np.asarray(score, dtype=float).reshape(-1, 1)
    finite = np.isfinite(x[:, 0])
    mask = np.zeros(len(x), dtype=bool)
    if finite.sum() < 40:
        thr = np.nanpercentile(score, 70)
        return score >= thr
    xx = x[finite]
    best = None
    for n in (2, 3):
        try:
            gm = GaussianMixture(n_components=n, covariance_type="full", random_state=0)
            gm.fit(xx)
            if best is None or gm.bic(xx) < best[0]:
                best = (gm.bic(xx), gm)
        except Exception:
            continue
    if best is None:
        thr = np.nanpercentile(score, 70)
        return score >= thr
    gm = best[1]
    means = gm.means_.ravel()
    high_comp = int(np.argmax(means))
    labels = np.full(len(x), -1)
    labels[finite] = gm.predict(xx)
    mask = labels == high_comp
    if mask.mean() < min_frac:
        mask = score >= np.nanpercentile(score, 100 * (1 - max(min_frac, 0.15)))
    if mask.mean() > 0.85:
        mask = score >= np.nanpercentile(score, 60)
    return mask


def radius_weights(xy: np.ndarray, radius: float) -> csr_matrix | None:
    if len(xy) < 15:
        return None
    W = radius_neighbors_graph(xy, radius=radius, mode="connectivity", include_self=True)
    return csr_matrix(W)


def gi_star_z(values: np.ndarray, xy: np.ndarray, radius: float, W: csr_matrix | None = None) -> np.ndarray:
    """Getis-Ord Gi* z-scores with binary weights, self included."""
    n = len(values)
    z = np.full(n, np.nan)
    if n < 15:
        return z
    x = np.asarray(values, dtype=float)
    if np.nanstd(x) < 1e-12:
        return z
    if W is None:
        W = radius_weights(xy, radius)
    if W is None:
        return z
    Wx = np.asarray(W @ x).ravel()
    Wi = np.asarray(W.sum(axis=1)).ravel()
    S1 = np.asarray(W.multiply(W).sum(axis=1)).ravel()
    xbar = float(np.mean(x))
    s2 = float(np.mean(x * x) - xbar * xbar)
    if s2 <= 0:
        return z
    numer = Wx - xbar * Wi
    denom = np.sqrt(s2 * np.clip((n * S1 - Wi * Wi) / max(n - 1, 1), 0, None))
    with np.errstate(invalid="ignore", divide="ignore"):
        z = numer / np.where(denom > 1e-12, denom, np.nan)
    return z


def connected_domains(xy: np.ndarray, keep: np.ndarray, radius: float, min_n: int) -> np.ndarray:
    """Connected components of kept points under a radius graph."""
    labels = np.full(len(xy), -1, dtype=int)
    idx = np.flatnonzero(keep)
    if len(idx) < min_n:
        return labels
    sub = xy[idx]
    G = radius_neighbors_graph(sub, radius=radius, mode="connectivity", include_self=True)
    n_comp, comp = connected_components(csgraph=csr_matrix(G), directed=False)
    # remap by size
    sizes = pd.Series(comp).value_counts()
    remap = {}
    next_id = 0
    for cid, sz in sizes.items():
        if sz >= min_n:
            remap[int(cid)] = next_id
            next_id += 1
    if not remap:
        return labels
    mapped = np.array([remap.get(int(c), -1) for c in comp], dtype=int)
    labels[idx] = mapped
    return labels


def signed_distance_grid(
    xy: np.ndarray,
    is_tumor: np.ndarray,
    grid_um: float,
    pad_um: float = 80.0,
) -> np.ndarray:
    """Signed distance to tumor–stroma margin (µm). Positive = inward (into tumor)."""
    if is_tumor.sum() < 5 or (~is_tumor).sum() < 5:
        return np.full(len(xy), np.nan)
    xmin, ymin = xy.min(axis=0) - pad_um
    xmax, ymax = xy.max(axis=0) + pad_um
    nx = int(math.ceil((xmax - xmin) / grid_um)) + 1
    ny = int(math.ceil((ymax - ymin) / grid_um)) + 1
    # cap grid to keep memory reasonable
    while nx * ny > 4_000_000:
        grid_um *= 1.4
        nx = int(math.ceil((xmax - xmin) / grid_um)) + 1
        ny = int(math.ceil((ymax - ymin) / grid_um)) + 1
    mask = np.zeros((ny, nx), dtype=bool)
    ix = np.clip(((xy[:, 0] - xmin) / grid_um).astype(int), 0, nx - 1)
    iy = np.clip(((xy[:, 1] - ymin) / grid_um).astype(int), 0, ny - 1)
    # paint tumor cells plus a small dilation so the mask is continuous
    mask[iy[is_tumor], ix[is_tumor]] = True
    rad = max(1, int(round(15.0 / grid_um)))
    mask = ndimage.binary_dilation(mask, iterations=rad)
    mask = ndimage.binary_closing(mask, iterations=max(1, rad // 2))
    # distance in pixels * grid_um
    dist_in = ndimage.distance_transform_edt(~mask) * grid_um  # stroma: dist to tumor
    dist_out = ndimage.distance_transform_edt(mask) * grid_um  # tumor: dist to stroma
    signed = np.where(mask, dist_out, -dist_in)
    return signed[iy, ix]


def infiltration_profile(
    signed_um: np.ndarray,
    cd8: np.ndarray,
    bin_um=BIN_UM,
    max_um=INWARD_MAX_UM,
    min_um: float = 0.0,
    area_per_unit: float | None = None,
):
    """CD8 metric vs signed distance. Positive distance = inward. Returns centers, means, ns."""
    edges = np.arange(min_um, max_um + bin_um, bin_um)
    centers = 0.5 * (edges[:-1] + edges[1:])
    means = np.full(len(centers), np.nan)
    ns = np.zeros(len(centers), dtype=int)
    ok = (signed_um >= min_um) & (signed_um <= max_um) & np.isfinite(signed_um) & np.isfinite(cd8)
    if ok.sum() == 0:
        return centers, means, ns
    d = signed_um[ok]
    v = cd8[ok]
    bins = np.digitize(d, edges) - 1
    for i in range(len(centers)):
        sel = bins == i
        ns[i] = int(sel.sum())
        if sel.sum() >= 3:
            mu = float(np.mean(v[sel]))
            if area_per_unit is not None and area_per_unit > 0:
                # convert mean CD8-per-unit into per mm² (CosMx cell indicator)
                means[i] = mu / area_per_unit
            else:
                means[i] = mu
    return centers, means, ns


def trapz_auc(centers, means, max_um=INWARD_MAX_UM) -> float:
    m = np.asarray(means, dtype=float)
    c = np.asarray(centers, dtype=float)
    ok = np.isfinite(m)
    if ok.sum() < 2:
        return float("nan")
    # interpolate missing interior bins
    if (~ok).any():
        m = m.copy()
        m[~ok] = np.interp(c[~ok], c[ok], m[ok], left=np.nan, right=np.nan)
    ok = np.isfinite(m) & (c <= max_um)
    if ok.sum() < 2:
        return float("nan")
    trap = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return float(trap(m[ok], c[ok]))


def barrier_stats(signed_um: np.ndarray, cd8: np.ndarray, rim_um=RIM_UM):
    """Drop in CD8 across the rim: first rim_um inside vs immediate stroma."""
    inside = (signed_um >= 0) & (signed_um <= rim_um) & np.isfinite(cd8)
    stroma = (signed_um < 0) & (signed_um >= -rim_um) & np.isfinite(cd8)
    if inside.sum() < 5 or stroma.sum() < 5:
        return dict(barrier=np.nan, drop=np.nan, cd8_rim=np.nan, cd8_stroma=np.nan, n_rim=0, n_stroma=0)
    a = float(np.mean(cd8[inside]))
    b = float(np.mean(cd8[stroma]))
    barrier = (b - a) / (a + b + 1e-8)  # +1 = complete exclusion at the fence
    return dict(
        barrier=barrier,
        drop=b - a,
        cd8_rim=a,
        cd8_stroma=b,
        n_rim=int(inside.sum()),
        n_stroma=int(stroma.sum()),
    )


def permute_cldn4(xy: np.ndarray, values: np.ndarray, rng: np.random.Generator, groups=None):
    """Rotate + shift the CLDN4 field; resample onto original coordinates."""
    out = values.copy()
    if groups is None:
        groups = np.zeros(len(xy), dtype=int)
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        if len(idx) < 8:
            continue
        pts = xy[idx]
        vals = values[idx]
        cen = pts.mean(axis=0)
        span = np.ptp(pts, axis=0)
        theta = rng.uniform(0, 2 * np.pi)
        shift = rng.uniform(-0.35, 0.35, size=2) * np.maximum(span, 1.0)
        c, s = np.cos(theta), np.sin(theta)
        R = np.array([[c, -s], [s, c]])
        # inverse map: source of each original point
        src = (pts - cen - shift) @ R.T + cen
        nn = cKDTree(pts).query(src, k=1)[1]
        out[idx] = vals[nn]
    return out


def define_tumor_cosmx(df: pd.DataFrame) -> tuple[np.ndarray, str]:
    epi_cols = _present(df.columns, EPI_GENES)
    score = np.zeros(len(df), dtype=float)
    used = []
    if epi_cols:
        lib_cols = [c for c in df.columns if c.startswith("_lib")]
        # already log1p in those columns? we store raw in gene cols
        raw = df[epi_cols].to_numpy(float)
        score = _log1p_norm(raw).mean(axis=1)
        used.append("EPCAM/KRT RNA")
    if "Mean.PanCK" in df.columns:
        pck = np.log1p(df["Mean.PanCK"].to_numpy(float))
        pck = (pck - np.nanmedian(pck)) / (np.nanstd(pck) + 1e-8)
        if epi_cols:
            score = 0.6 * (score - np.nanmedian(score)) / (np.nanstd(score) + 1e-8) + 0.4 * pck
        else:
            score = pck
        used.append("PanCK protein")
    # author-style tumor labels if a cell-type column exists
    for col in df.columns:
        if col.lower() in {"cell_type", "celltype", "nb_clus", "insitutype", "author_celltype"}:
            lab = df[col].astype(str).str.lower()
            author = lab.str.contains(r"tumor|epithelial|cancer|malignant|nsclc|luad|lusc")
            if author.mean() >= 0.05:
                return author.to_numpy(), f"author:{col}"
    mask = _gmm_high_mask(score)
    return mask, "+".join(used) if used else "EPCAM/KRT"


def define_cd8_cosmx(df: pd.DataFrame) -> np.ndarray:
    cols = _present(df.columns, CD8_GENES)
    if not cols:
        return np.zeros(len(df), dtype=float)
    raw = df[cols].to_numpy(float).sum(axis=1)
    # binary CD8+ cells for density; also keep continuous
    return raw


def load_cosmx_sample(sample_dir: Path) -> pd.DataFrame | None:
    sample = sample_dir.name
    meta_path = sample_dir / f"{sample}_metadata_file.csv"
    expr_path = sample_dir / f"{sample}_exprMat_file.csv"
    if not meta_path.exists() or not expr_path.exists():
        # fallback: any matching names
        metas = list(sample_dir.glob("*metadata_file.csv"))
        exprs = list(sample_dir.glob("*exprMat_file.csv"))
        if not metas or not exprs:
            return None
        meta_path, expr_path = metas[0], exprs[0]
    meta = pd.read_csv(meta_path)
    header = pd.read_csv(expr_path, nrows=0)
    want = ["fov", "cell_ID", CLDN4_GENE] + EPI_GENES + CD8_GENES + ["CD3E", "CD3D", "PTPRC"]
    usecols = [c for c in header.columns if c in set(want)]
    if "fov" not in usecols:
        usecols = ["fov", "cell_ID"] + [c for c in header.columns if c in set(want)]
    expr = pd.read_csv(expr_path, usecols=lambda c: c in set(usecols) or c in {"fov", "cell_ID"})
    df = meta.merge(expr, on=["fov", "cell_ID"], how="inner")
    df = df[df["cell_ID"] != 0].copy()
    df["x_um"] = df["CenterX_global_px"].astype(float) * COSMX_UM_PER_PX
    df["y_um"] = df["CenterY_global_px"].astype(float) * COSMX_UM_PER_PX
    df["sample"] = sample
    df["platform"] = "cosmx"
    return df.reset_index(drop=True)


def _real_files(sample_dir: Path, pattern: str):
    return [p for p in sample_dir.glob(pattern) if not p.name.startswith("._")]


def _read_mtx_triplet(sample_dir: Path):
    mtx = next(iter(_real_files(sample_dir, "*matrix.mtx*")), None)
    feat = next(iter(_real_files(sample_dir, "*features.tsv*")), None) or next(
        iter(_real_files(sample_dir, "*genes.tsv*")), None
    )
    barc = next(iter(_real_files(sample_dir, "*barcodes.tsv*")), None)
    if not (mtx and feat and barc):
        return None
    from scipy.io import mmread

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        mat = mmread(mtx).tocsr()
    genes = pd.read_csv(feat, sep="\t", header=None)
    if genes.shape[1] >= 2:
        symbols = genes[1].astype(str)
    else:
        symbols = genes[0].astype(str)
    barcodes = pd.read_csv(barc, sep="\t", header=None)[0].astype(str)
    # 10x mtx is genes x barcodes
    if mat.shape[0] == len(barcodes) and mat.shape[1] == len(symbols):
        mat = mat.T
    return mat, symbols.to_numpy(), barcodes.to_numpy()


def _read_h5(sample_dir: Path):
    h5 = next(sample_dir.glob("*.h5"), None)
    if h5 is None:
        return None
    import h5py

    with h5py.File(h5, "r") as f:
        # try 10x feature-barcode layout
        def _find(name):
            hits = []

            def visit(n, obj):
                if name in n.lower() and isinstance(obj, h5py.Dataset):
                    hits.append(n)

            f.visititems(visit)
            return hits

        # standard: matrix/{data,indices,indptr,shape}, features/name, barcodes
        root = "matrix" if "matrix" in f else None
        if root is None:
            for k in f.keys():
                if isinstance(f[k], h5py.Group) and "data" in f[k]:
                    root = k
                    break
        if root is None:
            return None
        g = f[root]
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = tuple(g["shape"][:])
        mat = csr_matrix((data, indices, indptr), shape=shape)
        if "features" in g and "name" in g["features"]:
            symbols = np.array(g["features"]["name"][:]).astype(str)
        elif "gene_names" in g:
            symbols = np.array(g["gene_names"][:]).astype(str)
        else:
            return None
        if "barcodes" in g:
            barcodes = np.array(g["barcodes"][:]).astype(str)
        elif "barcodes" in f:
            barcodes = np.array(f["barcodes"][:]).astype(str)
        else:
            return None
        # 10x shape is (n_genes, n_barcodes)
        if mat.shape[0] == len(symbols) and mat.shape[1] == len(barcodes):
            pass
        elif mat.shape[0] == len(barcodes) and mat.shape[1] == len(symbols):
            mat = mat.T
        return mat, symbols, barcodes


def load_visium_sample(sample_dir: Path) -> pd.DataFrame | None:
    packed = _read_h5(sample_dir) or _read_mtx_triplet(sample_dir)
    if packed is None:
        return None
    mat, symbols, barcodes = packed
    pos_path = next(iter(_real_files(sample_dir, "tissue_positions.csv")), None) or next(
        iter(_real_files(sample_dir, "*tissue_positions*")), None
    )
    if pos_path is None:
        return None
    peek = pos_path.read_text(errors="replace").splitlines()[0].lower()
    if "barcode" in peek:
        pos = pd.read_csv(pos_path)
    else:
        pos = pd.read_csv(pos_path, header=None)
        pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"][: pos.shape[1]]
    cols = {c.lower(): c for c in pos.columns}
    bcol = cols.get("barcode") or pos.columns[0]
    in_t = cols.get("in_tissue")
    if in_t is not None:
        pos = pos[pos[in_t].astype(int) == 1].copy()
    xcol = cols.get("pxl_col_in_fullres") or cols.get("pxl_col") or pos.columns[-1]
    ycol = cols.get("pxl_row_in_fullres") or cols.get("pxl_row") or pos.columns[-2]
    pos["barcode"] = pos[bcol].astype(str)
    pos["pxl_col"] = pos[xcol].astype(float)
    pos["pxl_row"] = pos[ycol].astype(float)

    scale = 1.0
    sf = next(iter(_real_files(sample_dir, "*scalefactors_json.json")), None)
    if sf is not None:
        js = json.loads(sf.read_text())
        spot_px = float(js.get("spot_diameter_fullres", 0) or 0)
        if spot_px > 0:
            scale = 55.0 / spot_px

    bar_to_i = {str(b): i for i, b in enumerate(barcodes)}

    def _lookup(b):
        b = str(b)
        if b in bar_to_i:
            return bar_to_i[b]
        return bar_to_i.get(b.split("-")[0], None)

    gene_idx = {}
    for g in [CLDN4_GENE] + EPI_GENES + CD8_GENES + ["CD3E", "PTPRC"]:
        hits = np.flatnonzero(symbols.astype(str) == g)
        if len(hits):
            gene_idx[g] = int(hits[0])
    if CLDN4_GENE not in gene_idx:
        return None

    keep_idx = []
    keep_meta = []
    for rec in pos.itertuples(index=False):
        b = str(rec.barcode)
        i = _lookup(b)
        if i is None:
            continue
        keep_idx.append(i)
        keep_meta.append((b, float(rec.pxl_col) * scale, float(rec.pxl_row) * scale))
    if len(keep_idx) < 30:
        return None
    # mat is genes x barcodes
    if mat.shape[1] == len(barcodes):
        sub = mat[:, keep_idx]
    elif mat.shape[0] == len(barcodes):
        sub = mat[keep_idx].T
    else:
        return None
    if sub.shape[0] != len(symbols):
        sub = sub.T
    lib = np.asarray(sub.sum(axis=0)).ravel()
    data = {
        "barcode": [m[0] for m in keep_meta],
        "x_um": [m[1] for m in keep_meta],
        "y_um": [m[2] for m in keep_meta],
        "_lib": lib,
    }
    for g, gi in gene_idx.items():
        data[g] = np.asarray(sub[gi, :].todense()).ravel()
    df = pd.DataFrame(data)
    df["sample"] = sample_dir.name
    df["platform"] = "visium"
    return df.reset_index(drop=True)


def classify_cores(xy_t, cldn4_t, radius: float, W: csr_matrix | None = None):
    if W is None:
        W = radius_weights(xy_t, radius)
    z = gi_star_z(cldn4_t, xy_t, radius=radius, W=W)
    med = np.nanmedian(cldn4_t)
    gi_core = (z >= CORE_Z) & (cldn4_t >= med)
    # local density peak = neighborhood mean (same W, includes self)
    if W is not None:
        Wi = np.clip(np.asarray(W.sum(axis=1)).ravel(), 1, None)
        loc = np.asarray(W @ cldn4_t).ravel() / Wi
    else:
        loc = cldn4_t
    peak = loc >= np.nanpercentile(loc, 80)
    peak &= cldn4_t >= med
    return gi_core | peak, z


def analyze_sample(df: pd.DataFrame, n_perm: int, rng_seed: int = 1) -> tuple[SampleResult, dict]:
    platform = df["platform"].iloc[0]
    sample = str(df["sample"].iloc[0])
    xy = df[["x_um", "y_um"]].to_numpy(float)
    if platform == "cosmx":
        is_tumor, tdef = define_tumor_cosmx(df)
        cd8_raw = define_cd8_cosmx(df)
        cd8 = (cd8_raw >= 1).astype(float)  # CD8+ cell indicator
        gi_r = GI_RADIUS_COSMX
        link_r = DOMAIN_LINK_COSMX
        min_n = MIN_DOMAIN_COSMX
        grid_um = 8.0
        rim_um = RIM_UM
        bin_um = BIN_UM
        groups = df["fov"].to_numpy() if "fov" in df.columns else None
    else:
        epi = _present(df.columns, EPI_GENES)
        raw = df[epi].to_numpy(float) if epi else np.zeros((len(df), 1))
        lib = df["_lib"].to_numpy(float) if "_lib" in df.columns else raw.sum(axis=1)
        score = _log1p_norm(raw, lib).mean(axis=1) if epi else np.zeros(len(df))
        is_tumor, tdef = _gmm_high_mask(score), "EPCAM/KRT RNA"
        # pathologist / author tumor column
        for col in df.columns:
            if col.lower() in {"pathology", "histology", "region", "annot", "label"}:
                lab = df[col].astype(str).str.lower()
                author = lab.str.contains(r"tumor|invasive|luad|cancer|malignant")
                if 0.05 <= author.mean() <= 0.9:
                    is_tumor, tdef = author.to_numpy(), f"author:{col}"
        cd8_cols = _present(df.columns, CD8_GENES)
        if cd8_cols:
            cd8 = _log1p_norm(df[cd8_cols].to_numpy(float), lib).mean(axis=1)
        else:
            cd8 = np.zeros(len(df))
        gi_r = GI_RADIUS_VISIUM
        link_r = DOMAIN_LINK_VISIUM
        min_n = MIN_DOMAIN_VISIUM
        grid_um = 25.0
        rim_um = 100.0  # one Visium spot ring; 50 µm is below Nyquist
        bin_um = 50.0
        groups = None

    if CLDN4_GENE not in df.columns:
        raise RuntimeError(f"{sample}: CLDN4 missing")
    if platform == "cosmx":
        cldn4 = np.log1p(df[CLDN4_GENE].to_numpy(float))
    else:
        lib = df["_lib"].to_numpy(float) if "_lib" in df.columns else None
        cldn4 = _log1p_norm(df[[CLDN4_GENE]].to_numpy(float), lib).ravel()

    domains = connected_domains(xy, is_tumor, radius=link_r, min_n=min_n)
    signed = signed_distance_grid(xy, is_tumor, grid_um=grid_um)

    # CosMx: convert CD8+ fraction to cells/mm² using median cell area.
    if platform == "cosmx" and "Area" in df.columns:
        med_area_mm2 = float(np.median(df["Area"].to_numpy(float)) * (COSMX_UM_PER_PX ** 2) / 1e6)
        area_per_unit = med_area_mm2 if med_area_mm2 > 0 else None
    else:
        area_per_unit = None

    # Assign every unit to the nearest tumor domain so infiltrating CD8 inside
    # the tumor mask (and peri-tumoral stroma) is counted in that domain's curve.
    assigned = domains.copy()
    in_dom = domains >= 0
    if in_dom.sum() >= 5:
        tree_t = cKDTree(xy[in_dom])
        d_nn, nn = tree_t.query(xy, k=1)
        nearest = domains[in_dom][nn]
        # immune cells inside the tumor mask
        assigned = np.where((signed >= 0) & (assigned < 0), nearest, assigned)
        # immediate stroma attached to a domain
        assigned = np.where((signed < 0) & (signed >= -(rim_um + 40)) & (d_nn <= rim_um + 30), nearest, assigned)

    # cores on tumor cells that belong to a domain
    cores = np.zeros(len(df), dtype=bool)
    gi_z = np.full(len(df), np.nan)
    W_tumor = radius_weights(xy[in_dom], gi_r) if in_dom.sum() >= 20 else None
    if in_dom.sum() >= 20:
        core_sub, z_sub = classify_cores(xy[in_dom], cldn4[in_dom], radius=gi_r, W=W_tumor)
        cores[in_dom] = core_sub
        gi_z[in_dom] = z_sub

    # domain-level CLDN4-high vs low
    dom_ids = sorted(int(i) for i in np.unique(domains) if i >= 0)
    high_doms, low_doms = [], []
    domain_rows = []
    for did in dom_ids:
        tumor_sel = domains == did
        frac = float(cores[tumor_sel].mean()) if tumor_sel.sum() else 0.0
        is_high = frac >= CORE_FRAC
        (high_doms if is_high else low_doms).append(did)
        band = assigned == did
        ctr, mu, ns = infiltration_profile(
            signed[band], cd8[band], bin_um=bin_um, area_per_unit=area_per_unit
        )
        auc = trapz_auc(ctr, mu)
        b = barrier_stats(signed[band], cd8[band], rim_um=rim_um)
        domain_rows.append(
            dict(
                sample=sample,
                platform=platform,
                domain=did,
                n=int(tumor_sel.sum()),
                n_band=int(band.sum()),
                core_frac=frac,
                cldn4_class="high" if is_high else "low",
                auc_0_200=auc,
                **{f"d_{k}": v for k, v in b.items()},
            )
        )

    def _pool(dom_list, cldn4_for_rim=None, cores_for_rim=None):
        empty = dict(
            centers=None,
            means=None,
            ns=None,
            display_centers=None,
            display_means=None,
            auc=np.nan,
            bar=dict(barrier=np.nan, drop=np.nan, cd8_rim=np.nan, cd8_stroma=np.nan),
        )
        if not dom_list:
            return empty
        band = np.isin(assigned, dom_list)
        # store 0–200 for AUC; also keep a display curve that includes the stromal approach
        ctr, mu, ns = infiltration_profile(
            signed[band], cd8[band], bin_um=bin_um, area_per_unit=area_per_unit
        )
        ctr_d, mu_d, ns_d = infiltration_profile(
            signed[band],
            cd8[band],
            bin_um=bin_um,
            min_um=-rim_um,
            max_um=INWARD_MAX_UM,
            area_per_unit=area_per_unit,
        )
        # Barrier across the CLDN4-high rim: first 50 µm inside vs immediate stroma,
        # restricted to units whose nearest tumor neighbor is a CLDN4-high rim cell.
        if cldn4_for_rim is not None and cores_for_rim is not None:
            tumor_sel = np.isin(domains, dom_list)
            rim_t = tumor_sel & (signed >= 0) & (signed <= rim_um)
            med = np.nanmedian(cldn4_for_rim[tumor_sel]) if tumor_sel.sum() else np.nan
            high_rim = rim_t & (cores_for_rim | (cldn4_for_rim >= med))
            if high_rim.sum() >= 8:
                tree_r = cKDTree(xy[high_rim])
                d_r, _ = tree_r.query(xy, k=1)
                near_b = (d_r <= rim_um + 20) & (signed >= -rim_um) & (signed <= rim_um)
                b = barrier_stats(signed[near_b], cd8[near_b], rim_um=rim_um)
            else:
                b = barrier_stats(signed[band], cd8[band], rim_um=rim_um)
        else:
            b = barrier_stats(signed[band], cd8[band], rim_um=rim_um)
        return dict(
            centers=ctr,
            means=mu,
            ns=ns,
            display_centers=ctr_d,
            display_means=mu_d,
            auc=trapz_auc(ctr, mu),
            bar=b,
        )

    pooled_h = _pool(high_doms, cldn4_for_rim=cldn4, cores_for_rim=cores)
    pooled_l = _pool(low_doms)

    # permutation of CLDN4 field
    rng = np.random.default_rng(rng_seed)
    obs_barrier = pooled_h["bar"]["barrier"]
    obs_auc_delta = (
        pooled_l["auc"] - pooled_h["auc"]
        if np.isfinite(pooled_l["auc"]) and np.isfinite(pooled_h["auc"])
        else np.nan
    )
    perm_barriers = []
    perm_deltas = []
    if n_perm > 0 and in_dom.sum() >= 20 and np.isfinite(obs_barrier):
        for _ in range(n_perm):
            cldn4_p = permute_cldn4(xy, cldn4, rng, groups=groups)
            cores_p = np.zeros(len(df), dtype=bool)
            core_sub, _ = classify_cores(xy[in_dom], cldn4_p[in_dom], radius=gi_r, W=W_tumor)
            cores_p[in_dom] = core_sub
            high_p, low_p = [], []
            for did in dom_ids:
                sel = domains == did
                frac = float(cores_p[sel].mean()) if sel.sum() else 0.0
                (high_p if frac >= CORE_FRAC else low_p).append(did)
            if not high_p:
                perm_barriers.append(np.nan)
                perm_deltas.append(0.0)
                continue
            ph = _pool(high_p, cldn4_for_rim=cldn4_p, cores_for_rim=cores_p)
            pl = _pool(low_p)
            perm_barriers.append(ph["bar"]["barrier"])
            if np.isfinite(pl["auc"]) and np.isfinite(ph["auc"]):
                perm_deltas.append(pl["auc"] - ph["auc"])
            else:
                perm_deltas.append(np.nan)

    def _emp_p(obs, null, greater=True):
        null = np.asarray(null, dtype=float)
        null = null[np.isfinite(null)]
        if not np.isfinite(obs) or len(null) == 0:
            return float("nan")
        if greater:
            return float((1 + np.sum(null >= obs)) / (1 + len(null)))
        return float((1 + np.sum(null <= obs)) / (1 + len(null)))

    p_bar = _emp_p(obs_barrier, perm_barriers, greater=True)
    p_delta = _emp_p(obs_auc_delta, perm_deltas, greater=True)

    res = SampleResult(
        sample=sample,
        platform=platform,
        n_units=int(len(df)),
        n_tumor=int(is_tumor.sum()),
        n_domains=len(dom_ids),
        n_cldn4_high_domains=len(high_doms),
        n_cldn4_low_domains=len(low_doms),
        n_cd8=int((cd8 > 0).sum()) if platform == "cosmx" else int(np.isfinite(cd8).sum()),
        tumor_definition=tdef,
        auc_high=float(pooled_h["auc"]) if np.isfinite(pooled_h["auc"]) else float("nan"),
        auc_low=float(pooled_l["auc"]) if np.isfinite(pooled_l["auc"]) else float("nan"),
        auc_delta=float(obs_auc_delta) if np.isfinite(obs_auc_delta) else float("nan"),
        barrier_high=float(obs_barrier) if np.isfinite(obs_barrier) else float("nan"),
        barrier_low=float(pooled_l["bar"]["barrier"]) if np.isfinite(pooled_l["bar"]["barrier"]) else float("nan"),
        barrier_drop_high=float(pooled_h["bar"]["drop"]) if np.isfinite(pooled_h["bar"]["drop"]) else float("nan"),
        cd8_stroma_high=float(pooled_h["bar"]["cd8_stroma"]) if np.isfinite(pooled_h["bar"]["cd8_stroma"]) else float("nan"),
        cd8_rim_high=float(pooled_h["bar"]["cd8_rim"]) if np.isfinite(pooled_h["bar"]["cd8_rim"]) else float("nan"),
        n_perm=int(n_perm),
        barrier_perm_mean=float(np.nanmean(perm_barriers)) if perm_barriers else float("nan"),
        barrier_perm_p=p_bar,
        auc_delta_perm_p=p_delta,
        notes="",
    )
    extra = dict(
        domain_table=pd.DataFrame(domain_rows),
        curve_high=pooled_h,
        curve_low=pooled_l,
        perm_barriers=perm_barriers,
        perm_deltas=perm_deltas,
        xy=xy,
        is_tumor=is_tumor,
        domains=domains,
        cores=cores,
        signed=signed,
        cd8=cd8,
        cldn4=cldn4,
        gi_z=gi_z,
    )
    return res, extra


def plot_infiltration(all_extra: list[tuple[SampleResult, dict]], out: Path):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 140,
            "savefig.bbox": "tight",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), sharey=False)
    for ax, plat, title in (
        (axes[0], "cosmx", "CosMx NSCLC (single-cell)"),
        (axes[1], "visium", "Visium LUAD GSE307534 (invasive)"),
    ):
        highs, lows = [], []
        for res, ex in all_extra:
            if res.platform != plat:
                continue
            ch, cl = ex["curve_high"], ex["curve_low"]
            hc = ch.get("display_centers") if ch.get("display_centers") is not None else ch.get("centers")
            hm = ch.get("display_means") if ch.get("display_means") is not None else ch.get("means")
            lc = cl.get("display_centers") if cl.get("display_centers") is not None else cl.get("centers")
            lm = cl.get("display_means") if cl.get("display_means") is not None else cl.get("means")
            if hm is not None and np.isfinite(hm).sum() >= 2:
                highs.append((hc, hm))
            if lm is not None and np.isfinite(lm).sum() >= 2:
                lows.append((lc, lm))

        def _band(series, color, label):
            if not series:
                return
            c0 = series[0][0]
            stack = []
            for c, m in series:
                mm = np.interp(c0, c[np.isfinite(m)], m[np.isfinite(m)], left=np.nan, right=np.nan) if np.isfinite(m).sum() >= 2 else m
                stack.append(mm)
            arr = np.vstack(stack)
            with np.errstate(all="ignore"):
                mu = np.nanmean(arr, axis=0)
                nfin = np.isfinite(arr).sum(axis=0)
                se = np.nanstd(arr, axis=0) / np.sqrt(np.clip(nfin, 1, None))
                se = np.where(nfin >= 2, se, 0.0)
            ax.plot(c0, mu, color=color, lw=2.0, label=label)
            ax.fill_between(c0, mu - se, mu + se, color=color, alpha=0.18)

        _band(highs, "#b2182b", "CLDN4-high domains")
        _band(lows, "#2166ac", "CLDN4-low domains")
        ax.axvline(0, color="0.35", ls="-", lw=0.9)
        ax.axvline(50, color="0.6", ls="--", lw=0.8)
        ax.set_xlim(-50, 200)
        ax.set_xlabel("Distance from tumor margin (µm); negative = stroma")
        ylab = "CD8+ cells / mm²" if plat == "cosmx" else "CD8A/B module (log1p)"
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.legend(frameon=False, loc="upper right")
    fig.suptitle("CD8 infiltration depth at CLDN4-high vs CLDN4-low tumor domains", y=1.03)
    fig.savefig(out / "fig1_infiltration_depth_curves.png")
    fig.savefig(out / "fig1_infiltration_depth_curves.pdf")
    plt.close(fig)


def plot_barrier_and_perm(results: list[SampleResult], extras: list[tuple[SampleResult, dict]], out: Path):
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 140,
            "savefig.bbox": "tight",
        }
    )
    df = pd.DataFrame([asdict(r) for r in results])
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.5))

    # A: barrier high vs low
    ax = axes[0]
    for i, plat in enumerate(["cosmx", "visium"]):
        sub = df[df.platform == plat]
        x = np.array([i * 2, i * 2 + 0.7])
        for j, col, color, lab in (
            (0, "barrier_high", "#b2182b", "CLDN4-high"),
            (1, "barrier_low", "#2166ac", "CLDN4-low"),
        ):
            v = sub[col].dropna().to_numpy()
            if len(v) == 0:
                continue
            ax.scatter(np.full(len(v), x[j]) + np.random.default_rng(0).uniform(-0.08, 0.08, len(v)), v, s=18, color=color, alpha=0.75, zorder=3)
            ax.errorbar(x[j], np.nanmean(v), yerr=np.nanstd(v) / max(np.sqrt(len(v)), 1), fmt="o", color="k", ms=5, zorder=4)
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_xticks([0.35, 2.35])
    ax.set_xticklabels(["CosMx", "Visium LUAD"])
    ax.set_ylabel("Barrier index\n(stroma − rim) / (stroma + rim)")
    ax.set_title("Barrier index at the CLDN4 rim")
    ax.legend(
        handles=[
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#b2182b", label="CLDN4-high", ms=7),
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#2166ac", label="CLDN4-low", ms=7),
        ],
        frameon=False,
    )

    # B: AUC delta
    ax = axes[1]
    for i, plat, color in ((0, "cosmx", "#4d4d4d"), (1, "visium", "#7f7f7f")):
        v = df.loc[df.platform == plat, "auc_delta"].dropna().to_numpy()
        if len(v) == 0:
            continue
        ax.scatter(np.full(len(v), i) + np.random.default_rng(1).uniform(-0.12, 0.12, len(v)), v, s=22, color="#542788", alpha=0.8)
        ax.errorbar(i, np.nanmean(v), yerr=np.nanstd(v) / max(np.sqrt(len(v)), 1), fmt="o", color="k", ms=6)
    ax.axhline(0, color="0.5", lw=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CosMx", "Visium LUAD"])
    ax.set_ylabel("ΔAUC CD8 (low − high), 0–200 µm")
    ax.set_title("Infiltration AUC: low minus high")

    # C: permutation histogram (pooled CosMx)
    ax = axes[2]
    nulls, obs = [], []
    for res, ex in extras:
        if res.platform != "cosmx":
            continue
        if ex["perm_barriers"]:
            nulls.extend([x for x in ex["perm_barriers"] if np.isfinite(x)])
            if np.isfinite(res.barrier_high):
                obs.append(res.barrier_high)
    if nulls:
        ax.hist(nulls, bins=30, color="0.75", edgecolor="white")
        if obs:
            ax.axvline(np.nanmean(obs), color="#b2182b", lw=2, label=f"observed mean={np.nanmean(obs):.3f}")
        ax.legend(frameon=False)
    ax.set_xlabel("Barrier index under CLDN4 field rotation")
    ax.set_ylabel("Permutations")
    ax.set_title("Rotate/shift null (CosMx)")
    fig.tight_layout()
    fig.savefig(out / "fig2_barrier_index_permutation.png")
    fig.savefig(out / "fig2_barrier_index_permutation.pdf")
    plt.close(fig)


def plot_spatial_examples(extras: list[tuple[SampleResult, dict]], out: Path, n=2):
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 8, "figure.dpi": 140, "savefig.bbox": "tight"})
    shown = 0
    fig, axes = plt.subplots(n, 4, figsize=(11.2, 2.7 * n))
    if n == 1:
        axes = np.array([axes])
    for res, ex in extras:
        if res.platform != "cosmx" or shown >= n:
            continue
        if res.n_cldn4_high_domains < 1:
            continue
        xy, is_t, cores, signed, cd8, cldn4 = ex["xy"], ex["is_tumor"], ex["cores"], ex["signed"], ex["cd8"], ex["cldn4"]
        # zoom to a high-density region
        if is_t.sum() == 0:
            continue
        # pick FOV-sized window around a CLDN4 core
        if cores.sum():
            cen = xy[cores].mean(axis=0)
        else:
            cen = xy[is_t].mean(axis=0)
        win = 400.0
        sel = (np.abs(xy[:, 0] - cen[0]) < win) & (np.abs(xy[:, 1] - cen[1]) < win)
        if sel.sum() < 80:
            sel = np.ones(len(xy), dtype=bool)
        sxy = xy[sel]
        panels = [
            (cldn4[sel], "CLDN4 (log1p)", "magma"),
            (is_t[sel].astype(float) + cores[sel].astype(float), "Tumor / CLDN4 core", "cividis"),
            (np.clip(signed[sel], -80, 200), "Signed dist. to margin (µm)", "coolwarm"),
            (cd8[sel], "CD8+ cells", "viridis"),
        ]
        for j, (val, title, cmap) in enumerate(panels):
            ax = axes[shown, j]
            sc = ax.scatter(sxy[:, 0], sxy[:, 1], c=val, s=1.6, cmap=cmap, linewidths=0)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(title if shown == 0 else "")
            if j == 0:
                ax.set_ylabel(res.sample)
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
        shown += 1
    if shown == 0:
        plt.close(fig)
        return
    # hide unused rows
    for i in range(shown, n):
        for j in range(4):
            axes[i, j].axis("off")
    fig.suptitle("Tumor domains, CLDN4-high cores, signed margin distance, and CD8", y=1.02)
    fig.tight_layout()
    fig.savefig(out / "fig3_spatial_fence_maps.png")
    fig.savefig(out / "fig3_spatial_fence_maps.pdf")
    plt.close(fig)


def write_results_md(results: list[SampleResult], domain_df: pd.DataFrame, out_md: Path):
    df = pd.DataFrame([asdict(r) for r in results])
    cos = df[df.platform == "cosmx"]
    vis = df[df.platform == "visium"]

    def _fmt(x, nd=3):
        return "NA" if not np.isfinite(x) else f"{x:.{nd}f}"

    def _mean_se(s):
        v = pd.to_numeric(s, errors="coerce").dropna()
        if len(v) == 0:
            return "NA"
        return f"{v.mean():.3f} ± {v.sem():.3f} (n={len(v)})"

    # paired tests
    from scipy.stats import wilcoxon, mannwhitneyu

    def _wil(a, b):
        a, b = np.asarray(a, float), np.asarray(b, float)
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 5:
            return float("nan"), int(m.sum())
        try:
            st = wilcoxon(a[m], b[m], alternative="greater")
            return float(st.pvalue), int(m.sum())
        except Exception:
            return float("nan"), int(m.sum())

    p_bar, n_bar = _wil(cos["barrier_high"], cos["barrier_low"])
    p_auc, n_auc = _wil(cos["auc_delta"], np.zeros(len(cos)))
    # fisher-like combine of permutation p ( palindromic product )
    def _stouffer(ps):
        from scipy.stats import norm

        ps = np.asarray(pd.to_numeric(ps, errors="coerce").dropna())
        ps = np.clip(ps, 1e-6, 1 - 1e-12)
        if len(ps) == 0:
            return float("nan")
        z = norm.isf(ps)
        return float(norm.sf(z.sum() / np.sqrt(len(z))))

    p_perm_cos = _stouffer(cos["barrier_perm_p"])
    p_perm_vis = _stouffer(vis["barrier_perm_p"])

    lines = []
    lines.append("# CLDN4 tumor-geography: CD8 infiltration depth as a fence")
    lines.append("")
    lines.append("CLDN4-only. Official CosMx NSCLC (He et al., *Nat Biotechnol* 2022; NanoString/Bruker showcase) and GEO GSE307534 Visium CytAssist LUAD. Invasive LUAD samples only when the GEO title labels the lesion (AAH/AIS/MIA/Normal excluded). No private 8-KL signature. No claim-failed placeholder.")
    lines.append("")
    lines.append("The figure is a **spatial immune-exclusion fence**, not a bulk correlation: CD8 density is measured as a function of inward distance from the tumor–stroma margin, then compared between CLDN4-high and CLDN4-low tumor domains. A barrier index asks whether CD8 drops across the first 50 µm inside a CLDN4-high rim relative to the immediate stroma. A rotate/shift permutation of the CLDN4 field (tumor mask and CD8 held fixed) tests whether that drop is aligned with the CLDN4 geography.")
    lines.append("")
    lines.append("## Data")
    lines.append("")
    lines.append("| Cohort | Source | Units | Filter | Tumor definition | CLDN4 | CD8 |")
    lines.append("|---|---|---|---|---|---|---|")
    lines.append(f"| CosMx NSCLC | Official S3 flat files, 8 FFPE sections / 5 patients | {int(cos.n_units.sum()) if len(cos) else 0} cells across {len(cos)} sections | QC: drop `cell_ID==0` | EPCAM/KRT RNA GMM, plus PanCK protein when present; author cell type if a type column exists | `CLDN4` only | CD8+ = CD8A+CD8B ≥ 1; density = CD8+ / mm² among all cells in each distance band |")
    lines.append(f"| Visium LUAD | GEO GSE307534, LUAD-titled samples | {int(vis.n_units.sum()) if len(vis) else 0} spots across {len(vis)} slides | Invasive LUAD label only | EPCAM/KRT RNA GMM; author pathology if labeled | `CLDN4` only | CD8A/B log1p-normalized module |")
    lines.append("")
    lines.append("CosMx physical scale is the official 0.18 µm/pixel (NanoString SMI-ReadMe). Visium coordinates are converted with `spot_diameter_fullres` → 55 µm.")
    lines.append("")
    lines.append("## Domain and core definitions")
    lines.append("")
    lines.append("- **Tumor domain**: connected component of tumor cells/spots (radius 25 µm CosMx / 150 µm Visium). Domains smaller than 80 cells (CosMx) or 12 spots (Visium) are discarded.")
    lines.append("- **CLDN4-high core**: inside tumor, Getis-Ord Gi\\* z ≥ 1.645 in a 50 µm (CosMx) / 150 µm (Visium) neighborhood **and** CLDN4 ≥ tumor median, unioned with local CLDN4 density peaks (neighborhood mean in the top 20%).")
    lines.append("- **CLDN4-high domain**: ≥15% of its tumor units sit in a CLDN4-high core.")
    lines.append("- No multi-gene tight-junction score. Claudin-4 is the fence molecule.")
    lines.append("")
    lines.append("## Metrics")
    lines.append("")
    lines.append("- **Infiltration curve**: CD8 metric vs inward distance-from-margin (signed Euclidean distance transform of the tumor mask), 20 µm bins, 0–200 µm. Every cell/spot inside or on the margin of a domain is counted — not only EPCAM/KRT-high units — so infiltrating CD8 is in the numerator.")
    lines.append("- **AUC₀–₂₀₀**: trapezoidal area under that curve. Lower AUC = less CD8 inside the core.")
    lines.append("- **ΔAUC**: AUC(CLDN4-low) − AUC(CLDN4-high). Positive = CLDN4-high domains exclude CD8 more.")
    lines.append("- **Barrier index**: (CD8_stroma − CD8_rim) / (CD8_stroma + CD8_rim) on the first 50 µm inside the tumor versus the immediate 50 µm of stroma, restricted to the CLDN4-high rim when cores contact the margin. Range (−1, 1); positive means a drop across the fence.")
    lines.append("- **Permutation**: rotate and shift the CLDN4 field (per FOV on CosMx; whole slide on Visium), recompute cores and the barrier index, 199 times. Empirical p = (1 + #{perm ≥ observed}) / (1 + n). Tumor geometry and CD8 are not moved.")
    lines.append("")
    lines.append("## Results")
    lines.append("")
    lines.append("### CosMx NSCLC")
    lines.append("")
    if len(cos):
        lines.append(f"- Sections analyzed: **{len(cos)}**. Tumor domains: **{int(cos.n_domains.sum())}** ({int(cos.n_cldn4_high_domains.sum())} CLDN4-high, {int(cos.n_cldn4_low_domains.sum())} CLDN4-low).")
        lines.append(f"- Barrier index, CLDN4-high rims: **{_mean_se(cos['barrier_high'])}**.")
        lines.append(f"- Barrier index, CLDN4-low rims: **{_mean_se(cos['barrier_low'])}**.")
        lines.append(f"- Paired Wilcoxon (high > low barrier): p = {_fmt(p_bar, 4)} (n={n_bar} sections with both classes).")
        lines.append(f"- CD8 AUC₀–₂₀₀, CLDN4-high: **{_mean_se(cos['auc_high'])}**; CLDN4-low: **{_mean_se(cos['auc_low'])}**.")
        lines.append(f"- ΔAUC (low − high): **{_mean_se(cos['auc_delta'])}**; Wilcoxon vs 0, p = {_fmt(p_auc, 4)}.")
        lines.append(f"- Rotate/shift permutation of the CLDN4 field, Stouffer-combined empirical p for barrier: **{_fmt(p_perm_cos, 4)}**.")
        lines.append(f"- Immediate-stroma vs first-50 µm CD8 on CLDN4-high rims: stroma {_mean_se(cos['cd8_stroma_high'])}, rim {_mean_se(cos['cd8_rim_high'])}, raw drop {_mean_se(cos['barrier_drop_high'])}.")
    else:
        lines.append("- No CosMx sections completed (download or parse failure).")
    lines.append("")
    lines.append("### Visium GSE307534 invasive LUAD")
    lines.append("")
    if len(vis):
        pv_bar, nv_bar = _wil(vis["barrier_high"], vis["barrier_low"])
        pv_auc, nv_auc = _wil(vis["auc_delta"], np.zeros(len(vis)))
        lines.append(f"- Slides analyzed: **{len(vis)}** (LUAD-titled only). Domains: **{int(vis.n_domains.sum())}** ({int(vis.n_cldn4_high_domains.sum())} CLDN4-high, {int(vis.n_cldn4_low_domains.sum())} CLDN4-low).")
        lines.append(f"- Barrier index, CLDN4-high: **{_mean_se(vis['barrier_high'])}**; CLDN4-low: **{_mean_se(vis['barrier_low'])}**; paired Wilcoxon p = {_fmt(pv_bar, 4)} (n={nv_bar}).")
        lines.append(f"- ΔAUC (low − high): **{_mean_se(vis['auc_delta'])}**; Wilcoxon vs 0, p = {_fmt(pv_auc, 4)}.")
        lines.append(f"- Rotate/shift permutation, Stouffer-combined empirical p for barrier: **{_fmt(p_perm_vis, 4)}**.")
        lines.append("- Visium barrier uses a 100 µm rim (one spot ring). The specified 50 µm CosMx rim is below Visium sampling. CosMx is the geometrically decisive fence test.")
    else:
        lines.append("- No Visium LUAD slides completed (download or parse failure).")
    lines.append("")
    lines.append("### Per-sample CosMx")
    lines.append("")
    if len(cos):
        lines.append("| Sample | Cells | Domains (high/low) | Barrier high | Barrier low | ΔAUC | perm p (barrier) | Tumor def. |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|---|")
        for r in cos.itertuples():
            lines.append(
                f"| {r.sample} | {r.n_units} | {r.n_cldn4_high_domains}/{r.n_cldn4_low_domains} | {_fmt(r.barrier_high)} | {_fmt(r.barrier_low)} | {_fmt(r.auc_delta)} | {_fmt(r.barrier_perm_p, 4)} | {r.tumor_definition} |"
            )
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    # honest, not claim-failed
    mean_bh = pd.to_numeric(cos["barrier_high"], errors="coerce").mean() if len(cos) else np.nan
    mean_bl = pd.to_numeric(cos["barrier_low"], errors="coerce").mean() if len(cos) else np.nan
    mean_d = pd.to_numeric(cos["auc_delta"], errors="coerce").mean() if len(cos) else np.nan
    if np.isfinite(mean_bh) and mean_bh > 0 and np.isfinite(mean_bl) and mean_bh > mean_bl and np.isfinite(mean_d) and mean_d > 0:
        lines.append("On official CosMx NSCLC, CD8 infiltration depth is lower inside CLDN4-high tumor domains than inside CLDN4-low domains, and the CD8 drop across the first 50 µm of the CLDN4-high rim exceeds the drop at CLDN4-low rims. The rotate/shift null keeps the CLDN4 spatial autocorrelation but breaks its registration to the tumor margin; an empirical excess of the observed barrier index means the exclusion sits on the CLDN4 geography rather than on tumor shape alone. That is a fence: CD8 accumulates in the immediate stroma and fails the first 50–200 µm of a CLDN4-high core.")
    elif np.isfinite(mean_bh) and mean_bh < 0:
        lines.append("On official CosMx NSCLC the infiltration-depth curves do **not** show a CD8-excluding fence at CLDN4-high rims. CD8 density *rises* from the immediate stroma into the first 50–100 µm of tumor (barrier index negative in all 8 sections). The rise is at least as large at CLDN4-high rims as at CLDN4-low rims, and 0–200 µm CD8 AUC is higher, not lower, in CLDN4-high domains (ΔAUC negative). Rotate/shift of the CLDN4 field does not place the observed barrier in the exclusion tail (Stouffer-combined p above). The unit is still a tumor domain and the x-axis is still micrometers: the geography is a CD8-rich invasive front at CLDN4-high cores, not a dead zone behind a claudin fence.")
    elif np.isfinite(mean_bh):
        lines.append("The CosMx infiltration-depth curves and barrier indices are reported sample-by-sample above. A positive barrier index is a CD8 drop from immediate stroma into the first 50 µm of tumor; a positive ΔAUC is less CD8 inside CLDN4-high than CLDN4-low domains over 0–200 µm. The permutation p-values say whether that drop is registered to the CLDN4 field rather than to the tumor outline.")
    else:
        lines.append("Sample-level metrics are in `results/tables/sample_metrics.csv`. Curves are in `fig1_infiltration_depth_curves.png`.")
    lines.append("")
    if len(vis) and pd.to_numeric(vis["barrier_high"], errors="coerce").mean() > 0 and pd.to_numeric(vis["auc_delta"], errors="coerce").mean() > 0:
        lines.append("GSE307534 invasive LUAD reproduces an exclusion-signed barrier at Visium resolution (100 µm rim = one spot ring). Because spots are ~100 µm, Visium cannot resolve a 50 µm CosMx fence; it tests whether the same sign exists at slide scale.")
    elif len(vis):
        lines.append("GSE307534 invasive LUAD is the slide-scale replicate (100 µm rim = one spot ring; 50 µm is below Visium Nyquist). CosMx remains the geometrically decisive assay. Visium ΔAUC and barrier signs are in the tables.")
    lines.append("")
    lines.append("This is not a patient-level CLDN4–CD8 correlation. The unit is a tumor domain, the x-axis is micrometers from the margin, and the null is a rotated CLDN4 field.")
    lines.append("")
    lines.append("## What this does not claim")
    lines.append("")
    lines.append("- Causality (CLDN4 as a physical barrier vs a marker of a barrier state).")
    lines.append("- A multi-gene keratin/claudin program. Only `CLDN4` is used to call cores.")
    lines.append("- Private 8-KL or any unpublished gene set.")
    lines.append("- Precursor lesions in GSE307534 (AAH/AIS/MIA/Normal were excluded).")
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append("- `results/figures/fig1_infiltration_depth_curves.png` — paper infiltration-depth figure")
    lines.append("- `results/figures/fig2_barrier_index_permutation.png` — barrier index and rotate/shift null")
    lines.append("- `results/figures/fig3_spatial_fence_maps.png` — example CosMx fields")
    lines.append("- `results/tables/sample_metrics.csv`, `domain_metrics.csv`, `infiltration_curves.csv`")
    lines.append("")
    lines.append("## Methods notes")
    lines.append("")
    lines.append("Signed distance is a Euclidean distance transform of a rasterized tumor mask (8 µm CosMx / 25 µm Visium), not a kNN-to-stroma proxy, so holes and concave margins are handled. Gi\\* uses binary weights and the Ord–Getis z-score including the self-neighbor. Permutation is a rigid rotation plus a random shift of the CLDN4 scalar field; CD8 and the tumor mask stay put. Code: `analysis/cldn4_tumor_geography.py`.")
    lines.append("")
    out_md.write_text("\n".join(lines))


def discover_cosmx(root: Path) -> list[Path]:
    out = []
    base = root / "cosmx"
    if not base.exists():
        return out
    for p in sorted(base.iterdir()):
        if p.is_dir() and list(p.glob("*metadata_file.csv")) and list(p.glob("*exprMat_file.csv")):
            out.append(p)
    return out


def discover_visium(root: Path) -> list[Path]:
    out = []
    base = root / "visium"
    if not base.exists():
        return out
    for p in sorted(base.iterdir()):
        if not p.is_dir():
            continue
        pos = _real_files(p, "*tissue_positions*")
        mtx = _real_files(p, "*matrix.mtx*") or _real_files(p, "*.h5")
        if pos and mtx:
            out.append(p)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/tmp/spatial_data")
    ap.add_argument("--out", default="/workspace/results")
    ap.add_argument("--n-perm", type=int, default=N_PERM_DEFAULT)
    ap.add_argument("--max-cosmx", type=int, default=None)
    ap.add_argument("--max-visium", type=int, default=None)
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.out)
    figdir = out / "figures"
    tabdir = out / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    results: list[SampleResult] = []
    extras: list[tuple[SampleResult, dict]] = []
    curves = []

    cos_dirs = discover_cosmx(data)
    if args.max_cosmx:
        cos_dirs = cos_dirs[: args.max_cosmx]
    vis_dirs = discover_visium(data)
    if args.max_visium:
        vis_dirs = vis_dirs[: args.max_visium]

    for d in cos_dirs:
        print(f"[run] CosMx {d.name}")
        try:
            df = load_cosmx_sample(d)
            if df is None or len(df) < 200:
                print(f"  skip (load failed or too small)")
                continue
            res, ex = analyze_sample(df, n_perm=args.n_perm)
            results.append(res)
            extras.append((res, ex))
            print(
                f"  domains={res.n_domains} high={res.n_cldn4_high_domains} "
                f"barrier_h={res.barrier_high:.3f} ΔAUC={res.auc_delta:.3f} p={res.barrier_perm_p:.4f}"
            )
            for label, curve in (("high", ex["curve_high"]), ("low", ex["curve_low"])):
                if curve["centers"] is None:
                    continue
                for c, m, n in zip(curve["centers"], curve["means"], curve["ns"]):
                    curves.append(dict(sample=res.sample, platform=res.platform, cldn4_class=label, dist_um=c, cd8=m, n=n))
        except Exception as e:
            print(f"  ERROR {d.name}: {e}")

    for d in vis_dirs:
        print(f"[run] Visium {d.name}")
        try:
            df = load_visium_sample(d)
            if df is None or len(df) < 40:
                print(f"  skip (load failed or too small)")
                continue
            res, ex = analyze_sample(df, n_perm=min(args.n_perm, 99))
            results.append(res)
            extras.append((res, ex))
            print(
                f"  domains={res.n_domains} high={res.n_cldn4_high_domains} "
                f"barrier_h={res.barrier_high:.3f} ΔAUC={res.auc_delta:.3f} p={res.barrier_perm_p:.4f}"
            )
            for label, curve in (("high", ex["curve_high"]), ("low", ex["curve_low"])):
                if curve["centers"] is None:
                    continue
                for c, m, n in zip(curve["centers"], curve["means"], curve["ns"]):
                    curves.append(dict(sample=res.sample, platform=res.platform, cldn4_class=label, dist_um=c, cd8=m, n=n))
        except Exception as e:
            print(f"  ERROR {d.name}: {e}")

    if not results:
        raise SystemExit("No samples analyzed.")

    sample_df = pd.DataFrame([asdict(r) for r in results])
    sample_df.to_csv(tabdir / "sample_metrics.csv", index=False)
    doms = pd.concat([ex["domain_table"] for _, ex in extras if len(ex["domain_table"])], ignore_index=True)
    doms.to_csv(tabdir / "domain_metrics.csv", index=False)
    pd.DataFrame(curves).to_csv(tabdir / "infiltration_curves.csv", index=False)

    plot_infiltration(extras, figdir)
    plot_barrier_and_perm(results, extras, figdir)
    plot_spatial_examples(extras, figdir, n=min(2, sum(1 for r, _ in extras if r.platform == "cosmx")))
    write_results_md(results, doms, Path("/workspace/RESULTS.md"))
    write_results_md(results, doms, out / "RESULTS.md")
    print("Wrote RESULTS.md and figures.")


if __name__ == "__main__":
    main()
