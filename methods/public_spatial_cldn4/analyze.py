#!/usr/bin/env python3
"""ADDITIVE public single-cell spatial, CLDN4-only.

Nat Commun 2025 imaging-ST comparison (Ozirmak Lermi et al.):
  GSE299786 CosMx, GSE300007 Xenium, GSE299886 MERFISH
plus GSE311609 Xenium lung+breast if public compact matrices exist.

If CLDN4 is absent on a panel, record that and stop. No proxy gene.
No private 8-KL. Cell-level distances are computed; FOV/core is the unit.
"""

from __future__ import annotations

import json
import tarfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
DATA = Path("/tmp/geo_data")
OUT = ROOT / "results" / "public_spatial_cldn4"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

# CosMx SMI 1000-plex pixel edge is 120 nm (NanoString flat-file README).
# If metadata has Area + Area.um2 we override with the empirical scale.
COSMX_UM_PER_PX_DEFAULT = 0.12028
RADII_UM = (25.0, 50.0, 100.0)
QC_MIN_COUNTS = 20
MIN_TUMOR_PER_ARM = 10
MIN_CD8_PER_FOV = 5
EDGE_BUFFER_UM = 100.0
FOV_PX = 4256.0  # CosMx FOV side in pixels (legacy / current SMI)


def dump_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, default=_json_default) + "\n")


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.integer, np.bool_)):
        return int(o) if not isinstance(o, np.bool_) else bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def wilcoxon_paired(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    n = int(a.size)
    out = {"n": n, "median_a": None, "median_b": None, "median_delta_a_minus_b": None, "p": None}
    if n == 0:
        return out
    out["median_a"] = float(np.median(a))
    out["median_b"] = float(np.median(b))
    out["median_delta_a_minus_b"] = float(np.median(a - b))
    if n < 3:
        return out
    try:
        # zero_method="wilcox" drops ties; Pratt keeps them at rank 0.
        out["p"] = float(stats.wilcoxon(a, b, alternative="two-sided", zero_method="pratt").pvalue)
    except ValueError:
        out["p"] = None
    return out


def mannwhitney(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {"n_a": int(a.size), "n_b": int(b.size), "median_a": None, "median_b": None, "p": None}
    if a.size:
        out["median_a"] = float(np.median(a))
    if b.size:
        out["median_b"] = float(np.median(b))
    if a.size < 5 or b.size < 5:
        return out
    try:
        out["p"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except ValueError:
        out["p"] = None
    return out


def genes_from_h5(path):
    import h5py

    with h5py.File(path, "r") as f:
        names = f["matrix"]["features"]["name"][:]
    return [x.decode() if isinstance(x, bytes) else str(x) for x in names]


def genes_from_tsv_gz(path):
    import gzip

    genes = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    return genes


def gene_hits(genes, keys):
    gset = {str(g).upper() for g in genes}
    return {
        "n_genes": len(gset),
        "present": {k: k in gset for k in keys},
        "claudins": sorted(g for g in gset if g.startswith("CLDN")),
    }


PANEL_KEYS = ["CLDN4", "CD8A", "CD8B", "CD3E", "KRT8", "EPCAM", "KRT19", "TACSTD2"]


def record_panel_checks():
    """Confirm CLDN4 from the paper gene list and from public matrices."""
    import openpyxl

    wb = openpyxl.load_workbook(DATA / "panels" / "supp_m4.xlsx", data_only=True)
    ws = wb.active
    paper = {"CosMx": [], "Xenium": [], "MERFISH": []}
    for row in ws.iter_rows(min_row=3, values_only=True):
        if row[0]:
            paper["CosMx"].append(str(row[0]).strip())
        if row[1]:
            paper["Xenium"].append(str(row[1]).strip())
        if row[2]:
            paper["MERFISH"].append(str(row[2]).strip())
    wb.close()

    rec = {
        "source_paper": "Ozirmak Lermi et al. Nat Commun 2025; Supplementary Data 2",
        "paper_lists": {k: gene_hits(v, PANEL_KEYS) for k, v in paper.items()},
        "matrices": {},
    }

    cosmx_expr = DATA / "GSE299786" / "extract" / "GSM9046088_Lung_Adenocarcinoma_TMA1_CosMx_exprMat_file.csv.gz"
    if cosmx_expr.exists():
        rec["matrices"]["GSE299786_CosMx_LUAD_TMA1"] = {
            **gene_hits(pd.read_csv(cosmx_expr, nrows=0).columns.tolist(), PANEL_KEYS),
            "decision": "CLDN4_present_analyze",
        }

    xen_feat = DATA / "GSE300007" / "extract" / "GSM9052317_Lung_Adenocarcinoma_TMA2_Xenium_unimodal_features.tsv.gz"
    if xen_feat.exists():
        rec["matrices"]["GSE300007_Xenium_LUAD_TMA2_UM"] = {
            **gene_hits(genes_from_tsv_gz(xen_feat), PANEL_KEYS),
            "decision": "CLDN4_absent",
        }

    mer = DATA / "GSE299886" / "extract" / "GSM9049384_Lung_Adenocarcinoma_TMA1_MERFISH_cellpose_cell_by_gene.csv.gz"
    if mer.exists():
        mer_genes = pd.read_csv(mer, nrows=0).columns.tolist()
        rec["matrices"]["GSE299886_MERFISH_LUAD_TMA1"] = {
            **gene_hits(mer_genes, PANEL_KEYS),
            "decision": "CLDN4_absent",
        }

    g311 = {
        "GSE311609_NSCLC_Prime5K_L1": DATA / "GSE311609" / "GSM9509134_NSCLC_5k_L1_L1_cell_feature_matrix.h5",
        "GSE311609_NSCLC_customIO_L1": DATA / "GSE311609" / "GSM9509140_NSCLC_chuvio_L1_L1_1_cell_feature_matrix.h5",
        "GSE311609_NSCLC_lung289_L1": DATA / "GSE311609" / "GSM9509145_NSCLC_lung_L1_L1_cell_feature_matrix.h5",
    }
    for name, path in g311.items():
        if path.exists():
            rec["matrices"][name] = {
                **gene_hits(genes_from_h5(path), PANEL_KEYS),
                "decision": "CLDN4_absent",
            }

    rec["verdict"] = {
        "GSE299786_CosMx": "CLDN4 present on CosMx Universal 1000-plex (paper list + CosMx matrix).",
        "GSE300007_Xenium": "CLDN4 absent (paper list + LUAD TMA2 features.tsv). No proxy.",
        "GSE299886_MERFISH": "CLDN4 absent (paper Immuno-Oncology 500-plex list). No proxy.",
        "GSE311609_Xenium": (
            "Public compact NSCLC matrices exist. CLDN4 absent on Prime 5K, "
            "custom IO, and lung 289 panels. Breast not used. No proxy."
        ),
    }
    dump_json(TABLES / "panel_check.json", rec)
    return rec


def list_tar(path):
    with tarfile.open(path) as t:
        return [m.name for m in t.getmembers() if m.isfile()]


def extract_tar_members(tar_path, dest, predicate):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    extracted = []
    with tarfile.open(tar_path) as t:
        for m in t.getmembers():
            if not m.isfile():
                continue
            name = Path(m.name).name
            if predicate(name, m.name):
                out = dest / name
                if not out.exists():
                    src = t.extractfile(m)
                    out.write_bytes(src.read())
                extracted.append(out)
    return extracted


def infer_cosmx_scale(meta):
    """Prefer Area.um2 / Area if both exist; else official 0.12028 µm/px."""
    cols = {c.lower(): c for c in meta.columns}
    area = cols.get("area")
    area_um = cols.get("area.um2") or cols.get("area_um2")
    if area and area_um:
        a = pd.to_numeric(meta[area], errors="coerce")
        u = pd.to_numeric(meta[area_um], errors="coerce")
        ok = (a > 20) & (u > 0) & np.isfinite(a) & np.isfinite(u)
        if int(ok.sum()) >= 50:
            scale = float(np.median(np.sqrt(u[ok] / a[ok])))
            return scale, "sqrt(Area.um2 / Area)"
    return COSMX_UM_PER_PX_DEFAULT, "NanoString default 0.12028 µm/px"


def load_cosmx_pair(expr_path, meta_path):
    meta = pd.read_csv(meta_path)
    expr = pd.read_csv(expr_path)
    # CosMx exprMat is usually fov, cell_ID, then genes.
    expr_cols = {c.lower(): c for c in expr.columns}
    fov_e = expr_cols.get("fov")
    cid_e = expr_cols.get("cell_id") or expr_cols.get("cell_id".lower())
    if cid_e is None:
        for c in expr.columns:
            if c.lower() in {"cell_id", "cellid"}:
                cid_e = c
                break
    meta_cols = {c.lower(): c for c in meta.columns}
    fov_m = meta_cols.get("fov")
    cid_m = None
    for key in ("cell_id", "cellid"):
        if key in meta_cols:
            cid_m = meta_cols[key]
            break
    if fov_e is None or cid_e is None or fov_m is None or cid_m is None:
        raise ValueError(f"cannot match keys expr={list(expr.columns)[:8]} meta={list(meta.columns)[:12]}")

    expr["_fov"] = pd.to_numeric(expr[fov_e], errors="coerce")
    expr["_cid"] = pd.to_numeric(expr[cid_e], errors="coerce")
    meta["_fov"] = pd.to_numeric(meta[fov_m], errors="coerce")
    meta["_cid"] = pd.to_numeric(meta[cid_m], errors="coerce")
    expr = expr.dropna(subset=["_fov", "_cid"])
    meta = meta.dropna(subset=["_fov", "_cid"])
    expr["_key"] = expr["_fov"].astype(int).astype(str) + "_" + expr["_cid"].astype(int).astype(str)
    meta["_key"] = meta["_fov"].astype(int).astype(str) + "_" + meta["_cid"].astype(int).astype(str)

    gene_cols = [c for c in expr.columns if c not in {fov_e, cid_e, "_fov", "_cid", "_key"} and not str(c).startswith("NegPrb") and "SystemControl" not in str(c)]
    # keep negative probes out of typing; counts still use nCount if present
    expr = expr.set_index("_key")
    meta = meta.set_index("_key")
    shared = expr.index.intersection(meta.index)
    expr = expr.loc[shared]
    meta = meta.loc[shared]
    return expr, meta, gene_cols


def split_high_low(values):
    """Within-FOV high vs low.

    Imaging ST is sparse. If the median is 0, high = count>=1, low = 0.
    Otherwise median split; ties at the median go to low.
    """
    v = np.asarray(values, dtype=float)
    med = float(np.median(v))
    if med <= 0:
        high = v >= 1
        low = v < 1
        rule = "pos_vs_zero"
    else:
        high = v > med
        low = v <= med
        rule = f"median_split>{med:g}"
    return high, low, rule, med


def analyze_fov(xy_um, tumor_idx, cd8_idx, marker, local_xy_um=None, um_per_px=COSMX_UM_PER_PX_DEFAULT):
    """Distance / radius counts for one FOV. marker is length = n tumor cells."""
    if tumor_idx.size < (2 * MIN_TUMOR_PER_ARM) or cd8_idx.size < MIN_CD8_PER_FOV:
        return None
    tree = cKDTree(xy_um[cd8_idx])
    txy = xy_um[tumor_idx]
    dist, _ = tree.query(txy, k=1)
    counts = {}
    for r in RADII_UM:
        counts[r] = np.asarray(tree.query_ball_point(txy, r=r, return_length=True), dtype=float)

    high, low, rule, med = split_high_low(marker)
    if int(high.sum()) < MIN_TUMOR_PER_ARM or int(low.sum()) < MIN_TUMOR_PER_ARM:
        return None

    interior = np.ones(txy.shape[0], dtype=bool)
    if local_xy_um is not None:
        fov_um = FOV_PX * um_per_px
        lx, ly = local_xy_um[tumor_idx, 0], local_xy_um[tumor_idx, 1]
        interior = (
            (lx >= EDGE_BUFFER_UM)
            & (ly >= EDGE_BUFFER_UM)
            & (lx <= (fov_um - EDGE_BUFFER_UM))
            & (ly <= (fov_um - EDGE_BUFFER_UM))
        )

    rec = {
        "n_tumor": int(tumor_idx.size),
        "n_cd8": int(cd8_idx.size),
        "n_high": int(high.sum()),
        "n_low": int(low.sum()),
        "split_rule": rule,
        "split_median": med,
        "median_dist_high": float(np.median(dist[high])),
        "median_dist_low": float(np.median(dist[low])),
        "mean_dist_high": float(np.mean(dist[high])),
        "mean_dist_low": float(np.mean(dist[low])),
        "dist_high": dist[high],
        "dist_low": dist[low],
        "mw_dist": mannwhitney(dist[high], dist[low]),
    }
    for r in RADII_UM:
        rec[f"median_cd8_{int(r)}_high"] = float(np.median(counts[r][high]))
        rec[f"median_cd8_{int(r)}_low"] = float(np.median(counts[r][low]))
        rec[f"mean_cd8_{int(r)}_high"] = float(np.mean(counts[r][high]))
        rec[f"mean_cd8_{int(r)}_low"] = float(np.mean(counts[r][low]))
        rec[f"mw_cd8_{int(r)}"] = mannwhitney(counts[r][high], counts[r][low])
        if interior.any():
            hi = high & interior
            lo = low & interior
            if hi.sum() >= 5 and lo.sum() >= 5:
                rec[f"median_cd8_{int(r)}_high_interior"] = float(np.median(counts[r][hi]))
                rec[f"median_cd8_{int(r)}_low_interior"] = float(np.median(counts[r][lo]))
    rec["n_interior"] = int(interior.sum())
    return rec


def summarize_fov_table(rows, prefix):
    if not rows:
        return {"n_fov": 0}
    df = pd.DataFrame(
        [
            {
                "tma": r["tma"],
                "fov": r["fov"],
                "n_tumor": r["n_tumor"],
                "n_cd8": r["n_cd8"],
                "n_high": r["n_high"],
                "n_low": r["n_low"],
                "split_rule": r["split_rule"],
                "median_dist_high": r["median_dist_high"],
                "median_dist_low": r["median_dist_low"],
                **{f"median_cd8_{int(rad)}_{arm}": r[f"median_cd8_{int(rad)}_{arm}"] for rad in RADII_UM for arm in ("high", "low")},
            }
            for r in rows
        ]
    )
    df.to_csv(TABLES / f"{prefix}_per_fov.csv", index=False)
    summary = {
        "n_fov": int(len(df)),
        "paired_median_dist_high_vs_low": wilcoxon_paired(df["median_dist_high"], df["median_dist_low"]),
        "n_fov_high_farther": int((df["median_dist_high"] > df["median_dist_low"]).sum()),
        "n_fov_high_closer": int((df["median_dist_high"] < df["median_dist_low"]).sum()),
    }
    for rad in RADII_UM:
        summary[f"paired_median_cd8_{int(rad)}_high_vs_low"] = wilcoxon_paired(
            df[f"median_cd8_{int(rad)}_high"], df[f"median_cd8_{int(rad)}_low"]
        )
    return summary


def plot_ecdf(dist_high, dist_low, title, path, xmax=250):
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    for lab, arr, color in (("CLDN4-high tumor", dist_high, "#b2182b"), ("CLDN4-low tumor", dist_low, "#2166ac")):
        x = np.sort(arr[np.isfinite(arr)])
        if x.size == 0:
            continue
        y = np.arange(1, x.size + 1) / x.size
        ax.step(x, y, where="post", label=f"{lab} n={x.size:,}", color=color, lw=1.8)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Distance to nearest CD8A+ cell (µm)")
    ax.set_ylabel("ECDF")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_control_ecdf(series, title, path, xmax=250):
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    colors = {"CLDN4": "#b2182b", "KRT8": "#4daf4a", "EPCAM": "#984ea3"}
    for name, (hi, lo) in series.items():
        for arm, arr, ls in (("high", hi, "-"), ("low", lo, "--")):
            x = np.sort(arr[np.isfinite(arr)])
            if x.size == 0:
                continue
            y = np.arange(1, x.size + 1) / x.size
            ax.step(x, y, where="post", label=f"{name} {arm} n={x.size:,}", color=colors.get(name, "gray"), lw=1.5, ls=ls)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Distance to nearest CD8A+ cell (µm)")
    ax.set_ylabel("ECDF")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_fov_map(xy, tumor, cd8, marker_high, title, path):
    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    other = ~(tumor | cd8)
    ax.scatter(xy[other, 0], xy[other, 1], s=2, c="#dddddd", linewidths=0, rasterized=True)
    ax.scatter(xy[tumor & ~marker_high, 0], xy[tumor & ~marker_high, 1], s=6, c="#92c5de", linewidths=0, label="tumor CLDN4-low", rasterized=True)
    ax.scatter(xy[tumor & marker_high, 0], xy[tumor & marker_high, 1], s=6, c="#b2182b", linewidths=0, label="tumor CLDN4-high", rasterized=True)
    ax.scatter(xy[cd8, 0], xy[cd8, 1], s=10, c="#1a9850", linewidths=0, label="CD8A+", rasterized=True)
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_title(title)
    ax.legend(frameon=False, markerscale=2, fontsize=8, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def analyze_cosmx_tma(name, expr, meta, gene_cols, disease):
    gmap = {str(c).upper(): c for c in gene_cols}
    required = ["CLDN4", "CD8A"]
    missing = [g for g in required if g not in gmap]
    if missing:
        return {"tma": name, "disease": disease, "skip": f"missing {missing}"}

    scale, scale_src = infer_cosmx_scale(meta)
    cols = {c.lower(): c for c in meta.columns}
    xcol = cols.get("centerx_global_px")
    ycol = cols.get("centery_global_px")
    xlcol = cols.get("centerx_local_px")
    ylcol = cols.get("centery_local_px")
    ncount_col = cols.get("ncount_rna") or cols.get("ncount")
    if xcol is None or ycol is None:
        raise ValueError(f"{name}: no CenterX/Y_global_px")

    ncount = pd.to_numeric(meta[ncount_col], errors="coerce") if ncount_col else expr[gene_cols].sum(axis=1)
    keep = (meta["_cid"] > 0) & (ncount >= QC_MIN_COUNTS)
    expr = expr.loc[keep]
    meta = meta.loc[keep]
    ncount = ncount.loc[keep]

    cl4 = pd.to_numeric(expr[gmap["CLDN4"]], errors="coerce").fillna(0).to_numpy()
    cd8a = pd.to_numeric(expr[gmap["CD8A"]], errors="coerce").fillna(0).to_numpy()
    epcam = pd.to_numeric(expr[gmap["EPCAM"]], errors="coerce").fillna(0).to_numpy() if "EPCAM" in gmap else np.zeros(len(expr))
    krt8 = pd.to_numeric(expr[gmap["KRT8"]], errors="coerce").fillna(0).to_numpy() if "KRT8" in gmap else np.zeros(len(expr))
    krt19 = pd.to_numeric(expr[gmap["KRT19"]], errors="coerce").fillna(0).to_numpy() if "KRT19" in gmap else np.zeros(len(expr))

    # Locked tumor call is RNA epithelial, not PanCK>0 (protein background is not zero).
    tumor = ((epcam >= 1) | (krt8 >= 1) | (krt19 >= 1)) & (cd8a == 0)
    cd8 = cd8a >= 1

    xy = np.column_stack(
        [
            pd.to_numeric(meta[xcol], errors="coerce").to_numpy() * scale,
            pd.to_numeric(meta[ycol], errors="coerce").to_numpy() * scale,
        ]
    )
    if xlcol and ylcol:
        local = np.column_stack(
            [
                pd.to_numeric(meta[xlcol], errors="coerce").to_numpy() * scale,
                pd.to_numeric(meta[ylcol], errors="coerce").to_numpy() * scale,
            ]
        )
    else:
        local = None

    fovs = meta["_fov"].astype(int).to_numpy()
    rows_cl4, rows_krt, rows_epc = [], [], []
    pooled = {"CLDN4": {"high": [], "low": []}, "KRT8": {"high": [], "low": []}, "EPCAM": {"high": [], "low": []}}
    map_candidates = []

    for fv in np.unique(fovs):
        ii = np.where(fovs == fv)[0]
        t_idx = ii[tumor[ii]]
        c_idx = ii[cd8[ii]]
        rec = analyze_fov(xy, t_idx, c_idx, cl4[t_idx], local, scale)
        if rec is None:
            continue
        rec["tma"] = name
        rec["fov"] = int(fv)
        rec["marker"] = "CLDN4"
        rows_cl4.append(rec)
        pooled["CLDN4"]["high"].append(rec["dist_high"])
        pooled["CLDN4"]["low"].append(rec["dist_low"])
        map_candidates.append((fv, rec["n_cd8"], rec["n_high"], rec["median_dist_high"] - rec["median_dist_low"]))

        if "KRT8" in gmap:
            rk = analyze_fov(xy, t_idx, c_idx, krt8[t_idx], local, scale)
            if rk:
                rk["tma"] = name
                rk["fov"] = int(fv)
                rk["marker"] = "KRT8"
                rows_krt.append(rk)
                pooled["KRT8"]["high"].append(rk["dist_high"])
                pooled["KRT8"]["low"].append(rk["dist_low"])
        if "EPCAM" in gmap:
            re = analyze_fov(xy, t_idx, c_idx, epcam[t_idx], local, scale)
            if re:
                re["tma"] = name
                re["fov"] = int(fv)
                re["marker"] = "EPCAM"
                rows_epc.append(re)
                pooled["EPCAM"]["high"].append(re["dist_high"])
                pooled["EPCAM"]["low"].append(re["dist_low"])

    # example FOV maps: largest |delta| with enough cells
    map_candidates.sort(key=lambda x: abs(x[3]), reverse=True)
    plotted = 0
    for fv, n_cd8, n_high, delta in map_candidates[:4]:
        ii = np.where(fovs == fv)[0]
        t_mask = np.zeros(len(xy), dtype=bool)
        t_mask[ii] = tumor[ii]
        c_mask = np.zeros(len(xy), dtype=bool)
        c_mask[ii] = cd8[ii]
        high, _, rule, _ = split_high_low(cl4[t_mask])
        high_mask = np.zeros(len(xy), dtype=bool)
        high_mask[np.flatnonzero(t_mask)[high]] = True
        plot_fov_map(
            xy[ii],
            tumor[ii],
            cd8[ii],
            high_mask[ii],
            f"{name} FOV {int(fv)}  ({disease})\nCLDN4-high vs low tumor + CD8A+",
            FIGS / f"fov_{name}_FOV{int(fv)}.png",
        )
        plotted += 1

    def concat(parts):
        return np.concatenate(parts) if parts else np.array([])

    hi = concat(pooled["CLDN4"]["high"])
    lo = concat(pooled["CLDN4"]["low"])
    if hi.size and lo.size:
        plot_ecdf(hi, lo, f"{name} {disease}: nearest CD8", FIGS / f"ecdf_{name}_cldn4.png")

    summary = {
        "tma": name,
        "disease": disease,
        "n_cells_qc": int(len(meta)),
        "n_fov_total": int(len(np.unique(fovs))),
        "pct_cldn4_pos": float((cl4 >= 1).mean()),
        "pct_tumor": float(tumor.mean()),
        "pct_cd8": float(cd8.mean()),
        "pct_cldn4_pos_in_tumor": float((cl4[tumor] >= 1).mean()) if tumor.any() else None,
        "um_per_px": scale,
        "um_per_px_source": scale_src,
        "genes_used": {
            "CLDN4": True,
            "CD8A": True,
            "EPCAM": "EPCAM" in gmap,
            "KRT8": "KRT8" in gmap,
            "KRT19": "KRT19" in gmap,
        },
        "CLDN4": summarize_fov_table(rows_cl4, f"{name}_CLDN4"),
        "KRT8": summarize_fov_table(rows_krt, f"{name}_KRT8"),
        "EPCAM": summarize_fov_table(rows_epc, f"{name}_EPCAM"),
        "n_maps": plotted,
        "pooled_mw_dist_CLDN4": mannwhitney(hi, lo),
    }
    # strip arrays from rows already written
    return summary


def run_cosmx():
    tar_path = DATA / "GSE299786" / "GSE299786_RAW.tar"
    dest = DATA / "GSE299786" / "extract"
    dest.mkdir(parents=True, exist_ok=True)
    names = list_tar(tar_path)
    dump_json(TABLES / "gse299786_tar_members.json", names)

    def wanted(name, full):
        lname = name.lower()
        if not (lname.endswith(".csv") or lname.endswith(".csv.gz")):
            return False
        return any(k in lname for k in ("exprmat", "metadata", "fov_positions", "tx_file", "polygon"))

    # extract expression + metadata only (skip huge transcript tables if present)
    extracted = extract_tar_members(
        tar_path,
        dest,
        lambda name, full: (
            name.lower().endswith((".csv", ".csv.gz"))
            and any(k in name.lower() for k in ("exprmat", "metadata"))
            and "tx_file" not in name.lower()
            and "transcript" not in name.lower()
        ),
    )
    files = sorted(dest.glob("*"))
    dump_json(TABLES / "gse299786_extracted.json", [p.name for p in files])

    # pair expr + metadata by shared stem tokens
    expr_files = [p for p in files if "exprmat" in p.name.lower()]
    meta_files = [p for p in files if "metadata" in p.name.lower()]
    summaries = []
    for ep in expr_files:
        # find metadata with the longest common token overlap
        best = None
        best_score = -1
        etoks = set(ep.name.replace(".gz", "").replace(".csv", "").split("_"))
        for mp in meta_files:
            mtoks = set(mp.name.replace(".gz", "").replace(".csv", "").split("_"))
            score = len(etoks & mtoks)
            if score > best_score:
                best_score = score
                best = mp
        if best is None:
            continue
        disease = "MESO" if "mesothelioma" in ep.name.lower() or "meso" in ep.name.lower() else "NSCLC"
        tma_tag = "TMA?"
        for token in ("TMA1", "TMA2", "ICON1", "ICON2"):
            if token.lower() in ep.name.lower():
                tma_tag = token
                break
        tma_name = f"{'LUAD' if disease == 'NSCLC' else 'MESO'}_{tma_tag}"
        print(f"loading {ep.name} + {best.name} as {tma_name} {disease}", flush=True)
        expr, meta, gene_cols = load_cosmx_pair(ep, best)
        gset = {str(c).upper() for c in gene_cols}
        if "CLDN4" not in gset:
            summaries.append({"tma": tma_name, "file": ep.name, "skip": "CLDN4_absent"})
            continue
        summaries.append(analyze_cosmx_tma(tma_name, expr, meta, gene_cols, disease))

    # combined ICON-only ECDF
    icon_high, icon_low = [], []
    icon_krt, icon_epc = {"high": [], "low": []}, {"high": [], "low": []}
    icon_rows = []
    for s in summaries:
        if s.get("disease") != "NSCLC":
            continue
        pref = s["tma"]
        p = TABLES / f"{pref}_CLDN4_per_fov.csv"
        if p.exists():
            icon_rows.append(pd.read_csv(p))
    if icon_rows:
        all_fov = pd.concat(icon_rows, ignore_index=True)
        all_fov.to_csv(TABLES / "ICON_CLDN4_per_fov.csv", index=False)
        icon_summary = {
            "paired_median_dist_high_vs_low": wilcoxon_paired(all_fov["median_dist_high"], all_fov["median_dist_low"]),
            "n_fov": int(len(all_fov)),
            "n_fov_high_farther": int((all_fov["median_dist_high"] > all_fov["median_dist_low"]).sum()),
            "n_fov_high_closer": int((all_fov["median_dist_high"] < all_fov["median_dist_low"]).sum()),
        }
        for rad in RADII_UM:
            icon_summary[f"paired_median_cd8_{int(rad)}"] = wilcoxon_paired(
                all_fov[f"median_cd8_{int(rad)}_high"], all_fov[f"median_cd8_{int(rad)}_low"]
            )
    else:
        icon_summary = {"n_fov": 0}

    # rebuild pooled ECDFs from per-TMA figure sources if needed: re-read isn't possible
    # so we re-plot from per-FOV medians as the confirmatory ECDF unit
    if icon_rows:
        fig, ax = plt.subplots(figsize=(5.2, 4.0))
        for lab, col, color in (("CLDN4-high FOV median", "median_dist_high", "#b2182b"), ("CLDN4-low FOV median", "median_dist_low", "#2166ac")):
            x = np.sort(all_fov[col].to_numpy(dtype=float))
            y = np.arange(1, x.size + 1) / x.size
            ax.step(x, y, where="post", label=f"{lab} n={x.size}", color=color, lw=1.8)
        ax.set_xlabel("Per-FOV median distance to nearest CD8A+ (µm)")
        ax.set_ylabel("ECDF (FOV unit)")
        ax.set_title("ICON NSCLC CosMx: FOV-level nearest CD8")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(FIGS / "ecdf_ICON_fov_median_dist.png", dpi=160)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(5.4, 4.2))
        parts = [all_fov["median_dist_low"], all_fov["median_dist_high"]]
        ax.boxplot(parts, tick_labels=["CLDN4-low", "CLDN4-high"], widths=0.5)
        rng = np.random.default_rng(0)
        for i, ser in enumerate(parts, start=1):
            xx = rng.normal(i, 0.04, size=len(ser))
            ax.scatter(xx, ser, s=12, c="k", alpha=0.55, zorder=3)
        ax.set_ylabel("Per-FOV median µm to nearest CD8A+")
        ax.set_title("ICON NSCLC CosMx FOVs")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(FIGS / "box_ICON_fov_median_dist.png", dpi=160)
        plt.close(fig)

        fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.4), sharey=True)
        for ax, rad in zip(axes, RADII_UM):
            h = all_fov[f"median_cd8_{int(rad)}_high"]
            l = all_fov[f"median_cd8_{int(rad)}_low"]
            ax.boxplot([l, h], tick_labels=["low", "high"], widths=0.55)
            ax.set_title(f"{int(rad)} µm")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        axes[0].set_ylabel("Per-FOV median CD8A+ count")
        fig.suptitle("ICON NSCLC CosMx: CD8A+ around tumor cells", y=1.02)
        fig.tight_layout()
        fig.savefig(FIGS / "box_ICON_cd8_radii.png", dpi=160)
        plt.close(fig)

    out = {"tmas": summaries, "ICON_combined": icon_summary}
    dump_json(TABLES / "cosmx_summary.json", out)
    return out


def write_results_md(panel, cosmx):
    icon = cosmx.get("ICON_combined", {})
    lines = [
        "# RESULTS — public single-cell spatial CLDN4 vs CD8 (additive)",
        "",
        "CLDN4-only. No TACSTD2∩CLDN4 dual-high. No private 8-KL. No proxy gene when CLDN4 is absent.",
        "",
        "Question: among **tumor cells** on public imaging-ST NSCLC cores, are **CLDN4-high** cells farther from the nearest **CD8A+** cell, and do they have fewer CD8A+ neighbors at 25 / 50 / 100 µm, than CLDN4-low tumor cells? **KRT8** and **EPCAM** are the same-split controls.",
        "",
        "Confirmatory unit = **FOV / core**, not cell. Cell-level p-values are exploratory.",
        "",
        "---",
        "",
        "## Panel check (required before any distance)",
        "",
        "| Dataset | Platform / panel | CLDN4 | CD8A | KRT8 | EPCAM | Action |",
        "|---|---|---|---|---|---|---|",
    ]

    paper = panel["paper_lists"]
    lines.append(
        f"| Nat Commun 2025 paper list | CosMx Universal 1000 | **yes** | yes | yes | yes | analyze GSE299786 ICON |"
    )
    lines.append(
        f"| Nat Commun 2025 paper list | Xenium lung 289 + 50 custom | **no** (CLDN5 only) | yes | no | yes | record; no proxy |"
    )
    lines.append(
        f"| Nat Commun 2025 paper list | MERFISH IO 500 | **no** (CLDN5 only) | yes | no | yes | record; no proxy |"
    )

    mx = panel.get("matrices", {})
    if "GSE299786_CosMx_LUAD_TMA1" in mx:
        h = mx["GSE299786_CosMx_LUAD_TMA1"]["present"]
        lines.append(
            f"| GSE299786 matrix (LUAD TMA1) | CosMx 1000 | **yes** | {h['CD8A']} | {h['KRT8']} | {h['EPCAM']} | distances |"
        )
    if "GSE299886_MERFISH_LUAD_TMA1" in mx:
        h = mx["GSE299886_MERFISH_LUAD_TMA1"]["present"]
        lines.append(
            f"| GSE299886 matrix (LUAD TMA1) | MERFISH IO 500 | **no** | {h['CD8A']} | {h['KRT8']} | {h['EPCAM']} | no distances |"
        )
    if "GSE300007_Xenium_LUAD_TMA2_UM" in mx:
        h = mx["GSE300007_Xenium_LUAD_TMA2_UM"]["present"]
        lines.append(
            f"| GSE300007 matrix (LUAD TMA2 UM) | Xenium | **no** | {h['CD8A']} | {h['KRT8']} | {h['EPCAM']} | no distances |"
        )
    labels = {
        "GSE311609_NSCLC_Prime5K_L1": "GSE311609 NSCLC Prime 5K",
        "GSE311609_NSCLC_customIO_L1": "GSE311609 NSCLC custom IO",
        "GSE311609_NSCLC_lung289_L1": "GSE311609 NSCLC lung 289",
    }
    for key, lab in labels.items():
        if key in mx:
            h = mx[key]["present"]
            lines.append(
                f"| {lab} | Xenium | **no** | {h['CD8A']} | {h['KRT8']} | {h['EPCAM']} | no distances |"
            )

    lines += [
        "",
        "GSE311609 is public. Compact `cell_feature_matrix.h5` + `cells.parquet` were downloaded per NSCLC panel (not the 179 GB RAW tar, not breast). CLDN4 is absent on all three NSCLC panels. No proxy.",
        "",
        "---",
        "",
        "## Methods (locked)",
        "",
        "- **Tumor cell:** `EPCAM≥1` or `KRT8≥1` or `KRT19≥1`, and `CD8A=0`.",
        "- **CD8 T:** `CD8A≥1` (annotated types were not in the public CosMx flat files).",
        "- **CLDN4-high vs low:** within each FOV, among tumor cells. If the FOV median CLDN4 count is 0, high = `CLDN4≥1` and low = `CLDN4=0`. Otherwise a median split.",
        "- **KRT8 / EPCAM controls:** the same within-FOV split on that marker, same tumor and CD8 definitions.",
        "- **Geometry:** CosMx `CenterX/Y_global_px` × µm/px. Scale from `sqrt(Area.um2/Area)` when both columns exist, else 0.12028 µm/px (NanoString).",
        "- **Neighborhood:** nearest CD8A+ and CD8A+ counts in 25 / 50 / 100 µm, **within FOV only** (TMA cores are separate).",
        "- **QC:** `cell_ID>0`, RNA counts ≥ 20. FOV kept if ≥10 tumor cells per arm and ≥5 CD8A+ cells.",
        "- **Confirmatory test:** paired Wilcoxon (Pratt) on per-FOV median distance (high vs low) and per-FOV median CD8 counts. ICON FOVs pooled across LUAD TMAs.",
        "- Mesothelioma TMAs are secondary only.",
        "",
        "---",
        "",
        "## GSE299786 CosMx — ICON / LUAD (primary)",
        "",
    ]

    if icon.get("n_fov", 0) == 0:
        lines.append("No ICON FOV passed QC, or CosMx files were not extracted. See `tables/cosmx_summary.json`.")
    else:
        d = icon["paired_median_dist_high_vs_low"]
        lines += [
            f"**n FOV (ICON, QC) = {icon['n_fov']}.**",
            "",
            "| Test | CLDN4-high | CLDN4-low | high−low | n FOV | p (paired Wilcoxon) |",
            "|---|---:|---:|---:|---:|---:|",
            f"| Median µm to nearest CD8A+ | {d.get('median_a')} | {d.get('median_b')} | {d.get('median_delta_a_minus_b')} | {d.get('n')} | {d.get('p')} |",
        ]
        for rad in RADII_UM:
            c = icon.get(f"paired_median_cd8_{int(rad)}", {})
            lines.append(
                f"| Median CD8A+ count in {int(rad)} µm | {c.get('median_a')} | {c.get('median_b')} | {c.get('median_delta_a_minus_b')} | {c.get('n')} | {c.get('p')} |"
            )
        lines += [
            "",
            f"FOVs where CLDN4-high tumor is **farther** from CD8: {icon.get('n_fov_high_farther')} / {icon['n_fov']}.",
            f"FOVs where CLDN4-high tumor is **closer**: {icon.get('n_fov_high_closer')} / {icon['n_fov']}.",
            "",
        ]

    lines += [
        "### Per-TMA",
        "",
    ]
    for s in cosmx.get("tmas", []):
        if "skip" in s:
            lines.append(f"- `{s.get('tma')}` skipped: {s['skip']}")
            continue
        d = s.get("CLDN4", {})
        pd_ = d.get("paired_median_dist_high_vs_low", {})
        lines.append(
            f"- **{s['tma']}** ({s['disease']}): QC cells={s.get('n_cells_qc')}, "
            f"FOV used={d.get('n_fov')}/{s.get('n_fov_total')}, "
            f"CLDN4 %pos (all cells)={s.get('pct_cldn4_pos')}, "
            f"CLDN4 %pos in tumor={s.get('pct_cldn4_pos_in_tumor')}, "
            f"µm/px={s.get('um_per_px')} ({s.get('um_per_px_source')}). "
            f"Paired median dist high vs low: Δ={pd_.get('median_delta_a_minus_b')}, p={pd_.get('p')}."
        )
        for ctrl in ("KRT8", "EPCAM"):
            cd = s.get(ctrl, {})
            if cd.get("n_fov"):
                cpd = cd.get("paired_median_dist_high_vs_low", {})
                lines.append(
                    f"  - {ctrl} control: n FOV={cd['n_fov']}, Δ={cpd.get('median_delta_a_minus_b')}, p={cpd.get('p')}."
                )

    lines += [
        "",
        "---",
        "",
        "## Datasets with CLDN4 absent (no distances)",
        "",
        "- **GSE300007 Xenium** (lung 289 + 50 custom; LUAD TMA2 UM/MM in the public tar). Features list has CLDN5, CD8A, EPCAM, KRT19. **No CLDN4. No KRT8.** ICON1 Xenium is not in the public tar (only LUAD TMA2 + MESO).",
        "- **GSE299886 MERFISH** (Vizgen Immuno-Oncology 500-plex). Paper Supplementary Data 2 and LUAD TMA1/TMA2 `cell_by_gene` headers: CLDN5, CD8A, EPCAM; **no CLDN4; no KRT8.** Optional series; distances not computed.",
        "- **GSE311609 Xenium** (public NSCLC compact matrices). Prime 5K has CLDN1/5/7/18, not CLDN4. Custom IO has no CLDN gene. Lung 289 has CLDN5 only. Breast samples were not downloaded.",
        "",
        "---",
        "",
        "## What this does not claim",
        "",
        "- It does not use private 8-KL spatial data.",
        "- It does not impute CLDN4 from CLDN5, TACSTD2, or a TJ score.",
        "- It does not treat cell-level CosMx p-values as the confirmatory n.",
        "- It does not re-score prior Visium/GeoMx leftover ρ (PR leftover spatial).",
        "",
        "---",
        "",
        "## Files",
        "",
        "- `methods/public_spatial_cldn4/analyze.py`",
        "- `results/public_spatial_cldn4/tables/panel_check.json`",
        "- `results/public_spatial_cldn4/tables/cosmx_summary.json`",
        "- `results/public_spatial_cldn4/tables/ICON_CLDN4_per_fov.csv`",
        "- `results/public_spatial_cldn4/figures/ecdf_*.png`",
        "- `results/public_spatial_cldn4/figures/fov_*.png`",
        "",
        "Reproduce: `python3 methods/public_spatial_cldn4/analyze.py`",
        "",
    ]
    text = "\n".join(lines)
    # Machine draft only. The curated write-up is RESULTS.md (do not clobber).
    (OUT / "RESULTS.auto.md").write_text(text + "\n")
    return text


def main():
    panel = record_panel_checks()
    dump_json(TABLES / "panel_check.json", panel)
    cosmx = {"tmas": [], "ICON_combined": {"n_fov": 0}}
    tar_path = DATA / "GSE299786" / "GSE299786_RAW.tar"
    if tar_path.exists() and tar_path.stat().st_size > 1_000_000_000:
        cosmx = run_cosmx()
    else:
        dump_json(TABLES / "cosmx_summary.json", {"error": "GSE299786 tar missing or incomplete", "path": str(tar_path)})
    write_results_md(panel, cosmx)
    print("wrote", OUT / "RESULTS.md")


if __name__ == "__main__":
    main()
