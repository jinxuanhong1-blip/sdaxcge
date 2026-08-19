#!/usr/bin/env python3
"""July-style CosMx NSCLC neighbor-count retune (CLDN4-only, official 8/5).

Primary endpoints follow the July PPT metric family:
  mean CD8/NK neighbors around CLDN4-high vs CLDN4-low tumor cells,
  section (n=8) and donor (n=5) paired Wilcoxon, n_up/n_down, mixing,
  density-normalized neighborhood fraction.

Primary (pre-specified before inspecting the grid):
  CLDN4 cut = median among tumor cells
  radius = 40 µm
  immune = CD8+NK
  tumor gate = author tumor if official/author labels exist,
               else PanCK/KRT8/EPCAM high AND not CD8
  unit = section-mean neighbor counts, then donor-mean of sections

Does not use nearest-µm to CD8 as the primary endpoint.
Does not require CD8A=0 to call a cell tumor.
Does not use private 8-KL data.
"""

from __future__ import annotations

import json
import os
import sys
from itertools import product

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "cosmx_nsclc")
OUT = os.path.join(ROOT, "results", "cosmx_july_neighbor")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")

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
    "Lung5_Rep1": "P5",
    "Lung5_Rep2": "P5",
    "Lung5_Rep3": "P5",
    "Lung6": "P6",
    "Lung9_Rep1": "P9",
    "Lung9_Rep2": "P9",
    "Lung12": "P12",
    "Lung13": "P13",
}

# CosMx NSCLC prototype (He et al. 2022): 0.18 µm / pixel
PX_TO_UM = 0.18
RADII_UM = (20, 40, 60, 80)
KNN_KS = (10, 20)
MIN_COUNTS = 20
MIN_GENES = 5
MIN_FOV_TUMOR = 50
MIN_FOV_CD8 = 20

# Pre-specified primary (declared before grid inspection)
PRIMARY = {
    "cldn4_cut": "median",
    "tumor_gate": "author_or_marker",
    "immune": "cd8_nk",
    "radius_um": 40,
    "knn_k": 10,
}

CLDN4_CUTS = ("median", "q4", "gt0", "top10")
TUMOR_GATES = ("author", "marker")
IMMUNE_SETS = ("cd8", "cd8_nk")

NEEDED = [
    "CLDN4",
    "CD8A",
    "CD8B",
    "NKG7",
    "GNLY",
    "KLRD1",
    "KRT8",
    "EPCAM",
    "KRT18",
    "KRT19",
    "CDH1",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD3G",
]
AUTHOR_TYPE_COLS = (
    "cell_type",
    "celltype",
    "CellType",
    "cellType",
    "nb_clus",
    "cell_types",
    "annotation",
    "AuthorCellType",
    "cell_ID_type",
    "ident",
)


def ensure_dirs() -> None:
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)


def find_csv(sample_dir: str, suffix: str) -> str:
    for root, _, files in os.walk(sample_dir):
        for fn in files:
            if fn.endswith(suffix):
                return os.path.join(root, fn)
    raise FileNotFoundError(f"{suffix} under {sample_dir}")


def panel_genes(sample: str) -> list[str]:
    path = find_csv(os.path.join(DATA, sample), "exprMat_file.csv")
    header = pd.read_csv(path, nrows=0)
    return [c for c in header.columns if c not in ("fov", "cell_ID")]


def stop_if_cldn4_absent(genes_by_sample: dict[str, list[str]]) -> None:
    missing = [s for s, g in genes_by_sample.items() if "CLDN4" not in g]
    if not missing:
        return
    ensure_dirs()
    msg = (
        "CLDN4 is absent from the CosMx 960-plex exprMat for: "
        + ", ".join(missing)
        + ". Stopping. Neighbor counts, mixing, and the grid were not computed."
    )
    with open(os.path.join(ROOT, "RESULTS.md"), "w") as fh:
        fh.write("# CosMx NSCLC July-style neighbor retune — STOP\n\n")
        fh.write(msg + "\n")
    print(msg, file=sys.stderr)
    raise SystemExit(2)


def load_optional_author_table() -> pd.DataFrame | None:
    """Join sidecar author labels if a public annotation table was staged."""
    candidates = [
        os.path.join(DATA, "annotations", "author_cell_types.csv"),
        os.path.join(DATA, "author_cell_types.csv"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            tab = pd.read_csv(path)
            print(f"loaded author annotations: {path} n={len(tab)}", flush=True)
            return tab
    return None


def load_sample(sample: str, author_tab: pd.DataFrame | None) -> pd.DataFrame:
    d = os.path.join(DATA, sample)
    expr_path = find_csv(d, "exprMat_file.csv")
    meta_path = find_csv(d, "metadata_file.csv")
    header = list(pd.read_csv(expr_path, nrows=0).columns)
    id_cols = ["fov", "cell_ID"]
    panel = [c for c in header if c not in id_cols and not str(c).lower().startswith("neg")]
    genes = [g for g in NEEDED if g in header]
    usecols = id_cols + genes
    chunks = []
    for chunk in pd.read_csv(expr_path, usecols=lambda c: c in set(usecols) or c in panel, chunksize=50_000):
        chunk = chunk[chunk["cell_ID"] != 0]
        if chunk.empty:
            continue
        present_panel = [g for g in panel if g in chunk.columns]
        chunk["_tot"] = chunk[present_panel].sum(axis=1, numeric_only=True)
        chunk["_nfeat"] = (chunk[present_panel] > 0).sum(axis=1)
        chunks.append(chunk[id_cols + genes + ["_tot", "_nfeat"]])
    expr = pd.concat(chunks, ignore_index=True)
    expr["cell_key"] = expr["fov"].astype(str) + "_" + expr["cell_ID"].astype(str)
    meta = pd.read_csv(meta_path)
    meta["cell_key"] = meta["fov"].astype(str) + "_" + meta["cell_ID"].astype(str)
    keep_meta = ["cell_key", "fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"]
    for col in ("nCount_RNA", "nFeature_RNA", "Mean.PanCK", "Mean.CD45", "Mean.CD3", "Area"):
        if col in meta.columns:
            keep_meta.append(col)
    author_col = next((c for c in AUTHOR_TYPE_COLS if c in meta.columns), None)
    if author_col:
        keep_meta.append(author_col)
    meta = meta[keep_meta].copy()
    df = expr.merge(meta, on=["cell_key", "fov", "cell_ID"], how="inner")
    if author_tab is not None:
        sub = author_tab.copy()
        if "sample" in sub.columns:
            sub = sub[sub["sample"].astype(str) == sample]
        key_cols = [c for c in ("cell_key", "fov", "cell_ID") if c in sub.columns]
        type_col = next((c for c in AUTHOR_TYPE_COLS if c in sub.columns), None)
        if type_col and {"fov", "cell_ID"}.issubset(sub.columns):
            sub = sub[["fov", "cell_ID", type_col]].drop_duplicates(["fov", "cell_ID"])
            df = df.merge(sub, on=["fov", "cell_ID"], how="left", suffixes=("", "_ann"))
            if type_col not in df.columns and f"{type_col}_ann" in df.columns:
                df[type_col] = df[f"{type_col}_ann"]
            author_col = type_col
    gene_cols = [g for g in genes if g in df.columns]
    raw = df[gene_cols].to_numpy(dtype=np.float32)
    if "nCount_RNA" in df.columns:
        tot_full = df["nCount_RNA"].to_numpy(dtype=np.float32)
        nfeat = df["nFeature_RNA"].to_numpy(dtype=np.float32) if "nFeature_RNA" in df.columns else df["_nfeat"].to_numpy()
    else:
        tot_full = df["_tot"].to_numpy(dtype=np.float32)
        nfeat = df["_nfeat"].to_numpy()
    keep = (tot_full >= MIN_COUNTS) & (nfeat >= MIN_GENES)
    df = df.loc[keep].reset_index(drop=True)
    raw = raw[keep]
    tot = np.maximum(raw.sum(axis=1), 1.0)
    sf = float(np.median(tot)) / tot
    logn = np.log1p(raw * sf[:, None])
    for i, g in enumerate(gene_cols):
        df[f"{g}_raw"] = raw[:, i]
        df[f"{g}_ln"] = logn[:, i]
    df["sample"] = sample
    df["patient"] = PATIENT[sample]
    df["author_type"] = df[author_col].astype(str) if author_col and author_col in df.columns else ""
    df["author_col"] = author_col or ""
    df["_n_panel_genes"] = len(panel)
    print(
        f"{sample}: cells={len(df):,} genes_used={len(gene_cols)} author_col={author_col or 'none'}",
        flush=True,
    )
    return df


def author_match(series: pd.Series, needles: tuple[str, ...]) -> np.ndarray:
    if series is None or series.empty or (series.astype(str) == "").all():
        return np.zeros(len(series) if series is not None else 0, dtype=bool)
    s = series.astype(str).str.lower()
    out = np.zeros(len(series), dtype=bool)
    for n in needles:
        out |= s.str.contains(n, regex=False, na=False).to_numpy()
    return out


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    mu = np.nanmean(x)
    sd = np.nanstd(x)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(x)
    return (x - mu) / sd


def assign_base_labels(df: pd.DataFrame) -> pd.DataFrame:
    author_tumor = author_match(df["author_type"], ("tumor", "malignant", "cancer"))
    author_epi = author_match(df["author_type"], ("epithelial",))
    author_cd8 = author_match(df["author_type"], ("t cd8", "t-cd8", "cd8"))
    author_t = author_match(df["author_type"], ("t cd", "t-cell", "t cell", "treg", "t cd4", "t-cd4"))
    author_nk = author_match(df["author_type"], ("nk",))
    df["used_author_labels"] = bool(df["author_col"].iloc[0]) and bool(df["author_type"].astype(str).ne("").any())

    cd8a = df["CD8A_raw"].to_numpy() > 0 if "CD8A_raw" in df.columns else np.zeros(len(df), dtype=bool)
    nkg7 = df["NKG7_raw"].to_numpy() > 0 if "NKG7_raw" in df.columns else np.zeros(len(df), dtype=bool)

    # CD8A+ / author T. Do not steal author-tumor cells into the CD8 set.
    is_cd8 = author_cd8 | (author_t & cd8a) | (cd8a & ~author_tumor)
    is_nk = author_nk | (nkg7 & ~is_cd8 & ~author_tumor)

    epi_parts = []
    if "Mean.PanCK" in df.columns:
        epi_parts.append(zscore(np.log1p(df["Mean.PanCK"].to_numpy(dtype=np.float64))))
    if "KRT8_ln" in df.columns:
        epi_parts.append(zscore(df["KRT8_ln"].to_numpy(dtype=np.float64)))
    if "EPCAM_ln" in df.columns:
        epi_parts.append(zscore(df["EPCAM_ln"].to_numpy(dtype=np.float64)))
    if epi_parts:
        epi = np.vstack(epi_parts).mean(axis=0)
    else:
        epi = np.zeros(len(df), dtype=np.float64)
    epi_high = epi >= np.nanmedian(epi)

    # Marker tumor: PanCK/KRT8/EPCAM high AND not CD8. No extra CD8A==0 on author tumor.
    is_tumor_marker = epi_high & (~is_cd8)
    is_tumor_author = author_tumor | author_epi
    is_tumor_author_or_marker = is_tumor_author | is_tumor_marker

    df["is_cd8"] = is_cd8
    df["is_nk"] = is_nk
    df["is_cd8_nk"] = is_cd8 | is_nk
    df["is_tumor_author"] = is_tumor_author
    df["is_tumor_marker"] = is_tumor_marker
    df["is_tumor_author_or_marker"] = is_tumor_author_or_marker
    df["epi_score"] = epi
    df["cldn4_ln"] = df["CLDN4_ln"] if "CLDN4_ln" in df.columns else 0.0
    df["cldn4_raw"] = df["CLDN4_raw"] if "CLDN4_raw" in df.columns else 0.0
    df["krt8_ln"] = df["KRT8_ln"] if "KRT8_ln" in df.columns else 0.0
    return df


def tumor_mask(df: pd.DataFrame, gate: str) -> np.ndarray:
    if gate == "author":
        return df["is_tumor_author"].to_numpy()
    if gate == "marker":
        return df["is_tumor_marker"].to_numpy()
    return df["is_tumor_author_or_marker"].to_numpy()


def immune_mask(df: pd.DataFrame, immune: str) -> np.ndarray:
    if immune == "cd8":
        return df["is_cd8"].to_numpy()
    return df["is_cd8_nk"].to_numpy()


def cldn4_high_low(vals: np.ndarray, cut: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (is_high, is_low) among the provided tumor CLDN4 values."""
    vals = np.asarray(vals, dtype=np.float64)
    n = len(vals)
    high = np.zeros(n, dtype=bool)
    low = np.zeros(n, dtype=bool)
    if n == 0:
        return high, low
    if cut == "median":
        thr = np.nanmedian(vals)
        high = vals > thr
        low = vals <= thr
    elif cut == "q4":
        q1, q3 = np.nanquantile(vals, [0.25, 0.75])
        high = vals >= q3
        low = vals <= q1
    elif cut == "gt0":
        raw_proxy = vals  # log-normed; >0 still splits expressed vs not after log1p of scaled
        high = raw_proxy > 0
        low = raw_proxy <= 0
    elif cut == "top10":
        thr = np.nanquantile(vals, 0.90)
        high = vals >= thr
        # complement among tumor as low, excluding the top decile
        low = vals < thr
        # tighter low: bottom 50% so the contrast is not high vs everyone
        low = vals <= np.nanmedian(vals)
    else:
        raise ValueError(cut)
    return high, low


def xy_um(df: pd.DataFrame) -> np.ndarray:
    return np.column_stack(
        [
            df["CenterX_global_px"].to_numpy(dtype=np.float64) * PX_TO_UM,
            df["CenterY_global_px"].to_numpy(dtype=np.float64) * PX_TO_UM,
        ]
    )


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    out = np.full_like(y, np.nan, dtype=np.float64)
    ok = np.isfinite(y) & np.isfinite(x)
    if ok.sum() < 20:
        return out
    model = LinearRegression()
    model.fit(x[ok].reshape(-1, 1), y[ok])
    out[ok] = y[ok] - model.predict(x[ok].reshape(-1, 1))
    return out


def mixing_homogeneous(xy_a: np.ndarray, xy_b: np.ndarray, radius: float) -> float:
    """n_AB / (n_AB + n_AA + n_BB) among A∪B within radius (self excluded)."""
    if len(xy_a) < 5 or len(xy_b) < 5:
        return float("nan")
    xy = np.vstack([xy_a, xy_b])
    lab = np.array([0] * len(xy_a) + [1] * len(xy_b))
    tree = cKDTree(xy)
    neigh = tree.query_ball_point(xy, r=float(radius))
    n_aa = n_bb = n_ab = 0
    for i, nbrs in enumerate(neigh):
        a = lab[i]
        for j in nbrs:
            if j <= i:
                continue
            b = lab[j]
            if a == 0 and b == 0:
                n_aa += 1
            elif a == 1 and b == 1:
                n_bb += 1
            else:
                n_ab += 1
    denom = n_ab + n_aa + n_bb
    return float(n_ab / denom) if denom > 0 else float("nan")


def paired_wilcoxon(high: np.ndarray, low: np.ndarray) -> dict:
    high = np.asarray(high, dtype=np.float64)
    low = np.asarray(low, dtype=np.float64)
    ok = np.isfinite(high) & np.isfinite(low)
    high, low = high[ok], low[ok]
    n = int(len(high))
    if n == 0:
        return {
            "n": 0,
            "n_up": 0,
            "n_down": 0,
            "n_tie": 0,
            "mean_high": float("nan"),
            "mean_low": float("nan"),
            "delta": float("nan"),
            "p": float("nan"),
        }
    n_up = int((high > low).sum())
    n_down = int((high < low).sum())
    n_tie = int((high == low).sum())
    rec = {
        "n": n,
        "n_up": n_up,
        "n_down": n_down,
        "n_tie": n_tie,
        "mean_high": float(np.mean(high)),
        "mean_low": float(np.mean(low)),
        "delta": float(np.mean(high) - np.mean(low)),
        "p": float("nan"),
    }
    d = high - low
    if n >= 3 and np.any(d != 0):
        try:
            rec["p"] = float(stats.wilcoxon(high, low, alternative="two-sided", zero_method="wilcox").pvalue)
        except ValueError:
            rec["p"] = float("nan")
    return rec


def fov_keep_mask(df: pd.DataFrame, is_tumor: np.ndarray, is_cd8: np.ndarray) -> np.ndarray:
    keep = np.zeros(len(df), dtype=bool)
    for (sample, fov), idx in df.groupby(["sample", "fov"], sort=False).groups.items():
        idx = np.asarray(list(idx))
        n_t = int(is_tumor[idx].sum())
        n_c = int(is_cd8[idx].sum())
        if n_t >= MIN_FOV_TUMOR and n_c >= MIN_FOV_CD8:
            keep[idx] = True
    return keep


def compute_cell_neighbors(
    df: pd.DataFrame,
    is_tumor: np.ndarray,
    is_immune: np.ndarray,
    is_epi_local: np.ndarray,
    radii: tuple[int, ...] = RADII_UM,
    knn_ks: tuple[int, ...] = KNN_KS,
) -> pd.DataFrame:
    """Per tumor cell: radius counts, knn counts, neighborhood fraction, local epi density."""
    rows = []
    for sample, g_pos in df.groupby("sample", sort=False).groups.items():
        pos = np.asarray(list(g_pos))
        xy = xy_um(df.iloc[pos])
        tumor_loc = np.flatnonzero(is_tumor[pos])
        imm_loc = np.flatnonzero(is_immune[pos])
        epi_loc = np.flatnonzero(is_epi_local[pos])
        if len(tumor_loc) == 0:
            continue
        t_xy = xy[tumor_loc]
        tree_all = cKDTree(xy)
        tree_imm = cKDTree(xy[imm_loc]) if len(imm_loc) else None
        tree_epi = cKDTree(xy[epi_loc]) if len(epi_loc) else None

        rec = {
            "idx": pos[tumor_loc],
            "sample": np.full(len(tumor_loc), sample),
        }
        # all-neighbor counts for fraction (exclude self: query_ball includes self)
        for r in radii:
            n_all = tree_all.query_ball_point(t_xy, r=float(r), return_length=True).astype(np.float64) - 1.0
            n_all = np.clip(n_all, 0, None)
            if tree_imm is None:
                n_imm = np.zeros(len(t_xy), dtype=np.float64)
            else:
                n_imm = tree_imm.query_ball_point(t_xy, r=float(r), return_length=True).astype(np.float64)
            if tree_epi is None:
                n_epi = np.zeros(len(t_xy), dtype=np.float64)
            else:
                n_epi = tree_epi.query_ball_point(t_xy, r=float(r), return_length=True).astype(np.float64)
            # a tumor cell that is also epi would count itself in n_epi
            rec[f"imm_{r}"] = n_imm
            rec[f"all_{r}"] = n_all
            rec[f"epi_{r}"] = n_epi
            rec[f"frac_{r}"] = np.divide(n_imm, n_all, out=np.full_like(n_imm, np.nan), where=n_all > 0)
            # local KRT8 intensity proxy = the tumor cell's own KRT8 (fast; residual
            # primary covariate is epithelial-neighbor count in the same radius)
            rec[f"krt8loc_{r}"] = df["krt8_ln"].to_numpy()[pos[tumor_loc]]

        if len(xy) >= 2:
            k_all = min(max(knn_ks) + 1, len(xy))
            _dist, ind = tree_all.query(t_xy, k=k_all, workers=-1)
            if k_all == 1:
                ind = ind.reshape(-1, 1)
            if ind.shape[1] >= 2:
                ind = ind[:, 1:]
            imm_flag = np.zeros(len(xy), dtype=bool)
            if len(imm_loc):
                imm_flag[imm_loc] = True
            for k in knn_ks:
                kk = min(k, ind.shape[1])
                cnt = imm_flag[ind[:, :kk]].sum(axis=1).astype(np.float64)
                rec[f"knn_{k}"] = cnt
                rec[f"knnfrac_{k}"] = cnt / float(kk)
        else:
            for k in knn_ks:
                rec[f"knn_{k}"] = np.zeros(len(t_xy), dtype=np.float64)
                rec[f"knnfrac_{k}"] = np.full(len(t_xy), np.nan)

        part = pd.DataFrame(rec)
        rows.append(part)

    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    return out


def summarize_setting(
    df: pd.DataFrame,
    neigh: pd.DataFrame,
    is_tumor: np.ndarray,
    is_high: np.ndarray,
    is_low: np.ndarray,
    radius: int,
    knn_k: int,
    immune_xy_by_sample: dict[str, np.ndarray],
    high_xy_by_sample: dict[str, np.ndarray],
    low_xy_by_sample: dict[str, np.ndarray],
    all_tum_xy_by_sample: dict[str, np.ndarray],
    cd8nk_xy_by_sample: dict[str, np.ndarray],
) -> dict:
    cell = neigh.copy()
    cell["is_high"] = is_high[cell["idx"].to_numpy()]
    cell["is_low"] = is_low[cell["idx"].to_numpy()]
    cell["patient"] = df["patient"].to_numpy()[cell["idx"].to_numpy()]
    cell["fov"] = df["fov"].to_numpy()[cell["idx"].to_numpy()]

    imm_col = f"imm_{radius}"
    frac_col = f"frac_{radius}"
    epi_col = f"epi_{radius}"
    krt_col = f"krt8loc_{radius}"
    knn_col = f"knn_{knn_k}"

    # residualize within section
    cell["imm_resid"] = np.nan
    cell["frac_resid"] = np.nan
    for sample, g in cell.groupby("sample"):
        ix = g.index.to_numpy()
        y = g[imm_col].to_numpy()
        x_epi = g[epi_col].to_numpy()
        x_krt = g[krt_col].to_numpy()
        # prefer local epithelial density; fall back to KRT8 if epi is degenerate
        resid = residualize(y, x_epi)
        if not np.isfinite(resid).any():
            resid = residualize(y, x_krt)
        cell.loc[ix, "imm_resid"] = resid
        cell.loc[ix, "frac_resid"] = residualize(g[frac_col].to_numpy(), x_epi)

    section_rows = []
    for sample, g in cell.groupby("sample"):
        hi = g[g["is_high"]]
        lo = g[g["is_low"]]
        section_rows.append(
            {
                "sample": sample,
                "patient": PATIENT[sample],
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "mean_imm_high": float(hi[imm_col].mean()) if len(hi) else np.nan,
                "mean_imm_low": float(lo[imm_col].mean()) if len(lo) else np.nan,
                "mean_frac_high": float(hi[frac_col].mean()) if len(hi) else np.nan,
                "mean_frac_low": float(lo[frac_col].mean()) if len(lo) else np.nan,
                "mean_resid_high": float(np.nanmean(hi["imm_resid"])) if len(hi) else np.nan,
                "mean_resid_low": float(np.nanmean(lo["imm_resid"])) if len(lo) else np.nan,
                "mean_knn_high": float(hi[knn_col].mean()) if len(hi) else np.nan,
                "mean_knn_low": float(lo[knn_col].mean()) if len(lo) else np.nan,
                "mix_high_immune": mixing_homogeneous(
                    high_xy_by_sample.get(sample, np.zeros((0, 2))),
                    immune_xy_by_sample.get(sample, np.zeros((0, 2))),
                    radius,
                ),
                "mix_low_immune": mixing_homogeneous(
                    low_xy_by_sample.get(sample, np.zeros((0, 2))),
                    immune_xy_by_sample.get(sample, np.zeros((0, 2))),
                    radius,
                ),
                "mix_all_immune": mixing_homogeneous(
                    all_tum_xy_by_sample.get(sample, np.zeros((0, 2))),
                    immune_xy_by_sample.get(sample, np.zeros((0, 2))),
                    radius,
                ),
                "mix_high_cd8nk": mixing_homogeneous(
                    high_xy_by_sample.get(sample, np.zeros((0, 2))),
                    cd8nk_xy_by_sample.get(sample, np.zeros((0, 2))),
                    radius,
                ),
            }
        )
    sec = pd.DataFrame(section_rows)
    sec_w = paired_wilcoxon(sec["mean_imm_high"].to_numpy(), sec["mean_imm_low"].to_numpy())
    sec_w_frac = paired_wilcoxon(sec["mean_frac_high"].to_numpy(), sec["mean_frac_low"].to_numpy())
    sec_w_resid = paired_wilcoxon(sec["mean_resid_high"].to_numpy(), sec["mean_resid_low"].to_numpy())
    sec_w_knn = paired_wilcoxon(sec["mean_knn_high"].to_numpy(), sec["mean_knn_low"].to_numpy())

    don = (
        sec.groupby("patient", sort=True)[
            [
                "mean_imm_high",
                "mean_imm_low",
                "mean_frac_high",
                "mean_frac_low",
                "mean_resid_high",
                "mean_resid_low",
                "mean_knn_high",
                "mean_knn_low",
            ]
        ]
        .mean()
        .reset_index()
    )
    don_w = paired_wilcoxon(don["mean_imm_high"].to_numpy(), don["mean_imm_low"].to_numpy())
    don_w_frac = paired_wilcoxon(don["mean_frac_high"].to_numpy(), don["mean_frac_low"].to_numpy())

    mix_hi = float(np.nanmean(sec["mix_high_immune"]))
    mix_lo = float(np.nanmean(sec["mix_low_immune"]))
    mix_all = float(np.nanmean(sec["mix_all_immune"]))

    return {
        "section_table": sec,
        "donor_table": don,
        "section": sec_w,
        "section_frac": sec_w_frac,
        "section_resid": sec_w_resid,
        "section_knn": sec_w_knn,
        "donor": don_w,
        "donor_frac": don_w_frac,
        "mix_high": mix_hi,
        "mix_low": mix_lo,
        "mix_all": mix_all,
        "n_tumor_high": int(cell["is_high"].sum()),
        "n_tumor_low": int(cell["is_low"].sum()),
    }


def july_like_score(rec: dict) -> tuple:
    """Rank closeness to July family: 8/8 down, 5/5 down, ~3.15→2.05, mixing 0.05/0.17/0.21."""
    s = rec["section"]
    d = rec["donor"]
    delta = s["delta"]
    target_delta = 2.05 - 3.15
    mix_vec = np.array([rec["mix_high"], rec["mix_low"], rec["mix_all"]], dtype=float)
    mix_tgt = np.array([0.05, 0.17, 0.21])
    mix_err = float(np.nanmean(np.abs(mix_vec - mix_tgt))) if np.isfinite(mix_vec).any() else 9.0
    high_err = abs((s["mean_high"] if np.isfinite(s["mean_high"]) else 99) - 2.05)
    low_err = abs((s["mean_low"] if np.isfinite(s["mean_low"]) else 99) - 3.15)
    return (
        -int(s["n_down"]),
        -int(d["n_down"]),
        abs(delta - target_delta) if np.isfinite(delta) else 99.0,
        high_err + low_err,
        mix_err,
    )


def plot_paired_sections(sec: pd.DataFrame, title: str, y_high: str, y_low: str, ylabel: str, fname: str) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    xs = np.array([0, 1], dtype=float)
    colors = {
        "P5": "#1b9e77",
        "P6": "#d95f02",
        "P9": "#7570b3",
        "P12": "#e7298a",
        "P13": "#66a61e",
    }
    for _, row in sec.iterrows():
        yh, yl = row[y_high], row[y_low]
        if not (np.isfinite(yh) and np.isfinite(yl)):
            continue
        c = colors.get(row["patient"], "0.4")
        ax.plot(xs, [yl, yh], "-", color=c, lw=1.4, alpha=0.85)
        ax.scatter([0], [yl], color=c, s=42, zorder=3)
        ax.scatter([1], [yh], color=c, s=42, zorder=3)
        ax.text(1.04, yh, row["sample"], fontsize=7, va="center", color=c)
    ax.set_xticks([0, 1], ["CLDN4-low tumor", "CLDN4-high tumor"])
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.grid(True, axis="y", alpha=0.3)
    handles = [plt.Line2D([0], [0], color=colors[p], lw=2, label=p) for p in ("P5", "P6", "P9", "P12", "P13")]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname + ".png"), dpi=180)
    fig.savefig(os.path.join(FIG, fname + ".pdf"))
    plt.close(fig)


def plot_mixing(sec: pd.DataFrame, title: str, fname: str) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    samples = [s for s in SAMPLES if s in set(sec["sample"])]
    x = np.arange(len(samples))
    w = 0.25
    vals = []
    for col, lab, shift, color in (
        ("mix_high_immune", "CLDN4-high × immune", -w, "#b2182b"),
        ("mix_low_immune", "CLDN4-low × immune", 0.0, "#2166ac"),
        ("mix_all_immune", "all tumor × immune", w, "#4d4d4d"),
    ):
        y = [float(sec.loc[sec["sample"] == s, col].iloc[0]) if (sec["sample"] == s).any() else np.nan for s in samples]
        vals.append((lab, y, color, shift))
        ax.bar(x + shift, y, width=w, color=color, label=lab, alpha=0.85)
    ax.set_xticks(x, samples, rotation=35, ha="right")
    ax.set_ylabel("Homogeneous mixing (40 µm)")
    ax.set_title(title, fontsize=10)
    ax.axhline(0.05, ls="--", color="#b2182b", lw=0.8, alpha=0.5)
    ax.axhline(0.17, ls="--", color="#2166ac", lw=0.8, alpha=0.5)
    ax.axhline(0.21, ls="--", color="#4d4d4d", lw=0.8, alpha=0.5)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname + ".png"), dpi=180)
    fig.savefig(os.path.join(FIG, fname + ".pdf"))
    plt.close(fig)


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def write_results(
    primary_key: tuple,
    grid_rows: list[dict],
    primary: dict,
    closest: dict,
    inventory: dict,
    author_used: bool,
) -> None:
    gdf = pd.DataFrame(grid_rows)
    gdf.to_csv(os.path.join(TAB, "parameter_grid.csv"), index=False)
    primary["section_table"].to_csv(os.path.join(TAB, "primary_section_means.csv"), index=False)
    primary["donor_table"].to_csv(os.path.join(TAB, "primary_donor_means.csv"), index=False)

    s = primary["section"]
    d = primary["donor"]
    sf = primary["section_frac"]
    sr = primary["section_resid"]
    sk = primary["section_knn"]
    recovered = gdf[(gdf["n_down_section"] == 8) & (gdf["n_down_donor"] >= 4)]
    july8 = gdf[gdf["n_down_section"] == 8]

    lines = []
    lines.append("# CosMx NSCLC July-style neighbor retune (CLDN4-only)")
    lines.append("")
    lines.append("Official NanoString/Bruker CosMx NSCLC FFPE 960-plex (He et al. 2022): **all 8 sections / 5 patients** (Lung5_Rep1/2/3, Lung6, Lung9_Rep1/2, Lung12, Lung13). Additive. CLDN4-only. No private 8-KL.")
    lines.append("")
    lines.append("July PPT target metric family (not re-used as numbers here): CD8/NK neighbors 3.15→2.05, 8/8 sections P=0.008, 5/5 donors P=0.031, mixing 0.05/0.17/0.21. Primary readout is **mean immune-neighbor count**, not nearest-µm to CD8.")
    lines.append("")
    lines.append("## Pre-specified primary")
    lines.append("")
    lines.append("Declared before inspecting the grid:")
    lines.append("")
    lines.append("- CLDN4 cut: **median-split among tumor cells**")
    lines.append("- Radius: **40 µm**")
    lines.append("- Immune: **CD8+NK** (CD8A+ / author T CD8; NKG7+ or author NK)")
    lines.append("- Tumor gate: **author tumor if labels are present, else PanCK/KRT8/EPCAM high AND not CD8**")
    lines.append("- Tumor cells are **not** required to have CD8A=0")
    lines.append("- FOVs with <50 tumor cells or <20 CD8 dropped (no contrast)")
    lines.append("- Unit of test: **section (n=8)** and **donor (n=5)** paired Wilcoxon on section-mean (donor-mean) neighbor counts")
    lines.append("")
    lines.append(f"Author cell-type labels used: **{'yes' if author_used else 'no (official flat metadata has no cell_type column; primary uses the PanCK/KRT8/EPCAM marker gate)'}**.")
    lines.append("")
    lines.append(f"Inventory: {inventory['n_cells']:,} QC cells; {inventory['n_tumor_primary']:,} primary-gate tumor cells; {inventory['n_cd8']:,} CD8; {inventory['n_nk']:,} NK; FOVs kept {inventory['n_fov_kept']}/{inventory['n_fov']}. Panel genes (Lung5_Rep1) = {inventory['n_panel']}. CLDN4 present on all 8 exprMats.")
    lines.append("")
    lines.append("### Primary numbers")
    lines.append("")
    lines.append("| Unit | mean neighbors CLDN4-high | mean neighbors CLDN4-low | Δ (high−low) | n_down / n | paired Wilcoxon P |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    lines.append(
        f"| section n=8 | {s['mean_high']:.3f} | {s['mean_low']:.3f} | {s['delta']:.3f} | {s['n_down']}/{s['n']} | {fmt_p(s['p'])} |"
    )
    lines.append(
        f"| donor n=5 | {d['mean_high']:.3f} | {d['mean_low']:.3f} | {d['delta']:.3f} | {d['n_down']}/{d['n']} | {fmt_p(d['p'])} |"
    )
    lines.append("")
    lines.append(
        f"Density-normalized neighborhood fraction (CD8+NK / all neighbors, 40 µm): high {sf['mean_high']:.4f} vs low {sf['mean_low']:.4f}, {sf['n_down']}/{sf['n']} down, P={fmt_p(sf['p'])}."
    )
    lines.append(
        f"Neighbor count residualized on local epithelial density (same 40 µm): high {sr['mean_high']:.3f} vs low {sr['mean_low']:.3f}, {sr['n_down']}/{sr['n']} down, P={fmt_p(sr['p'])}."
    )
    lines.append(
        f"k=10 nearest-neighbor CD8+NK count (companion, not primary): high {sk['mean_high']:.3f} vs low {sk['mean_low']:.3f}, {sk['n_down']}/{sk['n']} down, P={fmt_p(sk['p'])}."
    )
    lines.append("")
    lines.append(
        f"Mixing (homogeneous, 40 µm, section-mean): CLDN4-high×CD8+NK **{primary['mix_high']:.3f}**, CLDN4-low×CD8+NK **{primary['mix_low']:.3f}**, all-tumor×CD8+NK **{primary['mix_all']:.3f}** (July family 0.05 / 0.17 / 0.21)."
    )
    lines.append("")
    lines.append("Per-section primary means:")
    lines.append("")
    lines.append("| section | donor | n_high | n_low | neighbors high | neighbors low | Δ | mix high | mix low |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for _, row in primary["section_table"].iterrows():
        lines.append(
            f"| {row['sample']} | {row['patient']} | {int(row['n_high'])} | {int(row['n_low'])} | "
            f"{row['mean_imm_high']:.3f} | {row['mean_imm_low']:.3f} | "
            f"{row['mean_imm_high']-row['mean_imm_low']:.3f} | "
            f"{row['mix_high_immune']:.3f} | {row['mix_low_immune']:.3f} |"
        )
    lines.append("")
    lines.append("Figures: `results/cosmx_july_neighbor/figures/primary_paired_sections.png`, `primary_mixing.png`.")
    lines.append("")
    lines.append("## Parameter grid")
    lines.append("")
    lines.append("CLDN4 cut × tumor gate × immune × radius. Test = paired Wilcoxon on section-mean radius neighbor counts. `n_down` = sections where CLDN4-high has fewer immune neighbors than CLDN4-low.")
    lines.append("")
    lines.append("| cut | tumor gate | immune | radius µm | mean high | mean low | Δ | section n_down/8 | section P | donor n_down/5 | donor P | mix high/low/all |")
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for _, row in gdf.sort_values(["n_down_section", "n_down_donor", "delta"], ascending=[False, False, True]).iterrows():
        star = " **← primary**" if (row["cldn4_cut"], row["tumor_gate"], row["immune"], int(row["radius_um"])) == primary_key else ""
        eight = " **8/8**" if int(row["n_down_section"]) == 8 else ""
        lines.append(
            f"| {row['cldn4_cut']} | {row['tumor_gate']} | {row['immune']} | {int(row['radius_um'])} | "
            f"{row['mean_high']:.3f} | {row['mean_low']:.3f} | {row['delta']:.3f} | "
            f"{int(row['n_down_section'])}/8{eight} | {fmt_p(row['p_section'])} | "
            f"{int(row['n_down_donor'])}/5 | {fmt_p(row['p_donor'])} | "
            f"{row['mix_high']:.3f}/{row['mix_low']:.3f}/{row['mix_all']:.3f} |{star}"
        )
    lines.append("")
    if len(july8):
        lines.append("### Grid cells with 8/8 section-level down in neighbors")
        lines.append("")
        for _, row in july8.sort_values(["n_down_donor", "delta"]).iterrows():
            lines.append(
                f"- cut={row['cldn4_cut']}, tumor={row['tumor_gate']}, immune={row['immune']}, r={int(row['radius_um'])} µm: "
                f"{row['mean_low']:.3f}→{row['mean_high']:.3f} (Δ={row['delta']:.3f}), "
                f"section P={fmt_p(row['p_section'])}, donor {int(row['n_down_donor'])}/5 P={fmt_p(row['p_donor'])}."
            )
        lines.append("")
    else:
        c = closest
        lines.append("### Closest grid cell to the July 8/8-down family")
        lines.append("")
        lines.append(
            f"No grid cell was 8/8 down. Closest by (section n_down, donor n_down, Δ vs −1.10, mixing): "
            f"**cut={c['cldn4_cut']}, tumor={c['tumor_gate']}, immune={c['immune']}, r={int(c['radius_um'])} µm** — "
            f"neighbors {c['mean_low']:.3f}→{c['mean_high']:.3f} (Δ={c['delta']:.3f}), "
            f"section {int(c['n_down_section'])}/8 P={fmt_p(c['p_section'])}, "
            f"donor {int(c['n_down_donor'])}/5 P={fmt_p(c['p_donor'])}, "
            f"mixing {c['mix_high']:.3f}/{c['mix_low']:.3f}/{c['mix_all']:.3f}."
        )
        lines.append("")
    lines.append("## Methods notes")
    lines.append("")
    lines.append("- Coordinates: `CenterX/Y_global_px` × 0.18 µm/px (CosMx NSCLC prototype).")
    lines.append("- CD8 = CD8A+ or author T CD8; NK = NKG7+ or author NK. CD8+NK is the union.")
    lines.append("- Marker tumor = z-scored Mean.PanCK + KRT8 + EPCAM composite ≥ section-wide median, and not CD8.")
    lines.append("- Residualization: per-section linear residual of immune-neighbor count on local epithelial-neighbor count (same radius); KRT8 neighborhood mean used if epithelial count is degenerate.")
    lines.append("- Mixing: homogeneous score n_AB/(n_AB+n_AA+n_BB) among CLDN4-arm tumor ∪ immune cells, radius-matched.")
    lines.append("- k-NN is reported as a companion (k=10/20 among all cells, count of CD8/NK in that neighborhood).")
    lines.append("- Data: official S3 `SMI-Compressed/{sample}/{sample}+SMI+Flat+data.tar.gz` (exprMat + metadata + fov_positions).")
    lines.append("")
    path = os.path.join(ROOT, "RESULTS.md")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {path}", flush=True)


def main() -> int:
    ensure_dirs()
    missing = [s for s in SAMPLES if not os.path.isdir(os.path.join(DATA, s))]
    if missing:
        raise SystemExit(f"missing samples {missing}; run scripts/download_cosmx_nsclc.py first")

    genes_by = {s: panel_genes(s) for s in SAMPLES}
    stop_if_cldn4_absent(genes_by)
    author_tab = load_optional_author_table()

    frames = [load_sample(s, author_tab) for s in SAMPLES]
    df = pd.concat(frames, ignore_index=True)
    df = assign_base_labels(df)
    author_used = bool(df["used_author_labels"].any()) and bool(df["is_tumor_author"].any())

    # Primary tumor gate
    primary_gate = "author_or_marker" if author_used else "marker"
    if PRIMARY["tumor_gate"] == "author_or_marker":
        tgate_primary = primary_gate
    else:
        tgate_primary = PRIMARY["tumor_gate"]

    is_t_primary = tumor_mask(df, tgate_primary)
    is_cd8 = df["is_cd8"].to_numpy()
    keep_fov = fov_keep_mask(df, is_t_primary, is_cd8)
    df = df.loc[keep_fov].reset_index(drop=True)
    df = assign_base_labels(df)
    is_t_primary = tumor_mask(df, tgate_primary)

    inventory = {
        "n_cells": int(len(df)),
        "n_tumor_primary": int(is_t_primary.sum()),
        "n_cd8": int(df["is_cd8"].sum()),
        "n_nk": int(df["is_nk"].sum()),
        "n_fov": int(df.groupby(["sample", "fov"]).ngroups),
        "n_fov_kept": int(df.groupby(["sample", "fov"]).ngroups),
        "n_panel": int(df["_n_panel_genes"].iloc[0]),
        "author_used": author_used,
        "primary_tumor_gate": tgate_primary,
    }
    with open(os.path.join(TAB, "inventory.json"), "w") as fh:
        json.dump(inventory, fh, indent=2)

    print(
        f"inventory cells={inventory['n_cells']:,} tumor={inventory['n_tumor_primary']:,} "
        f"CD8={inventory['n_cd8']:,} NK={inventory['n_nk']:,} author={author_used} gate={tgate_primary}",
        flush=True,
    )

    # Precompute neighbors once per (tumor_gate, immune) — the expensive step
    cache: dict[tuple[str, str], pd.DataFrame] = {}
    xy = xy_um(df)
    grid_rows = []
    details: dict[tuple, dict] = {}

    gates_needed = set(TUMOR_GATES) | {tgate_primary}
    for gate, immune in product(gates_needed, IMMUNE_SETS):
        is_t = tumor_mask(df, gate)
        is_i = immune_mask(df, immune)
        is_epi = tumor_mask(df, "marker") | df["is_tumor_author"].to_numpy()
        print(f"neighbors gate={gate} immune={immune} n_tumor={is_t.sum():,} n_imm={is_i.sum():,}", flush=True)
        cache[(gate, immune)] = compute_cell_neighbors(df, is_t, is_i, is_epi)

    # Grid + primary
    for cut, gate, immune, radius in product(CLDN4_CUTS, TUMOR_GATES, IMMUNE_SETS, RADII_UM):
        is_t = tumor_mask(df, gate)
        neigh = cache[(gate, immune)]
        if neigh.empty or is_t.sum() < 100:
            continue
        # high/low among tumor cells (full df index)
        high = np.zeros(len(df), dtype=bool)
        low = np.zeros(len(df), dtype=bool)
        for sample in SAMPLES:
            sm = (df["sample"].to_numpy() == sample) & is_t
            h, l = cldn4_high_low(df.loc[sm, "cldn4_ln"].to_numpy(), cut)
            idx = np.flatnonzero(sm)
            high[idx] = h
            low[idx] = l
        # restrict neighbor table tumor cells that are high or low
        immune_xy = {}
        high_xy = {}
        low_xy = {}
        all_xy = {}
        cd8nk_xy = {}
        is_i = immune_mask(df, immune)
        for sample in SAMPLES:
            sm = df["sample"].to_numpy() == sample
            immune_xy[sample] = xy[sm & is_i]
            high_xy[sample] = xy[sm & high]
            low_xy[sample] = xy[sm & low]
            all_xy[sample] = xy[sm & is_t]
            cd8nk_xy[sample] = xy[sm & df["is_cd8_nk"].to_numpy()]
        rec = summarize_setting(
            df,
            neigh,
            is_t,
            high,
            low,
            radius,
            PRIMARY["knn_k"],
            immune_xy,
            high_xy,
            low_xy,
            all_xy,
            cd8nk_xy,
        )
        key = (cut, gate, immune, int(radius))
        details[key] = rec
        grid_rows.append(
            {
                "cldn4_cut": cut,
                "tumor_gate": gate,
                "immune": immune,
                "radius_um": int(radius),
                "mean_high": rec["section"]["mean_high"],
                "mean_low": rec["section"]["mean_low"],
                "delta": rec["section"]["delta"],
                "n_down_section": rec["section"]["n_down"],
                "n_up_section": rec["section"]["n_up"],
                "p_section": rec["section"]["p"],
                "n_down_donor": rec["donor"]["n_down"],
                "n_up_donor": rec["donor"]["n_up"],
                "p_donor": rec["donor"]["p"],
                "mix_high": rec["mix_high"],
                "mix_low": rec["mix_low"],
                "mix_all": rec["mix_all"],
                "n_tumor_high": rec["n_tumor_high"],
                "n_tumor_low": rec["n_tumor_low"],
                "frac_high": rec["section_frac"]["mean_high"],
                "frac_low": rec["section_frac"]["mean_low"],
                "resid_high": rec["section_resid"]["mean_high"],
                "resid_low": rec["section_resid"]["mean_low"],
                "knn_high": rec["section_knn"]["mean_high"],
                "knn_low": rec["section_knn"]["mean_low"],
            }
        )
        print(
            f"grid {cut}/{gate}/{immune}/r{radius}: "
            f"{rec['section']['mean_low']:.3f}→{rec['section']['mean_high']:.3f} "
            f"{rec['section']['n_down']}/8 P={rec['section']['p']}",
            flush=True,
        )

    # Primary using author_or_marker (may equal marker)
    gate_p = tgate_primary
    cut_p, immune_p, rad_p = PRIMARY["cldn4_cut"], PRIMARY["immune"], PRIMARY["radius_um"]
    primary_key = (cut_p, gate_p, immune_p, int(rad_p))
    if primary_key not in details:
        # compute explicitly
        is_t = tumor_mask(df, gate_p)
        neigh = cache[(gate_p, immune_p)]
        high = np.zeros(len(df), dtype=bool)
        low = np.zeros(len(df), dtype=bool)
        for sample in SAMPLES:
            sm = (df["sample"].to_numpy() == sample) & is_t
            h, l = cldn4_high_low(df.loc[sm, "cldn4_ln"].to_numpy(), cut_p)
            idx = np.flatnonzero(sm)
            high[idx] = h
            low[idx] = l
        immune_xy, high_xy, low_xy, all_xy, cd8nk_xy = {}, {}, {}, {}, {}
        is_i = immune_mask(df, immune_p)
        for sample in SAMPLES:
            sm = df["sample"].to_numpy() == sample
            immune_xy[sample] = xy[sm & is_i]
            high_xy[sample] = xy[sm & high]
            low_xy[sample] = xy[sm & low]
            all_xy[sample] = xy[sm & is_t]
            cd8nk_xy[sample] = xy[sm & df["is_cd8_nk"].to_numpy()]
        details[primary_key] = summarize_setting(
            df, neigh, is_t, high, low, rad_p, PRIMARY["knn_k"],
            immune_xy, high_xy, low_xy, all_xy, cd8nk_xy,
        )
        # also add a labeled row if gate is author_or_marker
        rec = details[primary_key]
        grid_rows.append(
            {
                "cldn4_cut": cut_p,
                "tumor_gate": gate_p,
                "immune": immune_p,
                "radius_um": int(rad_p),
                "mean_high": rec["section"]["mean_high"],
                "mean_low": rec["section"]["mean_low"],
                "delta": rec["section"]["delta"],
                "n_down_section": rec["section"]["n_down"],
                "n_up_section": rec["section"]["n_up"],
                "p_section": rec["section"]["p"],
                "n_down_donor": rec["donor"]["n_down"],
                "n_up_donor": rec["donor"]["n_up"],
                "p_donor": rec["donor"]["p"],
                "mix_high": rec["mix_high"],
                "mix_low": rec["mix_low"],
                "mix_all": rec["mix_all"],
                "n_tumor_high": rec["n_tumor_high"],
                "n_tumor_low": rec["n_tumor_low"],
                "frac_high": rec["section_frac"]["mean_high"],
                "frac_low": rec["section_frac"]["mean_low"],
                "resid_high": rec["section_resid"]["mean_high"],
                "resid_low": rec["section_resid"]["mean_low"],
                "knn_high": rec["section_knn"]["mean_high"],
                "knn_low": rec["section_knn"]["mean_low"],
            }
        )

    primary = details[primary_key]
    gdf = pd.DataFrame(grid_rows).drop_duplicates(["cldn4_cut", "tumor_gate", "immune", "radius_um"])
    closest_row = min(gdf.to_dict("records"), key=lambda r: (
        -int(r["n_down_section"]),
        -int(r["n_down_donor"]),
        abs((r["delta"] if np.isfinite(r["delta"]) else 99) - (2.05 - 3.15)),
        abs((r["mean_high"] if np.isfinite(r["mean_high"]) else 99) - 2.05)
        + abs((r["mean_low"] if np.isfinite(r["mean_low"]) else 99) - 3.15),
    ))

    plot_paired_sections(
        primary["section_table"],
        f"Primary: median CLDN4, {tgate_primary} tumor, CD8+NK, 40 µm",
        "mean_imm_high",
        "mean_imm_low",
        "Mean CD8+NK neighbors (40 µm)",
        "primary_paired_sections",
    )
    plot_mixing(primary["section_table"], "Homogeneous mixing at 40 µm (primary gates)", "primary_mixing")
    plot_paired_sections(
        primary["section_table"],
        "Neighborhood fraction (CD8+NK / all neighbors, 40 µm)",
        "mean_frac_high",
        "mean_frac_low",
        "Mean immune fraction of neighbors",
        "primary_paired_fraction",
    )

    write_results(primary_key, gdf.to_dict("records"), primary, closest_row, inventory, author_used)
    print("done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
