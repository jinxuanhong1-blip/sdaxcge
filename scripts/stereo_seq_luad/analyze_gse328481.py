#!/usr/bin/env python3
"""CLDN4-only Stereo-seq LUAD analysis on public GSE328481 cell-bin H5ADs.

Per sample, at cell-bin and bin50:
  1. Spearman(CLDN4, CD8A) on log1p(CP10K)
  2. KRT8 residual: OLS residual of log1p(CLDN4) ~ log1p(KRT8), then Spearman vs CD8A
  3. Nearest CD8A+ distance from CLDN4-high vs CLDN4-low epithelial units
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import linregress, mannwhitneyu, spearmanr, wilcoxon

SAMPLES = [
    ("GSM9684206", "LUAD_P1", "GSM9684206_D06053D2.h5ad"),
    ("GSM9684207", "LUAD_P2", "GSM9684207_D06047C3.h5ad"),
    ("GSM9684208", "LUAD_P3", "GSM9684208_D06047F6.h5ad"),
    ("GSM9684209", "LUAD_P4", "GSM9684209_D06047E1.h5ad"),
    ("GSM9684210", "LUAD_P5", "GSM9684210_D06050A2.h5ad"),
    ("GSM9684211", "LUAD_P6", "GSM9684211_D06047A2.h5ad"),
    ("GSM9684212", "LUAD_P7", "GSM9684212_D06050C2.h5ad"),
    ("GSM9684213", "LUAD_P8", "GSM9684213_D06047D4.h5ad"),
    ("GSM9684214", "LUAD_P9", "GSM9684214_D06047E2.h5ad"),
    ("GSM9684215", "LUAD_P10", "GSM9684215_D06050D4.h5ad"),
    ("GSM9684216", "LUAD_P11", "GSM9684216_D06050E4.h5ad"),
]

EPI_LABELS = {"cancer_cell", "epi"}
GENES = ("CLDN4", "CD8A", "KRT8")
MIN_GENES_CELL = 50
MIN_GENES_BIN50 = 100
MIN_UMI_CELL = 20
MIN_UMI_BIN50 = 50


def _decode(x):
    if isinstance(x, bytes):
        return x.decode()
    return str(x)


def _categorical(group):
    cats = group["categories"]
    if isinstance(cats, h5py.Group):
        labels = [_decode(v) for v in cats["values"][:]]
    else:
        labels = [_decode(v) for v in cats[:]]
    codes = group["codes"][:]
    out = np.array(["NA"] * len(codes), dtype=object)
    ok = (codes >= 0) & (codes < len(labels))
    out[ok] = np.array(labels, dtype=object)[codes[ok]]
    return out


def load_needed(path: Path) -> dict:
    with h5py.File(path, "r") as f:
        var = np.array([_decode(v) for v in f["var/_index"][:]])
        idx = {g: int(np.where(var == g)[0][0]) for g in GENES}
        X = sparse.csr_matrix(
            (f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]),
            shape=(f["X/indptr"].shape[0] - 1, var.size),
        )
        counts = np.vstack([np.asarray(X[:, idx[g]].todense()).ravel() for g in GENES]).T.astype(np.float64)
        n_umi = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        n_genes = np.diff(X.indptr).astype(np.int32)
        anno = _categorical(f["obs/anno"])
        x = f["obs/x"][:].astype(np.float64)
        y = f["obs/y"][:].astype(np.float64)
        bin50_id = _categorical(f["obs/bin50_location_id"])
        bin50_x = f["obs/bin50_x"][:].astype(np.float64)
        bin50_y = f["obs/bin50_y"][:].astype(np.float64)
        resolution_nm = int(f["uns/resolution"][()])
    return {
        "counts": counts,
        "n_umi": n_umi,
        "n_genes": n_genes,
        "anno": anno,
        "x": x,
        "y": y,
        "bin50_id": bin50_id,
        "bin50_x": bin50_x,
        "bin50_y": bin50_y,
        "resolution_nm": resolution_nm,
        "n_cells_raw": counts.shape[0],
    }


def log_cp10k(counts: np.ndarray, n_umi: np.ndarray) -> np.ndarray:
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    return np.log1p(counts * scale[:, None])


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    ok = np.isfinite(x) & np.isfinite(y)
    resid = np.full_like(y, np.nan, dtype=np.float64)
    if ok.sum() < 10:
        return resid
    fit = linregress(x[ok], y[ok])
    resid[ok] = y[ok] - (fit.intercept + fit.slope * x[ok])
    return resid


def spearman_safe(a: np.ndarray, b: np.ndarray) -> tuple[float, float, int]:
    ok = np.isfinite(a) & np.isfinite(b)
    n = int(ok.sum())
    if n < 10 or np.unique(a[ok]).size < 2 or np.unique(b[ok]).size < 2:
        return np.nan, np.nan, n
    r, p = spearmanr(a[ok], b[ok])
    return float(r), float(p), n


def nearest_cd8_stats(
    xy_um: np.ndarray,
    epi: np.ndarray,
    cldn4_raw: np.ndarray,
    cldn4_log: np.ndarray,
    cd8_ref: np.ndarray,
) -> dict:
    """CLDN4-high = epithelial with CLDN4>0 (top half of detected if n_detected>=40).

    CD8 reference excludes epithelial units so mixed bins cannot self-match.
    """
    out = {
        "n_epi": int(epi.sum()),
        "n_cd8_ref": int(cd8_ref.sum()),
        "n_cldn4_high_epi": 0,
        "n_cldn4_low_epi": 0,
        "median_nn_um_high": np.nan,
        "median_nn_um_low": np.nan,
        "delta_high_minus_low_um": np.nan,
        "mw_u": np.nan,
        "mw_p": np.nan,
    }
    detected = epi & (cldn4_raw > 0)
    low = epi & (cldn4_raw == 0)
    if detected.sum() >= 40:
        med = np.median(cldn4_log[detected])
        high = detected & (cldn4_log >= med)
    else:
        high = detected
    out["n_cldn4_high_epi"] = int(high.sum())
    out["n_cldn4_low_epi"] = int(low.sum())
    if high.sum() < 10 or low.sum() < 10 or cd8_ref.sum() < 5:
        return out
    tree = cKDTree(xy_um[cd8_ref])
    d_high, _ = tree.query(xy_um[high], k=1)
    d_low, _ = tree.query(xy_um[low], k=1)
    out["median_nn_um_high"] = float(np.median(d_high))
    out["median_nn_um_low"] = float(np.median(d_low))
    out["delta_high_minus_low_um"] = out["median_nn_um_high"] - out["median_nn_um_low"]
    u, p = mannwhitneyu(d_high, d_low, alternative="two-sided")
    out["mw_u"] = float(u)
    out["mw_p"] = float(p)
    return out


def metrics_for_table(expr: np.ndarray, n_umi: np.ndarray, n_genes: np.ndarray, anno, x, y, resolution_nm, min_genes, min_umi) -> dict:
    keep = (n_genes >= min_genes) & (n_umi >= min_umi)
    expr = expr[keep]
    n_umi = n_umi[keep]
    anno = np.asarray(anno)[keep]
    x = x[keep]
    y = y[keep]
    log = log_cp10k(expr, n_umi)
    cldn4, cd8a, krt8 = log[:, 0], log[:, 1], log[:, 2]
    resid = residualize(cldn4, krt8)
    epi = np.isin(anno, list(EPI_LABELS))
    xy_um = np.column_stack([x, y]) * (resolution_nm / 1000.0)
    cd8_pos = expr[:, 1] > 0
    cd8_ref = cd8_pos & (~epi)
    nonzero = (expr[:, 0] > 0) | (expr[:, 1] > 0)
    rho_all = spearman_safe(cldn4, cd8a)
    rho_epi = spearman_safe(cldn4[epi], cd8a[epi])
    rho_nz = spearman_safe(cldn4[nonzero], cd8a[nonzero])
    rho_resid = spearman_safe(resid, cd8a)
    rho_resid_epi = spearman_safe(resid[epi], cd8a[epi])
    nn = nearest_cd8_stats(xy_um, epi, expr[:, 0], cldn4, cd8_ref)
    return {
        "n_units": int(keep.sum()),
        "n_epithelial": int(epi.sum()),
        "n_cd8_pos": int(cd8_pos.sum()),
        "n_cd8_ref_non_epi": int(cd8_ref.sum()),
        "frac_cd8_pos": float(cd8_pos.mean()) if keep.sum() else np.nan,
        "frac_cldn4_pos": float((expr[:, 0] > 0).mean()) if keep.sum() else np.nan,
        "frac_cldn4_pos_epi": float((expr[epi, 0] > 0).mean()) if epi.any() else np.nan,
        "n_both_pos": int(((expr[:, 0] > 0) & cd8_pos).sum()),
        "frac_cd8_among_cldn4_pos": float(cd8_pos[expr[:, 0] > 0].mean()) if (expr[:, 0] > 0).any() else np.nan,
        "frac_cd8_among_cldn4_neg": float(cd8_pos[expr[:, 0] == 0].mean()) if (expr[:, 0] == 0).any() else np.nan,
        "rho_cldn4_cd8a": rho_all[0],
        "p_cldn4_cd8a": rho_all[1],
        "n_cldn4_cd8a": rho_all[2],
        "rho_cldn4_cd8a_epi": rho_epi[0],
        "p_cldn4_cd8a_epi": rho_epi[1],
        "n_cldn4_cd8a_epi": rho_epi[2],
        "rho_cldn4_cd8a_nonzero": rho_nz[0],
        "p_cldn4_cd8a_nonzero": rho_nz[1],
        "n_cldn4_cd8a_nonzero": rho_nz[2],
        "rho_krt8_residual_cd8a": rho_resid[0],
        "p_krt8_residual_cd8a": rho_resid[1],
        "n_krt8_residual_cd8a": rho_resid[2],
        "rho_krt8_residual_cd8a_epi": rho_resid_epi[0],
        "p_krt8_residual_cd8a_epi": rho_resid_epi[1],
        "n_krt8_residual_cd8a_epi": rho_resid_epi[2],
        **nn,
    }


def aggregate_bin50(d: dict) -> dict:
    df = pd.DataFrame(
        {
            "bin50": d["bin50_id"],
            "CLDN4": d["counts"][:, 0],
            "CD8A": d["counts"][:, 1],
            "KRT8": d["counts"][:, 2],
            "n_umi": d["n_umi"],
            "n_genes": d["n_genes"],
            "x": d["bin50_x"],
            "y": d["bin50_y"],
            "epi": np.isin(d["anno"], list(EPI_LABELS)).astype(np.int32),
        }
    )
    g = df.groupby("bin50", sort=False)
    agg = g.agg(
        CLDN4=("CLDN4", "sum"),
        CD8A=("CD8A", "sum"),
        KRT8=("KRT8", "sum"),
        n_umi=("n_umi", "sum"),
        n_genes=("n_genes", "sum"),
        x=("x", "first"),
        y=("y", "first"),
        n_cells=("epi", "size"),
        n_epi=("epi", "sum"),
    )
    counts = agg[["CLDN4", "CD8A", "KRT8"]].to_numpy()
    anno = np.where(agg["n_epi"] / agg["n_cells"] >= 0.5, "epi", "other")
    return {
        "counts": counts,
        "n_umi": agg["n_umi"].to_numpy(),
        "n_genes": agg["n_genes"].to_numpy(),
        "anno": anno,
        "x": agg["x"].to_numpy(),
        "y": agg["y"].to_numpy(),
    }


def summarize(df: pd.DataFrame, cols: list[str]) -> dict:
    out = {}
    for col in cols:
        vals = df[col].dropna().to_numpy()
        out[f"{col}_median"] = float(np.median(vals)) if vals.size else np.nan
        out[f"{col}_n"] = int(vals.size)
        if vals.size >= 6 and np.unique(np.round(vals, 12)).size > 1:
            w, p = wilcoxon(vals)
            out[f"{col}_wilcoxon_stat"] = float(w)
            out[f"{col}_wilcoxon_p"] = float(p)
        else:
            out[f"{col}_wilcoxon_stat"] = np.nan
            out[f"{col}_wilcoxon_p"] = np.nan
        out[f"{col}_n_negative"] = int((vals < 0).sum()) if vals.size else 0
        out[f"{col}_n_positive"] = int((vals > 0).sum()) if vals.size else 0
    return out


def main():
    data_dir = Path(os.environ.get("GSE328481_DIR", "/tmp/gse328481"))
    out_dir = Path("results/stereo_seq_luad")
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    rows_cell, rows_bin = [], []
    for gsm, patient, fn in SAMPLES:
        path = data_dir / fn
        print(f"analyzing {patient} {fn}", flush=True)
        d = load_needed(path)
        cell = metrics_for_table(
            d["counts"], d["n_umi"], d["n_genes"], d["anno"], d["x"], d["y"],
            d["resolution_nm"], MIN_GENES_CELL, MIN_UMI_CELL,
        )
        cell.update({"gsm": gsm, "patient": patient, "chip": fn.replace(".h5ad", ""), "unit": "cell_bin"})
        rows_cell.append(cell)

        b = aggregate_bin50(d)
        binm = metrics_for_table(
            b["counts"], b["n_umi"], b["n_genes"], b["anno"], b["x"], b["y"],
            d["resolution_nm"], MIN_GENES_BIN50, MIN_UMI_BIN50,
        )
        binm.update({"gsm": gsm, "patient": patient, "chip": fn.replace(".h5ad", ""), "unit": "bin50"})
        rows_bin.append(binm)
        del d, b

    cell_df = pd.DataFrame(rows_cell)
    bin_df = pd.DataFrame(rows_bin)
    both = pd.concat([cell_df, bin_df], ignore_index=True)
    both.to_csv(out_dir / "tables" / "per_sample_metrics.csv", index=False)

    rho_cols = [
        "rho_cldn4_cd8a",
        "rho_cldn4_cd8a_epi",
        "rho_cldn4_cd8a_nonzero",
        "rho_krt8_residual_cd8a",
        "rho_krt8_residual_cd8a_epi",
        "delta_high_minus_low_um",
    ]
    summary = {
        "cell_bin": summarize(cell_df, rho_cols),
        "bin50": summarize(bin_df, rho_cols),
        "n_samples": 11,
        "series": "GSE328481",
        "platform": "BGI Stereo-seq / Stereo-XCR-seq cell_bins",
        "epithelial_definition": "author anno cancer_cell or epi; bin50 if >=50% of cells in those labels",
        "cd8_definition": "CD8A raw count > 0; nearest-neighbor reference excludes epithelial units",
        "cldn4_high": "epithelial with CLDN4>0 (top half of detected if n_detected>=40); low = epithelial CLDN4==0",
        "krt8_residual": "OLS residual of log1p(CP10K) CLDN4 ~ log1p(CP10K) KRT8",
        "distance_unit": "micrometers (coordinate * resolution_nm/1000; resolution=500)",
        "no_private_8kl": True,
        "cldn4_only": True,
    }
    with open(out_dir / "tables" / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
