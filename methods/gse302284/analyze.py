#!/usr/bin/env python3
"""Score TACSTD2, CLDN4, TJ, and IFN/MHC-I/APM in public GSE302284.

Design (from GEO, not inferred from results):
  - 6 Cell Ranger filtered H5 libraries.
  - No sacituzumab / SKB264 / TROP2-ADC treated RNA is deposited.
  - Treatment contrast that exists: osimertinib vs vehicle in DFCI282 and PC9
    (n_library = 1 vs 1 each). Patient 357 tumor vs lymph node is residual
    neoadjuvant EGFR-TKI tissue, not a drug contrast.

SKB264 / CLDN4 thesis is taken as given. This script reports honest public
direction + n + p + log2FC on the deposited contrasts.

Cell-level p-values are exploratory: libraries are the biological unit, and
each arm has n_library = 1. Pseudobulk log2FC is the primary effect size.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from download import SAMPLES  # noqa: E402
from gene_sets import (  # noqa: E402
    EPI_MARKERS,
    FOCAL_GENES,
    IFN_ISG,
    IMMUNE_MARKERS,
    MHC1_APM,
    PRIORITY6,
    SETS,
    TARGETS,
    TJ_SIG,
)

try:
    import h5py
except ImportError as e:  # pragma: no cover
    raise SystemExit("h5py is required: pip install -r methods/gse302284/requirements.txt") from e


CONTRASTS = [
    {
        "id": "DFCI282_osi_vs_veh",
        "model": "DFCI282",
        "treat_gsm": "GSM9128162",
        "ctrl_gsm": "GSM9128163",
        "treat_label": "osimertinib",
        "ctrl_label": "vehicle",
        "kind": "treatment",
        "compartment": "all_qc",
        "note": "n_library = 1 vs 1. Not a TROP2-ADC contrast.",
    },
    {
        "id": "PC9_osi_vs_veh",
        "model": "PC9",
        "treat_gsm": "GSM9128164",
        "ctrl_gsm": "GSM9128165",
        "treat_label": "osimertinib",
        "ctrl_label": "vehicle",
        "kind": "treatment",
        "compartment": "all_qc",
        "note": (
            "n_library = 1 vs 1. GEO labels PC9_Veh as cell line PC10; "
            "paired by library name. Not a TROP2-ADC contrast."
        ),
    },
    {
        "id": "patient357_tumor_vs_LN",
        "model": "patient_357",
        "treat_gsm": "GSM9101265",
        "ctrl_gsm": "GSM9101264",
        "treat_label": "tumor_residual",
        "ctrl_label": "lymph_node_residual",
        "kind": "descriptive_site",
        "compartment": "epithelial",
        "note": (
            "Same-patient residual tissue after neoadjuvant EGFR-TKI. "
            "Not osimertinib-vs-vehicle and not TROP2-ADC."
        ),
    },
]


def read_10x_h5(path: Path) -> tuple[sparse.csc_matrix, np.ndarray, np.ndarray]:
    """Return genes x cells CSC counts, gene symbols, barcodes."""
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        data = m["data"][:]
        indices = m["indices"][:]
        indptr = m["indptr"][:]
        shape = tuple(int(x) for x in m["shape"][:])
        barcodes = np.array(
            [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in m["barcodes"][:]]
        )
        names = np.array(
            [
                x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x)
                for x in m["features"]["name"][:]
            ]
        )
        ftype = np.array(
            [
                x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x)
                for x in m["features"]["feature_type"][:]
            ]
        )
    X = sparse.csc_matrix((data, indices, indptr), shape=shape)
    keep = ftype == "Gene Expression"
    if keep.any() and not keep.all():
        X = X[keep, :]
        names = names[keep]
    # Sum duplicate symbols
    uniq, inv = np.unique(names, return_inverse=True)
    if len(uniq) != len(names):
        X = sparse.csc_matrix(X)
        rows = []
        for i in range(len(uniq)):
            idx = np.where(inv == i)[0]
            if len(idx) == 1:
                rows.append(X[idx[0], :])
            else:
                rows.append(X[idx, :].sum(axis=0))
        X = sparse.vstack(rows, format="csc")
        names = uniq
    return X.tocsc(), names, barcodes


def qc_mask(X: sparse.csc_matrix, names: np.ndarray) -> np.ndarray:
    """Author-like QC: drop bottom 10% UMI, bottom 10% genes, mito > 25%."""
    n_umi = np.asarray(X.sum(axis=0)).ravel()
    n_gene = np.asarray((X > 0).sum(axis=0)).ravel()
    mito = np.array([str(g).startswith("MT-") for g in names])
    mito_umi = np.asarray(X[mito, :].sum(axis=0)).ravel() if mito.any() else np.zeros_like(n_umi)
    mito_frac = np.divide(mito_umi, n_umi, out=np.zeros_like(n_umi, dtype=float), where=n_umi > 0)
    umi_cut = np.quantile(n_umi, 0.10)
    gene_cut = np.quantile(n_gene, 0.10)
    return (n_umi >= umi_cut) & (n_gene >= gene_cut) & (mito_frac < 0.25) & (n_umi > 0)


def log1p_cp10k(X: sparse.csc_matrix) -> sparse.csc_matrix:
    n_umi = np.asarray(X.sum(axis=0)).ravel().astype(float)
    n_umi[n_umi == 0] = 1.0
    # scale columns
    Xc = X.tocsc().astype(float)
    Xc = Xc.multiply(1.0 / n_umi)
    Xc = Xc * 1e4
    Xc.data = np.log1p(Xc.data)
    return Xc.tocsc()


def gene_index(names: np.ndarray) -> dict[str, int]:
    return {str(g): i for i, g in enumerate(names)}


def row_dense(X: sparse.spmatrix, i: int) -> np.ndarray:
    return np.asarray(X[i, :].todense()).ravel()


def set_score(logX: sparse.spmatrix, idx: dict[str, int], genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [g for g in genes if g in idx]
    if not present:
        return np.full(logX.shape[1], np.nan), []
    mat = np.vstack([row_dense(logX, idx[g]) for g in present])
    return mat.mean(axis=0), present


def epithelial_mask(logX: sparse.spmatrix, idx: dict[str, int]) -> np.ndarray:
    epi, _ = set_score(logX, idx, EPI_MARKERS)
    imm, _ = set_score(logX, idx, IMMUNE_MARKERS)
    epcam = row_dense(logX, idx["EPCAM"]) if "EPCAM" in idx else np.zeros_like(epi)
    krt = np.zeros_like(epi)
    for g in ("KRT8", "KRT18", "KRT19", "KRT7"):
        if g in idx:
            krt = np.maximum(krt, row_dense(logX, idx[g]))
    return (epi >= imm) & ((epcam > 0) | (krt > 0))


def rank_biserial(a: np.ndarray, b: np.ndarray, u: float) -> float:
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return float("nan")
    # positive => treat > ctrl (a is treat)
    return (2.0 * u) / (n1 * n2) - 1.0


def welch_safe(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return float(t), float(p)


def mwu_safe(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    if len(a) < 1 or len(b) < 1:
        return float("nan"), float("nan")
    try:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        return float("nan"), 1.0
    return float(u), float(p)


def spearman_safe(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 5:
        return float("nan"), float("nan"), n
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), n


def load_library(cache_raw: Path, rec: dict) -> dict:
    path = cache_raw / rec["h5"]
    X, names, barcodes = read_10x_h5(path)
    keep = qc_mask(X, names)
    X = X[:, keep]
    barcodes = barcodes[keep]
    logX = log1p_cp10k(X)
    idx = gene_index(names)
    n_umi = np.asarray(X.sum(axis=0)).ravel()
    return {
        "rec": rec,
        "X": X,
        "logX": logX,
        "names": names,
        "idx": idx,
        "barcodes": barcodes,
        "n_raw": int(keep.size),
        "n_qc": int(keep.sum()),
        "n_umi": n_umi,
    }


def gene_stats(treat: dict, ctrl: dict, gene: str, tmask: np.ndarray, cmask: np.ndarray) -> dict:
    out = {
        "gene": gene,
        "present_treat": gene in treat["idx"],
        "present_ctrl": gene in ctrl["idx"],
    }
    if gene not in treat["idx"] or gene not in ctrl["idx"]:
        out.update(
            {
                "n_treat": int(tmask.sum()),
                "n_ctrl": int(cmask.sum()),
                "mean_log1p_treat": float("nan"),
                "mean_log1p_ctrl": float("nan"),
                "pct_pos_treat": float("nan"),
                "pct_pos_ctrl": float("nan"),
                "cpm_treat": float("nan"),
                "cpm_ctrl": float("nan"),
                "log2FC_pseudobulk": float("nan"),
                "delta_mean_log1p": float("nan"),
                "mwu_U": float("nan"),
                "mwu_p": float("nan"),
                "rank_biserial": float("nan"),
                "welch_t": float("nan"),
                "welch_p": float("nan"),
            }
        )
        return out

    Xt = treat["X"]
    Xc = ctrl["X"]
    it, ic = treat["idx"][gene], ctrl["idx"][gene]
    ct = np.asarray(Xt[it, :].todense()).ravel()[tmask]
    cc = np.asarray(Xc[ic, :].todense()).ravel()[cmask]
    lt = row_dense(treat["logX"], it)[tmask]
    lc = row_dense(ctrl["logX"], ic)[cmask]
    umi_t = treat["n_umi"][tmask].sum()
    umi_c = ctrl["n_umi"][cmask].sum()
    cpm_t = 1e6 * ct.sum() / umi_t if umi_t > 0 else float("nan")
    cpm_c = 1e6 * cc.sum() / umi_c if umi_c > 0 else float("nan")
    u, p = mwu_safe(lt, lc)
    tstat, tp = welch_safe(lt, lc)
    out.update(
        {
            "n_treat": int(tmask.sum()),
            "n_ctrl": int(cmask.sum()),
            "mean_log1p_treat": float(lt.mean()) if len(lt) else float("nan"),
            "mean_log1p_ctrl": float(lc.mean()) if len(lc) else float("nan"),
            "pct_pos_treat": float((ct > 0).mean()) if len(ct) else float("nan"),
            "pct_pos_ctrl": float((cc > 0).mean()) if len(cc) else float("nan"),
            "cpm_treat": float(cpm_t),
            "cpm_ctrl": float(cpm_c),
            "log2FC_pseudobulk": float(math.log2((cpm_t + 1.0) / (cpm_c + 1.0))),
            "delta_mean_log1p": float(lt.mean() - lc.mean()) if len(lt) and len(lc) else float("nan"),
            "mwu_U": u,
            "mwu_p": p,
            "rank_biserial": rank_biserial(lt, lc, u) if np.isfinite(u) else float("nan"),
            "welch_t": tstat,
            "welch_p": tp,
        }
    )
    return out


def score_stats(treat: dict, ctrl: dict, name: str, genes: list[str], tmask: np.ndarray, cmask: np.ndarray) -> dict:
    st, pt = set_score(treat["logX"], treat["idx"], genes)
    sc, pc = set_score(ctrl["logX"], ctrl["idx"], genes)
    st, sc = st[tmask], sc[cmask]
    u, p = mwu_safe(st, sc)
    tstat, tp = welch_safe(st, sc)
    # pseudobulk set score: mean of per-gene log2FC
    gene_lfc = []
    for g in genes:
        if g in treat["idx"] and g in ctrl["idx"]:
            ct = np.asarray(treat["X"][treat["idx"][g], :].todense()).ravel()[tmask].sum()
            cc = np.asarray(ctrl["X"][ctrl["idx"][g], :].todense()).ravel()[cmask].sum()
            umi_t = treat["n_umi"][tmask].sum()
            umi_c = ctrl["n_umi"][cmask].sum()
            cpm_t = 1e6 * ct / umi_t
            cpm_c = 1e6 * cc / umi_c
            gene_lfc.append(math.log2((cpm_t + 1.0) / (cpm_c + 1.0)))
    return {
        "set": name,
        "n_genes_requested": len(genes),
        "n_genes_present": len(set(pt) | set(pc)),
        "genes_present": ",".join(sorted(set(pt) | set(pc))),
        "n_treat": int(tmask.sum()),
        "n_ctrl": int(cmask.sum()),
        "mean_score_treat": float(np.nanmean(st)),
        "mean_score_ctrl": float(np.nanmean(sc)),
        "delta_mean_score": float(np.nanmean(st) - np.nanmean(sc)),
        "median_gene_log2FC": float(np.median(gene_lfc)) if gene_lfc else float("nan"),
        "mean_gene_log2FC": float(np.mean(gene_lfc)) if gene_lfc else float("nan"),
        "mwu_U": u,
        "mwu_p": p,
        "rank_biserial": rank_biserial(st, sc, u) if np.isfinite(u) else float("nan"),
        "welch_t": tstat,
        "welch_p": tp,
    }


def within_library_corr(lib: dict, mask: np.ndarray) -> list[dict]:
    idx, logX = lib["idx"], lib["logX"]
    rows = []
    pairs = [
        ("TACSTD2", "CLDN4"),
        ("TACSTD2", "TJ_SIG"),
        ("TACSTD2", "IFN_ISG"),
        ("TACSTD2", "MHC1_APM"),
        ("CLDN4", "TJ_SIG"),
        ("CLDN4", "IFN_ISG"),
        ("CLDN4", "MHC1_APM"),
        ("TJ_SIG", "IFN_ISG"),
        ("TJ_SIG", "MHC1_APM"),
    ]
    scores = {}
    if "TACSTD2" in idx:
        scores["TACSTD2"] = row_dense(logX, idx["TACSTD2"])[mask]
    if "CLDN4" in idx:
        scores["CLDN4"] = row_dense(logX, idx["CLDN4"])[mask]
    scores["TJ_SIG"], _ = set_score(logX, idx, TJ_SIG)
    scores["IFN_ISG"], _ = set_score(logX, idx, IFN_ISG)
    scores["MHC1_APM"], _ = set_score(logX, idx, MHC1_APM)
    for k in ("TJ_SIG", "IFN_ISG", "MHC1_APM"):
        scores[k] = scores[k][mask]
    rec = lib["rec"]
    for a, b in pairs:
        if a not in scores or b not in scores:
            r, p, n = float("nan"), float("nan"), int(mask.sum())
        else:
            r, p, n = spearman_safe(scores[a], scores[b])
        rows.append(
            {
                "gsm": rec["gsm"],
                "library": rec["library"],
                "model": rec["model"],
                "treatment": rec["treatment"],
                "pair": f"{a}_vs_{b}",
                "n_cells": n,
                "spearman_rho": r,
                "spearman_p": p,
            }
        )
    return rows


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p == 0 or p < 1e-300:
        return "<1e-300"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_n(x: float, nd: int = 3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{nd}f}"


def direction(x: float, up_good: bool | None = None) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    if abs(x) < 1e-12:
        return "flat"
    return "up" if x > 0 else "down"


def run_contrast(libs: dict[str, dict], spec: dict) -> dict:
    treat = libs[spec["treat_gsm"]]
    ctrl = libs[spec["ctrl_gsm"]]
    if spec["compartment"] == "epithelial":
        tmask = epithelial_mask(treat["logX"], treat["idx"])
        cmask = epithelial_mask(ctrl["logX"], ctrl["idx"])
    else:
        tmask = np.ones(treat["n_qc"], dtype=bool)
        cmask = np.ones(ctrl["n_qc"], dtype=bool)

    genes_to_score = list(dict.fromkeys(FOCAL_GENES + MHC1_APM + IFN_ISG + TJ_SIG))
    gene_rows = [gene_stats(treat, ctrl, g, tmask, cmask) for g in genes_to_score]
    for r in gene_rows:
        r["contrast"] = spec["id"]
        r["model"] = spec["model"]
        r["kind"] = spec["kind"]
        r["n_library_treat"] = 1
        r["n_library_ctrl"] = 1
        r["compartment"] = spec["compartment"]

    set_rows = []
    for name, genes in SETS.items():
        row = score_stats(treat, ctrl, name, genes, tmask, cmask)
        row["contrast"] = spec["id"]
        row["model"] = spec["model"]
        row["kind"] = spec["kind"]
        row["n_library_treat"] = 1
        row["n_library_ctrl"] = 1
        row["compartment"] = spec["compartment"]
        set_rows.append(row)

    corr_rows = within_library_corr(treat, tmask) + within_library_corr(ctrl, cmask)
    return {
        "spec": spec,
        "n_treat_qc": treat["n_qc"],
        "n_ctrl_qc": ctrl["n_qc"],
        "n_treat_used": int(tmask.sum()),
        "n_ctrl_used": int(cmask.sum()),
        "n_treat_raw": treat["n_raw"],
        "n_ctrl_raw": ctrl["n_raw"],
        "gene_rows": gene_rows,
        "set_rows": set_rows,
        "corr_rows": corr_rows,
        "treat_scores": _cell_scores(treat, tmask, spec["treat_label"]),
        "ctrl_scores": _cell_scores(ctrl, cmask, spec["ctrl_label"]),
    }


def _cell_scores(lib: dict, mask: np.ndarray, label: str) -> pd.DataFrame:
    idx, logX = lib["idx"], lib["logX"]
    d = {
        "library": lib["rec"]["library"],
        "gsm": lib["rec"]["gsm"],
        "model": lib["rec"]["model"],
        "arm": label,
        "barcode": lib["barcodes"][mask],
    }
    for g in TARGETS + TJ_SIG:
        d[g] = row_dense(logX, idx[g])[mask] if g in idx else np.nan
    for name, genes in SETS.items():
        sc, _ = set_score(logX, idx, genes)
        d[name] = sc[mask]
    return pd.DataFrame(d)


def plot_score_boxes(results: list[dict], out: Path) -> None:
    specs = [r for r in results if r["spec"]["kind"] == "treatment"]
    metrics = [
        ("TACSTD2", "TACSTD2"),
        ("CLDN4", "CLDN4"),
        ("TJ_SIG", "TJ signature"),
        ("IFN_ISG", "IFN/ISG"),
        ("MHC1_APM", "MHC-I/APM"),
    ]
    fig, axes = plt.subplots(len(specs), len(metrics), figsize=(14, 6.2), sharey=False)
    if len(specs) == 1:
        axes = np.array([axes])
    for i, res in enumerate(specs):
        df = pd.concat([res["treat_scores"], res["ctrl_scores"]], ignore_index=True)
        for j, (col, title) in enumerate(metrics):
            ax = axes[i, j]
            arms = [res["spec"]["ctrl_label"], res["spec"]["treat_label"]]
            data = [df.loc[df["arm"] == a, col].dropna().to_numpy() for a in arms]
            bp = ax.boxplot(data, labels=["veh", "osi"], widths=0.55, patch_artist=True, showfliers=False)
            colors = ["#9aa4b2", "#2f6fed"]
            for patch, c in zip(bp["boxes"], colors):
                patch.set_facecolor(c)
                patch.set_alpha(0.55)
            ax.set_title(f"{res['spec']['model']}\n{title}", fontsize=9)
            ax.set_ylabel("mean log1p CP10k" if j == 0 else "")
            ax.tick_params(labelsize=8)
    fig.suptitle(
        "GSE302284 — osimertinib vs vehicle (n_library = 1 vs 1). Not a TROP2-ADC contrast.",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def plot_heatmap(gene_df: pd.DataFrame, out: Path) -> None:
    genes = [g for g in FOCAL_GENES if g in set(gene_df["gene"])]
    contrasts = ["DFCI282_osi_vs_veh", "PC9_osi_vs_veh", "patient357_tumor_vs_LN"]
    mat = []
    for g in genes:
        row = []
        for c in contrasts:
            sub = gene_df[(gene_df["gene"] == g) & (gene_df["contrast"] == c)]
            row.append(float(sub["log2FC_pseudobulk"].iloc[0]) if len(sub) else np.nan)
        mat.append(row)
    M = np.array(mat, dtype=float)
    fig, ax = plt.subplots(figsize=(6.4, 7.2))
    vmax = np.nanmax(np.abs(M)) if np.isfinite(M).any() else 1.0
    vmax = max(vmax, 0.5)
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(contrasts)))
    ax.set_xticklabels(["DFCI282\nosi/veh", "PC9\nosi/veh", "pt357\ntumor/LN"], fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=6.5)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="pseudobulk log2FC")
    ax.set_title("GSE302284 focal genes — pseudobulk log2FC")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def plot_scatter(results: list[dict], out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    for ax, res in zip(axes, results):
        df = pd.concat([res["treat_scores"], res["ctrl_scores"]], ignore_index=True)
        for arm, color, label in (
            (res["spec"]["ctrl_label"], "#9aa4b2", res["spec"]["ctrl_label"]),
            (res["spec"]["treat_label"], "#2f6fed", res["spec"]["treat_label"]),
        ):
            sub = df[df["arm"] == arm]
            ax.scatter(sub["TACSTD2"], sub["CLDN4"], s=6, alpha=0.25, c=color, label=label, linewidths=0)
        ax.set_xlabel("TACSTD2 log1p CP10k")
        ax.set_ylabel("CLDN4 log1p CP10k")
        ax.set_title(res["spec"]["model"], fontsize=10)
        ax.legend(fontsize=7, markerscale=2)
    fig.suptitle("GSE302284 — TACSTD2 vs CLDN4 (cell-level, QC cells / epithelium)")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def plot_set_bars(set_df: pd.DataFrame, out: Path) -> None:
    sub = set_df[set_df["kind"] == "treatment"].copy()
    sets = ["TJ_SIG", "IFN_ISG", "MHC1_APM", "PRIORITY6"]
    models = ["DFCI282", "PC9"]
    x = np.arange(len(sets))
    width = 0.35
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    for i, model in enumerate(models):
        vals = []
        for s in sets:
            row = sub[(sub["model"] == model) & (sub["set"] == s)]
            vals.append(float(row["mean_gene_log2FC"].iloc[0]) if len(row) else np.nan)
        ax.bar(x + (i - 0.5) * width, vals, width, label=model)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(["TJ", "IFN/ISG", "MHC-I/APM", "PRIORITY6"])
    ax.set_ylabel("mean per-gene pseudobulk log2FC\n(osi − vehicle)")
    ax.set_title("GSE302284 set direction — osimertinib vs vehicle (n_library = 1 vs 1)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def write_writeup(
    out: Path,
    gene_df: pd.DataFrame,
    set_df: pd.DataFrame,
    corr_df: pd.DataFrame,
    results: list[dict],
    lib_qc: list[dict],
) -> None:
    def g(contrast: str, gene: str) -> pd.Series:
        sub = gene_df[(gene_df["contrast"] == contrast) & (gene_df["gene"] == gene)]
        return sub.iloc[0] if len(sub) else pd.Series(dtype=float)

    def s(contrast: str, name: str) -> pd.Series:
        sub = set_df[(set_df["contrast"] == contrast) & (set_df["set"] == name)]
        return sub.iloc[0] if len(sub) else pd.Series(dtype=float)

    def line(row: pd.Series, label: str) -> str:
        return (
            f"| {label} | {int(row.get('n_treat', 0))} / {int(row.get('n_ctrl', 0))} | "
            f"{fmt_n(row.get('log2FC_pseudobulk', row.get('mean_gene_log2FC')))} | "
            f"{direction(row.get('log2FC_pseudobulk', row.get('mean_gene_log2FC')))} | "
            f"{fmt_p(row.get('mwu_p'))} | {fmt_n(row.get('rank_biserial'))} |"
        )

    d282_t2 = g("DFCI282_osi_vs_veh", "TACSTD2")
    d282_c4 = g("DFCI282_osi_vs_veh", "CLDN4")
    d282_tj = s("DFCI282_osi_vs_veh", "TJ_SIG")
    d282_ifn = s("DFCI282_osi_vs_veh", "IFN_ISG")
    d282_apm = s("DFCI282_osi_vs_veh", "MHC1_APM")
    pc9_t2 = g("PC9_osi_vs_veh", "TACSTD2")
    pc9_c4 = g("PC9_osi_vs_veh", "CLDN4")
    pc9_tj = s("PC9_osi_vs_veh", "TJ_SIG")
    pc9_ifn = s("PC9_osi_vs_veh", "IFN_ISG")
    pc9_apm = s("PC9_osi_vs_veh", "MHC1_APM")

    qc_rows = "\n".join(
        f"| {q['gsm']} | {q['library']} | {q['model']} | {q['treatment']} | {q['n_raw']} | {q['n_qc']} |"
        for q in lib_qc
    )

    def corr_line(gsm: str, pair: str) -> str:
        sub = corr_df[(corr_df["gsm"] == gsm) & (corr_df["pair"] == pair)]
        if not len(sub):
            return "NA"
        r = sub.iloc[0]
        return f"ρ={fmt_n(r['spearman_rho'])}, p={fmt_p(r['spearman_p'])}, n={int(r['n_cells'])}"

    md = f"""# C — GSE302284 public context (TROP2 ADC / EGFR DTP NSCLC)

**Slice:** `methods/gse302284/` + `results/C_GSE302284/` only.
**Thesis:** SKB264 / CLDN4 is taken as given. This slice adds a public GEO series; it does not audit the thesis.
**Accession:** [GSE302284](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302284) (Baldacci / Brea / Jänne; public 2025-07-10). Cell Ranger v7.1 filtered H5, GRCh38.

## What is (and is not) on GEO

| Question | Answer |
|---|---|
| Sacituzumab / SKB264 / TROP2-ADC treated RNA? | **No.** Six libraries only. The series abstract discusses sacituzumab govitecan (SG) *in vivo* efficacy vs TROP2 CAR-T, but **no SG / SKB264 / TROP2-ADC treated sample is deposited.** |
| Public analog of SKB264? | **Not as a treatment RNA contrast.** SG is the same antibody lineage / TROP2-ADC class as SKB264; that analog cannot be scored here because those tumors are not on GEO. |
| What *can* be scored? | **Osimertinib vs vehicle** in DFCI282 and PC9 (n_library = 1 vs 1 each), plus descriptive residual tumor vs lymph node from one neoadjuvant EGFR-TKI patient (357). |

This is extra public context around TROP2-high EGFR DTP NSCLC, **not** a SKB264 on-treatment analog.

## Design (pre-specified)

- **Primary public contrast:** osimertinib − vehicle, separately in DFCI282 and PC9.
- **Unit:** library is the biological replicate. Each arm has **n = 1 library**. Cell-level MWU *p* is exploratory (cells are not independent biological replicates). Primary effect size = **pseudobulk log2FC** = log2((CPM_treat+1)/(CPM_ctrl+1)).
- **QC:** Cell Ranger filtered H5, then author-like filters (drop bottom 10% UMI, bottom 10% genes, mito > 25%). Scrublet was not re-run.
- **Patient 357:** epithelial cells only (EPCAM/KRT vs PTPRC). Cell-line models: all QC cells.
- **TJ signature:** mean log1p(CP10k) of CLDN1, CLDN4, CLDN7, F11R, PARD3, OCLN, TJP1.
- **IFN / MHC-I / APM:** PRIORITY6 (IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A), IFN_ISG, MHC1_APM (fixed lists in `methods/gse302284/gene_sets.py`).

## Coverage

| GSM | Library | Model | Treatment | n_raw (H5) | n_QC |
|---|---|---|---|---:|---:|
{qc_rows}

GEO note: `PC9_Veh` (GSM9128165) has `cell line: PC10` / source_name PC10; library name is PC9_Veh. It is used as the paired vehicle for PC9_Osi.

## Osimertinib vs vehicle — direction (honest)

Positive log2FC = higher in osimertinib. n_cells = QC cells used; n_library = **1 vs 1**.

| Model / feature | n_cells osi / veh | log2FC | direction | cell MWU *p* | rank-biserial |
|---|---:|---:|---|---:|---:|
{line(d282_t2, "DFCI282 TACSTD2")}
{line(d282_c4, "DFCI282 CLDN4")}
{line(d282_tj, "DFCI282 TJ signature")}
{line(d282_ifn, "DFCI282 IFN/ISG")}
{line(d282_apm, "DFCI282 MHC-I/APM")}
{line(pc9_t2, "PC9 TACSTD2")}
{line(pc9_c4, "PC9 CLDN4")}
{line(pc9_tj, "PC9 TJ signature")}
{line(pc9_ifn, "PC9 IFN/ISG")}
{line(pc9_apm, "PC9 MHC-I/APM")}

Set rows use **mean per-gene pseudobulk log2FC** (not a single composite CPM). Cell MWU is on the per-cell mean log1p score.

### Analog question (TJ down / IFN up) — not testable for SKB264 here

The SKB264 thesis predicts TJ down and IFN/MHC-I up after TROP2-ADC. **GSE302284 has no TROP2-ADC RNA**, so that analog is not scored. On the deposited **osimertinib vs vehicle** contrast (extra context only):

| Model | TJ (osi vs veh) | IFN/ISG (osi vs veh) | MHC-I/APM (osi vs veh) |
|---|---|---|---|
| DFCI282 | {direction(d282_tj.get("mean_gene_log2FC"))} (mean gene log2FC {fmt_n(d282_tj.get("mean_gene_log2FC"))}; n_lib=1+1; n_cells {int(d282_tj.get("n_treat", 0))}/{int(d282_tj.get("n_ctrl", 0))}; cell *p*={fmt_p(d282_tj.get("mwu_p"))}) | {direction(d282_ifn.get("mean_gene_log2FC"))} (mean gene log2FC {fmt_n(d282_ifn.get("mean_gene_log2FC"))}; cell *p*={fmt_p(d282_ifn.get("mwu_p"))}) | {direction(d282_apm.get("mean_gene_log2FC"))} (mean gene log2FC {fmt_n(d282_apm.get("mean_gene_log2FC"))}; cell *p*={fmt_p(d282_apm.get("mwu_p"))}) |
| PC9 | {direction(pc9_tj.get("mean_gene_log2FC"))} (mean gene log2FC {fmt_n(pc9_tj.get("mean_gene_log2FC"))}; n_lib=1+1; n_cells {int(pc9_tj.get("n_treat", 0))}/{int(pc9_tj.get("n_ctrl", 0))}; cell *p*={fmt_p(pc9_tj.get("mwu_p"))}) | {direction(pc9_ifn.get("mean_gene_log2FC"))} (mean gene log2FC {fmt_n(pc9_ifn.get("mean_gene_log2FC"))}; cell *p*={fmt_p(pc9_ifn.get("mwu_p"))}) | {direction(pc9_apm.get("mean_gene_log2FC"))} (mean gene log2FC {fmt_n(pc9_apm.get("mean_gene_log2FC"))}; cell *p*={fmt_p(pc9_apm.get("mwu_p"))}) |

Do not read osi-vs-vehicle as an SKB264 / SG on-treatment result.

## Within-library Spearman (cell-level)

| Library | TACSTD2 vs CLDN4 | TACSTD2 vs TJ | CLDN4 vs IFN/ISG |
|---|---|---|---|
| DFCI282 osi | {corr_line("GSM9128162", "TACSTD2_vs_CLDN4")} | {corr_line("GSM9128162", "TACSTD2_vs_TJ_SIG")} | {corr_line("GSM9128162", "CLDN4_vs_IFN_ISG")} |
| DFCI282 veh | {corr_line("GSM9128163", "TACSTD2_vs_CLDN4")} | {corr_line("GSM9128163", "TACSTD2_vs_TJ_SIG")} | {corr_line("GSM9128163", "CLDN4_vs_IFN_ISG")} |
| PC9 osi | {corr_line("GSM9128164", "TACSTD2_vs_CLDN4")} | {corr_line("GSM9128164", "TACSTD2_vs_TJ_SIG")} | {corr_line("GSM9128164", "CLDN4_vs_IFN_ISG")} |
| PC9 veh | {corr_line("GSM9128165", "TACSTD2_vs_CLDN4")} | {corr_line("GSM9128165", "TACSTD2_vs_TJ_SIG")} | {corr_line("GSM9128165", "CLDN4_vs_IFN_ISG")} |
| pt357 tumor epi | {corr_line("GSM9101265", "TACSTD2_vs_CLDN4")} | {corr_line("GSM9101265", "TACSTD2_vs_TJ_SIG")} | {corr_line("GSM9101265", "CLDN4_vs_IFN_ISG")} |
| pt357 LN epi | {corr_line("GSM9101264", "TACSTD2_vs_CLDN4")} | {corr_line("GSM9101264", "TACSTD2_vs_TJ_SIG")} | {corr_line("GSM9101264", "CLDN4_vs_IFN_ISG")} |

## Patient 357 residual tumor vs LN (descriptive only)

Epithelial cells after neoadjuvant EGFR-TKI. Not a treatment contrast.

| Feature | n_cells tumor / LN | log2FC | direction | cell MWU *p* |
|---|---:|---:|---|---:|
{line(g("patient357_tumor_vs_LN", "TACSTD2"), "TACSTD2")}
{line(g("patient357_tumor_vs_LN", "CLDN4"), "CLDN4")}
{line(s("patient357_tumor_vs_LN", "TJ_SIG"), "TJ signature")}
{line(s("patient357_tumor_vs_LN", "IFN_ISG"), "IFN/ISG")}
{line(s("patient357_tumor_vs_LN", "MHC1_APM"), "MHC-I/APM")}

## Caveats

- **n_library = 1 vs 1.** No sample-level *p* exists. Cell-level *p* will be small whenever thousands of cells shift even slightly; use log2FC / rank-biserial.
- No doublet re-call (authors used scrublet).
- Patient 357 is one person, two sites, both residual on EGFR-TKI.
- mRNA ≠ TROP2 / claudin-4 protein.
- PC9 vehicle GEO annotation inconsistency (PC10 vs PC9).
- Public data only. SRA FASTQ were not re-aligned.

## Reproduce

```bash
python3 methods/gse302284/download.py --cache-dir /tmp/gse302284
python3 methods/gse302284/analyze.py --cache-dir /tmp/gse302284
```

Figures: `results/C_GSE302284/figures/`. Tables: `results/C_GSE302284/tables/`.
"""
    (out / "WRITEUP.md").write_text(md)
    (out / "README.md").write_text(
        "# GSE302284 — TROP2 / CLDN4 / TJ / IFN public context\n\n"
        "See `WRITEUP.md`. Scripts: `methods/gse302284/`.\n"
        "No sacituzumab / SKB264 treated libraries are deposited; "
        "scored contrast is osimertinib vs vehicle.\n"
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache-dir", default="/tmp/gse302284")
    p.add_argument("--out-dir", default=None)
    args = p.parse_args()

    root = Path(__file__).resolve().parents[2]
    out = Path(args.out_dir) if args.out_dir else root / "results" / "C_GSE302284"
    figs = out / "figures"
    tables = out / "tables"
    figs.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    cache_raw = Path(args.cache_dir) / "raw"
    missing = [s["h5"] for s in SAMPLES if not (cache_raw / s["h5"]).exists()]
    if missing:
        raise SystemExit(
            f"Missing H5 in {cache_raw}: {missing}. Run methods/gse302284/download.py first."
        )

    libs = {}
    lib_qc = []
    for rec in SAMPLES:
        lib = load_library(cache_raw, rec)
        libs[rec["gsm"]] = lib
        lib_qc.append(
            {
                "gsm": rec["gsm"],
                "library": rec["library"],
                "model": rec["model"],
                "treatment": rec["treatment"],
                "n_raw": lib["n_raw"],
                "n_qc": lib["n_qc"],
            }
        )
        print(f"{rec['gsm']} {rec['library']}: raw={lib['n_raw']} qc={lib['n_qc']}", flush=True)

    results = [run_contrast(libs, spec) for spec in CONTRASTS]
    gene_df = pd.DataFrame([r for res in results for r in res["gene_rows"]])
    set_df = pd.DataFrame([r for res in results for r in res["set_rows"]])
    corr_df = pd.DataFrame([r for res in results for r in res["corr_rows"]])
    # de-duplicate corr (patient contrast + treatment both emit library corrs)
    corr_df = corr_df.drop_duplicates(subset=["gsm", "pair"]).reset_index(drop=True)

    # BH-FDR within each contrast for focal genes only
    gene_df["mwu_q_focal"] = np.nan
    for cid, sub in gene_df.groupby("contrast"):
        mask = sub["gene"].isin(FOCAL_GENES) & sub["mwu_p"].notna()
        if mask.sum() >= 2:
            q = np.full(len(sub), np.nan)
            _, qvals, _, _ = multipletests(sub.loc[mask, "mwu_p"].to_numpy(), method="fdr_bh")
            q[np.where(mask.to_numpy())[0]] = qvals
            gene_df.loc[sub.index, "mwu_q_focal"] = q

    gene_df.to_csv(tables / "gene_stats.tsv", sep="\t", index=False)
    set_df.to_csv(tables / "geneset_stats.tsv", sep="\t", index=False)
    corr_df.to_csv(tables / "within_library_spearman.tsv", sep="\t", index=False)
    pd.DataFrame(lib_qc).to_csv(tables / "library_qc.tsv", sep="\t", index=False)
    pd.DataFrame(SAMPLES).to_csv(tables / "sample_table.tsv", sep="\t", index=False)

    plot_score_boxes(results, figs / "fig1_osi_vs_veh_scores.png")
    plot_heatmap(gene_df, figs / "fig2_focal_log2fc_heatmap.png")
    plot_scatter(results, figs / "fig3_tacstd2_vs_cldn4.png")
    plot_set_bars(set_df, figs / "fig4_set_direction.png")

    # compact key stats
    key = []
    for res in results:
        cid = res["spec"]["id"]
        for gene in TARGETS:
            row = gene_df[(gene_df["contrast"] == cid) & (gene_df["gene"] == gene)].iloc[0]
            key.append(
                {
                    "contrast": cid,
                    "feature": gene,
                    "n_library": "1 vs 1",
                    "n_cells_treat": int(row["n_treat"]),
                    "n_cells_ctrl": int(row["n_ctrl"]),
                    "log2FC_pseudobulk": row["log2FC_pseudobulk"],
                    "direction": direction(row["log2FC_pseudobulk"]),
                    "mwu_p_cell": row["mwu_p"],
                    "rank_biserial": row["rank_biserial"],
                }
            )
        for name in ("TJ_SIG", "IFN_ISG", "MHC1_APM", "PRIORITY6"):
            row = set_df[(set_df["contrast"] == cid) & (set_df["set"] == name)].iloc[0]
            key.append(
                {
                    "contrast": cid,
                    "feature": name,
                    "n_library": "1 vs 1",
                    "n_cells_treat": int(row["n_treat"]),
                    "n_cells_ctrl": int(row["n_ctrl"]),
                    "log2FC_pseudobulk": row["mean_gene_log2FC"],
                    "direction": direction(row["mean_gene_log2FC"]),
                    "mwu_p_cell": row["mwu_p"],
                    "rank_biserial": row["rank_biserial"],
                }
            )
    key_df = pd.DataFrame(key)
    key_df.to_csv(tables / "key_stats.tsv", sep="\t", index=False)

    summary = {
        "series": "GSE302284",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "sacituzumab_or_skb264_rna_deposited": False,
        "public_analog_of_skb264": (
            "Not as a treatment RNA contrast. SG is mentioned in the series "
            "abstract; deposited libraries are osimertinib vs vehicle plus one "
            "residual patient pair."
        ),
        "n_library_per_arm": 1,
        "libraries": lib_qc,
        "key": key_df.to_dict(orient="records"),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    write_writeup(out, gene_df, set_df, corr_df, results, lib_qc)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
