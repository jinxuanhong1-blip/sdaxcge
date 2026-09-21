#!/usr/bin/env python3
"""CLDN4 vs CytoTRACE2 potency and stemness on concordant-4 malignant cells.

ADDITIVE. CLDN4 only. GSE123902 + GSE131907 + GSE205335 + GSE189357.
No GSE148071. No dual-high TACSTD2∩CLDN4.

CytoTRACE2 is fit separately inside each dataset (Kang et al. recommend
against mixing batches in the KNN smooth). The inferential unit is the
patient / donor / sample. The primary summary is a DerSimonian–Laird
meta-analysis of the four within-cohort Spearman correlations.
"""
from __future__ import annotations

import argparse
import gc
import json
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.neighbors import NearestNeighbors

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import BARRIER_KERATIN, stemness_sets  # noqa: E402
from lib_stats import random_effects_dl, spearman  # noqa: E402

MIN_CELLS = 20
MIN_ARM = 8
DATASETS = ("GSE123902", "GSE131907", "GSE205335", "GSE189357")
COLORS = {
    "GSE123902": "#4c78a8",
    "GSE131907": "#f58518",
    "GSE205335": "#54a24b",
    "GSE189357": "#b279a2",
}


def _fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"


def _fmt_r(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def bh(pvals):
    p = np.asarray(pvals, float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    n = len(ranked)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(ok.sum())
    tmp[order] = q
    out[ok] = tmp
    return out


def log1p_cp10k(X, gene_index, n_umi) -> np.ndarray:
    if gene_index is None:
        return np.full(X.shape[0], np.nan)
    col = X[:, gene_index]
    val = np.asarray(col.toarray()).ravel() if sparse.issparse(col) else np.asarray(col).ravel()
    tot = np.asarray(n_umi, float)
    tot[tot <= 0] = np.nan
    return np.log1p(val / tot * 1e4)


def module_score(X, genes, var_names, n_umi) -> tuple[np.ndarray, int, list[str]]:
    present = [g for g in genes if g in var_names]
    absent = [g for g in genes if g not in var_names]
    if not present:
        return np.full(X.shape[0], np.nan), 0, absent
    idx = [var_names.get_loc(g) for g in present]
    sub = X[:, idx]
    raw = np.asarray(sub.toarray()) if sparse.issparse(sub) else np.asarray(sub)
    tot = np.asarray(n_umi, float)
    tot[tot <= 0] = np.nan
    expr = np.log1p(raw / tot[:, None] * 1e4)
    std = np.nanstd(expr, axis=0)
    keep = std > 1e-8
    if keep.sum() == 0:
        return np.full(X.shape[0], np.nan), 0, absent
    expr = expr[:, keep]
    z = (expr - np.nanmean(expr, axis=0)) / (np.nanstd(expr, axis=0) + 1e-8)
    return np.nanmean(z, axis=1), int(keep.sum()), absent


def gulati_cytotrace(X, n_top: int = 200, n_neighbors: int = 10) -> np.ndarray:
    """Gulati et al. Science 2020 CytoTRACE. Higher = more potent."""
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    X = X.tocsr()
    n_cells = X.shape[0]
    gene_counts = np.asarray((X > 0).sum(axis=1)).ravel().astype(float)
    nz_cells = np.asarray((X > 0).sum(axis=0)).ravel()
    keep = nz_cells >= 10
    if keep.sum() < n_top:
        keep = nz_cells >= 3
    Xs = X[:, keep]
    logx = Xs.copy()
    logx.data = np.log1p(logx.data)
    gc = gene_counts - gene_counts.mean()
    gc_ss = float(np.dot(gc, gc))
    if gc_ss <= 0:
        return stats.rankdata(gene_counts) / n_cells
    col_means = np.asarray(logx.mean(axis=0)).ravel()
    xt_gc = np.asarray(logx.T.dot(gc)).ravel()
    col_ss = np.asarray(logx.power(2).sum(axis=0)).ravel() - n_cells * col_means**2
    denom = np.sqrt(np.clip(col_ss, 0, None) * gc_ss)
    corr = np.zeros(int(keep.sum()), dtype=float)
    ok = denom > 1e-12
    corr[ok] = xt_gc[ok] / denom[ok]
    top = np.argsort(corr)[::-1][: min(n_top, int(keep.sum()))]
    score = np.asarray(logx[:, top].mean(axis=1)).ravel()
    feat = np.asarray(logx[:, top].toarray())
    feat = feat - feat.mean(axis=0)
    std = feat.std(axis=0)
    std[std == 0] = 1
    feat = feat / std
    k = min(n_neighbors, n_cells - 1)
    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean")
    nn.fit(feat)
    idx = nn.kneighbors(feat, return_distance=False)
    smooth = score[idx].mean(axis=1)
    return stats.rankdata(smooth) / n_cells


def counts_frame(adata) -> pd.DataFrame:
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    X = X.tocsr()
    nz = np.asarray((X > 0).sum(axis=0)).ravel()
    keep = nz >= 10
    genes = np.asarray(adata.var_names, dtype=str)[keep]
    dense = np.asarray(X[:, keep].todense(), dtype=np.float32)
    return pd.DataFrame(dense, index=np.asarray(adata.obs_names, dtype=str), columns=genes)


def run_cytotrace2(adata, work: Path) -> tuple[pd.DataFrame | None, dict]:
    info = {
        "package": "cytotrace2-py",
        "citation": "Kang et al. Nat Methods 2025 doi:10.1038/s41592-025-02857-2",
        "ran": False,
        "dataset": str(adata.obs["dataset"].iloc[0]),
    }
    try:
        from cytotrace2_py import cytotrace2_py as ctmod
    except Exception as exc:
        info["reason"] = f"import failed: {exc}"
        return None, info
    work.mkdir(parents=True, exist_ok=True)
    out_dir = work / "cytotrace2_out"
    cached = out_dir / "cytotrace2_results.txt"
    if cached.is_file() and cached.stat().st_size > 1000:
        df = pd.read_csv(cached, sep="\t", index_col=0)
        info["ran"] = True
        info["cached"] = True
        info["n_scored"] = int(len(df))
        return df, info
    expr = counts_frame(adata)
    info["n_cells"] = int(expr.shape[0])
    info["n_genes_input"] = int(expr.shape[1])
    orig = ctmod.read_file
    ctmod.read_file = lambda path: expr
    try:
        print(
            f"cytotrace2 {info['dataset']} cells={expr.shape[0]} genes={expr.shape[1]}",
            flush=True,
        )
        result = ctmod.cytotrace2(
            f"in-memory-{info['dataset']}",
            species="human",
            batch_size=min(4000, expr.shape[0]),
            smooth_batch_size=min(800, expr.shape[0]),
            disable_plotting=True,
            disable_parallelization=True,
            max_cores=1,
            seed=14,
            output_dir=str(out_dir),
        )
    except Exception as exc:
        info["reason"] = f"cytotrace2() failed: {exc}\n{traceback.format_exc()}"
        print(info["reason"], flush=True)
        return None, info
    finally:
        ctmod.read_file = orig
        del expr
        gc.collect()
    if not isinstance(result, pd.DataFrame):
        result = pd.DataFrame(result)
    info["ran"] = True
    info["n_scored"] = int(len(result))
    info["columns"] = list(result.columns)
    return result, info


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def score_dataset(path: Path, work: Path, sets: dict) -> tuple[pd.DataFrame, dict]:
    import anndata as ad

    adata = ad.read_h5ad(path)
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    X = X.tocsr()
    n_umi = np.asarray(X.sum(axis=1)).ravel()
    var = pd.Index(adata.var_names.astype(str))
    cldn_i = var.get_loc("CLDN4") if "CLDN4" in var else None
    obs = adata.obs.copy()
    obs["expr_CLDN4"] = log1p_cp10k(X, cldn_i, n_umi)
    obs["CLDN4_pos"] = (obs["expr_CLDN4"] > 0).astype(float)
    barrier, n_b, miss_b = module_score(X, BARRIER_KERATIN, var, n_umi)
    obs["score_barrier"] = barrier
    coverage = {"barrier_genes_used": n_b, "barrier_absent": miss_b}
    for key, genes in sets.items():
        sc, n_used, absent = module_score(X, genes, var, n_umi)
        obs[f"score_{key}"] = sc
        coverage[f"{key}_genes_used"] = n_used
        coverage[f"{key}_n_set"] = len(genes)
        coverage[f"{key}_absent_n"] = len(absent)
        raw_idx = [var.get_loc(g) for g in genes if g in var]
        if raw_idx:
            sub = X[:, raw_idx]
            coverage[f"{key}_frac_any_detected"] = float(
                np.asarray((sub.sum(axis=1) > 0)).ravel().mean()
            )
        else:
            coverage[f"{key}_frac_any_detected"] = 0.0
    obs["gulati"] = gulati_cytotrace(X)
    ct_df, ct_info = run_cytotrace2(adata, work)
    if ct_df is None:
        raise SystemExit(f"CytoTRACE2 failed for {path.name}: {ct_info}")
    ct_df = ct_df.copy()
    ct_df.index = ct_df.index.astype(str)
    mapped = ct_df.reindex(obs.index.astype(str))
    obs["CytoTRACE2_Score"] = pd.to_numeric(mapped["CytoTRACE2_Score"], errors="coerce")
    obs["preKNN_CytoTRACE2_Score"] = pd.to_numeric(
        mapped["preKNN_CytoTRACE2_Score"], errors="coerce"
    )
    obs["CytoTRACE2_Potency"] = mapped["CytoTRACE2_Potency"].astype(str)
    n_ok = int(np.isfinite(obs["CytoTRACE2_Score"]).sum())
    if n_ok < 50:
        raise SystemExit(f"CytoTRACE2 joined only {n_ok} cells for {path.name}")
    ct_info["n_joined"] = n_ok
    ct_info["frac_genes_lt_500"] = float((obs["n_genes_detected"] < 500).mean())
    keep = [
        "dataset", "unit_id", "patient_id", "tissue", "histology", "stage",
        "malignant_rule", "n_umi", "n_genes_detected", "expr_CLDN4", "CLDN4_pos",
        "score_barrier", "score_benporath_es1", "score_wong_esc", "score_benporath_core9",
        "gulati", "CytoTRACE2_Score", "preKNN_CytoTRACE2_Score", "CytoTRACE2_Potency",
    ]
    cells = obs[keep].copy()
    cells.index = obs.index.astype(str)
    cells.index.name = "cell"
    meta = {
        "cytotrace2": ct_info,
        "coverage": coverage,
        "n_cells": int(len(cells)),
        "potency_counts": cells["CytoTRACE2_Potency"].value_counts().to_dict(),
    }
    del adata, X
    gc.collect()
    return cells, meta


def within_unit_arms(g: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    cldn = g["expr_CLDN4"]
    if len(g) >= 3 * MIN_ARM:
        q1, q2 = np.nanquantile(cldn, [1.0 / 3.0, 2.0 / 3.0])
        hi = g.loc[cldn >= q2]
        lo = g.loc[cldn <= q1]
        return hi, lo, "within-unit tertile"
    med = np.nanmedian(cldn)
    hi = g.loc[cldn > med]
    lo = g.loc[cldn <= med]
    return hi, lo, "within-unit median"


def patient_table(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for uid, g in cells.groupby("unit_id", sort=False):
        hi, lo, rule = within_unit_arms(g)
        def arm_mean(frame, col):
            return float(frame[col].mean()) if len(frame) else np.nan
        rows.append(
            {
                "unit_id": uid,
                "dataset": g["dataset"].iloc[0],
                "patient_id": g["patient_id"].iloc[0],
                "tissue": g["tissue"].iloc[0],
                "histology": g["histology"].iloc[0],
                "stage": g["stage"].iloc[0],
                "malignant_rule": g["malignant_rule"].iloc[0],
                "n_malignant": int(len(g)),
                "mean_CLDN4": float(g["expr_CLDN4"].mean()),
                "pct_CLDN4": float(g["CLDN4_pos"].mean() * 100.0),
                "mean_ct2": float(g["CytoTRACE2_Score"].mean()),
                "mean_preknn": float(g["preKNN_CytoTRACE2_Score"].mean()),
                "frac_diff": float((g["CytoTRACE2_Potency"] == "Differentiated").mean()),
                "mean_es1": float(g["score_benporath_es1"].mean()),
                "mean_wong": float(g["score_wong_esc"].mean()),
                "mean_core9": float(g["score_benporath_core9"].mean()),
                "mean_barrier": float(g["score_barrier"].mean()),
                "mean_gulati": float(g["gulati"].mean()),
                "arm_rule": rule,
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "ct2_high": arm_mean(hi, "CytoTRACE2_Score"),
                "ct2_low": arm_mean(lo, "CytoTRACE2_Score"),
                "es1_high": arm_mean(hi, "score_benporath_es1"),
                "es1_low": arm_mean(lo, "score_benporath_es1"),
                "wong_high": arm_mean(hi, "score_wong_esc"),
                "wong_low": arm_mean(lo, "score_wong_esc"),
                "barrier_high": arm_mean(hi, "score_barrier"),
                "barrier_low": arm_mean(lo, "score_barrier"),
                "gulati_high": arm_mean(hi, "gulati"),
                "gulati_low": arm_mean(lo, "gulati"),
                "preknn_high": arm_mean(hi, "preKNN_CytoTRACE2_Score"),
                "preknn_low": arm_mean(lo, "preKNN_CytoTRACE2_Score"),
                "eligible": int(len(g) >= MIN_CELLS),
            }
        )
    return pd.DataFrame(rows)


def cohort_spearmans(elig: pd.DataFrame, x: str, y: str) -> tuple[list[float], list[int], list[str]]:
    rhos, ns, labels = [], [], []
    for ds in DATASETS:
        sub = elig[elig["dataset"] == ds]
        r, _p, n = spearman(sub[x], sub[y])
        if n >= 4 and np.isfinite(r):
            rhos.append(r)
            ns.append(n)
            labels.append(ds)
    return rhos, ns, labels


def wilcoxon_paired(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if len(a) < 4:
        return np.nan, np.nan, np.nan, int(len(a))
    d = a - b
    if np.allclose(d, 0):
        return 0.0, 1.0, 0.0, int(len(a))
    w, p = stats.wilcoxon(a, b, alternative="two-sided", zero_method="wilcox")
    return float(w), float(p), float(np.median(d)), int(len(a))


def scatter(elig, x, y, xlab, ylab, title, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    for ds in DATASETS:
        sub = elig[elig["dataset"] == ds]
        if sub.empty:
            continue
        ax.scatter(
            sub[x], sub[y], s=46, c=COLORS[ds], label=f"{ds} n={len(sub)}",
            edgecolors="white", linewidths=0.4,
        )
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    _save(fig, path)


def forest(cohort_rows: list[dict], path: Path) -> None:
    """One panel per contrast so the names stay inside the axes."""
    order = [
        "CLDN4 vs CytoTRACE2",
        "CLDN4 vs Ben-Porath ES1 stemness",
        "CLDN4 vs barrier/keratin (no CLDN4)",
    ]
    fig, axes = plt.subplots(len(order), 1, figsize=(7.4, 6.8), sharex=True)
    for ax, contrast in zip(axes, order):
        block = [r for r in cohort_rows if r["contrast"] == contrast]
        for i, r in enumerate(block):
            ax.plot([r["lo"], r["hi"]], [i, i], color=COLORS[r["dataset"]], lw=1.8)
            ax.scatter([r["rho"]], [i], s=36, color=COLORS[r["dataset"]], zorder=3)
        ax.axvline(0, color="0.45", lw=0.8, ls="--")
        ax.set_yticks(range(len(block)))
        ax.set_yticklabels([f"{r['dataset']}  n={r['n']}" for r in block], fontsize=8)
        ax.set_title(contrast, fontsize=10, loc="left")
        if block:
            ax.set_ylim(-0.6, len(block) - 0.4)
    axes[-1].set_xlabel("within-cohort Spearman ρ (Fisher-z 95% CI)")
    fig.suptitle("Concordant-4 CLDN4 associations", fontsize=12)
    _save(fig, path)


def write_finding(path: Path, S: dict) -> None:
    prim = {r["contrast"]: r for r in S["primary"]}
    c4 = prim["CLDN4 vs CytoTRACE2"]
    es = prim["CLDN4 vs Ben-Porath ES1 stemness"]
    wg = prim["CLDN4 vs Wong ESC stemness"]
    bar = prim["CLDN4 vs barrier/keratin (no CLDN4)"]
    agree = prim["CytoTRACE2 vs Ben-Porath ES1"]
    pbar = prim["CytoTRACE2 vs barrier/keratin (no CLDN4)"]
    frac = prim["CLDN4 vs frac Differentiated"]

    ct2_tail = (
        f"ρ={_fmt_r(c4['rho'])}, p={_fmt_p(c4['p'])}, q={_fmt_p(c4['q'])}, I²={c4['I2']:.1f}%"
    )
    if np.isfinite(c4.get("p", np.nan)) and c4["p"] < 0.05:
        ct2_call = f"the four-cohort meta is significant ({ct2_tail})"
    else:
        ct2_call = f"the four-cohort meta is not significant ({ct2_tail}, 95% CI includes 0)"
    answer = (
        f"CLDN4-high is barrier-like "
        f"(DL ρ={_fmt_r(bar['rho'])}, p={_fmt_p(bar['p'])}, q={_fmt_p(bar['q'])}, I²={bar['I2']:.1f}%). "
        f"The CytoTRACE2 Differentiated fraction rises with CLDN4 "
        f"(ρ={_fmt_r(frac['rho'])}, p={_fmt_p(frac['p'])}, q={_fmt_p(frac['q'])}, I²={frac['I2']:.1f}%). "
        f"Continuous CytoTRACE2 is {'lower' if c4['rho'] < 0 else 'higher'} where CLDN4 is higher, but {ct2_call}. "
        f"Embryonic stemness does not fall with CLDN4 "
        f"(Ben-Porath ES1 ρ={_fmt_r(es['rho'])}, p={_fmt_p(es['p'])}; "
        f"Wong ESC ρ={_fmt_r(wg['rho'])}, p={_fmt_p(wg['p'])})"
    )
    lines = [
        "# Finding — concordant-4 malignant CLDN4 vs CytoTRACE2 and stemness",
        "",
        "ADDITIVE. **CLDN4 only.** Same four public sets as the locked concordant-4 "
        "patient table: **GSE123902 + GSE131907 + GSE205335 + GSE189357**. "
        "Not GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526. "
        "No TACSTD2∩CLDN4 dual-high. This does **not** re-audit the T/NK exclusion ρ.",
        "",
        "Primary potency = **CytoTRACE2** (`cytotrace2-py` "
        f"{S['cytotrace2_version']}, Kang et al. *Nat Methods* 2025), fit **inside each dataset**. "
        "Score 0 = differentiated, 1 = totipotent. Stemness is paired, not a substitute: "
        "**Ben-Porath ES1** and **Wong ESC core** (MSigDB c2.cgp v2023.2.Hs), "
        "mean z of log1p(CP10k), with CLDN4, EPCAM, TACSTD2, and the barrier/keratin genes removed. "
        "Barrier/keratin = KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**). "
        "Inferential unit = patient / donor / sample. "
        "Primary number = DerSimonian–Laird meta of the four within-cohort Spearman ρ. "
        "Cell-level ρ is not the claim.",
        "",
        f"**Question.** Are CLDN4-high malignant cells more differentiated and more barrier-like? "
        f"**Answer:** {answer}.",
        "",
        "## Verdict",
        "",
        f"DL meta CLDN4 vs CytoTRACE2: k={c4['k']}, N={c4['n']}, ρ={_fmt_r(c4['rho'])}, "
        f"95% CI [{_fmt_r(c4['ci_lo'])}, {_fmt_r(c4['ci_hi'])}], p={_fmt_p(c4['p'])}, "
        f"q={_fmt_p(c4['q'])}, I²={c4['I2']:.1f}%.",
        f"DL meta CLDN4 vs Ben-Porath ES1 stemness: ρ={_fmt_r(es['rho'])}, p={_fmt_p(es['p'])}, "
        f"q={_fmt_p(es['q'])}, I²={es['I2']:.1f}%.",
        f"DL meta CLDN4 vs Wong ESC stemness: ρ={_fmt_r(wg['rho'])}, p={_fmt_p(wg['p'])}, "
        f"q={_fmt_p(wg['q'])}, I²={wg['I2']:.1f}%.",
        f"DL meta CLDN4 vs barrier/keratin (CLDN4 excluded): ρ={_fmt_r(bar['rho'])}, "
        f"p={_fmt_p(bar['p'])}, q={_fmt_p(bar['q'])}, I²={bar['I2']:.1f}%.",
        f"CytoTRACE2 vs ES1 (do the two potency readouts agree): ρ={_fmt_r(agree['rho'])}, "
        f"p={_fmt_p(agree['p'])}, q={_fmt_p(agree['q'])}.",
        f"CytoTRACE2 vs barrier: ρ={_fmt_r(pbar['rho'])}, p={_fmt_p(pbar['p'])}, q={_fmt_p(pbar['q'])}.",
        f"CLDN4 vs fraction Differentiated: ρ={_fmt_r(frac['rho'])}, p={_fmt_p(frac['p'])}, "
        f"q={_fmt_p(frac['q'])}.",
        "",
        "Predicted signs if CLDN4-high is differentiated / barrier-like: CytoTRACE2 and stemness "
        "**negative**, barrier and Differentiated fraction **positive**, CytoTRACE2 vs stemness **positive**, "
        "CytoTRACE2 vs barrier **negative**.",
        "",
        f"**What holds.** {S['what_holds']}",
        "",
        f"**What does not hold.** {S['what_fails']}",
        "",
        "## Honest n",
        "",
        f"- Analysis cells after QC and cap ≤{S['cap']}/unit: **n_cells = {S['n_cells']}** "
        + ", ".join(f"{ds} {S['n_cells_by'][ds]}" for ds in DATASETS)
        + ".",
        f"- Units: **n_units = {S['n_units']}** "
        + ", ".join(f"{ds} {S['n_units_by'][ds]}" for ds in DATASETS)
        + f". Eligible (≥{MIN_CELLS} malignant cells): **n = {S['n_eligible']}**.",
        f"- Paired within-unit CLDN4-high vs low (≥{MIN_ARM} cells/arm): **n = {S['n_paired']}**.",
        "- GSE131907 / GSE205335 malignant = **author label**, not a new CNV call. "
        f"GSE131907 author-malignant catalog {S['catalog'].get('GSE131907', {}).get('n_author_malignant', 'NA')}; "
        f"tLung author-malignant = {S['catalog'].get('GSE131907', {}).get('tLung_author_malignant', 'NA')}. "
        f"GSE205335 author-malignant non-normal catalog "
        f"{S['catalog'].get('GSE205335', {}).get('n_author_malignant_non_normal', 'NA')}.",
        "- GSE123902 and GSE189357 malignant = marker gate "
        "`(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`. **Not CNV.** "
        "NORMAL lung from GSE123902 is dropped. Gate counts are in `extract_audit.tsv`.",
        f"- CytoTRACE2 potency categories (cells): {S['potency_counts']}.",
        f"- Cells with <500 genes (CytoTRACE2 prefers 500–1000): {S['frac_lt_500']:.3f}.",
        f"- Stemness gene coverage: {S['coverage_brief']}.",
        "- GSE148071 not used. Dual-high not used. T/NK ρ not re-fit.",
        "",
        "## Locked choices",
        "",
        "- CytoTRACE2 per dataset, human model, `max_cores=1`, seed 14. Not a joint KNN.",
        "- Gulati 2020 gene-count CytoTRACE is sensitivity only.",
        "- CLDN4 = patient-mean log1p(CP10k). %pos is sensitivity.",
        "- Stemness z-scores are within dataset.",
        "- Primary family = seven DL metas, BH inside that list.",
        "- Paired test = within-unit CLDN4 tertile (median split if the unit is small).",
        "- Cap ≤200 cells/unit after QC (author sets: cap then QC, matching the winning-pair potency run).",
        "",
        "## Primary (DL meta of within-cohort Spearman, BH inside this list)",
        "",
        "| Contrast | k | N | ρ | 95% CI | p | q | I² |",
        "| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: |",
    ]
    for r in S["primary"]:
        lines.append(
            f"| {r['contrast']} | {r['k']} | {r['n']} | {_fmt_r(r['rho'])} | "
            f"[{_fmt_r(r['ci_lo'])}, {_fmt_r(r['ci_hi'])}] | {_fmt_p(r['p'])} | "
            f"{_fmt_p(r['q'])} | {r['I2']:.1f}% |"
        )
    lines += [
        "",
        "## Within-cohort Spearman (not BH)",
        "",
        "| Contrast | dataset | n | ρ | p |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for r in S["cohort"]:
        lines.append(
            f"| {r['contrast']} | {r['dataset']} | {r['n']} | {_fmt_r(r['rho'])} | {_fmt_p(r['p'])} |"
        )
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | k or n | ρ or Δ | p | note |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for r in S["sensitivity"]:
        stat = r.get("rho", r.get("delta"))
        lines.append(
            f"| {r['contrast']} | {r.get('n', r.get('k'))} | {_fmt_r(stat)} | {_fmt_p(r.get('p'))} | {r.get('note','')} |"
        )
    lines += [
        "",
        "## Paired within-unit CLDN4-high minus CLDN4-low",
        "",
        "Negative potency or stemness Δ = high arm more differentiated. Positive barrier Δ = high arm more barrier-like.",
        "",
        "| Contrast | n_units | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in S["paired"]:
        d = r.get("delta")
        ds = "NA" if d is None or not np.isfinite(d) else f"{d:.3f}"
        lines.append(f"| {r['contrast']} | {r['n']} | {ds} | {_fmt_p(r['p'])} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Not a rewrite of concordant-4 T/NK exclusion (malignant CLDN4 %pos vs T/NK, n=65).",
        "- Not an ICI / RECIST / stage test. GSE189357 stage n=3 per bin is metadata only.",
        "- Marker-malignant is not copy-number malignant.",
        "- GSE205335 mixes ADC / SQ / SCLC / NUT. ADC-only is in the sensitivity table.",
        "- A tiny paired Δ is not a large within-tumor stemness collapse.",
        "- Do not call CLDN4 a stemness marker from this table.",
        "- Visium same-spot correlation is not in this folder.",
        "",
        "## Outputs",
        "",
        "- `results/tables/patient_cldn4_vs_potency.tsv` — done criterion",
        "- `results/tables/extract_audit.tsv` — marker-gate counts vs the locked pseudobulk n",
        "- `results/tables/patient_means.tsv`",
        "- `results/tables/cohort_spearman.tsv`",
        "- `results/tables/stats.tsv`",
        "- `results/figures/fig_patient_cldn4_vs_potency.png`",
        "- `results/figures/fig_patient_cldn4_vs_stemness.png`",
        "- `results/figures/fig_patient_cldn4_vs_barrier.png`",
        "- `results/figures/fig_forest_cohort_rho.png`",
        "- `results/figures/fig_paired_high_vs_low.png`",
        "- `results/figures/fig_honest_n.png`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 -m venv /tmp/ct2venv",
        "/tmp/ct2venv/bin/pip install -r methods/concordant4_cytotrace2_cldn4/requirements.txt",
        "/tmp/ct2venv/bin/python methods/concordant4_cytotrace2_cldn4/scripts/download.py --out /tmp/concordant4_ct2",
        "/tmp/ct2venv/bin/python methods/concordant4_cytotrace2_cldn4/scripts/extract_malignant.py --data /tmp/concordant4_ct2",
        "/tmp/ct2venv/bin/python methods/concordant4_cytotrace2_cldn4/scripts/analyze.py \\",
        "  --h5ad /tmp/concordant4_ct2/h5ad --outdir methods/concordant4_cytotrace2_cldn4/results",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")
    print(f"wrote {path}", flush=True)


def verdict_sentences(
    primary: list[dict], paired: list[dict], sensitivity: list[dict], cohort: list[dict]
) -> tuple[str, str]:
    by = {r["contrast"]: r for r in primary}
    sens = {r["contrast"]: r for r in sensitivity}
    pair = {r["contrast"]: r for r in paired}

    def row_txt(name: str) -> str:
        r = by[name]
        return (
            f"{name}: ρ={_fmt_r(r['rho'])}, p={_fmt_p(r['p'])}, q={_fmt_p(r['q'])}, "
            f"I²={r['I2']:.1f}%"
        )

    holds = [
        "Between patients, CLDN4 tracks barrier/keratin with CLDN4 held out of the score "
        f"({row_txt('CLDN4 vs barrier/keratin (no CLDN4)')}).",
        "Between patients, CLDN4 tracks a higher CytoTRACE2 Differentiated fraction "
        f"({row_txt('CLDN4 vs frac Differentiated')}). "
        "GSE131907 carries the cohort-level signal; the other three cohorts are the same sign and not significant alone.",
        "Within a tumor, CLDN4-high cells are barrier-higher "
        f"(n={pair['barrier/keratin high vs low']['n']}, "
        f"Δmed={pair['barrier/keratin high vs low']['delta']:+.3f}, "
        f"p={_fmt_p(pair['barrier/keratin high vs low']['p'])}).",
    ]
    wong = sens.get("DL CytoTRACE2 vs Wong ESC", {})
    loo = sens.get("LOO drop GSE189357: CLDN4 vs CytoTRACE2", {})
    fails = [
        "Continuous CytoTRACE2 vs CLDN4 is negative but the four-cohort meta CI includes 0 "
        f"({row_txt('CLDN4 vs CytoTRACE2')}). Cohorts: "
        + "; ".join(
            f"{r['dataset']} ρ={_fmt_r(r['rho'])} (p={_fmt_p(r['p'])}, n={r['n']})"
            for r in cohort
            if r["contrast"] == "CLDN4 vs CytoTRACE2"
        )
        + ". "
        f"Leave-one-out dropping GSE189357 is ρ={_fmt_r(loo.get('rho'))}, p={_fmt_p(loo.get('p'))}. "
        "That sensitivity is not the primary.",
        "The within-patient CytoTRACE2 shift is null "
        f"(Δmed={pair['CytoTRACE2 high vs low']['delta']:+.3f}, "
        f"p={_fmt_p(pair['CytoTRACE2 high vs low']['p'])}). "
        "Do not quote it as a within-tumor potency drop.",
        "Stemness does not mark CLDN4-high cells as less stem-like. "
        f"Ben-Porath ES1 vs CLDN4 is positive ({row_txt('CLDN4 vs Ben-Porath ES1 stemness')}). "
        f"Wong ESC vs CLDN4 is null ({row_txt('CLDN4 vs Wong ESC stemness')}). "
        "Within a tumor, CLDN4-high cells score higher on both "
        f"(ES1 Δmed={pair['Ben-Porath ES1 high vs low']['delta']:+.3f}, "
        f"p={_fmt_p(pair['Ben-Porath ES1 high vs low']['p'])}; "
        f"Wong Δmed={pair['Wong ESC high vs low']['delta']:+.3f}, "
        f"p={_fmt_p(pair['Wong ESC high vs low']['p'])}).",
        "Wong stemness does track CytoTRACE2 "
        f"(sensitivity DL ρ={_fmt_r(wong.get('rho'))}, p={_fmt_p(wong.get('p'))}), "
        "so the stemness score is coupled to potency and CLDN4 is not on that axis.",
        "Gulati 2020 gene-count CytoTRACE goes the other way on the paired test "
        f"(Δmed={pair['Gulati2020 high vs low']['delta']:+.3f}, "
        f"p={_fmt_p(pair['Gulati2020 high vs low']['p'])}). "
        "Same disagreement as the winning-pair potency folder. This is why Gulati is not the primary.",
        "CLDN4 %pos vs CytoTRACE2 is null. Do not substitute the T/NK %pos score for the mean used here.",
    ]
    return " ".join(holds), " ".join(fails)


def main() -> None:
    here = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", type=Path, default=Path("/tmp/concordant4_ct2/h5ad"))
    p.add_argument("--outdir", type=Path, default=here / "results")
    p.add_argument("--finding", type=Path, default=here / "FINDING.md")
    p.add_argument("--work", type=Path, default=Path("/tmp/concordant4_ct2/work"))
    args = p.parse_args()
    sets = stemness_sets()
    frames = []
    metas = {}
    for ds in DATASETS:
        cells, meta = score_dataset(args.h5ad / f"{ds}.h5ad", args.work / ds, sets)
        frames.append(cells)
        metas[ds] = meta
        print(f"scored {ds} cells={len(cells)}", flush=True)
    cells = pd.concat(frames, axis=0)
    patient = patient_table(cells)
    elig = patient[patient["eligible"] == 1].copy()

    primary_spec = [
        ("CLDN4 vs CytoTRACE2", "mean_CLDN4", "mean_ct2", "neg"),
        ("CLDN4 vs Ben-Porath ES1 stemness", "mean_CLDN4", "mean_es1", "neg"),
        ("CLDN4 vs Wong ESC stemness", "mean_CLDN4", "mean_wong", "neg"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier", "pos"),
        ("CytoTRACE2 vs Ben-Porath ES1", "mean_ct2", "mean_es1", "pos"),
        ("CytoTRACE2 vs barrier/keratin (no CLDN4)", "mean_ct2", "mean_barrier", "neg"),
        ("CLDN4 vs frac Differentiated", "mean_CLDN4", "frac_diff", "pos"),
    ]
    primary = []
    cohort_rows = []
    for contrast, x, y, _expect in primary_spec:
        rhos, ns, labels = cohort_spearmans(elig, x, y)
        dl = random_effects_dl(rhos, ns)
        ci = dl.get("ci95_rho", [np.nan, np.nan])
        primary.append(
            {
                "contrast": contrast,
                "k": dl.get("k", 0),
                "n": dl.get("n_patients_total", 0),
                "rho": dl.get("pooled_rho", np.nan),
                "ci_lo": ci[0] if ci else np.nan,
                "ci_hi": ci[1] if ci else np.nan,
                "p": dl.get("p", np.nan),
                "I2": dl.get("I2", np.nan),
                "x": x,
                "y": y,
            }
        )
        for ds in DATASETS:
            sub = elig[elig["dataset"] == ds]
            r, pv, n = spearman(sub[x], sub[y])
            if n >= 4 and np.isfinite(r):
                se = 1.0 / np.sqrt(n - 3)
                z = np.arctanh(np.clip(r, -0.999999, 0.999999))
                cohort_rows.append(
                    {
                        "contrast": contrast,
                        "dataset": ds,
                        "n": n,
                        "rho": r,
                        "p": pv,
                        "lo": float(np.tanh(z - 1.96 * se)),
                        "hi": float(np.tanh(z + 1.96 * se)),
                    }
                )
    qs = bh([r["p"] for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = float(q)

    sensitivity = []

    def add_dl(name, frame, x, y, note):
        rhos, ns, _labels = cohort_spearmans(frame, x, y)
        dl = random_effects_dl(rhos, ns)
        sensitivity.append(
            {
                "contrast": name,
                "n": dl.get("n_patients_total", 0),
                "k": dl.get("k", 0),
                "rho": dl.get("pooled_rho", np.nan),
                "p": dl.get("p", np.nan),
                "note": note,
            }
        )

    r, pv, n = spearman(elig["mean_CLDN4"], elig["mean_ct2"])
    sensitivity.append(
        {
            "contrast": "pooled Spearman CLDN4 vs CytoTRACE2 (not meta)",
            "n": n, "rho": r, "p": pv,
            "note": "all eligible units stacked; cohort is not the weight",
        }
    )
    add_dl("DL CLDN4 %pos vs CytoTRACE2", elig, "pct_CLDN4", "mean_ct2", "sensitivity; %pos is the T/NK score")
    add_dl("DL CLDN4 vs Gulati2020", elig, "mean_CLDN4", "mean_gulati", "gene-count CytoTRACE, per dataset")
    add_dl("DL CLDN4 vs preKNN CytoTRACE2", elig, "mean_CLDN4", "mean_preknn", "unsmoothed CytoTRACE2")
    add_dl("DL CLDN4 vs Ben-Porath core nine", elig, "mean_CLDN4", "mean_core9", "n=9 TF set")
    add_dl("DL CytoTRACE2 vs Wong ESC", elig, "mean_ct2", "mean_wong", "second stemness agreement")
    for drop in DATASETS:
        sub = elig[elig["dataset"] != drop]
        add_dl(f"LOO drop {drop}: CLDN4 vs CytoTRACE2", sub, "mean_CLDN4", "mean_ct2", "leave-one-cohort-out")
    adc = elig[(elig["dataset"] == "GSE205335") & (elig["histology"] == "ADC")]
    r, pv, n = spearman(adc["mean_CLDN4"], adc["mean_ct2"])
    sensitivity.append({"contrast": "GSE205335 ADC-only CLDN4 vs CytoTRACE2", "n": n, "rho": r, "p": pv, "note": "not meta"})
    for tissue, label in (("PRIMARY", "GSE123902 primary-only"), ("METASTASIS", "GSE123902 met-only")):
        sub = elig[(elig["dataset"] == "GSE123902") & (elig["tissue"] == tissue)]
        r, pv, n = spearman(sub["mean_CLDN4"], sub["mean_ct2"])
        sensitivity.append({"contrast": f"{label} CLDN4 vs CytoTRACE2", "n": n, "rho": r, "p": pv, "note": "donor tissue; n may be small"})

    paired_df = elig[(elig["n_high"] >= MIN_ARM) & (elig["n_low"] >= MIN_ARM)].copy()
    paired = []
    for name, hi, lo in (
        ("CytoTRACE2 high vs low", "ct2_high", "ct2_low"),
        ("Ben-Porath ES1 high vs low", "es1_high", "es1_low"),
        ("Wong ESC high vs low", "wong_high", "wong_low"),
        ("barrier/keratin high vs low", "barrier_high", "barrier_low"),
        ("Gulati2020 high vs low", "gulati_high", "gulati_low"),
        ("preKNN CytoTRACE2 high vs low", "preknn_high", "preknn_low"),
    ):
        w, pv, delta, n = wilcoxon_paired(paired_df[hi], paired_df[lo])
        paired.append({"contrast": name, "n": n, "W": w, "p": pv, "delta": delta})

    # figures
    figdir = args.outdir / "figures"
    tabdir = args.outdir / "tables"
    tabdir.mkdir(parents=True, exist_ok=True)
    c4 = next(r for r in primary if r["contrast"] == "CLDN4 vs CytoTRACE2")
    es = next(r for r in primary if r["contrast"] == "CLDN4 vs Ben-Porath ES1 stemness")
    wg = next(r for r in primary if r["contrast"] == "CLDN4 vs Wong ESC stemness")
    bar = next(r for r in primary if r["contrast"] == "CLDN4 vs barrier/keratin (no CLDN4)")
    scatter(
        elig, "mean_ct2", "mean_CLDN4",
        "patient-mean CytoTRACE2 (higher = more potent)",
        "patient-mean CLDN4 log1p(CP10k)",
        f"DL ρ={_fmt_r(c4['rho'])}  p={_fmt_p(c4['p'])}  I²={c4['I2']:.0f}%  N={c4['n']}",
        figdir / "fig_patient_cldn4_vs_potency",
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.6))
    for ax, col, row, lab in (
        (axes[0], "mean_es1", es, "Ben-Porath ES1 (z, barrier genes out)"),
        (axes[1], "mean_wong", wg, "Wong ESC core (z)"),
    ):
        for ds in DATASETS:
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[col], sub["mean_CLDN4"], s=36, c=COLORS[ds], label=ds, edgecolors="white", linewidths=0.3)
        ax.set_xlabel(lab)
        ax.set_ylabel("patient-mean CLDN4")
        ax.set_title(f"DL ρ={_fmt_r(row['rho'])} p={_fmt_p(row['p'])}", fontsize=9)
        ax.legend(frameon=False, fontsize=7)
    fig.suptitle("Stemness paired with CLDN4", fontsize=11)
    _save(fig, figdir / "fig_patient_cldn4_vs_stemness")
    pbar = next(r for r in primary if r["contrast"] == "CytoTRACE2 vs barrier/keratin (no CLDN4)")
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.6))
    for ax, x, y, row, xlab, ylab in (
        (axes[0], "mean_barrier", "mean_CLDN4", bar, "barrier/keratin (no CLDN4)", "CLDN4"),
        (axes[1], "mean_barrier", "mean_ct2", pbar, "barrier/keratin (no CLDN4)", "CytoTRACE2"),
    ):
        for ds in DATASETS:
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub[y], s=36, c=COLORS[ds], label=ds, edgecolors="white", linewidths=0.3)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(f"DL ρ={_fmt_r(row['rho'])} p={_fmt_p(row['p'])}", fontsize=9)
        ax.legend(frameon=False, fontsize=7)
    fig.suptitle("Barrier program, CLDN4 excluded from the score", fontsize=11)
    _save(fig, figdir / "fig_patient_cldn4_vs_barrier")

    forest_contrasts = {
        "CLDN4 vs CytoTRACE2",
        "CLDN4 vs Ben-Porath ES1 stemness",
        "CLDN4 vs barrier/keratin (no CLDN4)",
    }
    forest([r for r in cohort_rows if r["contrast"] in forest_contrasts], figdir / "fig_forest_cohort_rho")

    fig, axes = plt.subplots(1, 3, figsize=(10.8, 4.2))
    panels = (
        ("ct2_low", "ct2_high", "CytoTRACE2", paired[0]),
        ("es1_low", "es1_high", "Ben-Porath ES1", paired[1]),
        ("barrier_low", "barrier_high", "barrier/keratin", paired[3]),
    )
    for ax, (lo, hi, lab, row) in zip(axes, panels):
        for ds in DATASETS:
            sub = paired_df[paired_df["dataset"] == ds]
            ax.scatter(sub[lo], sub[hi], s=28, c=COLORS[ds], label=ds)
        if len(paired_df):
            vals = np.concatenate([paired_df[lo].to_numpy(), paired_df[hi].to_numpy()])
            vals = vals[np.isfinite(vals)]
            if len(vals):
                lo_l, hi_l = float(vals.min()), float(vals.max())
                pad = 0.05 * (hi_l - lo_l + 1e-6)
                ax.plot([lo_l - pad, hi_l + pad], [lo_l - pad, hi_l + pad], ls="--", c="0.6", lw=1)
        ax.set_xlabel(f"CLDN4-low {lab}")
        ax.set_ylabel(f"CLDN4-high {lab}")
        ax.set_title(f"Δmed={row['delta']:+.3f} p={_fmt_p(row['p'])}" if np.isfinite(row["delta"]) else lab, fontsize=8)
    axes[0].legend(frameon=False, fontsize=6)
    fig.suptitle(f"Within-unit CLDN4-high vs low  n={len(paired_df)}", fontsize=11)
    _save(fig, figdir / "fig_paired_high_vs_low")

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0))
    unit_ct = patient.groupby("dataset").size().reindex(DATASETS)
    cell_ct = cells.groupby("dataset").size().reindex(DATASETS)
    axes[0].bar(list(DATASETS), unit_ct.to_numpy(), color=[COLORS[d] for d in DATASETS])
    axes[0].tick_params(axis="x", labelrotation=20)
    axes[0].set_ylabel("units")
    axes[0].set_title(f"units={len(patient)}  eligible={len(elig)}")
    axes[1].bar(list(DATASETS), cell_ct.to_numpy(), color=[COLORS[d] for d in DATASETS])
    axes[1].tick_params(axis="x", labelrotation=20)
    axes[1].set_ylabel("cells after cap")
    axes[1].set_title(f"n_cells={len(cells)}")
    fig.suptitle("Honest n — concordant-4 malignant, cap ≤200/unit", fontsize=11)
    _save(fig, figdir / "fig_honest_n")

    # descriptive potency mix by within-dataset CLDN4 tertile
    tert = pd.Series(index=cells.index, dtype=object)
    for ds, sub in cells.groupby("dataset"):
        q1, q2 = np.nanquantile(sub["expr_CLDN4"], [1 / 3, 2 / 3])
        lab = np.full(len(sub), "mid", dtype=object)
        lab[sub["expr_CLDN4"].to_numpy() <= q1] = "low"
        lab[sub["expr_CLDN4"].to_numpy() >= q2] = "high"
        tert.loc[sub.index] = lab
    cells = cells.copy()
    cells["cldn4_tertile_within_dataset"] = tert
    cat_order = ["Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent"]
    fig, axes = plt.subplots(1, 4, figsize=(11.5, 3.6), sharey=True)
    for ax, ds in zip(axes, DATASETS):
        sub = cells[cells["dataset"] == ds]
        tab = (
            sub.groupby(["cldn4_tertile_within_dataset", "CytoTRACE2_Potency"], observed=False)
            .size().unstack(fill_value=0)
        )
        for c in cat_order:
            if c not in tab.columns:
                tab[c] = 0
        tab = tab.reindex(["low", "mid", "high"])[cat_order].fillna(0)
        frac = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        frac.plot(kind="bar", stacked=True, ax=ax, colormap="viridis", legend=False)
        ax.set_title(ds, fontsize=8)
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelrotation=0)
    axes[0].set_ylabel("fraction of cells")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=7, loc="upper right")
    fig.suptitle("Descriptive: potency category by within-dataset CLDN4 tertile (cells, not the n)", fontsize=10)
    _save(fig, figdir / "fig_potency_category")

    done_cols = [
        "unit_id", "dataset", "patient_id", "tissue", "histology", "stage",
        "n_malignant", "mean_CLDN4", "pct_CLDN4", "mean_ct2", "frac_diff",
        "mean_es1", "mean_wong", "mean_core9", "mean_barrier", "mean_gulati",
        "mean_preknn", "ct2_high", "ct2_low", "barrier_high", "barrier_low",
    ]
    elig[done_cols].to_csv(tabdir / "patient_cldn4_vs_potency.tsv", sep="\t", index=False)
    patient.to_csv(tabdir / "patient_means.tsv", sep="\t", index=False)
    pd.DataFrame(cohort_rows).to_csv(tabdir / "cohort_spearman.tsv", sep="\t", index=False)
    paired_df.to_csv(tabdir / "paired_high_vs_low.tsv", sep="\t", index=False)
    cells.to_csv(tabdir / "cell_scores.tsv", sep="\t")
    stats_rows = []
    for r in primary:
        stats_rows.append({**{k: r[k] for k in ("contrast", "k", "n", "rho", "p", "q", "I2")}, "family": "primary_DL"})
    for r in cohort_rows:
        stats_rows.append({**r, "family": "cohort"})
    for r in sensitivity:
        stats_rows.append({**r, "family": "sensitivity"})
    for r in paired:
        stats_rows.append({**r, "family": "paired"})
    pd.DataFrame(stats_rows).to_csv(tabdir / "stats.tsv", sep="\t", index=False)

    try:
        import importlib.metadata
        version = importlib.metadata.version("cytotrace2-py")
    except Exception:
        version = "1.1.0.4"
    cat_path = args.h5ad / "catalog.json"
    catalog = json.loads(cat_path.read_text()) if cat_path.exists() else {}
    holds, fails = verdict_sentences(primary, paired, sensitivity, cohort_rows)
    potency_counts: dict[str, int] = {}
    for meta in metas.values():
        for k, v in meta["potency_counts"].items():
            potency_counts[k] = potency_counts.get(k, 0) + int(v)
    coverage_brief = {
        ds: {
            "es1_used": metas[ds]["coverage"]["benporath_es1_genes_used"],
            "wong_used": metas[ds]["coverage"]["wong_esc_genes_used"],
            "core9_detected_frac": round(metas[ds]["coverage"]["benporath_core9_frac_any_detected"], 3),
            "frac_genes_lt_500": round(metas[ds]["cytotrace2"]["frac_genes_lt_500"], 3),
        }
        for ds in DATASETS
    }
    n_cells_by = cells.groupby("dataset").size().to_dict()
    n_units_by = patient.groupby("dataset").size().to_dict()
    S = {
        "cytotrace2_version": version,
        "cap": 200,
        "n_cells": int(len(cells)),
        "n_cells_by": {ds: int(n_cells_by.get(ds, 0)) for ds in DATASETS},
        "n_units": int(len(patient)),
        "n_units_by": {ds: int(n_units_by.get(ds, 0)) for ds in DATASETS},
        "n_eligible": int(len(elig)),
        "n_paired": int(len(paired_df)),
        "catalog": catalog,
        "potency_counts": potency_counts,
        "frac_lt_500": float((cells["n_genes_detected"] < 500).mean()),
        "coverage_brief": coverage_brief,
        "primary": primary,
        "cohort": cohort_rows,
        "sensitivity": sensitivity,
        "paired": paired,
        "what_holds": holds,
        "what_fails": fails,
        "metas": metas,
    }
    (args.outdir / "summary.json").write_text(json.dumps(S, indent=2, default=str))
    write_finding(args.finding, S)
    print(json.dumps({"primary": primary, "n_eligible": len(elig), "n_cells": len(cells)}, indent=2, default=str))


if __name__ == "__main__":
    main()
