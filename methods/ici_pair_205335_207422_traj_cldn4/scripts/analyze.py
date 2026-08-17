#!/usr/bin/env python3
"""PAGA + DPT per ICI dataset, then stacked patient tables.

ADDITIVE. CLDN4 only. GSE205335 (RECIST) + GSE207422 (MPR).
Not a 7-pool. No GSE148071. No dual-high. Patient is the unit.
PAGA + DPT are run separately (platforms differ; Harmony is not forced).
Root is leftover AT2-like epithelium, never CLDN4-high.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import (  # noqa: E402
    AUTHOR_MALIGNANT_205335,
    AUTHOR_NONMALIGNANT_205335,
    COMPARATOR,
    CONTROLS,
    FOCAL,
    QC_NEG,
    STATES,
)

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e


LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_GENES = 200
MIN_UMI = 500
MIN_CELLS_PER_GENE = 10
MIN_CELLS_PER_SAMPLE = 10
MIN_CELLS_MALIGNANT = 10


def _bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n, dtype=float)
    running = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        k = n - rank + 1
        running = min(running, p[i] * n / k)
        q[i] = running
    return [float(min(1.0, x)) for x in q]


def _spearman(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 4:
        return {"n": n, "rho": None, "p": None, "note": "n<4"}
    rho, p = stats.spearmanr(x[mask], y[mask])
    if not np.isfinite(rho):
        return {"n": n, "rho": None, "p": None, "note": "undefined"}
    return {"n": n, "rho": float(rho), "p": float(p), "note": ""}


def _mwu(a, b) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    n_a, n_b = int(len(a)), int(len(b))
    if n_a < 2 or n_b < 2:
        return {
            "n_a": n_a,
            "n_b": n_b,
            "median_a": float(np.median(a)) if n_a else None,
            "median_b": float(np.median(b)) if n_b else None,
            "mean_a": float(np.mean(a)) if n_a else None,
            "mean_b": float(np.mean(b)) if n_b else None,
            "U": None,
            "p": None,
            "note": "n<2",
        }
    method = "exact" if (n_a + n_b) <= 20 else "asymptotic"
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided", method=method)
    return {
        "n_a": n_a,
        "n_b": n_b,
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "delta_median_a_minus_b": float(np.median(a) - np.median(b)),
        "U": float(u),
        "p": float(p),
        "method": method,
        "note": "",
    }


def fisher_z_dl(rows: list[dict]) -> dict:
    """DerSimonian–Laird on Fisher-z of Spearman ρ. Descriptive."""
    usable = [r for r in rows if r.get("rho") is not None and r.get("n", 0) >= 4]
    if len(usable) < 1:
        return {"k": 0, "N": 0, "rho": None, "p": None, "I2": None, "note": "no cohorts"}
    if len(usable) == 1:
        r = usable[0]
        return {
            "k": 1,
            "N": r["n"],
            "rho": r["rho"],
            "p": r["p"],
            "I2": None,
            "note": "single cohort",
        }
    zs, ws, ns = [], [], []
    for r in usable:
        n = r["n"]
        rho = max(min(r["rho"], 0.999), -0.999)
        z = np.arctanh(rho)
        var = 1.0 / (n - 3)
        zs.append(z)
        ws.append(1.0 / var)
        ns.append(n)
    zs = np.asarray(zs)
    ws = np.asarray(ws)
    zbar = float(np.sum(ws * zs) / np.sum(ws))
    q = float(np.sum(ws * (zs - zbar) ** 2))
    k = len(usable)
    c = float(np.sum(ws) - np.sum(ws**2) / np.sum(ws))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    w_star = 1.0 / (1.0 / ws + tau2)
    z_re = float(np.sum(w_star * zs) / np.sum(w_star))
    se = float(np.sqrt(1.0 / np.sum(w_star)))
    from scipy.stats import norm

    p = float(2 * norm.sf(abs(z_re / se)))
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
    return {
        "k": k,
        "N": int(sum(ns)),
        "rho": float(np.tanh(z_re)),
        "p": p,
        "I2": float(i2),
        "tau2": tau2,
        "note": "DL Fisher-z",
    }


def _score_mean(adata, genes: tuple[str, ...], key: str) -> list[str]:
    present = [g for g in genes if g in adata.var_names]
    absent = [g for g in genes if g not in adata.var_names]
    if not present:
        adata.obs[key] = np.nan
        return absent
    X = adata[:, present].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    adata.obs[key] = np.asarray(X, dtype=float).mean(axis=1)
    return absent


def _paga_components(connect: np.ndarray, thresh: float = 0.0) -> list[set[int]]:
    n = connect.shape[0]
    seen = [False] * n
    comps = []
    for i in range(n):
        if seen[i]:
            continue
        stack = [i]
        seen[i] = True
        cur = {i}
        while stack:
            u = stack.pop()
            for v in range(n):
                if not seen[v] and connect[u, v] > thresh:
                    seen[v] = True
                    stack.append(v)
                    cur.add(v)
        comps.append(cur)
    return comps


def _pick_root(adata, dataset: str) -> tuple[int, dict]:
    """Leftover AT2-like epithelium. Never CLDN4-high."""
    leftover = adata.obs["is_leftover"].astype(str).eq("True")
    info = {
        "n_leftover": int(leftover.sum()),
        "rule": None,
        "dataset": dataset,
        "never_cldn4_high": True,
    }
    if int(leftover.sum()) >= 20:
        idx = np.flatnonzero(leftover.to_numpy())
        scores = adata.obs.loc[leftover, "score_AT2"].to_numpy()
        med = np.nanmedian(scores)
        pick = idx[int(np.nanargmin(np.abs(scores - med)))]
        info["rule"] = f"{dataset} leftover epithelium (median AT2 score)"
        return int(pick), info
    if int(leftover.sum()) >= 5:
        idx = np.flatnonzero(leftover.to_numpy())
        scores = adata.obs.loc[leftover, "score_AT2"].to_numpy()
        pick = idx[int(np.nanargmax(scores))]
        info["rule"] = f"{dataset} leftover epithelium (max AT2; leftover n<20)"
        return int(pick), info
    if "leiden" not in adata.obs:
        raise SystemExit("leiden missing before root pick")
    means = (
        adata.obs.groupby("leiden", observed=True)["score_AT2"]
        .mean()
        .sort_values(ascending=False)
    )
    top = str(means.index[0])
    cand = adata.obs["leiden"].astype(str) == top
    # among that cluster, pick high AT2 AND low CLDN4
    scores = adata.obs.loc[cand, "score_AT2"].to_numpy() - 0.25 * adata.obs.loc[cand, "expr_CLDN4"].to_numpy()
    idx = np.flatnonzero(cand.to_numpy())
    pick = idx[int(np.nanargmax(scores))]
    info["rule"] = f"{dataset} Leiden {top} max AT2-minus-CLDN4 (leftover n<5)"
    info["fallback_cluster"] = top
    return int(pick), info


def _fmt_rho(row: dict) -> str:
    if row.get("rho") is None:
        return f"n={row.get('n')}, ρ=NA, p=NA"
    return f"n={row['n']}, ρ={row['rho']:.3f}, p={row['p']:.3g}"


def _fmt_mwu(row: dict, a_lab: str, b_lab: str) -> str:
    if row.get("p") is None:
        return f"{a_lab} n={row.get('n_a')} vs {b_lab} n={row.get('n_b')}, p=NA"
    return (
        f"{a_lab} n={row['n_a']} med={row['median_a']:.3f} vs "
        f"{b_lab} n={row['n_b']} med={row['median_b']:.3f}, p={row['p']:.3g}"
    )


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def run_paga_dpt(adata, dataset: str, figdir: Path, tabdir: Path) -> dict:
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()
    adata.var["n_cells"] = np.array((adata.X > 0).sum(axis=0)).ravel()
    adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).ravel()
    adata.obs["n_umi"] = np.array(adata.X.sum(axis=1)).ravel()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    keep = (adata.obs["n_genes"] >= MIN_GENES) & (adata.obs["n_umi"] >= MIN_UMI)
    adata = adata[keep].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    absent: dict[str, list[str]] = {}
    for name, genes in STATES.items():
        absent[name] = _score_mean(adata, genes, f"score_{name}")
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        if g in adata.var_names:
            x = adata[:, g].X
            if hasattr(x, "toarray"):
                x = x.toarray()
            adata.obs[f"expr_{g}"] = np.asarray(x, dtype=float).ravel()
        else:
            adata.obs[f"expr_{g}"] = np.nan
            absent.setdefault("single_genes", []).append(g)

    try:
        sc.pp.highly_variable_genes(adata, layer="counts", flavor="seurat_v3", n_top_genes=N_HVG)
    except ImportError:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=N_HVG)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    # Per-dataset graph. Do not Harmony-force.
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2)
    except TypeError:
        sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.paga(adata, groups="leiden")
    sc.tl.diffmap(adata, n_comps=15)
    sc.tl.umap(adata)

    root_i, root_info = _pick_root(adata, dataset)
    adata.uns["iroot"] = root_i
    sc.tl.dpt(adata, n_dcs=10)
    root_info["index"] = root_i
    root_info["root_patient"] = str(adata.obs.iloc[root_i].get("patient_id", ""))
    root_info["root_subtype"] = str(adata.obs.iloc[root_i].get("author_subtype", ""))
    root_info["root_CLDN4"] = float(adata.obs.iloc[root_i].get("expr_CLDN4", np.nan))
    root_info["root_AT2"] = float(adata.obs.iloc[root_i].get("score_AT2", np.nan))

    dpt = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)
    adata.obs["dpt_pseudotime"] = dpt

    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    comps = _paga_components(connect, thresh=0.0)
    leiden_ids = [str(x) for x in adata.obs["leiden"].cat.categories]
    cluster_tab = []
    for cl in leiden_ids:
        sub = adata.obs[adata.obs["leiden"].astype(str) == cl]
        cluster_tab.append(
            {
                "dataset": dataset,
                "leiden": cl,
                "n_cells": int(len(sub)),
                "n_patients": int(sub["patient_id"].nunique()),
                "n_leftover": int(sub["is_leftover"].astype(str).eq("True").sum()),
                "n_malignant": int(sub["is_malignant"].astype(str).eq("True").sum()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
            }
        )
    pd.DataFrame(cluster_tab).to_csv(tabdir / f"{dataset.lower()}_leiden_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / f"{dataset.lower()}_paga_connectivities.tsv", sep="\t"
    )

    rows = []
    for pid, sub in adata.obs.groupby("patient_id", observed=True):
        n_mal = int(sub["is_malignant"].astype(str).eq("True").sum())
        n_left = int(sub["is_leftover"].astype(str).eq("True").sum())
        mal = sub[sub["is_malignant"].astype(str).eq("True")]
        rows.append(
            {
                "dataset": dataset,
                "patient": str(pid),
                "histology": str(sub["histology"].iloc[0]),
                "timing": str(sub["timing"].iloc[0]),
                "benefit_raw": str(sub["benefit_raw"].iloc[0]),
                "recist": str(sub["recist"].iloc[0]) if "recist" in sub.columns else "NA",
                "platform": str(sub["platform"].iloc[0]) if "platform" in sub.columns else "NA",
                "n_epi_cells": int(len(sub)),
                "n_malignant": n_mal,
                "n_leftover": n_left,
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_CLDN4_malignant": float(mal["expr_CLDN4"].mean()) if n_mal else np.nan,
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_dpt_malignant": float(mal["dpt_pseudotime"].mean()) if n_mal else np.nan,
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_SFTPC": float(sub["expr_SFTPC"].mean()) if "expr_SFTPC" in sub.columns else np.nan,
            }
        )
    patient_df = pd.DataFrame(rows)
    patient_df["eligible_spearman"] = patient_df["n_epi_cells"] >= MIN_CELLS_PER_SAMPLE
    patient_df["eligible_malignant"] = patient_df["n_malignant"] >= MIN_CELLS_MALIGNANT

    if dataset == "GSE205335":
        rec = patient_df["recist"].astype(str)
        patient_df["benefit_class"] = np.where(
            rec.eq("PR"),
            "benefit",
            np.where(rec.eq("PD"), "no_benefit", "excluded"),
        )
        patient_df["benefit_label"] = rec
        patient_df["outcome_scheme"] = "RECIST_PR_vs_PD"
    else:
        raw = patient_df["benefit_raw"].astype(str)
        patient_df["benefit_class"] = np.where(
            raw.eq("MPR"),
            "benefit",
            np.where(raw.eq("NMPR"), "no_benefit", "excluded"),
        )
        patient_df["benefit_label"] = raw
        patient_df["outcome_scheme"] = "MPR_vs_NMPR"

    # within-dataset rank / z so stacked DPT is not raw-scale
    elig_mask = patient_df["eligible_spearman"]
    ranks = patient_df.loc[elig_mask, "mean_dpt"].rank(method="average")
    patient_df["dpt_rank"] = np.nan
    patient_df.loc[elig_mask, "dpt_rank"] = ranks
    mu = patient_df.loc[elig_mask, "mean_dpt"].mean()
    sd = patient_df.loc[elig_mask, "mean_dpt"].std(ddof=1)
    patient_df["dpt_z"] = np.where(elig_mask, (patient_df["mean_dpt"] - mu) / (sd if sd else 1.0), np.nan)
    cldn_ranks = patient_df.loc[elig_mask, "mean_CLDN4"].rank(method="average")
    patient_df["cldn4_rank"] = np.nan
    patient_df.loc[elig_mask, "cldn4_rank"] = cldn_ranks

    patient_df.to_csv(tabdir / f"{dataset.lower()}_patient_means.tsv", sep="\t", index=False)

    elig = patient_df[patient_df["eligible_spearman"]].copy()
    malig = patient_df[patient_df["eligible_malignant"]].copy()
    contrasts = [
        ("CLDN4 vs DPT", "mean_CLDN4", "mean_dpt"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs club score", "mean_CLDN4", "mean_club"),
        ("CLDN4 vs basal score", "mean_CLDN4", "mean_basal"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs malignant-like score", "mean_CLDN4", "mean_malignant_like"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("SFTPC vs DPT (control)", "mean_SFTPC", "mean_dpt"),
        ("AT2 score vs DPT (control)", "mean_AT2", "mean_dpt"),
    ]
    primary = []
    for name, a, b in contrasts:
        primary.append({"dataset": dataset, "contrast": name, **_spearman(elig[a], elig[b])})
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q

    sensitivity = [
        {
            "dataset": dataset,
            "contrast": "malignant-only CLDN4 vs DPT",
            **_spearman(malig["mean_CLDN4_malignant"], malig["mean_dpt_malignant"]),
        },
        {
            "dataset": dataset,
            "contrast": "malignant-only CLDN4 vs AT2",
            **_spearman(malig["mean_CLDN4_malignant"], malig["mean_AT2"]),
        },
    ]
    if dataset == "GSE205335":
        nsclc = elig[elig["histology"].isin(["ADC", "SQ"])]
        sensitivity.append(
            {
                "dataset": dataset,
                "contrast": "ADC+SQ CLDN4 vs DPT",
                **_spearman(nsclc["mean_CLDN4"], nsclc["mean_dpt"]),
            }
        )

    # outcome tests
    if dataset == "GSE205335":
        ben = elig[elig["benefit_label"] == "PR"]
        nben = elig[elig["benefit_label"] == "PD"]
        a_lab, b_lab = "PR", "PD"
        extra_sd = elig[elig["benefit_label"].isin(["PR", "PD", "SD"])]
        extra_rnr = (
            extra_sd.assign(
                rnr=np.where(extra_sd["benefit_label"].eq("PR"), "R", "NR")
            )
        )
    else:
        ben = elig[elig["benefit_label"] == "MPR"]
        nben = elig[elig["benefit_label"] == "NMPR"]
        a_lab, b_lab = "MPR", "NMPR"
        extra_rnr = None

    outcome = [
        {
            "dataset": dataset,
            "contrast": f"DPT {a_lab} vs {b_lab}",
            "scheme": patient_df["outcome_scheme"].iloc[0],
            **_mwu(ben["mean_dpt"], nben["mean_dpt"]),
            "label_a": a_lab,
            "label_b": b_lab,
        },
        {
            "dataset": dataset,
            "contrast": f"CLDN4 {a_lab} vs {b_lab}",
            "scheme": patient_df["outcome_scheme"].iloc[0],
            **_mwu(ben["mean_CLDN4"], nben["mean_CLDN4"]),
            "label_a": a_lab,
            "label_b": b_lab,
        },
    ]
    if extra_rnr is not None:
        r = extra_rnr[extra_rnr["rnr"] == "R"]
        nr = extra_rnr[extra_rnr["rnr"] == "NR"]
        outcome.append(
            {
                "dataset": dataset,
                "contrast": "DPT RECIST R(PR) vs NR(SD+PD) [sensitivity]",
                "scheme": "RECIST_R_vs_NR",
                **_mwu(r["mean_dpt"], nr["mean_dpt"]),
                "label_a": "R",
                "label_b": "NR",
            }
        )
        nsclc_ben = ben[ben["histology"].isin(["ADC", "SQ"])]
        nsclc_nben = nben[nben["histology"].isin(["ADC", "SQ"])]
        outcome.append(
            {
                "dataset": dataset,
                "contrast": "DPT PR vs PD ADC+SQ [sensitivity]",
                "scheme": "RECIST_PR_vs_PD_NSCLC",
                **_mwu(nsclc_ben["mean_dpt"], nsclc_nben["mean_dpt"]),
                "label_a": "PR",
                "label_b": "PD",
            }
        )

    # figures
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.0))
    sc.pl.paga(
        adata,
        color="expr_CLDN4",
        ax=axes[0, 0],
        show=False,
        frameon=False,
        cmap="viridis",
        title=f"{dataset} PAGA mean CLDN4",
    )
    sc.pl.umap(adata, color="expr_CLDN4", ax=axes[0, 1], show=False, frameon=False, cmap="viridis", title="UMAP CLDN4")
    sc.pl.umap(adata, color="dpt_pseudotime", ax=axes[1, 0], show=False, frameon=False, cmap="magma", title="UMAP DPT (leftover-AT2 root)")
    ax = axes[1, 1]
    c4 = next(r for r in primary if r["contrast"] == "CLDN4 vs DPT")
    colors = {"benefit": "#2a6f97", "no_benefit": "#b23a48", "excluded": "#8a8a8a"}
    for cls, col in colors.items():
        sub = elig[elig["benefit_class"] == cls]
        if sub.empty:
            continue
        ax.scatter(sub["mean_dpt"], sub["mean_CLDN4"], s=52, c=col, label=f"{cls} n={len(sub)}")
    ax.set_xlabel("patient-mean DPT")
    ax.set_ylabel("patient-mean CLDN4")
    ax.set_title(f"{dataset} CLDN4 vs DPT  {_fmt_rho(c4)}")
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(f"{dataset} epithelium  n_cells={adata.n_obs}  n_patients={len(patient_df)}", fontsize=11)
    _save(fig, figdir / f"fig_{dataset.lower()}_trajectory_cldn4")

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    axes[0].bar(
        patient_df["patient"].astype(str),
        patient_df["n_epi_cells"],
        color="#4c6a92",
    )
    axes[0].tick_params(axis="x", rotation=60, labelsize=7)
    axes[0].set_ylabel("epithelial cells in object")
    axes[0].set_title(f"{dataset} cells/patient (not n)")
    lab_order = sorted(elig["benefit_label"].unique())
    for lab in lab_order:
        sub = elig[elig["benefit_label"] == lab]
        axes[1].scatter(np.repeat(lab, len(sub)), sub["mean_dpt"], s=48)
    axes[1].set_ylabel("patient-mean DPT")
    axes[1].set_title("DPT by outcome label")
    _save(fig, figdir / f"fig_{dataset.lower()}_honest_n")

    for color, fname, cmap in (
        ("author_subtype", f"fig_{dataset.lower()}_umap_subtype", None),
        ("dpt_pseudotime", f"fig_{dataset.lower()}_umap_dpt", "magma"),
        ("score_AT2", f"fig_{dataset.lower()}_umap_AT2", "viridis"),
        ("score_barrier_keratin", f"fig_{dataset.lower()}_umap_barrier", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    return {
        "dataset": dataset,
        "n_cells": int(adata.n_obs),
        "n_patients": int(len(patient_df)),
        "n_eligible": int(len(elig)),
        "n_malignant_eligible": int(len(malig)),
        "lineage_counts": {
            "malignant": int(adata.obs["is_malignant"].astype(str).eq("True").sum()),
            "leftover": int(adata.obs["is_leftover"].astype(str).eq("True").sum()),
        },
        "subtype_counts": adata.obs["author_subtype"].astype(str).value_counts().to_dict(),
        "histology_patients": patient_df.groupby("histology")["patient"].nunique().to_dict(),
        "benefit_counts": elig["benefit_label"].value_counts().to_dict(),
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "harmony": {"used": False, "reason": "per-dataset graph; Harmony not forced"},
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "outcome": outcome,
        "patient_df": patient_df,
        "qc": {
            "min_genes": MIN_GENES,
            "min_umi": MIN_UMI,
            "EPCAM_mean": float(adata.obs["expr_EPCAM"].mean()) if "expr_EPCAM" in adata.obs else None,
            "PTPRC_mean": float(adata.obs["expr_PTPRC"].mean()) if "expr_PTPRC" in adata.obs else None,
        },
    }


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    a = s["GSE205335"]
    b = s["GSE207422"]
    stacked = s["stacked"]

    def rho_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r.get('dataset', '')} | {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    def mwu_md(r: dict) -> str:
        pa = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        ma = "NA" if r.get("median_a") is None else f"{r['median_a']:.3f}"
        mb = "NA" if r.get("median_b") is None else f"{r['median_b']:.3f}"
        return (
            f"| {r.get('dataset', 'stacked')} | {r['contrast']} | "
            f"{r.get('n_a')} | {r.get('n_b')} | {ma} | {mb} | {pa} |"
        )

    lines = [
        "# Finding — ICI pair GSE205335 + GSE207422, CLDN4-only PAGA/DPT vs outcome",
        "",
        "ADDITIVE. **CLDN4 only.** The two public ICI scRNA sets that carry outcome labels: "
        "**GSE205335** (Ahn/Lee eLife 2024; palliative ICI; RECIST) + **GSE207422** "
        "(Hu *Genome Med* 2023; neoadjuvant PD-1 + chemo; MPR). "
        "**Not a 7-pool. No GSE148071. No dual-high / TACSTD2 gate.** "
        "Patient is the unit. Honest n is not cell count.",
        "",
        "PAGA + leftover-AT2-rooted diffusion pseudotime are run **per dataset**. "
        "Platforms differ (GSE205335 mixed 3′/5′ 10x biopsies/effusions; GSE207422 10x post-surgery). "
        "Harmony was **not** forced. Raw DPT is not commensurate across datasets; "
        "the stacked DPT tests use within-dataset rank / z.",
        "",
        "Root is leftover (non-malignant / A3-normal-lung) epithelium at median AT2 score. "
        "**Root is not CLDN4-high.** Barrier/keratin excludes CLDN4.",
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        "### GSE205335 (RECIST)",
        "",
        f"- GEO patients / samples: **26 / 33**. Catalog, not the test n.",
        f"- Tumor-tissue epithelium in the object: **n_cells = {a['n_cells']}** in **{a['n_patients']}** patients "
        f"(normal LN/lung/brain dropped; 4 normal-only patients never enter).",
        f"- Patients with ≥{MIN_CELLS_PER_SAMPLE} epithelial cells (Spearman / DPT): **n = {a['n_eligible']}**.",
        f"- Patients with ≥{MIN_CELLS_MALIGNANT} author malignant cells: **n = {a['n_malignant_eligible']}**.",
        f"- RECIST on eligible patients: {a['benefit_counts']}.",
        f"- Primary outcome contrast is **PR vs PD** (SD and NE excluded). "
        f"R vs NR (PR vs SD+PD) is sensitivity only.",
        f"- Lineage on the object: {a['lineage_counts']}. Author subtypes: {a['subtype_counts']}.",
        f"- Histology (patients): {a['histology_patients']}. SCLC/NUT stay in the all-comer row.",
        f"- DPT root: {a['root']['rule']} (patient {a['root'].get('root_patient')}; "
        f"CLDN4={a['root'].get('root_CLDN4'):.3f}; AT2={a['root'].get('root_AT2'):.3f}).",
        f"- PAGA components at connectivity>0: **{a['n_paga_components']}** among {a['n_leiden']} Leiden vertices.",
        "- MPR is not a GEO field here. RECIST is not substituted for MPR.",
        "",
        "### GSE207422 (MPR)",
        "",
        f"- GEO samples: **15** (3 pre-treatment TN + 12 post-surgery). Catalog, not the test n.",
        f"- Post-treatment epithelium in the object: **n_cells = {b['n_cells']}** in **{b['n_patients']}** patients. "
        "The 3 pre-biopsies are out of the graph.",
        f"- Patients with ≥{MIN_CELLS_PER_SAMPLE} epithelial cells: **n = {b['n_eligible']}**.",
        f"- Patients with ≥{MIN_CELLS_MALIGNANT} A3-malignant-like cells: **n = {b['n_malignant_eligible']}**. "
        "A3 empties some MPR residuals (normal-lung program); that is reported, not patched.",
        f"- Pathologic response on eligible patients: {b['benefit_counts']} (pCR P06 = MPR).",
        f"- Lineage on the object: {b['lineage_counts']}.",
        f"- DPT root: {b['root']['rule']} (patient {b['root'].get('root_patient')}; "
        f"CLDN4={b['root'].get('root_CLDN4'):.3f}; AT2={b['root'].get('root_AT2'):.3f}).",
        f"- PAGA components at connectivity>0: **{b['n_paga_components']}** among {b['n_leiden']} Leiden vertices.",
        "- Author CopyKAT / DRMref barcodes are not on GEO. A3 marker malignant-like is used.",
        "",
        "### Stacked (not a joint embedding)",
        "",
        f"- Stacked patient rows: **N = {stacked['n_rows']}** "
        f"(GSE205335 {stacked['n_205335']} + GSE207422 {stacked['n_207422']}).",
        f"- Stacked CLDN4 vs within-dataset DPT rank: **n = {stacked['cldn4_vs_dpt_rank']['n']}**.",
        f"- Stacked benefit vs no-benefit (PR+MPR vs PD+NMPR) on DPT z: "
        f"**n = {stacked['dpt_vs_benefit']['n_a']} vs {stacked['dpt_vs_benefit']['n_b']}**. "
        "SD/NE/TN are out of the binary stack.",
        "- Raw DPT is **not** pooled. Harmony was not run.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        "- Batch: none (per-dataset PCA/neighbors). Harmony not forced.",
        "- DPT root: leftover epithelium, median AT2 score. Never CLDN4-high.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- GSE205335 malignant = author `Malignant cells`. GSE207422 malignant-like = epithelial AND zero UMI of SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3.",
        "",
        "## Per-dataset — CLDN4 vs DPT (patient Spearman, BH inside each dataset list)",
        "",
        "| Dataset | Contrast | n | ρ | p | q |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for r in a["primary_spearman"] + b["primary_spearman"]:
        lines.append(rho_md(r))
    lines += [
        "",
        "## Per-dataset — DPT vs ICI benefit (patient MWU)",
        "",
        "| Dataset | Contrast | n_a | n_b | med_a | med_b | p |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in a["outcome"] + b["outcome"]:
        lines.append(mwu_md(r))
    lines += [
        "",
        "## Stacked patient table — CLDN4 vs DPT",
        "",
        "DPT is rank- and z-transformed **within dataset** before stacking. "
        "Fisher-z DL is descriptive (k=2).",
        "",
        "| Contrast | k | N | ρ | p | I² | note |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for r in stacked["cldn4_vs_dpt_rows"]:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        i2 = "—" if r.get("I2") is None else f"{100 * r['I2']:.0f}%"
        lines.append(
            f"| {r['contrast']} | {r.get('k', 1)} | {r['n'] if 'n' in r else r.get('N')} | "
            f"{rho} | {pv} | {i2} | {r.get('note', '')} |"
        )
    lines += [
        "",
        "## Stacked — DPT vs ICI benefit (PR+MPR vs PD+NMPR)",
        "",
        "| Contrast | n_benefit | n_no_benefit | med_z_benefit | med_z_no | p |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in stacked["benefit_rows"]:
        pa = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        ma = "NA" if r.get("median_a") is None else f"{r['median_a']:.3f}"
        mb = "NA" if r.get("median_b") is None else f"{r['median_b']:.3f}"
        lines.append(
            f"| {r['contrast']} | {r.get('n_a')} | {r.get('n_b')} | {ma} | {mb} | {pa} |"
        )
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Dataset | Contrast | n | ρ | p |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for r in a["sensitivity_spearman"] + b["sensitivity_spearman"]:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['dataset']} | {r['contrast']} | {r['n']} | {rho} | {pv} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- This is not a 7-pool, not GSE148071, not a TACSTD2 redo, not dual-high.",
        "- Per-dataset DPT clocks are not a shared latent time. Do not interpret stacked raw DPT.",
        "- Do not write “AT2 differentiates into ICI-resistant NSCLC because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- GSE205335 RECIST PR vs PD is not MPR. GSE207422 MPR is not RECIST.",
        "- A3 malignant-like is not CopyKAT. Author malignant on GSE205335 is not CNV re-called here.",
        "- Harmony was not run. A null stacked test is not evidence that a joint embedding would be null.",
        "- DPT is an ordering, not a clock. Slingshot was not required.",
        "",
        "## Outputs (done criterion)",
        "",
        "- `results/tables/gse205335_patient_means.tsv` — per-dataset patient table",
        "- `results/tables/gse207422_patient_means.tsv` — per-dataset patient table",
        "- `results/tables/per_dataset_cldn4_vs_dpt.tsv`",
        "- `results/tables/per_dataset_dpt_vs_benefit.tsv`",
        "- `results/tables/stacked_patient_table.tsv`",
        "- `results/tables/stacked_cldn4_vs_dpt.tsv`",
        "- `results/tables/stacked_dpt_vs_benefit.tsv`",
        "- `results/tables/honest_n.tsv`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/ici_pair_205335_207422_traj_cldn4/requirements.txt",
        "python3 methods/ici_pair_205335_207422_traj_cldn4/scripts/download.py \\",
        "  --out /tmp/ici_pair_205335_207422",
        "python3 methods/ici_pair_205335_207422_traj_cldn4/scripts/extract.py \\",
        "  --data /tmp/ici_pair_205335_207422 \\",
        "  --out /tmp/ici_pair_205335_207422/extracted",
        "python3 methods/ici_pair_205335_207422_traj_cldn4/scripts/analyze.py \\",
        "  --extracted /tmp/ici_pair_205335_207422/extracted \\",
        "  --outdir methods/ici_pair_205335_207422_traj_cldn4/results \\",
        "  --finding methods/ici_pair_205335_207422_traj_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--extracted", type=Path, default=Path("/tmp/ici_pair_205335_207422/extracted"))
    p.add_argument("--outdir", type=Path, default=Path("methods/ici_pair_205335_207422_traj_cldn4/results"))
    p.add_argument("--finding", type=Path, default=Path("methods/ici_pair_205335_207422_traj_cldn4/FINDING.md"))
    args = p.parse_args()
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    runs = {}
    for ds, fname in (
        ("GSE205335", "GSE205335_epithelium.h5ad"),
        ("GSE207422", "GSE207422_epithelium.h5ad"),
    ):
        path = Path(args.extracted) / fname
        if not path.is_file():
            raise SystemExit(f"missing {path}")
        print(f"=== {ds} PAGA+DPT ===", flush=True)
        adata = sc.read_h5ad(path)
        runs[ds] = run_paga_dpt(adata, ds, figdir, tabdir)
        del adata

    # stacked patient table
    stacked_df = pd.concat(
        [runs["GSE205335"]["patient_df"], runs["GSE207422"]["patient_df"]],
        ignore_index=True,
    )
    stacked_df.to_csv(tabdir / "stacked_patient_table.tsv", sep="\t", index=False)

    per_dpt = pd.DataFrame(
        runs["GSE205335"]["primary_spearman"] + runs["GSE207422"]["primary_spearman"]
    )
    per_dpt.to_csv(tabdir / "per_dataset_cldn4_vs_dpt.tsv", sep="\t", index=False)
    per_out = pd.DataFrame(runs["GSE205335"]["outcome"] + runs["GSE207422"]["outcome"])
    per_out.to_csv(tabdir / "per_dataset_dpt_vs_benefit.tsv", sep="\t", index=False)

    elig = stacked_df[stacked_df["eligible_spearman"]].copy()
    c4_205 = next(r for r in runs["GSE205335"]["primary_spearman"] if r["contrast"] == "CLDN4 vs DPT")
    c4_207 = next(r for r in runs["GSE207422"]["primary_spearman"] if r["contrast"] == "CLDN4 vs DPT")
    stacked_rank = {"contrast": "stacked CLDN4 vs DPT rank (within dataset)", **_spearman(elig["mean_CLDN4"], elig["dpt_rank"])}
    stacked_z = {"contrast": "stacked CLDN4 vs DPT z (within dataset)", **_spearman(elig["mean_CLDN4"], elig["dpt_z"])}
    dl = fisher_z_dl(
        [
            {**c4_205, "contrast": "GSE205335 CLDN4 vs DPT"},
            {**c4_207, "contrast": "GSE207422 CLDN4 vs DPT"},
        ]
    )
    cldn4_vs_dpt_rows = [
        {**c4_205, "contrast": "GSE205335 CLDN4 vs DPT", "k": 1, "I2": None, "note": "per-dataset raw DPT"},
        {**c4_207, "contrast": "GSE207422 CLDN4 vs DPT", "k": 1, "I2": None, "note": "per-dataset raw DPT"},
        {**stacked_rank, "k": 2, "I2": None, "note": "pooled patients; DPT ranked within dataset"},
        {**stacked_z, "k": 2, "I2": None, "note": "pooled patients; DPT z within dataset"},
        {
            "contrast": "DL Fisher-z of per-dataset CLDN4 vs DPT",
            "k": dl["k"],
            "N": dl["N"],
            "n": dl["N"],
            "rho": dl["rho"],
            "p": dl["p"],
            "I2": dl["I2"],
            "note": dl["note"],
        },
    ]
    pd.DataFrame(cldn4_vs_dpt_rows).to_csv(tabdir / "stacked_cldn4_vs_dpt.tsv", sep="\t", index=False)

    ben = elig[elig["benefit_class"] == "benefit"]
    nben = elig[elig["benefit_class"] == "no_benefit"]
    dpt_ben = {
        "contrast": "stacked DPT z, benefit (PR+MPR) vs no-benefit (PD+NMPR)",
        **_mwu(ben["dpt_z"], nben["dpt_z"]),
    }
    cldn_ben = {
        "contrast": "stacked CLDN4, benefit (PR+MPR) vs no-benefit (PD+NMPR)",
        **_mwu(ben["mean_CLDN4"], nben["mean_CLDN4"]),
    }
    dpt_rank_ben = {
        "contrast": "stacked DPT rank, benefit vs no-benefit",
        **_mwu(ben["dpt_rank"], nben["dpt_rank"]),
    }
    benefit_rows = [dpt_ben, dpt_rank_ben, cldn_ben]
    # per-dataset binary already in outcome; add NSCLC-only stack
    nsclc = elig[~((elig["dataset"] == "GSE205335") & ~elig["histology"].isin(["ADC", "SQ"]))]
    nsclc_ben = nsclc[nsclc["benefit_class"] == "benefit"]
    nsclc_nben = nsclc[nsclc["benefit_class"] == "no_benefit"]
    benefit_rows.append(
        {
            "contrast": "stacked DPT z, benefit vs no-benefit, drop GSE205335 SCLC/NUT",
            **_mwu(nsclc_ben["dpt_z"], nsclc_nben["dpt_z"]),
        }
    )
    pd.DataFrame(benefit_rows).to_csv(tabdir / "stacked_dpt_vs_benefit.tsv", sep="\t", index=False)

    honest = pd.DataFrame(
        [
            {
                "item": "GSE205335 GEO patients",
                "n": 26,
                "note": "catalog; 33 GSM",
            },
            {
                "item": "GSE205335 patients in graph",
                "n": runs["GSE205335"]["n_patients"],
                "note": "tumor epithelium; normals dropped",
            },
            {
                "item": "GSE205335 Spearman n",
                "n": runs["GSE205335"]["n_eligible"],
                "note": f">={MIN_CELLS_PER_SAMPLE} epi cells",
            },
            {
                "item": "GSE205335 PR vs PD",
                "n": int((elig.dataset.eq("GSE205335") & elig.benefit_label.isin(["PR", "PD"])).sum()),
                "note": "SD/NE excluded",
            },
            {
                "item": "GSE207422 GEO samples",
                "n": 15,
                "note": "3 pre + 12 post; catalog",
            },
            {
                "item": "GSE207422 patients in graph",
                "n": runs["GSE207422"]["n_patients"],
                "note": "post-treatment only",
            },
            {
                "item": "GSE207422 Spearman n",
                "n": runs["GSE207422"]["n_eligible"],
                "note": f">={MIN_CELLS_PER_SAMPLE} epi cells",
            },
            {
                "item": "GSE207422 MPR vs NMPR",
                "n": int((elig.dataset.eq("GSE207422") & elig.benefit_label.isin(["MPR", "NMPR"])).sum()),
                "note": "pCR counted as MPR",
            },
            {
                "item": "stacked patient rows",
                "n": int(len(elig)),
                "note": "eligible Spearman patients",
            },
            {
                "item": "stacked benefit vs no-benefit",
                "n": int(len(ben) + len(nben)),
                "note": f"{len(ben)} vs {len(nben)}; PR+MPR vs PD+NMPR",
            },
        ]
    )
    honest.to_csv(tabdir / "honest_n.tsv", sep="\t", index=False)

    # stacked figure
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2))
    ds_col = {"GSE205335": "#2a6f97", "GSE207422": "#b23a48"}
    for ds, col in ds_col.items():
        sub = elig[elig["dataset"] == ds]
        axes[0].scatter(sub["dpt_rank"], sub["mean_CLDN4"], s=52, c=col, label=f"{ds} n={len(sub)}")
    axes[0].set_xlabel("within-dataset DPT rank")
    axes[0].set_ylabel("patient-mean CLDN4")
    axes[0].set_title(_fmt_rho(stacked_rank))
    axes[0].legend(fontsize=8, frameon=False)
    rng = np.random.default_rng(0)
    for cls, col, lab in (
        ("benefit", "#2a6f97", "PR+MPR"),
        ("no_benefit", "#b23a48", "PD+NMPR"),
    ):
        sub = elig[elig["benefit_class"] == cls]
        x = np.where(sub["benefit_class"].eq("benefit"), 0, 1).astype(float)
        x = x + rng.uniform(-0.08, 0.08, size=len(x))
        axes[1].scatter(x, sub["dpt_z"], s=52, c=col, label=f"{lab} n={len(sub)}")
    axes[1].set_xticks([0, 1], ["benefit (PR+MPR)", "no-benefit (PD+NMPR)"])
    axes[1].set_ylabel("within-dataset DPT z")
    axes[1].set_title(_fmt_mwu(dpt_ben, "benefit", "no-benefit"))
    axes[1].legend(fontsize=8, frameon=False)
    fig.suptitle(f"Stacked ICI pair  N={len(elig)} patients  (DPT not raw-pooled)", fontsize=11)
    _save(fig, figdir / "fig_stacked_cldn4_dpt_benefit")

    # verdict text
    o205 = next(r for r in runs["GSE205335"]["outcome"] if r["contrast"].startswith("DPT PR"))
    o207 = next(r for r in runs["GSE207422"]["outcome"] if r["contrast"].startswith("DPT MPR"))
    parts = [
        f"GSE205335 patient-level CLDN4 vs leftover-AT2-rooted DPT: {_fmt_rho(c4_205)}.",
        f"GSE205335 DPT RECIST PR vs PD: {_fmt_mwu(o205, 'PR', 'PD')}.",
        f"GSE207422 patient-level CLDN4 vs leftover-AT2-rooted DPT: {_fmt_rho(c4_207)}.",
        f"GSE207422 DPT MPR vs NMPR: {_fmt_mwu(o207, 'MPR', 'NMPR')}.",
        f"Stacked CLDN4 vs within-dataset DPT rank: {_fmt_rho(stacked_rank)}.",
        f"Stacked DPT z, PR+MPR vs PD+NMPR: {_fmt_mwu(dpt_ben, 'benefit', 'no-benefit')}.",
        f"DL Fisher-z of the two per-dataset CLDN4–DPT ρ: k={dl['k']}, N={dl['N']}, "
        f"ρ={'NA' if dl['rho'] is None else f'{dl['rho']:.3f}'}, "
        f"p={'NA' if dl['p'] is None else f'{dl['p']:.3g}'}.",
        "Harmony was not run. Not a 7-pool. No GSE148071. No dual-high. Root is not CLDN4-high.",
    ]
    verdict = " ".join(parts)

    # strip dataframes from runs for json
    summary = {
        "accessions": ["GSE205335", "GSE207422"],
        "not_a_7pool": True,
        "gse148071": False,
        "dual_high": False,
        "primary_gene": "CLDN4",
        "unit": "patient",
        "harmony_forced": False,
        "GSE205335": {k: v for k, v in runs["GSE205335"].items() if k != "patient_df"},
        "GSE207422": {k: v for k, v in runs["GSE207422"].items() if k != "patient_df"},
        "stacked": {
            "n_rows": int(len(elig)),
            "n_205335": int((elig["dataset"] == "GSE205335").sum()),
            "n_207422": int((elig["dataset"] == "GSE207422").sum()),
            "cldn4_vs_dpt_rank": stacked_rank,
            "cldn4_vs_dpt_z": stacked_z,
            "dl": dl,
            "dpt_vs_benefit": dpt_ben,
            "cldn4_vs_dpt_rows": cldn4_vs_dpt_rows,
            "benefit_rows": benefit_rows,
        },
        "verdict": verdict,
        "done_tables": [
            "gse205335_patient_means.tsv",
            "gse207422_patient_means.tsv",
            "per_dataset_cldn4_vs_dpt.tsv",
            "per_dataset_dpt_vs_benefit.tsv",
            "stacked_patient_table.tsv",
            "stacked_cldn4_vs_dpt.tsv",
            "stacked_dpt_vs_benefit.tsv",
        ],
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    print(
        json.dumps(
            {
                "ok": True,
                "n_205335": runs["GSE205335"]["n_eligible"],
                "n_207422": runs["GSE207422"]["n_eligible"],
                "stacked_n": int(len(elig)),
                "finding": str(args.finding),
                "tables": str(tabdir),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
