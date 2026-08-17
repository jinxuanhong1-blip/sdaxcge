#!/usr/bin/env python3
"""PAGA + DPT on the subsampled winning-pair malignant (and T/NK) slice.

CLDN4 only. No dual-high. Unit = GSE131907 sample / GSE205335 patient.
DPT root = lowest-CLDN4 cell among the top AT2-score tercile of malignant
cells — never a CLDN4-high cell, and not the PR #325 nLung AT2 root
(those cells are not in this malignant metastasis + ICI tumor merge).
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lib.gene_sets import COMPARATOR, FOCAL, STATES  # noqa: E402

LEIDEN_RES = 0.6
N_HVG = 2000
N_NEIGHBORS = 20
N_PCS = 20
MIN_CELLS = 8
MIN_TERTILE = 6


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
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4:
        return {"n": n, "rho": float("nan"), "p": float("nan")}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


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


def quartile_map(given: Path) -> dict[str, str]:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from run_sccoda import load_units

    comp = load_units(given)
    return dict(zip(comp["unit_id"], comp["quartile"]))


def score_and_embed(adata, batch_key: str = "dataset"):
    import scanpy as sc

    sc.pp.filter_cells(adata, min_genes=100)
    if adata.n_obs < 200:
        raise SystemExit(f"too few cells after QC: {adata.n_obs}")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for name, genes in STATES.items():
        _score_mean(adata, genes, f"score_{name}")
    _score_mean(adata, FOCAL, "score_CLDN4")
    _score_mean(adata, COMPARATOR, "score_TACSTD2")
    flavor = "seurat"
    try:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=N_HVG, batch_key=batch_key, flavor=flavor
        )
    except Exception:
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor=flavor)
    hvg = adata[:, adata.var["highly_variable"]].copy()
    # scale on the HVG slice only (already subsampled)
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=N_PCS, svd_solver="arpack")
    used_rep = "X_pca"
    harmony_ok = False
    try:
        sc.external.pp.harmony_integrate(hvg, key=batch_key)
        used_rep = "X_pca_harmony"
        harmony_ok = True
    except Exception as exc:  # noqa: BLE001
        print(f"Harmony failed, using PCA: {exc}", flush=True)
    sc.pp.neighbors(hvg, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, use_rep=used_rep)
    sc.tl.leiden(hvg, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2)
    sc.tl.paga(hvg)
    sc.tl.umap(hvg, init_pos="paga")
    return hvg, {"harmony": harmony_ok, "use_rep": used_rep, "n_cells": int(hvg.n_obs)}


def pick_malignant_root(adata) -> tuple[int, dict]:
    cldn = adata.obs["score_CLDN4"].to_numpy(dtype=float)
    at2 = adata.obs["score_AT2"].to_numpy(dtype=float)
    info = {"rule": None, "n": int(len(cldn))}
    ok = np.isfinite(cldn) & np.isfinite(at2)
    if ok.sum() < 20:
        pick = int(np.nanargmin(cldn))
        info["rule"] = "min CLDN4 (AT2 score thin)"
        return pick, info
    thr = np.nanpercentile(at2[ok], 66)
    cand = np.flatnonzero(ok & (at2 >= thr))
    if cand.size == 0:
        cand = np.flatnonzero(ok)
    pick = int(cand[int(np.nanargmin(cldn[cand]))])
    info["rule"] = "min CLDN4 among top AT2 tercile (never CLDN4-high)"
    info["root_CLDN4"] = float(cldn[pick])
    info["root_AT2"] = float(at2[pick])
    return pick, info


def unit_means(adata, min_cells: int = MIN_CELLS) -> pd.DataFrame:
    cols = [
        c
        for c in adata.obs.columns
        if c.startswith("score_") or c in {"dpt_pseudotime", "dataset", "unit_id", "unit_cldn4_pct", "unit_quartile"}
    ]
    rows = []
    for uid, sub in adata.obs.groupby("unit_id", observed=True):
        if len(sub) < min_cells:
            continue
        rec = {
            "unit_id": uid,
            "dataset": sub["dataset"].iloc[0],
            "n_cells": int(len(sub)),
            "unit_cldn4_pct": float(sub["unit_cldn4_pct"].iloc[0]) if "unit_cldn4_pct" in sub else np.nan,
            "unit_quartile": str(sub["unit_quartile"].iloc[0]) if "unit_quartile" in sub else "NA",
        }
        for c in [
            "score_CLDN4",
            "score_TACSTD2",
            "score_AT2",
            "score_barrier_keratin",
            "score_club",
            "score_basal",
            "score_malignant_like",
            "score_tnk_cytotoxic",
            "score_tnk_exh",
            "score_tnk_naive",
            "dpt_pseudotime",
        ]:
            if c in sub:
                rec[c] = float(np.nanmean(sub[c].to_numpy(dtype=float)))
        rows.append(rec)
    return pd.DataFrame(rows)


def paga_tests(means: pd.DataFrame) -> pd.DataFrame:
    contrasts = [
        ("CLDN4 vs DPT", "score_CLDN4", "dpt_pseudotime"),
        ("CLDN4 vs AT2", "score_CLDN4", "score_AT2"),
        ("CLDN4 vs barrier/keratin", "score_CLDN4", "score_barrier_keratin"),
        ("CLDN4 vs club", "score_CLDN4", "score_club"),
        ("CLDN4 vs basal", "score_CLDN4", "score_basal"),
        ("CLDN4 vs malignant-like", "score_CLDN4", "score_malignant_like"),
        ("CLDN4 vs TACSTD2 (comparator)", "score_CLDN4", "score_TACSTD2"),
    ]
    rows = []
    slices = {"merge": means}
    for ds in sorted(means["dataset"].unique()):
        slices[ds] = means[means["dataset"] == ds]
    for sl, frame in slices.items():
        for name, x, y in contrasts:
            if x not in frame or y not in frame:
                continue
            sp = _spearman(frame[x], frame[y])
            rows.append({"slice": sl, "contrast": name, **sp})
    tab = pd.DataFrame(rows)
    merge = tab["slice"] == "merge"
    tab["q_bh_merge"] = np.nan
    if merge.any():
        tab.loc[merge, "q_bh_merge"] = _bh(tab.loc[merge, "p"].fillna(1).tolist())
    return tab


def paired_tertile(adata) -> pd.DataFrame:
    rows = []
    for uid, sub in adata.obs.groupby("unit_id", observed=True):
        cldn = sub["score_CLDN4"].to_numpy(dtype=float)
        if np.isfinite(cldn).sum() < 3 * MIN_TERTILE:
            continue
        try:
            tert = pd.qcut(pd.Series(cldn, index=sub.index).rank(method="average"), 3, labels=["L", "M", "H"], duplicates="drop")
        except ValueError:
            continue
        hi = sub.loc[tert == "H"]
        lo = sub.loc[tert == "L"]
        if len(hi) < MIN_TERTILE or len(lo) < MIN_TERTILE:
            continue
        rec = {"unit_id": uid, "dataset": sub["dataset"].iloc[0], "n_high": int(len(hi)), "n_low": int(len(lo))}
        for c in ["score_barrier_keratin", "score_AT2", "dpt_pseudotime"]:
            if c in sub:
                rec[f"delta_{c}"] = float(np.nanmedian(hi[c]) - np.nanmedian(lo[c]))
        rows.append(rec)
    if not rows:
        return pd.DataFrame()
    frame = pd.DataFrame(rows)
    out_rows = []
    for c in ["delta_score_barrier_keratin", "delta_score_AT2", "delta_dpt_pseudotime"]:
        if c not in frame:
            continue
        vals = frame[c].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) < 4:
            continue
        try:
            w, p = stats.wilcoxon(vals, alternative="two-sided", zero_method="wilcox")
        except ValueError:
            w, p = np.nan, np.nan
        out_rows.append(
            {
                "contrast": c.replace("delta_", "paired tertile Δ "),
                "n_units": int(len(vals)),
                "median_delta": float(np.median(vals)),
                "W": float(w) if np.isfinite(w) else np.nan,
                "p": float(p) if np.isfinite(p) else np.nan,
            }
        )
    return pd.DataFrame(out_rows)


def fig_trajectory(adata, out: Path, title: str) -> None:
    umap = adata.obsm["X_umap"]
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.4), constrained_layout=True)
    panels = [
        (axes[0, 0], adata.obs["score_CLDN4"], "CLDN4", "viridis"),
        (axes[0, 1], adata.obs.get("dpt_pseudotime", pd.Series(index=adata.obs.index, dtype=float)), "DPT", "magma"),
        (axes[1, 0], adata.obs["dataset"].astype("category").cat.codes, "cohort", "tab10"),
        (axes[1, 1], adata.obs["unit_quartile"].map({"Q1": 0, "Q2": 1, "Q3": 2, "Q4": 3}).fillna(1.5), "unit CLDN4 quartile", "coolwarm"),
    ]
    for ax, val, lab, cmap in panels:
        sc = ax.scatter(umap[:, 0], umap[:, 1], c=val, s=4, cmap=cmap, linewidths=0, alpha=0.85)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(lab)
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(title, fontsize=11)
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def fig_paga(adata, out: Path) -> None:
    import scanpy as sc

    fig, ax = plt.subplots(figsize=(5.2, 4.6), constrained_layout=True)
    sc.pl.paga(adata, ax=ax, show=False, title="PAGA (Leiden)", node_size_scale=1.2)
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def fig_extra_programs(means: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.2), constrained_layout=True)
    pairs = [
        ("score_CLDN4", "dpt_pseudotime", "CLDN4 vs DPT"),
        ("score_CLDN4", "score_barrier_keratin", "CLDN4 vs barrier/keratin"),
        ("score_CLDN4", "score_AT2", "CLDN4 vs AT2"),
    ]
    for ax, (x, y, title) in zip(axes, pairs):
        if x not in means or y not in means:
            ax.set_axis_off()
            continue
        for ds, sub in means.groupby("dataset"):
            ax.scatter(sub[x], sub[y], s=28, alpha=0.85, label=ds)
        ax.set_xlabel("unit-mean CLDN4")
        ax.set_ylabel(y.replace("score_", "").replace("dpt_pseudotime", "DPT"))
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, fontsize=7)
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def write_paga_tables(adata, outdir: Path, prefix: str) -> dict:
    import scanpy as sc

    conn = adata.uns["paga"]["connectivities"]
    if hasattr(conn, "toarray"):
        conn = conn.toarray()
    clusters = sorted(adata.obs["leiden"].astype(str).unique(), key=lambda x: int(x) if x.isdigit() else x)
    paga_df = pd.DataFrame(conn, index=clusters, columns=clusters)
    paga_df.to_csv(outdir / "tables" / f"{prefix}_paga_connectivities.tsv", sep="\t")
    vert = (
        adata.obs.groupby("leiden", observed=True)
        .agg(
            n_cells=("leiden", "size"),
            mean_CLDN4=("score_CLDN4", "mean"),
            mean_AT2=("score_AT2", "mean") if "score_AT2" in adata.obs else ("score_CLDN4", "mean"),
        )
        .reset_index()
    )
    vert.to_csv(outdir / "tables" / f"{prefix}_leiden_vertices.tsv", sep="\t", index=False)
    return {"n_leiden": int(len(clusters))}


def annotate_units(adata, given: Path) -> None:
    a = pd.read_csv(given / "GSE131907_samples.tsv", sep="\t")
    b = pd.read_csv(given / "GSE205335_patients.tsv", sep="\t")
    pct = {}
    for _, row in a.iterrows():
        pct[f"GSE131907:{row['sample']}"] = float(row["mal_CLDN4_pct"]) if pd.notna(row["mal_CLDN4_pct"]) else np.nan
    for _, row in b.iterrows():
        pct[f"GSE205335:{row['patient']}"] = (
            float(row["mal_CLDN4_pct_pos"]) if pd.notna(row["mal_CLDN4_pct_pos"]) else np.nan
        )
    qmap = quartile_map(given)
    adata.obs["unit_cldn4_pct"] = adata.obs["unit_id"].map(pct)
    adata.obs["unit_quartile"] = adata.obs["unit_id"].map(qmap).fillna("NA")


def run_compartment(adata, given: Path, outdir: Path, name: str, do_dpt: bool) -> dict:
    import scanpy as sc

    annotate_units(adata, given)
    hvg, embed = score_and_embed(adata)
    annotate_units(hvg, given)
    root_info = {}
    if do_dpt:
        root, root_info = pick_malignant_root(hvg)
        hvg.uns["iroot"] = root
        sc.tl.diffmap(hvg)
        sc.tl.dpt(hvg)
        hvg.obs["dpt_pseudotime"] = hvg.obs["dpt_pseudotime"].astype(float)
    means = unit_means(hvg)
    tests = paga_tests(means) if do_dpt else pd.DataFrame()
    paired = paired_tertile(hvg) if do_dpt else pd.DataFrame()
    means.to_csv(outdir / "tables" / f"{name}_unit_means.tsv", sep="\t", index=False)
    if len(tests):
        tests.to_csv(outdir / "tables" / f"{name}_sample_level_spearman.tsv", sep="\t", index=False)
    if len(paired):
        paired.to_csv(outdir / "tables" / f"{name}_tertile_paired.tsv", sep="\t", index=False)
    write_paga_tables(hvg, outdir, name)
    fig_trajectory(hvg, outdir / "figures" / f"fig_{name}_trajectory_cldn4.png", f"{name} PAGA/UMAP · CLDN4")
    try:
        fig_paga(hvg, outdir / "figures" / f"fig_{name}_paga.png")
    except Exception as exc:  # noqa: BLE001
        print(f"PAGA plot failed: {exc}", flush=True)
    if do_dpt and len(means):
        fig_extra_programs(means, outdir / "figures" / f"fig_{name}_extra_programs.png")
    info = {
        "compartment": name,
        "n_cells": int(hvg.n_obs),
        "n_units": int(hvg.obs["unit_id"].nunique()),
        "n_units_in_means": int(len(means)),
        "embed": embed,
        "root": root_info,
        "subsampled": True,
    }
    (outdir / "tables" / f"{name}_paga_info.json").write_text(json.dumps(info, indent=2))
    return info


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--given", type=Path, default=ROOT / "data" / "given")
    p.add_argument("--outdir", type=Path, default=ROOT / "results")
    args = p.parse_args()
    import anndata as ad
    import scanpy as sc

    sc.settings.verbosity = 2
    sc.settings.figdir = str(args.outdir / "figures")
    (args.outdir / "tables").mkdir(parents=True, exist_ok=True)
    (args.outdir / "figures").mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input)
    print(f"loaded {adata.n_obs} cells × {adata.n_vars} genes", flush=True)
    mal = adata[adata.obs["compartment"].astype(str) == "malignant"].copy()
    tnk = adata[adata.obs["compartment"].astype(str) == "TNK"].copy()
    print(f"malignant={mal.n_obs} TNK={tnk.n_obs}", flush=True)
    info_m = run_compartment(mal, args.given, args.outdir, "malignant", do_dpt=True)
    info_t = run_compartment(tnk, args.given, args.outdir, "tnk", do_dpt=False)
    (args.outdir / "tables" / "paga_run_info.json").write_text(
        json.dumps({"malignant": info_m, "tnk": info_t}, indent=2)
    )
    print(json.dumps({"malignant": info_m, "tnk": info_t}, indent=2), flush=True)


if __name__ == "__main__":
    main()
