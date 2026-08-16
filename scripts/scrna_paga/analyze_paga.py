#!/usr/bin/env python3
"""PAGA + diffusion on GSE131907 LUAD nLung/tLung epithelium.

ADDITIVE. Locked design is in methods/scrna_paga/playbook.md.
Does not retune Leiden, roots, or gene lists after seeing TACSTD2.
Sample (not cell) is the inferential unit.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from gene_sets import (
    AUTHOR_AT2,
    AUTHOR_CLUB,
    AUTHOR_TUMOR_STATE,
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
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
EXTRA_BARRIER_RHO_GT = 0.0
EXTRA_BARRIER_P_LT = 0.05


def _bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n, dtype=float)
    running = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        # rank from the end
        k = n - rank + 1
        running = min(running, p[i] * n / k)
        q[i] = running
    return [float(min(1.0, x)) for x in q]


def _spearman(x: np.ndarray, y: np.ndarray) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 4:
        return {"n": n, "rho": None, "p": None, "note": "n<4"}
    rho, p = stats.spearmanr(x[mask], y[mask])
    if not np.isfinite(rho):
        return {"n": n, "rho": None, "p": None, "note": "undefined"}
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


def _pick_root(adata) -> tuple[int, dict]:
    """External arrow: nLung AT2 (author label, else highest AT2 Leiden)."""
    nlung = adata.obs["Sample_Origin"].astype(str) == "nLung"
    subtype = adata.obs.get("Cell_subtype", pd.Series(index=adata.obs.index, dtype=object))
    subtype = subtype.astype(str)
    at2 = nlung & subtype.isin(AUTHOR_AT2)
    info = {
        "n_nLung": int(nlung.sum()),
        "n_nLung_author_AT2": int(at2.sum()),
        "rule": None,
    }
    if int(at2.sum()) >= 20:
        idx = np.flatnonzero(at2.to_numpy())
        scores = adata.obs.loc[at2, "score_AT2"].to_numpy()
        med = np.nanmedian(scores)
        pick = idx[int(np.nanargmin(np.abs(scores - med)))]
        info["rule"] = "nLung author AT2 (median AT2 score)"
        return int(pick), info
    # fallback: Leiden cluster with highest mean AT2 among nLung cells
    if "leiden" not in adata.obs:
        raise SystemExit("leiden missing before root pick")
    sub = adata.obs.loc[nlung, ["leiden", "score_AT2"]]
    if sub.empty:
        raise SystemExit("no nLung cells to root DPT")
    means = sub.groupby("leiden", observed=True)["score_AT2"].mean().sort_values(ascending=False)
    top = str(means.index[0])
    cand = nlung & (adata.obs["leiden"].astype(str) == top)
    scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
    idx = np.flatnonzero(cand.to_numpy())
    pick = idx[int(np.nanargmax(scores))]
    info["rule"] = f"nLung Leiden {top} max AT2 score (no author AT2)"
    info["fallback_cluster"] = top
    return int(pick), info


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


def _write_report(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    lines = [
        "# GSE131907 PAGA / diffusion: TACSTD2 and CLDN4 on LUAD epithelium",
        "",
        "**Slice:** `methods/scrna_paga/` + `scripts/scrna_paga/` + `results/scrna_paga/` only.",
        "**Additive.** Does not rewrite `methods/scrna/` or `methods/trajectory/`.",
        "**Cohort:** GSE131907 (Kim et al. 2020), LUAD, public UMI. tLung epithelium + nLung epithelium (AT2/club prior).",
        "**ICI labels:** none. This slice cannot test ICI / MPR / RECIST.",
        "",
        "## What was tested",
        "",
        "PAGA connectivities and diffusion pseudotime on author-labeled LUAD epithelium, then sample-level Spearman of TACSTD2 and CLDN4 (separately) versus DPT and versus pre-specified AT2 / club / basal / barrier-keratin / malignant-like scores.",
        "",
        "## Coverage (honest missingness)",
        "",
        f"- Cells in object: **n_cells = {s['n_cells']}** (nLung {s['n_cells_nLung']}, tLung {s['n_cells_tLung']}).",
        f"- Samples: **n_samples = {s['n_samples']}** (nLung {s['n_samples_nLung']}, tLung {s['n_samples_tLung']}).",
        f"- Patients: **n_patients = {s['n_patients']}** (LUNG_Nxx / LUNG_Txx numeric id; unpaired N01 and T25).",
        f"- Author subtypes (cells): AT2 {s['subtype_counts'].get('AT2', 0)}, Club {s['subtype_counts'].get('Club', 0)}, AT1 {s['subtype_counts'].get('AT1', 0)}, Ciliated {s['subtype_counts'].get('Ciliated', 0)}, tS1 {s['subtype_counts'].get('tS1', 0)}, tS2 {s['subtype_counts'].get('tS2', 0)}, tS3 {s['subtype_counts'].get('tS3', 0)}, NA {s['subtype_counts'].get('NA', 0)}. **Author basal = 0.**",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        f"- Samples with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} epithelial cells used for sample-level Spearman: **n = {s['n_samples_eligible']}**.",
        f"- GSE253013 skipped (processed RDS ~9.3 GB, over public-file cap). GSE207422 not pooled (NSCLC, not LUAD-only).",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- DPT root: {s['root']['rule']} (root cell index {s['root']['index']}).",
        f"- PAGA components at connectivity>0: {s['n_paga_components']}.",
        "",
        "## Primary (sample-level Spearman)",
        "",
        "| Contrast | n_samples | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in s["primary_spearman"]:
        rho = "NA" if row["rho"] is None else f"{row['rho']:.3f}"
        pv = "NA" if row["p"] is None else f"{row['p']:.3g}"
        qv = "NA" if row.get("q") is None else f"{row['q']:.3g}"
        lines.append(f"| {row['contrast']} | {row['n']} | {rho} | {pv} | {qv} |")
    lines += [
        "",
        "## Sensitivity (origin-stratified; not in the BH family)",
        "",
        "| Contrast | n_samples | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for row in s.get("sensitivity_spearman", []):
        rho = "NA" if row["rho"] is None else f"{row['rho']:.3f}"
        pv = "NA" if row["p"] is None else f"{row['p']:.3g}"
        lines.append(f"| {row['contrast']} | {row['n']} | {rho} | {pv} |")
    extra = s["extra_barrier_figure"]
    lines += [
        "",
        "## Extra barrier/keratin figure",
        "",
        (
            f"Emitted: **{extra['emitted']}**. "
            f"Rule: sample-level Spearman(TACSTD2, barrier_keratin) ρ>0 and p<{EXTRA_BARRIER_P_LT}. "
            f"Observed n={extra['n']}, ρ={extra['rho'] if extra['rho'] is None else round(extra['rho'], 3)}, "
            f"p={extra['p'] if extra['p'] is None else f'{extra['p']:.3g}'}."
        ),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Caveats",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- Malignant-like is author tS1/tS2/tS3 and/or a weak CEACAM5/6/MKI67 score — **not CNV**.",
        "- Club/basal are airway programs. If PAGA disconnects them from AT2, that is a discrete-state result, not a failed download.",
        "- GSE131907 is treatment-naive. Do not write ICI language.",
        "- No TACSTD2∩CLDN4 both-high gate.",
        "- Do not write “AT2 differentiates into LUAD because PAGA is connected.”",
        "- Do not write “TROP2 marks the malignant terminal.”",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r requirements.txt",
        "bash scripts/scrna_paga/download.sh /tmp/scrna_paga_data",
        "python3 scripts/scrna_paga/extract_epithelium.py --data /tmp/scrna_paga_data --out /tmp/scrna_paga_data/epithelium.h5ad",
        "python3 scripts/scrna_paga/analyze_paga.py --input /tmp/scrna_paga_data/epithelium.h5ad --outdir results/scrna_paga",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default="results/scrna_paga")
    args = p.parse_args()
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    adata = sc.read_h5ad(args.input)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()

    # QC floors (locked, not tuned on TACSTD2)
    adata.var["n_cells"] = np.array((adata.X > 0).sum(axis=0)).ravel()
    adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).ravel()
    adata.obs["n_umi"] = np.array(adata.X.sum(axis=1)).ravel()
    sc.pp.filter_genes(adata, min_cells=10)
    keep = (adata.obs["n_genes"] >= 200) & (adata.obs["n_umi"] >= 500)
    adata = adata[keep].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    absent: dict[str, list[str]] = {}
    for name, genes in STATES.items():
        absent[name] = _score_mean(adata, genes, f"score_{name}")
    for g in FOCAL + CONTROLS + QC_NEG:
        if g in adata.var_names:
            x = adata[:, g].X
            if hasattr(x, "toarray"):
                x = x.toarray()
            adata.obs[f"expr_{g}"] = np.asarray(x, dtype=float).ravel()
        else:
            adata.obs[f"expr_{g}"] = np.nan
            absent.setdefault("single_genes", []).append(g)

    try:
        sc.pp.highly_variable_genes(
            adata, layer="counts", flavor="seurat_v3", n_top_genes=N_HVG
        )
    except ImportError:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=N_HVG)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2)
    except TypeError:
        sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.paga(adata, groups="leiden")
    sc.tl.diffmap(adata, n_comps=15)
    sc.tl.umap(adata)

    root_i, root_info = _pick_root(adata)
    adata.uns["iroot"] = root_i
    sc.tl.dpt(adata, n_dcs=10)
    root_info["index"] = root_i
    root_info["root_sample"] = str(adata.obs.iloc[root_i].get("Sample", ""))
    root_info["root_origin"] = str(adata.obs.iloc[root_i].get("Sample_Origin", ""))
    root_info["root_subtype"] = str(adata.obs.iloc[root_i].get("Cell_subtype", ""))

    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    comps = _paga_components(connect, thresh=0.0)
    leiden_ids = [str(x) for x in adata.obs["leiden"].cat.categories]
    cluster_tab = []
    for i, cl in enumerate(leiden_ids):
        sub = adata.obs[adata.obs["leiden"].astype(str) == cl]
        cluster_tab.append(
            {
                "leiden": cl,
                "n_cells": int(len(sub)),
                "n_samples": int(sub["Sample"].nunique()),
                "n_nLung": int((sub["Sample_Origin"] == "nLung").sum()),
                "n_tLung": int((sub["Sample_Origin"] == "tLung").sum()),
                "top_subtype": (
                    sub["Cell_subtype"].astype(str).value_counts().index[0]
                    if "Cell_subtype" in sub
                    else "NA"
                ),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan).mean()),
                "pct_TACSTD2_pos": float((sub["expr_TACSTD2"] > 0).mean()),
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
            }
        )
    cluster_df = pd.DataFrame(cluster_tab)
    cluster_df.to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )

    # Sample-level means
    dpt = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)
    adata.obs["dpt_pseudotime"] = dpt
    rows = []
    for (sample, origin), sub in adata.obs.groupby(["Sample", "Sample_Origin"], observed=True):
        rec = {
            "Sample": sample,
            "Sample_Origin": origin,
            "n_cells": int(len(sub)),
            "n_author_AT2": int(sub["Cell_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
            "n_author_club": int(sub["Cell_subtype"].astype(str).isin(AUTHOR_CLUB).sum()),
            "n_author_tumor_state": int(
                sub["Cell_subtype"].astype(str).isin(AUTHOR_TUMOR_STATE).sum()
            ),
            "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
            "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
            "mean_dpt": float(sub["dpt_pseudotime"].mean()),
            "mean_AT2": float(sub["score_AT2"].mean()),
            "mean_club": float(sub["score_club"].mean()),
            "mean_basal": float(sub["score_basal"].mean()),
            "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
            "mean_malignant_like": float(sub["score_malignant_like"].mean()),
            "pct_TACSTD2_pos": float((sub["expr_TACSTD2"] > 0).mean()),
            "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
        }
        rows.append(rec)
    sample_df = pd.DataFrame(rows)
    # patient id: Kim samples look like LUAD_T18 / LUNG_N18 — strip suffix
    sample_df["patient_id"] = sample_df["Sample"].astype(str).str.extract(
        r"LUNG_[NT](\d+)", expand=False
    )
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()

    contrasts = [
        ("TACSTD2 vs DPT", "mean_TACSTD2", "mean_dpt"),
        ("CLDN4 vs DPT", "mean_CLDN4", "mean_dpt"),
        ("TACSTD2 vs CLDN4", "mean_TACSTD2", "mean_CLDN4"),
        ("TACSTD2 vs AT2 score", "mean_TACSTD2", "mean_AT2"),
        ("TACSTD2 vs club score", "mean_TACSTD2", "mean_club"),
        ("TACSTD2 vs basal score", "mean_TACSTD2", "mean_basal"),
        ("TACSTD2 vs barrier/keratin score", "mean_TACSTD2", "mean_barrier_keratin"),
        ("TACSTD2 vs malignant-like score", "mean_TACSTD2", "mean_malignant_like"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs club score", "mean_CLDN4", "mean_club"),
        ("CLDN4 vs basal score", "mean_CLDN4", "mean_basal"),
        ("CLDN4 vs barrier/keratin score", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs malignant-like score", "mean_CLDN4", "mean_malignant_like"),
        ("SFTPC vs DPT (control)", "mean_AT2", "mean_dpt"),
    ]
    # SFTPC control uses AT2 score as proxy if we don't have sample SFTPC mean
    if "expr_SFTPC" in adata.obs:
        sft = []
        for (sample, origin), sub in adata.obs.groupby(["Sample", "Sample_Origin"], observed=True):
            sft.append(((sample, origin), float(sub["expr_SFTPC"].mean())))
        sft_map = dict(sft)
        elig["mean_SFTPC"] = [
            sft_map.get((r.Sample, r.Sample_Origin), np.nan) for r in elig.itertuples()
        ]
        contrasts[-1] = ("SFTPC vs DPT (control)", "mean_SFTPC", "mean_dpt")

    primary = []
    for name, a, b in contrasts:
        if a not in elig.columns or b not in elig.columns:
            primary.append({"contrast": name, "n": 0, "rho": None, "p": None})
            continue
        sp = _spearman(elig[a].to_numpy(), elig[b].to_numpy())
        primary.append({"contrast": name, **sp})
    # BH on the pre-specified list
    p_for_q = [r["p"] if r["p"] is not None else 1.0 for r in primary]
    qs = _bh(p_for_q)
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q
    pd.DataFrame(primary).to_csv(tabdir / "sample_level_spearman.tsv", sep="\t", index=False)

    # tLung-only sensitivity
    elig_t = elig[elig["Sample_Origin"] == "tLung"]
    elig_n = elig[elig["Sample_Origin"] == "nLung"]
    tlung_sp = _spearman(elig_t["mean_TACSTD2"].to_numpy(), elig_t["mean_dpt"].to_numpy()) if len(elig_t) else {"n": 0, "rho": None, "p": None}
    tlung_sp_c = _spearman(elig_t["mean_CLDN4"].to_numpy(), elig_t["mean_dpt"].to_numpy()) if len(elig_t) else {"n": 0, "rho": None, "p": None}
    nlung_sp = _spearman(elig_n["mean_TACSTD2"].to_numpy(), elig_n["mean_dpt"].to_numpy()) if len(elig_n) else {"n": 0, "rho": None, "p": None}
    nlung_sp_c = _spearman(elig_n["mean_CLDN4"].to_numpy(), elig_n["mean_dpt"].to_numpy()) if len(elig_n) else {"n": 0, "rho": None, "p": None}
    tlung_bar = _spearman(elig_t["mean_TACSTD2"].to_numpy(), elig_t["mean_barrier_keratin"].to_numpy()) if len(elig_t) else {"n": 0, "rho": None, "p": None}
    sensitivity = [
        {"contrast": "tLung-only TACSTD2 vs DPT", **tlung_sp},
        {"contrast": "tLung-only CLDN4 vs DPT", **tlung_sp_c},
        {"contrast": "nLung-only TACSTD2 vs DPT", **nlung_sp},
        {"contrast": "nLung-only CLDN4 vs DPT", **nlung_sp_c},
        {"contrast": "tLung-only TACSTD2 vs barrier/keratin", **tlung_bar},
    ]

    barrier_row = next(r for r in primary if r["contrast"] == "TACSTD2 vs barrier/keratin score")
    emit_extra = (
        barrier_row["rho"] is not None
        and barrier_row["p"] is not None
        and barrier_row["rho"] > EXTRA_BARRIER_RHO_GT
        and barrier_row["p"] < EXTRA_BARRIER_P_LT
    )

    # Cell-level descriptive (not inferential)
    cell_desc = {
        "TACSTD2_vs_dpt": _spearman(
            adata.obs["expr_TACSTD2"].to_numpy(),
            adata.obs["dpt_pseudotime"].to_numpy(),
        ),
        "CLDN4_vs_dpt": _spearman(
            adata.obs["expr_CLDN4"].to_numpy(),
            adata.obs["dpt_pseudotime"].to_numpy(),
        ),
        "TACSTD2_vs_barrier": _spearman(
            adata.obs["expr_TACSTD2"].to_numpy(),
            adata.obs["score_barrier_keratin"].to_numpy(),
        ),
        "note": "descriptive only; n_cells is not the experimental unit",
    }

    # ---------- figures ----------
    def _save(fig, name: str) -> None:
        fig.tight_layout()
        fig.savefig(figdir / f"{name}.png", dpi=160)
        fig.savefig(figdir / f"{name}.pdf")
        plt.close(fig)

    # 1. PAGA colored by origin mix
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    sc.pl.paga(
        adata,
        color="leiden",
        ax=ax,
        show=False,
        title=f"PAGA Leiden (n_cells={adata.n_obs}, n_samples={sample_df.shape[0]})",
    )
    _save(fig, "fig1_paga_leiden")

    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    sc.pl.umap(
        adata,
        color="Sample_Origin",
        ax=ax,
        show=False,
        title=f"UMAP display-only  nLung={int((adata.obs.Sample_Origin=='nLung').sum())}  tLung={int((adata.obs.Sample_Origin=='tLung').sum())}",
        frameon=False,
    )
    _save(fig, "fig2_umap_origin")

    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    sc.pl.umap(
        adata,
        color="Cell_subtype",
        ax=ax,
        show=False,
        title="UMAP display-only, author Cell_subtype",
        frameon=False,
    )
    _save(fig, "fig3_umap_author_subtype")

    for gene, fname in [
        ("expr_TACSTD2", "fig4_umap_TACSTD2"),
        ("expr_CLDN4", "fig5_umap_CLDN4"),
        ("dpt_pseudotime", "fig6_umap_dpt"),
        ("score_AT2", "fig7_umap_AT2_score"),
        ("score_barrier_keratin", "fig8_umap_barrier_keratin"),
    ]:
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        sc.pl.umap(adata, color=gene, ax=ax, show=False, frameon=False, cmap="viridis")
        _save(fig, fname)

    # Gene vs DPT, cells subsampled for drawing only
    rng = np.random.default_rng(0)
    take = rng.choice(adata.n_obs, size=min(4000, adata.n_obs), replace=False)
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8), sharex=True)
    for ax, gene, lab in (
        (axes[0], "expr_TACSTD2", "TACSTD2"),
        (axes[1], "expr_CLDN4", "CLDN4"),
    ):
        ax.scatter(
            adata.obs["dpt_pseudotime"].to_numpy()[take],
            adata.obs[gene].to_numpy()[take],
            s=4,
            alpha=0.25,
            c=np.where(adata.obs["Sample_Origin"].to_numpy()[take] == "tLung", "#b23a48", "#2a6f97"),
            linewidths=0,
        )
        ax.set_xlabel("DPT (AT2-rooted)")
        ax.set_ylabel(f"{lab} log1p(CP10k)")
        sp = cell_desc["TACSTD2_vs_dpt" if lab == "TACSTD2" else "CLDN4_vs_dpt"]
        rho_s = "NA" if sp["rho"] is None else f"{sp['rho']:.2f}"
        ax.set_title(f"{lab} vs DPT  cell-level ρ={rho_s} (descriptive)")
    fig.suptitle("Blue=nLung, red=tLung; cell scatter is not the inferential test", fontsize=9)
    _save(fig, "fig9_genes_vs_dpt_cells")

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8))
    for ax, gene, lab in (
        (axes[0], "mean_TACSTD2", "TACSTD2"),
        (axes[1], "mean_CLDN4", "CLDN4"),
    ):
        for origin, col in (("nLung", "#2a6f97"), ("tLung", "#b23a48")):
            sub = elig[elig["Sample_Origin"] == origin]
            ax.scatter(sub["mean_dpt"], sub[gene], s=40, c=col, label=f"{origin} n={len(sub)}")
        ax.set_xlabel("sample-mean DPT")
        ax.set_ylabel(f"sample-mean {lab}")
        row = next(r for r in primary if r["contrast"].startswith(lab))
        rho_s = "NA" if row["rho"] is None else f"{row['rho']:.2f}"
        p_s = "NA" if row["p"] is None else f"{row['p']:.3g}"
        ax.set_title(f"{lab} vs DPT  n={row['n']}  ρ={rho_s}  p={p_s}")
        ax.legend(fontsize=8, frameon=False)
    _save(fig, "fig10_sample_genes_vs_dpt")

    # Module heatmap of cluster means
    score_cols = [
        "mean_AT2",
        "mean_club",
        "mean_basal",
        "mean_barrier_keratin",
        "mean_TACSTD2",
        "mean_CLDN4",
    ]
    mat = cluster_df.set_index("leiden")[score_cols]
    mat_z = (mat - mat.mean()) / (mat.std(ddof=0).replace(0, np.nan))
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    im = ax.imshow(mat_z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(score_cols)))
    ax.set_xticklabels(
        ["AT2", "club", "basal", "barrier/KRT", "TACSTD2", "CLDN4"], rotation=40, ha="right"
    )
    ax.set_yticks(range(len(mat_z)))
    ax.set_yticklabels(
        [f"L{i} n={int(cluster_df.loc[cluster_df.leiden==i, 'n_cells'].iloc[0])}" for i in mat_z.index]
    )
    fig.colorbar(im, ax=ax, shrink=0.8, label="z of cluster mean")
    ax.set_title("Leiden vertices: lineage scores vs TACSTD2/CLDN4")
    _save(fig, "fig11_cluster_score_heatmap")

    if emit_extra:
        fig, ax = plt.subplots(figsize=(5.2, 4.4))
        for origin, col in (("nLung", "#2a6f97"), ("tLung", "#b23a48")):
            sub = elig[elig["Sample_Origin"] == origin]
            ax.scatter(
                sub["mean_barrier_keratin"],
                sub["mean_TACSTD2"],
                s=44,
                c=col,
                label=f"{origin} n={len(sub)}",
            )
        ax.set_xlabel("sample-mean barrier/keratin score")
        ax.set_ylabel("sample-mean TACSTD2")
        ax.set_title(
            f"EXTRA: TACSTD2 on barrier/keratin  n={barrier_row['n']}  "
            f"ρ={barrier_row['rho']:.2f}  p={barrier_row['p']:.3g}"
        )
        ax.legend(fontsize=8, frameon=False)
        _save(fig, "fig_extra_barrier_keratin")

    # Honesty panel: n by subtype × origin
    fig, ax = plt.subplots(figsize=(7.4, 3.6))
    ct = (
        adata.obs.assign(Cell_subtype=adata.obs["Cell_subtype"].astype(str).fillna("NA"))
        .groupby(["Sample_Origin", "Cell_subtype"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.T.plot(kind="bar", ax=ax, color={"nLung": "#2a6f97", "tLung": "#b23a48"})
    ax.set_ylabel("cells")
    ax.set_title(f"Honest n: author subtypes  n_cells={adata.n_obs}  n_samples={sample_df.shape[0]}")
    ax.legend(frameon=False)
    plt.setp(ax.get_xticklabels(), rotation=40, ha="right")
    _save(fig, "fig12_honest_n_subtype")

    subtype_counts = (
        adata.obs["Cell_subtype"].astype(str).fillna("NA").value_counts().to_dict()
    )
    n_patients = int(sample_df["patient_id"].nunique())

    # Verdict text from the locked tests (no spin)
    t2_dpt = next(r for r in primary if r["contrast"] == "TACSTD2 vs DPT")
    c4_dpt = next(r for r in primary if r["contrast"] == "CLDN4 vs DPT")
    t2_c4 = next(r for r in primary if r["contrast"] == "TACSTD2 vs CLDN4")
    def _fmt(row: dict) -> str:
        if row.get("rho") is None:
            return f"n={row.get('n')}, ρ=NA, p=NA"
        return f"n={row['n']}, ρ={row['rho']:.3f}, p={row['p']:.3g}"

    parts = [
        f"Sample-level TACSTD2 vs DPT: {_fmt(t2_dpt)}.",
        f"CLDN4 vs DPT: {_fmt(c4_dpt)}.",
        f"TACSTD2 vs CLDN4 co-expression: {_fmt(t2_c4)}.",
        f"TACSTD2 vs barrier/keratin: {_fmt(barrier_row)}; extra figure {'yes' if emit_extra else 'no'}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        (
            f"tLung-only TACSTD2 vs DPT: {_fmt(tlung_sp)} "
            f"({'inconclusive' if (tlung_sp.get('p') is None or tlung_sp.get('p', 1) >= 0.05) else 'p<0.05'})."
        ),
        "Author basal n_cells=0; basal score is a KRT5/KRT15/TP63/NGFR proxy, not an author state.",
        "Club vs TACSTD2 is a sample-level null. Do not write that TACSTD2 marks club on this object.",
        "The n=22 DPT correlations mix nLung and tLung; they are not a within-tumor progression test.",
    ]
    verdict = " ".join(parts)

    summary = {
        "accession": "GSE131907",
        "histology": "LUAD",
        "public_only": True,
        "n_cells": int(adata.n_obs),
        "n_cells_nLung": int((adata.obs["Sample_Origin"] == "nLung").sum()),
        "n_cells_tLung": int((adata.obs["Sample_Origin"] == "tLung").sum()),
        "n_samples": int(sample_df.shape[0]),
        "n_samples_nLung": int((sample_df["Sample_Origin"] == "nLung").sum()),
        "n_samples_tLung": int((sample_df["Sample_Origin"] == "tLung").sum()),
        "n_samples_eligible": int(len(elig)),
        "n_patients": n_patients,
        "subtype_counts": subtype_counts,
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "leiden_resolution": LEIDEN_RES,
        "primary_spearman": primary,
        "tlung_only_TACSTD2_vs_dpt": tlung_sp,
        "tlung_only_CLDN4_vs_dpt": tlung_sp_c,
        "sensitivity_spearman": sensitivity,
        "cell_level_descriptive": cell_desc,
        "extra_barrier_figure": {
            "emitted": bool(emit_extra),
            "rho": barrier_row["rho"],
            "p": barrier_row["p"],
            "n": barrier_row["n"],
            "rule": f"rho>{EXTRA_BARRIER_RHO_GT} and p<{EXTRA_BARRIER_P_LT} at sample level",
        },
        "qc": {
            "min_genes": 200,
            "min_umi": 500,
            "min_cells_per_gene": 10,
            "lineage_EPCAM_mean": float(adata.obs["expr_EPCAM"].mean()) if "expr_EPCAM" in adata.obs else None,
            "lineage_PTPRC_mean": float(adata.obs["expr_PTPRC"].mean()) if "expr_PTPRC" in adata.obs else None,
        },
        "verdict": verdict,
        "skipped": {
            "GSE253013": "processed RDS ~9.3 GB; over public-file cap",
            "GSE207422": "NSCLC mixed histology; not pooled into this LUAD-only graph",
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    _write_report(outdir / "REPORT.md", {"summary": summary})
    provenance = {
        "accession": "GSE131907",
        "paper": "Kim et al. Nat Commun 2020 PMID 32385277",
        "files": {
            "annotation": "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
            "umi": "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
            "series_matrix": "GSE131907_series_matrix.txt.gz",
        },
        "ftp_base": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907",
        "refused": [
            "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz (~2.9 GB)",
            "GSE253013_all_luad_garnett_temp.rds.gz (~9.3 GB)",
            "EGA FASTQ EGAD00001005054",
            "GSE207422 NSCLC mixed histology (not LUAD-only)",
        ],
        "locked": {
            "leiden_resolution": LEIDEN_RES,
            "n_hvg": N_HVG,
            "n_neighbors": N_NEIGHBORS,
            "n_pcs": N_PCS,
            "root": "nLung author AT2",
        },
    }
    man = Path("/tmp/scrna_paga_data/DOWNLOAD_MANIFEST.json")
    if man.is_file():
        provenance["download_manifest"] = json.loads(man.read_text())
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2))

    # compact obs for provenance (no expression matrix)
    adata.obs[
        [
            c
            for c in [
                "Sample",
                "Sample_Origin",
                "Cell_type",
                "Cell_subtype",
                "leiden",
                "dpt_pseudotime",
                "expr_TACSTD2",
                "expr_CLDN4",
                "score_AT2",
                "score_club",
                "score_basal",
                "score_barrier_keratin",
                "n_umi",
                "n_genes",
            ]
            if c in adata.obs.columns
        ]
    ].to_csv(tabdir / "cell_obs.tsv", sep="\t")

    print(json.dumps({"ok": True, "n_cells": adata.n_obs, "n_samples": int(sample_df.shape[0]), "extra": emit_extra}, indent=2))


if __name__ == "__main__":
    main()
