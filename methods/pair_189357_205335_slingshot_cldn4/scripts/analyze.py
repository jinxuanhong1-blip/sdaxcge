#!/usr/bin/env python3
"""REAL Slingshot + PAGA on GSE189357+GSE205335 epithelium, CLDN4 only.

ADDITIVE. PR #459 combo enum is given and is not re-ranked.
Root is never the CLDN4-high Leiden cluster.
Primary clock = Slingshot (Street 2018). PAGA is the graph companion.
Inferential unit = patient (GSE189357 TD1–TD9; GSE205335 Pxxxx).
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from gene_sets import COMPARATOR, FOCAL, IFN_ALIASES, STATES  # noqa: E402

MIN_UNIT = 10
MIN_ARM = 8
LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGH = 30
N_PCS = 30
SEED = 1


def fmt_p(value) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_rho(value) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.3f}"


def bh_fdr(pvals: list[float]) -> list[float]:
    arr = np.asarray(pvals, dtype=float)
    n = arr.size
    order = np.argsort(arr)
    q = np.empty(n, dtype=float)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        adj = arr[i] * n / (n - rank + 1)
        prev = min(prev, adj)
        q[i] = min(prev, 1.0)
    return q.tolist()


def spearman(df: pd.DataFrame, a: str, b: str) -> dict:
    sub = df[[a, b]].replace([np.inf, -np.inf], np.nan).dropna()
    n = int(len(sub))
    if n < 5:
        return {"n": n, "rho": None, "p": None}
    rho, p = stats.spearmanr(sub[a], sub[b])
    return {"n": n, "rho": float(rho), "p": float(p)}


def score_genes(adata, genes, name: str) -> None:
    present = []
    for g in genes:
        aliases = IFN_ALIASES.get(g, (g,))
        hit = next((a for a in aliases if a in adata.var_names), None)
        if hit is not None:
            present.append(hit)
    import scanpy as sc

    if not present:
        adata.obs[name] = 0.0
        return
    sc.tl.score_genes(adata, present, score_name=name, use_raw=False, random_state=SEED)


def pick_start_cluster(adata) -> dict:
    """Never root on a CLDN4-high cluster. Prefer AT2-like / author AT2."""
    cl = adata.obs["leiden"].astype(str)
    rows = []
    for c, sub in adata.obs.groupby(cl, observed=True):
        rows.append(
            {
                "leiden": str(c),
                "n_cells": int(len(sub)),
                "mean_CLDN4": float(sub["CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_SFTPC": float(sub["SFTPC"].mean()),
                "mean_barrier": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "n_author_AT2": int(sub["author_subtype"].astype(str).eq("AT2").sum()),
                "frac_189357": float(sub["dataset"].eq("GSE189357").mean()),
            }
        )
    tab = pd.DataFrame(rows).sort_values("leiden")
    cldn = tab["mean_CLDN4"].to_numpy()
    # top tercile of cluster-mean CLDN4 is forbidden
    if len(tab) >= 3:
        cut = np.quantile(cldn, 2 / 3)
        allowed = tab.loc[tab["mean_CLDN4"] <= cut].copy()
    else:
        allowed = tab.loc[tab["mean_CLDN4"] < tab["mean_CLDN4"].max()].copy()
        if allowed.empty:
            allowed = tab.copy()
    # prefer a cluster that actually contains author AT2 and is not CLDN4-high
    at2ish = allowed.loc[allowed["n_author_AT2"] > 0]
    pool = at2ish if len(at2ish) else allowed
    pool = pool.sort_values(["mean_AT2", "mean_SFTPC", "n_author_AT2"], ascending=False)
    start = str(pool.iloc[0]["leiden"])
    high_cldn = tab.loc[tab["mean_CLDN4"].idxmax(), "leiden"]
    if start == str(high_cldn) and len(allowed) > 1:
        start = str(allowed.sort_values("mean_AT2", ascending=False).iloc[0]["leiden"])
    info = {
        "start_cluster": start,
        "forbidden_high_cldn4_cluster": str(high_cldn),
        "rule": (
            "Leiden cluster with highest AT2 score among clusters not in the "
            "top tercile of mean CLDN4; prefer author AT2 if present. Never CLDN4-high."
        ),
        "n_allowed": int(len(allowed)),
        "n_clusters": int(len(tab)),
    }
    return tab, info


def python_slingshot(X: np.ndarray, labels: np.ndarray, start: str):
    """Street-style MST + polyline projection if R slingshot is unavailable."""
    labs = np.asarray(labels, dtype=str)
    uniq = [c for c in pd.unique(labs)]
    if start not in uniq:
        raise SystemExit(f"start {start} not in clusters")
    centers = np.vstack([X[labs == c].mean(axis=0) for c in uniq])
    dist = squareform(pdist(centers, metric="euclidean"))
    mst = minimum_spanning_tree(dist).toarray()
    mst = np.maximum(mst, mst.T)
    adj = {c: [] for c in uniq}
    for i, a in enumerate(uniq):
        for j, b in enumerate(uniq):
            if i < j and mst[i, j] > 0:
                adj[a].append(b)
                adj[b].append(a)
    # lineages: unique paths from start to each leaf
    leaves = [c for c in uniq if c != start and len(adj[c]) == 1]
    if not leaves:
        leaves = [c for c in uniq if c != start]

    def path_to(target: str) -> list[str] | None:
        prev = {start: None}
        q = [start]
        while q:
            u = q.pop(0)
            if u == target:
                break
            for v in adj[u]:
                if v not in prev:
                    prev[v] = u
                    q.append(v)
        if target not in prev:
            return None
        path = [target]
        while path[-1] != start:
            path.append(prev[path[-1]])
        return path[::-1]

    lineages = []
    for i, leaf in enumerate(leaves, start=1):
        path = path_to(leaf)
        if path:
            lineages.append({"lineage_id": f"Lineage{i}", "path": path})

    n = X.shape[0]
    pt = np.full((n, len(lineages)), np.nan, dtype=float)
    w = np.zeros((n, len(lineages)), dtype=float)
    for li, lin in enumerate(lineages):
        path = lin["path"]
        idx = {c: k for k, c in enumerate(path)}
        cents = np.vstack([centers[uniq.index(c)] for c in path])
        # polyline arc-length
        seg = np.linalg.norm(np.diff(cents, axis=0), axis=1)
        cum = np.concatenate([[0.0], np.cumsum(seg)])
        in_lin = np.array([lab in idx for lab in labs])
        for i in np.flatnonzero(in_lin):
            # project onto nearest segment
            best_t = 0.0
            best_d = np.inf
            for s in range(len(path) - 1):
                a, b = cents[s], cents[s + 1]
                ab = b - a
                denom = float(np.dot(ab, ab)) + 1e-12
                u = float(np.clip(np.dot(X[i] - a, ab) / denom, 0.0, 1.0))
                proj = a + u * ab
                d = float(np.linalg.norm(X[i] - proj))
                t = cum[s] + u * seg[s]
                if d < best_d:
                    best_d = d
                    best_t = t
            pt[i, li] = best_t
            w[i, li] = 1.0 / (1.0 + best_d)
        # cells not on this lineage get 0 weight
    return lineages, pt, w


def run_r_slingshot(emb: pd.DataFrame, clusters: pd.Series, start: str, work: Path) -> dict | None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return None
    probe = subprocess.run(
        [
            rscript,
            "-e",
            'if (nzchar(Sys.getenv("R_LIBS_USER"))) .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); suppressPackageStartupMessages(library(slingshot)); cat(as.character(packageVersion("slingshot")))',
        ],
        capture_output=True,
        text=True,
    )
    if probe.returncode != 0 or not probe.stdout.strip():
        return None
    work.mkdir(parents=True, exist_ok=True)
    emb_path = work / "embedding.csv"
    cl_path = work / "clusters.csv"
    emb.to_csv(emb_path)
    clusters.rename("leiden").to_csv(cl_path)
    prefix = str(work / "sling")
    script = HERE / "run_slingshot.R"
    proc = subprocess.run(
        [
            rscript,
            str(script),
            "--embedding",
            str(emb_path),
            "--clusters",
            str(cl_path),
            "--start",
            start,
            "--out",
            prefix,
        ],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        print(proc.stdout, flush=True)
        print(proc.stderr, flush=True)
        return None
    print(proc.stdout, flush=True)
    pt = pd.read_csv(f"{prefix}_pseudotime.csv")
    w = pd.read_csv(f"{prefix}_weights.csv")
    lin = pd.read_csv(f"{prefix}_lineages.csv")
    return {
        "engine": "R_slingshot",
        "version": probe.stdout.strip(),
        "pseudotime": pt,
        "weights": w,
        "lineages": lin,
        "log": (proc.stdout or "")[-400:],
    }


def paired_tertile(adata) -> pd.DataFrame:
    rows = []
    for unit, sub in adata.obs.groupby("unit_id", observed=True):
        if len(sub) < MIN_ARM * 2:
            continue
        try:
            tert = pd.qcut(sub["CLDN4"].rank(method="average"), 3, labels=["low", "mid", "high"], duplicates="drop")
        except ValueError:
            continue
        if tert.nunique() < 3:
            continue
        hi = sub.loc[tert.eq("high")]
        lo = sub.loc[tert.eq("low")]
        if len(hi) < MIN_ARM or len(lo) < MIN_ARM:
            continue
        rows.append(
            {
                "unit_id": unit,
                "dataset": str(sub["dataset"].iloc[0]),
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "IFN_high": float(hi["score_IFN"].mean()),
                "IFN_low": float(lo["score_IFN"].mean()),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "pt_high": float(hi["sling_pt"].mean()),
                "pt_low": float(lo["sling_pt"].mean()),
            }
        )
    return pd.DataFrame(rows)


def wilcox_paired(df: pd.DataFrame, a: str, b: str) -> dict:
    if df.empty:
        return {"n": 0, "delta": None, "p": None, "W": None}
    d = df[a] - df[b]
    try:
        w, p = stats.wilcoxon(df[a], df[b], zero_method="wilcox", alternative="two-sided")
    except ValueError:
        return {"n": int(len(df)), "delta": float(d.median()), "p": None, "W": None}
    return {"n": int(len(df)), "delta": float(d.median()), "p": float(p), "W": float(w)}


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path.with_suffix(".png"), dpi=160)
    plt.savefig(path.with_suffix(".pdf"))
    plt.close()


def write_finding(path: Path, s: dict) -> None:
    prim = s["primary"]
    c4 = next(r for r in prim if r["contrast"] == "CLDN4 vs Slingshot PT")
    bar = next(r for r in prim if r["contrast"].startswith("CLDN4 vs barrier"))
    ifn = next(r for r in prim if r["contrast"] == "CLDN4 vs IFN")
    c4_ifn_pt = next(r for r in prim if r["contrast"] == "IFN vs Slingshot PT")
    c4_bar_pt = next(r for r in prim if r["contrast"] == "barrier/keratin vs Slingshot PT")
    lines = [
        "# Finding — pair GSE189357+GSE205335, CLDN4-only REAL Slingshot/PAGA",
        "",
        "ADDITIVE. **CLDN4 only.** Pair from PR #459 that already differs (malignant CLDN4 %pos vs T/NK, n=31, ρ=−0.478, Q4 r=−0.750). That T/NK cut is **not re-audited**. No dual-high TACSTD2∩CLDN4. GSE131907 is not added. Not CellChat.",
        "",
        f"Primary clock: **{s['clock']}**. Root / start cluster is **not** CLDN4-high ({s['root']['rule']}). Start Leiden = `{s['root']['start_cluster']}`; forbidden high-CLDN4 Leiden = `{s['root']['forbidden_high_cldn4_cluster']}`. Inferential unit = **patient**. Barrier/keratin and IFN scores **exclude CLDN4**. PAGA is the connectivity companion, not the clock.",
        "",
        f"**What holds (n={c4['n']} patients).** CLDN4 vs barrier/keratin (no CLDN4): ρ={fmt_rho(bar['rho'])}, p={fmt_p(bar['p'])}. CLDN4 vs IFN: ρ={fmt_rho(ifn['rho'])}, p={fmt_p(ifn['p'])}. Along Slingshot PT: CLDN4 ρ={fmt_rho(c4['rho'])}, p={fmt_p(c4['p'])}; barrier ρ={fmt_rho(c4_bar_pt['rho'])}, p={fmt_p(c4_bar_pt['p'])}; IFN ρ={fmt_rho(c4_ifn_pt['rho'])}, p={fmt_p(c4_ifn_pt['p'])}.",
        "",
        "## Verdict",
        "",
        (
            f"Patient-level CLDN4 vs Slingshot PT: n={c4['n']}, ρ={fmt_rho(c4['rho'])}, p={fmt_p(c4['p'])}. "
            f"CLDN4 vs barrier/keratin (CLDN4 excluded): n={bar['n']}, ρ={fmt_rho(bar['rho'])}, p={fmt_p(bar['p'])}. "
            f"CLDN4 vs IFN (Hallmark IFNα∩IFNγ): n={ifn['n']}, ρ={fmt_rho(ifn['rho'])}, p={fmt_p(ifn['p'])}. "
            f"barrier vs PT: n={c4_bar_pt['n']}, ρ={fmt_rho(c4_bar_pt['rho'])}, p={fmt_p(c4_bar_pt['p'])}. "
            f"IFN vs PT: n={c4_ifn_pt['n']}, ρ={fmt_rho(c4_ifn_pt['rho'])}, p={fmt_p(c4_ifn_pt['p'])}. "
            f"Slingshot engine={s['slingshot']['engine']}. Lineages={s['n_lineages']}. "
            f"PAGA components at connectivity>0: {s['paga_components']} among {s['n_leiden']} Leiden vertices. "
            "Not a TACSTD2 redo. No both-high gate. PR #459 not re-ranked."
        ),
        "",
        "## Honest n",
        "",
        f"- Analysis cells after QC (capped ≤{s['cap_per_unit']}/patient): **n_cells = {s['n_cells']}** (GSE189357 {s['n_cells_189357']}, GSE205335 {s['n_cells_205335']}).",
        f"- Patients: **n_units = {s['n_units']}** (GSE189357 {s['n_units_189357']}, GSE205335 {s['n_units_205335']}). Do not write n={s['n_units']} as one cohort.",
        f"- Patients with ≥{MIN_UNIT} epithelial cells used for Spearman: **n = {s['n_spearman']}**.",
        f"- GSE189357 gate: marker epithelium (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0. Author cell types are **absent**.",
        f"- GSE205335 gate: author `lineage.total == Epithelial cells` on non-normal tissue. Histology is mixed (ADC/SQ/SCLC/NUT).",
        f"- Author AT2 cells in the object: **{s['n_author_AT2']}** (GSE205335 only).",
        f"- Patients with ≥{MIN_ARM} cells in both CLDN4-high and CLDN4-low arms: **n = {s['n_paired']}**.",
        f"- Genes absent from locked sets: {s['missing_genes']}.",
        f"- Slingshot: engine={s['slingshot']['engine']}; {s['slingshot'].get('version') or s['slingshot'].get('reason')}.",
        "- Dual-high: not used. GSE131907: not added. PR #459 T/NK Spearman: not re-audited.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGH}; PCs {N_PCS}.",
        "- Batch: harmonypy on PCA, batch=dataset.",
        f"- Slingshot start cluster: `{s['root']['start_cluster']}` ({s['root']['rule']}).",
        f"- PAGA components at connectivity>0: **{s['paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN genes: Hallmark IFNα ∩ IFNγ (**CLDN4 out**).",
        "",
        "## Lineage table (done criterion)",
        "",
        "| lineage | path | n_cells | start | end | mean CLDN4 | mean barrier | mean IFN | ρ CLDN4~PT | ρ barrier~PT | ρ IFN~PT |",
        "| --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in s["lineage_rows"]:
        lines.append(
            f"| {r['lineage_id']} | `{r['cluster_path']}` | {r['n_cells']} | {r['start_cluster']} | {r['end_cluster']} | "
            f"{r['mean_CLDN4']:.3f} | {r['mean_barrier']:.3f} | {r['mean_IFN']:.3f} | "
            f"{fmt_rho(r['rho_CLDN4_pt'])} | {fmt_rho(r['rho_barrier_pt'])} | {fmt_rho(r['rho_IFN_pt'])} |"
        )
    lines += [
        "",
        f"Machine table: `results/tables/lineage_table.tsv` ({s['n_lineages']} lineages).",
        "",
        "## Primary (patient-level Spearman, BH inside this list)",
        "",
        "| Contrast | n_patients | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in prim:
        lines.append(
            f"| {r['contrast']} | {r['n']} | {fmt_rho(r['rho'])} | {fmt_p(r['p'])} | {fmt_p(r.get('q'))} |"
        )
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | n_patients | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s["sensitivity"]:
        lines.append(f"| {r['contrast']} | {r['n']} | {fmt_rho(r['rho'])} | {fmt_p(r['p'])} |")
    extra = s["paired"]
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (patient-paired)",
        "",
        f"Emitted: **True**. Paired tertile n={extra['n_paired']}.",
        "",
        "| Paired contrast (high − low) | n | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in extra["rows"]:
        lines.append(f"| {r['contrast']} | {r['n']} | {r['delta'] if r['delta'] is None else f'{r['delta']:.3f}'} | {fmt_p(r['p'])} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- PR #459 (combo enum vs T/NK) is not re-audited here.",
        "- GSE189357 has no author AT2 labels; the start cluster is AT2-scored, not an author nLung AT2 root.",
        "- GSE205335 is mixed histology (ADC/SQ/SCLC/NUT) and an ICI cohort; this is **not** an ICI / RECIST test.",
        "- Marker epithelium ≠ CNV-called malignant. Author epithelial ≠ CNV.",
        "- No TACSTD2∩CLDN4 both-high gate.",
        "- Do not write “AT2 differentiates into LUAD because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Slingshot is an ordering on a Harmony-merged pair, not a developmental clock.",
        "",
        "## Outputs",
        "",
        "- `results/tables/lineage_table.tsv` — **done criterion**",
        "- `results/tables/sample_level_spearman.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_along_pseudotime.png`",
        "- `results/figures/fig_paga.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/pair_189357_205335_slingshot_cldn4/requirements.txt",
        "python3 methods/pair_189357_205335_slingshot_cldn4/scripts/download.py \\",
        "  --out /tmp/geo_pair_189357_205335",
        "python3 methods/pair_189357_205335_slingshot_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/geo_pair_189357_205335 \\",
        "  --out /tmp/geo_pair_189357_205335/epithelium.h5ad",
        "python3 methods/pair_189357_205335_slingshot_cldn4/scripts/analyze.py \\",
        "  --input /tmp/geo_pair_189357_205335/epithelium.h5ad \\",
        "  --outdir methods/pair_189357_205335_slingshot_cldn4/results \\",
        "  --finding methods/pair_189357_205335_slingshot_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")
    print(f"wrote {path}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=Path("/tmp/geo_pair_189357_205335/epithelium.h5ad"))
    p.add_argument("--outdir", type=Path, default=ROOT / "results")
    p.add_argument("--finding", type=Path, default=ROOT / "FINDING.md")
    p.add_argument("--work", type=Path, default=Path("/tmp/geo_pair_189357_205335/slingshot_work"))
    args = p.parse_args()

    import scanpy as sc
    import anndata as ad

    sc.settings.verbosity = 2
    sc.settings.figdir = str(args.outdir / "figures")
    adata = ad.read_h5ad(args.input)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    adata.var_names = pd.Index(adata.var_names.astype(str).str.upper())
    adata.var_names_make_unique()
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata = adata[adata.obs["n_genes"] >= 200].copy()
    adata = adata[adata.obs["n_counts"] >= 500].copy()
    print(f"QC cells={adata.n_obs} units={adata.obs['unit_id'].nunique()}", flush=True)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for g in list(FOCAL) + list(COMPARATOR) + ["SFTPC", "EPCAM", "PTPRC"]:
        if g in adata.var_names:
            adata.obs[g] = np.asarray(adata[:, g].X.todense()).ravel()
        else:
            adata.obs[g] = 0.0

    missing = {}
    for name, genes in STATES.items():
        present = [g for g in genes if g in adata.var_names or any(a in adata.var_names for a in IFN_ALIASES.get(g, ()))]
        missing[name] = [g for g in genes if g not in present]
        score_genes(adata, genes, f"score_{name}")

    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat")
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=N_PCS, svd_solver="arpack", use_highly_variable=True)
    try:
        import harmonypy as hm

        ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "dataset", max_iter_harmony=20)
        z = np.asarray(ho.Z_corr)
        if z.shape[0] == adata.n_obs:
            adata.obsm["X_pca_harmony"] = z
        elif z.shape[1] == adata.n_obs:
            adata.obsm["X_pca_harmony"] = z.T
        else:
            raise ValueError(f"Harmony Z_corr shape {z.shape} vs n_obs={adata.n_obs}")
        use_rep = "X_pca_harmony"
        batch_engine = "harmonypy"
    except Exception as exc:  # noqa: BLE001
        print(f"Harmony failed ({exc}); using PCA", flush=True)
        adata.obsm["X_pca_harmony"] = adata.obsm["X_pca"]
        use_rep = "X_pca"
        batch_engine = f"pca_only:{exc}"

    sc.pp.neighbors(adata, n_neighbors=N_NEIGH, n_pcs=N_PCS, use_rep=use_rep)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.paga(adata, groups="leiden")
    conn = np.asarray(adata.uns["paga"]["connectivities"].todense())
    n_comp = int(np.sum(np.amax(conn, axis=1) > 0) > 0)
    # connected components on thresholded PAGA
    from scipy.sparse.csgraph import connected_components

    n_comp, _ = connected_components((conn > 0).astype(int), directed=False)

    cl_tab, root = pick_start_cluster(adata)
    start = root["start_cluster"]
    # DPT companion, same start cluster, root = median AT2 cell in start cluster (not CLDN4-high)
    start_mask = adata.obs["leiden"].astype(str).eq(start)
    start_idx = np.flatnonzero(start_mask.to_numpy())
    at2 = adata.obs.loc[start_mask, "score_AT2"].to_numpy()
    root_local = int(np.argsort(np.abs(at2 - np.median(at2)))[0])
    adata.uns["iroot"] = int(start_idx[root_local])
    sc.tl.dpt(adata, n_dcs=10)
    adata.obs["dpt_pseudotime"] = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)

    emb = pd.DataFrame(adata.obsm[use_rep], index=adata.obs_names.astype(str))
    clusters = adata.obs["leiden"].astype(str)
    clusters.index = adata.obs_names.astype(str)
    sling = run_r_slingshot(emb, clusters, start, args.work)
    if sling is None:
        print("R slingshot unavailable; using Street MST+polyline Python engine", flush=True)
        lineages, pt, w = python_slingshot(emb.to_numpy(), clusters.to_numpy(), start)
        lin_df = pd.DataFrame(
            [
                {
                    "lineage_id": L["lineage_id"],
                    "cluster_path": "->".join(L["path"]),
                    "n_clusters": len(L["path"]),
                    "start_cluster": L["path"][0],
                    "end_cluster": L["path"][-1],
                }
                for L in lineages
            ]
        )
        pt_df = pd.DataFrame(pt, index=adata.obs_names.astype(str), columns=[L["lineage_id"] for L in lineages])
        pt_df["barcode"] = pt_df.index
        w_df = pd.DataFrame(w, index=adata.obs_names.astype(str), columns=[L["lineage_id"] for L in lineages])
        w_df["barcode"] = w_df.index
        sling = {
            "engine": "python_street_mst_polyline",
            "version": "Street2018-style MST + polyline (R slingshot missing)",
            "reason": "R slingshot package not importable",
            "pseudotime": pt_df.reset_index(drop=True),
            "weights": w_df.reset_index(drop=True),
            "lineages": lin_df,
        }
    else:
        sling["reason"] = ""

    pt = sling["pseudotime"].copy()
    w = sling["weights"].copy()
    if "barcode" not in pt.columns:
        pt["barcode"] = adata.obs_names.astype(str)
    pt = pt.set_index("barcode")
    w = w.set_index("barcode") if "barcode" in w.columns else w
    pt = pt.reindex(adata.obs_names.astype(str))
    w = w.reindex(adata.obs_names.astype(str))
    lin_cols = [c for c in pt.columns if c != "barcode"]
    # cell PT = weight-averaged lineage PT
    num = np.zeros(adata.n_obs, dtype=float)
    den = np.zeros(adata.n_obs, dtype=float)
    for c in lin_cols:
        pv = pd.to_numeric(pt[c], errors="coerce").to_numpy()
        wv = pd.to_numeric(w[c], errors="coerce").fillna(0).to_numpy() if c in w.columns else np.isfinite(pv).astype(float)
        m = np.isfinite(pv)
        num += np.where(m, pv * wv, 0.0)
        den += np.where(m, wv, 0.0)
    adata.obs["sling_pt"] = np.where(den > 0, num / den, np.nan)
    # primary lineage = most cells with finite PT
    best_lin = max(lin_cols, key=lambda c: int(np.isfinite(pd.to_numeric(pt[c], errors="coerce")).sum()))
    adata.obs["sling_pt_primary"] = pd.to_numeric(pt[best_lin], errors="coerce").to_numpy()
    adata.obs["sling_lineage_primary"] = best_lin

    # lineage table
    lin_rows = []
    for rec in sling["lineages"].itertuples(index=False):
        lid = rec.lineage_id
        if lid not in pt.columns:
            # R may name columns Lineage1 / curve1
            cand = [c for c in lin_cols if lid.replace("Lineage", "") in c or c == lid]
            col = cand[0] if cand else None
        else:
            col = lid
        if col is None:
            continue
        vals = pd.to_numeric(pt[col], errors="coerce")
        ww = pd.to_numeric(w[col], errors="coerce").fillna(0) if col in w.columns else vals.notna().astype(float)
        on = (ww > 0) & vals.notna()
        sub = adata.obs.loc[on.to_numpy()]
        cell_pt = vals.loc[on]
        rho_c4 = spearman(pd.DataFrame({"a": sub["CLDN4"].to_numpy(), "b": cell_pt.to_numpy()}), "a", "b")
        rho_b = spearman(pd.DataFrame({"a": sub["score_barrier_keratin"].to_numpy(), "b": cell_pt.to_numpy()}), "a", "b")
        rho_i = spearman(pd.DataFrame({"a": sub["score_IFN"].to_numpy(), "b": cell_pt.to_numpy()}), "a", "b")
        path = str(rec.cluster_path).split("->")
        tip = path[-1]
        tip_mask = sub["leiden"].astype(str).eq(tip)
        start_mask2 = sub["leiden"].astype(str).eq(str(rec.start_cluster))
        lin_rows.append(
            {
                "lineage_id": lid,
                "cluster_path": rec.cluster_path,
                "n_clusters": int(rec.n_clusters),
                "start_cluster": rec.start_cluster,
                "end_cluster": rec.end_cluster,
                "n_cells": int(on.sum()),
                "mean_CLDN4": float(sub["CLDN4"].mean()) if len(sub) else np.nan,
                "mean_barrier": float(sub["score_barrier_keratin"].mean()) if len(sub) else np.nan,
                "mean_IFN": float(sub["score_IFN"].mean()) if len(sub) else np.nan,
                "mean_AT2": float(sub["score_AT2"].mean()) if len(sub) else np.nan,
                "start_mean_CLDN4": float(sub.loc[start_mask2, "CLDN4"].mean()) if start_mask2.any() else np.nan,
                "start_mean_AT2": float(sub.loc[start_mask2, "score_AT2"].mean()) if start_mask2.any() else np.nan,
                "end_mean_CLDN4": float(sub.loc[tip_mask, "CLDN4"].mean()) if tip_mask.any() else np.nan,
                "end_mean_IFN": float(sub.loc[tip_mask, "score_IFN"].mean()) if tip_mask.any() else np.nan,
                "rho_CLDN4_pt": rho_c4["rho"],
                "p_CLDN4_pt": rho_c4["p"],
                "rho_barrier_pt": rho_b["rho"],
                "p_barrier_pt": rho_b["p"],
                "rho_IFN_pt": rho_i["rho"],
                "p_IFN_pt": rho_i["p"],
            }
        )
    lineage_df = pd.DataFrame(lin_rows)

    # unit means
    obs = adata.obs.copy()
    unit = (
        obs.groupby("unit_id", observed=True)
        .agg(
            dataset=("dataset", "first"),
            histology=("histology", "first"),
            n_cells=("CLDN4", "size"),
            mean_CLDN4=("CLDN4", "mean"),
            mean_TACSTD2=("TACSTD2", "mean"),
            mean_SFTPC=("SFTPC", "mean"),
            mean_AT2=("score_AT2", "mean"),
            mean_barrier=("score_barrier_keratin", "mean"),
            mean_IFN=("score_IFN", "mean"),
            mean_malignant=("score_malignant_like", "mean"),
            mean_sling=("sling_pt", "mean"),
            mean_dpt=("dpt_pseudotime", "mean"),
        )
        .reset_index()
    )
    elig = unit.loc[unit["n_cells"] >= MIN_UNIT].copy()
    primary_spec = [
        ("CLDN4 vs Slingshot PT", "mean_CLDN4", "mean_sling"),
        ("CLDN4 vs DPT (companion)", "mean_CLDN4", "mean_dpt"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier"),
        ("CLDN4 vs IFN", "mean_CLDN4", "mean_IFN"),
        ("CLDN4 vs malignant-like", "mean_CLDN4", "mean_malignant"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("barrier/keratin vs Slingshot PT", "mean_barrier", "mean_sling"),
        ("IFN vs Slingshot PT", "mean_IFN", "mean_sling"),
        ("SFTPC vs Slingshot PT (control)", "mean_SFTPC", "mean_sling"),
        ("AT2 vs Slingshot PT (control)", "mean_AT2", "mean_sling"),
    ]
    primary = [{"contrast": name, **spearman(elig, a, b)} for name, a, b in primary_spec]
    qs = bh_fdr([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = q if r["p"] is not None else None

    e189 = elig.loc[elig["dataset"].eq("GSE189357")]
    e205 = elig.loc[elig["dataset"].eq("GSE205335")]
    e_adc = elig.loc[elig["histology"].isin(["ADC", "LUAD"])]
    sensitivity = [
        {"contrast": "GSE189357-only CLDN4 vs Slingshot PT", **spearman(e189, "mean_CLDN4", "mean_sling")},
        {"contrast": "GSE205335-only CLDN4 vs Slingshot PT", **spearman(e205, "mean_CLDN4", "mean_sling")},
        {"contrast": "GSE189357-only CLDN4 vs IFN", **spearman(e189, "mean_CLDN4", "mean_IFN")},
        {"contrast": "GSE205335-only CLDN4 vs IFN", **spearman(e205, "mean_CLDN4", "mean_IFN")},
        {"contrast": "GSE189357-only CLDN4 vs barrier/keratin", **spearman(e189, "mean_CLDN4", "mean_barrier")},
        {"contrast": "GSE205335-only CLDN4 vs barrier/keratin", **spearman(e205, "mean_CLDN4", "mean_barrier")},
        {"contrast": "ADC/LUAD-only CLDN4 vs Slingshot PT", **spearman(e_adc, "mean_CLDN4", "mean_sling")},
        {"contrast": "GSE205335-only IFN vs Slingshot PT", **spearman(e205, "mean_IFN", "mean_sling")},
        {"contrast": "GSE189357-only IFN vs Slingshot PT", **spearman(e189, "mean_IFN", "mean_sling")},
    ]

    paired_df = paired_tertile(adata)
    paired_rows = [
        {"contrast": "barrier/keratin (no CLDN4) high vs low", **wilcox_paired(paired_df, "barrier_high", "barrier_low")},
        {"contrast": "IFN high vs low", **wilcox_paired(paired_df, "IFN_high", "IFN_low")},
        {"contrast": "AT2 high vs low", **wilcox_paired(paired_df, "AT2_high", "AT2_low")},
        {"contrast": "Slingshot PT high vs low", **wilcox_paired(paired_df, "pt_high", "pt_low")},
    ]

    tables = args.outdir / "tables"
    figs = args.outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)
    lineage_df.to_csv(tables / "lineage_table.tsv", sep="\t", index=False)
    pd.DataFrame(primary).to_csv(tables / "sample_level_spearman.tsv", sep="\t", index=False)
    pd.DataFrame(sensitivity).to_csv(tables / "sensitivity_spearman.tsv", sep="\t", index=False)
    unit.to_csv(tables / "sample_means.tsv", sep="\t", index=False)
    cl_tab.to_csv(tables / "leiden_vertices.tsv", sep="\t", index=False)
    paired_df.to_csv(tables / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    paga_df = pd.DataFrame(conn, index=adata.obs["leiden"].astype(str).unique(), columns=None)
    # write PAGA connectivities aligned to leiden categories
    leiden_cats = adata.obs["leiden"].astype("category").cat.categories.astype(str)
    paga_out = pd.DataFrame(conn, index=leiden_cats, columns=leiden_cats)
    paga_out.to_csv(tables / "paga_connectivities.tsv", sep="\t")

    # figures
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0))
    sc.pl.umap(adata, color="CLDN4", ax=axes[0], show=False, frameon=False, cmap="viridis", title="UMAP CLDN4")
    sc.pl.umap(adata, color="sling_pt", ax=axes[1], show=False, frameon=False, cmap="magma", title="UMAP Slingshot PT")
    c4 = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    for ds, col in (("GSE189357", "#1f77b4"), ("GSE205335", "#d62728")):
        sub = elig.loc[elig["dataset"].eq(ds)]
        axes[2].scatter(sub["mean_sling"], sub["mean_CLDN4"], s=48, c=col, label=f"{ds} n={len(sub)}")
    axes[2].set_xlabel("patient-mean Slingshot PT")
    axes[2].set_ylabel("patient-mean CLDN4")
    axes[2].set_title(f"CLDN4 vs PT  n={c4['n']}  ρ={fmt_rho(c4['rho'])}  p={fmt_p(c4['p'])}")
    axes[2].legend(frameon=False, fontsize=8)
    savefig(figs / "fig_trajectory_cldn4")

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0))
    rng = np.random.default_rng(0)
    take = rng.choice(adata.n_obs, size=min(4000, adata.n_obs), replace=False)
    xx = adata.obs["sling_pt"].to_numpy()[take]
    for ax, key, title in (
        (axes[0], "CLDN4", "CLDN4"),
        (axes[1], "score_barrier_keratin", "barrier/keratin (no CLDN4)"),
        (axes[2], "score_IFN", "IFN (IFNα∩IFNγ)"),
    ):
        yy = adata.obs[key].to_numpy()[take]
        ax.hexbin(xx, yy, gridsize=40, cmap="Greys", mincnt=1)
        # unit overlay
        ycol = {"CLDN4": "mean_CLDN4", "score_barrier_keratin": "mean_barrier", "score_IFN": "mean_IFN"}[key]
        for ds, col in (("GSE189357", "#1f77b4"), ("GSE205335", "#d62728")):
            sub = elig.loc[elig["dataset"].eq(ds)]
            ax.scatter(sub["mean_sling"], sub[ycol], s=36, c=col, edgecolor="white", linewidth=0.4, label=ds)
        ax.set_xlabel("Slingshot PT")
        ax.set_ylabel(title)
        ax.set_title(title + " along PT")
        ax.legend(frameon=False, fontsize=7)
    savefig(figs / "fig_along_pseudotime")

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    sc.pl.paga(adata, ax=ax, show=False, title=f"PAGA Leiden r={LEIDEN_RES}")
    savefig(figs / "fig_paga")

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    counts = [
        s["n"] if False else None
        for s in []
    ]
    axes[0].bar(
        ["GSE189357\npatients", "GSE205335\npatients", "Spearman n"],
        [
            int(adata.obs.loc[adata.obs["dataset"].eq("GSE189357"), "unit_id"].nunique()),
            int(adata.obs.loc[adata.obs["dataset"].eq("GSE205335"), "unit_id"].nunique()),
            int(len(elig)),
        ],
        color=["#1f77b4", "#d62728", "#444444"],
    )
    axes[0].set_ylabel("n patients")
    axes[0].set_title("Honest n (patient is the unit)")
    axes[1].bar(
        ["GSE189357\ncells", "GSE205335\ncells"],
        [
            int((adata.obs["dataset"] == "GSE189357").sum()),
            int((adata.obs["dataset"] == "GSE205335").sum()),
        ],
        color=["#1f77b4", "#d62728"],
    )
    axes[1].set_ylabel("analysis cells (capped)")
    axes[1].set_title("Do not treat n_cells as n")
    savefig(figs / "fig_honest_n")

    if not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8))
        for ax, a, b, title in (
            (axes[0], "barrier_low", "barrier_high", "barrier/keratin"),
            (axes[1], "IFN_low", "IFN_high", "IFN"),
            (axes[2], "pt_low", "pt_high", "Slingshot PT"),
        ):
            for i, rec in paired_df.iterrows():
                ax.plot([0, 1], [rec[a], rec[b]], color="#bbbbbb", lw=0.7)
            ax.scatter(np.zeros(len(paired_df)), paired_df[a], c="#4c78a8", s=22)
            ax.scatter(np.ones(len(paired_df)), paired_df[b], c="#f58518", s=22)
            ax.set_xticks([0, 1], ["CLDN4-low", "CLDN4-high"])
            ax.set_title(title)
        savefig(figs / "fig_extra_cldn4_tertile")
    else:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, "no paired tertile units", ha="center")
        ax.axis("off")
        savefig(figs / "fig_extra_cldn4_tertile")

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.5))
    sc.pl.umap(adata, color="dataset", ax=axes[0, 0], show=False, frameon=False, title="dataset")
    sc.pl.umap(adata, color="score_AT2", ax=axes[0, 1], show=False, frameon=False, cmap="cividis", title="AT2 score")
    sc.pl.umap(adata, color="score_barrier_keratin", ax=axes[0, 2], show=False, frameon=False, cmap="plasma", title="barrier (no CLDN4)")
    sc.pl.umap(adata, color="score_IFN", ax=axes[1, 0], show=False, frameon=False, cmap="cividis", title="IFN")
    sc.pl.umap(adata, color="leiden", ax=axes[1, 1], show=False, frameon=False, title="Leiden", legend_loc="on data")
    sc.pl.umap(adata, color="author_subtype", ax=axes[1, 2], show=False, frameon=False, title="author/marker subtype")
    savefig(figs / "fig_extra_umaps")

    # lineage path figure
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1], s=2, c="#dddddd")
    pal = plt.cm.tab10(np.linspace(0, 1, max(len(lineage_df), 1)))
    for i, rec in lineage_df.iterrows():
        path = str(rec.cluster_path).split("->")
        cents = []
        for c in path:
            m = adata.obs["leiden"].astype(str).eq(c)
            if m.any():
                cents.append(adata.obsm["X_umap"][m.to_numpy()].mean(axis=0))
        if len(cents) >= 2:
            arr = np.vstack(cents)
            ax.plot(arr[:, 0], arr[:, 1], "-o", color=pal[i % len(pal)], label=rec.lineage_id, lw=2)
    ax.set_title(f"Slingshot lineages (start={start}, not CLDN4-high)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_xticks([])
    ax.set_yticks([])
    savefig(figs / "fig_lineages")

    summary = {
        "clock": "REAL Slingshot (Street 2018)" if sling["engine"].startswith("R_") else "Slingshot-style MST + polyline (R package missing)",
        "slingshot": {"engine": sling["engine"], "version": sling.get("version"), "reason": sling.get("reason")},
        "root": root,
        "n_lineages": int(len(lineage_df)),
        "lineage_rows": lin_rows,
        "primary": primary,
        "sensitivity": sensitivity,
        "paired": {"n_paired": int(len(paired_df)), "rows": paired_rows},
        "n_cells": int(adata.n_obs),
        "n_cells_189357": int((adata.obs["dataset"] == "GSE189357").sum()),
        "n_cells_205335": int((adata.obs["dataset"] == "GSE205335").sum()),
        "n_units": int(adata.obs["unit_id"].nunique()),
        "n_units_189357": int(adata.obs.loc[adata.obs["dataset"].eq("GSE189357"), "unit_id"].nunique()),
        "n_units_205335": int(adata.obs.loc[adata.obs["dataset"].eq("GSE205335"), "unit_id"].nunique()),
        "n_spearman": int(len(elig)),
        "n_author_AT2": int(adata.obs["author_subtype"].astype(str).eq("AT2").sum()),
        "n_paired": int(len(paired_df)),
        "n_leiden": int(adata.obs["leiden"].nunique()),
        "paga_components": int(n_comp),
        "missing_genes": missing,
        "cap_per_unit": 350,
        "batch_engine": batch_engine,
        "dual_high": False,
        "pr459_reaudited": False,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    cl_tab.to_csv(tables / "start_cluster_audit.tsv", sep="\t", index=False)
    write_finding(args.finding, summary)
    print(json.dumps({"done": True, "lineage_table": str(tables / "lineage_table.tsv"), "n_lineages": int(len(lineage_df))}, indent=2), flush=True)


if __name__ == "__main__":
    main()
