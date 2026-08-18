#!/usr/bin/env python3
"""REAL Slingshot + PAGA on GSE131907+GSE189357 epithelium, CLDN4 only.

ADDITIVE. PR #459 pair that already differs (T/NK not re-audited).
Root = GSE131907 nLung author AT2, never CLDN4-high.
Scores CLDN4 + barrier (no CLDN4) + IFN along Slingshot pseudotime.
No dual-high. No GSE148071. Done when the lineage table exists.
"""
from __future__ import annotations

import argparse
import json
import os
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import (  # noqa: E402
    AUTHOR_AT2,
    AUTHOR_CLUB,
    AUTHOR_TUMOR_STATE,
    COMPARATOR,
    CONTROLS,
    FOCAL,
    IFN_CORE,
    QC_NEG,
    STATES,
    hallmark_ifn_union,
)

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e

HERE = Path(__file__).resolve().parents[1]
LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
R_LIBS_USER = str(Path.home() / "R" / "library")


def slingshot_status() -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return {"available": False, "reason": "Rscript not on PATH"}
    env = os.environ.copy()
    env["R_LIBS_USER"] = R_LIBS_USER
    try:
        proc = subprocess.run(
            [rscript, "-e", 'cat(as.character(packageVersion("slingshot")))'],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
            env=env,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"available": False, "reason": f"Rscript slingshot probe failed: {exc}"}
    if proc.returncode != 0:
        return {
            "available": False,
            "reason": (proc.stderr or proc.stdout or "slingshot not installed").strip()[:400],
        }
    return {"available": True, "version": (proc.stdout or "").strip()}


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


def _spearman(x: np.ndarray, y: np.ndarray) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 4:
        return {"n": n, "rho": None, "p": None, "note": "n<4"}
    rho, p = stats.spearmanr(x[mask], y[mask])
    if not np.isfinite(rho):
        return {"n": n, "rho": None, "p": None, "note": "undefined"}
    return {"n": n, "rho": float(rho), "p": float(p)}


def _wilcoxon_paired(a: np.ndarray, b: np.ndarray) -> dict:
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 4:
        return {"n": n, "W": None, "p": None, "delta_median": None, "note": "n<4"}
    aa, bb = a[mask], b[mask]
    try:
        w, p = stats.wilcoxon(aa, bb, alternative="two-sided", zero_method="wilcox")
    except ValueError:
        return {
            "n": n,
            "W": None,
            "p": None,
            "delta_median": float(np.median(aa - bb)),
            "note": "wilcoxon failed",
        }
    return {
        "n": n,
        "W": float(w),
        "p": float(p),
        "delta_median": float(np.median(aa - bb)),
        "median_high": float(np.median(aa)),
        "median_low": float(np.median(bb)),
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


def _fmt(row: dict, keys=("rho", "p")) -> str:
    n = row.get("n")
    if row.get("rho") is None and "rho" in keys:
        return f"n={n}, ρ=NA, p=NA"
    if "W" in keys:
        if row.get("W") is None:
            return f"n={n}, W=NA, p=NA"
        return (
            f"n={n}, W={row['W']:.1f}, Δmed={row.get('delta_median'):.3f}, "
            f"p={row['p']:.3g}"
        )
    return f"n={n}, ρ={row['rho']:.3f}, p={row['p']:.3g}"


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def maybe_harmony(adata) -> dict:
    info = {"used": False, "reason": None}
    try:
        import harmonypy
    except ImportError:
        info["reason"] = "harmonypy not installed; neighbors on PCA"
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
        return info
    ho = harmonypy.run_harmony(
        adata.obsm["X_pca"][:, :N_PCS],
        adata.obs,
        "dataset",
        max_iter_harmony=20,
    )
    z = np.asarray(ho.Z_corr)
    if z.shape[0] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = z
    elif z.shape[1] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = z.T
    else:
        raise ValueError(f"Harmony Z_corr shape {z.shape} vs n_obs={adata.n_obs}")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_pca_harmony")
    info["used"] = True
    info["reason"] = "harmonypy on PCA, batch=dataset"
    return info


def _pick_root_and_start(adata) -> tuple[int, str, dict]:
    """External arrow: GSE131907 nLung author AT2. Never root on CLDN4-high."""
    nlung = (adata.obs["dataset"].astype(str) == "GSE131907") & (
        adata.obs["Sample_Origin"].astype(str) == "nLung"
    )
    subtype = adata.obs["author_subtype"].astype(str)
    at2 = nlung & subtype.isin(AUTHOR_AT2)
    leiden = adata.obs["leiden"].astype(str)
    info: dict = {
        "n_nLung": int(nlung.sum()),
        "n_nLung_author_AT2": int(at2.sum()),
        "rule": None,
        "rejected_cldn4_high_start": False,
    }
    if int(at2.sum()) < 20:
        raise SystemExit("too few GSE131907 nLung author AT2 cells to root")
    idx = np.flatnonzero(at2.to_numpy())
    scores = adata.obs.loc[at2, "score_AT2"].to_numpy()
    pick = int(idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))])
    start = str(leiden.iloc[pick])
    cl_means = (
        adata.obs.groupby(leiden, observed=True)["expr_CLDN4"].mean().sort_values(ascending=False)
    )
    top_cldn4 = str(cl_means.index[0])
    info["highest_cldn4_cluster"] = top_cldn4
    info["start_cluster_mean_CLDN4"] = float(cl_means.get(start, np.nan))
    info["top_cldn4_cluster_mean"] = float(cl_means.iloc[0])
    if start == top_cldn4:
        info["rejected_cldn4_high_start"] = True
        at2_by_cl = (
            adata.obs.loc[at2].groupby(leiden.loc[at2], observed=True)["score_AT2"].mean()
        )
        candidates = [
            c
            for c in at2_by_cl.sort_values(ascending=False).index.astype(str)
            if c != top_cldn4
        ]
        if not candidates:
            raise SystemExit("all AT2-bearing clusters are the CLDN4-high cluster")
        start = str(candidates[0])
        cand = at2 & (leiden == start)
        scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
        idx = np.flatnonzero(cand.to_numpy())
        pick = int(idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))])
        info["rule"] = (
            f"GSE131907 nLung author AT2, start cluster reassigned from "
            f"CLDN4-high Leiden {top_cldn4} to Leiden {start}"
        )
    else:
        info["rule"] = "GSE131907 nLung author AT2 (median AT2 score); start ≠ CLDN4-high cluster"
    info["start_cluster"] = start
    info["index"] = pick
    root_obs = adata.obs.iloc[pick]
    info["root_unit"] = str(root_obs.get("unit_id", ""))
    info["root_origin"] = str(root_obs.get("Sample_Origin", ""))
    info["root_subtype"] = str(root_obs.get("author_subtype", ""))
    info["root_dataset"] = str(root_obs.get("dataset", ""))
    info["root_CLDN4"] = float(root_obs.get("expr_CLDN4", np.nan))
    info["root_AT2"] = float(root_obs.get("score_AT2", np.nan))
    if str(root_obs.get("cldn4_tertile", "")) == "high":
        # Keep the start cluster; swap the root cell to an AT2 cell in that
        # cluster that is not CLDN4-high.
        cand = at2 & (leiden == start) & (adata.obs["cldn4_tertile"].astype(str) != "high")
        if int(cand.sum()) == 0:
            cand = at2 & (leiden == start)
        scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
        idx = np.flatnonzero(cand.to_numpy())
        pick = int(idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))])
        info["index"] = pick
        info["root_cell_was_cldn4_high"] = True
        root_obs = adata.obs.iloc[pick]
        info["root_CLDN4"] = float(root_obs.get("expr_CLDN4", np.nan))
        info["root_tertile"] = str(root_obs.get("cldn4_tertile", ""))
    else:
        info["root_cell_was_cldn4_high"] = False
        info["root_tertile"] = str(root_obs.get("cldn4_tertile", ""))
    return pick, start, info


def run_slingshot(adata, start_cluster: str, scratch: Path, tabdir: Path) -> dict:
    scratch.mkdir(parents=True, exist_ok=True)
    rep = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"
    rd = np.asarray(adata.obsm[rep][:, :N_PCS], dtype=float)
    cells = adata.obs_names.astype(str).to_numpy()
    rd_df = pd.DataFrame(rd, index=cells, columns=[f"PC{i+1}" for i in range(rd.shape[1])])
    rd_path = scratch / "slingshot_rd.tsv"
    cl_path = scratch / "slingshot_clusters.tsv"
    rd_df.to_csv(rd_path, sep="\t")
    pd.DataFrame({"cell_id": cells, "cluster": adata.obs["leiden"].astype(str).to_numpy()}).to_csv(
        cl_path, sep="\t", index=False
    )
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise SystemExit("Rscript missing after install; cannot write lineage table")
    env = os.environ.copy()
    env["R_LIBS_USER"] = R_LIBS_USER
    cmd = [
        rscript,
        str(HERE / "scripts" / "run_slingshot.R"),
        "--rd",
        str(rd_path),
        "--clusters",
        str(cl_path),
        "--start",
        str(start_cluster),
        "--outdir",
        str(tabdir),
    ]
    print("RUN", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, check=False, text=True, env=env)
    if proc.returncode != 0:
        raise SystemExit(f"Slingshot R failed with code {proc.returncode}")
    lin_path = tabdir / "slingshot_lineages.tsv"
    if not lin_path.is_file():
        raise SystemExit("Slingshot finished but lineage table is missing")
    pt = pd.read_csv(tabdir / "slingshot_cell_pseudotime.tsv", sep="\t")
    pt = pt.set_index("cell_id").reindex(adata.obs_names.astype(str))
    lin_cols = [c for c in pt.columns if c not in {"cluster"}]
    for c in lin_cols:
        adata.obs[f"sling_{c}"] = pd.to_numeric(pt[c], errors="coerce").to_numpy()
    if lin_cols:
        adata.obs["sling_pt_primary"] = adata.obs[f"sling_{lin_cols[0]}"]
        adata.obs["sling_pt_mean"] = np.nanmean(
            np.vstack([adata.obs[f"sling_{c}"].to_numpy() for c in lin_cols]), axis=0
        )
    else:
        raise SystemExit("Slingshot returned no lineage columns")
    info_path = tabdir / "slingshot_info.json"
    info = json.loads(info_path.read_text()) if info_path.is_file() else {}
    info["available"] = True
    info["representation"] = rep
    info["primary_column"] = lin_cols[0]
    info["lineage_table"] = str(lin_path)
    return info


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    extra = s["extra_figure"]

    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    lines = [
        "# Finding — pair GSE131907+GSE189357 epithelium, CLDN4-only REAL Slingshot/PAGA",
        "",
        "ADDITIVE. **CLDN4 only.** Pair that already differs in PR #459 "
        "(malignant CLDN4 %pos vs T/NK, n=30, ρ=−0.542; Q4 r=−0.619). "
        "That T/NK rho is **not re-audited**. Does **not** redo GSE131907-only PAGA "
        "(PR #325) or the winning-pair Slingshot/DPT (PR #449). "
        "GSE148071 is not added. No TACSTD2∩CLDN4 dual-high gate.",
        "",
        "Primary clock: **real Slingshot** (Street et al. 2018) on Harmony PCA + Leiden, "
        "start cluster = GSE131907 nLung author AT2, **never CLDN4-high**. "
        "PAGA is the graph geometry. Inferential unit = **GSE131907 sample + GSE189357 patient**. "
        "Cell-level ρ is descriptive. Barrier/keratin **excludes CLDN4**. "
        "IFN = Hallmark IFNα ∪ IFNγ mean (CLDN4 not a member).",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- Analysis cells after QC (capped ≤{s['cap_per_unit']}/unit): **n_cells = {s['n_cells']}** "
        f"(GSE131907 {s['n_cells_gse131907']}, GSE189357 {s['n_cells_gse189357']}).",
        f"- Units (GSE131907 Sample + GSE189357 patient): **n_units = {s['n_units']}** "
        f"(GSE131907 {s['n_units_gse131907']}, GSE189357 {s['n_units_gse189357']}). "
        "Do not write these as one patient n.",
        f"- Units with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} epithelial cells used for Spearman: "
        f"**n = {s['n_units_eligible']}**.",
        f"- GSE131907 nLung cells / author AT2 in the object: {s['n_cells_nLung']} / {s['n_author_AT2']}.",
        f"- GSE189357 stages (cells): {s['stage_counts']}.",
        f"- Author / marker subtypes (cells): {s['subtype_counts']}.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, "
        f"mid {s['cldn4_tertile_counts'].get('mid', 0)}, "
        f"high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Units with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms: "
        f"**n = {s['n_units_paired_tertile']}**.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        f"- IFN Hallmark union genes present: {s['n_ifn_present']} / {s['n_ifn_locked']}.",
        "- GSE148071 not used. GSE131907 PE unlabeled epithelium dropped.",
        f"- Slingshot: available={s['slingshot'].get('available')}; "
        f"version={s['slingshot'].get('version') or s['slingshot'].get('slingshot_version')}; "
        f"n_lineages={s['slingshot'].get('n_lineages')}.",
        f"- Root: {s['root']['rule']} (cell index {s['root']['index']}, unit {s['root'].get('root_unit')}, "
        f"tertile {s['root'].get('root_tertile')}, start Leiden {s['root'].get('start_cluster')}).",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- Batch: {s['harmony']['reason']}.",
        f"- Slingshot start: Leiden {s['root']['start_cluster']} "
        f"(rejected CLDN4-high start={s['root'].get('rejected_cldn4_high_start')}).",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN: Hallmark IFNα ∪ IFNγ mean (**CLDN4 out**).",
        "",
        "## Lineage table (done criterion)",
        "",
        "| lineage | start | end | n_clusters | n_cells | path | ρ CLDN4 | ρ barrier | ρ IFN |",
        "| --- | --- | --- | ---: | ---: | --- | ---: | ---: | ---: |",
    ]
    for r in s["lineages"]:
        def f(v):
            return "NA" if v is None else f"{v:.3f}"

        lines.append(
            f"| {r['lineage_id']} | {r['start_cluster']} | {r['end_cluster']} | "
            f"{r['n_clusters']} | {r['n_cells_finite_pt']} | `{r['path']}` | "
            f"{f(r.get('rho_CLDN4'))} | {f(r.get('rho_barrier'))} | {f(r.get('rho_IFN'))} |"
        )
    lines += [
        "",
        "## Primary (sample-level Spearman vs Slingshot PT, BH inside this list)",
        "",
        "| Contrast | n_units | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(row_md(r))
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | n_units | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s.get("sensitivity_spearman", []):
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {rho} | {pv} |")
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)",
        "",
        (
            f"Emitted: **{extra['emitted']}**. "
            f"Observed barrier Spearman n={extra['spearman_n']}, "
            f"ρ={extra.get('spearman_rho')}, p={extra.get('spearman_p')}."
        ),
        "",
        "| Paired contrast (high − low) | n_units | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s.get("paired_tertile", []):
        d = "NA" if r.get("delta_median") is None else f"{r['delta_median']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {d} | {pv} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- This is not a redo of PR #325 (GSE131907-only PAGA) or PR #449 (winning-pair DPT).",
        "- The pooled Slingshot Spearman mixes cohorts and nLung vs tumor and is **not** a within-tumor progression test.",
        "- GSE189357 epithelium is **marker-gated**, not author-labeled. Malignant-like is CEACAM5/6/MKI67 — **not CNV**.",
        "- PR #459 T/NK ρ is given and was not re-audited.",
        "- No TACSTD2∩CLDN4 both-high gate. GSE148071 not added.",
        "- Do not write “AT2 differentiates into LUAD because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Slingshot is an ordering, not a developmental clock.",
        "",
        "## Outputs",
        "",
        "- `results/tables/slingshot_lineages.tsv` — **done criterion** (also copied to `lineage_table.tsv`)",
        "- `results/tables/sample_level_spearman.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_along_pseudotime.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/pair_131907_189357_slingshot_cldn4/requirements.txt",
        "# R 4.3 + Bioconductor slingshot (user library ~/R/library)",
        "python3 methods/pair_131907_189357_slingshot_cldn4/scripts/run_all.py",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default=str(HERE / "results"))
    p.add_argument("--finding", default=str(HERE / "FINDING.md"))
    p.add_argument("--scratch", default="/tmp/pair_131907_189357/slingshot_scratch")
    args = p.parse_args()
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    sling = slingshot_status()
    print(json.dumps({"slingshot": sling}, indent=2), flush=True)
    if not sling.get("available"):
        raise SystemExit(f"REAL Slingshot required: {sling}")

    adata = sc.read_h5ad(args.input)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata.X = adata.layers["counts"].copy()

    adata.var["n_cells"] = np.array((adata.X > 0).sum(axis=0)).ravel()
    adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).ravel()
    adata.obs["n_umi"] = np.array(adata.X.sum(axis=1)).ravel()
    sc.pp.filter_genes(adata, min_cells=10)
    keep = (adata.obs["n_genes"] >= 200) & (adata.obs["n_umi"] >= 500)
    adata = adata[keep].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    ifn_genes = hallmark_ifn_union()
    absent: dict[str, list[str]] = {}
    for name, genes in STATES.items():
        absent[name] = _score_mean(adata, genes, f"score_{name}")
    absent["IFN_hallmark"] = _score_mean(adata, ifn_genes, "score_IFN")
    absent["IFN_core"] = _score_mean(adata, IFN_CORE, "score_IFN_core")
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        if g in adata.var_names:
            x = adata[:, g].X
            if hasattr(x, "toarray"):
                x = x.toarray()
            adata.obs[f"expr_{g}"] = np.asarray(x, dtype=float).ravel()
        else:
            adata.obs[f"expr_{g}"] = np.nan
            absent.setdefault("single_genes", []).append(g)

    cldn = adata.obs["expr_CLDN4"].to_numpy()
    q1, q2 = np.nanquantile(cldn, [1 / 3, 2 / 3])
    tert = np.full(adata.n_obs, "mid", dtype=object)
    tert[cldn <= q1] = "low"
    tert[cldn > q2] = "high"
    adata.obs["cldn4_tertile"] = pd.Categorical(tert, categories=["low", "mid", "high"], ordered=True)

    try:
        sc.pp.highly_variable_genes(adata, layer="counts", flavor="seurat_v3", n_top_genes=N_HVG)
    except ImportError:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=N_HVG)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    harmony = maybe_harmony(adata)
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2)
    except TypeError:
        sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.paga(adata, groups="leiden")
    sc.tl.diffmap(adata, n_comps=15)
    sc.tl.umap(adata)

    root_i, start_cluster, root_info = _pick_root_and_start(adata)
    adata.uns["iroot"] = root_i
    sc.tl.dpt(adata, n_dcs=10)
    adata.obs["dpt_pseudotime"] = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)

    sling_info = run_slingshot(adata, start_cluster, Path(args.scratch), tabdir)
    sling.update(sling_info)

    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    comps = _paga_components(connect, thresh=0.0)
    leiden_ids = [str(x) for x in adata.obs["leiden"].cat.categories]
    cluster_tab = []
    for cl in leiden_ids:
        sub = adata.obs[adata.obs["leiden"].astype(str) == cl]
        cluster_tab.append(
            {
                "leiden": cl,
                "n_cells": int(len(sub)),
                "n_units": int(sub["unit_id"].nunique()),
                "n_GSE131907": int((sub["dataset"] == "GSE131907").sum()),
                "n_GSE189357": int((sub["dataset"] == "GSE189357").sum()),
                "top_subtype": sub["author_subtype"].astype(str).value_counts().index[0],
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_sling_pt": float(sub["sling_pt_primary"].mean()),
                "is_start": cl == start_cluster,
            }
        )
    pd.DataFrame(cluster_tab).to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )

    lin_df = pd.read_csv(tabdir / "slingshot_lineages.tsv", sep="\t")
    lineage_rows = []
    for _, rec in lin_df.iterrows():
        col = rec["slingshot_column"]
        pt = adata.obs[f"sling_{col}"]
        mask = np.isfinite(pt.to_numpy())
        sub = adata.obs.loc[mask]
        c4 = _spearman(sub["expr_CLDN4"].to_numpy(), pt.loc[mask].to_numpy())
        br = _spearman(sub["score_barrier_keratin"].to_numpy(), pt.loc[mask].to_numpy())
        inf = _spearman(sub["score_IFN"].to_numpy(), pt.loc[mask].to_numpy())
        path = str(rec["path"]).split("->")
        start_cells = adata.obs[adata.obs["leiden"].astype(str) == path[0]]
        end_cells = adata.obs[adata.obs["leiden"].astype(str) == path[-1]]
        row = rec.to_dict()
        row.update(
            {
                "rho_CLDN4": c4["rho"],
                "p_CLDN4": c4["p"],
                "rho_barrier": br["rho"],
                "p_barrier": br["p"],
                "rho_IFN": inf["rho"],
                "p_IFN": inf["p"],
                "start_mean_AT2": float(start_cells["score_AT2"].mean()),
                "start_mean_CLDN4": float(start_cells["expr_CLDN4"].mean()),
                "end_mean_CLDN4": float(end_cells["expr_CLDN4"].mean()),
                "end_mean_barrier": float(end_cells["score_barrier_keratin"].mean()),
                "end_mean_IFN": float(end_cells["score_IFN"].mean()),
                "n_units_on_lineage": int(sub["unit_id"].nunique()),
            }
        )
        lineage_rows.append(row)
    lineage_out = pd.DataFrame(lineage_rows)
    lineage_out.to_csv(tabdir / "slingshot_lineages.tsv", sep="\t", index=False)
    lineage_out.to_csv(tabdir / "lineage_table.tsv", sep="\t", index=False)

    rows = []
    for unit, sub in adata.obs.groupby("unit_id", observed=True):
        rows.append(
            {
                "unit_id": unit,
                "dataset": str(sub["dataset"].iloc[0]),
                "Sample": str(sub["Sample"].iloc[0]),
                "Sample_Origin": str(sub["Sample_Origin"].iloc[0]),
                "patient_id": str(sub["patient_id"].iloc[0]),
                "histology": str(sub["histology"].iloc[0]),
                "stage": str(sub["stage"].iloc[0]),
                "n_cells": int(len(sub)),
                "n_author_AT2": int(sub["author_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
                "n_author_club": int(sub["author_subtype"].astype(str).isin(AUTHOR_CLUB).sum()),
                "n_author_tumor_state": int(
                    sub["author_subtype"].astype(str).isin(AUTHOR_TUMOR_STATE).sum()
                ),
                "n_cldn4_low": int((sub["cldn4_tertile"] == "low").sum()),
                "n_cldn4_mid": int((sub["cldn4_tertile"] == "mid").sum()),
                "n_cldn4_high": int((sub["cldn4_tertile"] == "high").sum()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_sling_pt": float(sub["sling_pt_primary"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_IFN_core": float(sub["score_IFN_core"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
            }
        )
    sample_df = pd.DataFrame(rows)
    if "expr_SFTPC" in adata.obs:
        sft_map = {
            unit: float(sub["expr_SFTPC"].mean())
            for unit, sub in adata.obs.groupby("unit_id", observed=True)
        }
        sample_df["mean_SFTPC"] = sample_df["unit_id"].map(sft_map)
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()

    contrasts = [
        ("CLDN4 vs Slingshot PT", "mean_CLDN4", "mean_sling_pt"),
        ("barrier/keratin (no CLDN4) vs Slingshot PT", "mean_barrier_keratin", "mean_sling_pt"),
        ("IFN (Hallmark α∪γ) vs Slingshot PT", "mean_IFN", "mean_sling_pt"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs IFN (Hallmark α∪γ)", "mean_CLDN4", "mean_IFN"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs malignant-like score", "mean_CLDN4", "mean_malignant_like"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("SFTPC vs Slingshot PT (control)", "mean_SFTPC", "mean_sling_pt"),
        ("AT2 score vs Slingshot PT (control)", "mean_AT2", "mean_sling_pt"),
    ]
    primary = []
    for name, a, b in contrasts:
        if a not in elig.columns or b not in elig.columns:
            primary.append({"contrast": name, "n": 0, "rho": None, "p": None})
            continue
        primary.append({"contrast": name, **_spearman(elig[a].to_numpy(), elig[b].to_numpy())})
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q
    pd.DataFrame(primary).to_csv(tabdir / "sample_level_spearman.tsv", sep="\t", index=False)

    elig_131 = elig[elig["dataset"] == "GSE131907"]
    elig_189 = elig[elig["dataset"] == "GSE189357"]
    elig_t = elig[elig["Sample_Origin"] == "tLung"]
    elig_n = elig[elig["Sample_Origin"] == "nLung"]
    elig_tumor = elig[~elig["Sample_Origin"].isin(["nLung"])]
    elig_ais = elig[elig["stage"].isin(["AIS", "MIA", "IAC"])]

    def sp(frame, a, b):
        if frame is None or len(frame) == 0:
            return {"n": 0, "rho": None, "p": None}
        return _spearman(frame[a].to_numpy(), frame[b].to_numpy())

    sensitivity = [
        {"contrast": "GSE131907-only CLDN4 vs Slingshot PT", **sp(elig_131, "mean_CLDN4", "mean_sling_pt")},
        {"contrast": "GSE189357-only CLDN4 vs Slingshot PT", **sp(elig_189, "mean_CLDN4", "mean_sling_pt")},
        {"contrast": "GSE131907-only barrier vs Slingshot PT", **sp(elig_131, "mean_barrier_keratin", "mean_sling_pt")},
        {"contrast": "GSE189357-only barrier vs Slingshot PT", **sp(elig_189, "mean_barrier_keratin", "mean_sling_pt")},
        {"contrast": "GSE131907-only IFN vs Slingshot PT", **sp(elig_131, "mean_IFN", "mean_sling_pt")},
        {"contrast": "GSE189357-only IFN vs Slingshot PT", **sp(elig_189, "mean_IFN", "mean_sling_pt")},
        {"contrast": "tLung-only CLDN4 vs Slingshot PT", **sp(elig_t, "mean_CLDN4", "mean_sling_pt")},
        {"contrast": "nLung-only CLDN4 vs Slingshot PT", **sp(elig_n, "mean_CLDN4", "mean_sling_pt")},
        {"contrast": "tumor-only (drop nLung) CLDN4 vs Slingshot PT", **sp(elig_tumor, "mean_CLDN4", "mean_sling_pt")},
        {"contrast": "GSE189357 AIS/MIA/IAC CLDN4 vs Slingshot PT", **sp(elig_ais, "mean_CLDN4", "mean_sling_pt")},
        {"contrast": "CLDN4 vs DPT (companion, not the clock)", **sp(elig, "mean_CLDN4", "mean_dpt")},
        {"contrast": "IFN core vs Slingshot PT", **sp(elig, "mean_IFN_core", "mean_sling_pt")},
    ]
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    paired_recs = []
    for unit, sub in adata.obs.groupby("unit_id", observed=True):
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        paired_recs.append(
            {
                "unit_id": unit,
                "dataset": str(sub["dataset"].iloc[0]),
                "Sample_Origin": str(sub["Sample_Origin"].iloc[0]),
                "n_high": int(len(hi)),
                "n_low": int(len(lo)),
                "AT2_high": float(hi["score_AT2"].mean()),
                "AT2_low": float(lo["score_AT2"].mean()),
                "barrier_high": float(hi["score_barrier_keratin"].mean()),
                "barrier_low": float(lo["score_barrier_keratin"].mean()),
                "IFN_high": float(hi["score_IFN"].mean()),
                "IFN_low": float(lo["score_IFN"].mean()),
                "sling_high": float(hi["sling_pt_primary"].mean()),
                "sling_low": float(lo["sling_pt_primary"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_recs)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    paired_rows = []
    for label, a, b in (
        ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
        ("IFN (Hallmark a+g) high vs low", "IFN_high", "IFN_low"),
        ("AT2 high vs low", "AT2_high", "AT2_low"),
        ("Slingshot PT high vs low", "sling_high", "sling_low"),
    ):
        if paired_df.empty:
            paired_rows.append({"contrast": label, "n": 0, "W": None, "p": None, "delta_median": None})
        else:
            paired_rows.append({"contrast": label, **_wilcoxon_paired(paired_df[a].to_numpy(), paired_df[b].to_numpy())})

    barrier_row = next(r for r in primary if r["contrast"].startswith("CLDN4 vs barrier"))
    emit_extra = True

    colors = {"GSE131907": "#2a6f97", "GSE189357": "#c47b2b"}

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.0))
    ax = axes[0, 0]
    sc.pl.paga(
        adata,
        color="expr_CLDN4",
        ax=ax,
        show=False,
        frameon=False,
        cmap="viridis",
        title="PAGA (Leiden) colored by mean CLDN4",
    )
    ax = axes[0, 1]
    sc.pl.umap(adata, color="expr_CLDN4", ax=ax, show=False, frameon=False, cmap="viridis", title="UMAP CLDN4")
    ax = axes[1, 0]
    sc.pl.umap(
        adata,
        color="sling_pt_primary",
        ax=ax,
        show=False,
        frameon=False,
        cmap="magma",
        title=f"UMAP Slingshot PT (start Leiden {start_cluster})",
    )
    ax = axes[1, 1]
    c4_pt = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    for ds, col in colors.items():
        sub = elig[elig["dataset"] == ds]
        ax.scatter(sub["mean_sling_pt"], sub["mean_CLDN4"], s=48, c=col, label=f"{ds} n={len(sub)}")
    ax.set_xlabel("sample-mean Slingshot PT")
    ax.set_ylabel("sample-mean CLDN4")
    ax.set_title(_fmt(c4_pt))
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"Pair GSE131907+GSE189357 epithelium scored by CLDN4   n_cells={adata.n_obs}  n_units={sample_df.shape[0]}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    # CLDN4 + barrier + IFN along pseudotime (requested extra)
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    pt = adata.obs["sling_pt_primary"].to_numpy()
    finite = np.isfinite(pt)
    bins = np.linspace(np.nanmin(pt[finite]), np.nanmax(pt[finite]), 12)
    centers = 0.5 * (bins[:-1] + bins[1:])
    panels = (
        ("expr_CLDN4", "CLDN4", "#3d5a80"),
        ("score_barrier_keratin", "barrier/keratin (no CLDN4)", "#588157"),
        ("score_IFN", "IFN Hallmark α∪γ", "#9b2226"),
    )
    for ax, (key, lab, col) in zip(axes, panels):
        y = adata.obs[key].to_numpy()
        means, lo, hi = [], [], []
        for i in range(len(bins) - 1):
            m = finite & (pt >= bins[i]) & (pt < bins[i + 1] if i < len(bins) - 2 else pt <= bins[i + 1])
            vals = y[m]
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                means.append(np.nan)
                lo.append(np.nan)
                hi.append(np.nan)
            else:
                means.append(float(np.mean(vals)))
                lo.append(float(np.quantile(vals, 0.25)))
                hi.append(float(np.quantile(vals, 0.75)))
        ax.fill_between(centers, lo, hi, color=col, alpha=0.18)
        ax.plot(centers, means, color=col, lw=2)
        ax.set_xlabel("Slingshot PT (primary lineage)")
        ax.set_ylabel(lab)
        ax.set_title(lab)
    fig.suptitle("EXTRA: CLDN4 + barrier + IFN along Slingshot pseudotime (cell bins)", fontsize=11)
    _save(fig, figdir / "fig_along_pseudotime")

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    for ax, x, y, title in (
        (axes[0], "mean_sling_pt", "mean_CLDN4", next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")),
        (axes[1], "mean_sling_pt", "mean_barrier_keratin", next(r for r in primary if r["contrast"].startswith("barrier"))),
        (axes[2], "mean_sling_pt", "mean_IFN", next(r for r in primary if r["contrast"].startswith("IFN"))),
    ):
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub[y], s=44, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel("sample-mean Slingshot PT")
        ax.set_title(_fmt(title))
        ax.legend(fontsize=7, frameon=False)
    axes[0].set_ylabel("CLDN4")
    axes[1].set_ylabel("barrier/keratin (no CLDN4)")
    axes[2].set_ylabel("IFN Hallmark α∪γ")
    fig.suptitle("EXTRA: sample-level CLDN4 / barrier / IFN vs Slingshot PT", fontsize=11)
    _save(fig, figdir / "fig_extra_sample_along_pt")

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    ct = (
        adata.obs.assign(author_subtype=adata.obs["author_subtype"].astype(str).fillna("NA"))
        .groupby(["dataset", "author_subtype"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.T.plot(kind="bar", ax=axes[0], color={"GSE131907": "#2a6f97", "GSE189357": "#c47b2b"})
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Subtypes  n_cells={adata.n_obs}")
    axes[0].legend(frameon=False, fontsize=8)
    plt.setp(axes[0].get_xticklabels(), rotation=40, ha="right")
    unit_ct = sample_df.groupby("dataset").size()
    axes[1].bar(
        unit_ct.index.astype(str),
        unit_ct.to_numpy(),
        color=["#2a6f97", "#c47b2b"][: len(unit_ct)],
    )
    axes[1].set_ylabel("units")
    axes[1].set_title(f"Honest n units={sample_df.shape[0]} (eligible {len(elig)})")
    _save(fig, figdir / "fig_honest_n")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("IFN_low", "IFN_high", "IFN Hallmark α∪γ"),
            ("sling_low", "sling_high", "Slingshot PT"),
        )
        for ax, (lo, hi, lab) in zip(axes, panels):
            for ds, col in colors.items():
                sub = paired_df[paired_df["dataset"] == ds]
                ax.scatter(sub[lo], sub[hi], s=44, c=col, label=f"{ds} n={len(sub)}")
            lims = [
                min(paired_df[lo].min(), paired_df[hi].min()),
                max(paired_df[lo].max(), paired_df[hi].max()),
            ]
            pad = 0.05 * (lims[1] - lims[0] + 1e-6)
            ax.plot([lims[0] - pad, lims[1] + pad], [lims[0] - pad, lims[1] + pad], ls="--", c="0.6", lw=1)
            ax.set_xlabel(f"CLDN4-low {lab}")
            ax.set_ylabel(f"CLDN4-high {lab}")
        axes[0].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("barrier")), keys=("W", "p")))
        axes[1].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("IFN")), keys=("W", "p")))
        axes[2].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("Slingshot")), keys=("W", "p")))
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle(f"EXTRA: within-unit CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    for color, fname, cmap in (
        ("dataset", "fig_umap_dataset", None),
        ("Sample_Origin", "fig_umap_origin", None),
        ("stage", "fig_umap_stage", None),
        ("author_subtype", "fig_umap_author_subtype", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
        ("score_IFN", "fig_umap_IFN", "viridis"),
        ("leiden", "fig_umap_leiden", None),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    sc.pl.paga(adata, ax=ax, show=False, frameon=False, title="PAGA Leiden graph")
    _save(fig, figdir / "fig_paga")

    subtype_counts = adata.obs["author_subtype"].astype(str).fillna("NA").value_counts().to_dict()
    stage_counts = adata.obs["stage"].astype(str).fillna("NA").value_counts().to_dict()
    tertile_counts = adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict()
    c4_pt = next(r for r in primary if r["contrast"] == "CLDN4 vs Slingshot PT")
    bar_pt = next(r for r in primary if r["contrast"].startswith("barrier"))
    ifn_pt = next(r for r in primary if r["contrast"].startswith("IFN"))
    c4_bar = barrier_row
    c4_ifn = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN (Hallmark α∪γ)")
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_ifn = next(r for r in paired_rows if r["contrast"].startswith("IFN"))
    s131 = next(r for r in sensitivity if r["contrast"] == "GSE131907-only CLDN4 vs Slingshot PT")
    s189 = next(r for r in sensitivity if r["contrast"] == "GSE189357-only CLDN4 vs Slingshot PT")

    parts = [
        f"Sample-level CLDN4 vs AT2-rooted Slingshot PT: {_fmt(c4_pt)}.",
        f"Barrier/keratin (no CLDN4) vs Slingshot PT: {_fmt(bar_pt)}.",
        f"IFN (Hallmark α∪γ) vs Slingshot PT: {_fmt(ifn_pt)}.",
        f"CLDN4 vs barrier/keratin: {_fmt(c4_bar)}.",
        f"CLDN4 vs IFN: {_fmt(c4_ifn)}.",
        f"CLDN4 vs AT2: {_fmt(c4_at2)}.",
        f"GSE131907-only CLDN4 vs Slingshot PT: {_fmt(s131)}.",
        f"GSE189357-only CLDN4 vs Slingshot PT: {_fmt(s189)}.",
        f"Paired CLDN4-high vs low barrier: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}.",
        f"Slingshot n_lineages={sling.get('n_lineages')} start={start_cluster}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        "The pooled Slingshot correlation mixes cohorts and nLung vs tumor and is not a within-tumor progression test.",
        "Real Slingshot was run. Not a TACSTD2 redo. No both-high gate. GSE148071 not added.",
    ]
    verdict = " ".join(parts)
    what_holds = (
        f"**What holds (n={c4_bar['n']} units).** "
        f"CLDN4 vs Slingshot PT {_fmt(c4_pt)}. "
        f"Barrier vs Slingshot PT {_fmt(bar_pt)}. "
        f"IFN vs Slingshot PT {_fmt(ifn_pt)}. "
        f"CLDN4 vs barrier (CLDN4 excluded) {_fmt(c4_bar)}. "
        f"Within-unit CLDN4-high vs low barrier {_fmt(pair_bar, keys=('W', 'p'))}; "
        f"IFN {_fmt(pair_ifn, keys=('W', 'p'))}. "
        f"**Honest split.** GSE131907-only {_fmt(s131)}; GSE189357-only {_fmt(s189)} "
        f"(GSE189357 n is the 9-patient arm — do not overclaim)."
    )

    n_ifn_present = int(sum(g in adata.var_names for g in ifn_genes))
    summary = {
        "accessions": ["GSE131907", "GSE189357"],
        "pr459_pair_differs": True,
        "gse148071_added": False,
        "not_a_pr325_redo": True,
        "not_a_pr449_redo": True,
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "dual_high": False,
        "clock": "slingshot",
        "slingshot": sling,
        "harmony": harmony,
        "cap_per_unit": 350,
        "n_cells": int(adata.n_obs),
        "n_cells_gse131907": int((adata.obs["dataset"] == "GSE131907").sum()),
        "n_cells_gse189357": int((adata.obs["dataset"] == "GSE189357").sum()),
        "n_cells_nLung": int((adata.obs["Sample_Origin"] == "nLung").sum()),
        "n_author_AT2": int(adata.obs["author_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
        "n_units": int(sample_df.shape[0]),
        "n_units_gse131907": int((sample_df["dataset"] == "GSE131907").sum()),
        "n_units_gse189357": int((sample_df["dataset"] == "GSE189357").sum()),
        "n_units_eligible": int(len(elig)),
        "n_units_paired_tertile": int(len(paired_df)),
        "subtype_counts": subtype_counts,
        "stage_counts": stage_counts,
        "cldn4_tertile_counts": tertile_counts,
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": {k: v for k, v in absent.items() if k != "IFN_hallmark"},
        "n_ifn_locked": int(len(ifn_genes)),
        "n_ifn_present": n_ifn_present,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "leiden_resolution": LEIDEN_RES,
        "lineages": lineage_rows,
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "extra_figure": {
            "emitted": bool(emit_extra),
            "spearman_rho": barrier_row["rho"],
            "spearman_p": barrier_row["p"],
            "spearman_n": barrier_row["n"],
        },
        "verdict": verdict,
        "what_holds": what_holds,
        "skipped": {
            "GSE148071": "explicitly not added",
            "GSE131907 PE": "unlabeled epithelium",
            "dual_high": "TACSTD2∩CLDN4 not used",
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accessions": ["GSE131907", "GSE189357"],
                "papers": {
                    "GSE131907": "Kim et al. Nat Commun 2020 PMID 32385277",
                    "GSE189357": "Zhu / Fan / Jiang et al. GEO GSE189357 PMID 36434043",
                    "PR459": "CLDN4 %pos vs T/NK n=30 ρ=-0.542 (given, not re-audited)",
                },
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "ifn_excludes_CLDN4": True,
                "slingshot": sling,
                "locked": {
                    "leiden_resolution": LEIDEN_RES,
                    "n_hvg": N_HVG,
                    "n_neighbors": N_NEIGHBORS,
                    "n_pcs": N_PCS,
                    "root": "GSE131907 nLung author AT2, never CLDN4-high",
                    "cap_per_unit": 350,
                },
            },
            indent=2,
        )
    )
    print(
        json.dumps(
            {
                "ok": True,
                "n_cells": adata.n_obs,
                "n_units": int(sample_df.shape[0]),
                "n_eligible": int(len(elig)),
                "n_lineages": sling.get("n_lineages"),
                "lineage_table": str(tabdir / "lineage_table.tsv"),
                "finding": args.finding,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
