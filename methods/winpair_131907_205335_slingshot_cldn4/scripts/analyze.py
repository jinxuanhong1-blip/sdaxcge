#!/usr/bin/env python3
"""Slingshot/Palantir-style trajectory on winning-pair epithelium, CLDN4 only.

ADDITIVE. GSE131907 + GSE205335 epithelium. GSE207422 is not added.
Slingshot R is not required: if R/slingshot is missing, AT2-rooted diffusion
pseudotime (scanpy DPT) is the documented primary clock. Palantir is optional.
No TACSTD2∩CLDN4 dual-high gate. Inferential unit = sample/patient.
Barrier/keratin excludes CLDN4.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import (  # noqa: E402
    AUTHOR_AT2,
    AUTHOR_CLUB,
    AUTHOR_TUMOR_STATE,
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
MIN_CELLS_PER_SAMPLE_FOR_MEAN = 10
MIN_CELLS_PER_TERTILE_ARM = 8
EXTRA_BARRIER_RHO_GT = 0.0
EXTRA_BARRIER_P_LT = 0.05


def slingshot_status() -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return {
            "available": False,
            "reason": "Rscript not on PATH",
            "fallback": "scanpy diffusion pseudotime (Haghverdi et al. 2016)",
        }
    try:
        proc = subprocess.run(
            [rscript, "-e", 'packageVersion("slingshot")'],
            check=False,
            capture_output=True,
            text=True,
            timeout=20,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {
            "available": False,
            "reason": f"Rscript slingshot probe failed: {exc}",
            "fallback": "scanpy diffusion pseudotime (Haghverdi et al. 2016)",
        }
    if proc.returncode != 0:
        return {
            "available": False,
            "reason": (proc.stderr or proc.stdout or "slingshot not installed").strip()[:300],
            "fallback": "scanpy diffusion pseudotime (Haghverdi et al. 2016)",
        }
    return {"available": True, "version": (proc.stdout or "").strip()}


def palantir_status() -> dict:
    try:
        import palantir  # noqa: F401

        return {"available": True, "module": "palantir"}
    except ImportError:
        return {
            "available": False,
            "reason": "python package palantir not installed",
            "note": "optional companion; DPT remains the primary clock",
        }


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


def _pick_root(adata) -> tuple[int, dict]:
    """External arrow: GSE131907 nLung author AT2. Never root on CLDN4-high."""
    nlung = (adata.obs["dataset"].astype(str) == "GSE131907") & (
        adata.obs["Sample_Origin"].astype(str) == "nLung"
    )
    subtype = adata.obs["author_subtype"].astype(str)
    at2 = nlung & subtype.isin(AUTHOR_AT2)
    info = {
        "n_nLung": int(nlung.sum()),
        "n_nLung_author_AT2": int(at2.sum()),
        "rule": None,
    }
    if int(at2.sum()) >= 20:
        idx = np.flatnonzero(at2.to_numpy())
        scores = adata.obs.loc[at2, "score_AT2"].to_numpy()
        pick = idx[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))]
        info["rule"] = "GSE131907 nLung author AT2 (median AT2 score)"
        return int(pick), info
    if "leiden" not in adata.obs:
        raise SystemExit("leiden missing before root pick")
    sub = adata.obs.loc[nlung, ["leiden", "score_AT2"]]
    if sub.empty:
        raise SystemExit("no GSE131907 nLung cells to root DPT")
    means = sub.groupby("leiden", observed=True)["score_AT2"].mean().sort_values(ascending=False)
    top = str(means.index[0])
    cand = nlung & (adata.obs["leiden"].astype(str) == top)
    scores = adata.obs.loc[cand, "score_AT2"].to_numpy()
    idx = np.flatnonzero(cand.to_numpy())
    pick = idx[int(np.nanargmax(scores))]
    info["rule"] = f"GSE131907 nLung Leiden {top} max AT2 score (no author AT2)"
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
    adata.obsm["X_pca_harmony"] = np.asarray(ho.Z_corr).T
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_pca_harmony")
    info["used"] = True
    info["reason"] = "harmonypy on PCA, batch=dataset"
    return info


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    extra = s["extra_figure"]
    extra_rho = "NA" if extra.get("spearman_rho") is None else f"{round(extra['spearman_rho'], 3)}"
    extra_p = "NA" if extra.get("spearman_p") is None else f"{extra['spearman_p']:.3g}"

    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    sling = s["slingshot"]
    clock = (
        "Slingshot"
        if sling.get("available")
        else "documented AT2-rooted diffusion pseudotime (scanpy DPT; Slingshot R missing)"
    )
    lines = [
        "# Finding — winning-pair GSE131907+GSE205335 epithelium, CLDN4-only trajectory",
        "",
        "ADDITIVE. **CLDN4 only.** Winning pair from the CLDN4-first combinatorial search "
        "(PR #290: author %pos GSE131907+GSE205335, n=43, ρ=−0.479 vs T/NK). "
        "This folder does **not** redo GSE131907-only PAGA (PR #325). "
        "GSE207422 is not added. No TACSTD2∩CLDN4 dual-high gate.",
        "",
        f"Primary clock: **{clock}**. Inferential unit = **sample/patient**. "
        "Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. "
        "Root is GSE131907 nLung author AT2, never CLDN4-high.",
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- Analysis cells after QC (capped ≤{s['cap_per_unit']}/unit): **n_cells = {s['n_cells']}** "
        f"(GSE131907 {s['n_cells_gse131907']}, GSE205335 {s['n_cells_gse205335']}).",
        f"- Units (GSE131907 Sample + GSE205335 patient): **n_units = {s['n_units']}** "
        f"(GSE131907 {s['n_units_gse131907']}, GSE205335 {s['n_units_gse205335']}).",
        f"- Units with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} epithelial cells used for Spearman: "
        f"**n = {s['n_units_eligible']}**.",
        f"- GSE131907 nLung cells / author AT2 in the object: {s['n_cells_nLung']} / {s['n_author_AT2']}.",
        f"- Author subtypes (cells): {s['subtype_counts']}.",
        f"- CLDN4 tertile cells: low {s['cldn4_tertile_counts'].get('low', 0)}, "
        f"mid {s['cldn4_tertile_counts'].get('mid', 0)}, "
        f"high {s['cldn4_tertile_counts'].get('high', 0)}.",
        f"- Units with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms: "
        f"**n = {s['n_units_paired_tertile']}**.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        "- GSE207422 not used. GSE131907 PE unlabeled epithelium dropped. "
        "GSE205335 normal-tissue samples dropped.",
        f"- Slingshot: available={sling.get('available')}; {sling.get('reason') or sling.get('version')}.",
        f"- Palantir: available={s['palantir'].get('available')}; "
        f"{s['palantir'].get('reason') or s['palantir'].get('module')}.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- Batch: {s['harmony']['reason']}.",
        f"- DPT root: {s['root']['rule']} (root cell index {s['root']['index']}, unit {s['root'].get('root_unit')}).",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "",
        "## Primary (sample-level Spearman, BH inside this list)",
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
            f"Rule: sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) ρ>0 and p<{EXTRA_BARRIER_P_LT}, "
            f"or any paired tertile Wilcoxon p<{EXTRA_BARRIER_P_LT}, or n_paired≥4. "
            f"Observed Spearman n={extra['spearman_n']}, ρ={extra_rho}, p={extra_p}."
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
        "- This is not a redo of PR #325 (GSE131907-only PAGA).",
        "- The pooled DPT Spearman mixes cohorts, nLung and tumor, and is **not** a within-tumor progression test.",
        "- Malignant-like is author tS1/tS2/tS3 / Malignant cells and/or CEACAM5/6/MKI67 — **not CNV**.",
        "- GSE205335 is an ICI biopsy/effusion cohort, but this analysis is **not** an ICI / MPR / RECIST test.",
        "- No TACSTD2∩CLDN4 both-high gate.",
        "- Do not write “AT2 differentiates into LUAD because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Slingshot R was not run. DPT is an ordering, not a clock.",
        "",
        "## Outputs",
        "",
        "- `results/tables/sample_level_spearman.tsv` — **done criterion**",
        "- `results/tables/sample_means.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/winpair_131907_205335_slingshot_cldn4/requirements.txt",
        "python3 methods/winpair_131907_205335_slingshot_cldn4/scripts/download.py \\",
        "  --out /tmp/winpair_131907_205335",
        "python3 methods/winpair_131907_205335_slingshot_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/winpair_131907_205335 \\",
        "  --out /tmp/winpair_131907_205335/epithelium.h5ad",
        "python3 methods/winpair_131907_205335_slingshot_cldn4/scripts/analyze.py \\",
        "  --input /tmp/winpair_131907_205335/epithelium.h5ad \\",
        "  --outdir methods/winpair_131907_205335_slingshot_cldn4/results \\",
        "  --finding methods/winpair_131907_205335_slingshot_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default="methods/winpair_131907_205335_slingshot_cldn4/results")
    p.add_argument("--finding", default="methods/winpair_131907_205335_slingshot_cldn4/FINDING.md")
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
    pal = palantir_status()
    print(json.dumps({"slingshot": sling, "palantir": pal}, indent=2), flush=True)

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

    root_i, root_info = _pick_root(adata)
    adata.uns["iroot"] = root_i
    sc.tl.dpt(adata, n_dcs=10)
    root_info["index"] = root_i
    root_info["root_unit"] = str(adata.obs.iloc[root_i].get("unit_id", ""))
    root_info["root_origin"] = str(adata.obs.iloc[root_i].get("Sample_Origin", ""))
    root_info["root_subtype"] = str(adata.obs.iloc[root_i].get("author_subtype", ""))
    root_info["root_dataset"] = str(adata.obs.iloc[root_i].get("dataset", ""))

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
                "n_GSE205335": int((sub["dataset"] == "GSE205335").sum()),
                "top_subtype": sub["author_subtype"].astype(str).value_counts().index[0],
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan).mean()),
            }
        )
    pd.DataFrame(cluster_tab).to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )

    dpt = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)
    adata.obs["dpt_pseudotime"] = dpt
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
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
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
        if a not in elig.columns or b not in elig.columns:
            primary.append({"contrast": name, "n": 0, "rho": None, "p": None})
            continue
        primary.append({"contrast": name, **_spearman(elig[a].to_numpy(), elig[b].to_numpy())})
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = None if r["p"] is None else q
    pd.DataFrame(primary).to_csv(tabdir / "sample_level_spearman.tsv", sep="\t", index=False)

    elig_131 = elig[elig["dataset"] == "GSE131907"]
    elig_205 = elig[elig["dataset"] == "GSE205335"]
    elig_t = elig[elig["Sample_Origin"] == "tLung"]
    elig_n = elig[elig["Sample_Origin"] == "nLung"]
    elig_tumor = elig[~elig["Sample_Origin"].isin(["nLung"])]
    elig_adcsq = elig_205[elig_205["histology"].isin(["ADC", "SQ"])]

    def sp(frame, a, b):
        if frame is None or len(frame) == 0:
            return {"n": 0, "rho": None, "p": None}
        return _spearman(frame[a].to_numpy(), frame[b].to_numpy())

    sensitivity = [
        {"contrast": "GSE131907-only CLDN4 vs DPT", **sp(elig_131, "mean_CLDN4", "mean_dpt")},
        {"contrast": "GSE205335-only CLDN4 vs DPT", **sp(elig_205, "mean_CLDN4", "mean_dpt")},
        {"contrast": "GSE131907-only CLDN4 vs AT2", **sp(elig_131, "mean_CLDN4", "mean_AT2")},
        {"contrast": "GSE205335-only CLDN4 vs AT2", **sp(elig_205, "mean_CLDN4", "mean_AT2")},
        {"contrast": "GSE131907-only CLDN4 vs barrier/keratin (no CLDN4)", **sp(elig_131, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "GSE205335-only CLDN4 vs barrier/keratin (no CLDN4)", **sp(elig_205, "mean_CLDN4", "mean_barrier_keratin")},
        {"contrast": "tLung-only CLDN4 vs DPT", **sp(elig_t, "mean_CLDN4", "mean_dpt")},
        {"contrast": "nLung-only CLDN4 vs DPT", **sp(elig_n, "mean_CLDN4", "mean_dpt")},
        {"contrast": "tumor-only (drop nLung) CLDN4 vs DPT", **sp(elig_tumor, "mean_CLDN4", "mean_dpt")},
        {"contrast": "GSE205335 ADC+SQ CLDN4 vs DPT", **sp(elig_adcsq, "mean_CLDN4", "mean_dpt")},
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
                "malignant_high": float(hi["score_malignant_like"].mean()),
                "malignant_low": float(lo["score_malignant_like"].mean()),
                "dpt_high": float(hi["dpt_pseudotime"].mean()),
                "dpt_low": float(lo["dpt_pseudotime"].mean()),
                "club_high": float(hi["score_club"].mean()),
                "club_low": float(lo["score_club"].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_recs)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    paired_rows = []
    for label, a, b in (
        ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
        ("AT2 high vs low", "AT2_high", "AT2_low"),
        ("malignant-like high vs low", "malignant_high", "malignant_low"),
        ("DPT high vs low", "dpt_high", "dpt_low"),
        ("club high vs low", "club_high", "club_low"),
    ):
        if paired_df.empty:
            paired_rows.append({"contrast": label, "n": 0, "W": None, "p": None, "delta_median": None})
        else:
            paired_rows.append({"contrast": label, **_wilcoxon_paired(paired_df[a].to_numpy(), paired_df[b].to_numpy())})

    barrier_row = next(r for r in primary if r["contrast"].startswith("CLDN4 vs barrier"))
    paired_sig = any(r.get("p") is not None and r["p"] < EXTRA_BARRIER_P_LT for r in paired_rows)
    emit_extra = (
        (
            barrier_row["rho"] is not None
            and barrier_row["p"] is not None
            and barrier_row["rho"] > EXTRA_BARRIER_RHO_GT
            and barrier_row["p"] < EXTRA_BARRIER_P_LT
        )
        or paired_sig
        or len(paired_df) >= 4
    )

    # ---- figures ----
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
    sc.pl.umap(adata, color="dpt_pseudotime", ax=ax, show=False, frameon=False, cmap="magma", title="UMAP DPT (AT2-rooted)")
    ax = axes[1, 1]
    c4_dpt = next(r for r in primary if r["contrast"] == "CLDN4 vs DPT")
    colors = {"GSE131907": "#2a6f97", "GSE205335": "#b23a48"}
    for ds, col in colors.items():
        sub = elig[elig["dataset"] == ds]
        ax.scatter(sub["mean_dpt"], sub["mean_CLDN4"], s=48, c=col, label=f"{ds} n={len(sub)}")
    ax.set_xlabel("sample-mean DPT")
    ax.set_ylabel("sample-mean CLDN4")
    rho_s = "NA" if c4_dpt["rho"] is None else f"{c4_dpt['rho']:.2f}"
    p_s = "NA" if c4_dpt["p"] is None else f"{c4_dpt['p']:.3g}"
    ax.set_title(f"CLDN4 vs DPT  n={c4_dpt['n']}  ρ={rho_s}  p={p_s}")
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"Winning pair epithelium scored by CLDN4   n_cells={adata.n_obs}  n_units={sample_df.shape[0]}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    ct = (
        adata.obs.assign(author_subtype=adata.obs["author_subtype"].astype(str).fillna("NA"))
        .groupby(["dataset", "author_subtype"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.T.plot(kind="bar", ax=axes[0], color={"GSE131907": "#2a6f97", "GSE205335": "#b23a48"})
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Author subtypes  n_cells={adata.n_obs}")
    axes[0].legend(frameon=False, fontsize=8)
    plt.setp(axes[0].get_xticklabels(), rotation=40, ha="right")
    unit_ct = sample_df.groupby("dataset").size()
    axes[1].bar(unit_ct.index.astype(str), unit_ct.to_numpy(), color=["#2a6f97", "#b23a48"][: len(unit_ct)])
    axes[1].set_ylabel("units")
    axes[1].set_title(f"Honest n units={sample_df.shape[0]} (eligible {len(elig)})")
    _save(fig, figdir / "fig_honest_n")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
            ("AT2_low", "AT2_high", "AT2 score"),
            ("dpt_low", "dpt_high", "DPT"),
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
        axes[1].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("AT2")), keys=("W", "p")))
        axes[2].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("DPT")), keys=("W", "p")))
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle(f"EXTRA: within-unit CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    for color, fname, cmap in (
        ("dataset", "fig_umap_dataset", None),
        ("Sample_Origin", "fig_umap_origin", None),
        ("author_subtype", "fig_umap_author_subtype", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    # extra: sample-level CLDN4 vs AT2 and vs barrier
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    for ax, x, y, row, xlab, ylab in (
        (axes[0], "mean_AT2", "mean_CLDN4", c4_at2, "sample-mean AT2", "sample-mean CLDN4"),
        (axes[1], "mean_barrier_keratin", "mean_CLDN4", barrier_row, "sample-mean barrier/keratin (no CLDN4)", "sample-mean CLDN4"),
    ):
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub[y], s=48, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: sample-level CLDN4 vs AT2 / barrier (CLDN4 excluded)", fontsize=11)
    _save(fig, figdir / "fig_extra_sample_cldn4_programs")

    subtype_counts = adata.obs["author_subtype"].astype(str).fillna("NA").value_counts().to_dict()
    tertile_counts = adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict()
    c4_dpt = next(r for r in primary if r["contrast"] == "CLDN4 vs DPT")
    c4_at2 = next(r for r in primary if r["contrast"] == "CLDN4 vs AT2 score")
    c4_bar = barrier_row
    c4_mal = next(r for r in primary if r["contrast"] == "CLDN4 vs malignant-like score")
    c4_t2 = next(r for r in primary if r["contrast"].startswith("CLDN4 vs TACSTD2"))
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_at2 = next(r for r in paired_rows if r["contrast"].startswith("AT2"))
    s131 = next(r for r in sensitivity if r["contrast"] == "GSE131907-only CLDN4 vs DPT")
    s205 = next(r for r in sensitivity if r["contrast"] == "GSE205335-only CLDN4 vs DPT")

    parts = [
        f"Sample-level CLDN4 vs AT2-rooted DPT: {_fmt(c4_dpt)}.",
        f"CLDN4 vs AT2 score: {_fmt(c4_at2)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded from the score): {_fmt(c4_bar)}.",
        f"CLDN4 vs malignant-like: {_fmt(c4_mal)}.",
        f"CLDN4 vs TACSTD2 (comparator only): {_fmt(c4_t2)}.",
        f"GSE131907-only CLDN4 vs DPT: {_fmt(s131)}.",
        f"GSE205335-only CLDN4 vs DPT: {_fmt(s205)}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low AT2: {_fmt(pair_at2, keys=('W', 'p'))}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        "The pooled DPT correlation mixes cohorts and nLung vs tumor and is not a within-tumor progression test.",
        "Slingshot R was not run; DPT is the documented clock. Not a TACSTD2 redo. No both-high gate. GSE207422 not added.",
    ]
    verdict = " ".join(parts)

    summary = {
        "accessions": ["GSE131907", "GSE205335"],
        "winning_pair": True,
        "gse207422_added": False,
        "not_a_pr325_redo": True,
        "primary_gene": "CLDN4",
        "not_a_tacstd2_redo": True,
        "dual_high": False,
        "clock": "scanpy_dpt" if not sling.get("available") else "slingshot",
        "slingshot": sling,
        "palantir": pal,
        "harmony": harmony,
        "cap_per_unit": 350,
        "n_cells": int(adata.n_obs),
        "n_cells_gse131907": int((adata.obs["dataset"] == "GSE131907").sum()),
        "n_cells_gse205335": int((adata.obs["dataset"] == "GSE205335").sum()),
        "n_cells_nLung": int((adata.obs["Sample_Origin"] == "nLung").sum()),
        "n_author_AT2": int(adata.obs["author_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
        "n_units": int(sample_df.shape[0]),
        "n_units_gse131907": int((sample_df["dataset"] == "GSE131907").sum()),
        "n_units_gse205335": int((sample_df["dataset"] == "GSE205335").sum()),
        "n_units_eligible": int(len(elig)),
        "n_units_paired_tertile": int(len(paired_df)),
        "subtype_counts": subtype_counts,
        "cldn4_tertile_counts": tertile_counts,
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(len(comps)),
        "n_leiden": int(len(leiden_ids)),
        "leiden_resolution": LEIDEN_RES,
        "primary_spearman": primary,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "extra_figure": {
            "emitted": bool(emit_extra),
            "spearman_rho": barrier_row["rho"],
            "spearman_p": barrier_row["p"],
            "spearman_n": barrier_row["n"],
        },
        "qc": {
            "min_genes": 200,
            "min_umi": 500,
            "lineage_EPCAM_mean": float(adata.obs["expr_EPCAM"].mean()) if "expr_EPCAM" in adata.obs else None,
            "lineage_PTPRC_mean": float(adata.obs["expr_PTPRC"].mean()) if "expr_PTPRC" in adata.obs else None,
        },
        "verdict": verdict,
        "skipped": {
            "GSE207422": "explicitly not added",
            "GSE131907 PE": "unlabeled epithelium",
            "GSE205335 normal tissues": "Normal Lung / LN / Brain dropped",
            "slingshot": sling.get("reason"),
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accessions": ["GSE131907", "GSE205335"],
                "papers": {
                    "GSE131907": "Kim et al. Nat Commun 2020 PMID 32385277",
                    "GSE205335": "Ahn / Lee et al. eLife 2024 (GEO GSE205335)",
                },
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "slingshot": sling,
                "palantir": pal,
                "locked": {
                    "leiden_resolution": LEIDEN_RES,
                    "n_hvg": N_HVG,
                    "n_neighbors": N_NEIGHBORS,
                    "n_pcs": N_PCS,
                    "root": "GSE131907 nLung author AT2",
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
                "extra": emit_extra,
                "finding": args.finding,
                "sample_table": str(tabdir / "sample_level_spearman.tsv"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
