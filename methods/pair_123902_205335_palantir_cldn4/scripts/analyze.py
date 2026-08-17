#!/usr/bin/env python3
"""REAL Palantir on GSE123902+GSE205335 epithelium, CLDN4-only.

ADDITIVE. Install and run the palantir package (Setty 2019). Not DPT-as-Palantir.
Root is not CLDN4-high. Donor/patient is the unit.
Track CLDN4, barrier (CLDN4 held out), and IFN along destinies.
No dual-high. No GSE148071. Pair IFN DE (−1.05) is taken as given.
Done when destiny tables exist.
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

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))
from gene_sets import COMPARATOR, CONTROLS, FOCAL, IFN, MHC_I_APM, QC_NEG, STATES  # noqa: E402

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e

LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
N_DM = 10
N_WAYPOINTS = 500
MIN_CELLS_PER_UNIT = 10
MIN_CELLS_PER_TERTILE_ARM = 8
N_TREND_BINS = 8
SEED = 1


def require_palantir() -> dict:
    try:
        import palantir
    except ImportError as e:
        raise SystemExit(
            "REAL Palantir required. pip install palantir. "
            "Do not fall back to DPT."
        ) from e
    return {
        "available": True,
        "module": "palantir",
        "version": getattr(palantir, "__version__", "unknown"),
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


def _fmt(row: dict, keys=("rho", "p")) -> str:
    n = row.get("n")
    if "W" in keys:
        if row.get("W") is None:
            return f"n={n}, W=NA, p=NA"
        return (
            f"n={n}, W={row['W']:.1f}, Δmed={row.get('delta_median'):.3f}, "
            f"p={row['p']:.3g}"
        )
    if row.get("rho") is None:
        return f"n={n}, ρ=NA, p=NA"
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


def pick_early_cell(adata) -> tuple[str, dict]:
    """Root is never CLDN4-high. Prefer GSE123902 NORMAL AT2-like, CLDN4-low."""
    c4 = adata.obs["expr_CLDN4"].to_numpy(dtype=float)
    q66 = float(np.nanquantile(c4, 2.0 / 3.0))
    not_high = adata.obs["expr_CLDN4"] <= q66
    is_norm = (adata.obs["dataset"].astype(str) == "GSE123902") & (
        adata.obs["Sample_Origin"].astype(str) == "NORMAL"
    )
    info = {
        "cldn4_q66": q66,
        "n_not_cldn4_high": int(not_high.sum()),
        "n_gse123902_normal": int(is_norm.sum()),
        "rule": None,
        "forbidden": "CLDN4-high (top tertile) never used as root",
    }
    pool = is_norm & not_high
    if int(pool.sum()) >= 10:
        c4_pool = adata.obs.loc[pool, "expr_CLDN4"]
        low_cut = float(np.nanquantile(c4_pool, 1.0 / 3.0))
        low = pool & (adata.obs["expr_CLDN4"] <= low_cut)
        if int(low.sum()) < 8:
            low = pool
        # farthest from malignant-like in PCA/harmony, among CLDN4-low NORMAL
        rep = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"
        coords = np.asarray(adata.obsm[rep][:, : min(10, adata.obsm[rep].shape[1])], dtype=float)
        mal = adata.obs["author_subtype"].astype(str).eq("Malignant cells") | (
            (adata.obs["dataset"].astype(str) == "GSE123902")
            & (adata.obs["role"].astype(str) == "tumor")
        )
        if int(mal.sum()) >= 20:
            mal_cent = coords[mal.to_numpy()].mean(axis=0)
            idx = np.flatnonzero(low.to_numpy())
            dist = np.linalg.norm(coords[idx] - mal_cent, axis=1)
            # break ties toward higher AT2
            at2 = adata.obs["score_AT2"].to_numpy(dtype=float)[idx]
            rank = dist + 0.05 * (at2 - np.nanmin(at2)) / (np.nanmax(at2) - np.nanmin(at2) + 1e-9)
            pick_i = idx[int(np.nanargmax(rank))]
            info["rule"] = (
                "GSE123902 NORMAL, CLDN4-low tertile, farthest from malignant "
                "centroid in Harmony/PCA (AT2 tie-break). Not CLDN4-high."
            )
        else:
            scores = adata.obs.loc[low, "score_AT2"]
            pick_i = int(np.flatnonzero(adata.obs_names == scores.idxmax())[0])
            info["rule"] = "GSE123902 NORMAL CLDN4-low, max AT2 score. Not CLDN4-high."
        name = str(adata.obs_names[pick_i])
        info["early_cell"] = name
        info["early_cldn4"] = float(adata.obs.loc[name, "expr_CLDN4"])
        info["early_AT2"] = float(adata.obs.loc[name, "score_AT2"])
        info["early_unit"] = str(adata.obs.loc[name, "unit_id"])
        if info["early_cldn4"] > q66:
            raise SystemExit("root picker violated CLDN4-high ban")
        return name, info

    # Fallback: AT2-high, CLDN4-low, still not CLDN4-high
    cand = not_high
    if int(cand.sum()) < 10:
        raise SystemExit("not enough non-CLDN4-high cells to pick a root")
    scores = adata.obs.loc[cand, "score_AT2"]
    name = str(scores.idxmax())
    info["rule"] = "fallback: max AT2 among cells not in CLDN4-high tertile"
    info["early_cell"] = name
    info["early_cldn4"] = float(adata.obs.loc[name, "expr_CLDN4"])
    info["early_AT2"] = float(adata.obs.loc[name, "score_AT2"])
    info["early_unit"] = str(adata.obs.loc[name, "unit_id"])
    if info["early_cldn4"] > q66:
        raise SystemExit("root picker violated CLDN4-high ban")
    return name, info


def run_real_palantir(adata, early_cell: str, pal_info: dict) -> tuple[object, dict]:
    import palantir

    run = {
        "engine": "palantir",
        "version": pal_info["version"],
        "n_components": N_DM,
        "knn": N_NEIGHBORS,
        "num_waypoints": N_WAYPOINTS,
        "early_cell": early_cell,
        "terminals_user_specified": False,
        "terminals_defined_by_CLDN4": False,
        "note": "terminals auto-detected by Palantir; not CLDN4-high by construction",
    }
    pca = np.asarray(adata.obsm.get("X_pca_harmony", adata.obsm["X_pca"]), dtype=float)
    pca_df = pd.DataFrame(pca[:, :N_PCS], index=adata.obs_names.astype(str))

    # Prefer AnnData API (palantir ≥1.3); fall back to DataFrame API.
    used = None
    try:
        adata.obsm["X_pca"] = pca[:, :N_PCS]
        palantir.utils.run_diffusion_maps(adata, n_components=N_DM)
        palantir.utils.determine_multiscale_space(adata)
        pr_res = palantir.core.run_palantir(
            adata,
            early_cell,
            knn=N_NEIGHBORS,
            num_waypoints=N_WAYPOINTS,
            n_jobs=4,
            use_early_cell_as_start=True,
        )
        used = "anndata"
    except TypeError:
        used = None
    except Exception as exc:
        print(f"AnnData Palantir API failed ({exc}); trying DataFrame API", flush=True)
        used = None

    if used is None:
        dm_res = palantir.utils.run_diffusion_maps(pca_df, n_components=N_DM)
        ms_data = palantir.utils.determine_multiscale_space(dm_res)
        if "DM_EigenVectors" in getattr(dm_res, "keys", lambda: [])():
            adata.obsm["DM_EigenVectors"] = np.asarray(dm_res["DM_EigenVectors"])
        elif isinstance(dm_res, dict) and "EigenVectors" in dm_res:
            adata.obsm["DM_EigenVectors"] = np.asarray(dm_res["EigenVectors"])
        pr_res = palantir.core.run_palantir(
            ms_data,
            early_cell,
            knn=N_NEIGHBORS,
            num_waypoints=N_WAYPOINTS,
            n_jobs=4,
            use_early_cell_as_start=True,
        )
        used = "dataframe"

    run["api"] = used
    adata.obs["palantir_pseudotime"] = pr_res.pseudotime.reindex(adata.obs_names).to_numpy()
    adata.obs["palantir_entropy"] = pr_res.entropy.reindex(adata.obs_names).to_numpy()
    branch = pr_res.branch_probs.reindex(adata.obs_names)
    branch.columns = [str(c) for c in branch.columns]
    for col in branch.columns:
        safe = "fate_" + _safe_name(col)
        adata.obs[safe] = branch[col].to_numpy(dtype=float)
    adata.obs["destiny"] = branch.idxmax(axis=1).astype(str)
    adata.obs["destiny_prob"] = branch.max(axis=1).to_numpy(dtype=float)
    run["n_terminals"] = int(branch.shape[1])
    run["terminal_ids"] = list(branch.columns)
    return pr_res, run


def _safe_name(x: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in str(x))[:80]


def label_terminals(adata, branch: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for term in branch.columns:
        name = str(term)
        if name in adata.obs_names:
            cell = name
        else:
            # waypoint / cluster label: take cells with this destiny and extreme pseudotime
            mask = adata.obs["destiny"].astype(str) == name
            if int(mask.sum()) == 0:
                cell = name
            else:
                cell = str(adata.obs.loc[mask, "palantir_pseudotime"].idxmax())
        row = {
            "terminal_id": name,
            "terminal_cell": cell if cell in adata.obs_names else "NA",
            "n_cells_maxfate": int((adata.obs["destiny"].astype(str) == name).sum()),
        }
        if cell in adata.obs_names:
            o = adata.obs.loc[cell]
            row.update(
                {
                    "dataset": str(o.get("dataset", "NA")),
                    "unit_id": str(o.get("unit_id", "NA")),
                    "tissue": str(o.get("Sample_Origin", "NA")),
                    "author_subtype": str(o.get("author_subtype", "NA")),
                    "role": str(o.get("role", "NA")),
                    "expr_CLDN4": float(o.get("expr_CLDN4", np.nan)),
                    "score_barrier": float(o.get("score_barrier", np.nan)),
                    "score_IFN": float(o.get("score_IFN", np.nan)),
                    "score_AT2": float(o.get("score_AT2", np.nan)),
                    "score_malignant_like": float(o.get("score_malignant_like", np.nan)),
                    "palantir_pseudotime": float(o.get("palantir_pseudotime", np.nan)),
                }
            )
        # post-hoc label from terminal-cell scores — not used to pick the terminal
        scores = {
            "AT2": row.get("score_AT2", np.nan),
            "malignant_like": row.get("score_malignant_like", np.nan),
            "IFN": row.get("score_IFN", np.nan),
            "barrier": row.get("score_barrier", np.nan),
        }
        finite = {k: v for k, v in scores.items() if v is not None and np.isfinite(v)}
        row["posthoc_label"] = max(finite, key=finite.get) if finite else "unlabeled"
        row["defined_by_CLDN4"] = False
        rows.append(row)
    return pd.DataFrame(rows)


def destiny_trends(adata, terminals: list[str]) -> pd.DataFrame:
    rows = []
    for term in terminals:
        mask = adata.obs["destiny"].astype(str) == str(term)
        sub = adata.obs.loc[mask, :].copy()
        if len(sub) < N_TREND_BINS:
            continue
        try:
            sub["pt_bin"] = pd.qcut(sub["palantir_pseudotime"], N_TREND_BINS, duplicates="drop")
        except ValueError:
            continue
        for b, g in sub.groupby("pt_bin", observed=True):
            rows.append(
                {
                    "destiny": str(term),
                    "pt_bin": str(b),
                    "pt_mean": float(g["palantir_pseudotime"].mean()),
                    "n_cells": int(len(g)),
                    "n_units": int(g["unit_id"].nunique()),
                    "CLDN4": float(g["expr_CLDN4"].mean()),
                    "barrier": float(g["score_barrier"].mean()),
                    "IFN": float(g["score_IFN"].mean()),
                    "AT2": float(g["score_AT2"].mean()),
                    "entropy": float(g["palantir_entropy"].mean()),
                }
            )
    return pd.DataFrame(rows)


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    pal = s["palantir"]
    root = s["root"]

    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    lines = [
        "# Finding — pair GSE123902+GSE205335, CLDN4-only REAL Palantir destinies",
        "",
        "ADDITIVE. **CLDN4 only.** No TACSTD2∩CLDN4 dual-high gate. No GSE148071.",
        "This is the pair that **differs**: malignant-cell IFN family DE is already",
        "**logFC = −1.05** (Q4 vs Q1, n=9/9, p=3.16e-4; pair IFN DE folder). That DE",
        "is taken as given and is not re-run. This folder asks whether CLDN4, a",
        "CLDN4-excluded barrier score, and IFN travel together along **real Palantir",
        "destinies** (Setty et al. 2019).",
        "",
        f"Engine: **palantir {pal.get('version')}** (`{pal.get('api')}` API). "
        "Not DPT-as-Palantir. Inferential unit = **donor (GSE123902) / patient (GSE205335)**. "
        "Cell-level ρ is descriptive. Barrier score **excludes CLDN4**. "
        f"Root rule: {root.get('rule')} "
        f"(cell `{root.get('early_cell')}`, CLDN4={root.get('early_cldn4'):.3f}, "
        f"AT2={root.get('early_AT2'):.3f}). Root is **not** CLDN4-high.",
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
        f"(GSE123902 {s['n_cells_gse123902']}, GSE205335 {s['n_cells_gse205335']}).",
        f"- Units (GSE123902 donor + GSE205335 patient; NORMAL root-pool units separate): "
        f"**n_units = {s['n_units']}** (tumor {s['n_units_tumor']}, "
        f"NORMAL-root {s['n_units_normal']}; GSE123902 {s['n_units_gse123902']}, "
        f"GSE205335 {s['n_units_gse205335']}).",
        f"- Tumor units with ≥{MIN_CELLS_PER_UNIT} epithelial cells used for Spearman: "
        f"**n = {s['n_units_eligible']}**.",
        f"- GSE123902 NORMAL epithelial cells in the object (root pool): **{s['n_cells_normal']}**.",
        f"- Root cell CLDN4 = {root.get('early_cldn4'):.3f} vs object CLDN4-high tertile cut "
        f"{root.get('cldn4_q66'):.3f} (root below the cut).",
        f"- Palantir terminals (auto): **{s['n_terminals']}** — {s['terminal_ids']}.",
        f"- Units with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms: "
        f"**n = {s['n_units_paired_tertile']}**.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        "- GSE148071 not used. Dual-high not used. GSE205335 normal-tissue epithelium dropped.",
        f"- Palantir: available=True; version={pal.get('version')}; waypoints={N_WAYPOINTS}.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}; "
        f"diffusion components {N_DM}; waypoints {N_WAYPOINTS}.",
        f"- Batch: {s['harmony'].get('reason')}.",
        f"- DPT is **not** the clock. Palantir destinies are.",
        "- Barrier genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        "- IFN score: Hallmark IFNα ∪ IFNγ (same family as the given −1.05 DE).",
        "- Terminals are Palantir-auto, then labeled post-hoc. **Not** the CLDN4-high quantile.",
        "",
        "## Primary (tumor-unit Spearman, BH inside this list)",
        "",
        "| Contrast | n_units | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for r in s["primary_spearman"]:
        lines.append(row_md(r))
    lines += [
        "",
        "## Destinies (CLDN4, barrier, IFN along each terminal)",
        "",
        "Cell assignment = argmax Palantir fate probability. Trends are bin means of",
        "unsmoothed log-normalized scores vs Palantir pseudotime on that destiny.",
        "Unit-level Spearman is the claim; cell-level is descriptive.",
        "",
    ]
    for r in s.get("destiny_unit_spearman", []):
        lines.append(f"- {r['contrast']}: {_fmt(r)}")
    lines += [
        "",
        "## Sensitivity (not in the BH family)",
        "",
        "| Contrast | n_units | ρ | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s["sensitivity_spearman"]:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {rho} | {pv} |")
    extra = s["extra_figure"]
    extra_rho = "NA" if extra.get("spearman_rho") is None else f"{round(extra['spearman_rho'], 3)}"
    extra_p = "NA" if extra.get("spearman_p") is None else f"{extra['spearman_p']:.3g}"
    lines += [
        "",
        "## Extra figure — CLDN4-high vs CLDN4-low (unit-paired)",
        "",
        f"Emitted: **{extra.get('emitted')}**. Observed tumor-unit Spearman(CLDN4, IFN) "
        f"n={extra.get('spearman_n')}, ρ={extra_rho}, p={extra_p}.",
        "",
        "| Paired contrast (high − low) | n_units | Δ median | p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in s["paired_tertile"]:
        dm = "NA" if r.get("delta_median") is None else f"{r['delta_median']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {dm} | {pv} |")
    lines += [
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- This is not a redo of the pair IFN DE (−1.05 is given).",
        "- Palantir destinies on a two-cohort merge mix protocol and tissue; they are not a clock.",
        "- Terminals were not defined as CLDN4-high. Post-hoc labels are descriptive.",
        "- Barrier score excludes CLDN4. No TACSTD2∩CLDN4 both-high gate.",
        "- GSE148071 is not in this object.",
        "- Not evidence that CLDN4 *causes* IFN or barrier change.",
        "- GSE123902 malignant/epithelium is marker-gated (no author AT2/malignant labels).",
        "",
        "## Outputs",
        "",
        "- `results/tables/destiny_cell.tsv` — **done criterion (cell destinies)**",
        "- `results/tables/destiny_unit.tsv` — **done criterion (donor/patient destinies)**",
        "- `results/tables/destiny_terminals.tsv`",
        "- `results/tables/destiny_trends.tsv` — CLDN4 / barrier / IFN along destinies",
        "- `results/tables/destiny_unit_spearman.tsv`",
        "- `results/tables/honest_n.tsv`",
        "- `results/figures/fig_destiny_trends.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/pair_123902_205335_palantir_cldn4/requirements.txt",
        "python3 methods/pair_123902_205335_palantir_cldn4/scripts/download.py \\",
        "  --outdir /tmp/pair_123902_205335_palantir",
        "python3 methods/pair_123902_205335_palantir_cldn4/scripts/extract.py \\",
        "  --data /tmp/pair_123902_205335_palantir \\",
        "  --out /tmp/pair_123902_205335_palantir/epithelium.h5ad",
        "python3 methods/pair_123902_205335_palantir_cldn4/scripts/analyze.py \\",
        "  --input /tmp/pair_123902_205335_palantir/epithelium.h5ad \\",
        "  --outdir methods/pair_123902_205335_palantir_cldn4/results \\",
        "  --finding methods/pair_123902_205335_palantir_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    ap.add_argument("--finding", type=Path, default=HERE / "FINDING.md")
    args = ap.parse_args()

    pal_info = require_palantir()
    print(json.dumps({"palantir": pal_info}, indent=2), flush=True)

    adata = sc.read_h5ad(args.input)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, min_counts=500)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    absent: dict[str, list[str]] = {}
    for key, genes in STATES.items():
        absent[key] = _score_mean(adata, genes, f"score_{key}")
    absent["IFN"] = _score_mean(adata, IFN, "score_IFN")
    absent["MHC_I_APM"] = _score_mean(adata, MHC_I_APM, "score_MHC_I")
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        if g in adata.var_names:
            X = adata[:, g].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            adata.obs[f"expr_{g}"] = np.asarray(X, dtype=float).ravel()
        else:
            adata.obs[f"expr_{g}"] = np.nan
            absent.setdefault("focal", []).append(g)

    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts")
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=N_PCS, svd_solver="arpack")
    harmony = maybe_harmony(adata)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata)

    early, root_info = pick_early_cell(adata)
    print(json.dumps({"root": root_info}, indent=2), flush=True)
    pr_res, pal_run = run_real_palantir(adata, early, pal_info)
    pal_run.update(pal_info)
    branch = pr_res.branch_probs.reindex(adata.obs_names)
    branch.columns = [str(c) for c in branch.columns]
    terminals = list(branch.columns)
    term_df = label_terminals(adata, branch)
    trends = destiny_trends(adata, terminals)

    # tertiles on CLDN4 (object-wide, readout only)
    c4 = adata.obs["expr_CLDN4"].to_numpy(dtype=float)
    q1, q2 = np.nanquantile(c4, [1.0 / 3.0, 2.0 / 3.0])
    tert = np.full(adata.n_obs, "mid", dtype=object)
    tert[c4 <= q1] = "low"
    tert[c4 > q2] = "high"
    adata.obs["cldn4_tertile"] = tert

    fate_cols = [c for c in adata.obs.columns if c.startswith("fate_")]
    cell_df = adata.obs[
        [
            "dataset",
            "unit_id",
            "patient_id",
            "Sample_Origin",
            "role",
            "author_subtype",
            "is_normal_tissue",
            "expr_CLDN4",
            "score_barrier",
            "score_IFN",
            "score_AT2",
            "score_malignant_like",
            "palantir_pseudotime",
            "palantir_entropy",
            "destiny",
            "destiny_prob",
            "cldn4_tertile",
            "leiden",
        ]
        + fate_cols
    ].copy()
    cell_df.index.name = "cell"

    agg = {
        "dataset": "first",
        "role": lambda s: "root_pool" if (s == "root_pool").any() and (s == "tumor").sum() == 0 else "tumor",
        "is_normal_tissue": lambda s: "True" if (s.astype(str) == "True").all() else "False",
        "expr_CLDN4": "mean",
        "score_barrier": "mean",
        "score_IFN": "mean",
        "score_AT2": "mean",
        "score_malignant_like": "mean",
        "palantir_pseudotime": "mean",
        "palantir_entropy": "mean",
    }
    for c in fate_cols:
        agg[c] = "mean"
    unit_df = (
        adata.obs.groupby("unit_id", observed=True)
        .agg(agg | {"patient_id": "first", "Sample_Origin": "first"})
        .rename(
            columns={
                "expr_CLDN4": "mean_CLDN4",
                "score_barrier": "mean_barrier",
                "score_IFN": "mean_IFN",
                "score_AT2": "mean_AT2",
                "score_malignant_like": "mean_malignant_like",
                "palantir_pseudotime": "mean_pseudotime",
                "palantir_entropy": "mean_entropy",
            }
        )
    )
    n_per = adata.obs.groupby("unit_id", observed=True).size().rename("n_cells")
    unit_df = unit_df.join(n_per)
    dest_mode = (
        adata.obs.groupby("unit_id", observed=True)["destiny"]
        .agg(lambda s: s.value_counts().index[0])
        .rename("destiny_mode")
    )
    unit_df = unit_df.join(dest_mode)
    unit_df["is_tumor_unit"] = unit_df["role"].astype(str).ne("root_pool")

    elig = unit_df[(unit_df["n_cells"] >= MIN_CELLS_PER_UNIT) & unit_df["is_tumor_unit"]].copy()
    all_elig = unit_df[unit_df["n_cells"] >= MIN_CELLS_PER_UNIT].copy()

    primary = []
    contrasts = [
        ("CLDN4 vs Palantir pseudotime", "mean_CLDN4", "mean_pseudotime"),
        ("CLDN4 vs barrier (no CLDN4)", "mean_CLDN4", "mean_barrier"),
        ("CLDN4 vs IFN", "mean_CLDN4", "mean_IFN"),
        ("barrier vs IFN", "mean_barrier", "mean_IFN"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("IFN vs Palantir pseudotime", "mean_IFN", "mean_pseudotime"),
        ("barrier vs Palantir pseudotime", "mean_barrier", "mean_pseudotime"),
        ("CLDN4 vs entropy", "mean_CLDN4", "mean_entropy"),
    ]
    for fate in fate_cols:
        contrasts.append((f"CLDN4 vs {fate}", "mean_CLDN4", fate))
        contrasts.append((f"IFN vs {fate}", "mean_IFN", fate))
        contrasts.append((f"barrier vs {fate}", "mean_barrier", fate))
    for lab, x, y in contrasts:
        primary.append({"contrast": lab, **_spearman(elig[x].to_numpy(), elig[y].to_numpy())})
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary])
    for r, q in zip(primary, qs):
        r["q"] = q if r["p"] is not None else None

    destiny_unit_rows = [r for r in primary if "fate_" in r["contrast"] or "IFN" in r["contrast"] or "barrier" in r["contrast"]]

    sensitivity = []
    for ds in ("GSE123902", "GSE205335"):
        sub = elig[elig["dataset"] == ds]
        for lab, x, y in (
            (f"{ds}-only CLDN4 vs pseudotime", "mean_CLDN4", "mean_pseudotime"),
            (f"{ds}-only CLDN4 vs IFN", "mean_CLDN4", "mean_IFN"),
            (f"{ds}-only CLDN4 vs barrier", "mean_CLDN4", "mean_barrier"),
            (f"{ds}-only IFN vs pseudotime", "mean_IFN", "mean_pseudotime"),
        ):
            sensitivity.append({"contrast": lab, **_spearman(sub[x].to_numpy(), sub[y].to_numpy())})
    sensitivity.append(
        {
            "contrast": "all-units incl NORMAL-root CLDN4 vs IFN",
            **_spearman(all_elig["mean_CLDN4"].to_numpy(), all_elig["mean_IFN"].to_numpy()),
        }
    )

    # paired tertile within tumor units
    paired_rows = []
    paired_rec = []
    for uid, sub in adata.obs.groupby("unit_id", observed=True):
        if str(unit_df.loc[uid, "role"]) == "root_pool":
            continue
        hi = sub[sub["cldn4_tertile"] == "high"]
        lo = sub[sub["cldn4_tertile"] == "low"]
        if len(hi) < MIN_CELLS_PER_TERTILE_ARM or len(lo) < MIN_CELLS_PER_TERTILE_ARM:
            continue
        rec = {
            "unit_id": uid,
            "dataset": str(sub["dataset"].iloc[0]),
            "barrier_high": float(hi["score_barrier"].mean()),
            "barrier_low": float(lo["score_barrier"].mean()),
            "IFN_high": float(hi["score_IFN"].mean()),
            "IFN_low": float(lo["score_IFN"].mean()),
            "pt_high": float(hi["palantir_pseudotime"].mean()),
            "pt_low": float(lo["palantir_pseudotime"].mean()),
            "AT2_high": float(hi["score_AT2"].mean()),
            "AT2_low": float(lo["score_AT2"].mean()),
        }
        paired_rec.append(rec)
    paired_df = pd.DataFrame(paired_rec)
    if not paired_df.empty:
        for lab, hi, lo in (
            ("barrier (no CLDN4) high vs low", "barrier_high", "barrier_low"),
            ("IFN high vs low", "IFN_high", "IFN_low"),
            ("Palantir PT high vs low", "pt_high", "pt_low"),
            ("AT2 high vs low", "AT2_high", "AT2_low"),
        ):
            w = _wilcoxon_paired(paired_df[hi].to_numpy(), paired_df[lo].to_numpy())
            paired_rows.append({"contrast": lab, **w})
    else:
        paired_rows = [
            {"contrast": lab, "n": 0, "W": None, "p": None, "delta_median": None}
            for lab in (
                "barrier (no CLDN4) high vs low",
                "IFN high vs low",
                "Palantir PT high vs low",
                "AT2 high vs low",
            )
        ]

    ifn_row = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN")
    emit_extra = bool(
        (ifn_row.get("rho") is not None and ifn_row["rho"] != 0 and ifn_row.get("p") is not None)
        or any((r.get("p") is not None and r["p"] < 0.05) for r in paired_rows)
        or len(paired_df) >= 4
    )

    outdir = args.outdir
    tabdir = outdir / "tables"
    figdir = outdir / "figures"
    tabdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    cell_df.to_csv(tabdir / "destiny_cell.tsv", sep="\t")
    unit_df.to_csv(tabdir / "destiny_unit.tsv", sep="\t")
    term_df.to_csv(tabdir / "destiny_terminals.tsv", sep="\t", index=False)
    trends.to_csv(tabdir / "destiny_trends.tsv", sep="\t", index=False)
    pd.DataFrame(destiny_unit_rows).to_csv(tabdir / "destiny_unit_spearman.tsv", sep="\t", index=False)
    pd.DataFrame(primary).to_csv(tabdir / "sample_level_spearman.tsv", sep="\t", index=False)
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)

    honest = pd.DataFrame(
        [
            {"item": "n_cells", "n": int(adata.n_obs)},
            {"item": "n_cells_GSE123902", "n": int((adata.obs["dataset"] == "GSE123902").sum())},
            {"item": "n_cells_GSE205335", "n": int((adata.obs["dataset"] == "GSE205335").sum())},
            {"item": "n_cells_NORMAL_root_pool", "n": int((adata.obs["role"] == "root_pool").sum())},
            {"item": "n_units", "n": int(unit_df.shape[0])},
            {"item": "n_units_tumor", "n": int(unit_df["is_tumor_unit"].sum())},
            {"item": "n_units_eligible_tumor", "n": int(len(elig))},
            {"item": "n_units_paired_tertile", "n": int(len(paired_df))},
            {"item": "n_terminals", "n": int(len(terminals))},
            {"item": "root_is_CLDN4_high", "n": 0},
        ]
    )
    honest.to_csv(tabdir / "honest_n.tsv", sep="\t", index=False)

    # figures
    colors = {"GSE123902": "#2a6f97", "GSE205335": "#b23a48"}
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))
    if not trends.empty:
        for ax, score, lab in (
            (axes[0], "CLDN4", "CLDN4"),
            (axes[1], "barrier", "barrier (no CLDN4)"),
            (axes[2], "IFN", "IFN (Hallmark α∪γ)"),
        ):
            for dest, g in trends.groupby("destiny"):
                g = g.sort_values("pt_mean")
                ax.plot(g["pt_mean"], g[score], marker="o", label=str(dest)[:18])
            ax.set_xlabel("Palantir pseudotime (bin mean)")
            ax.set_ylabel(lab)
            ax.set_title(lab)
        axes[0].legend(fontsize=6, frameon=False)
    else:
        axes[0].text(0.5, 0.5, "no destiny trends", ha="center")
    fig.suptitle("CLDN4, barrier, IFN along Palantir destinies", fontsize=11)
    _save(fig, figdir / "fig_destiny_trends")

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))
    for ax, x, y, title in (
        (axes[0], "mean_pseudotime", "mean_CLDN4", "CLDN4 vs PT"),
        (axes[1], "mean_CLDN4", "mean_barrier", "CLDN4 vs barrier"),
        (axes[2], "mean_CLDN4", "mean_IFN", "CLDN4 vs IFN"),
    ):
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub[y], s=44, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel(x.replace("mean_", ""))
        ax.set_ylabel(y.replace("mean_", ""))
        ax.set_title(title)
        ax.legend(fontsize=6, frameon=False)
    fig.suptitle(f"Tumor-unit destinies  n={len(elig)}", fontsize=11)
    _save(fig, figdir / "fig_unit_destiny")

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ct = adata.obs.groupby(["dataset", "destiny"], observed=True).size().unstack(fill_value=0)
    ct.T.plot(kind="bar", ax=axes[0], color=colors)
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Destiny membership  n_cells={adata.n_obs}")
    axes[0].legend(frameon=False, fontsize=7)
    plt.setp(axes[0].get_xticklabels(), rotation=40, ha="right")
    unit_ct = elig.groupby("dataset").size()
    axes[1].bar(unit_ct.index.astype(str), unit_ct.to_numpy(), color=["#2a6f97", "#b23a48"][: len(unit_ct)])
    axes[1].set_ylabel("tumor units")
    axes[1].set_title(f"Honest n tumor units={len(elig)} / {unit_df.shape[0]}")
    _save(fig, figdir / "fig_honest_n")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier (no CLDN4)"),
            ("IFN_low", "IFN_high", "IFN"),
            ("pt_low", "pt_high", "Palantir PT"),
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
        fig.suptitle(f"EXTRA: within-unit CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        axes[0].legend(fontsize=7, frameon=False)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    # extra: root is not CLDN4-high
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.hist(adata.obs["expr_CLDN4"], bins=40, color="0.75", edgecolor="white")
    ax.axvline(root_info["cldn4_q66"], color="#b23a48", ls="--", label="CLDN4-high tertile cut")
    ax.axvline(root_info["early_cldn4"], color="#2a6f97", lw=2, label="root cell")
    ax.set_xlabel("CLDN4 (log-normalized)")
    ax.set_ylabel("cells")
    ax.set_title("EXTRA: root is not CLDN4-high")
    ax.legend(fontsize=8, frameon=False)
    _save(fig, figdir / "fig_extra_root_not_cldn4high")

    # UMAP extras
    for color, fname, cmap in (
        ("dataset", "fig_umap_dataset", None),
        ("destiny", "fig_umap_destiny", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("expr_CLDN4", "fig_umap_cldn4", "viridis"),
        ("score_barrier", "fig_umap_barrier_no_cldn4", "viridis"),
        ("score_IFN", "fig_umap_ifn", "viridis"),
        ("palantir_pseudotime", "fig_umap_pseudotime", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    n_cells_123 = int((adata.obs["dataset"] == "GSE123902").sum())
    n_cells_205 = int((adata.obs["dataset"] == "GSE205335").sum())
    n_norm = int((adata.obs["role"].astype(str) == "root_pool").sum())
    n_u_123 = int((unit_df["dataset"] == "GSE123902").sum())
    n_u_205 = int((unit_df["dataset"] == "GSE205335").sum())
    n_u_tumor = int(unit_df["is_tumor_unit"].sum())
    n_u_norm = int((~unit_df["is_tumor_unit"]).sum())

    c4_pt = next(r for r in primary if r["contrast"] == "CLDN4 vs Palantir pseudotime")
    c4_bar = next(r for r in primary if r["contrast"] == "CLDN4 vs barrier (no CLDN4)")
    c4_ifn = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN")
    bar_ifn = next(r for r in primary if r["contrast"] == "barrier vs IFN")
    ifn_pt = next(r for r in primary if r["contrast"] == "IFN vs Palantir pseudotime")
    pair_ifn = next(r for r in paired_rows if r["contrast"].startswith("IFN"))
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))

    parts = [
        f"REAL Palantir v{pal_run.get('version')} on GSE123902+GSE205335 epithelium.",
        f"Root not CLDN4-high ({root_info.get('rule')}; CLDN4={root_info.get('early_cldn4'):.3f}).",
        f"Tumor-unit CLDN4 vs Palantir PT: {_fmt(c4_pt)}.",
        f"CLDN4 vs barrier (CLDN4 excluded): {_fmt(c4_bar)}.",
        f"CLDN4 vs IFN: {_fmt(c4_ifn)}.",
        f"barrier vs IFN: {_fmt(bar_ifn)}.",
        f"IFN vs Palantir PT: {_fmt(ifn_pt)}.",
        f"Paired CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low barrier: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"{len(terminals)} auto terminal(s): {terminals}.",
        "Pair IFN DE −1.05 is given and not re-run. No dual-high. No GSE148071.",
        "Donor/patient is the unit. Destiny tables written.",
    ]
    verdict = " ".join(parts)
    what_holds = (
        f"**What this folder tests (tumor units n={c4_ifn['n']}).** "
        f"CLDN4 vs IFN {_fmt(c4_ifn)}; CLDN4 vs barrier {_fmt(c4_bar)}; "
        f"IFN vs Palantir PT {_fmt(ifn_pt)}. "
        f"Within-unit CLDN4-high vs low IFN {_fmt(pair_ifn, keys=('W', 'p'))}. "
        "The given pair IFN DE (−1.05, 9 vs 9) is a different contrast "
        "(malignant pseudobulk Q4 vs Q1) and is not re-audited here."
    )

    summary = {
        "accessions": ["GSE123902", "GSE205335"],
        "pair_that_differs": True,
        "given_ifn_de_logfc": -1.05,
        "gse148071_added": False,
        "primary_gene": "CLDN4",
        "dual_high": False,
        "clock": "palantir",
        "palantir": pal_run,
        "harmony": harmony,
        "cap_per_unit": 250,
        "n_cells": int(adata.n_obs),
        "n_cells_gse123902": n_cells_123,
        "n_cells_gse205335": n_cells_205,
        "n_cells_normal": n_norm,
        "n_units": int(unit_df.shape[0]),
        "n_units_gse123902": n_u_123,
        "n_units_gse205335": n_u_205,
        "n_units_tumor": n_u_tumor,
        "n_units_normal": n_u_norm,
        "n_units_eligible": int(len(elig)),
        "n_units_paired_tertile": int(len(paired_df)),
        "n_terminals": int(len(terminals)),
        "terminal_ids": terminals,
        "genes_absent": {k: v for k, v in absent.items() if v},
        "root": root_info,
        "primary_spearman": primary,
        "destiny_unit_spearman": destiny_unit_rows,
        "sensitivity_spearman": sensitivity,
        "paired_tertile": paired_rows,
        "extra_figure": {
            "emitted": bool(emit_extra),
            "spearman_rho": ifn_row.get("rho"),
            "spearman_p": ifn_row.get("p"),
            "spearman_n": ifn_row.get("n"),
        },
        "verdict": verdict,
        "what_holds": what_holds,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accessions": ["GSE123902", "GSE205335"],
                "papers": {
                    "GSE123902": "Laughney et al. Nat Med 2020 PMID 32066974",
                    "GSE205335": "Ahn / Lee et al. eLife 2024 (GEO GSE205335)",
                    "palantir": "Setty et al. Nat Biotechnol 2019 PMID 30899105",
                },
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "ifn_de_given": -1.05,
                "palantir": pal_run,
                "root_not_cldn4_high": True,
                "gse148071": False,
                "dual_high": False,
            },
            indent=2,
        )
    )
    print(
        json.dumps(
            {
                "ok": True,
                "n_cells": adata.n_obs,
                "n_units": int(unit_df.shape[0]),
                "n_eligible": int(len(elig)),
                "n_terminals": len(terminals),
                "destiny_tables": [
                    str(tabdir / "destiny_cell.tsv"),
                    str(tabdir / "destiny_unit.tsv"),
                    str(tabdir / "destiny_terminals.tsv"),
                    str(tabdir / "destiny_trends.tsv"),
                ],
                "finding": str(args.finding),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
