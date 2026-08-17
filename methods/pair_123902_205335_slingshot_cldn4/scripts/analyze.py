#!/usr/bin/env python3
"""Real Palantir/PAGA trajectory on GSE123902+GSE205335 epithelium, CLDN4 only.

ADDITIVE. Pair that differs in PR #459 (T/NK not re-audited).
PR #473 malignant IFN DE (IFN −1.05) is a different contrast and is not re-run.
Root is GSE123902 NORMAL AT2-like, never CLDN4-high.
Inferential unit = GSE123902 donor + GSE205335 patient (tumor units).
Barrier/keratin excludes CLDN4. No dual-high. No GSE148071.
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
    COMPARATOR,
    CONTROLS,
    FOCAL,
    QC_NEG,
    STATES,
    ifn_genes,
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
GIVEN = {
    "combo": "GSE123902+GSE205335",
    "score": "%pos",
    "n": 35,
    "rho": -0.522,
    "rho_p": 0.002,
    "q4_r": -0.802,
    "q4_p": 0.005,
    "n_q1": 9,
    "n_q4": 9,
    "source": "PR459_given_not_re_audited",
}
IFN_DE = {
    "source": "PR473_given",
    "family": "IFN",
    "n": 18,
    "n_q1": 9,
    "n_q4": 9,
    "logFC": -1.046,
    "p": 0.0003161,
}


def slingshot_status() -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return {
            "available": False,
            "reason": "Rscript not on PATH",
            "note": "Palantir is the installed real clock",
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
        return {"available": False, "reason": f"Rscript slingshot probe failed: {exc}"}
    if proc.returncode != 0:
        return {
            "available": False,
            "reason": (proc.stderr or proc.stdout or "slingshot not installed").strip()[:300],
        }
    return {"available": True, "version": (proc.stdout or "").strip()}


def palantir_status() -> dict:
    try:
        import palantir

        return {"available": True, "module": "palantir", "version": getattr(palantir, "__version__", "unknown")}
    except ImportError:
        return {"available": False, "reason": "python package palantir not installed"}


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
    """External arrow: GSE123902 NORMAL AT2-like. Never root on CLDN4-high."""
    normal = (adata.obs["dataset"].astype(str) == "GSE123902") & (
        adata.obs["Sample_Origin"].astype(str) == "NORMAL"
    )
    tert = adata.obs["cldn4_tertile"].astype(str)
    not_high = tert != "high"
    at2 = adata.obs["score_AT2"].to_numpy()
    info = {
        "n_normal": int(normal.sum()),
        "n_normal_not_cldn4_high": int((normal & not_high).sum()),
        "rule": None,
    }
    cand = normal & not_high
    if int(cand.sum()) >= 10:
        scores = at2[cand.to_numpy()]
        # AT2-like: top AT2 tercile among NORMAL, not CLDN4-high
        cut = np.nanquantile(scores, 2 / 3)
        at2_like = cand & (adata.obs["score_AT2"] >= cut)
        if int(at2_like.sum()) < 5:
            at2_like = cand
        idx = np.flatnonzero(at2_like.to_numpy())
        pick_scores = at2[idx]
        pick = idx[int(np.nanargmin(np.abs(pick_scores - np.nanmedian(pick_scores))))]
        # refuse if the pick is somehow CLDN4-high
        if str(adata.obs.iloc[pick]["cldn4_tertile"]) == "high":
            raise SystemExit("root pick landed on CLDN4-high; refusing")
        info["rule"] = "GSE123902 NORMAL AT2-like (median AT2; CLDN4-high excluded)"
        info["n_candidate"] = int(at2_like.sum())
        return int(pick), info
    # fallback: lowest CLDN4 in top AT2 tercile of the whole object, still not CLDN4-high
    cut = np.nanquantile(at2, 2 / 3)
    fb = (adata.obs["score_AT2"] >= cut) & not_high
    if int(fb.sum()) < 5:
        fb = not_high
    idx = np.flatnonzero(fb.to_numpy())
    cldn = adata.obs["expr_CLDN4"].to_numpy()[idx]
    pick = idx[int(np.nanargmin(cldn))]
    if str(adata.obs.iloc[pick]["cldn4_tertile"]) == "high":
        raise SystemExit("fallback root pick landed on CLDN4-high; refusing")
    info["rule"] = "fallback: min CLDN4 in top AT2 tercile (not CLDN4-high)"
    info["n_candidate"] = int(fb.sum())
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


def run_palantir(adata, root_i: int, pca_key: str) -> dict:
    import palantir

    root_name = str(adata.obs_names[root_i])
    palantir.utils.run_diffusion_maps(adata, n_components=10, knn=N_NEIGHBORS, pca_key=pca_key)
    palantir.utils.determine_multiscale_space(adata)
    palantir.core.run_palantir(
        adata,
        early_cell=root_name,
        knn=N_NEIGHBORS,
        num_waypoints=500,
        seed=20,
    )
    if "palantir_pseudotime" not in adata.obs:
        raise SystemExit("Palantir ran but palantir_pseudotime missing")
    pt = adata.obs["palantir_pseudotime"].replace([np.inf, -np.inf], np.nan)
    return {
        "available": True,
        "root_cell": root_name,
        "n_finite": int(np.isfinite(pt).sum()),
        "min": float(np.nanmin(pt)),
        "max": float(np.nanmax(pt)),
        "pca_key": pca_key,
    }


def write_finding(path: Path, ctx: dict) -> None:
    s = ctx["summary"]
    extra = s["extra_figure"]

    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    clock = s["clock_label"]
    lines = [
        "# Finding — pair GSE123902+GSE205335 epithelium, CLDN4-only Palantir/PAGA",
        "",
        "ADDITIVE. **CLDN4 only.** Pair that **differs** in PR #459 "
        "(GSE123902+GSE205335 %pos vs T/NK, n=35, ρ=−0.522, Q4 vs Q1 r=−0.802). "
        "That T/NK Spearman is **taken as given** and is not re-audited. "
        "PR #473 malignant IFN DE (IFN logFC −1.05, n=9/9) is a **different** contrast "
        "and is not re-run. No TACSTD2∩CLDN4 dual-high gate. GSE148071 is not added.",
        "",
        f"Primary clock: **{clock}**. Inferential unit = **GSE123902 donor + GSE205335 patient** "
        "(tumor units only). Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. "
        "IFN score = Hallmark IFNα ∪ IFNγ (same family as PR #473). "
        "Root is GSE123902 NORMAL AT2-like, **never CLDN4-high**.",
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
        f"- Tumor units (GSE123902 donor + GSE205335 patient): **n_tumor_units = {s['n_tumor_units']}** "
        f"(GSE123902 {s['n_tumor_gse123902']}, GSE205335 {s['n_tumor_gse205335']}).",
        f"- Tumor units with ≥{MIN_CELLS_PER_SAMPLE_FOR_MEAN} epithelial cells used for Spearman: **n = {s['n_units_eligible']}**.",
        f"- GSE123902 NORMAL cells in the object (root material, not primary units): **{s['n_cells_normal']}** "
        f"({s['n_normal_units']} libraries).",
        f"- Author / marker subtypes (cells): {s['subtype_counts']}.",
        f"- CLDN4 tertile cells: {s['cldn4_tertile_counts']}.",
        f"- Units with ≥{MIN_CELLS_PER_TERTILE_ARM} cells in both CLDN4-high and CLDN4-low arms: **n = {s['n_units_paired_tertile']}**.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        "- Given PR #459 T/NK n=35 is **not** this n.",
        "- Given PR #473 IFN DE n=18 (9/9) is **not** this n.",
        "- GSE148071 not used. Dual-high not used. T/NK infiltrate not re-scored.",
        f"- Slingshot: available={s['slingshot'].get('available')}; {s['slingshot'].get('reason', s['slingshot'].get('version', ''))}.",
        f"- Palantir: available={s['palantir'].get('available')}; version={s['palantir'].get('version', s['palantir'].get('reason', ''))}.",
        "",
        "## Locked choices",
        "",
        f"- Leiden resolution {LEIDEN_RES}; HVG {N_HVG}; neighbors {N_NEIGHBORS}; PCs {N_PCS}.",
        f"- Batch: {s['harmony'].get('reason')}.",
        f"- Root: {s['root'].get('rule')} (root cell index {s['root'].get('index')}, unit {s['root'].get('root_unit')}).",
        f"- PAGA components at connectivity>0: **{s['n_paga_components']}** among {s['n_leiden']} Leiden vertices.",
        "- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).",
        f"- IFN genes: Hallmark IFNα ∪ IFNγ, n_present={s['n_ifn_present']} / {s['n_ifn_locked']}.",
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
    for r in s["sensitivity_spearman"]:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        lines.append(f"| {r['contrast']} | {r['n']} | {rho} | {pv} |")
    extra_rho = "NA" if extra.get("spearman_rho") is None else f"{round(extra['spearman_rho'], 3)}"
    extra_p = "NA" if extra.get("spearman_p") is None else f"{extra['spearman_p']:.3g}"
    lines += [
        "",
        "## Extra figure — CLDN4 / IFN / barrier along Palantir PT",
        "",
        f"Emitted: **{extra.get('emitted')}**. "
        f"Sample-level CLDN4 vs barrier ρ={extra_rho}, p={extra_p}.",
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
        "- This is not a redo of PR #459 T/NK infiltrate or PR #473 malignant IFN DE.",
        "- The pooled PT Spearman mixes two cohorts and is **not** a within-tumor progression test.",
        "- Marker-malignant on GSE123902 is EPCAM|KRT>0 and PTPRC==0 — **not CNV**.",
        "- GSE205335 is an ICI biopsy/effusion cohort, but this analysis is **not** an ICI / MPR / RECIST test.",
        "- No TACSTD2∩CLDN4 both-high gate.",
        "- Do not write “NORMAL AT2 differentiates into LUAD because PAGA is connected.”",
        "- Do not write “CLDN4 marks the malignant terminal.”",
        "- Palantir/DPT is an ordering, not a clock of real time.",
        "- GSE148071 was not added.",
        "",
        "## Outputs",
        "",
        "- `results/tables/lineage_leiden.tsv` — **done criterion (lineage)**",
        "- `results/tables/sample_level_spearman.tsv` — **done criterion (sample-level)**",
        "- `results/tables/sample_means.tsv`",
        "- `results/tables/paga_connectivities.tsv`",
        "- `results/figures/fig_trajectory_cldn4.png`",
        "- `results/figures/fig_along_pt_programs.png`",
        "- `results/figures/fig_extra_cldn4_tertile.png`",
        "- `results/figures/fig_honest_n.png`",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/pair_123902_205335_slingshot_cldn4/requirements.txt",
        "python3 methods/pair_123902_205335_slingshot_cldn4/scripts/download.py \\",
        "  --outdir /tmp/pair_123902_205335_traj",
        "python3 methods/pair_123902_205335_slingshot_cldn4/scripts/extract_epithelium.py \\",
        "  --data /tmp/pair_123902_205335_traj \\",
        "  --out /tmp/pair_123902_205335_traj/epithelium.h5ad",
        "python3 methods/pair_123902_205335_slingshot_cldn4/scripts/analyze.py \\",
        "  --input /tmp/pair_123902_205335_traj/epithelium.h5ad \\",
        "  --outdir methods/pair_123902_205335_slingshot_cldn4/results \\",
        "  --finding methods/pair_123902_205335_slingshot_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", default="methods/pair_123902_205335_slingshot_cldn4/results")
    p.add_argument("--finding", default="methods/pair_123902_205335_slingshot_cldn4/FINDING.md")
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
    if not pal.get("available"):
        raise SystemExit("Palantir must be installed for this pair (real trajectory).")

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
    ifn = ifn_genes()
    absent["IFN"] = _score_mean(adata, ifn, "score_IFN")
    n_ifn_present = len(ifn) - len(absent["IFN"])
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
    root_info["root_cldn4_tertile"] = str(adata.obs.iloc[root_i].get("cldn4_tertile", ""))
    root_info["root_CLDN4"] = float(adata.obs.iloc[root_i]["expr_CLDN4"])
    root_info["root_AT2"] = float(adata.obs.iloc[root_i]["score_AT2"])
    if root_info["root_cldn4_tertile"] == "high":
        raise SystemExit("refusing CLDN4-high root")

    pca_key = "X_pca_harmony" if harmony.get("used") and "X_pca_harmony" in adata.obsm else "X_pca"
    pal_run = run_palantir(adata, root_i, pca_key)
    pal.update(pal_run)
    adata.obs["pseudotime"] = adata.obs["palantir_pseudotime"].replace([np.inf, -np.inf], np.nan)
    dpt = adata.obs["dpt_pseudotime"].replace([np.inf, -np.inf], np.nan)
    adata.obs["dpt_pseudotime"] = dpt
    clock_label = "Palantir (Setty et al. 2019); PAGA geometry; DPT companion"
    clock_col = "pseudotime"

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
                "n_GSE123902": int((sub["dataset"] == "GSE123902").sum()),
                "n_GSE205335": int((sub["dataset"] == "GSE205335").sum()),
                "n_normal": int((sub["Sample_Origin"] == "NORMAL").sum()),
                "top_subtype": sub["author_subtype"].astype(str).value_counts().index[0],
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_palantir": float(sub[clock_col].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
            }
        )
    lineage = pd.DataFrame(cluster_tab).sort_values("mean_palantir")
    lineage.to_csv(tabdir / "lineage_leiden.tsv", sep="\t", index=False)
    lineage.to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    pd.DataFrame(connect, index=leiden_ids, columns=leiden_ids).to_csv(
        tabdir / "paga_connectivities.tsv", sep="\t"
    )

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
                "is_tumor_unit": str(sub["is_tumor_unit"].iloc[0]),
                "n_cells": int(len(sub)),
                "n_cldn4_low": int((sub["cldn4_tertile"] == "low").sum()),
                "n_cldn4_mid": int((sub["cldn4_tertile"] == "mid").sum()),
                "n_cldn4_high": int((sub["cldn4_tertile"] == "high").sum()),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_TACSTD2": float(sub["expr_TACSTD2"].mean()),
                "mean_pt": float(sub[clock_col].mean()),
                "mean_dpt": float(sub["dpt_pseudotime"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_club": float(sub["score_club"].mean()),
                "mean_basal": float(sub["score_basal"].mean()),
                "mean_barrier_keratin": float(sub["score_barrier_keratin"].mean()),
                "mean_malignant_like": float(sub["score_malignant_like"].mean()),
                "pct_CLDN4_pos": float((sub["expr_CLDN4"] > 0).mean()),
                "mean_SFTPC": float(sub["expr_SFTPC"].mean()) if "expr_SFTPC" in sub else np.nan,
            }
        )
    sample_df = pd.DataFrame(rows)
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)
    tumor = sample_df[sample_df["is_tumor_unit"] == "True"].copy()
    elig = tumor[tumor["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN].copy()

    contrasts = [
        ("CLDN4 vs Palantir PT", "mean_CLDN4", "mean_pt"),
        ("IFN vs Palantir PT", "mean_IFN", "mean_pt"),
        ("barrier/keratin (no CLDN4) vs Palantir PT", "mean_barrier_keratin", "mean_pt"),
        ("CLDN4 vs IFN", "mean_CLDN4", "mean_IFN"),
        ("CLDN4 vs barrier/keratin (no CLDN4)", "mean_CLDN4", "mean_barrier_keratin"),
        ("CLDN4 vs AT2 score", "mean_CLDN4", "mean_AT2"),
        ("CLDN4 vs TACSTD2 (comparator)", "mean_CLDN4", "mean_TACSTD2"),
        ("SFTPC vs Palantir PT (control)", "mean_SFTPC", "mean_pt"),
        ("AT2 score vs Palantir PT (control)", "mean_AT2", "mean_pt"),
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

    elig_123 = elig[elig["dataset"] == "GSE123902"]
    elig_205 = elig[elig["dataset"] == "GSE205335"]
    elig_adcsq = elig_205[elig_205["histology"].isin(["ADC", "SQ"])]
    all_elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE_FOR_MEAN]

    def sp(frame, a, b):
        if frame is None or len(frame) == 0:
            return {"n": 0, "rho": None, "p": None}
        return _spearman(frame[a].to_numpy(), frame[b].to_numpy())

    sensitivity = [
        {"contrast": "GSE123902-only CLDN4 vs Palantir PT", **sp(elig_123, "mean_CLDN4", "mean_pt")},
        {"contrast": "GSE205335-only CLDN4 vs Palantir PT", **sp(elig_205, "mean_CLDN4", "mean_pt")},
        {"contrast": "GSE123902-only IFN vs Palantir PT", **sp(elig_123, "mean_IFN", "mean_pt")},
        {"contrast": "GSE205335-only IFN vs Palantir PT", **sp(elig_205, "mean_IFN", "mean_pt")},
        {"contrast": "GSE123902-only barrier vs Palantir PT", **sp(elig_123, "mean_barrier_keratin", "mean_pt")},
        {"contrast": "GSE205335-only barrier vs Palantir PT", **sp(elig_205, "mean_barrier_keratin", "mean_pt")},
        {"contrast": "GSE123902-only CLDN4 vs IFN", **sp(elig_123, "mean_CLDN4", "mean_IFN")},
        {"contrast": "GSE205335-only CLDN4 vs IFN", **sp(elig_205, "mean_CLDN4", "mean_IFN")},
        {"contrast": "include-NORMAL-units CLDN4 vs Palantir PT", **sp(all_elig, "mean_CLDN4", "mean_pt")},
        {"contrast": "GSE205335 ADC+SQ CLDN4 vs Palantir PT", **sp(elig_adcsq, "mean_CLDN4", "mean_pt")},
        {"contrast": "tumor CLDN4 vs DPT (companion)", **sp(elig, "mean_CLDN4", "mean_dpt")},
        {"contrast": "tumor IFN vs DPT (companion)", **sp(elig, "mean_IFN", "mean_dpt")},
    ]
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    paired_recs = []
    for unit, sub in adata.obs.groupby("unit_id", observed=True):
        if str(sub["is_tumor_unit"].iloc[0]) != "True":
            continue
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
                "pt_high": float(hi[clock_col].mean()),
                "pt_low": float(lo[clock_col].mean()),
            }
        )
    paired_df = pd.DataFrame(paired_recs)
    if not paired_df.empty:
        paired_df.to_csv(tabdir / "cldn4_tertile_paired.tsv", sep="\t", index=False)
    paired_rows = []
    for label, a, b in (
        ("barrier/keratin (no CLDN4) high vs low", "barrier_high", "barrier_low"),
        ("IFN high vs low", "IFN_high", "IFN_low"),
        ("AT2 high vs low", "AT2_high", "AT2_low"),
        ("Palantir PT high vs low", "pt_high", "pt_low"),
    ):
        if paired_df.empty:
            paired_rows.append({"contrast": label, "n": 0, "W": None, "p": None, "delta_median": None})
        else:
            paired_rows.append({"contrast": label, **_wilcoxon_paired(paired_df[a].to_numpy(), paired_df[b].to_numpy())})

    barrier_row = next(r for r in primary if "barrier/keratin (no CLDN4)" in r["contrast"] and "vs Palantir" not in r["contrast"])
    emit_extra = True

    colors = {"GSE123902": "#3d5a80", "GSE205335": "#b23a48"}

    # ---- figures ----
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.0))
    sc.pl.paga(
        adata,
        color="expr_CLDN4",
        ax=axes[0, 0],
        show=False,
        frameon=False,
        cmap="viridis",
        title="PAGA (Leiden) colored by mean CLDN4",
    )
    sc.pl.umap(adata, color="expr_CLDN4", ax=axes[0, 1], show=False, frameon=False, cmap="viridis", title="UMAP CLDN4")
    sc.pl.umap(adata, color=clock_col, ax=axes[1, 0], show=False, frameon=False, cmap="magma", title="UMAP Palantir PT")
    c4_pt = next(r for r in primary if r["contrast"] == "CLDN4 vs Palantir PT")
    for ds, col in colors.items():
        sub = elig[elig["dataset"] == ds]
        axes[1, 1].scatter(sub["mean_pt"], sub["mean_CLDN4"], s=48, c=col, label=f"{ds} n={len(sub)}")
    axes[1, 1].set_xlabel("sample-mean Palantir PT")
    axes[1, 1].set_ylabel("sample-mean CLDN4")
    axes[1, 1].set_title(_fmt(c4_pt))
    axes[1, 1].legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"Pair GSE123902+GSE205335 epithelium scored by CLDN4   n_cells={adata.n_obs}  n_tumor={len(elig)}",
        fontsize=11,
    )
    _save(fig, figdir / "fig_trajectory_cldn4")

    # extra: CLDN4 + IFN + barrier along PT (cell bins + sample scatter)
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    pt = adata.obs[clock_col].to_numpy()
    finite = np.isfinite(pt)
    bins = pd.qcut(pt[finite], 12, duplicates="drop")
    for ax, key, lab in (
        (axes[0], "expr_CLDN4", "CLDN4"),
        (axes[1], "score_IFN", "IFN (Hallmark α∪γ)"),
        (axes[2], "score_barrier_keratin", "barrier/keratin (no CLDN4)"),
    ):
        y = adata.obs[key].to_numpy()[finite]
        tmp = pd.DataFrame({"bin": bins, "y": y})
        g = tmp.groupby("bin", observed=True)["y"].mean()
        ax.plot(np.arange(len(g)), g.to_numpy(), marker="o", c="#1b4965")
        ax.set_xlabel("Palantir PT bin")
        ax.set_ylabel(f"mean {lab}")
        ax.set_title(lab)
    fig.suptitle("EXTRA: CLDN4 + IFN + barrier along Palantir PT (cell bins)", fontsize=11)
    _save(fig, figdir / "fig_along_pt_programs")

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    ifn_pt = next(r for r in primary if r["contrast"] == "IFN vs Palantir PT")
    bar_pt = next(r for r in primary if r["contrast"].startswith("barrier/keratin (no CLDN4) vs"))
    for ax, x, y, row, xlab, ylab in (
        (axes[0], "mean_pt", "mean_CLDN4", c4_pt, "sample-mean Palantir PT", "sample-mean CLDN4"),
        (axes[1], "mean_pt", "mean_IFN", ifn_pt, "sample-mean Palantir PT", "sample-mean IFN"),
        (axes[2], "mean_pt", "mean_barrier_keratin", bar_pt, "sample-mean Palantir PT", "sample-mean barrier"),
    ):
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub[y], s=48, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: sample-level CLDN4 / IFN / barrier vs Palantir PT", fontsize=11)
    _save(fig, figdir / "fig_extra_sample_along_pt")

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    ct = (
        adata.obs.assign(author_subtype=adata.obs["author_subtype"].astype(str).fillna("NA"))
        .groupby(["dataset", "author_subtype"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct.T.plot(kind="bar", ax=axes[0], color={"GSE123902": "#3d5a80", "GSE205335": "#b23a48"})
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Subtypes  n_cells={adata.n_obs}")
    axes[0].legend(frameon=False, fontsize=8)
    plt.setp(axes[0].get_xticklabels(), rotation=40, ha="right")
    unit_ct = elig.groupby("dataset").size()
    axes[1].bar(unit_ct.index.astype(str), unit_ct.to_numpy(), color=["#3d5a80", "#b23a48"][: len(unit_ct)])
    axes[1].set_ylabel("tumor units")
    axes[1].set_title(f"Honest n tumor units={len(elig)} (catalog tumor {len(tumor)})")
    _save(fig, figdir / "fig_honest_n")

    if emit_extra and not paired_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8))
        panels = (
            ("barrier_low", "barrier_high", "barrier/keratin (no CLDN4)"),
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
        axes[0].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("barrier")), keys=("W", "p")))
        axes[1].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("IFN")), keys=("W", "p")))
        axes[2].set_title(_fmt(next(r for r in paired_rows if r["contrast"].startswith("Palantir")), keys=("W", "p")))
        axes[0].legend(fontsize=7, frameon=False)
        fig.suptitle(f"EXTRA: within-unit CLDN4-high vs low  paired n={len(paired_df)}", fontsize=11)
        _save(fig, figdir / "fig_extra_cldn4_tertile")

    for color, fname, cmap in (
        ("dataset", "fig_umap_dataset", None),
        ("Sample_Origin", "fig_umap_origin", None),
        ("author_subtype", "fig_umap_author_subtype", None),
        ("cldn4_tertile", "fig_umap_cldn4_tertile", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("score_IFN", "fig_umap_IFN", "viridis"),
        ("score_barrier_keratin", "fig_umap_barrier_no_cldn4", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    sc.pl.paga(adata, color="score_IFN", ax=ax, show=False, frameon=False, cmap="viridis", title="PAGA IFN")
    _save(fig, figdir / "fig_paga_ifn")

    subtype_counts = adata.obs["author_subtype"].astype(str).fillna("NA").value_counts().to_dict()
    tertile_counts = adata.obs["cldn4_tertile"].astype(str).value_counts().to_dict()
    c4_ifn = next(r for r in primary if r["contrast"] == "CLDN4 vs IFN")
    c4_bar = barrier_row
    pair_bar = next(r for r in paired_rows if r["contrast"].startswith("barrier"))
    pair_ifn = next(r for r in paired_rows if r["contrast"].startswith("IFN"))
    s123 = next(r for r in sensitivity if r["contrast"] == "GSE123902-only CLDN4 vs Palantir PT")
    s205 = next(r for r in sensitivity if r["contrast"] == "GSE205335-only CLDN4 vs Palantir PT")

    parts = [
        f"Sample-level CLDN4 vs Palantir PT: {_fmt(c4_pt)}.",
        f"IFN vs Palantir PT: {_fmt(ifn_pt)}.",
        f"barrier/keratin (CLDN4 excluded) vs Palantir PT: {_fmt(bar_pt)}.",
        f"CLDN4 vs IFN: {_fmt(c4_ifn)}.",
        f"CLDN4 vs barrier/keratin (CLDN4 excluded): {_fmt(c4_bar)}.",
        f"GSE123902-only CLDN4 vs PT: {_fmt(s123)}.",
        f"GSE205335-only CLDN4 vs PT: {_fmt(s205)}.",
        f"Paired CLDN4-high vs low barrier/keratin: {_fmt(pair_bar, keys=('W', 'p'))}.",
        f"Paired CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}.",
        f"PAGA has {len(comps)} component(s) at connectivity>0 among {len(leiden_ids)} Leiden vertices.",
        "Root is GSE123902 NORMAL AT2-like and is not CLDN4-high.",
        "Given PR #459 T/NK and PR #473 IFN DE were not re-audited. No both-high gate. GSE148071 not added.",
    ]
    verdict = " ".join(parts)
    what_holds = (
        f"**What holds (n={c4_bar['n']} tumor units).** "
        f"CLDN4 vs Palantir PT is {_fmt(c4_pt)}. "
        f"IFN vs Palantir PT is {_fmt(ifn_pt)}. "
        f"barrier/keratin (no CLDN4) vs Palantir PT is {_fmt(bar_pt)}. "
        f"CLDN4 vs IFN is {_fmt(c4_ifn)}. "
        f"CLDN4 vs barrier is {_fmt(c4_bar)}. "
        f"Within-unit CLDN4-high vs low IFN: {_fmt(pair_ifn, keys=('W', 'p'))}; "
        f"barrier: {_fmt(pair_bar, keys=('W', 'p'))}. "
        f"**What this is not.** Given T/NK n=35 and IFN DE n=18 are different contrasts. "
        "Pooled PT mixes cohorts and is not a within-tumor progression test."
    )

    summary = {
        "accessions": ["GSE123902", "GSE205335"],
        "gse148071_added": False,
        "dual_high": False,
        "tnk_re_audited": False,
        "given_pr459": GIVEN,
        "given_pr473_ifn_de": IFN_DE,
        "primary_gene": "CLDN4",
        "clock": "palantir",
        "clock_label": clock_label,
        "slingshot": sling,
        "palantir": pal,
        "harmony": harmony,
        "cap_per_unit": 350,
        "n_cells": int(adata.n_obs),
        "n_cells_gse123902": int((adata.obs["dataset"] == "GSE123902").sum()),
        "n_cells_gse205335": int((adata.obs["dataset"] == "GSE205335").sum()),
        "n_cells_normal": int((adata.obs["Sample_Origin"] == "NORMAL").sum()),
        "n_normal_units": int(adata.obs.loc[adata.obs["Sample_Origin"] == "NORMAL", "unit_id"].nunique()),
        "n_tumor_units": int(len(tumor)),
        "n_tumor_gse123902": int((tumor["dataset"] == "GSE123902").sum()),
        "n_tumor_gse205335": int((tumor["dataset"] == "GSE205335").sum()),
        "n_units_eligible": int(len(elig)),
        "n_units_paired_tertile": int(len(paired_df)),
        "subtype_counts": subtype_counts,
        "cldn4_tertile_counts": tertile_counts,
        "cldn4_tertile_cuts": {"q1": float(q1), "q2": float(q2)},
        "genes_absent": absent,
        "n_ifn_locked": int(len(ifn)),
        "n_ifn_present": int(n_ifn_present),
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
        "verdict": verdict,
        "what_holds": what_holds,
        "skipped": {
            "GSE148071": "explicitly not added",
            "dual_high": "TACSTD2∩CLDN4 not used",
            "T/NK": "PR #459 given; not re-audited",
            "IFN_DE": "PR #473 given; not re-run",
            "slingshot": sling.get("reason"),
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(Path(args.finding), {"summary": summary})
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accessions": ["GSE123902", "GSE205335"],
                "papers": {
                    "GSE123902": "Laughney et al. Nat Med 2020 PMID 32066974 (GEO GSE123902)",
                    "GSE205335": "Ahn / Lee et al. eLife 2024 (GEO GSE205335)",
                },
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "ifn_set": "Hallmark IFNα ∪ IFNγ",
                "slingshot": sling,
                "palantir": pal,
                "locked": {
                    "leiden_resolution": LEIDEN_RES,
                    "n_hvg": N_HVG,
                    "n_neighbors": N_NEIGHBORS,
                    "n_pcs": N_PCS,
                    "root": "GSE123902 NORMAL AT2-like, not CLDN4-high",
                    "cap_per_unit": 350,
                    "unit": "GSE123902 donor + GSE205335 patient",
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
                "n_tumor_units": int(len(tumor)),
                "n_eligible": int(len(elig)),
                "root": root_info,
                "finding": args.finding,
                "lineage_table": str(tabdir / "lineage_leiden.tsv"),
                "sample_table": str(tabdir / "sample_level_spearman.tsv"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
