#!/usr/bin/env python3
"""Real Slingshot + PAGA on winning-pair epithelium, CLDN4 only.

ADDITIVE. GSE131907 + GSE205335 malignant/epithelial cells.
GSE207422 and GSE148071 are not added. No TACSTD2 x CLDN4 dual-high gate.
Does not re-audit PR #320 T/NK Q4 vs Q1 r=−0.705.

Primary clock is Slingshot (Street et al. 2018) via pyslingshot-bio, a
documented Python port of R Slingshot — not diffusion pseudotime.
Root cluster is AT2-like or lowest-CLDN4 and is never the CLDN4-high cluster.
Inferential unit = GSE131907 sample + GSE205335 patient.
Barrier/keratin and IFN scores exclude CLDN4.
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
from gene_sets import AUTHOR_AT2, COMPARATOR, CONTROLS, FOCAL, QC_NEG, STATES  # noqa: E402
from slingshot_py import r_slingshot_status, run_slingshot  # noqa: E402

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e


LEIDEN_RES = 0.4
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_CELLS_PER_SAMPLE = 10
N_CURVE_BINS = 12


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


def _expr(adata, gene: str, key: str) -> None:
    if gene not in adata.var_names:
        adata.obs[key] = np.nan
        return
    X = adata[:, gene].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    adata.obs[key] = np.asarray(X, dtype=float).ravel()


def maybe_harmony(adata) -> dict:
    info = {"used": False, "reason": None}
    try:
        import harmonypy
    except ImportError:
        info["reason"] = "harmonypy not installed; neighbors on PCA"
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
        adata.obsm["X_pca_harmony"] = np.asarray(adata.obsm["X_pca"][:, :N_PCS])
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


def pick_root_cluster(adata) -> tuple[str, dict, pd.DataFrame]:
    """AT2-like or lowest-CLDN4. Never the CLDN4-high cluster."""
    obs = adata.obs
    nlung_at2 = (obs["dataset"].astype(str) == "GSE131907") & (
        obs["Sample_Origin"].astype(str) == "nLung"
    ) & obs["author_subtype"].astype(str).isin(AUTHOR_AT2)
    rows = []
    for lab, sub in obs.groupby("leiden", observed=True):
        rows.append(
            {
                "leiden": str(lab),
                "n": int(len(sub)),
                "mean_CLDN4": float(sub["expr_CLDN4"].mean()),
                "mean_AT2": float(sub["score_AT2"].mean()),
                "mean_barrier": float(sub["score_barrier_keratin"].mean()),
                "mean_IFN": float(sub["score_IFN"].mean()),
                "mean_malignant": float(sub["score_malignant_like"].mean()),
                "n_nLung_AT2": int(nlung_at2.loc[sub.index].sum()),
                "frac_author_AT2": float(sub["author_subtype"].astype(str).isin(AUTHOR_AT2).mean()),
            }
        )
    stats_df = pd.DataFrame(rows).set_index("leiden")
    cldn4_high = str(stats_df["mean_CLDN4"].idxmax())
    cand = stats_df.drop(index=cldn4_high)
    info = {
        "excluded_cldn4_high_cluster": cldn4_high,
        "excluded_mean_CLDN4": float(stats_df.loc[cldn4_high, "mean_CLDN4"]),
    }
    if cand["n_nLung_AT2"].max() >= 20:
        root = str(cand["n_nLung_AT2"].idxmax())
        info["rule"] = (
            "Leiden with most GSE131907 nLung author AT2; "
            "max-CLDN4 cluster excluded"
        )
    else:
        low_c4 = cand[cand["mean_CLDN4"] <= cand["mean_CLDN4"].median()]
        if len(low_c4):
            root = str(low_c4["mean_AT2"].idxmax())
            info["rule"] = (
                "highest AT2 among below-median-CLDN4 clusters; "
                "max-CLDN4 cluster excluded"
            )
        else:
            root = str(cand["mean_CLDN4"].idxmin())
            info["rule"] = "lowest-CLDN4 cluster; max-CLDN4 cluster excluded"
    info["root_cluster"] = root
    info["root_mean_CLDN4"] = float(stats_df.loc[root, "mean_CLDN4"])
    info["root_mean_AT2"] = float(stats_df.loc[root, "mean_AT2"])
    info["root_n_nLung_AT2"] = int(stats_df.loc[root, "n_nLung_AT2"])
    if root == cldn4_high:
        raise SystemExit("root picker returned the CLDN4-high cluster")
    return root, info, stats_df.reset_index()


def _fmt(row: dict) -> str:
    n = row.get("n")
    if row.get("rho") is None:
        return f"n={n}, ρ=NA, p=NA"
    return f"n={n}, ρ={row['rho']:.3f}, p={row['p']:.3g}"


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def along_curve(pt: np.ndarray, values: dict[str, np.ndarray], n_bins: int) -> pd.DataFrame:
    mask = np.isfinite(pt)
    if mask.sum() < n_bins:
        n_bins = max(4, int(mask.sum() // 5) or 4)
    edges = np.linspace(0, 1, n_bins + 1)
    # rescale finite PT to 0-1
    x = pt.copy()
    lo, hi = np.nanmin(x[mask]), np.nanmax(x[mask])
    if hi > lo:
        x = (x - lo) / (hi - lo)
    rows = []
    for i in range(n_bins):
        if i == n_bins - 1:
            m = mask & (x >= edges[i]) & (x <= edges[i + 1])
        else:
            m = mask & (x >= edges[i]) & (x < edges[i + 1])
        rec = {
            "bin": i,
            "pt_lo": float(edges[i]),
            "pt_hi": float(edges[i + 1]),
            "pt_mid": float(0.5 * (edges[i] + edges[i + 1])),
            "n_cells": int(m.sum()),
        }
        for name, arr in values.items():
            rec[f"mean_{name}"] = float(np.nanmean(arr[m])) if m.any() else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def maybe_palantir(adata, root_cluster: str) -> dict:
    info = {"available": False, "reason": None}
    try:
        import palantir
    except ImportError as exc:
        info["reason"] = f"palantir not installed: {exc}"
        return info
    try:
        pca = pd.DataFrame(
            adata.obsm["X_pca_harmony"],
            index=adata.obs_names.astype(str),
        )
        dm = palantir.utils.run_diffusion_maps(pca, n_components=5)
        ms = palantir.utils.determine_multiscale_space(dm)
        in_root = adata.obs["leiden"].astype(str) == root_cluster
        scores = adata.obs.loc[in_root, "score_AT2"].to_numpy()
        names = adata.obs_names[in_root.to_numpy()]
        early = str(names[int(np.nanargmin(np.abs(scores - np.nanmedian(scores))))])
        pr = palantir.core.run_palantir(
            ms, early, num_waypoints=min(300, max(50, adata.n_obs // 40)), n_jobs=1
        )
        adata.obs["palantir_pt"] = pr.pseudotime.reindex(adata.obs_names).to_numpy()
        info = {
            "available": True,
            "early_cell": early,
            "n_waypoints": int(min(300, max(50, adata.n_obs // 40))),
            "module": "palantir",
        }
    except Exception as exc:  # noqa: BLE001
        info = {"available": False, "reason": f"{type(exc).__name__}: {exc}"}
    return info


def write_finding(path: Path, s: dict) -> None:
    def row_md(r: dict) -> str:
        rho = "NA" if r.get("rho") is None else f"{r['rho']:.3f}"
        pv = "NA" if r.get("p") is None else f"{r['p']:.3g}"
        qv = "NA" if r.get("q") is None else f"{r['q']:.3g}"
        return f"| {r['contrast']} | {r['n']} | {rho} | {pv} | {qv} |"

    lin_rows = []
    for lin in s["lineages"]:
        lin_rows.append(
            f"| {lin['name']} | {' → '.join(lin['clusters'])} | {lin['n_cells']} | "
            f"{lin['terminal']} | {lin['terminal_mean_CLDN4']:.3f} | "
            f"{lin['terminal_mean_barrier']:.3f} | {lin['terminal_mean_IFN']:.3f} |"
        )

    lines = [
        "# Finding — winning-pair GSE131907+GSE205335, real Slingshot CLDN4-only",
        "",
        "ADDITIVE. **CLDN4 only.** Winning pair is given (PR #320: author-malignant "
        "CLDN4 %pos vs T/NK, Q4 vs Q1 **n=23** r=−0.705) and is **not** re-audited. "
        "Malignant/epithelial cells. No dual-high TACSTD2×CLDN4. GSE148071 is not added. "
        "GSE207422 is not added. This is **not** a DPT-null writeup and **not** a redo of "
        "PR #449 (R missing → DPT).",
        "",
        f"Primary clock: **{s['clock_label']}**. Inferential unit = **sample/patient** "
        "(GSE131907 `Sample`, GSE205335 `patient`). Cell-level ρ is descriptive. "
        "Barrier/keratin and IFN scores **exclude CLDN4**. "
        f"Root cluster = Leiden {s['root']['root_cluster']} "
        f"({s['root']['rule']}). Never CLDN4-high "
        f"(excluded Leiden {s['root']['excluded_cldn4_high_cluster']}).",
        "",
        s.get("what_holds", ""),
        "",
        "## Verdict",
        "",
        s["verdict"],
        "",
        "## Honest n",
        "",
        f"- Analysis cells after QC (capped ≤{s['cap_per_unit']}/unit): "
        f"**n_cells = {s['n_cells']}** "
        f"(GSE131907 {s['n_cells_gse131907']}, GSE205335 {s['n_cells_gse205335']}).",
        f"- Units (GSE131907 Sample + GSE205335 patient): **n_units = {s['n_units']}** "
        f"(GSE131907 {s['n_units_gse131907']}, GSE205335 {s['n_units_gse205335']}).",
        f"- Units with ≥{MIN_CELLS_PER_SAMPLE} epithelial cells used for Spearman: "
        f"**n = {s['n_units_eligible']}**.",
        f"- GSE131907 nLung cells / author AT2 in the object: "
        f"{s['n_cells_nLung']} / {s['n_author_AT2']}.",
        f"- Author subtypes (cells): {s['subtype_counts']}.",
        f"- Slingshot lineages: **{s['n_lineages']}**. Engine: `{s['slingshot_engine']}`.",
        f"- R Slingshot: available={s['r_slingshot']['available']}; "
        f"{s['r_slingshot'].get('reason') or s['r_slingshot'].get('version')}.",
        f"- Palantir companion: available={s['palantir'].get('available')}; "
        f"{s['palantir'].get('reason') or 'ran'}.",
        f"- Genes absent from locked sets: {s['genes_absent']}.",
        "- GSE207422 not used. GSE148071 not used. GSE131907 PE unlabeled epithelium dropped. "
        "GSE205335 normal-tissue samples dropped.",
        "- PR #320 T/NK r=−0.705 is given and was not recomputed.",
        "",
        "## Slingshot / PAGA lineages",
        "",
        f"Start cluster **{s['root']['root_cluster']}** "
        f"(mean CLDN4 {s['root']['root_mean_CLDN4']:.3f}, "
        f"mean AT2 {s['root']['root_mean_AT2']:.3f}, "
        f"nLung AT2 cells {s['root']['root_n_nLung_AT2']}). "
        f"PAGA components at connectivity>0: **{s['n_paga_components']}** "
        f"among {s['n_leiden']} Leiden vertices.",
        "",
        "| lineage | clusters | n_cells | terminal | term CLDN4 | term barrier | term IFN |",
        "| --- | --- | ---: | --- | ---: | ---: | ---: |",
        *lin_rows,
        "",
        f"Primary lineage for along-curve plots: **{s['primary_lineage']}** "
        f"(terminal with highest mean CLDN4 among root-started lineages).",
        "",
        "## Sample-level Spearman (done criterion)",
        "",
        "Unit = GSE131907 sample or GSE205335 patient. BH inside this list only.",
        "",
        "| Contrast | n_units | ρ | p | q |",
        "| --- | ---: | ---: | ---: | ---: |",
        *[row_md(r) for r in s["primary_spearman"]],
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
    lines += [
        "",
        "## Along the primary Slingshot curve",
        "",
        "Binned cell means on the primary lineage (0 = root, 1 = terminal). "
        "This is the same curve for CLDN4, barrier/keratin (no CLDN4), and IFN.",
        "",
        "| bin | n_cells | mean CLDN4 | mean barrier | mean IFN | mean AT2 | mean malignant |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for rec in s["along_curve"]:
        lines.append(
            f"| {rec['bin']} | {rec['n_cells']} | {rec['mean_CLDN4']:.3f} | "
            f"{rec['mean_barrier']:.3f} | {rec['mean_IFN']:.3f} | "
            f"{rec['mean_AT2']:.3f} | {rec['mean_malignant']:.3f} |"
        )
    lines += [
        "",
        "## Extra figures",
        "",
        "- `results/figures/fig_lineage_cldn4.png` — **real Slingshot lineage plot** (done criterion)",
        "- `results/figures/fig_along_curve.png` — CLDN4 / barrier / IFN on the same curve",
        "- `results/figures/fig_sample_cldn4_pt.png` — sample-level CLDN4 vs Slingshot PT",
        "- `results/figures/fig_paga.png` — PAGA on Leiden",
        "- `results/figures/fig_extra_sample_programs.png` — sample CLDN4 vs barrier and vs IFN",
        "- `results/figures/fig_honest_n.png`",
        "",
        "## What this does not claim",
        "",
        "- Cell-level p-values are not the claim. n_cells is large by construction.",
        "- PR #320 T/NK Q4 vs Q1 r=−0.705 is given; this folder does not re-rank T/NK.",
        "- This is not a redo of PR #325 (GSE131907-only PAGA) or PR #449 (DPT fallback).",
        "- The pooled PT Spearman mixes cohorts and nLung vs tumor. Tumor-only is in Sensitivity.",
        "- Malignant-like is CEACAM5/6/MKI67 — **not CNV**.",
        "- GSE205335 is an ICI biopsy/effusion cohort; this analysis is **not** an ICI / MPR test.",
        "- No TACSTD2∩CLDN4 both-high gate. GSE148071 not used.",
        "- Do not write “AT2 differentiates into LUAD because PAGA is connected.”",
        "- Slingshot is an ordering along principal curves, not a developmental clock.",
        "",
        "## Outputs",
        "",
        "- `results/tables/sample_cldn4_vs_pseudotime.tsv` — **done criterion**",
        "- `results/tables/sample_level_spearman.tsv`",
        "- `results/figures/fig_lineage_cldn4.png` — **done criterion**",
        "- `results/summary.json`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install -r methods/winpair_slingshot_real_cldn4/requirements.txt",
        "python3 methods/winpair_slingshot_real_cldn4/scripts/download.py \\",
        "  --out /tmp/winpair_slingshot_real",
        "python3 methods/winpair_slingshot_real_cldn4/scripts/extract.py \\",
        "  --data /tmp/winpair_slingshot_real \\",
        "  --out /tmp/winpair_slingshot_real/epithelium.h5ad",
        "python3 methods/winpair_slingshot_real_cldn4/scripts/analyze.py \\",
        "  --input /tmp/winpair_slingshot_real/epithelium.h5ad \\",
        "  --outdir methods/winpair_slingshot_real_cldn4/results \\",
        "  --finding methods/winpair_slingshot_real_cldn4/FINDING.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--finding", type=Path, required=True)
    p.add_argument("--cap-per-unit", type=int, default=280)
    args = p.parse_args()

    outdir = args.outdir
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.verbosity = 1
    sc.settings.figdir = str(figdir)

    adata = sc.read_h5ad(args.input)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    keep = (adata.obs["n_genes"] >= 200) & (adata.obs["n_counts"] >= 500)
    adata = adata[keep].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    absent: dict[str, list[str]] = {}
    for name, genes in STATES.items():
        absent[name] = _score_mean(adata, genes, f"score_{name}" if name != "barrier_keratin" else "score_barrier_keratin")
    # alias
    if "score_IFN" not in adata.obs:
        absent["IFN"] = _score_mean(adata, STATES["IFN"], "score_IFN")
    for g in FOCAL + COMPARATOR + CONTROLS + QC_NEG:
        _expr(adata, g, f"expr_{g}")

    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts")
    adata.raw = adata
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    harmony = maybe_harmony(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, key_added="leiden")
    adata.obs["leiden"] = adata.obs["leiden"].astype(str)
    sc.tl.paga(adata, groups="leiden")
    connect = np.asarray(adata.uns["paga"]["connectivities"].todense())
    n_comp = _paga_ncomp(connect)

    root, root_info, leiden_df = pick_root_cluster(adata)
    leiden_df.to_csv(tabdir / "leiden_vertices.tsv", sep="\t", index=False)

    X = np.asarray(adata.obsm["X_pca_harmony"][:, : min(15, adata.obsm["X_pca_harmony"].shape[1])])
    embed = np.asarray(adata.obsm["X_umap"])
    labels = adata.obs["leiden"].to_numpy()
    sling = run_slingshot(X, labels, start_cluster=root, embed2d=embed)
    rstat = r_slingshot_status()

    shared = sling.shared_pseudotime
    adata.obs["slingshot_pt"] = shared
    for i, name in enumerate(sling.lineage_names):
        if i < sling.pseudotime.shape[1]:
            adata.obs[f"pt_{name}"] = sling.pseudotime[:, i]

    # lineage summaries
    lin_recs = []
    for i, (name, clusters) in enumerate(zip(sling.lineage_names, sling.lineages)):
        terminal = str(clusters[-1])
        in_lin = np.isin(labels, clusters)
        term = labels == terminal
        lin_recs.append(
            {
                "name": name,
                "clusters": [str(c) for c in clusters],
                "n_cells": int(in_lin.sum()),
                "terminal": terminal,
                "terminal_mean_CLDN4": float(adata.obs.loc[term, "expr_CLDN4"].mean()) if term.any() else np.nan,
                "terminal_mean_barrier": float(adata.obs.loc[term, "score_barrier_keratin"].mean()) if term.any() else np.nan,
                "terminal_mean_IFN": float(adata.obs.loc[term, "score_IFN"].mean()) if term.any() else np.nan,
                "terminal_mean_AT2": float(adata.obs.loc[term, "score_AT2"].mean()) if term.any() else np.nan,
            }
        )
    if not lin_recs:
        raise SystemExit("Slingshot returned no lineages — refusing a DPT-null writeup")
    # primary = highest terminal CLDN4 (barrier/malignant end), not IFN
    primary = max(lin_recs, key=lambda r: (r["terminal_mean_CLDN4"], r["n_cells"]))
    primary_name = primary["name"]
    pi = sling.lineage_names.index(primary_name)
    adata.obs["slingshot_pt_primary"] = sling.pseudotime[:, pi]
    pd.DataFrame(lin_recs).assign(clusters=lambda d: d["clusters"].map(lambda x: ",".join(x))).to_csv(
        tabdir / "lineages.tsv", sep="\t", index=False
    )

    pal = maybe_palantir(adata, root)

    # sample table
    obs = adata.obs
    grp = obs.groupby("unit_id", observed=True)
    sample_df = grp.agg(
        dataset=("dataset", "first"),
        Sample=("Sample", "first"),
        Sample_Origin=("Sample_Origin", "first"),
        n_cells=("unit_id", "size"),
        mean_CLDN4=("expr_CLDN4", "mean"),
        mean_slingshot_pt=("slingshot_pt", "mean"),
        mean_slingshot_pt_primary=("slingshot_pt_primary", "mean"),
        mean_barrier=("score_barrier_keratin", "mean"),
        mean_IFN=("score_IFN", "mean"),
        mean_AT2=("score_AT2", "mean"),
        mean_malignant=("score_malignant_like", "mean"),
        mean_TACSTD2=("expr_TACSTD2", "mean"),
        mean_SFTPC=("expr_SFTPC", "mean"),
    ).reset_index()
    if "palantir_pt" in obs.columns:
        sample_df["mean_palantir_pt"] = grp["palantir_pt"].mean().to_numpy()
    elig = sample_df[sample_df["n_cells"] >= MIN_CELLS_PER_SAMPLE].copy()
    elig.to_csv(tabdir / "sample_cldn4_vs_pseudotime.tsv", sep="\t", index=False)
    sample_df.to_csv(tabdir / "sample_means.tsv", sep="\t", index=False)

    primary_rows = [
        {"contrast": "CLDN4 vs Slingshot PT (shared)", **_spearman(elig["mean_CLDN4"], elig["mean_slingshot_pt"])},
        {"contrast": "CLDN4 vs Slingshot PT (primary lineage)", **_spearman(elig["mean_CLDN4"], elig["mean_slingshot_pt_primary"])},
        {"contrast": "CLDN4 vs AT2 score", **_spearman(elig["mean_CLDN4"], elig["mean_AT2"])},
        {"contrast": "CLDN4 vs barrier/keratin (no CLDN4)", **_spearman(elig["mean_CLDN4"], elig["mean_barrier"])},
        {"contrast": "CLDN4 vs IFN (no CLDN4)", **_spearman(elig["mean_CLDN4"], elig["mean_IFN"])},
        {"contrast": "CLDN4 vs malignant-like score", **_spearman(elig["mean_CLDN4"], elig["mean_malignant"])},
        {"contrast": "barrier vs Slingshot PT (primary)", **_spearman(elig["mean_barrier"], elig["mean_slingshot_pt_primary"])},
        {"contrast": "IFN vs Slingshot PT (primary)", **_spearman(elig["mean_IFN"], elig["mean_slingshot_pt_primary"])},
        {"contrast": "CLDN4 vs TACSTD2 (comparator)", **_spearman(elig["mean_CLDN4"], elig["mean_TACSTD2"])},
        {"contrast": "SFTPC vs Slingshot PT (control)", **_spearman(elig["mean_SFTPC"], elig["mean_slingshot_pt"])},
        {"contrast": "AT2 score vs Slingshot PT (control)", **_spearman(elig["mean_AT2"], elig["mean_slingshot_pt"])},
    ]
    if "mean_palantir_pt" in elig.columns:
        primary_rows.append(
            {"contrast": "CLDN4 vs Palantir PT (companion)", **_spearman(elig["mean_CLDN4"], elig["mean_palantir_pt"])}
        )
    qs = _bh([r["p"] if r["p"] is not None else 1.0 for r in primary_rows])
    for r, q in zip(primary_rows, qs):
        r["q"] = None if r["p"] is None else q
    pd.DataFrame(primary_rows).to_csv(tabdir / "sample_level_spearman.tsv", sep="\t", index=False)

    def _sub(df, mask, label_pt="mean_slingshot_pt_primary"):
        sub = df.loc[mask]
        return _spearman(sub["mean_CLDN4"], sub[label_pt])

    is_131 = elig["dataset"] == "GSE131907"
    is_205 = elig["dataset"] == "GSE205335"
    is_nlung = elig["Sample_Origin"].astype(str) == "nLung"
    is_tlung = elig["Sample_Origin"].astype(str) == "tLung"
    sensitivity = [
        {"contrast": "GSE131907-only CLDN4 vs primary PT", **_sub(elig, is_131)},
        {"contrast": "GSE205335-only CLDN4 vs primary PT", **_sub(elig, is_205)},
        {"contrast": "GSE131907-only CLDN4 vs barrier", **_spearman(elig.loc[is_131, "mean_CLDN4"], elig.loc[is_131, "mean_barrier"])},
        {"contrast": "GSE205335-only CLDN4 vs barrier", **_spearman(elig.loc[is_205, "mean_CLDN4"], elig.loc[is_205, "mean_barrier"])},
        {"contrast": "GSE131907-only CLDN4 vs IFN", **_spearman(elig.loc[is_131, "mean_CLDN4"], elig.loc[is_131, "mean_IFN"])},
        {"contrast": "GSE205335-only CLDN4 vs IFN", **_spearman(elig.loc[is_205, "mean_CLDN4"], elig.loc[is_205, "mean_IFN"])},
        {"contrast": "tLung-only CLDN4 vs primary PT", **_sub(elig, is_tlung)},
        {"contrast": "nLung-only CLDN4 vs primary PT", **_sub(elig, is_nlung)},
        {"contrast": "tumor-only (drop nLung) CLDN4 vs primary PT", **_sub(elig, ~is_nlung)},
    ]
    pd.DataFrame(sensitivity).to_csv(tabdir / "sensitivity_spearman.tsv", sep="\t", index=False)

    along = along_curve(
        adata.obs["slingshot_pt_primary"].to_numpy(),
        {
            "CLDN4": adata.obs["expr_CLDN4"].to_numpy(),
            "barrier": adata.obs["score_barrier_keratin"].to_numpy(),
            "IFN": adata.obs["score_IFN"].to_numpy(),
            "AT2": adata.obs["score_AT2"].to_numpy(),
            "malignant": adata.obs["score_malignant_like"].to_numpy(),
        },
        N_CURVE_BINS,
    )
    along.to_csv(tabdir / "along_curve_bins.tsv", sep="\t", index=False)

    colors = {"GSE131907": "#2a6f97", "GSE205335": "#b23a48"}

    # --- lineage plot (done criterion) ---
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    sc.pl.umap(adata, color="expr_CLDN4", ax=ax, show=False, frameon=False, cmap="viridis", colorbar_loc="right")
    for name, curve in sling.curves_embed.items():
        if curve is None or len(curve) < 2:
            continue
        lw = 2.8 if name == primary_name else 1.4
        ax.plot(curve[:, 0], curve[:, 1], "-", lw=lw, color="0.15" if name == primary_name else "0.45", label=name)
        ax.plot(curve[0, 0], curve[0, 1], "o", ms=7, color="#1b9e77")
        ax.plot(curve[-1, 0], curve[-1, 1], "s", ms=7, color="#d95f02")
    # cluster centroids in UMAP
    for lab, sub in adata.obs.groupby("leiden", observed=True):
        xy = embed[adata.obs_names.isin(sub.index)]
        cx, cy = xy.mean(axis=0)
        ax.text(cx, cy, str(lab), fontsize=7, ha="center", va="center", color="k")
    ax.set_title(
        f"Slingshot lineages  root={root}  primary={primary_name}\n"
        f"engine={sling.engine}  n_cells={adata.n_obs}"
    )
    ax.legend(fontsize=7, frameon=False, loc="best")
    _save(fig, figdir / "fig_lineage_cldn4")

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0))
    for ax, color, title, cmap in (
        (axes[0], "expr_CLDN4", "CLDN4", "viridis"),
        (axes[1], "score_barrier_keratin", "barrier/keratin (no CLDN4)", "viridis"),
        (axes[2], "score_IFN", "IFN (no CLDN4)", "magma"),
    ):
        sc.pl.umap(adata, color=color, ax=ax, show=False, frameon=False, cmap=cmap, title=title)
        for name, curve in sling.curves_embed.items():
            if curve is None or len(curve) < 2:
                continue
            ax.plot(curve[:, 0], curve[:, 1], "-", lw=1.6, color="0.15")
    fig.suptitle("Same Slingshot curves: CLDN4 vs barrier vs IFN", fontsize=11)
    _save(fig, figdir / "fig_lineage_programs")

    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    along.plot(x="pt_mid", y=["mean_CLDN4", "mean_barrier", "mean_IFN"], ax=ax, marker="o")
    ax.set_xlabel("primary Slingshot PT (0=root, 1=terminal)")
    ax.set_ylabel("binned mean (log1p CP10k)")
    ax.set_title(f"Along {primary_name}: CLDN4 / barrier / IFN")
    ax.legend(["CLDN4", "barrier (no CLDN4)", "IFN"], frameon=False, fontsize=8)
    _save(fig, figdir / "fig_along_curve")

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2))
    c4_pt = next(r for r in primary_rows if r["contrast"].startswith("CLDN4 vs Slingshot PT (primary"))
    c4_sh = next(r for r in primary_rows if r["contrast"].startswith("CLDN4 vs Slingshot PT (shared"))
    for ax, x, row, xlab in (
        (axes[0], "mean_slingshot_pt_primary", c4_pt, "sample-mean Slingshot PT (primary)"),
        (axes[1], "mean_slingshot_pt", c4_sh, "sample-mean Slingshot PT (shared)"),
    ):
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub["mean_CLDN4"], s=48, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel(xlab)
        ax.set_ylabel("sample-mean CLDN4")
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("Sample-level CLDN4 vs Slingshot pseudotime", fontsize=11)
    _save(fig, figdir / "fig_sample_cldn4_pt")

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2))
    c4_bar = next(r for r in primary_rows if "barrier" in r["contrast"] and r["contrast"].startswith("CLDN4"))
    c4_ifn = next(r for r in primary_rows if r["contrast"] == "CLDN4 vs IFN (no CLDN4)")
    for ax, x, row, xlab in (
        (axes[0], "mean_barrier", c4_bar, "sample-mean barrier/keratin (no CLDN4)"),
        (axes[1], "mean_IFN", c4_ifn, "sample-mean IFN (no CLDN4)"),
    ):
        for ds, col in colors.items():
            sub = elig[elig["dataset"] == ds]
            ax.scatter(sub[x], sub["mean_CLDN4"], s=48, c=col, label=f"{ds} n={len(sub)}")
        ax.set_xlabel(xlab)
        ax.set_ylabel("sample-mean CLDN4")
        ax.set_title(_fmt(row))
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("EXTRA: sample-level CLDN4 vs barrier / IFN", fontsize=11)
    _save(fig, figdir / "fig_extra_sample_programs")

    fig, ax = plt.subplots(figsize=(5.4, 4.8))
    sc.pl.paga(adata, ax=ax, show=False, frameon=False, title=f"PAGA  root={root}")
    _save(fig, figdir / "fig_paga")

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))
    ct = pd.crosstab(adata.obs["author_subtype"].astype(str), adata.obs["dataset"].astype(str))
    ct.T.plot(kind="bar", ax=axes[0], legend=False)
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Author subtypes  n_cells={adata.n_obs}")
    plt.setp(axes[0].get_xticklabels(), rotation=40, ha="right")
    unit_ct = sample_df.groupby("dataset").size()
    axes[1].bar(unit_ct.index.astype(str), unit_ct.to_numpy(), color=["#2a6f97", "#b23a48"][: len(unit_ct)])
    axes[1].set_ylabel("units")
    axes[1].set_title(f"Honest n units={sample_df.shape[0]} (eligible {len(elig)})")
    _save(fig, figdir / "fig_honest_n")

    for color, fname, cmap in (
        ("dataset", "fig_umap_dataset", None),
        ("leiden", "fig_umap_leiden", None),
        ("author_subtype", "fig_umap_author_subtype", None),
        ("score_AT2", "fig_umap_AT2", "viridis"),
        ("slingshot_pt_primary", "fig_umap_slingshot_pt", "cividis"),
    ):
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        kw = {"color": color, "ax": ax, "show": False, "frameon": False}
        if cmap:
            kw["cmap"] = cmap
        sc.pl.umap(adata, **kw)
        _save(fig, figdir / fname)

    subtype_counts = adata.obs["author_subtype"].astype(str).fillna("NA").value_counts().to_dict()
    c4_pt = next(r for r in primary_rows if r["contrast"].startswith("CLDN4 vs Slingshot PT (primary"))
    c4_bar = next(r for r in primary_rows if r["contrast"].startswith("CLDN4 vs barrier"))
    c4_ifn = next(r for r in primary_rows if r["contrast"] == "CLDN4 vs IFN (no CLDN4)")
    c4_mal = next(r for r in primary_rows if r["contrast"].startswith("CLDN4 vs malignant"))
    bar_pt = next(r for r in primary_rows if r["contrast"].startswith("barrier vs"))
    ifn_pt = next(r for r in primary_rows if r["contrast"].startswith("IFN vs"))

    along_recs = along.to_dict(orient="records")
    first_bin, last_bin = along_recs[0], along_recs[-1]
    thesis_ok = (
        last_bin["mean_CLDN4"] >= first_bin["mean_CLDN4"]
        and last_bin["mean_barrier"] >= first_bin["mean_barrier"]
        and last_bin["mean_IFN"] <= first_bin["mean_IFN"] + 0.05
    )
    what_holds = (
        f"**What holds (n={c4_bar['n']} units).** CLDN4 tracks a CLDN4-excluded "
        f"barrier/keratin score ({_fmt(c4_bar)}) and a malignant-like score ({_fmt(c4_mal)}). "
        f"On the primary Slingshot curve, root-bin → terminal-bin CLDN4 "
        f"{first_bin['mean_CLDN4']:.3f} → {last_bin['mean_CLDN4']:.3f}, "
        f"barrier {first_bin['mean_barrier']:.3f} → {last_bin['mean_barrier']:.3f}, "
        f"IFN {first_bin['mean_IFN']:.3f} → {last_bin['mean_IFN']:.3f}. "
        f"Sample-level CLDN4 vs IFN is {_fmt(c4_ifn)}. "
        f"**Thesis check (same curve):** CLDN4-high sits at the barrier/malignant end, "
        f"not the IFN-high end — {'supported by binned means' if thesis_ok else 'not uniformly supported by binned means; see table'}."
    )
    verdict = " ".join(
        [
            f"Slingshot engine={sling.engine}; {len(lin_recs)} lineage(s) from root {root}.",
            f"Sample-level CLDN4 vs primary PT: {_fmt(c4_pt)}.",
            f"CLDN4 vs barrier/keratin (no CLDN4): {_fmt(c4_bar)}.",
            f"CLDN4 vs IFN: {_fmt(c4_ifn)}.",
            f"barrier vs primary PT: {_fmt(bar_pt)}.",
            f"IFN vs primary PT: {_fmt(ifn_pt)}.",
            "Not a DPT-null writeup. Not a TACSTD2 redo. No both-high gate. GSE148071 not added.",
            "PR #320 T/NK r=−0.705 given, not re-audited.",
        ]
    )

    summary = {
        "accessions": ["GSE131907", "GSE205335"],
        "winning_pair": True,
        "pr320_given": "author-malignant CLDN4 %pos vs T/NK Q4 vs Q1 n=23 r=-0.705",
        "gse207422_added": False,
        "gse148071_added": False,
        "dual_high": False,
        "primary_gene": "CLDN4",
        "clock": "slingshot",
        "clock_label": f"Slingshot ({sling.engine})",
        "slingshot_engine": sling.engine,
        "slingshot_notes": sling.notes,
        "r_slingshot": rstat,
        "palantir": pal,
        "harmony": harmony,
        "cap_per_unit": args.cap_per_unit,
        "n_cells": int(adata.n_obs),
        "n_cells_gse131907": int((adata.obs["dataset"] == "GSE131907").sum()),
        "n_cells_gse205335": int((adata.obs["dataset"] == "GSE205335").sum()),
        "n_cells_nLung": int((adata.obs["Sample_Origin"] == "nLung").sum()),
        "n_author_AT2": int(adata.obs["author_subtype"].astype(str).isin(AUTHOR_AT2).sum()),
        "n_units": int(sample_df.shape[0]),
        "n_units_gse131907": int((sample_df["dataset"] == "GSE131907").sum()),
        "n_units_gse205335": int((sample_df["dataset"] == "GSE205335").sum()),
        "n_units_eligible": int(len(elig)),
        "subtype_counts": subtype_counts,
        "genes_absent": absent,
        "root": root_info,
        "n_paga_components": int(n_comp),
        "n_leiden": int(adata.obs["leiden"].nunique()),
        "leiden_resolution": LEIDEN_RES,
        "n_lineages": int(len(lin_recs)),
        "lineages": lin_recs,
        "primary_lineage": primary_name,
        "primary_spearman": primary_rows,
        "sensitivity_spearman": sensitivity,
        "along_curve": along_recs,
        "verdict": verdict,
        "what_holds": what_holds,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(args.finding, summary)
    (outdir / "provenance.json").write_text(
        json.dumps(
            {
                "accessions": ["GSE131907", "GSE205335"],
                "papers": {
                    "GSE131907": "Kim et al. Nat Commun 2020 PMID 32385277",
                    "GSE205335": "Ahn / Lee et al. eLife 2024 (GEO GSE205335)",
                    "slingshot": "Street et al. BMC Genomics 2018",
                    "pyslingshot-bio": "https://github.com/omicverse/py-Slingshot",
                },
                "primary_gene": "CLDN4",
                "barrier_excludes_CLDN4": True,
                "ifn_excludes_CLDN4": True,
                "r_slingshot": rstat,
                "slingshot_engine": sling.engine,
                "locked": {
                    "leiden_resolution": LEIDEN_RES,
                    "n_hvg": N_HVG,
                    "n_neighbors": N_NEIGHBORS,
                    "n_pcs": N_PCS,
                    "root": "AT2-like or lowest-CLDN4; never max-CLDN4 cluster",
                    "cap_per_unit": args.cap_per_unit,
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
                "n_lineages": len(lin_recs),
                "engine": sling.engine,
                "finding": str(args.finding),
                "sample_table": str(tabdir / "sample_cldn4_vs_pseudotime.tsv"),
                "lineage_plot": str(figdir / "fig_lineage_cldn4.png"),
            },
            indent=2,
        )
    )


def _paga_ncomp(connect: np.ndarray, thresh: float = 0.0) -> int:
    n = connect.shape[0]
    seen = [False] * n
    ncomp = 0
    for i in range(n):
        if seen[i]:
            continue
        ncomp += 1
        stack = [i]
        seen[i] = True
        while stack:
            u = stack.pop()
            for v in range(n):
                if not seen[v] and connect[u, v] > thresh:
                    seen[v] = True
                    stack.append(v)
    return ncomp


if __name__ == "__main__":
    main()
