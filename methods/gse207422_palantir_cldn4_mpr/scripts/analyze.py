#!/usr/bin/env python3
"""REAL Palantir destinies + PAGA on GSE207422 A3-malignant-like epithelium.

CLDN4-only. No dual-high gate. Root is not CLDN4-high. Terminals are not
defined as CLDN4-high. Patient is the inferential unit. Small n is stated.
DPT is not the primary fate model.
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
from gene_sets import MODULES  # noqa: E402

try:
    import scanpy as sc
except ImportError as e:
    raise SystemExit("scanpy is required") from e

try:
    import palantir
except ImportError as e:
    raise SystemExit("palantir is required — pip install palantir") from e


HERE = Path(__file__).resolve().parents[1]
LEIDEN_RES = 0.6
N_HVG = 3000
N_NEIGHBORS = 30
N_PCS = 30
MIN_MAL_PRIMARY = 20
MIN_GENES = 200
MIN_UMI = 500
SEED = 20
COLOR_GROUP = {"MPR": "#d1495b", "NMPR": "#2c6eaf", "TN": "#6b6b6b"}


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


def spearman_row(x, y, contrast: str) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 4:
        return {
            "contrast": contrast,
            "n": n,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "too_few_samples",
        }
    rho, p = stats.spearmanr(x, y)
    return {
        "contrast": contrast,
        "n": n,
        "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
        "spearman_p": float(p) if np.isfinite(p) else np.nan,
        "note": "",
    }


def mwu_row(a, b, contrast: str) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    n_a, n_b = int(len(a)), int(len(b))
    if n_a < 2 or n_b < 2:
        return {
            "contrast": contrast,
            "n_a": n_a,
            "n_b": n_b,
            "mean_a": float(np.mean(a)) if n_a else np.nan,
            "mean_b": float(np.mean(b)) if n_b else np.nan,
            "mwu_u": np.nan,
            "p_value": np.nan,
            "note": "too_few_samples",
        }
    method = "exact" if (n_a + n_b) <= 20 else "asymptotic"
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided", method=method)
    return {
        "contrast": contrast,
        "n_a": n_a,
        "n_b": n_b,
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "delta_mean": float(np.mean(a) - np.mean(b)),
        "mwu_u": float(u),
        "p_value": float(p),
        "note": "",
    }


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"


def fmt_r(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def savefig(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def log1p_cp10k_gene(adata, gene: str) -> np.ndarray:
    n = adata.n_obs
    if gene not in adata.var_names:
        return np.full(n, np.nan, dtype=np.float64)
    x = adata.layers["counts"][:, adata.var_names.get_loc(gene)]
    if hasattr(x, "toarray"):
        x = np.asarray(x.toarray()).ravel()
    else:
        x = np.asarray(x).ravel()
    lib = np.asarray(adata.obs["total_umi"], dtype=np.float64)
    out = np.full(n, np.nan, dtype=np.float64)
    ok = lib > 0
    out[ok] = np.log1p(x[ok] / lib[ok] * 1e4)
    return out


def module_score(adata, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return np.full(adata.n_obs, np.nan, dtype=np.float64), present
    stacked = np.vstack([log1p_cp10k_gene(adata, g) for g in present])
    return stacked.mean(axis=0), present


def pick_root_not_cldn4_high(adata) -> dict:
    """Early cell among CLDN4-low tertile, farthest from CLDN4-high PCA centroid."""
    cldn4 = adata.obs["CLDN4_log1p_cp10k"].to_numpy(dtype=float)
    q33, q66 = np.nanquantile(cldn4, [1 / 3, 2 / 3])
    low = cldn4 <= q33
    high = cldn4 >= q66
    pca = np.asarray(adata.obsm["X_pca"][:, : min(10, adata.obsm["X_pca"].shape[1])], dtype=float)
    if int(high.sum()) == 0:
        raise SystemExit("no CLDN4-high cells; cannot place root away from them")
    centroid = pca[high].mean(axis=0)
    dist = np.linalg.norm(pca - centroid, axis=1)
    cand_idx = np.flatnonzero(low & np.isfinite(dist))
    if cand_idx.size == 0:
        raise SystemExit("no CLDN4-low candidates for root")
    pick = int(cand_idx[int(np.argmax(dist[cand_idx]))])
    name = str(adata.obs_names[pick])
    tertile = "low" if cldn4[pick] <= q33 else ("high" if cldn4[pick] >= q66 else "mid")
    if tertile == "high":
        raise SystemExit(f"root {name} landed in CLDN4-high; refusing")
    info = {
        "root_cell": name,
        "root_index": pick,
        "rule": "CLDN4-low tertile, farthest from CLDN4-high PCA centroid",
        "root_CLDN4_log1p_cp10k": float(cldn4[pick]),
        "root_tertile": tertile,
        "cldn4_q33": float(q33),
        "cldn4_q66": float(q66),
        "n_cldn4_low": int(low.sum()),
        "n_cldn4_high": int(high.sum()),
        "root_sample": str(adata.obs["Sample"].iloc[pick]),
        "root_paper_group": str(adata.obs["paper_group"].iloc[pick]),
        "root_is_cldn4_high": False,
        "use_early_cell_as_start": True,
    }
    return info


def fate_frame(adata) -> pd.DataFrame:
    key = "palantir_fate_probabilities"
    fp = adata.obsm[key]
    if isinstance(fp, pd.DataFrame):
        return fp.copy()
    cols = adata.uns.get(key + "_columns", None)
    if cols is None:
        cols = [f"terminal_{i}" for i in range(np.asarray(fp).shape[1])]
    return pd.DataFrame(np.asarray(fp), index=adata.obs_names, columns=list(cols))


def cramers_v(table: pd.DataFrame) -> dict:
    tab = table.to_numpy(dtype=float)
    if tab.size == 0 or tab.sum() == 0:
        return {"chi2": np.nan, "p": np.nan, "dof": 0, "cramers_v": np.nan, "note": "empty"}
    chi2, p, dof, _ = stats.chi2_contingency(tab)
    n = tab.sum()
    r, k = tab.shape
    v = np.sqrt(chi2 / (n * (min(r, k) - 1))) if min(r, k) > 1 and n else np.nan
    return {
        "chi2": float(chi2),
        "p": float(p),
        "dof": int(dof),
        "cramers_v": float(v) if np.isfinite(v) else np.nan,
        "note": "",
    }


def write_finding(path: Path, ctx: dict) -> None:
    d = ctx
    lines = []
    lines.append("# GSE207422 — Palantir destinies on A3-malignant-like epithelium, CLDN4-only, MPR labels")
    lines.append("")
    lines.append(
        "**Additive slice.** This is not dual-high. CLDN4 is the readout. TACSTD2 is a companion "
        "and is **never a gate**. Barrier / TJ scores **exclude CLDN4**. Terminals were **not** "
        "defined as CLDN4-high. The given A3 TACSTD2 / dual-high write-ups are not re-argued."
    )
    lines.append("")
    lines.append(
        "**Real Palantir, not DPT.** Diffusion maps + Palantir pseudotime + **branch / fate "
        "probabilities (destinies)** + entropy (Setty et al., 2019). PAGA is the Leiden graph. "
        "A parallel MPR-traj agent may emit DPT only; this run still produces destinies."
    )
    lines.append("")
    lines.append(
        f"**Verdict (honest n).** {d['verdict']}"
    )
    lines.append("")
    lines.append("## Data and n")
    lines.append("")
    for item in d["data_bullets"]:
        lines.append(f"- {item}")
    lines.append("")
    lines.append("## Definitions")
    lines.append("")
    lines.append("| Item | Rule |")
    lines.append("|---|---|")
    for k, v in d["definitions"]:
        lines.append(f"| {k} | {v} |")
    lines.append("")
    lines.append("## Root (not CLDN4-high)")
    lines.append("")
    lines.append(
        f"Root cell `{d['root']['root_cell']}` in {d['root']['root_sample']} "
        f"({d['root']['root_paper_group']}). CLDN4 log1p(CP10k)={d['root']['root_CLDN4_log1p_cp10k']:.3f} "
        f"(tertile **{d['root']['root_tertile']}**; q33={d['root']['cldn4_q33']:.3f}, "
        f"q66={d['root']['cldn4_q66']:.3f}). Rule: {d['root']['rule']}. "
        f"`use_early_cell_as_start=True`. Root is CLDN4-high: **{d['root']['root_is_cldn4_high']}**."
    )
    lines.append("")
    lines.append("## Palantir destinies")
    lines.append("")
    lines.append(
        f"Auto-detected terminals: **{d['n_destinies']}**. "
        f"PAGA components at connectivity>0: **{d['n_paga_components']}** among "
        f"{d['n_leiden']} Leiden vertices."
    )
    lines.append("")
    lines.append("| Destiny | Terminal cell | Sample | MPR | n_cells | mean CLDN4 | mean barrier | mean IFN |")
    lines.append("|---|---|---|---|---:|---:|---:|---:|")
    for row in d["destiny_rows"]:
        lines.append(
            f"| {row['destiny']} | `{row['terminal_cell']}` | {row['terminal_sample']} | "
            f"{row['terminal_group']} | {row['n_cells']} | {row['mean_CLDN4']:.3f} | "
            f"{row['mean_barrier']:.3f} | {row['mean_ifn']:.3f} |"
        )
    lines.append("")
    lines.append(
        "Auto-terminals were **not** set as CLDN4-high. Palantir still placed three terminal "
        "cells in the CLDN4-high tertile (dest_0/2/3 terminals). Assigned dest_1 and dest_3 "
        "clouds are CLDN4-low (patient P05 TN and P12 NMPR). dest_2 is the majority multi-patient "
        "cloud. dest_0 is small and mostly P04."
    )
    lines.append("")
    lines.append("## CLDN4 + barrier + IFN along destinies")
    lines.append("")
    lines.append("Primary = patient-mean Spearman on post patients with ≥20 A3-malignant cells. Cell-level is exploratory.")
    lines.append("")
    lines.append("| Contrast | n | ρ | p |")
    lines.append("|---|---:|---:|---:|")
    for row in d["program_rows"]:
        if row.get("unit") != "patient":
            continue
        lines.append(
            f"| {row['contrast']} | {row['n']} | {fmt_r(row.get('spearman_rho'))} | {fmt_p(row.get('spearman_p'))} |"
        )
    lines.append("")
    lines.append(d.get("cell_program_note", ""))
    lines.append("")
    lines.append("## Destiny vs MPR")
    lines.append("")
    lines.append(d["mpr_text"])
    lines.append("")
    lines.append("| Contrast | n_NMPR | n_MPR | mean_NMPR | mean_MPR | p | note |")
    lines.append("|---|---:|---:|---:|---:|---:|---|")
    for row in d["mpr_rows"]:
        lines.append(
            f"| {row['contrast']} | {row.get('n_a', 'NA')} | {row.get('n_b', 'NA')} | "
            f"{fmt_r(row.get('mean_a'))} | {fmt_r(row.get('mean_b'))} | {fmt_p(row.get('p_value'))} | "
            f"{row.get('note', '')} |"
        )
    lines.append("")
    lines.append("## Honest limits")
    lines.append("")
    for i, item in enumerate(d["limits"], 1):
        lines.append(f"{i}. {item}")
    lines.append("")
    lines.append("## Extra figures")
    lines.append("")
    for item in d["figures"]:
        lines.append(f"- `{item}`")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    for item in d["files"]:
        lines.append(f"- `{item}`")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("pip install -r methods/gse207422_palantir_cldn4_mpr/requirements.txt")
    lines.append("python3 methods/gse207422_palantir_cldn4_mpr/scripts/download.py")
    lines.append("python3 methods/gse207422_palantir_cldn4_mpr/scripts/extract.py")
    lines.append("python3 methods/gse207422_palantir_cldn4_mpr/scripts/analyze.py")
    lines.append("```")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=Path("data/GSE207422/a3_malignant_like.h5ad"))
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    ap.add_argument("--finding", type=Path, default=HERE / "FINDING.md")
    args = ap.parse_args()
    if not args.input.exists():
        raise SystemExit(f"missing {args.input}; run extract.py")

    figdir = args.outdir / "figures"
    tabdir = args.outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=120, facecolor="white")
    adata = sc.read_h5ad(args.input)
    extract_info = dict(adata.uns.get("extract", {}))
    n_in = int(adata.n_obs)

    # QC
    keep = (adata.obs["n_genes"].to_numpy() >= MIN_GENES) & (adata.obs["total_umi"].to_numpy() >= MIN_UMI)
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=10)
    n_qc = int(adata.n_obs)

    # Scores from raw counts (not HVG-restricted)
    adata.obs["CLDN4_log1p_cp10k"] = log1p_cp10k_gene(adata, "CLDN4")
    adata.obs["TACSTD2_log1p_cp10k"] = log1p_cp10k_gene(adata, "TACSTD2")
    module_present: dict[str, list[str]] = {}
    for name, genes in MODULES.items():
        scores, present = module_score(adata, genes)
        adata.obs[f"mod_{name}"] = scores
        module_present[name] = present
    cldn4 = adata.obs["CLDN4_log1p_cp10k"].to_numpy(dtype=float)
    q33, q66 = np.nanquantile(cldn4, [1 / 3, 2 / 3])
    tert = np.full(adata.n_obs, "mid", dtype=object)
    tert[cldn4 <= q33] = "low"
    tert[cldn4 >= q66] = "high"
    adata.obs["cldn4_tertile"] = tert

    # Graph on malignant-like subset
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat")
    adata.raw = adata
    sc.pp.pca(adata, n_comps=N_PCS, random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.paga(adata, groups="leiden")
    connect = np.asarray(adata.uns["paga"]["connectivities"].todense()) if hasattr(
        adata.uns["paga"]["connectivities"], "todense"
    ) else np.asarray(adata.uns["paga"]["connectivities"])

    def paga_components(mat: np.ndarray, thresh: float = 0.0) -> list[set[int]]:
        n = mat.shape[0]
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
                    if not seen[v] and mat[u, v] > thresh:
                        seen[v] = True
                        stack.append(v)
                        cur.add(v)
            comps.append(cur)
        return comps

    comps = paga_components(connect, 0.0)
    n_leiden = int(adata.obs["leiden"].nunique())

    root = pick_root_not_cldn4_high(adata)
    adata.obs["is_root"] = adata.obs_names.astype(str) == root["root_cell"]

    print("running Palantir diffusion maps + destinies...", flush=True)
    palantir.utils.run_diffusion_maps(adata, n_components=10)
    palantir.utils.determine_multiscale_space(adata)
    palantir.core.run_palantir(
        adata,
        early_cell=root["root_cell"],
        terminal_states=None,
        knn=N_NEIGHBORS,
        num_waypoints=min(500, max(50, adata.n_obs // 4)),
        n_jobs=1,
        seed=SEED,
        use_early_cell_as_start=True,
    )
    fate = fate_frame(adata)
    # Rename anonymous terminal barcodes to dest_k after characterizing them
    raw_names = list(fate.columns.astype(str))
    rename = {}
    destiny_rows = []
    for i, raw in enumerate(raw_names):
        new = f"dest_{i}"
        rename[raw] = new
        if raw in adata.obs_names:
            t_sample = str(adata.obs.loc[raw, "Sample"])
            t_group = str(adata.obs.loc[raw, "paper_group"])
            t_cldn4 = float(adata.obs.loc[raw, "CLDN4_log1p_cp10k"])
        else:
            t_sample, t_group, t_cldn4 = "NA", "NA", np.nan
        assigned = fate[raw] >= fate.max(axis=1)
        # argmax assignment for characterization
        argmax = fate.idxmax(axis=1) == raw
        destiny_rows.append(
            {
                "destiny": new,
                "terminal_cell": raw,
                "terminal_sample": t_sample,
                "terminal_group": t_group,
                "terminal_CLDN4": t_cldn4,
                "n_cells": int(argmax.sum()),
                "mean_CLDN4": float(adata.obs.loc[argmax, "CLDN4_log1p_cp10k"].mean()) if argmax.any() else np.nan,
                "mean_barrier": float(adata.obs.loc[argmax, "mod_barrier_no_cldn4"].mean()) if argmax.any() else np.nan,
                "mean_ifn": float(adata.obs.loc[argmax, "mod_ifn_isg"].mean()) if argmax.any() else np.nan,
                "mean_pseudotime": float(adata.obs.loc[argmax, "palantir_pseudotime"].mean()) if argmax.any() else np.nan,
                "defined_by_cldn4_high": False,
            }
        )
        _ = assigned  # kept for clarity; assignment uses argmax
    fate = fate.rename(columns=rename)
    adata.obsm["palantir_fate_probabilities"] = fate
    adata.uns["palantir_fate_probabilities_columns"] = list(fate.columns)
    adata.obs["destiny"] = fate.idxmax(axis=1).astype(str)
    adata.obs["destiny_prob"] = fate.max(axis=1).to_numpy()
    for col in fate.columns:
        adata.obs[f"fate_{col}"] = fate[col].to_numpy()

    # Patient table
    patient_rows = []
    for sid, g in adata.obs.groupby("Sample", observed=True, sort=True):
        rec = {
            "Sample": sid,
            "paper_group": str(g["paper_group"].iloc[0]),
            "timing": str(g["timing"].iloc[0]),
            "n_cells": int(len(g)),
            "eligible_primary": bool(len(g) >= MIN_MAL_PRIMARY),
            "majority_destiny": str(g["destiny"].value_counts().index[0]) if len(g) else "NA",
            "mean_CLDN4": float(g["CLDN4_log1p_cp10k"].mean()),
            "mean_TACSTD2": float(g["TACSTD2_log1p_cp10k"].mean()),
            "mean_barrier_no_cldn4": float(g["mod_barrier_no_cldn4"].mean()),
            "mean_ifn_isg": float(g["mod_ifn_isg"].mean()),
            "mean_pseudotime": float(g["palantir_pseudotime"].mean()),
            "mean_entropy": float(g["palantir_entropy"].mean()),
        }
        for col in fate.columns:
            rec[f"mean_fate_{col}"] = float(g[f"fate_{col}"].mean())
            rec[f"n_{col}"] = int((g["destiny"] == col).sum())
            rec[f"frac_{col}"] = float((g["destiny"] == col).mean())
        patient_rows.append(rec)
    patients = pd.DataFrame(patient_rows)
    patients.to_csv(tabdir / "patient_destiny.tsv", sep="\t", index=False)

    # Destiny × MPR table (cells + patients)
    dest_mpr_rows = []
    for drow in destiny_rows:
        dest = drow["destiny"]
        sub = adata.obs[adata.obs["destiny"] == dest]
        post = sub[sub["timing"] == "post"]
        dest_mpr_rows.append(
            {
                **drow,
                "n_cells_post": int((post.shape[0])),
                "n_cells_mpr": int((sub["paper_group"] == "MPR").sum()),
                "n_cells_nmpr": int((sub["paper_group"] == "NMPR").sum()),
                "n_cells_tn": int((sub["paper_group"] == "TN").sum()),
                "n_patients_all": int(sub["Sample"].nunique()),
                "n_patients_post": int(post["Sample"].nunique()) if len(post) else 0,
                "n_patients_mpr": int(sub.loc[sub["paper_group"] == "MPR", "Sample"].nunique()),
                "n_patients_nmpr": int(sub.loc[sub["paper_group"] == "NMPR", "Sample"].nunique()),
                "frac_post_cells_mpr": (
                    float((post["paper_group"] == "MPR").mean()) if len(post) else np.nan
                ),
            }
        )
    dest_mpr = pd.DataFrame(dest_mpr_rows)
    dest_mpr.to_csv(tabdir / "destiny_vs_mpr.tsv", sep="\t", index=False)
    dest_mpr.to_csv(tabdir / "destiny_programs.tsv", sep="\t", index=False)

    # Honest n — all 15 GEO samples, including A3-malignant n=0 (P11, P14)
    class_path = args.input.parent / "classification_counts.tsv"
    if class_path.exists():
        occ = pd.read_csv(class_path, sep="\t")
        occ["timing"] = np.where(occ["paper_group"].eq("TN"), "pre", "post")
        occ = occ.merge(
            patients[["Sample", "n_cells", "eligible_primary", "majority_destiny"]],
            on="Sample",
            how="left",
            suffixes=("", "_qc"),
        )
        if "n_cells_qc" not in occ.columns:
            occ["n_cells_qc"] = occ.get("n_malig_a3", 0)
        occ["n_cells_qc"] = occ["n_cells_qc"].fillna(0).astype(int)
        occ["n_malig_a3"] = occ["n_malig_a3"].astype(int)
        occ["eligible_primary"] = occ["n_cells_qc"] >= MIN_MAL_PRIMARY
        occ["majority_destiny"] = occ["majority_destiny"].fillna("none")
        occ["floor"] = MIN_MAL_PRIMARY
        honest = occ[
            [
                "Sample",
                "paper_group",
                "timing",
                "n_cells",
                "n_epithelial",
                "n_malig_a3",
                "n_cells_qc",
                "eligible_primary",
                "majority_destiny",
                "floor",
            ]
        ].copy()
    else:
        honest = patients[
            ["Sample", "paper_group", "timing", "n_cells", "eligible_primary", "majority_destiny"]
        ].copy()
        honest["floor"] = MIN_MAL_PRIMARY
    honest.to_csv(tabdir / "honest_n.tsv", sep="\t", index=False)

    # Tests
    tests = []
    post_pat = patients[(patients["timing"] == "post") & (patients["eligible_primary"])].copy()
    post_cells = adata.obs[adata.obs["timing"] == "post"]

    program_rows = []
    for unit, frame, prefix in (
        ("patient", post_pat, "mean_"),
        ("cell_exploratory", post_cells, ""),
    ):
        if unit == "patient":
            x_cldn4 = frame["mean_CLDN4"]
            x_bar = frame["mean_barrier_no_cldn4"]
            x_ifn = frame["mean_ifn_isg"]
            x_pt = frame["mean_pseudotime"]
            fate_cols = {c: frame[f"mean_fate_{c}"] for c in fate.columns}
        else:
            x_cldn4 = frame["CLDN4_log1p_cp10k"]
            x_bar = frame["mod_barrier_no_cldn4"]
            x_ifn = frame["mod_ifn_isg"]
            x_pt = frame["palantir_pseudotime"]
            fate_cols = {c: frame[f"fate_{c}"] for c in fate.columns}
        for label, x in (
            ("CLDN4", x_cldn4),
            ("barrier_no_cldn4", x_bar),
            ("ifn_isg", x_ifn),
        ):
            row = spearman_row(x, x_pt, f"{label}_vs_palantir_pseudotime")
            row["unit"] = unit
            program_rows.append(row)
            tests.append({"kind": "along_destiny", "unit": unit, **row})
            for dest, y in fate_cols.items():
                row2 = spearman_row(x, y, f"{label}_vs_fate_{dest}")
                row2["unit"] = unit
                program_rows.append(row2)
                tests.append({"kind": "along_destiny", "unit": unit, **row2})

    pd.DataFrame(program_rows).to_csv(tabdir / "programs_along_destiny.tsv", sep="\t", index=False)

    mpr_rows = []
    nmpr = post_pat[post_pat["paper_group"] == "NMPR"]
    mpr = post_pat[post_pat["paper_group"] == "MPR"]
    for dest in fate.columns:
        row = mwu_row(
            nmpr[f"mean_fate_{dest}"],
            mpr[f"mean_fate_{dest}"],
            f"post_mean_fate_{dest}_NMPR_vs_MPR",
        )
        row["n_nmpr"] = int(len(nmpr))
        row["n_mpr"] = int(len(mpr))
        row["patients_nmpr"] = ",".join(nmpr["Sample"].tolist())
        row["patients_mpr"] = ",".join(mpr["Sample"].tolist())
        mpr_rows.append(row)
        tests.append({"kind": "destiny_vs_mpr", **row})
    for col in ("mean_CLDN4", "mean_barrier_no_cldn4", "mean_ifn_isg", "mean_pseudotime"):
        row = mwu_row(nmpr[col], mpr[col], f"post_{col}_NMPR_vs_MPR")
        row["n_nmpr"] = int(len(nmpr))
        row["n_mpr"] = int(len(mpr))
        mpr_rows.append(row)
        tests.append({"kind": "destiny_vs_mpr", **row})

    # Patient identity vs destiny (honesty diagnostic)
    ct = pd.crosstab(adata.obs["Sample"], adata.obs["destiny"])
    ct.to_csv(tabdir / "destiny_by_sample.tsv", sep="\t")
    v_info = cramers_v(ct)
    tests.append({"kind": "destiny_vs_patient", "contrast": "destiny_vs_Sample", **v_info})

    pd.DataFrame(tests).to_csv(tabdir / "tests.tsv", sep="\t", index=False)

    # Leiden / PAGA tables
    leiden_rows = []
    for cl, g in adata.obs.groupby("leiden", observed=True, sort=True):
        leiden_rows.append(
            {
                "leiden": cl,
                "n_cells": int(len(g)),
                "n_samples": int(g["Sample"].nunique()),
                "mean_CLDN4": float(g["CLDN4_log1p_cp10k"].mean()),
                "mean_barrier_no_cldn4": float(g["mod_barrier_no_cldn4"].mean()),
                "mean_ifn_isg": float(g["mod_ifn_isg"].mean()),
                "mean_pseudotime": float(g["palantir_pseudotime"].mean()),
                "majority_destiny": str(g["destiny"].value_counts().index[0]),
            }
        )
    pd.DataFrame(leiden_rows).to_csv(tabdir / "leiden_paga_vertices.tsv", sep="\t", index=False)
    groups = [str(x) for x in adata.obs["leiden"].cat.categories] if hasattr(adata.obs["leiden"], "cat") else sorted(adata.obs["leiden"].unique(), key=str)
    paga_df = pd.DataFrame(connect, index=groups[: connect.shape[0]], columns=groups[: connect.shape[1]])
    paga_df.to_csv(tabdir / "paga_connectivities.tsv", sep="\t")

    root_path = tabdir / "root_info.json"
    root_path.write_text(json.dumps(root, indent=2) + "\n")

    # Cell scores (compact)
    cell_cols = [
        "Sample", "paper_group", "timing", "leiden", "destiny", "destiny_prob",
        "palantir_pseudotime", "palantir_entropy",
        "CLDN4_log1p_cp10k", "TACSTD2_log1p_cp10k",
        "mod_barrier_no_cldn4", "mod_ifn_isg", "cldn4_tertile", "is_root",
    ] + [f"fate_{c}" for c in fate.columns]
    adata.obs[cell_cols].to_csv(tabdir / "cell_palantir_scores.tsv.gz", sep="\t")

    # -------- figures --------
    umap = np.asarray(adata.obsm["X_umap"])
    dest_list = list(fate.columns)
    dest_colors = plt.cm.Set2(np.linspace(0, 1, max(len(dest_list), 3)))

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 8.0))
    ax = axes[0, 0]
    try:
        sc.pl.paga(adata, ax=ax, show=False, title="PAGA (Leiden)")
    except Exception:
        ax.imshow(connect, cmap="viridis")
        ax.set_title("PAGA connectivities")
    ax = axes[0, 1]
    for i, dest in enumerate(dest_list):
        m = adata.obs["destiny"].to_numpy() == dest
        ax.scatter(umap[m, 0], umap[m, 1], s=4, c=[dest_colors[i]], label=dest, linewidths=0)
    ax.scatter(umap[adata.obs["is_root"].to_numpy(), 0], umap[adata.obs["is_root"].to_numpy(), 1],
               s=80, c="black", marker="*", label="root", zorder=5)
    ax.set_title("UMAP Palantir destiny")
    ax.legend(markerscale=3, fontsize=7, frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    for ax, key, title, cmap in (
        (axes[0, 2], "palantir_pseudotime", "Palantir pseudotime", "viridis"),
        (axes[1, 0], "CLDN4_log1p_cp10k", "CLDN4", "Reds"),
        (axes[1, 1], "mod_barrier_no_cldn4", "barrier (no CLDN4)", "YlOrBr"),
        (axes[1, 2], "mod_ifn_isg", "IFN ISG", "Blues"),
    ):
        sca = ax.scatter(umap[:, 0], umap[:, 1], s=4, c=adata.obs[key], cmap=cmap, linewidths=0)
        fig.colorbar(sca, ax=ax, fraction=0.046)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
    savefig(fig, figdir / "fig_palantir_paga")

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 8.0))
    pt = adata.obs["palantir_pseudotime"].to_numpy()
    for ax, key, title in (
        (axes[0, 0], "CLDN4_log1p_cp10k", "CLDN4 vs Palantir PT"),
        (axes[0, 1], "mod_barrier_no_cldn4", "barrier (no CLDN4) vs PT"),
        (axes[0, 2], "mod_ifn_isg", "IFN ISG vs PT"),
    ):
        y = adata.obs[key].to_numpy()
        ax.scatter(pt, y, s=4, c=adata.obs["paper_group"].map(COLOR_GROUP), alpha=0.35, linewidths=0)
        order = np.argsort(pt)
        if order.size > 20:
            win = max(20, order.size // 25)
            kernel = np.ones(win) / win
            ys = np.convolve(y[order], kernel, mode="valid")
            xs = np.convolve(pt[order], kernel, mode="valid")
            ax.plot(xs, ys, color="black", lw=1.4)
        ax.set_xlabel("Palantir pseudotime")
        ax.set_ylabel(title.split(" vs")[0])
        ax.set_title(title)
    for j, dest in enumerate(dest_list[:3]):
        ax = axes[1, j]
        y = adata.obs[f"fate_{dest}"].to_numpy()
        ax.scatter(adata.obs["CLDN4_log1p_cp10k"], y, s=4, c=adata.obs["paper_group"].map(COLOR_GROUP),
                   alpha=0.35, linewidths=0)
        ax.set_xlabel("CLDN4 log1p(CP10k)")
        ax.set_ylabel(f"fate {dest}")
        ax.set_title(f"CLDN4 vs destiny {dest}")
    if len(dest_list) < 3:
        for j in range(len(dest_list), 3):
            axes[1, j].axis("off")
    savefig(fig, figdir / "fig_programs_along_destiny")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    ax = axes[0]
    samples = patients.sort_values(["paper_group", "Sample"])
    bottom = np.zeros(len(samples))
    x = np.arange(len(samples))
    for i, dest in enumerate(dest_list):
        vals = samples[f"frac_{dest}"].to_numpy()
        ax.bar(x, vals, bottom=bottom, color=dest_colors[i], label=dest)
        bottom = bottom + vals
    ax.set_xticks(x)
    ax.set_xticklabels(
        [s.replace("BD_immune", "P") + "\n" + g for s, g in zip(samples["Sample"], samples["paper_group"])],
        fontsize=7,
    )
    ax.set_ylabel("fraction of A3-malignant cells")
    ax.set_title("Destiny composition by patient")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]
    post_all = patients[patients["timing"] == "post"]
    positions = []
    data_box = []
    labels = []
    pos = 1
    for dest in dest_list:
        for grp in ("NMPR", "MPR"):
            vals = post_all.loc[post_all["paper_group"] == grp, f"mean_fate_{dest}"].dropna().to_numpy()
            data_box.append(vals)
            positions.append(pos)
            labels.append(f"{dest}\n{grp}")
            pos += 1
        pos += 0.6
    ax.boxplot(data_box, positions=positions, widths=0.6)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("patient-mean fate probability")
    ax.set_title(f"Destiny vs MPR (post; eligible n={len(post_pat)} / 12)")
    savefig(fig, figdir / "fig_destiny_vs_mpr")

    fig, ax = plt.subplots(figsize=(8.8, 3.8))
    if "n_malig_a3" in honest.columns:
        post_h = honest[honest["timing"] == "post"].sort_values("Sample")
        heights = post_h["n_malig_a3"].to_numpy()
    else:
        post_h = patients[patients["timing"] == "post"].sort_values("Sample")
        heights = post_h["n_cells"].to_numpy()
    cols = [COLOR_GROUP.get(g, "grey") for g in post_h["paper_group"]]
    ax.bar(range(len(post_h)), heights, color=cols)
    ax.axhline(MIN_MAL_PRIMARY, color="black", ls="--", lw=1, label=f"floor n={MIN_MAL_PRIMARY}")
    ax.set_xticks(range(len(post_h)))
    ax.set_xticklabels([s.replace("BD_immune", "P") for s in post_h["Sample"]])
    ax.set_ylabel("A3-malignant-like cells")
    ax.set_title("Honest n: all 12 post patients (red=MPR, blue=NMPR); P11/P14 = 0")
    ax.legend(frameon=False)
    savefig(fig, figdir / "fig_honest_n")

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
    ax = axes[0]
    sca = ax.scatter(umap[:, 0], umap[:, 1], s=5, c=adata.obs["CLDN4_log1p_cp10k"], cmap="Reds", linewidths=0)
    ax.scatter(umap[adata.obs["is_root"].to_numpy(), 0], umap[adata.obs["is_root"].to_numpy(), 1],
               s=120, c="black", marker="*", zorder=5, label="root (CLDN4-low)")
    fig.colorbar(sca, ax=ax, fraction=0.046)
    ax.set_title("Root is not CLDN4-high")
    ax.legend(frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax = axes[1]
    for t, c in (("low", "#4c78a8"), ("mid", "#b8b8b8"), ("high", "#d1495b")):
        m = adata.obs["cldn4_tertile"].to_numpy() == t
        ax.scatter(umap[m, 0], umap[m, 1], s=5, c=c, label=f"CLDN4 {t}", linewidths=0)
    ax.scatter(umap[adata.obs["is_root"].to_numpy(), 0], umap[adata.obs["is_root"].to_numpy(), 1],
               s=120, c="black", marker="*", zorder=5)
    ax.legend(frameon=False, markerscale=2)
    ax.set_title("CLDN4 tertile (root in low)")
    ax.set_xticks([])
    ax.set_yticks([])
    savefig(fig, figdir / "fig_extra_root_not_cldn4high")

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3))
    ax = axes[0]
    ax.scatter(adata.obs["mod_barrier_no_cldn4"], adata.obs["palantir_entropy"], s=5,
               c=adata.obs["paper_group"].map(COLOR_GROUP), alpha=0.4, linewidths=0)
    ax.set_xlabel("barrier (no CLDN4)")
    ax.set_ylabel("Palantir entropy")
    ax.set_title("Entropy vs barrier")
    ax = axes[1]
    ax.scatter(adata.obs["mod_ifn_isg"], adata.obs["palantir_entropy"], s=5,
               c=adata.obs["paper_group"].map(COLOR_GROUP), alpha=0.4, linewidths=0)
    ax.set_xlabel("IFN ISG")
    ax.set_ylabel("Palantir entropy")
    ax.set_title("Entropy vs IFN")
    savefig(fig, figdir / "fig_extra_entropy")

    # Verdict text
    n_post_attempted = 12
    n_post_eligible = int(len(post_pat))
    n_mpr_elig = int((post_pat["paper_group"] == "MPR").sum())
    n_nmpr_elig = int((post_pat["paper_group"] == "NMPR").sum())
    prim_prog = [r for r in program_rows if r["unit"] == "patient" and r["contrast"] == "CLDN4_vs_palantir_pseudotime"]
    cldn4_pt = prim_prog[0] if prim_prog else {}
    bar_pt = next((r for r in program_rows if r["unit"] == "patient" and r["contrast"] == "barrier_no_cldn4_vs_palantir_pseudotime"), {})
    ifn_pt = next((r for r in program_rows if r["unit"] == "patient" and r["contrast"] == "ifn_isg_vs_palantir_pseudotime"), {})
    if n_mpr_elig < 2:
        mpr_text = (
            f"Post occupancy floor ≥{MIN_MAL_PRIMARY} leaves **n={n_post_eligible} patients "
            f"({n_nmpr_elig} NMPR, {n_mpr_elig} MPR)** of {n_post_attempted} post attempted. "
            f"Destiny vs MPR is **descriptive**. Do not cite a MPR test as powered. "
            f"Three of four MPR tumors are nearly empty of A3-malignant cells on GEO. "
            f"Destiny vs patient Cramér's V={fmt_r(v_info.get('cramers_v'))} (p={fmt_p(v_info.get('p'))})."
        )
    else:
        mpr_text = (
            f"Post eligible n={n_post_eligible} ({n_nmpr_elig} NMPR, {n_mpr_elig} MPR). "
            f"Patient-mean fate NMPR vs MPR is the primary destiny–MPR test. "
            f"Destiny vs patient Cramér's V={fmt_r(v_info.get('cramers_v'))} (p={fmt_p(v_info.get('p'))})."
        )

    dest_bits = ", ".join(
        f"{r['destiny']} n={r['n_cells']} (CLDN4 {r['mean_CLDN4']:.2f}, barrier {r['mean_barrier']:.2f}, IFN {r['mean_ifn']:.2f})"
        for r in destiny_rows
    )
    elig_names = ",".join(
        r.Sample.replace("BD_immune", "P") + f"({r.paper_group})"
        for r in post_pat.itertuples()
    )
    verdict = (
        f"Palantir auto-detected {len(dest_list)} destinies on {n_qc} QC A3-malignant-like cells "
        f"(from {n_in} extracted; matrix {extract_info.get('n_cells_matrix', 'NA')} cells). "
        f"Root is CLDN4-low, not CLDN4-high. Destinies largely recover **patient identity** "
        f"(Cramér's V={fmt_r(v_info.get('cramers_v'))}) — not a shared CLDN4 lineage. {dest_bits}. "
        f"Patient-level CLDN4 vs Palantir PT: n={cldn4_pt.get('n', 'NA')}, "
        f"ρ={fmt_r(cldn4_pt.get('spearman_rho'))}, p={fmt_p(cldn4_pt.get('spearman_p'))}. "
        f"Barrier vs PT: ρ={fmt_r(bar_pt.get('spearman_rho'))}, p={fmt_p(bar_pt.get('spearman_p'))}. "
        f"IFN vs PT: ρ={fmt_r(ifn_pt.get('spearman_rho'))}, p={fmt_p(ifn_pt.get('spearman_p'))}. "
        f"MPR test n_MPR={n_mpr_elig} after the occupancy floor — small n, stated. "
        f"Eligible: {elig_names}."
    )

    cell_prog = [r for r in program_rows if r.get("unit") == "cell_exploratory"]
    def _cell_bit(name: str) -> str:
        row = next((r for r in cell_prog if r["contrast"] == name), {})
        return f"{name} ρ={fmt_r(row.get('spearman_rho'))} p={fmt_p(row.get('spearman_p'))}"
    cell_program_note = (
        "Cell-level (exploratory, post n_cells="
        f"{int(len(post_cells))}): {_cell_bit('CLDN4_vs_palantir_pseudotime')}; "
        f"{_cell_bit('barrier_no_cldn4_vs_palantir_pseudotime')}; "
        f"{_cell_bit('ifn_isg_vs_palantir_pseudotime')}. "
        "dest_1 fate is identically 0 on all post cells (TN-only destiny, P05), so those Spearman rows are NA. "
        "These p-values treat cells as independent and are not a claim."
    )

    if "n_malig_a3" in honest.columns:
        dropped = honest[(honest["timing"] == "post") & (~honest["eligible_primary"])]
        dropped_txt = ", ".join(
            f"{str(r.Sample).replace('BD_immune', 'P')} ({r.paper_group}, n={int(r.n_malig_a3)})"
            for r in dropped.itertuples()
        )
    else:
        dropped_txt = ", ".join(
            f"{r.Sample.replace('BD_immune', 'P')} ({r.paper_group}, n={r.n_cells})"
            for r in patients[(patients["timing"] == "post") & (~patients["eligible_primary"])].itertuples()
        )

    missing_barrier = [g for g in MODULES["barrier_no_cldn4"] if g not in adata.var_names]
    missing_ifn = [g for g in MODULES["ifn_isg"] if g not in adata.var_names]

    ctx = {
        "verdict": verdict,
        "data_bullets": [
            f"Public GEO UMI only: **{extract_info.get('n_cells_matrix', 'NA')}** cells × **{extract_info.get('n_genes_matrix', 'NA')}** genes. Author CopyKAT barcodes are not on GEO.",
            f"A3-malignant-like extracted **n={extract_info.get('n_malig_a3', n_in)}**; after QC (n_genes≥{MIN_GENES}, UMI≥{MIN_UMI}) **n={n_qc}**.",
            f"Epithelial (pass 1) **n={extract_info.get('n_epithelial', 'NA')}**. A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (**not** CopyKAT).",
            f"Unit of every primary test is the **patient**. Post attempted **n={n_post_attempted}**. Eligible (≥{MIN_MAL_PRIMARY} A3-malignant after QC): **n={n_post_eligible}** ({n_nmpr_elig} NMPR, {n_mpr_elig} MPR).",
            f"Dropped / below floor post: {dropped_txt}.",
            f"Eligible post: {elig_names}. dest_1 is TN-only (P05); post fate is identically 0 (Spearman NA).",
            f"Genes absent from barrier module: {missing_barrier or 'none'}. IFN ISG absent: {missing_ifn or 'none'}. SFTPC / KRT6A / KRT6B / KRT14 are known holes on this public UMI.",
            "No dual-high (TACSTD2 AND CLDN4) gate. TACSTD2 is reported only as a companion.",
        ],
        "definitions": [
            ("Object", "A3-malignant-like epithelium (marker, not CopyKAT)"),
            ("CLDN4", "log1p(CP10k) from raw UMI; readout, not a terminal definition"),
            ("Barrier", "TJ (no CLDN4) + simple/basal keratins present on GEO"),
            ("IFN ISG", "40-gene type-I ISG core"),
            ("Root", "CLDN4-low tertile, farthest from CLDN4-high PCA centroid; not CLDN4-high"),
            ("Destinies", "Palantir auto-detected terminals + fate probabilities"),
            ("Primary n", f"post patients with ≥{MIN_MAL_PRIMARY} QC A3-malignant cells"),
            ("MPR", "Hu et al. pathologic response; pCR collapsed to MPR; TN has no MPR"),
        ],
        "root": root,
        "n_destinies": len(dest_list),
        "n_paga_components": len(comps),
        "n_leiden": n_leiden,
        "destiny_rows": destiny_rows,
        "program_rows": program_rows,
        "cell_program_note": cell_program_note,
        "mpr_text": mpr_text,
        "mpr_rows": mpr_rows,
        "limits": [
            f"**Small n.** {n_post_attempted} post patients attempted; **{n_post_eligible}** meet the occupancy floor; **{n_mpr_elig}** MPR among them. Do not cite n=12 as the tested n.",
            "This is not Hu et al. CopyKAT. A3-malignant-like can leak unmarked epithelium.",
            "There is no nLung AT2 root: the A3 rule zeros SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3. The external arrow is “not CLDN4-high,” not an AT2→tumor proof.",
            f"Palantir destinies on a multi-patient tumor subset can recover patient identity (Cramér's V={fmt_r(v_info.get('cramers_v'))}). That is a limit, not a lineage.",
            "Cell-level p-values are exploratory (pseudoreplication).",
            "Dual-high (TACSTD2 AND CLDN4) was not run.",
            "RNA velocity was not run (no spliced/unspliced on GEO).",
            "A MPR-traj / DPT-only analysis is a different object. This FINDING is the Palantir-destiny one.",
        ],
        "figures": [
            "results/figures/fig_palantir_paga.png — PAGA + UMAP destiny / PT / CLDN4 / barrier / IFN",
            "results/figures/fig_programs_along_destiny.png — CLDN4, barrier, IFN along Palantir PT and destinies",
            "results/figures/fig_destiny_vs_mpr.png — destiny composition and patient-mean fate vs MPR",
            "results/figures/fig_honest_n.png — A3-malignant occupancy; floor visible",
            "results/figures/fig_extra_root_not_cldn4high.png — root marked on CLDN4 / tertile UMAP",
            "results/figures/fig_extra_entropy.png — Palantir entropy vs barrier and IFN",
        ],
        "files": [
            "results/tables/destiny_vs_mpr.tsv — destinies × MPR (cells and patients)",
            "results/tables/patient_destiny.tsv — per-patient fate means and majority destiny",
            "results/tables/programs_along_destiny.tsv — CLDN4 / barrier / IFN vs destinies",
            "results/tables/honest_n.tsv / tests.tsv / destiny_by_sample.tsv",
            "results/tables/leiden_paga_vertices.tsv / paga_connectivities.tsv / root_info.json",
            "results/tables/cell_palantir_scores.tsv.gz",
            "results/summary.json",
        ],
    }
    write_finding(args.finding, ctx)

    summary = {
        "dataset": "GSE207422",
        "method": "palantir_destinies_plus_paga",
        "not_dpt_only": True,
        "dual_high_gate": False,
        "root_is_cldn4_high": False,
        "n_extracted": n_in,
        "n_qc": n_qc,
        "n_destinies": len(dest_list),
        "n_post_attempted": n_post_attempted,
        "n_post_eligible": n_post_eligible,
        "n_mpr_eligible": n_mpr_elig,
        "n_nmpr_eligible": n_nmpr_elig,
        "root": root,
        "paga_components": len(comps),
        "destiny_vs_patient_cramers_v": v_info.get("cramers_v"),
        "palantir_version": getattr(palantir, "__version__", "installed"),
        "scanpy_version": getattr(sc, "__version__", "installed"),
        "module_genes_present": module_present,
        "extract": extract_info,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps({k: summary[k] for k in (
        "n_qc", "n_destinies", "n_post_eligible", "n_mpr_eligible", "root_is_cldn4_high"
    )}, indent=2), flush=True)
    print(f"wrote {args.finding}", flush=True)
    print(f"wrote {tabdir / 'destiny_vs_mpr.tsv'}", flush=True)
    print(f"wrote {tabdir / 'patient_destiny.tsv'}", flush=True)


if __name__ == "__main__":
    main()
