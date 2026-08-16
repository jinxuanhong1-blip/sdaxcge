#!/usr/bin/env python3
"""TACSTD2/CLDN4 vs immune on HTAN/HCA Visium slides (<2GB)."""

from __future__ import annotations

import json
import traceback
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import cKDTree

from lib_htan import (
    CLDN4_IDS,
    TACSTD2_IDS,
    bh_fdr,
    ensure_dirs,
    extract_vector,
    mwu_safe,
    resolve_gene_index,
    signature_scores,
    spearman_safe,
    write_json,
)

VISIUM_KEYS = [
    "htan_visium_Ctrl1_A1",
    "htan_visium_Ctrl1_B1",
    "htan_visium_Ctrl2_C1",
    "htan_visium_Ctrl2_D1",
    "htan_visium_DT1_A1",
    "htan_visium_DT1_B1",
    "htan_visium_DT2_C1",
    "htan_visium_DT2_D1",
    "hca_visium_LngSP10193347",
]


def spatial_coords(adata):
    for key in ("spatial", "X_spatial", "spatial_coords"):
        if key in adata.obsm:
            xy = np.asarray(adata.obsm[key])
            if xy.ndim == 2 and xy.shape[1] >= 2:
                return xy[:, :2]
    for cols in (("array_row", "array_col"), ("x", "y"), ("pxl_col_in_fullres", "pxl_row_in_fullres")):
        if all(c in adata.obs.columns for c in cols):
            return adata.obs[list(cols)].to_numpy(dtype=float)
    return None


def neighbor_mean(xy, values, k=8):
    tree = cKDTree(xy)
    # k+1 includes self; drop self
    _, idx = tree.query(xy, k=min(k + 1, len(xy)))
    if idx.ndim == 1:
        idx = idx[:, None]
    neigh = idx[:, 1:] if idx.shape[1] > 1 else idx
    vals = np.asarray(values, dtype=float)
    return np.nanmean(vals[neigh], axis=1)


def analyze_slide(key: str, path: Path, tables: Path, figs: Path):
    adata = ad.read_h5ad(path, backed="r")
    t_idx, t_name = resolve_gene_index(adata.var, TACSTD2_IDS)
    c_idx, c_name = resolve_gene_index(adata.var, CLDN4_IDS)
    info = {
        "key": key,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "tacstd2_name": t_name,
        "cldn4_name": c_name,
    }
    if t_idx is None or c_idx is None:
        info["error"] = "TACSTD2 or CLDN4 not found"
        write_json(tables / f"{key}_info.json", info)
        if hasattr(adata, "file") and adata.file is not None:
            adata.file.close()
        return info, None

    tac = extract_vector(adata, t_idx)
    cld = extract_vector(adata, c_idx)
    scores, coverage = signature_scores(adata)
    info["signature_coverage"] = coverage
    xy = spatial_coords(adata)
    info["has_spatial"] = xy is not None

    df = pd.DataFrame(
        {
            "TACSTD2": tac,
            "CLDN4": cld,
            **{c: scores[c].to_numpy() for c in scores.columns},
        },
        index=adata.obs_names,
    )
    # in_tissue filter if present
    if "in_tissue" in adata.obs.columns:
        keep = adata.obs["in_tissue"].astype(int).to_numpy() == 1
        df = df.loc[keep]
        if xy is not None:
            xy = xy[keep]
        info["n_in_tissue"] = int(keep.sum())

    rows = []
    r = spearman_safe(df["TACSTD2"], df["CLDN4"])
    r.update({"dataset": key, "x": "TACSTD2", "y": "CLDN4", "kind": "spot"})
    rows.append(r)
    for gene in ("TACSTD2", "CLDN4"):
        for sig in scores.columns:
            rr = spearman_safe(df[gene], df[sig])
            rr.update({"dataset": key, "x": gene, "y": sig, "kind": "spot"})
            rows.append(rr)

    if xy is not None and len(df) >= 20:
        for sig in ("CD8_Tcell", "CYT", "Immune_general", "IFNG_6gene"):
            if sig not in df.columns:
                continue
            neigh = neighbor_mean(xy, df[sig].to_numpy(), k=8)
            df[f"neigh8_{sig}"] = neigh
            for gene in ("TACSTD2", "CLDN4"):
                rr = spearman_safe(df[gene], neigh)
                rr.update({"dataset": key, "x": gene, "y": f"neigh8_{sig}", "kind": "spatial_neighbor"})
                rows.append(rr)

        # high vs low TACSTD2 (median split among expressed-or-all spots)
        hi = df["TACSTD2"] >= np.nanmedian(df["TACSTD2"])
        for sig in ("CD8_Tcell", "CYT", "Immune_general", "IFNG_6gene", "TLS_12chemokine"):
            if sig not in df.columns:
                continue
            mw = mwu_safe(df.loc[hi, sig], df.loc[~hi, sig])
            mw.update({
                "dataset": key, "x": "TACSTD2_high_vs_low", "y": sig,
                "kind": "mwu_spot", "n": int(len(df)),
            })
            rows.append(mw)
            if f"neigh8_{sig}" in df.columns:
                mw2 = mwu_safe(df.loc[hi, f"neigh8_{sig}"], df.loc[~hi, f"neigh8_{sig}"])
                mw2.update({
                    "dataset": key, "x": "TACSTD2_high_vs_low", "y": f"neigh8_{sig}",
                    "kind": "mwu_neighbor", "n": int(len(df)),
                })
                rows.append(mw2)

        # spatial scatter
        sns.set_theme(style="white", font_scale=0.85)
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
        for ax, col, cmap in zip(
            axes,
            ("TACSTD2", "CLDN4", "CD8_Tcell"),
            ("Reds", "Oranges", "Blues"),
        ):
            if col not in df.columns:
                ax.axis("off")
                continue
            sc = ax.scatter(xy[:, 0], xy[:, 1], c=df[col], s=4, cmap=cmap, linewidths=0)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(col)
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        fig.suptitle(key)
        fig.tight_layout()
        fig.savefig(figs / f"{key}_spatial.png", dpi=130)
        plt.close(fig)

    assoc = pd.DataFrame(rows)
    if len(assoc):
        assoc.to_csv(tables / f"{key}_associations.tsv", sep="\t", index=False)
    write_json(tables / f"{key}_info.json", info)
    if hasattr(adata, "file") and adata.file is not None:
        adata.file.close()
    return info, assoc


def stouffer_combine(pvals, effects):
    """Stouffer Z from two-sided p and effect sign. Returns z, p."""
    from scipy.stats import norm

    p = np.asarray(pvals, dtype=float)
    e = np.asarray(effects, dtype=float)
    m = np.isfinite(p) & np.isfinite(e) & (p > 0)
    if m.sum() < 2:
        return np.nan, np.nan, int(m.sum())
    z = norm.isf(np.clip(p[m], 1e-300, 1) / 2.0) * np.sign(e[m])
    zc = z.sum() / np.sqrt(len(z))
    pc = 2 * norm.sf(abs(zc))
    return float(zc), float(pc), int(m.sum())


def main():
    paths = ensure_dirs()
    tables, figs = paths["tables"], paths["figures"]
    all_assoc = []
    runlog = []
    for key in VISIUM_KEYS:
        p = paths["data"] / f"{key}.h5ad"
        if not p.exists():
            runlog.append({"key": key, "status": "missing"})
            continue
        print(f"analyze {key}", flush=True)
        try:
            info, assoc = analyze_slide(key, p, tables, figs)
            runlog.append({"key": key, "status": "ok", "n_obs": info.get("n_obs")})
            if assoc is not None and len(assoc):
                all_assoc.append(assoc)
        except Exception as e:
            runlog.append({"key": key, "status": "error", "error": str(e),
                           "trace": traceback.format_exc()[-800:]})
            print("ERROR", key, e, flush=True)

    if all_assoc:
        cat = pd.concat(all_assoc, ignore_index=True)
        # unify p column (mwu uses p, spearman uses p)
        if "p" not in cat.columns and "pvalue" in cat.columns:
            cat["p"] = cat["pvalue"]
        cat["q"] = bh_fdr(cat["p"].to_numpy() if "p" in cat.columns else np.array([]))
        cat.to_csv(tables / "spatial_associations_all.tsv", sep="\t", index=False)

        # meta across HTAN visium slides only
        htan = cat[cat["dataset"].astype(str).str.startswith("htan_visium")]
        meta = []
        for (kind, x, y), sub in htan.groupby(["kind", "x", "y"], dropna=False):
            effect = sub["rho"] if "rho" in sub.columns else sub.get("rank_biserial")
            z, p, n = stouffer_combine(sub["p"], effect)
            meta.append({
                "kind": kind, "x": x, "y": y, "n_slides": n,
                "median_effect": float(np.nanmedian(effect)) if effect is not None else np.nan,
                "stouffer_z": z, "stouffer_p": p,
            })
        meta_df = pd.DataFrame(meta)
        if len(meta_df):
            meta_df["stouffer_q"] = bh_fdr(meta_df["stouffer_p"].to_numpy())
            meta_df.to_csv(tables / "spatial_htan_visium_meta.tsv", sep="\t", index=False)

        # summary heatmap of spot-level Spearman TACSTD2 vs signatures
        spot = cat[(cat["kind"] == "spot") & (cat["x"].isin(["TACSTD2", "CLDN4"]))]
        if len(spot):
            pv = spot.pivot_table(index="dataset", columns=["x", "y"], values="rho", aggfunc="first")
            fig, ax = plt.subplots(figsize=(12, 4.8))
            sns.heatmap(pv, cmap="RdBu_r", center=0, ax=ax, linewidths=0.2)
            ax.set_title("Visium spot-level Spearman ρ")
            fig.tight_layout()
            fig.savefig(figs / "spatial_spot_spearman_heatmap.png", dpi=140)
            plt.close(fig)

    write_json(paths["notes"] / "spatial_runlog.json", runlog)
    print(json.dumps(runlog, indent=2))


if __name__ == "__main__":
    main()
