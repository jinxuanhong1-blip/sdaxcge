#!/usr/bin/env python3
"""PCA + Harmony Milo DA on GSE207422 epithelium+immune cells.

Contrasts (sample/patient unit):
  1. patient-level malignant TACSTD2 high vs low (median split)
  2. MPR vs NMPR on post-treatment samples if both arms have ≥3 samples

Kept neighborhoods: T/NK-depleted interface nhoods next to TACSTD2-high
epithelium, reported with honest SpatialFDR and n.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from milo_core import count_nhoods, da_nhoods_qlf, make_nhoods

import anndata as ad
import scanpy as sc
import harmonypy as hm

LINEAGES = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T": ["CD3D", "CD3E", "CD2"],
    "NK": ["NKG7", "GNLY", "FGFBP2"],
    "B": ["CD79A", "MS4A1"],
    "plasma": ["IGHG1", "MZB1"],
    "myeloid": ["LYZ", "CD68", "CD14"],
    "neutrophil": ["CSF3R"],
    "mast": ["KIT"],
    "pDC": ["LILRA4"],
    "fibroblast": ["COL1A1", "DCN"],
    "endothelial": ["VWF", "PECAM1"],
}
NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]
IMMUNE = {"T", "NK", "B", "plasma", "myeloid", "neutrophil", "mast", "pDC"}
KEEP_LINEAGES = IMMUNE | {"epithelial"}

# Hu et al. Fig. 1 scRNA groups: pre-biopsy = TN even if later surgery was NMPR.
PAPER_GROUP = {
    "BD_immune01": "TN",
    "BD_immune02": "NMPR",
    "BD_immune03": "MPR",
    "BD_immune04": "NMPR",
    "BD_immune05": "TN",
    "BD_immune06": "MPR",
    "BD_immune07": "NMPR",
    "BD_immune08": "TN",
    "BD_immune09": "NMPR",
    "BD_immune10": "NMPR",
    "BD_immune11": "MPR",
    "BD_immune12": "NMPR",
    "BD_immune13": "NMPR",
    "BD_immune14": "MPR",
    "BD_immune15": "NMPR",
}


def score(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGES)
    scores = np.vstack([score(expr, LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def load_sample_meta(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path)
    raw = raw.dropna(subset=["Sample"]).copy()
    raw = raw[~raw["Sample"].astype(str).str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)]
    raw["Sample"] = raw["Sample"].astype(str)
    raw["paper_group"] = raw["Sample"].map(PAPER_GROUP)
    raw["path_response"] = raw["Pathologic Response"].replace({"pCR": "MPR"})
    raw["timing"] = np.where(
        raw["Resource"].astype(str).str.contains("Pre", case=False, na=False),
        "pre",
        "post",
    )
    return raw


def median_split(values: pd.Series) -> pd.Series:
    med = float(values.median())
    out = pd.Series(index=values.index, dtype=object)
    out[values < med] = "low"
    out[values > med] = "high"
    out[values == med] = np.nan
    return out, med


def nhood_composition(
    nhoods: sparse.csr_matrix,
    lineage: np.ndarray,
    tacstd2_log1p: np.ndarray,
    malignant: np.ndarray,
) -> pd.DataFrame:
    n_h = nhoods.shape[1]
    rows = []
    lin = np.asarray(lineage)
    is_tnk = (lin == "T") | (lin == "NK")
    is_epi = lin == "epithelial"
    for j in range(n_h):
        cells = nhoods[:, j].nonzero()[0]
        n = len(cells)
        n_tnk = int(is_tnk[cells].sum())
        n_epi = int(is_epi[cells].sum())
        n_mal = int(malignant[cells].sum())
        epi_ix = cells[is_epi[cells]]
        tnk_ix = cells[is_tnk[cells]]
        rows.append(
            {
                "nhood": j,
                "n_cells": n,
                "n_TNK": n_tnk,
                "n_epithelial": n_epi,
                "n_malignant_like": n_mal,
                "frac_TNK": n_tnk / n if n else np.nan,
                "frac_epithelial": n_epi / n if n else np.nan,
                "mean_TACSTD2_log1p_all": float(tacstd2_log1p[cells].mean()) if n else np.nan,
                "mean_TACSTD2_log1p_epi": float(tacstd2_log1p[epi_ix].mean()) if n_epi else np.nan,
                "mean_TACSTD2_log1p_tnk": float(tacstd2_log1p[tnk_ix].mean()) if n_tnk else np.nan,
            }
        )
    return pd.DataFrame(rows)


def run_one_embedding(
    adata: ad.AnnData,
    embedding_key: str,
    n_neighbors: int,
    prop: float,
    seed: int,
) -> dict:
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep=embedding_key,
        key_added=None,
    )
    nhoods, index_ixs, kth = make_nhoods(
        adata.obsp["connectivities"],
        adata.obsp["distances"],
        np.asarray(adata.obsm[embedding_key]),
        n_neighbors=n_neighbors,
        prop=prop,
        seed=seed,
    )
    counts, samples = count_nhoods(nhoods, adata.obs["sample"].to_numpy())
    lib = pd.Series(adata.obs["sample"]).value_counts().reindex(samples).fillna(0).to_numpy()
    return {
        "nhoods": nhoods,
        "index_ixs": index_ixs,
        "kth": kth,
        "counts": counts,
        "samples": samples,
        "lib": lib,
        "n_neighbors_used": n_neighbors,
        "n_nhoods": int(nhoods.shape[1]),
        "median_nhood_size": float(np.median(np.asarray(nhoods.sum(axis=0)).ravel())),
    }


def da_for_contrast(bundle: dict, sample_meta: pd.DataFrame, cond_col: str, alt_level: str) -> pd.DataFrame:
    sm = sample_meta.set_index("Sample")
    cond = sm.reindex(bundle["samples"])[cond_col].to_numpy()
    keep = pd.notna(cond)
    if keep.sum() < 4 or len(pd.unique(cond[keep])) < 2:
        return pd.DataFrame()
    return da_nhoods_qlf(
        bundle["counts"][:, keep],
        cond[keep],
        bundle["lib"][keep],
        bundle["kth"],
        alt_level=alt_level,
    )


def select_kept(comp: pd.DataFrame, da: pd.DataFrame, label: str) -> pd.DataFrame:
    """Interface nhoods with TACSTD2-high epithelium and depleted T/NK."""
    df = comp.merge(da, on="nhood", how="left", suffixes=("", "_da"))
    interface = (df["n_epithelial"] >= 3) & (df["n_TNK"] >= 3)
    epi_t = df.loc[interface, "mean_TACSTD2_log1p_epi"].median()
    tnk_q25 = df.loc[interface, "frac_TNK"].quantile(0.25)
    df["is_interface"] = interface
    df["tacstd2_high_epi_nhood"] = interface & (df["mean_TACSTD2_log1p_epi"] >= epi_t)
    df["tnk_depleted_nhood"] = interface & (df["frac_TNK"] <= tnk_q25)
    df["kept_composition"] = df["tacstd2_high_epi_nhood"] & df["tnk_depleted_nhood"]
    # DA support: this depleted-T/NK + high-TACSTD2-epi state enriched in alt (high / MPR)
    df["kept_da_enriched"] = (
        df["kept_composition"]
        & df["testable"].fillna(False)
        & (df["logFC"] > 0)
        & (df["SpatialFDR"] < 0.1)
    )
    # T/NK-rich nhoods depleted in TACSTD2-high patients, adjacent via interface rule
    tnk_rich = df["frac_TNK"] >= 0.5
    df["tnk_rich_depleted_in_alt"] = (
        tnk_rich
        & df["testable"].fillna(False)
        & (df["logFC"] < 0)
        & (df["SpatialFDR"] < 0.1)
    )
    df["keep_reason"] = np.where(
        df["kept_composition"],
        "interface_TACSTD2high_epi_TNKdepleted",
        "",
    )
    df["contrast"] = label
    df["epi_tacstd2_threshold"] = epi_t
    df["tnk_depleted_threshold"] = tnk_q25
    return df


def volcano(df: pd.DataFrame, title: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    x = df["logFC"].to_numpy()
    y = -np.log10(np.clip(df["SpatialFDR"].to_numpy(), 1e-12, 1))
    ax.scatter(x, y, s=8, c="#9aa0a6", alpha=0.6, linewidths=0, label="all testable")
    kept = df["kept_composition"].fillna(False)
    if kept.any():
        ax.scatter(
            x[kept],
            y[kept],
            s=22,
            c="#c0392b",
            alpha=0.9,
            linewidths=0,
            label="kept (T/NK↓ next to TACSTD2-high epi)",
        )
    ax.axhline(-np.log10(0.1), color="#555", ls="--", lw=0.8)
    ax.axvline(0, color="#555", ls=":", lw=0.8)
    ax.set_xlabel("log2 FC (alt vs ref)")
    ax.set_ylabel("−log10 SpatialFDR")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def composition_scatter(df: pd.DataFrame, title: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    ax.scatter(
        df["mean_TACSTD2_log1p_epi"],
        df["frac_TNK"],
        s=10,
        c="#9aa0a6",
        alpha=0.5,
        linewidths=0,
        label="interface or epi-containing",
    )
    kept = df["kept_composition"].fillna(False)
    if kept.any():
        ax.scatter(
            df.loc[kept, "mean_TACSTD2_log1p_epi"],
            df.loc[kept, "frac_TNK"],
            s=24,
            c="#c0392b",
            alpha=0.9,
            linewidths=0,
            label="kept",
        )
    ax.set_xlabel("neighborhood mean log1p TACSTD2 (epithelial cells)")
    ax.set_ylabel("neighborhood T/NK fraction")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument("--outdir", type=Path, default=Path("results/scrna_milo"))
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--n-neighbors", type=int, default=30)
    ap.add_argument("--prop", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    counts = sparse.load_npz(args.workdir / "counts_csr.npz").tocsr()
    pack = np.load(args.workdir / "hvg_marker_counts.npz", allow_pickle=True)
    markers = np.load(args.workdir / "markers.npz", allow_pickle=True)
    cells = pack["cells"]
    genes = [str(g) for g in pack["genes"]]
    total = pack["total"].astype(np.float64)
    expr = {str(g): markers["mat"][i] for i, g in enumerate(markers["genes"])}
    n = len(cells)
    sample_ids = np.array([str(c).rsplit("_", 1)[0] for c in cells])
    lineage = assign_lineage(expr, n)
    normal_umi = np.zeros(n, dtype=np.int32)
    for g in NORMAL_LUNG:
        if g in expr:
            normal_umi += expr[g]
    malignant = (lineage == "epithelial") & (normal_umi == 0)
    tacstd2 = expr["TACSTD2"].astype(np.float64)
    cp10k = np.where(total > 0, tacstd2 / total * 1e4, np.nan)
    tacstd2_log1p = np.log1p(tacstd2)
    tacstd2_log1p_cp10k = np.log1p(cp10k)

    meta = load_sample_meta(args.workdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    keep_cell = np.isin(lineage, list(KEEP_LINEAGES))
    # drop empty-library leftovers
    keep_cell &= total > 0

    adata = ad.AnnData(
        X=counts[keep_cell],
        obs=pd.DataFrame(
            {
                "barcode": cells[keep_cell],
                "sample": sample_ids[keep_cell],
                "lineage": lineage[keep_cell],
                "malignant_like": malignant[keep_cell],
                "total": total[keep_cell],
                "TACSTD2": tacstd2[keep_cell],
                "TACSTD2_log1p": tacstd2_log1p[keep_cell],
                "TACSTD2_log1p_cp10k": tacstd2_log1p_cp10k[keep_cell],
            },
            index=cells[keep_cell].astype(str),
        ),
        var=pd.DataFrame(index=genes),
    )
    adata.obs["sample"] = adata.obs["sample"].astype(str)
    adata.obs = adata.obs.join(meta.set_index("Sample"), on="sample")

    # patient-level malignant TACSTD2 (CP10k log1p), then median split
    mal = adata.obs[adata.obs["malignant_like"]]
    per_sample = (
        mal.groupby("sample")
        .agg(
            n_malignant_like=("malignant_like", "size"),
            mal_TACSTD2_log1p_cp10k=("TACSTD2_log1p_cp10k", "mean"),
            mal_TACSTD2_log1p=("TACSTD2_log1p", "mean"),
            mal_TACSTD2_pct_pos=("TACSTD2", lambda s: float((s > 0).mean() * 100)),
        )
        .reset_index()
        .rename(columns={"sample": "Sample"})
    )
    all_comp = (
        adata.obs.groupby("sample")
        .agg(
            n_epi_immune=("sample", "size"),
            n_TNK=("lineage", lambda s: int(s.isin(["T", "NK"]).sum())),
            n_epithelial=("lineage", lambda s: int((s == "epithelial").sum())),
            frac_TNK=("lineage", lambda s: float(s.isin(["T", "NK"]).mean())),
        )
        .reset_index()
        .rename(columns={"sample": "Sample"})
    )
    sample_tbl = meta.merge(per_sample, on="Sample", how="left").merge(all_comp, on="Sample", how="left")
    eligible_tac = sample_tbl["n_malignant_like"].fillna(0) >= 10
    split, med = median_split(sample_tbl.loc[eligible_tac, "mal_TACSTD2_log1p_cp10k"])
    sample_tbl["tacstd2_class"] = pd.Series(pd.NA, index=sample_tbl.index, dtype="object")
    sample_tbl.loc[eligible_tac, "tacstd2_class"] = split.to_numpy()
    sample_tbl["tacstd2_split_median"] = med
    sample_tbl["mpr_class"] = sample_tbl["paper_group"].where(
        sample_tbl["paper_group"].isin(["MPR", "NMPR"])
    )
    sample_tbl.to_csv(args.outdir / "sample_table.tsv", sep="\t", index=False)

    # normalize + HVG already selected; PCA on the kept genes
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=args.n_pcs, svd_solver="arpack")
    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        "sample",
        max_iter_harmony=20,
        verbose=True,
    )
    z_corr = np.asarray(ho.Z_corr, dtype=np.float32)
    if z_corr.shape == (adata.n_obs, args.n_pcs):
        adata.obsm["X_harmony"] = z_corr
    elif z_corr.shape == (args.n_pcs, adata.n_obs):
        adata.obsm["X_harmony"] = z_corr.T
    else:
        raise ValueError(f"unexpected Harmony Z_corr shape {z_corr.shape}")

    lineage_counts = (
        adata.obs.groupby(["sample", "lineage"]).size().unstack(fill_value=0)
    )
    lineage_counts.to_csv(args.outdir / "lineage_counts_by_sample.tsv", sep="\t")

    summaries = {
        "dataset": "GSE207422",
        "source": "GEO processed UMI matrix (public)",
        "n_cells_matrix": int(n),
        "n_cells_epi_immune": int(adata.n_obs),
        "n_genes_embedding": int(adata.n_vars),
        "n_samples": int(adata.obs["sample"].nunique()),
        "lineage_rule": "argmax mean log1p canonical markers; other if weak/ambiguous",
        "malignant_like_rule": "epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (not CopyKAT)",
        "tacstd2_split": {
            "metric": "sample-level mean log1p(CP10k) TACSTD2 in malignant-like cells",
            "median": med,
            "n_high": int((sample_tbl["tacstd2_class"] == "high").sum()),
            "n_low": int((sample_tbl["tacstd2_class"] == "low").sum()),
            "n_undefined": int(sample_tbl["tacstd2_class"].isna().sum()),
            "eligible_min_malignant_like": 10,
        },
        "mpr": {
            "n_MPR": int((sample_tbl["mpr_class"] == "MPR").sum()),
            "n_NMPR": int((sample_tbl["mpr_class"] == "NMPR").sum()),
            "includes_pCR_as_MPR": True,
            "excludes_treatment_naive_biopsies": True,
            "ran": bool(
                (sample_tbl["mpr_class"] == "MPR").sum() >= 3
                and (sample_tbl["mpr_class"] == "NMPR").sum() >= 3
            ),
        },
        "embeddings": ["X_pca", "X_harmony"],
        "n_pcs": args.n_pcs,
        "n_neighbors": args.n_neighbors,
        "nhood_prop": args.prop,
        "seed": args.seed,
        "da_method": "Python NB-QLF stand-in for miloR/edgeR glmQLFTest; SpatialFDR = cydar weighted BH on 1/kth distance",
        "neighborhood_meaning": "transcriptomic KNN neighborhood in epithelium+immune embedding, not physical space",
    }

    n_high = summaries["tacstd2_split"]["n_high"]
    n_low = summaries["tacstd2_split"]["n_low"]
    run_mpr = summaries["mpr"]["ran"]

    all_kept_rows = []
    da_counts = {}
    for emb_name, emb_key in [("pca", "X_pca"), ("harmony", "X_harmony")]:
        print(f"=== embedding {emb_name} ===", flush=True)
        bundle = run_one_embedding(adata, emb_key, args.n_neighbors, args.prop, args.seed)
        comp = nhood_composition(
            bundle["nhoods"],
            adata.obs["lineage"].to_numpy(),
            adata.obs["TACSTD2_log1p"].to_numpy(),
            adata.obs["malignant_like"].to_numpy(),
        )
        index_lin = adata.obs["lineage"].to_numpy()[bundle["index_ixs"]]
        comp["index_lineage"] = index_lin
        comp["index_cell"] = adata.obs_names.to_numpy()[bundle["index_ixs"]]
        comp["kth_distance"] = bundle["kth"]
        comp["embedding"] = emb_name
        summaries[f"{emb_name}_nhoods"] = {
            "n_nhoods": bundle["n_nhoods"],
            "median_nhood_size": bundle["median_nhood_size"],
            "n_neighbors": bundle["n_neighbors_used"],
        }

        contrasts = []
        if n_high >= 3 and n_low >= 3:
            contrasts.append(("tacstd2_high_vs_low", "tacstd2_class", "high"))
        if run_mpr:
            contrasts.append(("mpr_vs_nmpr", "mpr_class", "MPR"))

        for contrast, col, alt in contrasts:
            print(f"  DA {contrast} on {emb_name}", flush=True)
            da = da_for_contrast(bundle, sample_tbl, col, alt)
            if da.empty:
                summaries.setdefault("skipped", []).append(f"{emb_name}:{contrast}:empty")
                continue
            labeled = select_kept(comp, da, contrast)
            labeled["embedding"] = emb_name
            labeled.to_csv(
                args.outdir / f"nhoods_{emb_name}_{contrast}.tsv",
                sep="\t",
                index=False,
            )
            kept = labeled[labeled["kept_composition"]].copy()
            kept.to_csv(
                args.outdir / f"kept_{emb_name}_{contrast}.tsv",
                sep="\t",
                index=False,
            )
            all_kept_rows.append(kept)
            testable = labeled[labeled["testable"].fillna(False)]
            da_counts[f"{emb_name}:{contrast}"] = {
                "n_nhoods": int(len(labeled)),
                "n_testable": int(len(testable)),
                "n_SpatialFDR_lt_0.1": int((testable["SpatialFDR"] < 0.1).sum()),
                "n_SpatialFDR_lt_0.2": int((testable["SpatialFDR"] < 0.2).sum()),
                "n_FDR_BH_lt_0.1": int((testable["FDR_BH"] < 0.1).sum()),
                "n_kept_composition": int(labeled["kept_composition"].sum()),
                "n_kept_SpatialFDR_lt_0.1": int(
                    (labeled["kept_composition"] & (labeled["SpatialFDR"] < 0.1)).sum()
                ),
                "n_kept_SpatialFDR_lt_0.2": int(
                    (labeled["kept_composition"] & (labeled["SpatialFDR"] < 0.2)).sum()
                ),
                "n_kept_da_enriched_SpatialFDR_lt_0.1": int(labeled["kept_da_enriched"].sum()),
                "n_tnk_rich_depleted_in_alt_SpatialFDR_lt_0.1": int(
                    labeled["tnk_rich_depleted_in_alt"].sum()
                ),
                "n_ref": int(da["n_ref"].iloc[0]) if len(da) else 0,
                "n_alt": int(da["n_alt"].iloc[0]) if len(da) else 0,
                "ref_level": str(da["ref_level"].iloc[0]) if len(da) else None,
                "alt_level": str(da["alt_level"].iloc[0]) if len(da) else None,
            }
            volcano(
                testable,
                f"{emb_name} {contrast} (n_alt={da_counts[f'{emb_name}:{contrast}']['n_alt']}, "
                f"n_ref={da_counts[f'{emb_name}:{contrast}']['n_ref']})",
                args.outdir / f"fig_volcano_{emb_name}_{contrast}.png",
            )
            composition_scatter(
                labeled[labeled["is_interface"]],
                f"{emb_name} {contrast}: interface nhoods",
                args.outdir / f"fig_interface_{emb_name}_{contrast}.png",
            )

    nhood_files = sorted(args.outdir.glob("nhoods_*.tsv"))
    if nhood_files:
        fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.4), sharex=True, sharey=True)
        for ax, path in zip(axes.ravel(), nhood_files):
            nd = pd.read_csv(path, sep="\t")
            tt = nd[nd["testable"].fillna(False)]
            ax.hist(tt["PValue"].dropna(), bins=25, color="#4c6a92", alpha=0.85)
            ax.set_title(path.stem.replace("nhoods_", ""), fontsize=9)
            ax.set_xlabel("PValue")
            ax.set_ylabel("nhoods")
        fig.suptitle("Nominal neighborhood P (not SpatialFDR)", fontsize=11)
        fig.tight_layout()
        fig.savefig(args.outdir / "fig_pvalue_histograms.png", dpi=160)
        plt.close(fig)

    if all_kept_rows:
        kept_all = pd.concat(all_kept_rows, ignore_index=True)
        kept_all.to_csv(args.outdir / "kept_neighborhoods.tsv", sep="\t", index=False)
        # replicate across embeddings for the TACSTD2 contrast
        tac = kept_all[kept_all["contrast"] == "tacstd2_high_vs_low"]
        if not tac.empty:
            both = (
                tac.groupby("index_cell")["embedding"]
                .nunique()
                .reset_index(name="n_embeddings")
            )
            both.to_csv(args.outdir / "kept_index_cells_by_embedding.tsv", sep="\t", index=False)
            summaries["kept_replicated_index_cells"] = int((both["n_embeddings"] == 2).sum())

    summaries["da_counts"] = da_counts
    summaries["lineage_totals"] = {
        k: int(v) for k, v in adata.obs["lineage"].value_counts().to_dict().items()
    }
    summaries["n_malignant_like"] = int(adata.obs["malignant_like"].sum())
    any_s01 = any(v["n_SpatialFDR_lt_0.1"] for v in da_counts.values())
    any_s02 = any(v["n_SpatialFDR_lt_0.2"] for v in da_counts.values())
    any_bh = any(v["n_FDR_BH_lt_0.1"] for v in da_counts.values())
    summaries["verdict"] = {
        "any_SpatialFDR_lt_0.1": bool(any_s01),
        "any_SpatialFDR_lt_0.2": bool(any_s02),
        "any_FDR_BH_lt_0.1": bool(any_bh),
        "kept_neighborhoods_are_compositional_only": True,
        "note": (
            "Nominal P<0.05 can exist; SpatialFDR 0.1/0.2 and BH 0.1 are the "
            "quantities that decide a DA claim. n is samples per arm."
        ),
    }
    da_rows = []
    for key, v in da_counts.items():
        emb, contrast = key.split(":", 1)
        da_rows.append({"embedding": emb, "contrast": contrast, **v})
    if da_rows:
        pd.DataFrame(da_rows).to_csv(args.outdir / "da_counts.tsv", sep="\t", index=False)

    # sample-level sanity (not the Milo test)
    elig = sample_tbl.dropna(subset=["mal_TACSTD2_log1p_cp10k", "frac_TNK"])
    if len(elig) >= 4:
        rho, p = stats.spearmanr(elig["mal_TACSTD2_log1p_cp10k"], elig["frac_TNK"])
        summaries["sanity_sample_spearman_malTACSTD2_vs_TNK"] = {
            "rho": float(rho),
            "p": float(p),
            "n": int(len(elig)),
            "note": "patient unit; not neighborhood DA",
        }

    (args.outdir / "summary.json").write_text(json.dumps(summaries, indent=2))
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
