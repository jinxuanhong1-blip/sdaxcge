#!/usr/bin/env python3
"""Squidpy FOV graphs for He 2022 CosMx NSCLC (figshare 25976224).

Author tumor cells (tumor 5/6/9/12/13) are split at the within-sample median
of normalized CLDN4. CD8 is author T CD8 memory + T CD8 naive. NK is author NK.
Graphs are Squidpy radius graphs at 50 and 100 µm (0.18 µm per spatial unit).

Does not overwrite the locked exclusion-not-muzzling neighbor ratios.
"""

from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import squidpy as sq
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import binomtest, wilcoxon

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

H5AD = Path("/workspace/data/cosmx_human_nsclc_clustered.h5ad")
OUT = Path("/workspace/results/cosmx_squidpy_fov")
UM_PER_PX = 0.18  # CosMx SMI NSCLC image scale; h5ad spatial is global pixels
RADII_UM = (50.0, 100.0)
SEED = 20260921
N_PERMS_NHOOD = 1000
N_PERMS_LABEL = 999
MIN_ARM = 30
MIN_EFF = 10
TUMOR_TYPES = ("tumor 5", "tumor 6", "tumor 9", "tumor 12", "tumor 13")
CD8_TYPES = ("T CD8 memory", "T CD8 naive")
NK_TYPES = ("NK",)
CATS = ("cldn4_high", "cldn4_low", "cd8", "nk", "other")
EFFECTOR_GENES = ("GZMB", "PRF1", "NKG7", "IFNG")
SAMPLE_ORDER = (
    "LUAD-5 R1",
    "LUAD-5 R2",
    "LUAD-5 R3",
    "LUSC-6",
    "LUAD-9 R1",
    "LUAD-9 R2",
    "LUAD-12",
    "LUAD-13",
)


def _col(mat, idx: int) -> np.ndarray:
    block = mat[:, idx]
    if sparse.issparse(block):
        block = block.toarray()
    return np.asarray(block).ravel()


def load_inputs(path: Path) -> dict:
    a = ad.read_h5ad(path, backed="r")
    if a.n_obs != 765771 or a.n_vars != 960:
        raise SystemExit(f"Unexpected shape {a.shape}; expected 765771 x 960")
    var = list(a.var_names)
    need = ["CLDN4", *EFFECTOR_GENES]
    missing = [g for g in need if g not in var]
    if missing:
        raise SystemExit(f"Missing genes: {missing}")
    obs = a.obs[["sample", "patient", "fov", "cell_type", "n_counts"]].copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["patient"] = obs["patient"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["fov"] = obs["fov"].astype(int)
    spatial_px = np.asarray(a.obsm["spatial"], dtype=np.float64)
    genes = {}
    counts = {}
    for g in need:
        i = var.index(g)
        genes[g] = _col(a.X, i).astype(np.float32)
        counts[g] = _col(a.layers["counts"], i).astype(np.float32)
    n_counts = obs["n_counts"].to_numpy(dtype=np.float64)
    return {
        "obs": obs,
        "spatial_um": spatial_px * UM_PER_PX,
        "cldn4": genes["CLDN4"],
        "expr": genes,
        "counts": counts,
        "n_counts": n_counts,
    }


def assign_labels(obs: pd.DataFrame, cldn4: np.ndarray) -> tuple[np.ndarray, pd.DataFrame]:
    ct = obs["cell_type"].to_numpy()
    sample = obs["sample"].to_numpy()
    is_tumor = np.isin(ct, TUMOR_TYPES)
    is_cd8 = np.isin(ct, CD8_TYPES)
    is_nk = np.isin(ct, NK_TYPES)
    labels = np.full(len(obs), "other", dtype=object)
    labels[is_cd8] = "cd8"
    labels[is_nk] = "nk"
    rows = []
    for s in SAMPLE_ORDER:
        m = (sample == s) & is_tumor
        vals = cldn4[m]
        thr = float(np.median(vals))
        high = m & (cldn4 > thr)
        low = m & (cldn4 <= thr)
        labels[high] = "cldn4_high"
        labels[low] = "cldn4_low"
        rows.append(
            {
                "sample": s,
                "patient": obs.loc[m, "patient"].iloc[0],
                "n_tumor": int(m.sum()),
                "cldn4_median_norm": thr,
                "n_high": int(high.sum()),
                "n_low": int(low.sum()),
                "frac_cldn4_zero_tumor": float(np.mean(vals == 0)),
                "high_rule": "CLDN4 > within-sample tumor median"
                + (" (median is 0, so high = detected)" if thr == 0 else ""),
            }
        )
    if int(is_tumor.sum()) != 302313:
        raise SystemExit(f"Tumor count {int(is_tumor.sum())} != 302313")
    if int(is_cd8.sum()) != 15955 or int(is_nk.sum()) != 7814:
        raise SystemExit(f"CD8/NK counts {int(is_cd8.sum())}/{int(is_nk.sum())}")
    # tumor labels must cover every tumor cell exactly once
    if int(((labels == "cldn4_high") | (labels == "cldn4_low")).sum()) != int(is_tumor.sum()):
        raise SystemExit("CLDN4 split did not cover tumor cells")
    return labels, pd.DataFrame(rows)


def _present_cats(labs: np.ndarray) -> list[str]:
    return [c for c in CATS if np.any(labs == c)]


def _at(mat: np.ndarray, cats: list[str], a: str, b: str, missing: float = float("nan")) -> float:
    if a not in cats or b not in cats:
        return missing
    return float(mat[cats.index(a), cats.index(b)])


def _bin_conn(conn) -> sparse.csr_matrix:
    conn = conn.tocsr()
    return sparse.csr_matrix(
        (np.ones(conn.nnz, dtype=np.float64), conn.indices, conn.indptr),
        shape=conn.shape,
    )


def squidpy_block(
    xy: np.ndarray, labs: np.ndarray, seed: int, do_cooccur: bool, n_jobs: int = 1
) -> list[dict]:
    cats = _present_cats(labs)
    obs = pd.DataFrame({"lab": pd.Categorical(labs, categories=cats)})
    adata = ad.AnnData(X=np.zeros((len(labs), 1), dtype=np.float32), obs=obs)
    adata.obsm["spatial"] = np.asarray(xy, dtype=np.float64)
    adata.obs_names = [str(i) for i in range(adata.n_obs)]
    out = []
    if do_cooccur:
        interval = np.array([0.0, *RADII_UM], dtype=np.float64)
        occ, iv = sq.gr.co_occurrence(adata, cluster_key="lab", interval=interval, copy=True)
        # Squidpy 1.8.3 uses interval[1:] as cumulative distance thresholds.
        thr = [float(x) for x in np.asarray(iv)[1:]]
    else:
        occ, thr = None, []
    for radius in RADII_UM:
        sq.gr.spatial_neighbors_radius(adata, radius=float(radius))
        nh = sq.gr.nhood_enrichment(
            adata,
            cluster_key="lab",
            n_perms=N_PERMS_NHOOD,
            seed=seed,
            copy=True,
            n_jobs=n_jobs,
            show_progress_bar=False,
        )
        counts = sq.gr.interaction_matrix(adata, cluster_key="lab", normalized=False, copy=True)
        fracs = sq.gr.interaction_matrix(adata, cluster_key="lab", normalized=True, copy=True)
        z = np.asarray(nh.zscore, dtype=np.float64)
        counts = np.asarray(counts, dtype=np.float64)
        fracs = np.asarray(fracs, dtype=np.float64)
        n = {c: int(np.sum(labs == c)) for c in CATS}
        row = {
            "radius_um": radius,
            "n_cells": int(len(labs)),
            **{f"n_{c}": n[c] for c in CATS},
        }
        for arm in ("cldn4_high", "cldn4_low"):
            for eff in ("cd8", "nk"):
                row[f"z_{arm}_{eff}"] = _at(z, cats, arm, eff)
                raw = _at(counts, cats, arm, eff, missing=0.0)
                row[f"edges_{arm}_{eff}"] = raw
                row[f"mean_n_{arm}_{eff}"] = raw / n[arm] if n[arm] else float("nan")
                row[f"frac_{arm}_{eff}"] = _at(fracs, cats, arm, eff, missing=0.0)
            row[f"mean_n_{arm}_cd8nk"] = row[f"mean_n_{arm}_cd8"] + row[f"mean_n_{arm}_nk"]
            row[f"frac_{arm}_cd8nk"] = row[f"frac_{arm}_cd8"] + row[f"frac_{arm}_nk"]
        for eff in ("cd8", "nk", "cd8nk"):
            hi = row[f"mean_n_cldn4_high_{eff}"]
            lo = row[f"mean_n_cldn4_low_{eff}"]
            row[f"delta_mean_n_{eff}"] = hi - lo
            row[f"ratio_mean_n_{eff}"] = hi / lo if lo and np.isfinite(lo) and lo != 0 else float("nan")
            row[f"delta_z_{eff}"] = row[f"z_cldn4_high_{eff}"] - row[f"z_cldn4_low_{eff}"] if eff != "cd8nk" else float("nan")
            row[f"delta_frac_{eff}"] = row[f"frac_cldn4_high_{eff}"] - row[f"frac_cldn4_low_{eff}"]
        if occ is not None:
            # occ axis 0 is the anchor category in `cats` order; last axis matches thr
            for j, t in enumerate(thr):
                if abs(t - radius) > 1e-6:
                    continue
                for arm in ("cldn4_high", "cldn4_low"):
                    for eff in ("cd8", "nk"):
                        row[f"occ_{arm}_{eff}"] = _at(occ[:, :, j], cats, arm, eff)
                for eff in ("cd8", "nk"):
                    hi = row[f"occ_cldn4_high_{eff}"]
                    lo = row[f"occ_cldn4_low_{eff}"]
                    row[f"delta_occ_{eff}"] = hi - lo
                    row[f"ratio_occ_{eff}"] = hi / lo if lo and np.isfinite(lo) and lo != 0 else float("nan")
        # CLDN4-label shuffle on this fixed radius graph (tumor cells only).
        bin_conn = _bin_conn(adata.obsp["spatial_connectivities"])
        is_eff = np.isin(labs, ("cd8", "nk"))
        eff_n = np.asarray(bin_conn.dot(is_eff.astype(np.float64))).ravel()
        tumor = np.isin(labs, ("cldn4_high", "cldn4_low"))
        high = labs == "cldn4_high"
        ct = eff_n[tumor]
        ht = high[tumor]
        n_high = int(ht.sum())
        n_low = int((~ht).sum())
        if n_high >= MIN_ARM and n_low >= MIN_ARM and int(is_eff.sum()) >= MIN_EFF:
            obs_delta = float(ct[ht].mean() - ct[~ht].mean())
            rng = np.random.default_rng(seed + int(radius))
            null = np.empty(N_PERMS_LABEL, dtype=np.float64)
            for i in range(N_PERMS_LABEL):
                perm = rng.permutation(ct.shape[0])
                shuffled = ct[perm]
                null[i] = shuffled[:n_high].mean() - shuffled[n_high:].mean()
            # one-sided: high has fewer CD8+NK neighbors than low
            p = float((np.sum(null <= obs_delta) + 1) / (N_PERMS_LABEL + 1))
            row["perm_delta_mean_n_cd8nk"] = obs_delta
            row["perm_p_high_lt_low_cd8nk"] = p
        else:
            row["perm_delta_mean_n_cd8nk"] = float("nan")
            row["perm_p_high_lt_low_cd8nk"] = float("nan")
        out.append(row)
        del bin_conn
    del adata
    return out


def fov_passes(row: dict, eff: str) -> bool:
    if row["n_cldn4_high"] < MIN_ARM or row["n_cldn4_low"] < MIN_ARM:
        return False
    if eff == "cd8nk":
        return (row["n_cd8"] + row["n_nk"]) >= MIN_EFF
    return row[f"n_{eff}"] >= MIN_EFF


def summarize_units(df: pd.DataFrame, unit_cols: list[str]) -> pd.DataFrame:
    """Equal-weight mean of rows inside each unit. One row per unit x radius."""
    rows = []
    for keys, g in df.groupby(unit_cols + ["radius_um"], sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        rec = dict(zip(unit_cols + ["radius_um"], keys))
        rec["n_graphs"] = int(len(g))
        for col in g.columns:
            if col in rec or col in unit_cols or col == "radius_um":
                continue
            if pd.api.types.is_numeric_dtype(g[col]):
                rec[col] = float(np.nanmean(g[col].to_numpy(dtype=float)))
        rows.append(rec)
    return pd.DataFrame(rows)


def sign_summary(deltas: np.ndarray) -> dict:
    d = np.asarray(deltas, dtype=float)
    d = d[np.isfinite(d)]
    n = int(d.size)
    n_neg = int(np.sum(d < 0))
    n_pos = int(np.sum(d > 0))
    n_zero = int(np.sum(d == 0))
    nonzero = d[d != 0]
    if nonzero.size == 0:
        p = float("nan")
    else:
        p = float(binomtest(int(np.sum(nonzero < 0)), int(nonzero.size), 0.5, alternative="greater").pvalue)
    return {
        "n": n,
        "n_high_lt_low": n_neg,
        "n_high_gt_low": n_pos,
        "n_tie": n_zero,
        "median_delta": float(np.median(d)) if n else float("nan"),
        "mean_delta": float(np.mean(d)) if n else float("nan"),
        "sign_p_one_sided_high_lt_low": p,
    }


def wilcoxon_less(high: np.ndarray, low: np.ndarray) -> dict:
    h = np.asarray(high, dtype=float)
    lo = np.asarray(low, dtype=float)
    m = np.isfinite(h) & np.isfinite(lo)
    h, lo = h[m], lo[m]
    out = {"n": int(h.size)}
    if h.size < 5 or np.allclose(h, lo):
        out["wilcoxon_p_less"] = float("nan")
        out["wilcoxon_p_two_sided"] = float("nan")
        return out
    w_less = wilcoxon(h, lo, alternative="less", zero_method="wilcox", method="auto")
    w_two = wilcoxon(h, lo, alternative="two-sided", zero_method="wilcox", method="auto")
    out["wilcoxon_stat_less"] = float(w_less.statistic)
    out["wilcoxon_p_less"] = float(w_less.pvalue)
    out["wilcoxon_p_two_sided"] = float(w_two.pvalue)
    return out


def muzzling_table(bundle: dict, labels: np.ndarray) -> pd.DataFrame:
    obs = bundle["obs"]
    xy = bundle["spatial_um"]
    sample = obs["sample"].to_numpy()
    patient = obs["patient"].to_numpy()
    is_tumor = np.isin(labels, ("cldn4_high", "cldn4_low"))
    is_eff = np.isin(labels, ("cd8", "nk"))
    rows = []
    for s in SAMPLE_ORDER:
        st = (sample == s) & is_tumor
        se = (sample == s) & is_eff
        if st.sum() == 0 or se.sum() == 0:
            continue
        tree = cKDTree(xy[st])
        dist, ix = tree.query(xy[se], k=1, workers=1)
        tumor_lab = labels[st][ix]
        eff_index = np.where(se)[0]
        pat = patient[se][0]
        for radius in RADII_UM:
            near = dist <= radius
            for gene in EFFECTOR_GENES:
                counts = bundle["counts"][gene][eff_index]
                nct = bundle["n_counts"][eff_index]
                cpm = 1e4 * counts / np.maximum(nct, 1.0)
                for arm, bit in (("cldn4_high", near & (tumor_lab == "cldn4_high")), ("cldn4_low", near & (tumor_lab == "cldn4_low"))):
                    rows.append(
                        {
                            "sample": s,
                            "patient": pat,
                            "radius_um": radius,
                            "gene": gene,
                            "arm": arm,
                            "n_effector": int(bit.sum()),
                            "mean_count": float(counts[bit].mean()) if bit.any() else float("nan"),
                            "mean_cpm": float(cpm[bit].mean()) if bit.any() else float("nan"),
                            "frac_detected": float(np.mean(counts[bit] > 0)) if bit.any() else float("nan"),
                        }
                    )
    long = pd.DataFrame(rows)
    # wide hi/lo ratios per sample
    pieces = []
    for (s, radius, gene), g in long.groupby(["sample", "radius_um", "gene"], sort=False):
        hi = g[g.arm == "cldn4_high"].iloc[0]
        lo = g[g.arm == "cldn4_low"].iloc[0]
        pieces.append(
            {
                "sample": s,
                "patient": hi["patient"],
                "radius_um": radius,
                "gene": gene,
                "n_near_high": int(hi["n_effector"]),
                "n_near_low": int(lo["n_effector"]),
                "mean_count_high": hi["mean_count"],
                "mean_count_low": lo["mean_count"],
                "mean_cpm_high": hi["mean_cpm"],
                "mean_cpm_low": lo["mean_cpm"],
                "frac_det_high": hi["frac_detected"],
                "frac_det_low": lo["frac_detected"],
                "ratio_count": _safe_ratio(hi["mean_count"], lo["mean_count"]),
                "ratio_cpm": _safe_ratio(hi["mean_cpm"], lo["mean_cpm"]),
            }
        )
    return long, pd.DataFrame(pieces)


def _safe_ratio(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or den == 0:
        return float("nan")
    return float(num / den)


def donor_from_section(section: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    """Equal weight per section, then one row per patient x radius."""
    rows = []
    for (patient, radius), g in section.groupby(["patient", "radius_um"], sort=False):
        rec = {"patient": patient, "radius_um": radius, "n_sections": int(g["sample"].nunique())}
        for col in value_cols:
            rec[col] = float(np.nanmean(g[col].to_numpy(dtype=float)))
        rows.append(rec)
    return pd.DataFrame(rows)


def make_figures(section: pd.DataFrame, fov: pd.DataFrame, muzz: pd.DataFrame, figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    # Paired section mean CD8+NK neighbor counts
    samples = [s for s in SAMPLE_ORDER if s in set(section["sample"])]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=False)
    for ax, radius in zip(axes, RADII_UM):
        sub = section[section.radius_um == radius].set_index("sample").loc[samples]
        x = np.arange(len(samples))
        ax.plot(x, sub["mean_n_cldn4_low_cd8nk"], "o", color="#4C78A8", label="CLDN4-low tumor")
        ax.plot(x, sub["mean_n_cldn4_high_cd8nk"], "o", color="#F58518", label="CLDN4-high tumor")
        for i in x:
            ax.plot(
                [i, i],
                [sub["mean_n_cldn4_low_cd8nk"].iloc[i], sub["mean_n_cldn4_high_cd8nk"].iloc[i]],
                color="0.6",
                lw=1,
                zorder=0,
            )
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, ha="right")
        ax.set_title(f"{int(radius)} µm radius graph")
        ax.set_ylabel("Mean CD8+NK neighbors per tumor cell")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Squidpy interaction matrix: CD8+NK neighbors of author tumor", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "section_paired_cd8nk_neighbors.png", dpi=160)
    plt.close(fig)

    # Co-occurrence at 50 µm, FOV means per section
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
    use = fov[(fov.radius_um == 50) & (fov.n_cldn4_high >= MIN_ARM) & (fov.n_cldn4_low >= MIN_ARM)].copy()
    for ax, eff, title in zip(axes, ("cd8", "nk"), ("CD8", "NK")):
        rows_h, rows_l = [], []
        samples_c = [s for s in SAMPLE_ORDER if s in set(use["sample"])] or list(SAMPLE_ORDER)
        for s in samples_c:
            g = use[(use["sample"] == s) & (use[f"n_{eff}"] >= MIN_EFF)]
            rows_h.append(np.nanmean(g[f"occ_cldn4_high_{eff}"]) if len(g) else np.nan)
            rows_l.append(np.nanmean(g[f"occ_cldn4_low_{eff}"]) if len(g) else np.nan)
        x = np.arange(len(samples_c))
        ax.plot(x, rows_l, "o", color="#4C78A8", label="around CLDN4-low")
        ax.plot(x, rows_h, "o", color="#F58518", label="around CLDN4-high")
        for i in x:
            if np.isfinite(rows_h[i]) and np.isfinite(rows_l[i]):
                ax.plot([i, i], [rows_l[i], rows_h[i]], color="0.6", lw=1, zorder=0)
        ax.axhline(1.0, color="0.3", ls="--", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(samples_c, rotation=45, ha="right")
        ax.set_title(title)
        ax.set_ylabel("Co-occurrence ratio ≤ 50 µm")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Squidpy co-occurrence (FOV means). 1 = distance-matched base rate", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "section_cooccurrence_50um.png", dpi=160)
    plt.close(fig)

    # Muzzling CPM ratios at 50 µm
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    sub = muzz[muzz.radius_um == 50]
    width = 0.18
    for j, gene in enumerate(EFFECTOR_GENES):
        samples_m = [s for s in SAMPLE_ORDER if s in set(gg_all := sub["sample"])]
        gg = sub[sub.gene == gene].set_index("sample").reindex(samples_m)
        x = np.arange(len(samples_m)) + (j - 1.5) * width
        ax.bar(x, gg["ratio_cpm"], width=width, label=gene)
    ax.axhline(1.0, color="0.2", lw=0.8)
    ax.set_xticks(np.arange(len(samples_m)))
    ax.set_xticklabels(samples_m, rotation=45, ha="right")
    ax.set_ylabel("CPM ratio (near CLDN4-high / near CLDN4-low)")
    ax.set_title("Effector genes in CD8/NK, nearest tumor ≤ 50 µm")
    ax.legend(frameon=False, ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "muzzling_cpm_ratio_50um.png", dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=H5AD)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--limit-fovs", type=int, default=0)
    parser.add_argument("--skip-fov", action="store_true")
    parser.add_argument("--skip-section", action="store_true")
    parser.add_argument("--max-sections", type=int, default=0)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    tab = args.out / "tables"
    fig = args.out / "figures"
    tab.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print("loading", flush=True)
    bundle = load_inputs(args.h5ad)
    labels, thr_df = assign_labels(bundle["obs"], bundle["cldn4"])
    thr_df.to_csv(tab / "cldn4_thresholds.csv", index=False)
    print(thr_df.to_string(index=False), flush=True)
    obs = bundle["obs"]
    xy = bundle["spatial_um"]
    labels_s = pd.Series(labels, index=obs.index)

    fov_rows = []
    if not args.skip_fov:
        groups = list(obs.groupby(["sample", "fov"], sort=True).groups.items())
        if args.limit_fovs:
            groups = groups[: args.limit_fovs]
        print(f"FOV graphs: {len(groups)}", flush=True)
        for i, ((sample, fov), idx) in enumerate(groups, start=1):
            labs = labels_s.loc[idx].to_numpy()
            n_high = int(np.sum(labs == "cldn4_high"))
            n_low = int(np.sum(labs == "cldn4_low"))
            n_cd8 = int(np.sum(labs == "cd8"))
            n_nk = int(np.sum(labs == "nk"))
            patient = obs.loc[idx, "patient"].iloc[0]
            meta = {
                "sample": sample,
                "patient": patient,
                "fov": int(fov),
                "n_high_pre": n_high,
                "n_low_pre": n_low,
                "n_cd8_pre": n_cd8,
                "n_nk_pre": n_nk,
            }
            if n_high < MIN_ARM or n_low < MIN_ARM or (n_cd8 + n_nk) < MIN_EFF:
                meta.update({"skipped": True, "radius_um": np.nan})
                fov_rows.append(meta)
                continue
            t = time.time()
            blocks = squidpy_block(xy[obs.index.get_indexer(idx)], labs, SEED + i, do_cooccur=True)
            for b in blocks:
                fov_rows.append({**meta, "skipped": False, **b})
            print(
                f"FOV {i}/{len(groups)} {sample} fov{fov} {time.time()-t:.1f}s",
                flush=True,
            )
            if i % 25 == 0:
                pd.DataFrame(fov_rows).to_csv(tab / "fov_squidpy.partial.csv", index=False)
        fov_df = pd.DataFrame(fov_rows)
        fov_df.to_csv(tab / "fov_squidpy.csv", index=False)
    else:
        fov_df = pd.read_csv(tab / "fov_squidpy.csv")

    sec_rows = []
    if not args.skip_section:
        print("section graphs", flush=True)
        samples = list(SAMPLE_ORDER)
        if args.max_sections:
            samples = samples[: args.max_sections]
        for sample in samples:
            idx = obs.index[obs["sample"] == sample]
            labs = labels_s.loc[idx].to_numpy()
            t = time.time()
            blocks = squidpy_block(
                xy[obs.index.get_indexer(idx)], labs, SEED, do_cooccur=False, n_jobs=1
            )
            patient = obs.loc[idx, "patient"].iloc[0]
            for b in blocks:
                sec_rows.append({"sample": sample, "patient": patient, "level": "section", **b})
            print(f"section {sample} {time.time()-t:.1f}s n={len(idx)}", flush=True)
        sec_df = pd.DataFrame(sec_rows)
        sec_df.to_csv(tab / "section_squidpy.csv", index=False)
    else:
        sec_df = pd.read_csv(tab / "section_squidpy.csv")

    print("muzzling", flush=True)
    muzz_long, muzz = muzzling_table(bundle, labels)
    muzz_long.to_csv(tab / "muzzling_long.csv", index=False)
    muzz.to_csv(tab / "muzzling_ratios.csv", index=False)

    # FOV filter flags and section means of passing FOVs
    ran = fov_df[fov_df["skipped"] == False].copy() if "skipped" in fov_df.columns else fov_df.copy()
    if "skipped" in ran.columns:
        ran = ran[ran["skipped"] == False]
    fov_pass_rows = []
    for eff in ("cd8", "nk", "cd8nk"):
        for radius in RADII_UM:
            sub = ran[ran.radius_um == radius]
            keep = sub[sub.apply(lambda r: fov_passes(r, eff), axis=1)]
            fov_pass_rows.append(
                {
                    "eff": eff,
                    "radius_um": radius,
                    "n_fov_pass": int(len(keep)),
                    "n_fov_delta_neg": int(np.sum(keep[f"delta_mean_n_{eff}"] < 0)),
                    **{f"fov_{k}": v for k, v in sign_summary(keep[f"delta_mean_n_{eff}"].to_numpy()).items()},
                }
            )
    fov_sign = pd.DataFrame(fov_pass_rows)
    fov_sign.to_csv(tab / "fov_sign_summary.csv", index=False)

    value_cols = [
        c
        for c in sec_df.columns
        if c.startswith(("mean_n_", "frac_", "z_", "delta_", "ratio_", "occ_", "perm_", "n_"))
    ]
    donor = donor_from_section(sec_df, value_cols)
    donor.to_csv(tab / "donor_from_sections.csv", index=False)

    # tests
    tests = {"section_interaction": [], "donor_interaction": [], "muzzling_cpm": [], "cooccurrence_fov": []}
    for radius in RADII_UM:
        sub = sec_df[sec_df.radius_um == radius]
        for eff in ("cd8", "nk", "cd8nk"):
            rec = {
                "radius_um": radius,
                "eff": eff,
                "mean_high": float(sub[f"mean_n_cldn4_high_{eff}"].mean()),
                "mean_low": float(sub[f"mean_n_cldn4_low_{eff}"].mean()),
                "median_ratio": float(sub[f"ratio_mean_n_{eff}"].median()),
                "ratios": sub[f"ratio_mean_n_{eff}"].round(6).tolist(),
                "samples": sub["sample"].tolist(),
                **{f"sign_{k}": v for k, v in sign_summary(sub[f"delta_mean_n_{eff}"].to_numpy()).items()},
                **wilcoxon_less(sub[f"mean_n_cldn4_high_{eff}"], sub[f"mean_n_cldn4_low_{eff}"]),
            }
            tests["section_interaction"].append(rec)
            dsub = donor[donor.radius_um == radius]
            drec = {
                "radius_um": radius,
                "eff": eff,
                "ratios": (dsub[f"mean_n_cldn4_high_{eff}"] / dsub[f"mean_n_cldn4_low_{eff}"]).round(6).tolist(),
                "patients": dsub["patient"].tolist(),
                **sign_summary(
                    (dsub[f"mean_n_cldn4_high_{eff}"] - dsub[f"mean_n_cldn4_low_{eff}"]).to_numpy()
                ),
            }
            tests["donor_interaction"].append(drec)
        gsub = muzz[muzz.radius_um == radius]
        for gene in EFFECTOR_GENES:
            gg = gsub[gsub.gene == gene]
            tests["muzzling_cpm"].append(
                {
                    "radius_um": radius,
                    "gene": gene,
                    "ratios": gg["ratio_cpm"].round(6).tolist(),
                    "samples": gg["sample"].tolist(),
                    "n_ratio_lt_1": int(np.sum(gg["ratio_cpm"] < 1)),
                    "n_ratio_gt_1": int(np.sum(gg["ratio_cpm"] > 1)),
                    "min_ratio": float(gg["ratio_cpm"].min()),
                    "max_ratio": float(gg["ratio_cpm"].max()),
                    "median_ratio": float(gg["ratio_cpm"].median()),
                }
            )
    # co-occurrence sign on FOV means per section (passing FOVs)
    if "occ_cldn4_high_cd8" in ran.columns:
        for radius in RADII_UM:
            sub = ran[ran.radius_um == radius]
            for eff in ("cd8", "nk"):
                sec_means = []
                for s in SAMPLE_ORDER:
                    g = sub[(sub["sample"] == s) & (sub[f"n_{eff}"] >= MIN_EFF) & (sub.n_cldn4_high >= MIN_ARM) & (sub.n_cldn4_low >= MIN_ARM)]
                    if g.empty:
                        continue
                    sec_means.append(
                        {
                            "sample": s,
                            "patient": g["patient"].iloc[0],
                            "delta": float(np.nanmean(g[f"delta_occ_{eff}"])),
                            "occ_high": float(np.nanmean(g[f"occ_cldn4_high_{eff}"])),
                            "occ_low": float(np.nanmean(g[f"occ_cldn4_low_{eff}"])),
                            "n_fov": int(len(g)),
                        }
                    )
                sm = pd.DataFrame(sec_means)
                # donor equal-weight of sections
                dlt = sm.groupby("patient")["delta"].mean()
                tests["cooccurrence_fov"].append(
                    {
                        "radius_um": radius,
                        "eff": eff,
                        "section_occ_high": sm["occ_high"].round(6).tolist(),
                        "section_occ_low": sm["occ_low"].round(6).tolist(),
                        "samples": sm["sample"].tolist(),
                        "section_sign": sign_summary(sm["delta"].to_numpy()),
                        "donor_sign": sign_summary(dlt.to_numpy()),
                        "donor_delta": {k: float(v) for k, v in dlt.items()},
                    }
                )

    summary = {
        "squidpy": sq.__version__,
        "anndata": ad.__version__,
        "n_obs": 765771,
        "um_per_px": UM_PER_PX,
        "radii_um": list(RADII_UM),
        "n_perms_nhood": N_PERMS_NHOOD,
        "n_perms_label_shuffle": N_PERMS_LABEL,
        "seed": SEED,
        "min_arm": MIN_ARM,
        "min_effector": MIN_EFF,
        "tumor_types": list(TUMOR_TYPES),
        "cd8_types": list(CD8_TYPES),
        "nk_types": list(NK_TYPES),
        "n_fov_total": int(obs.groupby(["sample", "fov"]).ngroups),
        "n_fov_ran": int(ran.groupby(["sample", "fov"]).ngroups) if len(ran) else 0,
        "elapsed_sec": time.time() - t0,
        "tests": tests,
        "thresholds": thr_df.to_dict(orient="records"),
    }
    # section permutation p values
    summary["section_perm_p"] = sec_df[
        ["sample", "patient", "radius_um", "perm_delta_mean_n_cd8nk", "perm_p_high_lt_low_cd8nk", "ratio_mean_n_cd8nk", "mean_n_cldn4_high_cd8nk", "mean_n_cldn4_low_cd8nk"]
    ].to_dict(orient="records")
    (tab / "summary.json").write_text(json.dumps(summary, indent=2))
    make_figures(sec_df, ran, muzz, fig)
    print("done", round(time.time() - t0, 1), "s", flush=True)
    print(json.dumps(tests["section_interaction"], indent=2)[:4000])


if __name__ == "__main__":
    main()
