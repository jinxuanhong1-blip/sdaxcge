#!/usr/bin/env python3
"""Sensitivity grid for concordant-4 Milo.

Pre-specified search, not a replacement of the primary model
(k=30, d=30, Q3+Q4 vs Q1+Q2).

Objective, fixed before reading the grid: among full-graph specs with
n_high >= 10 and n_low >= 10, maximize the number of T/NK neighbourhoods
with SpatialFDR < 0.05 and log2FC < 0 (down in CLDN4-high-rich).
Ties: more negative median log2FC of those hits, then more negative
median log2FC of all tested T/NK neighbourhoods.

The malignant-only graph is reported separately. It contains no T/NK
cells, so it cannot increase the T/NK-down count.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_milo import (  # noqa: E402
    class_fractions,
    count_by_sample,
    refined_indices,
    spatial_fdr_kdistance,
    undirected_members,
)

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "results" / "cache"
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "figures"
DEEP = CACHE / "deepen"
MIN_PATIENTS = 5
PROP = 0.1
SEED = 1
KS = (15, 30, 50)
DS = (10, 20, 30)
SFDR = 0.05

CLASS_NAMES = ["B", "T", "NK", "myeloid", "malignant", "other"]


def class_group_from_counts(cc: np.ndarray) -> np.ndarray:
    labels = []
    for row in cc:
        total = row.sum()
        fr = row / total if total else row
        order = np.argsort(-fr)
        top = CLASS_NAMES[int(order[0])]
        if total == 0 or fr[order[0]] < 0.5:
            labels.append("mixed")
        elif top in ("T", "NK"):
            labels.append("T/NK")
        else:
            labels.append(top)
    return np.array(labels)


def within_rank_groups(sm: pd.DataFrame) -> pd.DataFrame:
    out = sm.copy()
    out["q4_vs_rest"] = (out["cldn4_quartile"] == "Q4").astype(int)
    out["keep_q4_vs_rest"] = 1
    out["q34"] = out["cldn4_high"].astype(int)
    out["keep_q34"] = 1
    out["q4q1"] = (out["cldn4_quartile"] == "Q4").astype(int)
    out["keep_q4q1"] = out["cldn4_quartile"].isin(["Q1", "Q4"]).astype(int)

    half = pd.Series(0, index=out.index, dtype=int)
    keep_half = pd.Series(1, index=out.index, dtype=int)
    tert_high = pd.Series(0, index=out.index, dtype=int)
    keep_tert = pd.Series(0, index=out.index, dtype=int)
    for ds, idx in out.groupby("dataset").groups.items():
        pct = out.loc[idx, "mal_CLDN4_pct"].astype(float)
        # Upper half vs lower half within cohort. Middle rank of an odd n is low.
        r = pct.rank(method="first")
        half.loc[idx] = (r > (len(idx) / 2.0)).astype(int)
        # Tertiles on the within-cohort rank. Drop the middle tertile.
        labels = pd.qcut(r, 3, labels=["T1", "T2", "T3"])
        tert_high.loc[idx] = (labels == "T3").astype(int)
        keep_tert.loc[idx] = labels.isin(["T1", "T3"]).astype(int)
    out["half"] = half
    out["keep_half"] = keep_half
    out["tert"] = tert_high
    out["keep_tert"] = keep_tert
    return out


CONTRASTS = [
    ("q34", "q34", "keep_q34"),
    ("q4_vs_q1", "q4q1", "keep_q4q1"),
    ("q4_vs_rest", "q4_vs_rest", "keep_q4_vs_rest"),
    ("upper_half", "half", "keep_half"),
    ("tertile", "tert", "keep_tert"),
]


def build_graph(X: np.ndarray, k: int):
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", metric="euclidean")
    nn.fit(X)
    dist, idx = nn.kneighbors(X)
    knn_idx = idx[:, 1:].astype(np.int32)
    knn_dist = dist[:, 1:]
    indices = refined_indices(X, knn_idx, PROP, SEED)
    members = undirected_members(knn_idx, indices)
    kdist = knn_dist[indices, -1]
    return members, kdist


def summarize_job(da: pd.DataFrame, groups: np.ndarray, kdist: np.ndarray) -> dict:
    da = da.copy()
    da["SpatialFDR"] = spatial_fdr_kdistance(da["PValue"].to_numpy(), kdist)
    da["class_group"] = groups
    tnk = da[da["class_group"] == "T/NK"]
    mal = da[da["class_group"] == "malignant"]
    down = tnk[(tnk["SpatialFDR"] < SFDR) & (tnk["logFC"] < 0)]
    up = tnk[(tnk["SpatialFDR"] < SFDR) & (tnk["logFC"] > 0)]
    mal_up = mal[(mal["SpatialFDR"] < SFDR) & (mal["logFC"] > 0)]
    mal_down = mal[(mal["SpatialFDR"] < SFDR) & (mal["logFC"] < 0)]
    return {
        "n_tested": int(len(da)),
        "n_tnk_tested": int(len(tnk)),
        "n_tnk_down_sfdr05": int(len(down)),
        "n_tnk_up_sfdr05": int(len(up)),
        "frac_tnk_down_sfdr05": float(len(down) / len(tnk)) if len(tnk) else np.nan,
        "median_logfc_tnk": float(tnk["logFC"].median()) if len(tnk) else np.nan,
        "median_logfc_tnk_down_hits": float(down["logFC"].median()) if len(down) else np.nan,
        "n_mal_tested": int(len(mal)),
        "n_mal_up_sfdr05": int(len(mal_up)),
        "n_mal_down_sfdr05": int(len(mal_down)),
        "median_logfc_mal": float(mal["logFC"].median()) if len(mal) else np.nan,
        "min_spatial_fdr": float(da["SpatialFDR"].min()) if len(da) else np.nan,
        "da": da,
    }


def main() -> None:
    DEEP.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(CACHE / "meta.tsv", sep="\t", index_col=0)
    X_all = pd.read_csv(CACHE / "harmony.tsv.gz", sep="\t", index_col=0)
    X_all = X_all.loc[meta.index].to_numpy(dtype=np.float64)
    sm = pd.read_csv(TABLES / "sample_meta.tsv", sep="\t")
    sm = sm.set_index("unit_key")
    sm = within_rank_groups(sm)
    meta = meta.copy()
    meta["unit_key"] = meta["dataset"].astype(str) + "|" + meta["unit_id"].astype(str)
    samples = [s for s in sm.index if s in set(meta["unit_key"])]
    sm = sm.loc[samples]
    sm_path = DEEP / "sample_meta_contrasts.tsv"
    sm.reset_index().to_csv(sm_path, sep="\t", index=False)

    # Balance table
    bal_rows = []
    for name, col, keep in CONTRASTS:
        sub = sm[sm[keep] == 1]
        bal_rows.append(
            {
                "contrast": name,
                "n_high": int(sub[col].sum()),
                "n_low": int((sub[col] == 0).sum()),
                **{
                    f"{ds}_high": int(((sub.dataset == ds) & (sub[col] == 1)).sum())
                    for ds in sorted(sub.dataset.unique())
                },
                **{
                    f"{ds}_low": int(((sub.dataset == ds) & (sub[col] == 0)).sum())
                    for ds in sorted(sub.dataset.unique())
                },
            }
        )
    pd.DataFrame(bal_rows).to_csv(TABLES / "contrast_balance.tsv", sep="\t", index=False)

    jobs = []
    graph_info = {}

    def add_graph(tag: str, cell_mask: np.ndarray, ks, ds, malignant_only: bool):
        sub_meta = meta.loc[cell_mask]
        X = X_all[cell_mask]
        # sample index over the full locked sample list; empty samples stay zero
        sample_index = {s: i for i, s in enumerate(samples)}
        codes = sub_meta["unit_key"].map(sample_index).to_numpy(dtype=np.int32)
        class_map = {c: i for i, c in enumerate(CLASS_NAMES)}
        class_codes = sub_meta["cell_class"].map(class_map).to_numpy(dtype=np.int32)
        for d in ds:
            Xd = X[:, :d]
            for k in ks:
                print(f"graph {tag} k={k} d={d} cells={len(sub_meta)}", flush=True)
                members, kdist = build_graph(Xd, k)
                counts = count_by_sample(members, codes, len(samples))
                cc = class_fractions(members, class_codes, len(CLASS_NAMES))
                groups = class_group_from_counts(cc)
                key = f"{tag}_k{k}_d{d}"
                nhood_ids = [f"{key}_nh{i}" for i in range(len(members))]
                cdf = pd.DataFrame(counts, index=nhood_ids, columns=samples)
                cdf.index.name = "nhood_id"
                cpath = DEEP / f"counts_{key}.tsv"
                cdf.to_csv(cpath, sep="\t")
                graph_info[key] = {
                    "kdist": dict(zip(nhood_ids, kdist.tolist())),
                    "group": dict(zip(nhood_ids, groups.tolist())),
                    "malignant_only": malignant_only,
                    "k": k,
                    "d": d,
                    "tag": tag,
                    "n_cells": int(len(sub_meta)),
                    "n_nhoods": len(nhood_ids),
                }
                for cname, col, keep in CONTRASTS:
                    if malignant_only and cname not in ("q34", "q4_vs_q1"):
                        continue
                    jobs.append(
                        {
                            "name": f"{key}__{cname}",
                            "counts": str(cpath),
                            "meta": str(sm_path),
                            "formula": f"~ dataset + {col}",
                            "coef": col,
                            "norm": "TMM",
                            "lib_mode": "colsum",
                            "lib_column": "",
                            "sample_column": "unit_key",
                            "keep_column": keep,
                            "key": key,
                            "contrast": cname,
                        }
                    )

    add_graph("full", np.ones(len(meta), dtype=bool), KS, DS, False)
    mal_mask = (meta["cell_class"] == "malignant").to_numpy()
    add_graph("mal", mal_mask, KS, DS, True)

    # Row filter: >=5 patients inside the contrast, applied by rewriting counts
    # to a tested subset so edgeR dispersion is fit on that set.
    filtered_jobs = []
    for job in jobs:
        counts = pd.read_csv(job["counts"], sep="\t", index_col=0)
        meta_c = pd.read_csv(job["meta"], sep="\t")
        meta_c = meta_c[meta_c[job["keep_column"]].astype(int) == 1]
        use_samples = [s for s in meta_c["unit_key"] if s in counts.columns]
        sub = counts[use_samples]
        n_pat = (sub > 0).sum(axis=1)
        keep_rows = n_pat[n_pat >= MIN_PATIENTS].index
        # also drop all-zero (already implied)
        dest = DEEP / f"tested_{job['name']}.tsv"
        sub.loc[keep_rows].to_csv(dest, sep="\t")
        job = dict(job)
        job["counts"] = str(dest)
        job["n_prefilter"] = int(len(counts))
        job["n_tested_rows"] = int(len(keep_rows))
        filtered_jobs.append(job)

    spec_path = DEEP / "models.json"
    # R script does not need key/contrast/n_* extras; extra keys are ok if
    # fromJSON returns lists. Strip to the fields edger_qlf.R reads.
    r_jobs = []
    for job in filtered_jobs:
        r_jobs.append({k: job[k] for k in (
            "name", "counts", "meta", "formula", "coef", "norm", "lib_mode",
            "lib_column", "sample_column", "keep_column",
        )})
    spec_path.write_text(json.dumps(r_jobs))
    edger_out = DEEP / "edger"
    edger_out.mkdir(parents=True, exist_ok=True)
    print(f"edgeR jobs {len(r_jobs)}", flush=True)
    subprocess.check_call(
        ["Rscript", str(ROOT / "scripts" / "edger_qlf.R"), str(spec_path), str(edger_out)]
    )

    rows = []
    saved = {}
    for job in filtered_jobs:
        design_path = edger_out / f"{job['name']}.design.tsv"
        design = pd.read_csv(design_path, sep="\t").iloc[0].to_dict() if design_path.exists() else {}
        da_path = edger_out / f"{job['name']}.tsv"
        info = graph_info[job["key"]]
        base = {
            "graph": info["tag"],
            "k": info["k"],
            "d": info["d"],
            "contrast": job["contrast"],
            "malignant_only": info["malignant_only"],
            "n_cells": info["n_cells"],
            "n_nhoods_sampled": info["n_nhoods"],
            "n_high": int(design.get("n_samples", 0) and 0),
            "singular": bool(design.get("singular", True)),
            "n_samples": int(design.get("n_samples", 0)) if design else 0,
            "rank": int(design.get("rank", -1)) if design else -1,
            "n_coef": int(design.get("n_coef", -1)) if design else -1,
        }
        # fill n_high/n_low from balance
        if not da_path.exists():
            base.update({
                "n_tested": 0,
                "n_tnk_tested": 0,
                "n_tnk_down_sfdr05": 0,
                "n_tnk_up_sfdr05": 0,
                "frac_tnk_down_sfdr05": np.nan,
                "median_logfc_tnk": np.nan,
                "median_logfc_tnk_down_hits": np.nan,
                "n_mal_tested": 0,
                "n_mal_up_sfdr05": 0,
                "n_mal_down_sfdr05": 0,
                "median_logfc_mal": np.nan,
                "min_spatial_fdr": np.nan,
                "status": "not_fit",
            })
            rows.append(base)
            continue
        da = pd.read_csv(da_path, sep="\t")
        ids = da["nhood_id"].tolist()
        groups = np.array([info["group"][i] for i in ids])
        kdist = np.array([info["kdist"][i] for i in ids], dtype=float)
        smry = summarize_job(da, groups, kdist)
        da_out = smry.pop("da")
        base.update(smry)
        base["status"] = "fit"
        rows.append(base)
        saved[job["name"]] = da_out

    grid = pd.DataFrame(rows)
    # n_high n_low
    bal = pd.read_csv(TABLES / "contrast_balance.tsv", sep="\t")
    grid = grid.drop(columns=["n_high"]).merge(bal[["contrast", "n_high", "n_low"]], on="contrast", how="left")

    full = grid[(grid["graph"] == "full") & (grid["status"] == "fit")].copy()
    eligible = full[(full["n_high"] >= 10) & (full["n_low"] >= 10) & (~full["singular"])].copy()
    eligible = eligible.sort_values(
        by=["n_tnk_down_sfdr05", "median_logfc_tnk_down_hits", "median_logfc_tnk"],
        ascending=[False, True, True],
    )
    eligible["selected"] = False
    if len(eligible):
        eligible.iloc[0, eligible.columns.get_loc("selected")] = True
        winner = eligible.iloc[0]
        win_name = f"full_k{int(winner.k)}_d{int(winner.d)}__{winner.contrast}"
        saved[win_name].to_csv(TABLES / "da_max_tnk_sfdr05.tsv", sep="\t", index=False)
        tnk_hits = saved[win_name]
        tnk_hits = tnk_hits[(tnk_hits.class_group == "T/NK") & (tnk_hits.SpatialFDR < SFDR) & (tnk_hits.logFC < 0)]
        tnk_hits.to_csv(TABLES / "tnk_down_sfdr05_maxspec.tsv", sep="\t", index=False)
    # primary reference slice
    prim_name = "full_k30_d30__q34"
    if prim_name in saved:
        saved[prim_name].to_csv(TABLES / "da_grid_primary_k30_d30_q34.tsv", sep="\t", index=False)
        p = saved[prim_name]
        p[(p.class_group == "T/NK") & (p.SpatialFDR < SFDR) & (p.logFC < 0)].to_csv(
            TABLES / "tnk_down_sfdr05_primary.tsv", sep="\t", index=False
        )

    grid.to_csv(TABLES / "grid_all.tsv", sep="\t", index=False)
    eligible.to_csv(TABLES / "grid_full_ranked.tsv", sep="\t", index=False)

    # Heatmap of T/NK-down counts for full graph
    contrasts = [c[0] for c in CONTRASTS]
    fig, axes = plt.subplots(1, len(contrasts), figsize=(2.4 * len(contrasts), 3.4), sharey=True)
    for ax, cname in zip(axes, contrasts):
        sub = full[full.contrast == cname]
        mat = np.full((len(KS), len(DS)), np.nan)
        for i, k in enumerate(KS):
            for j, d in enumerate(DS):
                hit = sub[(sub.k == k) & (sub.d == d)]
                if len(hit):
                    mat[i, j] = hit.iloc[0]["n_tnk_down_sfdr05"]
        im = ax.imshow(mat, cmap="Blues", vmin=0, vmax=np.nanmax(full["n_tnk_down_sfdr05"]))
        ax.set_xticks(range(len(DS)), [str(d) for d in DS])
        ax.set_yticks(range(len(KS)), [str(k) for k in KS])
        ax.set_xlabel("d")
        ax.set_title(cname, fontsize=9)
        for i in range(len(KS)):
            for j in range(len(DS)):
                if np.isfinite(mat[i, j]):
                    ax.text(j, i, f"{int(mat[i, j])}", ha="center", va="center", fontsize=8)
    axes[0].set_ylabel("k")
    fig.suptitle("T/NK neighbourhoods down in CLDN4-high-rich, SpatialFDR < 0.05", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "grid_tnk_down_sfdr05.png", dpi=160)
    fig.savefig(FIGS / "grid_tnk_down_sfdr05.pdf")
    plt.close(fig)

    # Malignant-only summary at SFDR 0.05
    mal = grid[grid.graph == "mal"].copy()
    mal.to_csv(TABLES / "grid_malignant_only.tsv", sep="\t", index=False)

    print(eligible.head(8).to_string(index=False))
    print("--- malignant q34 ---")
    print(mal[mal.contrast == "q34"][
        ["k", "d", "n_tested", "n_mal_up_sfdr05", "n_mal_down_sfdr05", "median_logfc_mal", "min_spatial_fdr", "status"]
    ].to_string(index=False))


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        sys.exit(exc.returncode)
