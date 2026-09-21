#!/usr/bin/env python3
"""Wider Milo sweep: maximize T/NK neighbourhoods down in CLDN4-high-rich.

The locked primary model is unchanged (k=30, d=30, Q3+Q4 vs Q1+Q2,
n=65). This script does not replace that fit. It searches a wider
grid and reports the specifications that increase the T/NK-down count
and the median |log2FC| of those neighbourhoods.

Cohort and unit, fixed before the fits:
- The same 65 locked units (GSE123902 13, GSE131907 21, GSE205335 22,
  GSE189357 9). No cohort is added or dropped.
- Every contrast keeps all 65 units in the model. A quantile cut
  labels units high or low; it does not delete the middle.
- The sample in the GLM is the locked unit, not the cell.

Graph, same rules as the primary:
- Euclidean kNN on the first d Harmony dimensions, undirected,
  Milo median refinement, prop=0.1, seed=1.
- A neighbourhood is tested when cells from at least 5 of the 65
  units fall in it.
- edgeR glmQLFit(robust=TRUE) + glmQLFTest, TMM, library size =
  column sums, design ~ dataset + high. SpatialFDR is k-distance.

Objective, fixed before reading this grid. Eligible rows are
full-rank fits with n_samples=65, n_high>=10 and n_low>=10.
A hit is a T/NK-majority neighbourhood with SpatialFDR < 0.05 and
log2FC < 0. Two quantities are maximised together:
- the count of those hits
- the median |log2FC| of those hits
The Pareto front of that pair is the result. Three labeled rows are
also written: the count maximum, the |log2FC| maximum among rows
with at least 20 hits, and the maximum of count times that median.
Ties break toward the other quantity, then toward a more negative
median log2FC of all tested T/NK neighbourhoods.

SpatialFDR cuts 0.01, 0.10 and 0.20 are recorded on the same fits.
They are not the selection rule. Loosening the cut increases the
count by definition.
"""

from __future__ import annotations

import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
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
DEEP = CACHE / "max_effect"
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "figures"
EDGE = ROOT / "scripts" / "edger_qlf.R"

PROP = 0.1
SEED = 1
MIN_PATIENTS = 5
SFDR_SELECT = 0.05
SFDR_EXTRA = (0.01, 0.10, 0.20)
MIN_HITS_FOR_LFC = 20
MIN_ARM = 10
N_WORKERS = 3

# Searched range. The previous grid was k in {15, 30, 50}, d in {10, 20, 30}.
KS = (5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60)
DS = (3, 5, 8, 10, 15, 20, 30)

CLASS_NAMES = ["B", "T", "NK", "myeloid", "malignant", "other"]


def class_group_from_counts(cc: np.ndarray) -> np.ndarray:
    labels = []
    for row in cc:
        total = float(row.sum())
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


def within_fraction_high(sm: pd.DataFrame, frac: float) -> pd.Series:
    """Within each dataset, the top `frac` by malignant CLDN4 %pos are high.

    Rank uses method='first' so ties do not move two units together.
    Every unit stays in the contrast.
    """
    high = pd.Series(0, index=sm.index, dtype=int)
    for _, idx in sm.groupby("dataset").groups.items():
        pct = sm.loc[idx, "mal_CLDN4_pct"].astype(float)
        rank = pct.rank(method="first")
        high.loc[idx] = (rank > (len(idx) * (1.0 - frac))).astype(int)
    return high


def build_contrasts(sm: pd.DataFrame) -> list[dict]:
    """Quantile cuts that keep all 65 units. Names are stable."""
    specs = [
        ("q4_vs_rest", (sm["cldn4_quartile"] == "Q4").astype(int), "locked Q4 vs Q1+Q2+Q3"),
        ("q34", sm["cldn4_high"].astype(int), "locked Q3+Q4 vs Q1+Q2"),
    ]
    for frac, name in (
        (0.15, "top15"),
        (0.20, "top20"),
        (0.25, "top25"),
        (1.0 / 3.0, "top33"),
        (0.40, "top40"),
        (0.50, "half"),
    ):
        specs.append((name, within_fraction_high(sm, frac), f"within-cohort top {frac:.3f} vs rest"))
    out = []
    for name, high, rule in specs:
        high = high.astype(int)
        n_high = int(high.sum())
        n_low = int((high == 0).sum())
        out.append(
            {
                "contrast": name,
                "high": high,
                "rule": rule,
                "n_high": n_high,
                "n_low": n_low,
                "n_samples_labeled": int(len(high)),
            }
        )
    return out


def build_graph(X: np.ndarray, k: int):
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", metric="euclidean", n_jobs=-1)
    nn.fit(X)
    dist, idx = nn.kneighbors(X)
    knn_idx = idx[:, 1:].astype(np.int32)
    knn_dist = dist[:, 1:]
    indices = refined_indices(X, knn_idx, PROP, SEED)
    members = undirected_members(knn_idx, indices)
    kdist = knn_dist[indices, -1].astype(np.float64)
    sizes = np.array([len(m) for m in members], dtype=np.int32)
    return members, kdist, sizes


def summarize(da: pd.DataFrame, groups: np.ndarray, kdist: np.ndarray, sizes: np.ndarray) -> dict:
    da = da.copy()
    da["SpatialFDR"] = spatial_fdr_kdistance(da["PValue"].to_numpy(), kdist)
    da["class_group"] = groups
    da["n_cells"] = sizes
    tnk = da[da["class_group"] == "T/NK"]
    mal = da[da["class_group"] == "malignant"]
    down = tnk[(tnk["SpatialFDR"] < SFDR_SELECT) & (tnk["logFC"] < 0)]
    up = tnk[(tnk["SpatialFDR"] < SFDR_SELECT) & (tnk["logFC"] > 0)]
    row = {
        "n_tested": int(len(da)),
        "n_tnk_tested": int(len(tnk)),
        "n_tnk_down_sfdr05": int(len(down)),
        "n_tnk_up_sfdr05": int(len(up)),
        "frac_tnk_down_sfdr05": float(len(down) / len(tnk)) if len(tnk) else np.nan,
        "median_logfc_tnk": float(tnk["logFC"].median()) if len(tnk) else np.nan,
        "median_abs_logfc_tnk_down": float(down["logFC"].abs().median()) if len(down) else np.nan,
        "median_logfc_tnk_down": float(down["logFC"].median()) if len(down) else np.nan,
        "median_ncells_tnk_down": float(down["n_cells"].median()) if len(down) else np.nan,
        "median_ncells_tnk": float(tnk["n_cells"].median()) if len(tnk) else np.nan,
        "n_mal_tested": int(len(mal)),
        "n_mal_up_sfdr05": int(((mal["SpatialFDR"] < SFDR_SELECT) & (mal["logFC"] > 0)).sum()) if len(mal) else 0,
        "n_mal_down_sfdr05": int(((mal["SpatialFDR"] < SFDR_SELECT) & (mal["logFC"] < 0)).sum()) if len(mal) else 0,
        "median_logfc_mal": float(mal["logFC"].median()) if len(mal) else np.nan,
        "min_spatial_fdr": float(np.nanmin(da["SpatialFDR"])) if len(da) else np.nan,
    }
    for cut in SFDR_EXTRA:
        tag = f"{cut:.2f}".replace(".", "")
        dcut = tnk[(tnk["SpatialFDR"] < cut) & (tnk["logFC"] < 0)]
        row[f"n_tnk_down_sfdr{tag}"] = int(len(dcut))
        row[f"median_abs_logfc_down_sfdr{tag}"] = float(dcut["logFC"].abs().median()) if len(dcut) else np.nan
    row["da"] = da
    return row


def pareto_mask(count: np.ndarray, effect: np.ndarray) -> np.ndarray:
    """True where no other row has count >= and effect >=, one of them strict."""
    n = len(count)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        if not np.isfinite(count[i]) or not np.isfinite(effect[i]):
            keep[i] = False
            continue
        dominates = (
            (count >= count[i])
            & (effect >= effect[i])
            & ((count > count[i]) | (effect > effect[i]))
        )
        dominates[i] = False
        if dominates.any():
            keep[i] = False
    return keep


def main() -> None:
    DEEP.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    (DEEP / "edger").mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(CACHE / "meta.tsv", sep="\t", index_col=0)
    X_all = pd.read_csv(CACHE / "harmony.tsv.gz", sep="\t", index_col=0)
    X_all = X_all.loc[meta.index].to_numpy(dtype=np.float64)
    if X_all.shape[1] < max(DS):
        raise SystemExit(f"Harmony has {X_all.shape[1]} dimensions, grid asks for {max(DS)}")

    sm = pd.read_csv(TABLES / "sample_meta.tsv", sep="\t").set_index("unit_key")
    meta = meta.copy()
    meta["unit_key"] = meta["dataset"].astype(str) + "|" + meta["unit_id"].astype(str)
    samples = [s for s in sm.index if s in set(meta["unit_key"])]
    if len(samples) != 65:
        raise SystemExit(f"expected 65 locked units in the embedding, found {len(samples)}")
    sm = sm.loc[samples]
    contrasts = build_contrasts(sm)
    bal_rows = []
    for c in contrasts:
        sub = sm.copy()
        sub["high"] = c["high"].to_numpy()
        row = {"contrast": c["contrast"], "rule": c["rule"], "n_high": c["n_high"], "n_low": c["n_low"], "n_units": 65}
        for ds in sorted(sub.dataset.unique()):
            row[f"{ds}_high"] = int(((sub.dataset == ds) & (sub.high == 1)).sum())
            row[f"{ds}_low"] = int(((sub.dataset == ds) & (sub.high == 0)).sum())
        bal_rows.append(row)
    pd.DataFrame(bal_rows).to_csv(TABLES / "max_effect_contrast_balance.tsv", sep="\t", index=False)

    sample_index = {s: i for i, s in enumerate(samples)}
    codes = meta["unit_key"].map(sample_index).to_numpy(dtype=np.int32)
    class_map = {c: i for i, c in enumerate(CLASS_NAMES)}
    class_codes = meta["cell_class"].map(class_map).to_numpy(dtype=np.int32)

    meta_path = DEEP / "sample_meta.tsv"
    meta_out = sm.reset_index()
    for c in contrasts:
        meta_out[c["contrast"]] = c["high"].to_numpy()
    meta_out.to_csv(meta_path, sep="\t", index=False)

    jobs = []
    graph_info = {}
    for d in DS:
        Xd = X_all[:, :d]
        for k in KS:
            print(f"graph k={k} d={d}", flush=True)
            members, kdist, sizes = build_graph(Xd, k)
            counts = count_by_sample(members, codes, len(samples))
            cc = class_fractions(members, class_codes, len(CLASS_NAMES))
            groups = class_group_from_counts(cc)
            key = f"k{k}_d{d}"
            n_pat = (counts > 0).sum(axis=1)
            keep = n_pat >= MIN_PATIENTS
            ids = np.array([f"{key}_nh{i}" for i in range(counts.shape[0])])
            tested_ids = ids[keep]
            cdf = pd.DataFrame(counts[keep], index=tested_ids, columns=samples)
            cdf.index.name = "nhood_id"
            cpath = DEEP / f"counts_{key}.tsv"
            cdf.to_csv(cpath, sep="\t")
            graph_info[key] = {
                "k": k,
                "d": d,
                "kdist": dict(zip(tested_ids.tolist(), kdist[keep].tolist())),
                "group": dict(zip(tested_ids.tolist(), groups[keep].tolist())),
                "size": dict(zip(tested_ids.tolist(), sizes[keep].tolist())),
                "n_nhoods": int(counts.shape[0]),
                "n_tested_rows": int(keep.sum()),
                "median_ncells": float(np.median(sizes)),
            }
            for c in contrasts:
                jobs.append(
                    {
                        "name": f"{key}__{c['contrast']}",
                        "counts": str(cpath),
                        "meta": str(meta_path),
                        "formula": f"~ dataset + {c['contrast']}",
                        "coef": c["contrast"],
                        "norm": "TMM",
                        "lib_mode": "colsum",
                        "lib_column": "",
                        "sample_column": "unit_key",
                        "keep_column": "",
                        "key": key,
                        "contrast": c["contrast"],
                    }
                )
            del members, counts, cc

    # Shard edgeR jobs so a few R processes run at once.
    shards = [jobs[i::N_WORKERS] for i in range(N_WORKERS)]
    futs = []
    with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
        for i, shard in enumerate(shards):
            if not shard:
                continue
            spec = [
                {kk: job[kk] for kk in (
                    "name", "counts", "meta", "formula", "coef", "norm", "lib_mode",
                    "lib_column", "sample_column", "keep_column",
                )}
                for job in shard
            ]
            spec_path = DEEP / f"models_{i}.json"
            spec_path.write_text(json.dumps(spec))
            futs.append(
                pool.submit(
                    subprocess.check_call,
                    ["Rscript", str(EDGE), str(spec_path), str(DEEP / "edger")],
                )
            )
        for fut in as_completed(futs):
            fut.result()

    rows = []
    for job in jobs:
        design_path = DEEP / "edger" / f"{job['name']}.design.tsv"
        da_path = DEEP / "edger" / f"{job['name']}.tsv"
        info = graph_info[job["key"]]
        design = pd.read_csv(design_path, sep="\t").iloc[0].to_dict() if design_path.exists() else {}
        cinfo = next(c for c in contrasts if c["contrast"] == job["contrast"])
        base = {
            "k": info["k"],
            "d": info["d"],
            "contrast": job["contrast"],
            "rule": cinfo["rule"],
            "n_high": cinfo["n_high"],
            "n_low": cinfo["n_low"],
            "n_nhoods_sampled": info["n_nhoods"],
            "median_ncells_graph": info["median_ncells"],
            "n_samples": int(design.get("n_samples", 0) or 0),
            "rank": int(design.get("rank", -1) or -1),
            "n_coef": int(design.get("n_coef", -1) or -1),
            "singular": bool(design.get("singular", True)),
        }
        if not da_path.exists():
            base.update({"status": "not_fit", "n_tested": 0, "n_tnk_down_sfdr05": 0, "median_abs_logfc_tnk_down": np.nan})
            rows.append(base)
            continue
        da = pd.read_csv(da_path, sep="\t")
        ids = da["nhood_id"].tolist()
        groups = np.array([info["group"][i] for i in ids])
        kdist = np.array([info["kdist"][i] for i in ids], dtype=float)
        sizes = np.array([info["size"][i] for i in ids], dtype=float)
        smry = summarize(da, groups, kdist, sizes)
        da_out = smry.pop("da")
        base.update(smry)
        base["status"] = "fit"
        base["n65"] = bool(base["n_samples"] == 65 and not base["singular"])
        rows.append(base)
        del da_out

    grid = pd.DataFrame(rows)
    eligible = grid[
        (grid["status"] == "fit")
        & (~grid["singular"])
        & (grid["n_samples"] == 65)
        & (grid["n_high"] >= MIN_ARM)
        & (grid["n_low"] >= MIN_ARM)
    ].copy()
    eligible["joint_score"] = eligible["n_tnk_down_sfdr05"] * eligible["median_abs_logfc_tnk_down"]
    eligible["pareto"] = pareto_mask(
        eligible["n_tnk_down_sfdr05"].to_numpy(dtype=float),
        eligible["median_abs_logfc_tnk_down"].to_numpy(dtype=float),
    )
    eligible = eligible.sort_values(
        by=["n_tnk_down_sfdr05", "median_abs_logfc_tnk_down", "median_logfc_tnk"],
        ascending=[False, False, True],
    )
    eligible["selected_count"] = False
    eligible["selected_abs_logfc"] = False
    eligible["selected_joint"] = False
    if len(eligible):
        eligible.iloc[0, eligible.columns.get_loc("selected_count")] = True
        lfc_pool = eligible[eligible["n_tnk_down_sfdr05"] >= MIN_HITS_FOR_LFC]
        if len(lfc_pool):
            lfc_idx = lfc_pool["median_abs_logfc_tnk_down"].idxmax()
            # idxmax on a filtered frame returns the original index
            eligible.loc[lfc_idx, "selected_abs_logfc"] = True
        joint_idx = eligible["joint_score"].idxmax()
        eligible.loc[joint_idx, "selected_joint"] = True

    written = set()

    def dump_selected(flag: str, fname: str) -> None:
        hit = eligible[eligible[flag]]
        if not len(hit):
            return
        r = hit.iloc[0]
        name = f"k{int(r.k)}_d{int(r.d)}__{r.contrast}"
        if name in written:
            return
        written.add(name)
        da_path = DEEP / "edger" / f"{name}.tsv"
        info = graph_info[f"k{int(r.k)}_d{int(r.d)}"]
        da = pd.read_csv(da_path, sep="\t")
        ids = da["nhood_id"].tolist()
        smry = summarize(
            da,
            np.array([info["group"][i] for i in ids]),
            np.array([info["kdist"][i] for i in ids], dtype=float),
            np.array([info["size"][i] for i in ids], dtype=float),
        )
        da = smry["da"]
        da.to_csv(TABLES / fname, sep="\t", index=False)
        down = da[(da.class_group == "T/NK") & (da.SpatialFDR < SFDR_SELECT) & (da.logFC < 0)]
        down.to_csv(TABLES / fname.replace("da_", "tnk_down_"), sep="\t", index=False)

    dump_selected("selected_count", "da_max_effect_count.tsv")
    dump_selected("selected_abs_logfc", "da_max_effect_abslogfc.tsv")
    dump_selected("selected_joint", "da_max_effect_joint.tsv")

    grid.to_csv(TABLES / "max_effect_grid_all.tsv", sep="\t", index=False)
    eligible.to_csv(TABLES / "max_effect_grid_eligible.tsv", sep="\t", index=False)
    front = eligible[eligible["pareto"]].sort_values(
        ["n_tnk_down_sfdr05", "median_abs_logfc_tnk_down"], ascending=[False, False]
    )
    front.to_csv(TABLES / "max_effect_pareto.tsv", sep="\t", index=False)

    # Scatter: count vs median |log2FC|. Pareto points marked.
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    contrasts_order = [c["contrast"] for c in contrasts]
    cmap = plt.get_cmap("tab10")
    for i, name in enumerate(contrasts_order):
        sub = eligible[eligible.contrast == name]
        ax.scatter(
            sub["n_tnk_down_sfdr05"],
            sub["median_abs_logfc_tnk_down"],
            s=22,
            color=cmap(i % 10),
            label=name,
            alpha=0.85,
            linewidths=0,
        )
    if len(front):
        ax.scatter(
            front["n_tnk_down_sfdr05"],
            front["median_abs_logfc_tnk_down"],
            s=90,
            facecolors="none",
            edgecolors="black",
            linewidths=1.1,
            label="Pareto",
            zorder=3,
        )
    for flag, text in (
        ("selected_count", "count and joint"),
        ("selected_abs_logfc", "median |log2FC|"),
    ):
        hit = eligible[eligible[flag]]
        if not len(hit):
            continue
        r = hit.iloc[0]
        ax.annotate(
            f"{text}\nk={int(r.k)}, d={int(r.d)}, {r.contrast}",
            xy=(r.n_tnk_down_sfdr05, r.median_abs_logfc_tnk_down),
            xytext=(8, 8),
            textcoords="offset points",
            fontsize=7,
        )
    ax.set_xlabel("T/NK neighbourhoods down (SpatialFDR < 0.05)")
    ax.set_ylabel("median |log2FC| of those neighbourhoods")
    ax.set_title("Concordant-4 Milo sweep, n = 65 units")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGS / "max_effect_pareto.png", dpi=160)
    fig.savefig(FIGS / "max_effect_pareto.pdf")
    plt.close(fig)

    # Count heatmaps, one panel per contrast.
    fig, axes = plt.subplots(2, 4, figsize=(16.5, 8.4), sharex=True, sharey=True)
    vmax = float(eligible["n_tnk_down_sfdr05"].max()) if len(eligible) else 1
    im = None
    for ax, name in zip(axes.ravel(), contrasts_order):
        sub = eligible[eligible.contrast == name]
        mat = np.full((len(KS), len(DS)), np.nan)
        for i, k in enumerate(KS):
            for j, d in enumerate(DS):
                hit = sub[(sub.k == k) & (sub.d == d)]
                if len(hit):
                    mat[i, j] = hit.iloc[0]["n_tnk_down_sfdr05"]
        im = ax.imshow(mat, cmap="Blues", vmin=0, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(DS)), [str(d) for d in DS], fontsize=8)
        ax.set_yticks(range(len(KS)), [str(k) for k in KS], fontsize=8)
        ax.set_title(name, fontsize=11)
        for i in range(len(KS)):
            for j in range(len(DS)):
                if np.isfinite(mat[i, j]):
                    color = "white" if mat[i, j] > 0.55 * vmax else "black"
                    ax.text(j, i, f"{int(mat[i, j])}", ha="center", va="center", fontsize=7, color=color)
    for ax in axes[1, :]:
        ax.set_xlabel("d")
    for ax in axes[:, 0]:
        ax.set_ylabel("k")
    fig.suptitle("T/NK neighbourhoods down in CLDN4-high-rich, SpatialFDR < 0.05, n = 65", fontsize=12)
    fig.tight_layout()
    if im is not None:
        fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02, label="T/NK down")
    fig.savefig(FIGS / "max_effect_count_heatmaps.png", dpi=160)
    fig.savefig(FIGS / "max_effect_count_heatmaps.pdf")
    plt.close(fig)

    summary = {
        "n_units": 65,
        "n_jobs": len(jobs),
        "n_eligible": int(len(eligible)),
        "n_pareto": int(eligible["pareto"].sum()) if len(eligible) else 0,
        "sfdr": SFDR_SELECT,
        "ks": list(KS),
        "ds": list(DS),
        "contrasts": bal_rows,
    }
    if len(eligible):
        for flag, label in (
            ("selected_count", "count"),
            ("selected_abs_logfc", "abs_logfc"),
            ("selected_joint", "joint"),
        ):
            hit = eligible[eligible[flag]]
            if len(hit):
                summary[label] = json.loads(hit.iloc[0].drop(labels=["pareto"]).to_json())
    (TABLES / "max_effect_summary.json").write_text(json.dumps(summary, indent=2))
    print(front[
        ["k", "d", "contrast", "n_high", "n_low", "n_tnk_tested", "n_tnk_down_sfdr05",
         "median_abs_logfc_tnk_down", "median_logfc_tnk", "joint_score",
         "selected_count", "selected_abs_logfc", "selected_joint"]
    ].to_string(index=False))
    print("--- count / |logFC| / joint winners ---")
    cols = ["k", "d", "contrast", "n_tnk_down_sfdr05", "median_abs_logfc_tnk_down", "median_logfc_tnk", "n_tnk_up_sfdr05", "joint_score"]
    print(eligible[eligible.selected_count | eligible.selected_abs_logfc | eligible.selected_joint][cols].to_string(index=False))


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        sys.exit(exc.returncode)
