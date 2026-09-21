#!/usr/bin/env python3
"""Milo on the concordant-4 scVI latent: TACSTD2 vs CLDN4, then TACSTD2 | CLDN4.

Pre-specified before the fit:

Graph. Undirected k=30 nearest neighbours on the 20-d patient-batch scVI latent.
Index cells are the Milo reduced-dimension refinement at prop=0.1, seed=1.
A neighbourhood is tested when it contains cells from at least 5 units.
DA is edgeR glmQLFit(robust=TRUE) + glmQLFTest with TMM, the test miloR uses.
SpatialFDR uses k-distance weights.

Covariates are patient-level, computed on every malignant cell, not on the 10%
sample. CLDN4 high-rich is the locked within-cohort Q3+Q4. TACSTD2 high-rich is
the same quartile rule on malignant TACSTD2 % positive. Dataset is the batch term.

Models, all on the same neighbourhoods:
  tac        ~ dataset + tacstd2_high
  cld        ~ dataset + cldn4_high
  tac|cld    ~ dataset + cldn4_high + tacstd2_high     (coef tacstd2_high)
  cld|tac    ~ dataset + cldn4_high + tacstd2_high     (coef cldn4_high)
Continuous sensitivity replaces the binaries with % positive / 10.

Question. Among T/NK neighbourhoods with SpatialFDR<0.05 and log2FC<0 in `tac`,
does the TACSTD2 log2FC move toward 0 in `tac|cld`? The same summary is reported
for every tested T/NK neighbourhood, and symmetrically for CLDN4 conditional on
TACSTD2. The hit-set Wilcoxon is a description of that selected set, not a
second discovery p-value.

This graph is a transcriptional neighbourhood. It is not a tissue neighbourhood
and it is not spatial exclusion.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import coo_matrix
from sklearn.neighbors import NearestNeighbors

HERE = Path(__file__).resolve().parent
PREP = HERE.parent / "scanvi_concordant4_cldn4"
sys.path.insert(0, str(PREP))
from common import dl_spearman, spearman, within_quartile  # noqa: E402

CACHE = HERE / "results" / "cache"
TABLES = HERE / "results" / "tables"
FIGS = HERE / "figures"
LOCKED = PREP / "data" / "locked_patient_units.tsv"
LOCKED_RHO = {
    "GSE123902": -0.659,
    "GSE131907": -0.522,
    "GSE205335": -0.435,
    "GSE189357": -0.600,
}

K = 30
PROP = 0.1
SEED = 1
MIN_PATIENTS = 5
SFDR = 0.05


def spatial_fdr_kdistance(pvalues: np.ndarray, k_distance: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    d = np.asarray(k_distance, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p) & np.isfinite(d) & (d > 0)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    w = 1.0 / d[ok]
    w[~np.isfinite(w)] = 1.0
    order = np.argsort(pv, kind="mergesort")
    pv_o = pv[order]
    w_o = w[order]
    raw = w_o.sum() * pv_o / np.cumsum(w_o)
    adj_o = np.minimum(1.0, np.minimum.accumulate(raw[::-1])[::-1])
    adj = np.empty_like(adj_o)
    adj[order] = adj_o
    out[ok] = adj
    return out


def refined_indices(X: np.ndarray, knn_idx: np.ndarray, prop: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    n_seed = int(math.floor(prop * n))
    seeds = rng.choice(n, size=n_seed, replace=False)
    med = np.median(X[knn_idx[seeds]], axis=1)
    nn = NearestNeighbors(n_neighbors=1, algorithm="auto", metric="euclidean")
    nn.fit(X)
    idx = nn.kneighbors(med, return_distance=False)[:, 0]
    _, first = np.unique(idx, return_index=True)
    return idx[np.sort(first)].astype(np.int32)


def undirected_members(knn_idx: np.ndarray, indices: np.ndarray) -> list[np.ndarray]:
    n, k = knn_idx.shape
    src = np.repeat(np.arange(n, dtype=np.int32), k)
    dst = knn_idx.ravel()
    mask = src != dst
    src, dst = src[mask], dst[mask]
    rows = np.concatenate([src, dst])
    cols = np.concatenate([dst, src])
    data = np.ones(rows.shape[0], dtype=np.float32)
    A = coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    A.sum_duplicates()
    A.data[:] = 1.0
    members = []
    for i in indices:
        nbrs = A.indices[A.indptr[i] : A.indptr[i + 1]]
        cells = np.unique(np.concatenate(([np.int32(i)], nbrs)))
        members.append(cells)
    return members


def count_by_sample(members: list[np.ndarray], sample_codes: np.ndarray, n_samples: int) -> np.ndarray:
    n_cells = sample_codes.shape[0]
    indptr = [0]
    indices = []
    for cells in members:
        indices.extend(int(c) for c in cells)
        indptr.append(len(indices))
    data = np.ones(len(indices), dtype=np.float64)
    M = coo_matrix(
        (data, (np.repeat(np.arange(len(members)), np.diff(indptr)), indices)),
        shape=(len(members), n_cells),
    ).tocsr()
    S = coo_matrix(
        (np.ones(n_cells), (np.arange(n_cells), sample_codes)),
        shape=(n_cells, n_samples),
    ).tocsr()
    return np.asarray((M @ S).toarray(), dtype=np.float64)


def load_units() -> pd.DataFrame:
    parts = []
    for ds in ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]:
        parts.append(pd.read_csv(PREP / "results" / "cache" / ds / "units.tsv", sep="\t"))
    units = pd.concat(parts, ignore_index=True)
    units["unit_id"] = units["unit_id"].astype(str)
    locked = pd.read_csv(LOCKED, sep="\t")
    locked["unit_id"] = locked["unit_id"].astype(str)
    m = locked.merge(units, on=["dataset", "unit_id"], suffixes=("_locked", "_new"))
    if len(m) != 65:
        raise SystemExit(f"expected 65 locked units, joined {len(m)}")
    d_pct = (m["mal_CLDN4_pct_new"] - m["mal_CLDN4_pct_locked"]).abs()
    if float(d_pct.max()) > 1.0:
        raise SystemExit(f"CLDN4 %pos drifted from the locked table, max |d|={float(d_pct.max())}")
    for ds, rho_locked in LOCKED_RHO.items():
        d = m.loc[m["dataset"] == ds]
        rho = stats.spearmanr(d["mal_CLDN4_pct_new"], d["frac_tnk_new"]).statistic
        if abs(rho - rho_locked) > 0.02:
            raise SystemExit(f"{ds} CLDN4 vs T/NK rho={rho} locked={rho_locked}")
    m["cldn4_high"] = m["cldn4_quartile"].isin(["Q3", "Q4"]).astype(int)
    m["cldn4_per10"] = m["mal_CLDN4_pct_new"].astype(float) / 10.0
    m["tacstd2_per10"] = m["mal_TACSTD2_pct"].astype(float) / 10.0
    q = []
    for _, g in m.groupby("dataset", sort=False):
        lab = within_quartile(g["mal_TACSTD2_pct"])
        q.append(lab)
    m["tacstd2_quartile"] = pd.concat(q).astype(str)
    m["tacstd2_high"] = m["tacstd2_quartile"].isin(["Q3", "Q4"]).astype(int)
    m["unit_key"] = m["dataset"] + "|" + m["unit_id"]
    m["frac_tnk"] = m["frac_tnk_new"]
    m["n_tnk"] = m["n_tnk_new"]
    m["n_cells"] = m["n_cells_new"]
    return m


def partial_spearman(x, y, z) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(ok.sum())
    if n < 5:
        return float("nan"), float("nan"), n
    rx = stats.rankdata(x[ok])
    ry = stats.rankdata(y[ok])
    rz = stats.rankdata(z[ok])

    def resid(a, b):
        design = np.column_stack([np.ones(len(b)), b])
        coef, _, _, _ = np.linalg.lstsq(design, a, rcond=None)
        return a - design @ coef

    r, p = stats.spearmanr(resid(rx, rz), resid(ry, rz))
    return float(r), float(p), n


def patient_level(units: pd.DataFrame) -> dict:
    rows = []
    for xcol in ["mal_CLDN4_pct_new", "mal_TACSTD2_pct"]:
        cohort = []
        for ds, g in units.groupby("dataset"):
            r, p, n = spearman(g[xcol], g["frac_tnk"])
            rows.append({"score": xcol, "dataset": ds, "n": n, "rho": r, "p": p, "level": "cohort"})
            cohort.append((r, n))
        meta = dl_spearman([r for r, _ in cohort], [n for _, n in cohort])
        rows.append({
            "score": xcol, "dataset": "DL_meta", "n": meta["N"], "rho": meta["rho"], "p": meta["p"],
            "level": "meta", "I2": meta["I2"], "ci_lo": meta["ci_lo"], "ci_hi": meta["ci_hi"],
        })
    partial_cohort = []
    for ds, g in units.groupby("dataset"):
        r, p, n = partial_spearman(g["mal_TACSTD2_pct"], g["frac_tnk"], g["mal_CLDN4_pct_new"])
        rows.append({"score": "TACSTD2_partial_CLDN4", "dataset": ds, "n": n, "rho": r, "p": p, "level": "cohort"})
        partial_cohort.append((r, n))
        r2, p2, n2 = partial_spearman(g["mal_CLDN4_pct_new"], g["frac_tnk"], g["mal_TACSTD2_pct"])
        rows.append({"score": "CLDN4_partial_TACSTD2", "dataset": ds, "n": n2, "rho": r2, "p": p2, "level": "cohort"})
    meta_p = dl_spearman([r for r, _ in partial_cohort], [n for _, n in partial_cohort])
    rows.append({
        "score": "TACSTD2_partial_CLDN4", "dataset": "DL_meta", "n": meta_p["N"], "rho": meta_p["rho"],
        "p": meta_p["p"], "level": "meta", "I2": meta_p["I2"], "ci_lo": meta_p["ci_lo"], "ci_hi": meta_p["ci_hi"],
    })
    # symmetric partials were stored per cohort; pool them too
    cld_partial = [r for r in rows if r["score"] == "CLDN4_partial_TACSTD2" and r["level"] == "cohort"]
    meta_c = dl_spearman([r["rho"] for r in cld_partial], [r["n"] for r in cld_partial])
    rows.append({
        "score": "CLDN4_partial_TACSTD2", "dataset": "DL_meta", "n": meta_c["N"], "rho": meta_c["rho"],
        "p": meta_c["p"], "level": "meta", "I2": meta_c["I2"], "ci_lo": meta_c["ci_lo"], "ci_hi": meta_c["ci_hi"],
    })
    r_genes, p_genes, n_genes = spearman(units["mal_TACSTD2_pct"], units["mal_CLDN4_pct_new"])
    # within-cohort DL of the gene-gene Spearman
    gg = []
    for ds, g in units.groupby("dataset"):
        r, p, n = spearman(g["mal_TACSTD2_pct"], g["mal_CLDN4_pct_new"])
        rows.append({"score": "TACSTD2_vs_CLDN4", "dataset": ds, "n": n, "rho": r, "p": p, "level": "cohort"})
        gg.append((r, n))
    meta_g = dl_spearman([r for r, _ in gg], [n for _, n in gg])
    rows.append({
        "score": "TACSTD2_vs_CLDN4", "dataset": "DL_meta", "n": meta_g["N"], "rho": meta_g["rho"],
        "p": meta_g["p"], "level": "meta", "I2": meta_g["I2"], "ci_lo": meta_g["ci_lo"], "ci_hi": meta_g["ci_hi"],
    })
    tab = pd.DataFrame(rows)
    tab.to_csv(TABLES / "patient_level_spearman.tsv", sep="\t", index=False)
    ct = pd.crosstab(units["tacstd2_high"], units["cldn4_high"])
    both = int(((units["tacstd2_high"] == 1) & (units["cldn4_high"] == 1)).sum())
    tac_only = int(((units["tacstd2_high"] == 1) & (units["cldn4_high"] == 0)).sum())
    cld_only = int(((units["tacstd2_high"] == 0) & (units["cldn4_high"] == 1)).sum())
    neither = int(((units["tacstd2_high"] == 0) & (units["cldn4_high"] == 0)).sum())
    # VIF of tacstd2_high on cldn4_high + dataset
    dummies = pd.get_dummies(units["dataset"], drop_first=True, dtype=float)
    design = np.column_stack([
        np.ones(len(units)),
        units["cldn4_high"].to_numpy(dtype=float),
        dummies.to_numpy(dtype=float),
    ])
    y = units["tacstd2_high"].to_numpy(dtype=float)
    coef, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ coef
    ss_res = float(np.sum((y - fitted) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    vif = 1 / (1 - r2) if np.isfinite(r2) and r2 < 1 else float("inf")
    out = {
        "n_units": int(len(units)),
        "n_tacstd2_high": int(units["tacstd2_high"].sum()),
        "n_cldn4_high": int(units["cldn4_high"].sum()),
        "n_both_high": both,
        "n_tacstd2_only": tac_only,
        "n_cldn4_only": cld_only,
        "n_neither": neither,
        "spearman_genes_pooled": {"rho": r_genes, "p": p_genes, "n": n_genes},
        "spearman_genes_dl": meta_g,
        "tacstd2_vs_tnk_dl": tab[(tab["score"] == "mal_TACSTD2_pct") & (tab["dataset"] == "DL_meta")].iloc[0].to_dict(),
        "cldn4_vs_tnk_dl": tab[(tab["score"] == "mal_CLDN4_pct_new") & (tab["dataset"] == "DL_meta")].iloc[0].to_dict(),
        "tacstd2_partial_dl": meta_p,
        "cldn4_partial_dl": meta_c,
        "vif_tacstd2_on_cldn4_and_dataset": float(vif),
        "r2_tacstd2_on_cldn4_and_dataset": float(r2),
        "crosstab": ct.to_dict(),
    }
    return out


def annotate(members, indices, knn_dist, meta, sample_codes, class_names, class_counts) -> pd.DataFrame:
    rows = []
    for i, cells in enumerate(members):
        cc = class_counts[i]
        total = cc.sum()
        fr = cc / total if total else np.zeros_like(cc)
        order = np.argsort(-fr)
        top = class_names[int(order[0])]
        if total == 0 or fr[order[0]] < 0.5 or (len(order) > 1 and fr[order[0]] == fr[order[1]]):
            group = "mixed"
        else:
            group = top
        present = np.unique(sample_codes[cells])
        rows.append({
            "nhood_id": f"nh{i}",
            "index_row": int(indices[i]),
            "k_distance": float(knn_dist[int(indices[i]), -1]),
            "n_cells": int(total),
            "n_patients": int(present.size),
            "majority_class": top,
            "majority_frac": float(fr[order[0]]) if total else float("nan"),
            "class_group": group,
            **{f"n_{c}": int(cc[j]) for j, c in enumerate(class_names)},
            **{f"frac_{c}": float(fr[j]) for j, c in enumerate(class_names)},
        })
    return pd.DataFrame(rows)


def run_edger(jobs: list[dict]) -> None:
    spec = CACHE / "models.json"
    spec.write_text(json.dumps(jobs))
    subprocess.check_call(["Rscript", str(HERE / "edger_qlf.R"), str(spec), str(CACHE / "edger")])


def read_da(name: str, annot: pd.DataFrame) -> pd.DataFrame | None:
    design_path = CACHE / "edger" / f"{name}.design.tsv"
    path = CACHE / "edger" / f"{name}.tsv"
    design = pd.read_csv(design_path, sep="\t")
    if bool(design["singular"].iloc[0]) or not path.exists():
        return None
    da = pd.read_csv(path, sep="\t")
    kd = annot.set_index("nhood_id").loc[da["nhood_id"], "k_distance"].to_numpy()
    da["SpatialFDR"] = spatial_fdr_kdistance(da["PValue"].to_numpy(), kd)
    da["n_samples"] = int(design["n_samples"].iloc[0])
    da["rank"] = int(design["rank"].iloc[0])
    da["n_coef"] = int(design["n_coef"].iloc[0])
    return da


def class_summary(da: pd.DataFrame, annot: pd.DataFrame, label: str) -> pd.DataFrame:
    m = da.merge(annot[["nhood_id", "class_group"]], on="nhood_id", how="left")
    rows = []
    for g, sub in m.groupby("class_group"):
        down = (sub["SpatialFDR"] < SFDR) & (sub["logFC"] < 0)
        up = (sub["SpatialFDR"] < SFDR) & (sub["logFC"] > 0)
        rows.append({
            "model": label,
            "class_group": g,
            "n_tested": int(len(sub)),
            "n_down": int(down.sum()),
            "n_up": int(up.sum()),
            "median_logFC": float(sub["logFC"].median()),
            "median_logFC_down": float(sub.loc[down, "logFC"].median()) if down.any() else float("nan"),
        })
    return pd.DataFrame(rows)


def shrinkage(alone: pd.DataFrame, cond: pd.DataFrame, annot: pd.DataFrame, cls: str) -> dict:
    a = alone.merge(annot[["nhood_id", "class_group"]], on="nhood_id")
    c = cond[["nhood_id", "logFC", "SpatialFDR", "PValue"]].rename(
        columns={"logFC": "logFC_cond", "SpatialFDR": "SpatialFDR_cond", "PValue": "PValue_cond"}
    )
    m = a.merge(c, on="nhood_id")
    m = m.loc[m["class_group"] == cls].copy()
    hits = m[(m["SpatialFDR"] < SFDR) & (m["logFC"] < 0)].copy()

    def block(df: pd.DataFrame) -> dict:
        if len(df) == 0:
            return {"n": 0}
        delta = df["logFC_cond"].to_numpy() - df["logFC"].to_numpy()
        med_a = float(np.median(df["logFC"]))
        med_c = float(np.median(df["logFC_cond"]))
        try:
            w = stats.wilcoxon(df["logFC_cond"], df["logFC"], alternative="two-sided", zero_method="wilcox")
            wp, ws = float(w.pvalue), float(w.statistic)
        except ValueError:
            wp, ws = float("nan"), float("nan")
        still = (df["SpatialFDR_cond"] < SFDR) & (df["logFC_cond"] < 0)
        return {
            "n": int(len(df)),
            "median_logFC_alone": med_a,
            "median_logFC_cond": med_c,
            "median_delta_cond_minus_alone": float(np.median(delta)),
            "frac_abs_smaller": float(np.mean(np.abs(df["logFC_cond"]) < np.abs(df["logFC"]))),
            "shrink_frac_of_medians": float(1 - med_c / med_a) if med_a != 0 else float("nan"),
            "n_still_down": int(still.sum()),
            "wilcoxon_p": wp,
            "wilcoxon_stat": ws,
        }

    return {"class": cls, "hits_down": block(hits), "all_tested": block(m)}


def fmt_rho(d: dict) -> str:
    return (
        f"ρ={d['rho']:.3f} (p={d['p']:.3g}, I²={100 * float(d['I2']):.1f}%, "
        f"{d['ci_lo']:.3f} to {d['ci_hi']:.3f}, N={int(d['n'])})"
    )


def interpret(tac_hits: dict, cld_hits: dict) -> str:
    """Reporting rule fixed in this file: a coefficient shrinks when its median
    absolute value falls below 80% of the unadjusted median. The numbers are
    the result; the sentence only names that comparison."""
    def shrunk(block: dict) -> bool | None:
        if block.get("n", 0) == 0:
            return None
        a = block["median_logFC_alone"]
        c = block["median_logFC_cond"]
        if not np.isfinite(a) or a == 0 or not np.isfinite(c):
            return None
        return abs(c) < 0.8 * abs(a)

    t = shrunk(tac_hits)
    c = shrunk(cld_hits)
    if t is None:
        return "No TACSTD2 T/NK-down neighbourhoods at SpatialFDR < 0.05, so there is no set whose coefficient can shrink."
    if t and c is False:
        return (
            "TACSTD2-associated T/NK-down neighbourhoods shrink conditional on CLDN4 "
            "(median |log2FC| below 80% of the unadjusted median). "
            "The CLDN4-associated T/NK-down set does not shrink by that rule conditional on TACSTD2."
        )
    if t and c:
        return (
            "Both coefficients shrink by the 80% rule. Conditioning does not isolate TACSTD2 from CLDN4; "
            "the shared patient-level variation is large enough that each term moves toward 0 when the other is included."
        )
    if t is False and c:
        return (
            "The TACSTD2 coefficient on its T/NK-down neighbourhoods does not shrink conditional on CLDN4. "
            "The CLDN4 coefficient does shrink conditional on TACSTD2."
        )
    return (
        "The TACSTD2 coefficient on its T/NK-down neighbourhoods does not shrink conditional on CLDN4 "
        "(median |log2FC| stays at or above 80% of the unadjusted median)."
    )


def write_finding(summary: dict, class_tab: pd.DataFrame, units: pd.DataFrame) -> None:
    pl = summary["patient_level"]
    sh = summary["shrink_binary_tnk"]
    hits = sh["hits_down"]
    allc = sh["all_tested"]
    sym = summary["shrink_binary_cldn4_on_tnk"]["hits_down"]
    lines = []
    lines.append("# Milo/scVI concordant-4: TACSTD2 T/NK neighbourhoods conditional on CLDN4")
    lines.append("")
    lines.append(
        "ADDITIVE. Same 65 locked units (GSE123902, GSE131907, GSE205335, GSE189357). "
        "The locked malignant CLDN4 % positive vs T/NK Spearman is unchanged. "
        "This folder asks whether T/NK neighbourhoods that are down in TACSTD2-high units "
        "stay down after CLDN4 is in the same model. The graph is a transcriptional kNN "
        "on the scVI latent, not a tissue distance and not spatial exclusion."
    )
    lines.append("")
    lines.append("## Design")
    lines.append("")
    lines.append(
        f"scVI {summary['train']['scvi_version']}, n_latent={summary['train']['n_latent']}, "
        f"n_layers=2, negative binomial, patient batch, seed 1, {summary['train']['epochs_scvi']} epochs. "
        f"Cells in the graph: {summary['n_cells']}. "
        "Each QC-pass cell was kept with probability 0.10 so embedded T/NK counts stay proportional "
        "to the unit. The 220/140/40 cap used for the earlier scVI integration was not used: that cap "
        "binds for almost every unit and would erase T/NK abundance. "
        f"Embedded T/NK count vs full-unit T/NK count Spearman ρ={summary['embed_vs_full_tnk_rho']:.3f} "
        f"(n={summary['n_units_in_graph']} units present in the graph)."
    )
    lines.append("")
    lines.append(
        f"Milo k=30 on all 20 latent dimensions, prop=0.1, seed=1. "
        f"Tested neighbourhoods (cells from ≥5 units): {summary['n_tested']}. "
        f"edgeR {summary['edger_version']}, TMM, quasi-likelihood F, k-distance SpatialFDR. "
        f"Design rank of the joint model: {summary['joint_rank']}/{summary['joint_n_coef']} "
        f"on {summary['joint_n_samples']} samples."
    )
    lines.append("")
    lines.append(
        f"CLDN4 high-rich is the locked Q3+Q4 ({int(units['cldn4_high'].sum())} vs "
        f"{int((1 - units['cldn4_high']).sum())}). "
        f"TACSTD2 high-rich is the within-cohort quartile Q3+Q4 of malignant TACSTD2 % positive "
        f"({int(units['tacstd2_high'].sum())} vs {int((1 - units['tacstd2_high']).sum())}). "
        f"Discordant units: TACSTD2-only {pl['n_tacstd2_only']}, CLDN4-only {pl['n_cldn4_only']}, "
        f"both high {pl['n_both_high']}, neither {pl['n_neither']}. "
        f"VIF of TACSTD2-high on CLDN4-high plus dataset = {pl['vif_tacstd2_on_cldn4_and_dataset']:.2f}."
    )
    lines.append("")
    lines.append("## Patient-level companion (not Milo)")
    lines.append("")
    lines.append(f"- Malignant CLDN4 %pos vs T/NK: {fmt_rho(pl['cldn4_vs_tnk_dl'])}.")
    lines.append(f"- Malignant TACSTD2 %pos vs T/NK: {fmt_rho(pl['tacstd2_vs_tnk_dl'])}.")
    lines.append(f"- TACSTD2 %pos vs CLDN4 %pos, within-cohort DL: {fmt_rho(pl['spearman_genes_dl'])}.")
    lines.append(
        f"- Partial Spearman, TACSTD2 vs T/NK given CLDN4, within-cohort DL: {fmt_rho(pl['tacstd2_partial_dl'])}."
    )
    lines.append(
        f"- Partial Spearman, CLDN4 vs T/NK given TACSTD2, within-cohort DL: {fmt_rho(pl['cldn4_partial_dl'])}."
    )
    lines.append("")
    lines.append("## Neighbourhood DA")
    lines.append("")
    lines.append("SpatialFDR < 0.05. log2FC is the high arm versus the low arm.")
    lines.append("")
    lines.append("| model | class | tested | down | up | median log2FC |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for rec in class_tab.itertuples(index=False):
        lines.append(
            f"| {rec.model} | {rec.class_group} | {rec.n_tested} | {rec.n_down} | {rec.n_up} | {rec.median_logFC:.3f} |"
        )
    lines.append("")
    lines.append("## Does the TACSTD2 T/NK-down set shrink conditional on CLDN4?")
    lines.append("")
    if hits.get("n", 0) == 0:
        lines.append("No T/NK neighbourhood was down at SpatialFDR < 0.05 in the TACSTD2-only model.")
    else:
        lines.append(
            f"TACSTD2-only T/NK-down neighbourhoods: n={hits['n']}, "
            f"median log2FC {hits['median_logFC_alone']:.3f}. "
            f"The same coefficient with CLDN4 in the model: median log2FC {hits['median_logFC_cond']:.3f} "
            f"(median paired change {hits['median_delta_cond_minus_alone']:.3f}). "
            f"{hits['frac_abs_smaller'] * 100:.1f}% have a smaller absolute coefficient. "
            f"{hits['n_still_down']} remain down at SpatialFDR < 0.05. "
            f"Wilcoxon signed-rank on the paired log2FC, two-sided p={hits['wilcoxon_p']:.3g}. "
            "That p describes this selected set. It is not an independent test."
        )
    lines.append("")
    if allc.get("n", 0):
        lines.append(
            f"All tested T/NK neighbourhoods, not only the hits: n={allc['n']}, "
            f"median TACSTD2 log2FC {allc['median_logFC_alone']:.3f} alone and "
            f"{allc['median_logFC_cond']:.3f} conditional on CLDN4 "
            f"(median paired change {allc['median_delta_cond_minus_alone']:.3f}, "
            f"Wilcoxon p={allc['wilcoxon_p']:.3g})."
        )
        lines.append("")
    if sym.get("n", 0):
        lines.append(
            f"Symmetric set, CLDN4-only T/NK-down neighbourhoods: n={sym['n']}, "
            f"median log2FC {sym['median_logFC_alone']:.3f} alone and "
            f"{sym['median_logFC_cond']:.3f} conditional on TACSTD2 "
            f"({sym['n_still_down']} still down)."
        )
        lines.append("")
    cont = summary.get("shrink_continuous_tnk", {}).get("hits_down", {})
    if cont.get("n", 0):
        lines.append(
            f"Continuous sensitivity (log2FC per 10 percentage points of TACSTD2 %pos). "
            f"T/NK-down set n={cont['n']}, median {cont['median_logFC_alone']:.3f} alone and "
            f"{cont['median_logFC_cond']:.3f} with CLDN4 %pos in the model "
            f"({cont['n_still_down']} still down)."
        )
        lines.append("")
    lines.append(interpret(hits, sym))
    lines.append("")
    lines.append("## What this does not say")
    lines.append("")
    lines.append("- A neighbourhood on the scVI latent is not a radius around a tumour cell.")
    lines.append("- The 10% sample is not the full tissue. Library sizes are about a tenth of the unit.")
    lines.append("- GSE131907 is sampled from the annotation and then QC-filtered. The other three cohorts are QC-filtered and then sampled. Both use the same 10% hash.")
    lines.append("- No fifth cohort. No private 8KL matrix. No Visium co-localisation written as exclusion.")
    lines.append("")
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    (CACHE / "edger").mkdir(parents=True, exist_ok=True)
    units = load_units()
    units.to_csv(TABLES / "patient_units.tsv", sep="\t", index=False)
    pl = patient_level(units)
    print(
        "patient-level",
        f"tac_rho={pl['tacstd2_vs_tnk_dl']['rho']:.3f}",
        f"partial={pl['tacstd2_partial_dl']['rho']:.3f}",
        f"gene_rho={pl['spearman_genes_dl']['rho']:.3f}",
        f"vif={pl['vif_tacstd2_on_cldn4_and_dataset']:.2f}",
        f"discordant={pl['n_tacstd2_only']}/{pl['n_cldn4_only']}",
        flush=True,
    )

    obs = pd.read_csv(CACHE / "obs.tsv.gz", sep="\t")
    X = np.load(CACHE / "latent.npy")
    if len(obs) != len(X):
        raise SystemExit(f"obs {len(obs)} != latent {len(X)}")
    obs["unit_key"] = obs["dataset"].astype(str) + "|" + obs["unit_id"].astype(str)
    keep = obs["unit_key"].isin(set(units["unit_key"])).to_numpy()
    obs = obs.loc[keep].reset_index(drop=True)
    X = X[keep]
    if X.shape[1] != 20:
        raise SystemExit(f"expected 20 latent dimensions, got {X.shape[1]}")

    present = obs.groupby("unit_key").size().rename("n_embedded")
    n_tnk_emb = obs.loc[obs["seed_class"] == "T/NK"].groupby("unit_key").size().rename("n_tnk_embedded")
    diag = units.set_index("unit_key").join(present).join(n_tnk_emb)
    diag["n_embedded"] = diag["n_embedded"].fillna(0)
    diag["n_tnk_embedded"] = diag["n_tnk_embedded"].fillna(0)
    rho_emb = stats.spearmanr(diag["n_tnk"], diag["n_tnk_embedded"]).statistic
    diag.reset_index().to_csv(TABLES / "embedded_counts.tsv", sep="\t", index=False)
    print(f"embed vs full T/NK rho={rho_emb:.3f} cells={len(obs)}", flush=True)

    samples = list(units["unit_key"])
    sample_index = {s: i for i, s in enumerate(samples)}
    sample_codes = obs["unit_key"].map(sample_index).to_numpy(dtype=np.int32)
    print(f"kNN k={K} d=20 cells={len(obs)}", flush=True)
    nn = NearestNeighbors(n_neighbors=K + 1, algorithm="auto", metric="euclidean")
    nn.fit(X)
    dist, idx = nn.kneighbors(X)
    knn_idx = idx[:, 1:].astype(np.int32)
    knn_dist = dist[:, 1:]
    indices = refined_indices(X, knn_idx, PROP, SEED)
    members = undirected_members(knn_idx, indices)
    class_names = ["malignant", "T/NK", "other"]
    class_map = {c: i for i, c in enumerate(class_names)}
    unknown = set(obs["seed_class"]) - set(class_names)
    if unknown:
        raise SystemExit(f"unexpected seed_class: {unknown}")
    class_codes = obs["seed_class"].map(class_map).to_numpy(dtype=np.int32)
    counts = count_by_sample(members, sample_codes, len(samples))
    class_counts = _class_counts(members, class_codes, len(class_names))
    annot = annotate(members, indices, knn_dist, obs, sample_codes, class_names, class_counts)
    tested = annot["n_patients"] >= MIN_PATIENTS
    annot["tested"] = tested.to_numpy()
    annot.to_csv(TABLES / "nhood_annotation.tsv", sep="\t", index=False)
    print(
        f"nhoods {len(annot)} tested {int(tested.sum())} "
        f"T/NK {(annot.loc[tested, 'class_group'] == 'T/NK').sum()}",
        flush=True,
    )
    use = annot.loc[tested, "nhood_id"].tolist()
    # counts row order is nh0.. 
    row_index = [int(n[2:]) for n in use]
    count_df = pd.DataFrame(counts[row_index], index=use, columns=samples)
    count_df.index.name = "nhood_id"
    count_path = CACHE / "nhood_counts.tsv"
    count_df.to_csv(count_path, sep="\t")
    meta = units[["unit_key", "dataset", "tacstd2_high", "cldn4_high", "tacstd2_per10", "cldn4_per10"]].copy()
    meta_path = CACHE / "sample_meta.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)

    jobs = []
    specs = [
        ("tac", "~ dataset + tacstd2_high", "tacstd2_high"),
        ("cld", "~ dataset + cldn4_high", "cldn4_high"),
        ("tac_given_cld", "~ dataset + cldn4_high + tacstd2_high", "tacstd2_high"),
        ("cld_given_tac", "~ dataset + cldn4_high + tacstd2_high", "cldn4_high"),
        ("tac_cont", "~ dataset + tacstd2_per10", "tacstd2_per10"),
        ("cld_cont", "~ dataset + cldn4_per10", "cldn4_per10"),
        ("tac_cont_given_cld", "~ dataset + cldn4_per10 + tacstd2_per10", "tacstd2_per10"),
        ("cld_cont_given_tac", "~ dataset + cldn4_per10 + tacstd2_per10", "cldn4_per10"),
    ]
    for name, formula, coef in specs:
        jobs.append({
            "name": name,
            "counts": str(count_path),
            "meta": str(meta_path),
            "formula": formula,
            "coef": coef,
            "norm": "TMM",
            "lib_mode": "colsum",
            "lib_column": "",
            "sample_column": "unit_key",
            "keep_column": "",
        })
    print("edgeR", flush=True)
    run_edger(jobs)
    das = {}
    for name, _, _ in specs:
        da = read_da(name, annot.loc[tested])
        if da is None:
            print(f"{name} singular or missing", flush=True)
            continue
        da.to_csv(TABLES / f"da_{name}.tsv", sep="\t", index=False)
        das[name] = da
        print(f"{name} min SpatialFDR {da['SpatialFDR'].min():.3g} n={len(da)}", flush=True)

    required = ["tac", "cld", "tac_given_cld", "cld_given_tac"]
    missing = [n for n in required if n not in das]
    if missing:
        raise SystemExit(f"required models missing: {missing}")

    class_tab = pd.concat([
        class_summary(das["tac"], annot, "TACSTD2 high"),
        class_summary(das["cld"], annot, "CLDN4 high"),
        class_summary(das["tac_given_cld"], annot, "TACSTD2 | CLDN4"),
        class_summary(das["cld_given_tac"], annot, "CLDN4 | TACSTD2"),
    ], ignore_index=True)
    class_tab.to_csv(TABLES / "class_summary.tsv", sep="\t", index=False)

    sh_tac = shrinkage(das["tac"], das["tac_given_cld"], annot, "T/NK")
    sh_cld = shrinkage(das["cld"], das["cld_given_tac"], annot, "T/NK")
    sh_cont = {}
    if "tac_cont" in das and "tac_cont_given_cld" in das:
        sh_cont = shrinkage(das["tac_cont"], das["tac_cont_given_cld"], annot, "T/NK")
    # hit-level table
    a = das["tac"].merge(annot[["nhood_id", "class_group", "n_patients"]], on="nhood_id")
    c = das["tac_given_cld"][["nhood_id", "logFC", "SpatialFDR"]].rename(
        columns={"logFC": "logFC_cond", "SpatialFDR": "SpatialFDR_cond"}
    )
    hits = a.merge(c, on="nhood_id")
    hits = hits[(hits["class_group"] == "T/NK") & (hits["SpatialFDR"] < SFDR) & (hits["logFC"] < 0)]
    hits.to_csv(TABLES / "tnk_down_tacstd2.tsv", sep="\t", index=False)

    edger_v = subprocess.check_output(
        ["Rscript", "-e", '.libPaths(c(path.expand("~/R/library"), .libPaths())); cat(as.character(packageVersion("edgeR")))'],
        text=True,
    ).strip()
    train = json.loads((CACHE / "train_info.json").read_text())
    joint = das["tac_given_cld"]
    summary = {
        "train": train,
        "n_cells": int(len(obs)),
        "n_tested": int(tested.sum()),
        "n_units_in_graph": int((diag["n_embedded"] > 0).sum()),
        "embed_vs_full_tnk_rho": float(rho_emb),
        "edger_version": edger_v,
        "joint_rank": int(joint["rank"].iloc[0]),
        "joint_n_coef": int(joint["n_coef"].iloc[0]),
        "joint_n_samples": int(joint["n_samples"].iloc[0]),
        "patient_level": pl,
        "shrink_binary_tnk": sh_tac,
        "shrink_binary_cldn4_on_tnk": sh_cld,
        "shrink_continuous_tnk": sh_cont,
        "interpretation": interpret(sh_tac["hits_down"], sh_cld["hits_down"]),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=_json))
    write_finding(summary, class_tab, units)
    plot(das, annot)
    print(summary["interpretation"], flush=True)


def _class_counts(members, class_codes, n_classes):
    n_cells = class_codes.shape[0]
    indptr = [0]
    indices = []
    for cells in members:
        indices.extend(int(c) for c in cells)
        indptr.append(len(indices))
    data = np.ones(len(indices), dtype=np.float64)
    M = coo_matrix(
        (data, (np.repeat(np.arange(len(members)), np.diff(indptr)), indices)),
        shape=(len(members), n_cells),
    ).tocsr()
    C = coo_matrix(
        (np.ones(n_cells), (np.arange(n_cells), class_codes)),
        shape=(n_cells, n_classes),
    ).tocsr()
    return np.asarray((M @ C).toarray(), dtype=np.float64)


def _json(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, dict):
        return {str(k): _json(v) for k, v in obj.items()}
    raise TypeError(type(obj))


def plot(das, annot) -> None:
    a = das["tac"].merge(annot[["nhood_id", "class_group"]], on="nhood_id")
    c = das["tac_given_cld"][["nhood_id", "logFC"]].rename(columns={"logFC": "logFC_cond"})
    m = a.merge(c, on="nhood_id")
    tnk = m[m["class_group"] == "T/NK"]
    hits = (tnk["SpatialFDR"] < SFDR) & (tnk["logFC"] < 0)
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    ax.scatter(tnk.loc[~hits, "logFC"], tnk.loc[~hits, "logFC_cond"], s=12, c="#9e9e9e", label="other T/NK", alpha=0.8)
    ax.scatter(tnk.loc[hits, "logFC"], tnk.loc[hits, "logFC_cond"], s=18, c="#2166ac", label="TACSTD2 down, SpatialFDR<0.05")
    lims = [
        np.nanmin([tnk["logFC"].min(), tnk["logFC_cond"].min(), -0.5]),
        np.nanmax([tnk["logFC"].max(), tnk["logFC_cond"].max(), 0.5]),
    ]
    ax.plot(lims, lims, color="black", lw=0.8)
    ax.axhline(0, color="black", lw=0.4)
    ax.axvline(0, color="black", lw=0.4)
    ax.set_xlabel("TACSTD2 log2FC, unadjusted")
    ax.set_ylabel("TACSTD2 log2FC | CLDN4")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "tacstd2_lfc_conditional.png", dpi=150)
    fig.savefig(FIGS / "tacstd2_lfc_conditional.pdf")
    plt.close()


if __name__ == "__main__":
    main()
