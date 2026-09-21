#!/usr/bin/env python3
"""CLDN4 in HTAN lung single cells vs DNA-damage and STING scores.

Analyzed objects are the CELLxGENE copies of Synapse Level-4 h5ads
(anonymous download). Per-sample Synapse matrices that return HTTP 403
without a Synapse session are inventoried and not scored.

Inferential unit is the donor. Cell-level Spearman is saved as a
descriptive table only.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = Path("/tmp/htan_sting/data")
GENE_SETS = json.loads((Path(__file__).parent / "gene_sets.json").read_text())["sets"]
OUT_TAB = ROOT / "results" / "htan_cldn4_sting" / "tables"
OUT_FIG = ROOT / "results" / "htan_cldn4_sting" / "figures"
NOTES = ROOT / "notes" / "htan_cldn4_sting"

# Alias groups collapse to one column if both symbols exist.
ALIAS = {
    "H2AX": ("H2AX", "H2AFX"),
    "STING1": ("STING1", "TMEM173"),
}

SCORE_ORDER = [
    "HALLMARK_DNA_REPAIR",
    "DDR_DSB",
    "PROLIF",
    "STING_SENSOR",
    "STING_ISG",
    "REACTOME_STING",
    "CGAS",
    "STING1",
    "H2AX",
]

MIN_DONOR_CELLS = 40
MIN_SIDE = 10
MIN_SPEARMAN_N = 8


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    mask = np.isfinite(p)
    if mask.sum() == 0:
        return out
    pv = p[mask]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(n)
    tmp[order] = q
    out[mask] = tmp
    return out


def spearman_safe(x, y, min_n=MIN_SPEARMAN_N):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < min_n or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return n, np.nan, np.nan
    rho, p = stats.spearmanr(x[m], y[m])
    return n, float(rho), float(p)


def partial_spearman(x, y, z, min_n=MIN_SPEARMAN_N):
    """Pearson correlation of rank-residuals of x,y after z."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < min_n:
        return n, np.nan, np.nan
    rx = stats.rankdata(x[m])
    ry = stats.rankdata(y[m])
    rz = stats.rankdata(z[m])
    if np.std(rx) == 0 or np.std(ry) == 0 or np.std(rz) == 0:
        return n, np.nan, np.nan

    def resid(a, b):
        design = np.column_stack([np.ones(len(b)), b])
        coef, _, _, _ = np.linalg.lstsq(design, a, rcond=None)
        return a - design @ coef

    xr, yr = resid(rx, rz), resid(ry, rz)
    if np.std(xr) == 0 or np.std(yr) == 0:
        return n, np.nan, np.nan
    rho, p = stats.pearsonr(xr, yr)
    return n, float(rho), float(p)


def symbol_index(var):
    sym = var["feature_name"].astype(str).str.upper()
    mapping = {}
    for i, g in enumerate(sym):
        mapping.setdefault(g, i)
    return mapping


def canonical_genes(gene_list, mapping):
    """Return unique (canonical_name, index, matched_symbol) in list order."""
    found = []
    missing = []
    used = set()
    for g in gene_list:
        cands = ALIAS.get(g, (g,))
        hit = None
        for c in cands:
            if c.upper() in mapping:
                hit = c.upper()
                break
        if hit is None:
            missing.append(g)
            continue
        idx = mapping[hit]
        if idx in used:
            continue
        used.add(idx)
        found.append((g, idx, hit))
    return found, missing


def mean_z(mat):
    mat = np.asarray(mat, dtype=float)
    if mat.ndim == 1:
        mat = mat[:, None]
    mu = np.nanmean(mat, axis=0)
    sd = np.nanstd(mat, axis=0)
    sd = np.where(sd == 0, np.nan, sd)
    z = (mat - mu) / sd
    z = np.where(np.isfinite(z), z, 0.0)
    return z.mean(axis=1)


def load_matrix(path):
    a = ad.read_h5ad(path, backed="r")
    mapping = symbol_index(a.var)
    needed = set()
    for genes in GENE_SETS.values():
        for g in genes:
            for c in ALIAS.get(g, (g,)):
                if c.upper() in mapping:
                    needed.add(mapping[c.upper()])
                    break
    for g in ("CLDN4", "TACSTD2"):
        if g in mapping:
            needed.add(mapping[g])
    idx = sorted(needed)
    sub = a[:, idx].to_memory()
    obs = a.obs.copy()
    a.file.close()
    # Rebuild a dense symbol matrix in original idx order.
    X = sub.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)
    p99 = float(np.nanpercentile(X, 99))
    logged = False
    if p99 > 30:
        X = np.log1p(X)
        logged = True
    colnames = sub.var["feature_name"].astype(str).str.upper().tolist()
    expr = pd.DataFrame(X, index=obs.index, columns=colnames)
    # Collapse alias duplicate columns if both survived (should not).
    expr = expr.loc[:, ~expr.columns.duplicated()]
    return obs, expr, {"p99_before_optional_log1p": p99, "applied_log1p": logged, "n_obs": int(obs.shape[0]), "n_genes_pulled": int(expr.shape[1])}


def score_frame(expr, mask):
    """Z-scores are computed inside ``mask`` so immune cells do not set the scale."""
    sub = expr.loc[mask]
    coverage = {}
    scores = {}
    for name, genes in GENE_SETS.items():
        found, missing = canonical_genes(genes, {c: i for i, c in enumerate(sub.columns)})
        cols = [sub.columns[i] for _, i, _ in found]
        coverage[name] = {
            "n_requested": len(genes),
            "n_used": len(cols),
            "used": [m for _, _, m in found],
            "missing": missing,
        }
        if not cols:
            scores[name] = np.full(sub.shape[0], np.nan)
        else:
            scores[name] = mean_z(sub[cols].to_numpy())
    out = pd.DataFrame(scores, index=sub.index)
    if "CLDN4" not in sub.columns:
        raise RuntimeError("CLDN4 absent")
    out.insert(0, "CLDN4", sub["CLDN4"].to_numpy())
    return out, coverage


def donor_tests(df, stratum):
    """df columns: donor, CLDN4, score columns."""
    paired_rows = []
    within_rows = []
    for score in SCORE_ORDER:
        if score not in df.columns:
            continue
        deltas = []
        rhos = []
        n_used = 0
        n_pos = 0
        split_modes = set()
        for donor, g in df.groupby("donor"):
            if len(g) < MIN_DONOR_CELLS:
                continue
            x = g["CLDN4"].to_numpy(dtype=float)
            y = g[score].to_numpy(dtype=float)
            n, rho, p = spearman_safe(x, y, min_n=MIN_DONOR_CELLS)
            if np.isfinite(rho):
                rhos.append(rho)
                within_rows.append({
                    "stratum": stratum,
                    "score": score,
                    "donor": donor,
                    "n_cells": n,
                    "rho": rho,
                    "p": p,
                })
            med = float(np.median(x))
            if med == 0:
                hi = y[x > 0]
                lo = y[x == 0]
                mode = "gt0_vs_0"
            else:
                hi = y[x >= med]
                lo = y[x < med]
                mode = "median"
            split_modes.add(mode)
            if hi.size < MIN_SIDE or lo.size < MIN_SIDE:
                continue
            if np.nanstd(hi) == 0 and np.nanstd(lo) == 0:
                continue
            deltas.append(float(np.nanmean(hi) - np.nanmean(lo)))
            n_used += 1
            if deltas[-1] > 0:
                n_pos += 1
        deltas = np.asarray(deltas, dtype=float)
        rhos = np.asarray(rhos, dtype=float)
        if deltas.size >= 6 and np.any(deltas != 0):
            w = stats.wilcoxon(deltas, alternative="two-sided", zero_method="wilcox")
            wp, wstat = float(w.pvalue), float(w.statistic)
        else:
            wp, wstat = np.nan, np.nan
        if rhos.size >= 6 and np.any(rhos != 0):
            wr = stats.wilcoxon(rhos, alternative="two-sided", zero_method="wilcox")
            rp, rstat = float(wr.pvalue), float(wr.statistic)
        else:
            rp, rstat = np.nan, np.nan
        paired_rows.append({
            "stratum": stratum,
            "score": score,
            "n_donors_delta": int(deltas.size),
            "median_delta_hi_minus_lo": float(np.median(deltas)) if deltas.size else np.nan,
            "n_donors_delta_pos": int(n_pos),
            "wilcoxon_delta_stat": wstat,
            "wilcoxon_delta_p": wp,
            "split": ",".join(sorted(split_modes)) if split_modes else "",
            "n_donors_rho": int(rhos.size),
            "median_within_donor_rho": float(np.median(rhos)) if rhos.size else np.nan,
            "n_donors_rho_pos": int(np.sum(rhos > 0)) if rhos.size else 0,
            "wilcoxon_rho_stat": rstat,
            "wilcoxon_rho_p": rp,
        })
    # Within-donor partial Spearman: module vs CLDN4 after the proliferation score.
    if "PROLIF" in df.columns:
        for score in ("HALLMARK_DNA_REPAIR", "DDR_DSB", "STING_SENSOR", "STING_ISG"):
            if score not in df.columns:
                continue
            prhos = []
            for donor, g in df.groupby("donor"):
                if len(g) < MIN_DONOR_CELLS:
                    continue
                n, rho, _p = partial_spearman(
                    g["CLDN4"], g[score], g["PROLIF"], min_n=MIN_DONOR_CELLS
                )
                if np.isfinite(rho):
                    prhos.append(rho)
            prhos = np.asarray(prhos, dtype=float)
            if prhos.size >= 6 and np.any(prhos != 0):
                wr = stats.wilcoxon(prhos, alternative="two-sided", zero_method="wilcox")
                rp, rstat = float(wr.pvalue), float(wr.statistic)
            else:
                rp, rstat = np.nan, np.nan
            paired_rows.append({
                "stratum": stratum,
                "score": f"{score}_partial_PROLIF",
                "n_donors_delta": 0,
                "median_delta_hi_minus_lo": np.nan,
                "n_donors_delta_pos": 0,
                "wilcoxon_delta_stat": np.nan,
                "wilcoxon_delta_p": np.nan,
                "split": "partial_rank_residual",
                "n_donors_rho": int(prhos.size),
                "median_within_donor_rho": float(np.median(prhos)) if prhos.size else np.nan,
                "n_donors_rho_pos": int(np.sum(prhos > 0)) if prhos.size else 0,
                "wilcoxon_rho_stat": rstat,
                "wilcoxon_rho_p": rp,
            })
    # Donor-mean Spearman.
    means = df.groupby("donor").mean(numeric_only=True)
    # Keep donors that had enough cells.
    big = df.groupby("donor").size()
    means = means.loc[big[big >= MIN_DONOR_CELLS].index]
    mean_rows = []
    for score in SCORE_ORDER:
        if score not in means.columns:
            continue
        n, rho, p = spearman_safe(means["CLDN4"], means[score], min_n=MIN_SPEARMAN_N)
        mean_rows.append({
            "stratum": stratum,
            "score": score,
            "n_donors": n,
            "rho": rho,
            "p": p,
        })
    for score, control in (
        ("HALLMARK_DNA_REPAIR", "PROLIF"),
        ("DDR_DSB", "PROLIF"),
        ("STING_ISG", "PROLIF"),
        ("STING_SENSOR", "PROLIF"),
    ):
        if score not in means.columns or control not in means.columns:
            continue
        n, rho, p = partial_spearman(means["CLDN4"], means[score], means[control])
        mean_rows.append({
            "stratum": stratum,
            "score": f"{score}_partial_{control}",
            "n_donors": n,
            "rho": rho,
            "p": p,
        })
    return pd.DataFrame(paired_rows), pd.DataFrame(within_rows), pd.DataFrame(mean_rows), means


def cell_level(df, stratum):
    rows = []
    for score in SCORE_ORDER:
        if score not in df.columns:
            continue
        n, rho, p = spearman_safe(df["CLDN4"], df[score], min_n=MIN_SPEARMAN_N)
        rows.append({"stratum": stratum, "score": score, "n_cells": n, "rho": rho, "p": p, "unit": "cell_descriptive"})
    return pd.DataFrame(rows)


def stratum_masks(name, obs):
    if name == "chan":
        fine = obs["cell_type_fine"].astype(str)
        histo = obs["histo"].astype(str)
        return {
            "chan_sclc_malignant": fine.isin(["SCLC-A", "SCLC-N", "SCLC-P"]),
            "chan_nsclc_malignant": (fine == "NSCLC") & (histo == "LUAD"),
            "chan_nsclc_label_in_sclc": (fine == "NSCLC") & (histo == "SCLC"),
            "chan_normal_epithelial": fine.isin(
                ["AE1", "AEP", "Basal", "Ciliated", "Club", "Ionocyte", "Mucinous", "Tuft"]
            ) & (histo == "normal"),
        }, "donor_id"
    if name == "glasner":
        epi = obs["cell_lineage"].astype(str).eq("Epithelial")
        ttype = obs["Tissue Type"].astype(str)
        return {
            "glasner_luad_epithelium": epi,
            "glasner_met_epithelium": epi & ttype.eq("Metastasis"),
            "glasner_primary_epithelium": epi & ttype.eq("Primary"),
        }, "donor_id"
    raise KeyError(name)


def analyze_one(key, path):
    obs, expr, info = load_matrix(path)
    masks, donor_col = stratum_masks(key, obs)
    info["donor_col"] = donor_col
    info["file"] = str(path)
    cover_rows = []
    cell_rows = []
    paired_rows = []
    within_rows = []
    mean_rows = []
    desc_rows = []
    donor_mean_frames = []
    for stratum, mask in masks.items():
        mask = mask.fillna(False).to_numpy()
        n = int(mask.sum())
        donors = obs.loc[mask, donor_col].astype(str)
        desc_rows.append({
            "stratum": stratum,
            "n_cells": n,
            "n_donors": int(donors.nunique()),
            "cldn4_mean": float(expr.loc[mask, "CLDN4"].mean()) if n and "CLDN4" in expr.columns else np.nan,
            "cldn4_pct_gt0": float((expr.loc[mask, "CLDN4"] > 0).mean()) if n else np.nan,
            "cldn4_median": float(expr.loc[mask, "CLDN4"].median()) if n else np.nan,
        })
        if n < 50:
            continue
        scored, coverage = score_frame(expr, mask)
        for sname, cov in coverage.items():
            cover_rows.append({"stratum": stratum, "score": sname, **{k: cov[k] if not isinstance(cov[k], list) else ",".join(cov[k]) for k in cov}})
        scored = scored.copy()
        scored["donor"] = donors.to_numpy()
        # Single-gene columns already in GENE_SETS? No. Add from expr inside the mask using same z? 
        # Single genes are listed in SCORE_ORDER and also as sets of one via explicit columns.
        for gene in ("CGAS", "STING1", "H2AX"):
            cands = ALIAS.get(gene, (gene,))
            col = next((c for c in cands if c in expr.columns), None)
            if col is None:
                continue
            # Use the raw (optionally log1p) expression, not a z of one gene.
            # For Spearman this is monotone with a 1-gene z.
            scored[gene] = expr.loc[mask, col].to_numpy()
        cell_rows.append(cell_level(scored, stratum))
        paired, within, means, mean_df = donor_tests(scored, stratum)
        paired_rows.append(paired)
        within_rows.append(within)
        mean_rows.append(means)
        mean_df = mean_df.copy()
        mean_df["stratum"] = stratum
        mean_df["donor"] = mean_df.index
        donor_mean_frames.append(mean_df.reset_index(drop=True))
    return {
        "info": info,
        "desc": pd.DataFrame(desc_rows),
        "coverage": pd.DataFrame(cover_rows),
        "cell": pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame(),
        "paired": pd.concat(paired_rows, ignore_index=True) if paired_rows else pd.DataFrame(),
        "within": pd.concat(within_rows, ignore_index=True) if within_rows else pd.DataFrame(),
        "means": pd.concat(mean_rows, ignore_index=True) if mean_rows else pd.DataFrame(),
        "donor_means": pd.concat(donor_mean_frames, ignore_index=True) if donor_mean_frames else pd.DataFrame(),
    }


def add_fdr(df, pcol, family_mask):
    df = df.copy()
    df["q"] = np.nan
    if family_mask.any():
        df.loc[family_mask, "q"] = bh_fdr(df.loc[family_mask, pcol].to_numpy())
    return df


def plot_donor_means(donor_means):
    specs = [
        ("chan_sclc_malignant", "HALLMARK_DNA_REPAIR", "Chan SCLC malignant"),
        ("chan_sclc_malignant", "STING_ISG", "Chan SCLC malignant"),
        ("chan_nsclc_malignant", "HALLMARK_DNA_REPAIR", "Chan NSCLC (LUAD histo)"),
        ("chan_nsclc_malignant", "STING_ISG", "Chan NSCLC (LUAD histo)"),
        ("glasner_luad_epithelium", "HALLMARK_DNA_REPAIR", "Glasner LUAD epithelium"),
        ("glasner_luad_epithelium", "STING_ISG", "Glasner LUAD epithelium"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.6), constrained_layout=True)
    for ax, (stratum, score, title) in zip(axes.ravel(), specs):
        sub = donor_means[(donor_means["stratum"] == stratum)]
        if sub.empty or score not in sub.columns:
            ax.set_axis_off()
            continue
        x = sub["CLDN4"].to_numpy()
        y = sub[score].to_numpy()
        ax.scatter(x, y, s=28, c="#1f4e79", alpha=0.85, linewidths=0)
        n, rho, p = spearman_safe(x, y, min_n=5)
        ax.set_title(f"{title}\n{score}\nn={n}  ρ={rho:.2f}  p={p:.2g}" if np.isfinite(rho) else title, fontsize=9)
        ax.set_xlabel("donor-mean CLDN4")
        ax.set_ylabel(score.replace("_", " "), fontsize=8)
    fig.savefig(OUT_FIG / "donor_mean_cldn4_vs_scores.png", dpi=140)
    plt.close(fig)


def plot_median_rho(paired):
    keep_strata = [
        "chan_sclc_malignant",
        "chan_nsclc_malignant",
        "glasner_luad_epithelium",
    ]
    keep_scores = [
        "HALLMARK_DNA_REPAIR",
        "DDR_DSB",
        "PROLIF",
        "STING_SENSOR",
        "STING_ISG",
        "REACTOME_STING",
    ]
    sub = paired[paired["stratum"].isin(keep_strata) & paired["score"].isin(keep_scores)].copy()
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(8.4, 4.8), constrained_layout=True)
    ypos = np.arange(len(sub))
    colors = {
        "chan_sclc_malignant": "#1f4e79",
        "chan_nsclc_malignant": "#b85c38",
        "glasner_luad_epithelium": "#2f6f4e",
    }
    ax.axvline(0, color="#888", lw=0.8)
    ax.scatter(
        sub["median_within_donor_rho"],
        ypos,
        c=[colors.get(s, "#333") for s in sub["stratum"]],
        s=36,
        zorder=3,
    )
    labels = [f"{r.stratum} | {r.score} (n={int(r.n_donors_rho)})" for r in sub.itertuples()]
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("median within-donor Spearman ρ (CLDN4 vs score)")
    ax.set_title("Donor-level CLDN4 associations")
    fig.savefig(OUT_FIG / "median_within_donor_rho.png", dpi=140)
    plt.close(fig)


def main():
    OUT_TAB.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    NOTES.mkdir(parents=True, exist_ok=True)
    jobs = [
        ("chan", DATA / "chan_combined.h5ad"),
        ("glasner", DATA / "glasner_luad.h5ad"),
    ]
    parts = []
    for key, path in jobs:
        print(f"analyzing {key} {path}", flush=True)
        parts.append(analyze_one(key, path))
    desc = pd.concat([p["desc"] for p in parts], ignore_index=True)
    coverage = pd.concat([p["coverage"] for p in parts], ignore_index=True)
    cell = pd.concat([p["cell"] for p in parts], ignore_index=True)
    paired = pd.concat([p["paired"] for p in parts], ignore_index=True)
    within = pd.concat([p["within"] for p in parts], ignore_index=True)
    means = pd.concat([p["means"] for p in parts], ignore_index=True)
    donor_means = pd.concat([p["donor_means"] for p in parts], ignore_index=True)

    # Primary family: Wilcoxon on within-donor rho, author-malignant Chan strata
    # plus Glasner epithelium labeled as not CNV-called. Scores are the six modules.
    primary_scores = {
        "HALLMARK_DNA_REPAIR", "DDR_DSB", "PROLIF",
        "STING_SENSOR", "STING_ISG", "REACTOME_STING",
    }
    primary_strata = {
        "chan_sclc_malignant", "chan_nsclc_malignant", "glasner_luad_epithelium",
    }
    fam = paired["score"].isin(primary_scores) & paired["stratum"].isin(primary_strata)
    paired = add_fdr(paired, "wilcoxon_rho_p", fam)
    means_fam = means["score"].isin(primary_scores) & means["stratum"].isin(primary_strata)
    means = add_fdr(means, "p", means_fam)

    desc.to_csv(OUT_TAB / "stratum_descriptives.tsv", sep="\t", index=False)
    coverage.to_csv(OUT_TAB / "signature_coverage.tsv", sep="\t", index=False)
    cell.to_csv(OUT_TAB / "cell_level_spearman_descriptive.tsv", sep="\t", index=False)
    paired.to_csv(OUT_TAB / "donor_paired_and_within_rho.tsv", sep="\t", index=False)
    within.to_csv(OUT_TAB / "within_donor_spearman.tsv", sep="\t", index=False)
    means.to_csv(OUT_TAB / "donor_mean_spearman.tsv", sep="\t", index=False)
    donor_means.to_csv(OUT_TAB / "donor_means.tsv", sep="\t", index=False)
    info = {p["info"]["file"]: p["info"] for p in parts}
    (NOTES / "analysis_info.json").write_text(json.dumps(info, indent=2) + "\n")
    plot_donor_means(donor_means)
    plot_median_rho(paired)
    print(desc.to_string(index=False))
    print("--- paired primary ---")
    cols = ["stratum", "score", "n_donors_rho", "median_within_donor_rho", "wilcoxon_rho_p", "q",
            "n_donors_delta", "median_delta_hi_minus_lo", "wilcoxon_delta_p", "split"]
    show = paired.loc[fam, cols]
    print(show.to_string(index=False))
    print("--- donor means primary ---")
    print(means.loc[means_fam, ["stratum", "score", "n_donors", "rho", "p", "q"]].to_string(index=False))
    print("--- partial ---")
    print(means.loc[means["score"].astype(str).str.contains("partial"), ["stratum", "score", "n_donors", "rho", "p"]].to_string(index=False))


if __name__ == "__main__":
    main()
