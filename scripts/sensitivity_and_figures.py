#!/usr/bin/env python3
"""
B4 robustness sweep + immune-signature sanity check + figures.

Goal: honestly characterize how fragile the B4 result (NR > R for the TJ
signature in GSE126044) is, and whether the user's p=0.019 is reachable under
reasonable-but-different pre-processing / scoring choices. We do NOT tune to
hit 0.019; we enumerate a grid of defensible options and report the full range.
"""
import json
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_B3B4 import (
    DATA, OUT, TJ_CORE, TJ_KEGG, CD8, GEP, resolve_genes,
    tcga_t, tj_tcga, cd8_tcga, gep_tcga,
)

# ---------------------------------------------------------------------------
# GSE126044 load (repeat here so this script is standalone-ish)
# ---------------------------------------------------------------------------
gse = pd.read_csv(DATA / "GSE126044_counts.txt", sep="\t", index_col=0)
gse.index = gse.index.astype(str)
RESPONDERS = {"Dis_02", "Dis_15", "Dis_04", "Dis_17", "Dis_10"}
FFPE = {"Dis_06", "Dis_18", "Dis_09", "Dis_05", "Dis_08"}
labels = pd.Series(["R" if c in RESPONDERS else "NR" for c in gse.columns], index=gse.columns)

TJ_CLAUDIN = [g for g in TJ_CORE if g.startswith("CLDN")]
TJ_SETS = {"TJ_core": TJ_CORE, "TJ_kegg": TJ_KEGG, "TJ_claudin": TJ_CLAUDIN}


def norm_matrix(counts, how):
    if how == "logCPM":
        cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
        return np.log2(cpm + 1.0)
    if how == "logcounts":
        return np.log2(counts + 1.0)
    raise ValueError(how)


def score(expr, genes, method):
    matched, _ = resolve_genes(genes, expr.index)
    sub = expr.loc[matched]
    if method == "zmean":
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
        return z.mean(axis=0)
    if method == "meanexpr":
        return sub.mean(axis=0)
    if method == "rankmean":
        r = sub.rank(axis=1)
        r = (r - 1) / (r.shape[1] - 1)
        return r.mean(axis=0)
    raise ValueError(method)


def tests(nr, r):
    out = {}
    out["MWU_two"] = float(stats.mannwhitneyu(nr, r, alternative="two-sided")[1])
    out["MWU_NRgtR"] = float(stats.mannwhitneyu(nr, r, alternative="greater")[1])
    out["Welch_two"] = float(stats.ttest_ind(nr, r, equal_var=False)[1])
    out["Welch_NRgtR"] = float(stats.ttest_ind(nr, r, equal_var=False, alternative="greater")[1])
    return out


# ---------------------------------------------------------------------------
# B4 sensitivity grid
# ---------------------------------------------------------------------------
rows = []
for tjname, sample_set, norm, method in itertools.product(
    TJ_SETS, ["all16", "fresh_only"], ["logCPM", "logcounts"], ["zmean", "meanexpr", "rankmean"]
):
    cols = gse.columns if sample_set == "all16" else [c for c in gse.columns if c not in FFPE]
    counts = gse[cols]
    lab = labels[cols]
    expr = norm_matrix(counts, norm)
    s = score(expr, TJ_SETS[tjname], method)
    nr, r = s[lab == "NR"].values, s[lab == "R"].values
    row = {
        "TJ_set": tjname, "samples": sample_set, "norm": norm, "method": method,
        "n_NR": int((lab == "NR").sum()), "n_R": int((lab == "R").sum()),
        "median_NR": float(np.median(nr)), "median_R": float(np.median(r)),
        "direction": "NR>R" if np.median(nr) > np.median(r) else "NR<=R",
    }
    row.update(tests(nr, r))
    rows.append(row)

sens = pd.DataFrame(rows)
sens.to_csv(OUT / "B4_sensitivity_grid.csv", index=False)

# ---------------------------------------------------------------------------
# Sanity: do canonical immune signatures separate R/NR in this cohort?
# (uses the primary pre-specified config: logCPM + zmean, all 16)
# ---------------------------------------------------------------------------
expr_main = norm_matrix(gse, "logCPM")
sanity = {}
for name, genes in [("CD8", CD8), ("GEP", GEP), ("TJ_core", TJ_CORE)]:
    s = score(expr_main, genes, "zmean")
    nr, r = s[labels == "NR"].values, s[labels == "R"].values
    sanity[name] = {
        "median_NR": float(np.median(nr)), "median_R": float(np.median(r)),
        "direction_R_vs_NR": "R>NR" if np.median(r) > np.median(nr) else "R<=NR",
        **tests(nr, r),
    }

summary = {
    "B4_sensitivity": {
        "n_configs": len(sens),
        "n_direction_NR_gt_R": int((sens["direction"] == "NR>R").sum()),
        "MWU_two_min": float(sens["MWU_two"].min()),
        "MWU_two_median": float(sens["MWU_two"].median()),
        "MWU_two_max": float(sens["MWU_two"].max()),
        "MWU_NRgtR_min": float(sens["MWU_NRgtR"].min()),
        "any_config_p_two_below_0.05": bool((sens["MWU_two"] < 0.05).any()),
        "any_config_p_one_below_0.05": bool((sens["MWU_NRgtR"] < 0.05).any()),
        "configs_near_0.019 (|p-0.019|<0.01, one-sided)": sens.loc[
            (sens["MWU_NRgtR"] - 0.019).abs() < 0.01,
            ["TJ_set", "samples", "norm", "method", "MWU_NRgtR"],
        ].to_dict("records"),
    },
    "GSE126044_immune_sanity_check": sanity,
}
with open(OUT / "B4_sensitivity_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
# B3 scatter (TCGA-LUAD, TJ_core vs CD8 and GEP)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, (yname, yvals) in zip(axes, [("CD8 signature", cd8_tcga), ("GEP signature", gep_tcga)]):
    rho, p = stats.spearmanr(tj_tcga, yvals)
    ax.scatter(tj_tcga, yvals, s=10, alpha=0.5, edgecolor="none", color="#2b6cb0")
    m, b = np.polyfit(tj_tcga, yvals, 1)
    xs = np.linspace(tj_tcga.min(), tj_tcga.max(), 50)
    ax.plot(xs, m * xs + b, color="#c53030", lw=2)
    ax.set_xlabel("TJ_core signature (z-mean)")
    ax.set_ylabel(yname + " (z-mean)")
    ax.set_title(f"TCGA-LUAD (n={len(tj_tcga)})\nSpearman rho={rho:.2f}, p={p:.1e}")
fig.suptitle("Claim B3: Tight-junction signature vs CD8 / GEP in TCGA-LUAD", fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "B3_TCGA_LUAD_scatter.png", dpi=150)
plt.close(fig)

# B4 boxplot (GSE126044, TJ_core, primary config)
s_tj = score(expr_main, TJ_CORE, "zmean")
nr, r = s_tj[labels == "NR"].values, s_tj[labels == "R"].values
p_two = stats.mannwhitneyu(nr, r, alternative="two-sided")[1]
p_one = stats.mannwhitneyu(nr, r, alternative="greater")[1]
fig, ax = plt.subplots(figsize=(5.2, 5))
ax.boxplot([r, nr], tick_labels=[f"R (n={len(r)})", f"NR (n={len(nr)})"], widths=0.5,
           showfliers=False)
for i, (vals, col) in enumerate([(r, "#2f855a"), (nr, "#c53030")], start=1):
    ax.scatter(np.random.normal(i, 0.05, len(vals)), vals, color=col, s=35, zorder=3, alpha=0.8)
ax.set_ylabel("TJ_core signature (z-mean)")
ax.set_title("Claim B4: GSE126044 TJ signature by response\n"
             f"MWU two-sided p={p_two:.3f}; one-sided NR>R p={p_one:.3f}")
fig.tight_layout()
fig.savefig(OUT / "B4_GSE126044_boxplot.png", dpi=150)
plt.close(fig)
print("figures + tables written to", OUT)
