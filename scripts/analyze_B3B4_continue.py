#!/usr/bin/env python3
"""
Claims B3 + B4 — continuation, aligned to the official claim pages.

Source of the claims
--------------------
claims/B3.md, claims/B4.md (user PPT 2026-08-17):
  B3  TCGA-LUAD: TJ-high tumors have low CD8 / low GEP, p < 1e-6.
      Test plan: score TJ by ssGSEA/GSVA (KEGG tight junction / CLDN1/4/7,
      F11R, TJP1/2, OCLN); CD8 = CD8A/CD8B (+ cytolytic); GEP = Ayers 18-gene;
      Spearman + TJ-high vs TJ-low Wilcoxon.
  B4  GSE126044: NR have higher TJ than R, p = 0.019.
      Same TJ definition as B3; Wilcoxon/t-test; also check CLDN4 alone.

This script does NOT invent a gene list to hit the claimed p-values. It uses
the lists named in the claim pages plus the 15-gene structural set already
used in this repo (PR #69), scores them two ways (z-mean and ssGSEA), and
reports the actual numbers.

Primary TJ set: the 15-gene structural signature from PR #69
  CLDN1/3/4/7, OCLN, TJP1/2/3, F11R, JAM2/3, MARVELD2/3, CGN, CGNL1.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
OUT = REPO / "results" / "claim_B3B4"
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260816)

# ---------------------------------------------------------------------------
# Gene sets
# ---------------------------------------------------------------------------
# Official 15-gene structural TJ used in this repo (PR #69). PRIMARY.
TJ_15 = [
    "CLDN1", "CLDN3", "CLDN4", "CLDN7",
    "OCLN",
    "TJP1", "TJP2", "TJP3",
    "F11R", "JAM2", "JAM3",
    "MARVELD2", "MARVELD3",
    "CGN", "CGNL1",
]

# Minimal module named on the B3 claim page.
TJ_CLAIM = ["CLDN1", "CLDN4", "CLDN7", "F11R", "TJP1", "TJP2", "OCLN"]

CD8 = ["CD8A", "CD8B"]
CYT = ["GZMA", "PRF1"]  # Rooney et al. cytolytic score
GEP = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]

SYNONYMS = {
    "PDCD1LG2": ["PDL2"],
    "F11R": ["JAM1", "JAMA"],
    "HLA-DQA1": ["HLA.DQA1"],
    "HLA-DRB1": ["HLA.DRB1"],
    "HLA-E": ["HLA.E"],
}


def msigdb_symbols(path: Path, key: str) -> list[str]:
    obj = json.loads(path.read_text())
    return list(obj[key]["geneSymbols"])


KEGG_TJ = msigdb_symbols(DATA / "KEGG_TIGHT_JUNCTION.json", "KEGG_TIGHT_JUNCTION")
GOCC_TJ = msigdb_symbols(DATA / "GOCC_TIGHT_JUNCTION.json", "GOCC_TIGHT_JUNCTION")

TJ_SETS = {
    "TJ_15": TJ_15,
    "TJ_claim_module": TJ_CLAIM,
    "KEGG_TIGHT_JUNCTION": KEGG_TJ,
    "GOCC_TIGHT_JUNCTION": GOCC_TJ,
}

# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------
def resolve(wanted, available) -> tuple[list[str], list[str]]:
    avail = set(available)
    matched, missing = [], []
    for g in wanted:
        hit = next((c for c in [g] + SYNONYMS.get(g, []) if c in avail), None)
        if hit is None:
            missing.append(g)
        else:
            matched.append(hit)
    return matched, missing


def zmean(expr: pd.DataFrame, genes) -> tuple[pd.Series, list, list]:
    matched, missing = resolve(genes, expr.index)
    sub = expr.loc[matched]
    z = sub.sub(sub.mean(axis=1), axis=0).div(
        sub.std(axis=1, ddof=0).replace(0, np.nan), axis=0
    )
    return z.mean(axis=0), matched, missing


def ssgsea(expr: pd.DataFrame, genes, tau: float = 0.25) -> tuple[pd.Series, list, list]:
    """GSVA-style ssGSEA (Barbie 2009): area under the KS-like running sum.

    For each sample, genes are ranked high-expression → rank 1. Set members
    contribute (rank^tau) / sum(rank^tau); non-members contribute −1/(N−n).
    The score is the sum of the cumulative walk (not just the maximum).
    Rank-based tests (Spearman, Mann-Whitney) are invariant to a subsequent
    per-set min-max normalization, so we leave scores unnormalized.
    """
    matched, missing = resolve(genes, expr.index)
    ranks = expr.rank(axis=0, method="average", ascending=False)
    N = ranks.shape[0]
    in_set = np.asarray(ranks.index.isin(matched))
    n = int(in_set.sum())
    if n == 0 or n == N:
        return pd.Series(np.nan, index=expr.columns), matched, missing
    scores = {}
    in_set_vals = in_set
    for col in ranks.columns:
        r = ranks[col].to_numpy()
        order = np.argsort(r, kind="mergesort")
        in_ord = in_set_vals[order]
        r_ord = r[order]
        r_alpha = np.power(r_ord, tau)
        denom = r_alpha[in_ord].sum()
        if denom == 0:
            scores[col] = np.nan
            continue
        step = np.where(in_ord, r_alpha / denom, -1.0 / (N - n))
        scores[col] = float(np.cumsum(step).sum())
    return pd.Series(scores), matched, missing


def spearman_ci(x, y, n_boot=2000):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    rho, p = stats.spearmanr(x, y)
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = RNG.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    lo, hi = np.quantile(boots, [0.025, 0.975])
    return {"rho": float(rho), "p": float(p), "n": int(n),
            "ci95_lo": float(lo), "ci95_hi": float(hi)}


def partial_spearman(x, y, z):
    """Rank-residual partial Spearman: Spearman of (x|z) vs (y|z) residuals."""
    xr = stats.rankdata(x)
    yr = stats.rankdata(y)
    zr = stats.rankdata(z)
    bx = np.polyfit(zr, xr, 1)
    by = np.polyfit(zr, yr, 1)
    rx = xr - (bx[0] * zr + bx[1])
    ry = yr - (by[0] * zr + by[1])
    return spearman_ci(rx, ry)


def high_low(score, outcome):
    """Median-split of `score`; Wilcoxon of `outcome` in high vs low."""
    med = np.median(score)
    hi = outcome[score > med]
    lo = outcome[score <= med]
    # expect TJ-high → lower CD8/GEP, so hi < lo
    u, p_two = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    _, p_less = stats.mannwhitneyu(hi, lo, alternative="less")
    return {
        "n_high": int(len(hi)), "n_low": int(len(lo)),
        "median_outcome_TJhigh": float(np.median(hi)),
        "median_outcome_TJlow": float(np.median(lo)),
        "direction": "TJhigh_lower" if np.median(hi) < np.median(lo) else "TJhigh_higher_or_eq",
        "MWU_p_two": float(p_two),
        "MWU_p_TJhigh_lower": float(p_less),
    }


def group_nr_r(score, labels):
    nr = score[labels == "NR"].to_numpy()
    r = score[labels == "R"].to_numpy()
    return {
        "n_NR": int(len(nr)), "n_R": int(len(r)),
        "median_NR": float(np.median(nr)), "median_R": float(np.median(r)),
        "mean_NR": float(np.mean(nr)), "mean_R": float(np.mean(r)),
        "direction": "NR>R" if np.median(nr) > np.median(r) else "NR<=R",
        "MWU_p_two": float(stats.mannwhitneyu(nr, r, alternative="two-sided")[1]),
        "MWU_p_NRgtR": float(stats.mannwhitneyu(nr, r, alternative="greater")[1]),
        "Welch_p_two": float(stats.ttest_ind(nr, r, equal_var=False)[1]),
        "Welch_p_NRgtR": float(stats.ttest_ind(nr, r, equal_var=False, alternative="greater")[1]),
    }


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("loading TCGA-LUAD...")
tcga = pd.read_csv(DATA / "TCGA_LUAD_HiSeqV2", sep="\t", index_col=0)
tcga.index = tcga.index.astype(str)
sample_type = pd.Series(tcga.columns).str.split("-").str[3].str[:2].to_numpy()
tcga_t = tcga.loc[:, tcga.columns[sample_type == "01"]]

print("loading Aran ESTIMATE purity...")
purity_raw = pd.read_excel(DATA / "aran_purity.xlsx", header=3)
purity_raw = purity_raw.rename(columns={"Sample ID": "sample", "Cancer type": "cancer",
                                        "ESTIMATE": "ESTIMATE_purity", "CPE": "CPE"})
luad_p = purity_raw.loc[purity_raw["cancer"] == "LUAD", ["sample", "ESTIMATE_purity", "CPE"]].copy()
# Xena barcodes are TCGA-XX-XXXX-01 (no vial letter). Aran uses -01A.
luad_p["barcode15"] = luad_p["sample"].astype(str).str[:15]
purity_map = luad_p.dropna(subset=["ESTIMATE_purity"]).drop_duplicates("barcode15").set_index("barcode15")["ESTIMATE_purity"]

print("loading GSE126044...")
gse = pd.read_csv(DATA / "GSE126044_counts.txt", sep="\t", index_col=0)
gse.index = gse.index.astype(str)
logcpm = np.log2(gse.div(gse.sum(axis=0), axis=1) * 1e6 + 1.0)
RESPONDERS = {"Dis_02", "Dis_15", "Dis_04", "Dis_17", "Dis_10"}
labels = pd.Series(["R" if c in RESPONDERS else "NR" for c in gse.columns], index=gse.columns)

# ---------------------------------------------------------------------------
# Score everything
# ---------------------------------------------------------------------------
print("scoring TCGA...")
tcga_scores = {}
tcga_meta = {}
for name, genes in TJ_SETS.items():
    for method, fn in [("zmean", zmean), ("ssgsea", ssgsea)]:
        s, matched, missing = fn(tcga_t, genes)
        key = f"{name}__{method}"
        tcga_scores[key] = s
        tcga_meta[key] = {"n_requested": len(genes), "n_matched": len(matched),
                          "missing": missing}

cd8_z, cd8_m, cd8_miss = zmean(tcga_t, CD8)
cyt_z, cyt_m, cyt_miss = zmean(tcga_t, CYT)
gep_z, gep_m, gep_miss = zmean(tcga_t, GEP)
cd8_ss, _, _ = ssgsea(tcga_t, CD8)
gep_ss, _, _ = ssgsea(tcga_t, GEP)

# ---------------------------------------------------------------------------
# B3
# ---------------------------------------------------------------------------
print("B3 correlations...")
B3 = {"n_primary_tumor": int(tcga_t.shape[1]), "signatures": tcga_meta,
      "immune_genes": {
          "CD8": {"matched": cd8_m, "missing": cd8_miss},
          "CYT": {"matched": cyt_m, "missing": cyt_miss},
          "GEP": {"matched": gep_m, "missing": gep_miss},
      },
      "correlations": {}, "high_low": {}, "purity_partial": {}}

immune = {"CD8_zmean": cd8_z, "CYT_zmean": cyt_z, "GEP_zmean": gep_z,
          "CD8_ssgsea": cd8_ss, "GEP_ssgsea": gep_ss}

for tj_key, tj_s in tcga_scores.items():
    B3["correlations"][tj_key] = {iname: spearman_ci(tj_s, iv) for iname, iv in immune.items()}
    B3["high_low"][tj_key] = {iname: high_low(tj_s.to_numpy(), iv.to_numpy())
                              for iname, iv in immune.items()}

# purity-adjusted (intersect samples with Aran ESTIMATE)
common = [c for c in tcga_t.columns if c in purity_map.index]
pur = purity_map.loc[common]
B3["purity"] = {
    "source": "Aran et al. Nat Commun 2015 Supp Data 1, ESTIMATE column",
    "n_with_purity": int(len(common)),
}
for tj_key, tj_s in tcga_scores.items():
    x = tj_s.loc[common]
    B3["purity_partial"][tj_key] = {
        "vs_CD8_zmean": partial_spearman(x, cd8_z.loc[common], pur),
        "vs_GEP_zmean": partial_spearman(x, gep_z.loc[common], pur),
        "vs_CYT_zmean": partial_spearman(x, cyt_z.loc[common], pur),
        "TJ_vs_purity": spearman_ci(x, pur),
    }
B3["purity_partial"]["CD8_vs_purity"] = spearman_ci(cd8_z.loc[common], pur)
B3["purity_partial"]["GEP_vs_purity"] = spearman_ci(gep_z.loc[common], pur)

# ---------------------------------------------------------------------------
# B4
# ---------------------------------------------------------------------------
print("scoring GSE126044...")
B4 = {"n_R": int((labels == "R").sum()), "n_NR": int((labels == "NR").sum()),
      "tests": {}, "CLDN4_alone": {}, "immune_sanity": {}}

gse_scores = {}
for name, genes in TJ_SETS.items():
    for method, fn in [("zmean", zmean), ("ssgsea", ssgsea)]:
        s, matched, missing = fn(logcpm, genes)
        key = f"{name}__{method}"
        gse_scores[key] = s
        B4["tests"][key] = {**group_nr_r(s, labels),
                            "n_matched": len(matched), "missing": missing}

# CLDN4 alone (claim page explicitly asks)
if "CLDN4" in logcpm.index:
    B4["CLDN4_alone"]["logCPM"] = group_nr_r(logcpm.loc["CLDN4"], labels)
    # also z-scored single gene (same ranks as logCPM → identical MWU)
    B4["CLDN4_alone"]["note"] = "single-gene MWU is invariant to z-score / rank transform of CLDN4"

# immune sanity: do CD8/GEP separate R vs NR?
for iname, genes in [("CD8", CD8), ("CYT", CYT), ("GEP", GEP)]:
    s, _, _ = zmean(logcpm, genes)
    B4["immune_sanity"][iname] = group_nr_r(s, labels)

# ---------------------------------------------------------------------------
# Persist
# ---------------------------------------------------------------------------
payload = {
    "B3": B3,
    "B4": B4,
    "methods": {
        "TJ_15": TJ_15,
        "TJ_claim_module": TJ_CLAIM,
        "KEGG_TIGHT_JUNCTION_n": len(KEGG_TJ),
        "GOCC_TIGHT_JUNCTION_n": len(GOCC_TJ),
        "scoring": ["zmean = mean per-gene z-score",
                    "ssgsea = GSVA-style Barbie 2009 area-under-walk, tau=0.25, unnormalized"],
        "B3_expression": "UCSC Xena TCGA.LUAD.sampleMap/HiSeqV2 log2(norm_count+1); primary tumors only",
        "B4_expression": "GEO GSE126044 raw counts → log2(CPM+1); labels from series matrix",
        "purity": "Aran et al. 2015 ESTIMATE purity; partial Spearman via rank residuals",
    },
}
# drop huge nested matched lists from the main JSON? keep missing only — already done.

with open(OUT / "results_B3B4_continue.json", "w") as f:
    json.dump(payload, f, indent=2)

# per-sample tables
tcga_tbl = pd.DataFrame(tcga_scores)
tcga_tbl["CD8_zmean"] = cd8_z
tcga_tbl["CYT_zmean"] = cyt_z
tcga_tbl["GEP_zmean"] = gep_z
tcga_tbl["ESTIMATE_purity"] = purity_map.reindex(tcga_tbl.index)
tcga_tbl.to_csv(OUT / "TCGA_LUAD_scores_continue.csv")

gse_tbl = pd.DataFrame(gse_scores)
gse_tbl["response"] = labels
gse_tbl["CLDN4_logCPM"] = logcpm.loc["CLDN4"] if "CLDN4" in logcpm.index else np.nan
gse_tbl.to_csv(OUT / "GSE126044_scores_continue.csv")

# compact summary table for README
rows = []
for tj_key, block in B3["correlations"].items():
    for iname, st in block.items():
        rows.append({"claim": "B3", "TJ": tj_key, "vs": iname,
                     "rho": st["rho"], "p": st["p"],
                     "ci95_lo": st["ci95_lo"], "ci95_hi": st["ci95_hi"]})
pd.DataFrame(rows).to_csv(OUT / "B3_correlation_table.csv", index=False)

rows = []
for tj_key, st in B4["tests"].items():
    rows.append({"TJ": tj_key, **{k: st[k] for k in
                 ["n_NR", "n_R", "median_NR", "median_R", "direction",
                  "MWU_p_two", "MWU_p_NRgtR", "Welch_p_two", "Welch_p_NRgtR"]}})
if B4["CLDN4_alone"]:
    st = B4["CLDN4_alone"]["logCPM"]
    rows.append({"TJ": "CLDN4_alone", **st})
pd.DataFrame(rows).to_csv(OUT / "B4_group_table.csv", index=False)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
print("figures...")

def scatter(ax, x, y, title, xlabel, ylabel):
    st = spearman_ci(x, y, n_boot=400)
    ax.scatter(x, y, s=10, alpha=0.45, edgecolor="none", color="#2b6cb0")
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(np.min(x), np.max(x), 50)
    ax.plot(xs, m * xs + b, color="#c53030", lw=2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\nSpearman ρ={st['rho']:.2f}, p={st['p']:.1e}")


fig, axes = plt.subplots(2, 2, figsize=(11, 9.2))
pairs = [
    (axes[0, 0], tcga_scores["TJ_15__zmean"], cd8_z, "TJ_15 z-mean vs CD8", "TJ_15 (z-mean)", "CD8 (z-mean)"),
    (axes[0, 1], tcga_scores["TJ_15__zmean"], gep_z, "TJ_15 z-mean vs GEP", "TJ_15 (z-mean)", "GEP (z-mean)"),
    (axes[1, 0], tcga_scores["TJ_15__ssgsea"], cd8_z, "TJ_15 ssGSEA vs CD8", "TJ_15 (ssGSEA)", "CD8 (z-mean)"),
    (axes[1, 1], tcga_scores["KEGG_TIGHT_JUNCTION__ssgsea"], cd8_z,
     "KEGG_TIGHT_JUNCTION ssGSEA vs CD8", "KEGG TJ (ssGSEA)", "CD8 (z-mean)"),
]
for ax, x, y, t, xl, yl in pairs:
    scatter(ax, x.to_numpy(), y.to_numpy(), t, xl, yl)
fig.suptitle("Claim B3 (continued): TCGA-LUAD TJ vs CD8 / GEP — official 15-gene + ssGSEA",
             fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "B3_continue_scatter.png", dpi=150)
plt.close(fig)

# purity residual scatter for primary TJ_15 z-mean vs CD8
x = tcga_scores["TJ_15__zmean"].loc[common]
y = cd8_z.loc[common]
xr = stats.rankdata(x)
yr = stats.rankdata(y)
zr = stats.rankdata(pur)
rx = xr - np.polyval(np.polyfit(zr, xr, 1), zr)
ry = yr - np.polyval(np.polyfit(zr, yr, 1), zr)
st = spearman_ci(rx, ry, n_boot=400)
fig, ax = plt.subplots(figsize=(5.6, 5.2))
ax.scatter(rx, ry, s=10, alpha=0.45, edgecolor="none", color="#2b6cb0")
m, b = np.polyfit(rx, ry, 1)
xs = np.linspace(rx.min(), rx.max(), 50)
ax.plot(xs, m * xs + b, color="#c53030", lw=2)
ax.set_xlabel("TJ_15 z-mean rank residual | ESTIMATE purity")
ax.set_ylabel("CD8 z-mean rank residual | ESTIMATE purity")
ax.set_title(f"B3 purity-adjusted (n={len(common)})\npartial Spearman ρ={st['rho']:.2f}, p={st['p']:.1e}")
fig.tight_layout()
fig.savefig(OUT / "B3_continue_purity_partial.png", dpi=150)
plt.close(fig)

# B4 boxplots: claim-page 7-gene module (the one that hits p=0.019), plus alternatives
fig, axes = plt.subplots(1, 4, figsize=(14, 4.6))
panels = [
    ("TJ_claim 7-gene z-mean", gse_scores["TJ_claim_module__zmean"]),
    ("TJ_15 z-mean", gse_scores["TJ_15__zmean"]),
    ("KEGG TJ ssGSEA", gse_scores["KEGG_TIGHT_JUNCTION__ssgsea"]),
    ("CLDN4 logCPM", logcpm.loc["CLDN4"]),
]
for ax, (title, s) in zip(axes, panels):
    nr = s[labels == "NR"].to_numpy()
    r = s[labels == "R"].to_numpy()
    p_two = stats.mannwhitneyu(nr, r, alternative="two-sided")[1]
    p_one = stats.mannwhitneyu(nr, r, alternative="greater")[1]
    ax.boxplot([r, nr], tick_labels=[f"R\nn={len(r)}", f"NR\nn={len(nr)}"],
               widths=0.55, showfliers=False)
    ax.scatter(np.full(len(r), 1) + RNG.normal(0, 0.04, len(r)), r,
               color="#2f855a", s=28, zorder=3)
    ax.scatter(np.full(len(nr), 2) + RNG.normal(0, 0.04, len(nr)), nr,
               color="#c53030", s=28, zorder=3)
    ax.set_title(f"{title}\n2s p={p_two:.3f}; NR>R p={p_one:.3f}")
fig.suptitle("Claim B4 (continued): GSE126044 NR vs R — official TJ defs + CLDN4",
             fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "B4_continue_boxplots.png", dpi=150)
plt.close(fig)

# print a compact human summary
print("\n======== B3 (Spearman, selected) ========")
for key in ["TJ_15__zmean", "TJ_15__ssgsea", "TJ_claim_module__zmean",
            "KEGG_TIGHT_JUNCTION__ssgsea", "KEGG_TIGHT_JUNCTION__zmean"]:
    c = B3["correlations"][key]
    print(f"{key:32s}  CD8 ρ={c['CD8_zmean']['rho']:+.3f} p={c['CD8_zmean']['p']:.2e}   "
          f"GEP ρ={c['GEP_zmean']['rho']:+.3f} p={c['GEP_zmean']['p']:.2e}")
print("\n======== B3 purity-partial TJ_15 zmean ========")
pp = B3["purity_partial"]["TJ_15__zmean"]
print("vs CD8", pp["vs_CD8_zmean"])
print("vs GEP", pp["vs_GEP_zmean"])
print("TJ vs purity", pp["TJ_vs_purity"])
print("CD8 vs purity", B3["purity_partial"]["CD8_vs_purity"])
print("\n======== B4 ========")
for key, st in B4["tests"].items():
    print(f"{key:36s} dir={st['direction']:6s}  MWU2={st['MWU_p_two']:.3f}  "
          f"NRgtR={st['MWU_p_NRgtR']:.3f}")
print("CLDN4", B4["CLDN4_alone"].get("logCPM"))
print("immune sanity", {k: (v["direction"], v["MWU_p_two"]) for k, v in B4["immune_sanity"].items()})
print("wrote", OUT)
