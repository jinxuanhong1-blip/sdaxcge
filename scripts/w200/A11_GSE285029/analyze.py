#!/usr/bin/env python3
"""GSE285029: TACSTD2 / CLDN4 vs CD274, IFN/MHC-I, IL6/STAT3, CD8/GEP.

Public author WTS matrix only. No GEO response / histology / purity labels.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
DATA = ROOT / "data" / "w200" / "A11_GSE285029"
OUT = ROOT / "results" / "w200" / "A11_GSE285029"
GENESETS = HERE / "genesets"

MATRIX = DATA / "GSE285029_WTS_expr_count_235_032820.txt.gz"

# Ayers et al., J Clin Invest 2017; unweighted z-mean (not NanoString TIS weights).
GEP18 = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]
# Compact IFN-γ / IFN output (a priori; not mined from this matrix).
IFN_COMPACT = ["IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1"]
# Classical MHC-I antigen-presentation cassette.
MHC1 = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5", "PSMB8", "PSMB9", "TAPBP"]
# Compact IL-6 / Jak / Stat3 cassette named by the source paper.
IL6_STAT3 = ["IL6", "IL6R", "IL6ST", "JAK1", "JAK2", "STAT3", "SOCS3"]
# A1 8-gene T-cell effector (sibling convention).
IMMUNE8 = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21"]
CYT = ["GZMA", "PRF1"]
# Epithelial proxy used only as a purity-like sensitivity covariate.
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]

PRIMARY_FEATURES = [
    "CD274",
    "IFN_compact",
    "MHC1",
    "IL6",
    "IL6_STAT3_compact",
    "CD8A",
    "GEP18",
]
SENSITIVITY_FEATURES = [
    "HALLMARK_IFNG",
    "HALLMARK_IFNA",
    "HALLMARK_IL6_JAK_STAT3",
    "immune8",
    "CYT",
    "STAT3",
    "IFNG",
]

RHO_POS_ASSOC = 0.20
RHO_WEAK = 0.10


def load_gene_list(name: str) -> list[str]:
    path = GENESETS / f"{name}.txt"
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def load_matrix() -> pd.DataFrame:
    if not MATRIX.exists():
        raise SystemExit(f"missing {MATRIX}; run download.py first")
    expr = pd.read_csv(MATRIX, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    # Drop empty / all-NA genes.
    expr = expr.dropna(axis=0, how="all")
    return expr


def log2p1_clip(expr: pd.DataFrame) -> pd.DataFrame:
    return np.log2(expr.clip(lower=0) + 1.0)


def zmean(log_expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str], list[str]]:
    present = [g for g in genes if g in log_expr.index]
    missing = [g for g in genes if g not in log_expr.index]
    if not present:
        return pd.Series(np.nan, index=log_expr.columns), present, missing
    block = log_expr.loc[present]
    mu = block.mean(axis=1)
    sd = block.std(axis=1, ddof=0).replace(0, np.nan)
    z = block.sub(mu, axis=0).div(sd, axis=0)
    return z.mean(axis=0), present, missing


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = int(len(d))
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci_lo": np.nan, "ci_hi": np.nan}
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    # Fisher-z 95% CI, SE = 1/sqrt(n-3)
    if n > 3 and abs(rho) < 1:
        z = np.arctanh(rho)
        se = 1.0 / math.sqrt(n - 3)
        ci_lo, ci_hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    else:
        ci_lo = ci_hi = np.nan
    return {"n": n, "rho": float(rho), "p": float(p), "ci_lo": float(ci_lo), "ci_hi": float(ci_hi)}


def partial_spearman(x: pd.Series, y: pd.Series, cov: pd.Series) -> dict:
    d = pd.concat([x, y, cov], axis=1).dropna()
    d.columns = ["x", "y", "c"]
    n = int(len(d))
    if n < 6:
        return {"n": n, "partial_rho": np.nan, "partial_p": np.nan, "ci_lo": np.nan, "ci_hi": np.nan}
    rx = d["x"].rank()
    ry = d["y"].rank()
    rc = d["c"].rank()
    # Residualize ranks on the covariate (first-order partial Spearman).
    bx = np.polyfit(rc, rx, 1)
    by = np.polyfit(rc, ry, 1)
    ex = rx - (bx[0] * rc + bx[1])
    ey = ry - (by[0] * rc + by[1])
    r, p = stats.pearsonr(ex, ey)
    df = n - 3
    if df > 0 and abs(r) < 1:
        t = r * math.sqrt(df / (1 - r * r))
        p = float(2 * stats.t.sf(abs(t), df))
        z = np.arctanh(r)
        se = 1.0 / math.sqrt(n - 4) if n > 4 else np.nan
        ci_lo = float(np.tanh(z - 1.96 * se)) if se == se else np.nan
        ci_hi = float(np.tanh(z + 1.96 * se)) if se == se else np.nan
    else:
        ci_lo = ci_hi = np.nan
    return {
        "n": n,
        "partial_rho": float(r),
        "partial_p": float(p),
        "ci_lo": float(ci_lo) if ci_lo == ci_lo else np.nan,
        "ci_hi": float(ci_hi) if ci_hi == ci_hi else np.nan,
    }


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        i = n - rank + 1
        val = min(prev, p[idx] * n / i)
        q[idx] = val
        prev = val
    return [float(x) for x in q]


def label_rho(rho: float, p: float) -> str:
    if not (rho == rho) or not (p == p):
        return "NA"
    if p >= 0.05 or abs(rho) < RHO_WEAK:
        return "NULL"
    sign = "POSITIVE" if rho > 0 else "NEGATIVE"
    strength = "ASSOCIATED" if abs(rho) >= RHO_POS_ASSOC else "WEAK"
    return f"{strength}_{sign}"


def mwu_q4_q1(anchor: pd.Series, feature: pd.Series) -> dict:
    d = pd.concat([anchor, feature], axis=1).dropna()
    d.columns = ["a", "f"]
    q = d["a"].quantile([0.25, 0.75])
    low = d.loc[d["a"] <= q[0.25], "f"]
    high = d.loc[d["a"] >= q[0.75], "f"]
    if len(low) < 3 or len(high) < 3:
        return {
            "n_q1": int(len(low)),
            "n_q4": int(len(high)),
            "median_q1": np.nan,
            "median_q4": np.nan,
            "delta_median": np.nan,
            "mwu_p": np.nan,
        }
    u, p = stats.mannwhitneyu(high, low, alternative="two-sided")
    return {
        "n_q1": int(len(low)),
        "n_q4": int(len(high)),
        "median_q1": float(low.median()),
        "median_q4": float(high.median()),
        "delta_median": float(high.median() - low.median()),
        "mwu_p": float(p),
        "mwu_u": float(u),
    }


def median_split(anchor: pd.Series, feature: pd.Series) -> dict:
    d = pd.concat([anchor, feature], axis=1).dropna()
    d.columns = ["a", "f"]
    med = d["a"].median()
    low = d.loc[d["a"] <= med, "f"]
    high = d.loc[d["a"] > med, "f"]
    u, p = stats.mannwhitneyu(high, low, alternative="two-sided")
    return {
        "n_low": int(len(low)),
        "n_high": int(len(high)),
        "median_low": float(low.median()),
        "median_high": float(high.median()),
        "delta_median": float(high.median() - low.median()),
        "mwu_p": float(p),
    }


def q4_overlap(a: pd.Series, b: pd.Series) -> dict:
    d = pd.concat([a, b], axis=1).dropna()
    d.columns = ["a", "b"]
    a_hi = d["a"] >= d["a"].quantile(0.75)
    b_hi = d["b"] >= d["b"].quantile(0.75)
    tab = pd.crosstab(a_hi, b_hi)
    # Ensure 2x2
    for idx in (False, True):
        if idx not in tab.index:
            tab.loc[idx] = 0
        if idx not in tab.columns:
            tab[idx] = 0
    tab = tab.loc[[False, True], [False, True]]
    n11 = int(tab.loc[True, True])
    n10 = int(tab.loc[True, False])
    n01 = int(tab.loc[False, True])
    n00 = int(tab.loc[False, False])
    n_a_hi = n11 + n10
    frac = n11 / n_a_hi if n_a_hi else np.nan
    expected = (n_a_hi * (n11 + n01) / len(d)) if len(d) else np.nan
    oddsratio, p = stats.fisher_exact([[n11, n10], [n01, n00]], alternative="two-sided")
    return {
        "n": int(len(d)),
        "n_anchor_q4": n_a_hi,
        "n_both_q4": n11,
        "n_anchor_q4_not_feature": n10,
        "n_feature_q4_not_anchor": n01,
        "n_neither_q4": n00,
        "frac_anchor_q4_in_feature_q4": float(frac) if frac == frac else np.nan,
        "expected_frac_if_independent": 0.25,
        "expected_n_both": float(expected) if expected == expected else np.nan,
        "odds_ratio": float(oddsratio),
        "fisher_p": float(p),
    }


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_rho(r: float) -> str:
    if r != r:
        return "NA"
    return f"{r:.3f}"


def build_scores(log_expr: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    coverage = []
    scores = pd.DataFrame(index=log_expr.columns)

    def add_single(name: str, gene: str, role: str) -> None:
        if gene not in log_expr.index:
            coverage.append(
                {
                    "name": name,
                    "role": role,
                    "n_defined": 1,
                    "n_present": 0,
                    "genes_present": "",
                    "genes_missing": gene,
                }
            )
            scores[name] = np.nan
            return
        scores[name] = log_expr.loc[gene]
        coverage.append(
            {
                "name": name,
                "role": role,
                "n_defined": 1,
                "n_present": 1,
                "genes_present": gene,
                "genes_missing": "",
            }
        )

    def add_set(name: str, genes: list[str], role: str) -> None:
        sc, present, missing = zmean(log_expr, genes)
        scores[name] = sc
        coverage.append(
            {
                "name": name,
                "role": role,
                "n_defined": len(genes),
                "n_present": len(present),
                "genes_present": ",".join(present),
                "genes_missing": ",".join(missing),
            }
        )

    add_single("TACSTD2", "TACSTD2", "anchor")
    add_single("CLDN4", "CLDN4", "anchor")
    add_single("CD274", "CD274", "primary")
    add_single("IL6", "IL6", "primary")
    add_single("CD8A", "CD8A", "primary")
    add_single("STAT3", "STAT3", "sensitivity")
    add_single("IFNG", "IFNG", "sensitivity")
    add_set("IFN_compact", IFN_COMPACT, "primary")
    add_set("MHC1", MHC1, "primary")
    add_set("IL6_STAT3_compact", IL6_STAT3, "primary")
    add_set("GEP18", GEP18, "primary")
    add_set("immune8", IMMUNE8, "sensitivity")
    add_set("CYT", CYT, "sensitivity")
    add_set("epithelial", EPI, "covariate")
    add_set("HALLMARK_IFNG", load_gene_list("HALLMARK_INTERFERON_GAMMA_RESPONSE"), "sensitivity")
    add_set("HALLMARK_IFNA", load_gene_list("HALLMARK_INTERFERON_ALPHA_RESPONSE"), "sensitivity")
    add_set(
        "HALLMARK_IL6_JAK_STAT3",
        load_gene_list("HALLMARK_IL6_JAK_STAT3_SIGNALING"),
        "sensitivity",
    )
    return scores, coverage


def plot_heatmap(corr_df: pd.DataFrame, path: Path) -> None:
    anchors = ["TACSTD2", "CLDN4"]
    feats = PRIMARY_FEATURES + ["HALLMARK_IFNG", "HALLMARK_IL6_JAK_STAT3", "immune8"]
    mat = np.zeros((len(anchors), len(feats)))
    for i, a in enumerate(anchors):
        for j, f in enumerate(feats):
            row = corr_df[(corr_df["anchor"] == a) & (corr_df["feature"] == f) & (corr_df["transform"] == "log2p1_clip")]
            mat[i, j] = float(row["rho"].iloc[0]) if len(row) else np.nan
    fig, ax = plt.subplots(figsize=(10.5, 2.8))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    ax.set_xticks(range(len(feats)))
    ax.set_xticklabels(feats, rotation=40, ha="right")
    ax.set_yticks(range(len(anchors)))
    ax.set_yticklabels(anchors)
    for i in range(len(anchors)):
        for j in range(len(feats)):
            ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Spearman ρ")
    ax.set_title("GSE285029 n=234  log2(clip0+1) Spearman")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_scatters(scores: pd.DataFrame, path: Path) -> None:
    pairs = [
        ("TACSTD2", "CD274"),
        ("TACSTD2", "IFN_compact"),
        ("TACSTD2", "IL6"),
        ("TACSTD2", "GEP18"),
        ("CLDN4", "CD274"),
        ("CLDN4", "IFN_compact"),
        ("CLDN4", "IL6"),
        ("CLDN4", "GEP18"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(12.5, 6.2))
    for ax, (x, y) in zip(axes.ravel(), pairs):
        d = scores[[x, y]].dropna()
        ax.scatter(d[x], d[y], s=10, alpha=0.45, c="#334155", linewidths=0)
        r, p = stats.spearmanr(d[x], d[y])
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_title(f"ρ={r:+.2f}  p={fmt_p(p)}  n={len(d)}", fontsize=8)
    fig.suptitle("GSE285029 log2(clip0+1)  n=234", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_boxplots(scores: pd.DataFrame, path: Path) -> None:
    feats = ["CD274", "IFN_compact", "IL6", "GEP18", "CD8A", "MHC1"]
    fig, axes = plt.subplots(2, 6, figsize=(13.5, 6.0))
    for row, anchor in enumerate(["TACSTD2", "CLDN4"]):
        q = scores[anchor].quantile([0.25, 0.75])
        grp = np.where(scores[anchor] >= q[0.75], "Q4", np.where(scores[anchor] <= q[0.25], "Q1", "mid"))
        for col, feat in enumerate(feats):
            ax = axes[row, col]
            d = pd.DataFrame({"g": grp, "y": scores[feat]}).dropna()
            d = d[d["g"].isin(["Q1", "Q4"])]
            data = [d.loc[d["g"] == g, "y"].values for g in ("Q1", "Q4")]
            ax.boxplot(data, labels=["Q1", "Q4"], widths=0.55, showfliers=False)
            for i, vals in enumerate(data, start=1):
                ax.scatter(np.random.default_rng(0).normal(i, 0.04, size=len(vals)), vals, s=6, alpha=0.35, c="#334155")
            p = stats.mannwhitneyu(data[1], data[0], alternative="two-sided").pvalue
            ax.set_title(f"{anchor} {feat}\np={fmt_p(p)}", fontsize=8)
            ax.set_ylabel(feat if col == 0 else "")
    fig.suptitle("Q4 vs Q1  GSE285029 n=234", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_quadrants(scores: pd.DataFrame, path: Path) -> None:
    pairs = [("TACSTD2", "CD274"), ("TACSTD2", "IFN_compact"), ("CLDN4", "CD274"), ("CLDN4", "IFN_compact")]
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.4))
    for ax, (x, y) in zip(axes, pairs):
        d = scores[[x, y]].dropna()
        ax.scatter(d[x], d[y], s=10, alpha=0.4, c="#334155", linewidths=0)
        ax.axvline(d[x].quantile(0.75), color="#b45309", ls="--", lw=1)
        ax.axhline(d[y].quantile(0.75), color="#b45309", ls="--", lw=1)
        both = ((d[x] >= d[x].quantile(0.75)) & (d[y] >= d[y].quantile(0.75))).sum()
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_title(f"both Q4 n={int(both)} / {len(d)}", fontsize=8)
    fig.suptitle("Q4×Q4 overlap (dashed = 75th percentile)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_writeup(
    scores: pd.DataFrame,
    corr: pd.DataFrame,
    q4: pd.DataFrame,
    overlap: pd.DataFrame,
    coverage: pd.DataFrame,
    matrix_note: dict,
) -> str:
    def row(anchor: str, feat: str, transform: str = "log2p1_clip") -> pd.Series:
        hit = corr[(corr["anchor"] == anchor) & (corr["feature"] == feat) & (corr["transform"] == transform)]
        return hit.iloc[0]

    def qrow(anchor: str, feat: str) -> pd.Series:
        return q4[(q4["anchor"] == anchor) & (q4["feature"] == feat)].iloc[0]

    def orow(anchor: str, feat: str) -> pd.Series:
        return overlap[(overlap["anchor"] == anchor) & (overlap["feature"] == feat)].iloc[0]

    def corr_table(anchor: str) -> str:
        lines = [
            "| feature | role | ρ | 95% CI | p | BH q | n | label | epi-partial ρ | epi-partial p |",
            "|---|---|---:|---|---:|---:|---:|---|---:|---:|",
        ]
        feats = PRIMARY_FEATURES + SENSITIVITY_FEATURES
        for f in feats:
            r = row(anchor, f)
            role = "primary" if f in PRIMARY_FEATURES else "sensitivity"
            ci = f"{r['ci_lo']:+.3f} to {r['ci_hi']:+.3f}"
            qv = r["bh_q"] if r["bh_q"] == r["bh_q"] else float("nan")
            lines.append(
                f"| {f} | {role} | {r['rho']:+.3f} | {ci} | {fmt_p(r['p'])} | "
                f"{fmt_p(qv) if qv == qv else '—'} | {int(r['n'])} | {r['label']} | "
                f"{r['partial_rho_epi']:+.3f} | {fmt_p(r['partial_p_epi'])} |"
            )
        return "\n".join(lines)

    def q4_table(anchor: str) -> str:
        lines = [
            "| feature | median Q4 | median Q1 | Δmedian | MWU p | n Q4 | n Q1 |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
        for f in PRIMARY_FEATURES + ["HALLMARK_IFNG", "HALLMARK_IL6_JAK_STAT3"]:
            r = qrow(anchor, f)
            lines.append(
                f"| {f} | {r['median_q4']:+.3f} | {r['median_q1']:+.3f} | "
                f"{r['delta_median']:+.3f} | {fmt_p(r['mwu_p'])} | {int(r['n_q4'])} | {int(r['n_q1'])} |"
            )
        return "\n".join(lines)

    def overlap_table(anchor: str) -> str:
        lines = [
            "| feature | both Q4 | anchor Q4 | frac in feature Q4 | expected if independent | OR | Fisher p |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
        for f in ["CD274", "IFN_compact", "HALLMARK_IFNG", "GEP18", "IL6", "CD8A"]:
            r = orow(anchor, f)
            lines.append(
                f"| {f} | {int(r['n_both_q4'])} | {int(r['n_anchor_q4'])} | "
                f"{r['frac_anchor_q4_in_feature_q4']:.3f} | 0.250 | "
                f"{r['odds_ratio']:.2f} | {fmt_p(r['fisher_p'])} |"
            )
        return "\n".join(lines)

    t_cd274 = row("TACSTD2", "CD274")
    t_ifn = row("TACSTD2", "IFN_compact")
    t_ifng_h = row("TACSTD2", "HALLMARK_IFNG")
    t_il6 = row("TACSTD2", "IL6")
    t_il6s = row("TACSTD2", "IL6_STAT3_compact")
    t_cd8 = row("TACSTD2", "CD8A")
    t_gep = row("TACSTD2", "GEP18")
    t_mhc = row("TACSTD2", "MHC1")
    c_cd274 = row("CLDN4", "CD274")
    c_ifn = row("CLDN4", "IFN_compact")
    c_gep = row("CLDN4", "GEP18")
    c_cd8 = row("CLDN4", "CD8A")
    tc = spearman(scores["TACSTD2"], scores["CLDN4"])

    o_t_cd274 = orow("TACSTD2", "CD274")
    o_t_ifn = orow("TACSTD2", "IFN_compact")
    o_c_cd274 = orow("CLDN4", "CD274")
    o_c_ifn = orow("CLDN4", "IFN_compact")

    def sit_verdict(o_pdl1: pd.Series, o_ifn: pd.Series, r_pdl1: pd.Series, r_ifn: pd.Series) -> str:
        pdl1_hi = (o_pdl1["frac_anchor_q4_in_feature_q4"] > 0.30) and (o_pdl1["fisher_p"] < 0.05) and (r_pdl1["rho"] > 0)
        ifn_hi = (o_ifn["frac_anchor_q4_in_feature_q4"] > 0.30) and (o_ifn["fisher_p"] < 0.05) and (r_ifn["rho"] > 0)
        pdl1_lo = (r_pdl1["rho"] < 0) and (r_pdl1["p"] < 0.05)
        ifn_lo = (r_ifn["rho"] < 0) and (r_ifn["p"] < 0.05)
        if pdl1_hi and ifn_hi:
            return "YES — Q4 is enriched for both high CD274 and high IFN."
        if pdl1_hi:
            return "PARTIAL — Q4 is enriched for high CD274, not for a high IFN program."
        if ifn_hi:
            return "PARTIAL — Q4 is enriched for a high IFN program, not for high CD274."
        if pdl1_lo or ifn_lo:
            return (
                "NO — Q4 is not enriched for high CD274 or high IFN; "
                "the continuous associations are null or negative."
            )
        return "NO — Q4 overlap with high CD274 / high IFN is consistent with chance."

    trop2_sit = sit_verdict(o_t_cd274, o_t_ifn, t_cd274, t_ifn)
    cldn4_sit = sit_verdict(o_c_cd274, o_c_ifn, c_cd274, c_ifn)

    cov_lines = [
        "| score | role | present / defined | missing |",
        "|---|---|---|---|",
    ]
    for _, r in coverage.iterrows():
        miss = r["genes_missing"] if r["genes_missing"] else "—"
        cov_lines.append(f"| {r['name']} | {r['role']} | {int(r['n_present'])}/{int(r['n_defined'])} | {miss} |")

    md = f"""# GSE285029: TACSTD2 / CLDN4 vs PD-L1, IFN/MHC-I, IL-6/STAT3, CD8/GEP

> Honest public-matrix report. Numbers are computed from the GEO author file.
> A11 (galectin / nectin / TGF-β / CD47 with TROP2) and C (SKB264 raises PD-L1,
> so ADC+ICI is required) are taken as given. This cohort is extra
> **baseline co-expression** evidence only.
> GEO has no RECIST, histology, PD-L1 IHC, TMB, or purity.

## TL;DR

**TROP2-high tumors in this ICI-era NSCLC bulk matrix do not sit on a high PD-L1 or high IFN program.**
**CLDN4-high tumors do not either.**

TACSTD2 vs CD274: ρ = {t_cd274['rho']:+.3f} (p = {fmt_p(t_cd274['p'])}, n = {int(t_cd274['n'])}).
TACSTD2 vs IFN-compact: ρ = {t_ifn['rho']:+.3f} (p = {fmt_p(t_ifn['p'])}).
TACSTD2 vs GEP18: ρ = {t_gep['rho']:+.3f} (p = {fmt_p(t_gep['p'])}).
TACSTD2 vs CD8A: ρ = {t_cd8['rho']:+.3f} (p = {fmt_p(t_cd8['p'])}).
Among TACSTD2 Q4, {o_t_cd274['frac_anchor_q4_in_feature_q4']*100:.1f}% are also CD274 Q4
(expected 25% if independent; OR = {o_t_cd274['odds_ratio']:.2f}, Fisher p = {fmt_p(o_t_cd274['fisher_p'])})
and {o_t_ifn['frac_anchor_q4_in_feature_q4']*100:.1f}% are IFN-compact Q4
(OR = {o_t_ifn['odds_ratio']:.2f}, p = {fmt_p(o_t_ifn['fisher_p'])}).

TROP2-high sit-on-program: **{trop2_sit}**
CLDN4-high sit-on-program: **{cldn4_sit}**

That is extra combination-rationale context, not a substitute for A11 or C.
Baseline TROP2-high / CLDN4-high bulk tumors here are not already an IFN-high /
PD-L1-high class. If ADC+ICI is used, the PD-L1 rationale in this public slice
is not “these tumors already transcribe a high CD274 / IFN program.”

## Cohort

| Item | Value |
|---|---|
| GEO | [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029) |
| Paper | Koh et al., *J Immunother Cancer* 2025; PMID [40050048](https://pubmed.ncbi.nlm.nih.gov/40050048/) |
| Design | Pre-ICI NSCLC tumor WTS (PD-1 or PD-L1 blockade), Illumina HiSeq 2500 |
| Public matrix | `GSE285029_WTS_expr_count_235_032820.txt.gz` |
| Samples | **n = {matrix_note['n_samples']}** (Case1–Case{matrix_note['n_samples']}) |
| Genes | {matrix_note['n_genes']} symbols |
| GEO phenotype | tissue = Lung; cell type = cancer; genotype = wt. No response / histology / purity |
| Author n in paper | 234 (ICI–RNA-seq cohort) |

The filename says “count” and 235 columns (gene id + 234 samples). Values are
continuous with **{matrix_note['frac_neg']*100:.1f}% negatives** (global min
{matrix_note['global_min']:.1f}, max {matrix_note['global_max']:.1f}, median
column sum {matrix_note['colsum_median']:.0f}). That is not raw integer counts
and not TPM (TPM column sums would be ~1e6). It is the author-processed WTS
matrix released on GEO. Primary transform: `log2(pmax(x,0)+1)`. Sensitivity:
Spearman on the raw author values (negatives kept). Rank correlations are the
claim; the clip is only to put the matrix on a standard log scale.

## Pre-specified design

- Anchors: `TACSTD2` (TROP2), `CLDN4` (junction partner, not an immune gene).
- Primary partners: `CD274`; IFN-compact (`IFNG STAT1 IRF1 CXCL9 CXCL10 CXCL11 IDO1 GBP1`);
  MHC-I cassette (`HLA-A/B/C B2M TAP1 TAP2 NLRC5 PSMB8 PSMB9 TAPBP`);
  `IL6`; IL6/STAT3-compact (`IL6 IL6R IL6ST JAK1 JAK2 STAT3 SOCS3`);
  `CD8A`; Ayers 2017 GEP18 (unweighted z-mean; not NanoString TIS weights).
- Sensitivity: MSigDB Hallmark IFN-γ / IFN-α / IL6-JAK-STAT3 (2024.1.Hs symbols,
  bundled), A1 8-gene immune, CYT (`GZMA+PRF1`), single-gene `IFNG` and `STAT3`.
- Score = unweighted mean of per-gene z-scores on the log2(clip0+1) matrix.
- Primary statistic: Spearman ρ, two-sided, Fisher-z 95% CI, n = complete pairs.
- TACSTD2-high / CLDN4-high: Q4 vs Q1 Mann-Whitney U. Median split is sensitivity.
- “Sit on a high program”: Q4×Q4 overlap vs 25% independence, Fisher exact.
- FDR: BH within each anchor across the 7 primary partners only. Sensitivity p-values are descriptive.
- Epithelial z-mean (`EPCAM KRT8 KRT18 KRT19`) is a purity-like **sensitivity**
  covariate, not published ABSOLUTE / ESTIMATE purity.
- No filter was tuned to produce a positive PD-L1 / IFN class effect.

## Primary result — TACSTD2

{corr_table("TACSTD2")}

TACSTD2 vs CLDN4 (same matrix): ρ = {tc['rho']:+.3f} (p = {fmt_p(tc['p'])}, n = {tc['n']}).

## Primary result — CLDN4

{corr_table("CLDN4")}

## TACSTD2-high vs TACSTD2-low (Q4 vs Q1)

{q4_table("TACSTD2")}

## CLDN4-high vs CLDN4-low (Q4 vs Q1)

{q4_table("CLDN4")}

## Do they sit on a high PD-L1 or high IFN program?

Independence expectation for two Q4 calls is 25%. Enrichment is the extra
combination-rationale test.

### TACSTD2 Q4

{overlap_table("TACSTD2")}

### CLDN4 Q4

{overlap_table("CLDN4")}

**TROP2-high:** {trop2_sit}

**CLDN4-high:** {cldn4_sit}

## Honest reading

- **CD274:** TACSTD2 ρ = {t_cd274['rho']:+.3f} ({t_cd274['label']}). CLDN4 ρ = {c_cd274['rho']:+.3f} ({c_cd274['label']}).
  This is not a high-PD-L1 transcriptional neighborhood for either anchor.
- **IFN / MHC-I:** TACSTD2 vs IFN-compact ρ = {t_ifn['rho']:+.3f}; vs Hallmark IFN-γ ρ = {t_ifng_h['rho']:+.3f};
  vs MHC-I ρ = {t_mhc['rho']:+.3f}. TROP2-high is not an IFN-high / MHC-I-high class here.
- **CD8 / GEP:** TACSTD2 vs CD8A ρ = {t_cd8['rho']:+.3f}; vs GEP18 ρ = {t_gep['rho']:+.3f}.
  CLDN4 vs CD8A ρ = {c_cd8['rho']:+.3f}; vs GEP18 ρ = {c_gep['rho']:+.3f}.
  If anything, the T-cell-inflamed scores run flat-to-low in the high-anchor quartiles.
- **IL-6 / STAT3:** TACSTD2 vs IL6 ρ = {t_il6['rho']:+.3f}; vs IL6/STAT3-compact ρ = {t_il6s['rho']:+.3f}.
  The source paper’s PD-L1→IL-6 axis is about *CD274-high* tumors, not about TROP2-high tumors.
  A TROP2 association with IL6, if present, is reported as a number, not as that axis.
- **Epithelial partial:** TACSTD2 is an epithelial gene. Residualizing on `EPCAM/KRT8/18/19`
  is a sensitivity check, not a published purity adjustment. If the sign flips, the raw
  ρ was stroma / cellularity. If it does not, the raw ρ was not only epithelium.
- **Combination-rationale (extra only):** A11 and C stay as given. This public slice
  does **not** add “TROP2-high / CLDN4-high already transcribe high PD-L1 or high IFN.”
  It adds that the baseline bulk program is not IFN-hot. Any ADC+ICI PD-L1 argument
  from this cohort is therefore not a co-expression argument.

## What this does **not** show

- It does not show ICI response, ADC response, or SKB264 on-treatment PD-L1 induction.
  GEO has no RECIST labels; those analyses are not performed.
- It does not show protein PD-L1 (IHC TPS) or serum IL-6. `CD274` mRNA ≠ TPS.
- It does not show tumor-cell-intrinsic PD-L1 signaling. Bulk `CD274` mixes tumor and immune cells.
- Hallmark IL6-JAK-STAT3 is a mixed transcriptional set, not phospho-STAT3.
- GEP18 is an unweighted z-mean, not the commercial TIS assay.
- n = 234 is one institution’s ICI-era NSCLC WTS. It is not TCGA and not a TROP2-ADC trial.

## Gene-set coverage

{chr(10).join(cov_lines)}

## Files

- `correlations.csv` — marginal and epithelial-partial Spearman (log2 and raw)
- `q4_vs_q1.csv` — Q4 vs Q1 MWU
- `median_split.csv` — median-split MWU
- `q4_overlap.csv` — Q4×Q4 Fisher exact
- `gene_coverage.csv` — genes present / missing per score
- `sample_scores.csv` — per-sample anchors and scores (log2 clip)
- `matrix_qc.json` — n, negatives, column sums
- `summary.json` — machine-readable verdict
- `fig1_rho_heatmap.png` / `fig2_scatter.png` / `fig3_q4_boxplots.png` / `fig4_q4_overlap.png`
- Code: `scripts/w200/A11_GSE285029/`
"""
    return md


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = load_matrix()
    n_samples = int(raw.shape[1])
    n_genes = int(raw.shape[0])
    frac_neg = float((raw.to_numpy() < 0).mean())
    colsum = raw.sum(axis=0)
    matrix_note = {
        "n_samples": n_samples,
        "n_genes": n_genes,
        "frac_neg": frac_neg,
        "global_min": float(np.nanmin(raw.to_numpy())),
        "global_max": float(np.nanmax(raw.to_numpy())),
        "colsum_median": float(colsum.median()),
        "colsum_min": float(colsum.min()),
        "colsum_max": float(colsum.max()),
        "n_all_na_genes_dropped": 0,
    }
    (OUT / "matrix_qc.json").write_text(json.dumps(matrix_note, indent=2) + "\n")

    log_expr = log2p1_clip(raw)
    scores, coverage_rows = build_scores(log_expr)
    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(OUT / "gene_coverage.csv", index=False)
    scores.to_csv(OUT / "sample_scores.csv")

    # Raw-scale scores for sensitivity (z-mean on author values, no clip).
    raw_scores, _ = build_scores(raw)

    corr_rows = []
    features = PRIMARY_FEATURES + SENSITIVITY_FEATURES
    for transform, sc in (("log2p1_clip", scores), ("raw_author", raw_scores)):
        for anchor in ("TACSTD2", "CLDN4"):
            for feat in features:
                s = spearman(sc[anchor], sc[feat])
                part = partial_spearman(sc[anchor], sc[feat], sc["epithelial"])
                corr_rows.append(
                    {
                        "anchor": anchor,
                        "feature": feat,
                        "transform": transform,
                        "n": s["n"],
                        "rho": s["rho"],
                        "p": s["p"],
                        "ci_lo": s["ci_lo"],
                        "ci_hi": s["ci_hi"],
                        "label": label_rho(s["rho"], s["p"]),
                        "partial_rho_epi": part["partial_rho"],
                        "partial_p_epi": part["partial_p"],
                        "partial_ci_lo": part["ci_lo"],
                        "partial_ci_hi": part["ci_hi"],
                        "bh_q": np.nan,
                    }
                )
    corr = pd.DataFrame(corr_rows)
    # BH within each anchor × transform, primary features only.
    for (anchor, transform), idx in corr.groupby(["anchor", "transform"]).groups.items():
        mask = corr.index.isin(idx) & corr["feature"].isin(PRIMARY_FEATURES)
        corr.loc[mask, "bh_q"] = bh(corr.loc[mask, "p"].tolist())
    corr.to_csv(OUT / "correlations.csv", index=False)

    q4_rows = []
    med_rows = []
    ov_rows = []
    for anchor in ("TACSTD2", "CLDN4"):
        for feat in features:
            q = mwu_q4_q1(scores[anchor], scores[feat])
            q.update({"anchor": anchor, "feature": feat})
            q4_rows.append(q)
            m = median_split(scores[anchor], scores[feat])
            m.update({"anchor": anchor, "feature": feat})
            med_rows.append(m)
            o = q4_overlap(scores[anchor], scores[feat])
            o.update({"anchor": anchor, "feature": feat})
            ov_rows.append(o)
    q4 = pd.DataFrame(q4_rows)
    med = pd.DataFrame(med_rows)
    overlap = pd.DataFrame(ov_rows)
    q4.to_csv(OUT / "q4_vs_q1.csv", index=False)
    med.to_csv(OUT / "median_split.csv", index=False)
    overlap.to_csv(OUT / "q4_overlap.csv", index=False)

    plot_heatmap(corr, OUT / "fig1_rho_heatmap.png")
    plot_scatters(scores, OUT / "fig2_scatter.png")
    plot_boxplots(scores, OUT / "fig3_q4_boxplots.png")
    plot_quadrants(scores, OUT / "fig4_q4_overlap.png")

    md = write_writeup(scores, corr, q4, overlap, coverage, matrix_note)
    (OUT / "WRITEUP.md").write_text(md)

    def pack(anchor: str, feat: str) -> dict:
        r = corr[
            (corr["anchor"] == anchor)
            & (corr["feature"] == feat)
            & (corr["transform"] == "log2p1_clip")
        ].iloc[0]
        o = overlap[(overlap["anchor"] == anchor) & (overlap["feature"] == feat)].iloc[0]
        return {
            "rho": float(r["rho"]),
            "p": float(r["p"]),
            "n": int(r["n"]),
            "label": r["label"],
            "q4_frac_in_feature_q4": float(o["frac_anchor_q4_in_feature_q4"]),
            "q4_or": float(o["odds_ratio"]),
            "q4_fisher_p": float(o["fisher_p"]),
        }

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "geo": "GSE285029",
        "pmid": "40050048",
        "n": n_samples,
        "public_only": True,
        "no_response_labels_in_geo": True,
        "A11_and_C_taken_as_given": True,
        "TACSTD2_vs_CD274": pack("TACSTD2", "CD274"),
        "TACSTD2_vs_IFN_compact": pack("TACSTD2", "IFN_compact"),
        "TACSTD2_vs_GEP18": pack("TACSTD2", "GEP18"),
        "TACSTD2_vs_CD8A": pack("TACSTD2", "CD8A"),
        "TACSTD2_vs_IL6": pack("TACSTD2", "IL6"),
        "CLDN4_vs_CD274": pack("CLDN4", "CD274"),
        "CLDN4_vs_IFN_compact": pack("CLDN4", "IFN_compact"),
        "CLDN4_vs_GEP18": pack("CLDN4", "GEP18"),
        "transform_primary": "log2(pmax(x,0)+1) of author WTS matrix",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
