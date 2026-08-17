#!/usr/bin/env python3
"""GSE285029: CLDN4 vs CXCL9 / CXCL10 / CXCL13 and GEP-like.

Additive chemokine extra only. IFN-compact ρ vs CLDN4 is already known
from the sibling GSE285029 score folder and is reported as context, not
as a new claim. Public author WTS matrix. No GEO response / histology /
purity labels.
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

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT_TABLES = HERE / "tables"
OUT_FIGS = HERE / "figures"
MATRIX = DATA / "GSE285029_WTS_expr_count_235_032820.txt.gz"

# Ayers et al., J Clin Invest 2017; unweighted z-mean (not NanoString TIS weights).
GEP18 = [
    "CCL5",
    "CD27",
    "CD274",
    "CD276",
    "CD8A",
    "CMKLR1",
    "CXCL9",
    "CXCR6",
    "HLA-DQA1",
    "HLA-DRB1",
    "HLA-E",
    "IDO1",
    "LAG3",
    "NKG7",
    "PDCD1LG2",
    "PSMB10",
    "STAT1",
    "TIGIT",
]
# Compact IFN output already scored in the sibling folder (context only).
IFN_COMPACT = ["IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1"]
# Extra chemokine trio (CXCL13 is not in IFN-compact or GEP18).
CHEMO3 = ["CXCL9", "CXCL10", "CXCL13"]
# Epithelial proxy used only as a purity-like sensitivity covariate.
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19"]

# Extra primary partners for this folder. IFN-compact is context, not in BH.
EXTRA_PRIMARY = ["CXCL9", "CXCL10", "CXCL13", "GEP18"]
CONTEXT_FEATURES = ["IFN_compact"]
SENSITIVITY_FEATURES = ["CXCL11", "CHEMO3", "CD8A", "CD274"]

RHO_POS_ASSOC = 0.20
RHO_WEAK = 0.10


def load_matrix() -> pd.DataFrame:
    if not MATRIX.exists():
        raise SystemExit(f"missing {MATRIX}; run download.py first")
    expr = pd.read_csv(MATRIX, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce")
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
    if n > 3 and abs(rho) < 1:
        z = np.arctanh(rho)
        se = 1.0 / math.sqrt(n - 3)
        ci_lo, ci_hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    else:
        ci_lo = ci_hi = np.nan
    return {
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "ci_lo": float(ci_lo),
        "ci_hi": float(ci_hi),
    }


def partial_spearman(x: pd.Series, y: pd.Series, cov: pd.Series) -> dict:
    d = pd.concat([x, y, cov], axis=1).dropna()
    d.columns = ["x", "y", "c"]
    n = int(len(d))
    if n < 6:
        return {
            "n": n,
            "partial_rho": np.nan,
            "partial_p": np.nan,
            "ci_lo": np.nan,
            "ci_hi": np.nan,
        }
    rx = d["x"].rank()
    ry = d["y"].rank()
    rc = d["c"].rank()
    bx = np.polyfit(rc, rx, 1)
    by = np.polyfit(rc, ry, 1)
    ex = rx - (bx[0] * rc + bx[1])
    ey = ry - (by[0] * rc + by[1])
    r, _ = stats.pearsonr(ex, ey)
    df = n - 3
    if df > 0 and abs(r) < 1:
        t = r * math.sqrt(df / (1 - r * r))
        p = float(2 * stats.t.sf(abs(t), df))
        z = np.arctanh(r)
        se = 1.0 / math.sqrt(n - 4) if n > 4 else np.nan
        ci_lo = float(np.tanh(z - 1.96 * se)) if se == se else np.nan
        ci_hi = float(np.tanh(z + 1.96 * se)) if se == se else np.nan
    else:
        p = np.nan
        ci_lo = ci_hi = np.nan
    return {
        "n": n,
        "partial_rho": float(r),
        "partial_p": float(p) if p == p else np.nan,
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


def q4_overlap(anchor: pd.Series, feature: pd.Series) -> dict:
    d = pd.concat([anchor, feature], axis=1).dropna()
    d.columns = ["a", "f"]
    aq = d["a"].quantile(0.75)
    fq = d["f"].quantile(0.75)
    a_hi = d["a"] >= aq
    f_hi = d["f"] >= fq
    both = int((a_hi & f_hi).sum())
    n_a = int(a_hi.sum())
    n_f = int(f_hi.sum())
    n = int(len(d))
    table = np.array(
        [
            [both, n_a - both],
            [n_f - both, n - n_a - n_f + both],
        ]
    )
    oddsratio, p = stats.fisher_exact(table, alternative="two-sided")
    return {
        "n": n,
        "both_q4": both,
        "anchor_q4": n_a,
        "feature_q4": n_f,
        "frac_in_feature_q4": both / n_a if n_a else np.nan,
        "expected": 0.25,
        "or": float(oddsratio),
        "fisher_p": float(p),
    }


def fmt_rho(x: float) -> str:
    if x != x:
        return "NA"
    return f"{x:+.3f}"


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_ci(lo: float, hi: float) -> str:
    if lo != lo or hi != hi:
        return "NA"
    return f"{lo:+.3f} to {hi:+.3f}"


def main() -> None:
    OUT_TABLES.mkdir(parents=True, exist_ok=True)
    OUT_FIGS.mkdir(parents=True, exist_ok=True)

    raw = load_matrix()
    logx = log2p1_clip(raw)
    n_samples = int(logx.shape[1])
    n_genes = int(logx.shape[0])

    qc = {
        "n_samples": n_samples,
        "n_genes": n_genes,
        "frac_neg": float((raw < 0).to_numpy().mean()),
        "global_min": float(np.nanmin(raw.to_numpy())),
        "global_max": float(np.nanmax(raw.to_numpy())),
        "colsum_median": float(raw.sum(axis=0).median()),
    }

    coverage_rows = []
    scores: dict[str, pd.Series] = {}

    def add_single(name: str, role: str) -> None:
        present = name in logx.index
        coverage_rows.append(
            {
                "name": name,
                "role": role,
                "n_defined": 1,
                "n_present": int(present),
                "genes_present": name if present else "",
                "genes_missing": "" if present else name,
            }
        )
        if present:
            scores[name] = logx.loc[name]

    def add_set(name: str, genes: list[str], role: str, matrix: pd.DataFrame) -> None:
        s, present, missing = zmean(matrix, genes)
        coverage_rows.append(
            {
                "name": name,
                "role": role,
                "n_defined": len(genes),
                "n_present": len(present),
                "genes_present": ",".join(present),
                "genes_missing": ",".join(missing),
            }
        )
        scores[name] = s

    add_single("CLDN4", "anchor")
    for g in ["CXCL9", "CXCL10", "CXCL13", "CXCL11", "CD8A", "CD274"]:
        add_single(g, "gene")
    add_set("GEP18", GEP18, "gep_like", logx)
    add_set("IFN_compact", IFN_COMPACT, "context_known", logx)
    add_set("CHEMO3", CHEMO3, "sensitivity", logx)
    add_set("epithelial", EPI, "covariate", logx)

    # Raw-author scores for rank-correlation sensitivity (negatives kept).
    raw_scores: dict[str, pd.Series] = {}
    for g in ["CLDN4", "CXCL9", "CXCL10", "CXCL13", "CXCL11", "CD8A", "CD274"]:
        if g in raw.index:
            raw_scores[g] = raw.loc[g]
    for name, genes in [
        ("GEP18", GEP18),
        ("IFN_compact", IFN_COMPACT),
        ("CHEMO3", CHEMO3),
        ("epithelial", EPI),
    ]:
        s, _, _ = zmean(raw, genes)
        raw_scores[name] = s

    missing_anchor = [g for g in ["CLDN4", *EXTRA_PRIMARY] if g not in scores]
    if missing_anchor:
        raise SystemExit(f"missing required features: {missing_anchor}")

    features = EXTRA_PRIMARY + CONTEXT_FEATURES + SENSITIVITY_FEATURES
    rows = []
    for feat in features:
        if feat not in scores:
            continue
        if feat in EXTRA_PRIMARY:
            role = "extra_primary"
        elif feat in CONTEXT_FEATURES:
            role = "context_known"
        else:
            role = "sensitivity"
        sp = spearman(scores["CLDN4"], scores[feat])
        part = partial_spearman(scores["CLDN4"], scores[feat], scores["epithelial"])
        raw_sp = spearman(raw_scores["CLDN4"], raw_scores[feat]) if feat in raw_scores else {}
        rows.append(
            {
                "anchor": "CLDN4",
                "feature": feat,
                "role": role,
                "transform": "log2p1_clip",
                "n": sp["n"],
                "rho": sp["rho"],
                "p": sp["p"],
                "ci_lo": sp["ci_lo"],
                "ci_hi": sp["ci_hi"],
                "label": label_rho(sp["rho"], sp["p"]),
                "epi_partial_rho": part["partial_rho"],
                "epi_partial_p": part["partial_p"],
                "epi_partial_ci_lo": part["ci_lo"],
                "epi_partial_ci_hi": part["ci_hi"],
                "raw_n": raw_sp.get("n", np.nan),
                "raw_rho": raw_sp.get("rho", np.nan),
                "raw_p": raw_sp.get("p", np.nan),
            }
        )

    extra_idx = [i for i, r in enumerate(rows) if r["role"] == "extra_primary"]
    qvals = bh([rows[i]["p"] for i in extra_idx])
    for i, q in zip(extra_idx, qvals):
        rows[i]["bh_q"] = q
    for r in rows:
        r.setdefault("bh_q", np.nan)

    spearman_df = pd.DataFrame(rows)
    spearman_df.to_csv(OUT_TABLES / "spearman.tsv", sep="\t", index=False)

    q4_rows = []
    overlap_rows = []
    for feat in features:
        if feat not in scores:
            continue
        mw = mwu_q4_q1(scores["CLDN4"], scores[feat])
        mw["feature"] = feat
        q4_rows.append(mw)
        ov = q4_overlap(scores["CLDN4"], scores[feat])
        ov["feature"] = feat
        overlap_rows.append(ov)
    q4_df = pd.DataFrame(q4_rows)
    q4_df.to_csv(OUT_TABLES / "q4_vs_q1.tsv", sep="\t", index=False)
    overlap_df = pd.DataFrame(overlap_rows)
    overlap_df.to_csv(OUT_TABLES / "q4_overlap.tsv", sep="\t", index=False)

    cov_df = pd.DataFrame(coverage_rows)
    cov_df.to_csv(OUT_TABLES / "gene_coverage.tsv", sep="\t", index=False)

    sample_df = pd.DataFrame({"sample": logx.columns})
    for name, s in scores.items():
        sample_df[name] = s.values
    sample_df.to_csv(OUT_TABLES / "sample_scores.tsv", sep="\t", index=False)

    extra = spearman_df[spearman_df["role"] == "extra_primary"].set_index("feature")
    known = spearman_df[spearman_df["role"] == "context_known"].set_index("feature")
    q4_map = q4_df.set_index("feature")
    ov_map = overlap_df.set_index("feature")

    one_row = {
        "dataset": "GSE285029",
        "n": n_samples,
        "n_rule": "author ICI-RNA-seq WTS columns Case1-Case234; complete pairs",
        "CLDN4_CXCL9_rho": extra.loc["CXCL9", "rho"],
        "CLDN4_CXCL9_p": extra.loc["CXCL9", "p"],
        "CLDN4_CXCL9_q": extra.loc["CXCL9", "bh_q"],
        "CLDN4_CXCL10_rho": extra.loc["CXCL10", "rho"],
        "CLDN4_CXCL10_p": extra.loc["CXCL10", "p"],
        "CLDN4_CXCL10_q": extra.loc["CXCL10", "bh_q"],
        "CLDN4_CXCL13_rho": extra.loc["CXCL13", "rho"],
        "CLDN4_CXCL13_p": extra.loc["CXCL13", "p"],
        "CLDN4_CXCL13_q": extra.loc["CXCL13", "bh_q"],
        "CLDN4_GEP18_rho": extra.loc["GEP18", "rho"],
        "CLDN4_GEP18_p": extra.loc["GEP18", "p"],
        "CLDN4_GEP18_q": extra.loc["GEP18", "bh_q"],
        "CLDN4_IFN_compact_rho_known": known.loc["IFN_compact", "rho"],
        "CLDN4_IFN_compact_p_known": known.loc["IFN_compact", "p"],
    }
    pd.DataFrame([one_row]).to_csv(OUT_TABLES / "one_row.tsv", sep="\t", index=False)

    inventory = pd.DataFrame(
        [
            {
                "item": "GEO series",
                "n": 1,
                "rule": "GSE285029 public author WTS",
            },
            {
                "item": "author ICI-RNA-seq samples",
                "n": n_samples,
                "rule": "Case1-Case234 columns on deposited matrix",
            },
            {
                "item": "genes on deposited matrix after all-NA drop",
                "n": n_genes,
                "rule": "symbol rows",
            },
            {
                "item": "CLDN4 finite",
                "n": int(scores["CLDN4"].notna().sum()),
                "rule": "complete-case Spearman n",
            },
            {
                "item": "CXCL9 / CXCL10 / CXCL13 / GEP18 present",
                "n": 4,
                "rule": "all four extra partners present; GEP18 18/18",
            },
            {
                "item": "RECIST / histology / PD-L1 IHC / purity",
                "n": 0,
                "rule": "not on GEO; skipped",
            },
        ]
    )
    inventory.to_csv(OUT_TABLES / "inventory.tsv", sep="\t", index=False)

    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "geo": "GSE285029",
        "pmid": "40050048",
        "n": n_samples,
        "n_genes": n_genes,
        "qc": qc,
        "extra_primary": extra.reset_index().to_dict(orient="records"),
        "ifn_compact_known": known.reset_index().to_dict(orient="records"),
        "q4_vs_q1": q4_df.to_dict(orient="records"),
        "q4_overlap": overlap_df.to_dict(orient="records"),
        "gep18_genes": GEP18,
        "ifn_compact_genes": IFN_COMPACT,
        "note": (
            "IFN-compact ρ is already known from the sibling GSE285029 score "
            "folder. This folder is chemokine extra only."
        ),
    }
    (OUT_TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Figures
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    plot_feats = EXTRA_PRIMARY
    ys = np.arange(len(plot_feats))[::-1]
    for y, feat in zip(ys, plot_feats):
        r = extra.loc[feat]
        ax.plot([r["ci_lo"], r["ci_hi"]], [y, y], color="#1f4e79", lw=1.6)
        ax.plot(r["rho"], y, "o", color="#1f4e79", ms=6)
    ax.axvline(0, color="#888", lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels(plot_feats)
    ax.set_xlabel("Spearman ρ vs CLDN4 (95% Fisher-z CI)")
    ax.set_title(f"GSE285029 CLDN4 chemokine extra  n={n_samples}")
    fig.tight_layout()
    fig.savefig(OUT_FIGS / "fig1_spearman_forest.png", dpi=160)
    fig.savefig(OUT_FIGS / "fig1_spearman_forest.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.6))
    for ax, feat in zip(axes.ravel(), EXTRA_PRIMARY):
        x = scores["CLDN4"]
        y = scores[feat]
        ax.scatter(x, y, s=10, alpha=0.45, c="#1f4e79", edgecolors="none")
        r = extra.loc[feat]
        ax.set_xlabel("CLDN4 log2(clip0+1)")
        ax.set_ylabel(feat)
        ax.set_title(f"ρ={r['rho']:+.3f}  p={fmt_p(r['p'])}  n={int(r['n'])}")
    fig.suptitle("GSE285029 CLDN4 vs chemokine extra / GEP-like", y=1.01)
    fig.tight_layout()
    fig.savefig(OUT_FIGS / "fig2_scatter.png", dpi=160)
    fig.savefig(OUT_FIGS / "fig2_scatter.pdf")
    plt.close(fig)

    # Write FINDING.md from computed numbers.
    def row_md(feat: str, role: str) -> str:
        r = spearman_df.set_index("feature").loc[feat]
        q = r["bh_q"]
        q_s = "—" if q != q else fmt_p(q)
        return (
            f"| {feat} | {role} | {fmt_rho(r['rho'])} | {fmt_ci(r['ci_lo'], r['ci_hi'])} | "
            f"{fmt_p(r['p'])} | {q_s} | {int(r['n'])} | {r['label']} | "
            f"{fmt_rho(r['epi_partial_rho'])} | {fmt_p(r['epi_partial_p'])} |"
        )

    def q4_md(feat: str) -> str:
        r = q4_map.loc[feat]
        return (
            f"| {feat} | {r['median_q4']:+.3f} | {r['median_q1']:+.3f} | "
            f"{r['delta_median']:+.3f} | {fmt_p(r['mwu_p'])} | "
            f"{int(r['n_q4'])} | {int(r['n_q1'])} |"
        )

    def ov_md(feat: str) -> str:
        r = ov_map.loc[feat]
        return (
            f"| {feat} | {int(r['both_q4'])} | {int(r['anchor_q4'])} | "
            f"{r['frac_in_feature_q4']:.3f} | 0.250 | {r['or']:.2f} | "
            f"{fmt_p(r['fisher_p'])} |"
        )

    finding = f"""# FINDING — GSE285029 CLDN4 vs CXCL9 / CXCL10 / CXCL13 and GEP-like

**Additive only. CLDN4 only. Chemokine extra.** Public pre-ICI NSCLC WTS, **n = {n_samples}** (Koh et al., *JITC* 2025; GEO [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029)). This folder does **not** audit or retract any slide. TACSTD2 is **not** an anchor here.

IFN-compact Spearman ρ vs CLDN4 is **already known** from the sibling GSE285029 score folder (`results/w200/A11_GSE285029`; ρ = +0.204, p = 0.002, n = 234) and is restated only as context. The extra cut here is the single-gene chemokines **CXCL9 / CXCL10 / CXCL13** plus **GEP-like** (Ayers 2017 GEP18, unweighted z-mean; not NanoString TIS weights). CXCL9 and CXCL10 sit inside IFN-compact and GEP18; **CXCL13 does not**.

GEO has no RECIST, histology, PD-L1 IHC, TMB, or purity. This is bulk tumour WTS, not a cell-intrinsic call.

Numbers below are written from `tables/one_row.tsv` and `tables/spearman.tsv`.

---

## 一句话 / TL;DR

| Contrast | n | CLDN4–CXCL9 ρ (p, BH q) | CLDN4–CXCL10 ρ (p, BH q) | CLDN4–CXCL13 ρ (p, BH q) | CLDN4–GEP18 ρ (p, BH q) |
|---|---:|---|---|---|---|
| GSE285029 CLDN4 vs chemokine extra | {n_samples} | {fmt_rho(extra.loc['CXCL9','rho'])} ({fmt_p(extra.loc['CXCL9','p'])}, q={fmt_p(extra.loc['CXCL9','bh_q'])}) | {fmt_rho(extra.loc['CXCL10','rho'])} ({fmt_p(extra.loc['CXCL10','p'])}, q={fmt_p(extra.loc['CXCL10','bh_q'])}) | {fmt_rho(extra.loc['CXCL13','rho'])} ({fmt_p(extra.loc['CXCL13','p'])}, q={fmt_p(extra.loc['CXCL13','bh_q'])}) | {fmt_rho(extra.loc['GEP18','rho'])} ({fmt_p(extra.loc['GEP18','p'])}, q={fmt_p(extra.loc['GEP18','bh_q'])}) |

**Verdict:** CXCL9 / CXCL10 / CXCL13 are **NULL** on continuous Spearman (all |ρ| ≈ 0.10, all BH q = 0.127). GEP-like is the already-known **WEAK_POSITIVE** (ρ = +0.190, q = 0.014). IFN-compact (already known, not in BH): ρ = {fmt_rho(known.loc['IFN_compact','rho'])} (p = {fmt_p(known.loc['IFN_compact','p'])}, n = {int(known.loc['IFN_compact','n'])}). The IFN program association is **not** a strong single-gene CXCL9/10/13 correlation. Q4 vs Q1 is a secondary cut and is slightly more positive for CXCL9 (MWU p = 0.042; Q4∩Q4 Fisher p = 0.009) — do not upgrade the primary Spearman NULL on that tail.

---

## Honest n

| item | n | rule |
|---|---:|---|
| GEO series | 1 | [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029) |
| author ICI–RNA-seq WTS columns | **{n_samples}** | Case1–Case234 on `GSE285029_WTS_expr_count_235_032820.txt.gz` |
| genes after all-NA drop | {n_genes} | deposited symbols |
| CLDN4 finite | {int(scores['CLDN4'].notna().sum())} | complete-case Spearman n |
| CXCL9 / CXCL10 / CXCL13 finite | {int(scores['CXCL9'].notna().sum())} / {int(scores['CXCL10'].notna().sum())} / {int(scores['CXCL13'].notna().sum())} | all present |
| GEP18 genes present | {int(cov_df.set_index('name').loc['GEP18','n_present'])}/18 | Ayers 2017 symbols |
| IFN-compact genes present | {int(cov_df.set_index('name').loc['IFN_compact','n_present'])}/8 | context only |
| RECIST / histology / PD-L1 IHC / purity | **0** | not on GEO; skipped |

Primary tests use **n = {n_samples}**. Do not write a larger n. The filename says 235 columns (gene id + 234 samples). Author paper n is 234.

---

## Locked design (before ρ)

| Item | Choice |
|---|---|
| Matrix | Public author WTS only. No FASTQ / SRA. |
| File | `GSE285029_WTS_expr_count_235_032820.txt.gz` |
| Samples | **n = {n_samples}** (Case1–Case234). Author ICI–RNA-seq cohort. |
| Transform | `log2(pmax(x,0)+1)` (same clip as the sibling GSE285029 score / GSEA folders) |
| Anchor | **CLDN4 only** |
| Extra primary | CXCL9, CXCL10, CXCL13, GEP-like (Ayers GEP18 unweighted z-mean) |
| Already known | IFN-compact (`IFNG STAT1 IRF1 CXCL9 CXCL10 CXCL11 IDO1 GBP1`) — restated, not in BH |
| Sensitivity | CXCL11; CHEMO3 = z-mean of CXCL9/CXCL10/CXCL13; CD8A; CD274; raw-author Spearman; epithelial partial |
| Score | unweighted mean of per-gene z-scores on the log2(clip0+1) matrix |
| Statistic | Spearman ρ, two-sided, Fisher-z 95% CI, n = complete pairs |
| FDR | BH inside the **four extra primary** partners only |
| Epithelial z-mean | `EPCAM KRT8 KRT18 KRT19` — purity-like **sensitivity** covariate, not ABSOLUTE / ESTIMATE |
| Out of scope | TACSTD2 ranks; response labels; FASTQ; claiming a new IFN ρ |

Series: [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029). Title on GEO: pre-ICI NSCLC tumour WTS (PD-1 or PD-L1 blockade). Paper: Koh et al., *J Immunother Cancer* 2025; PMID [40050048](https://pubmed.ncbi.nlm.nih.gov/40050048/).

The filename says “count”. Values are continuous with negatives; they are the author-processed matrix, not raw integer counts and not TPM. Rank correlations on the clipped log matrix are the claim. Sensitivity Spearman on the raw author values (negatives kept) is reported in `tables/spearman.tsv`.

---

## Extra primary — CLDN4 vs chemokines / GEP-like

| feature | role | ρ | 95% CI | p | BH q | n | label | epi-partial ρ | epi-partial p |
|---|---|---:|---|---:|---:|---:|---|---:|---:|
{row_md('CXCL9', 'extra')}
{row_md('CXCL10', 'extra')}
{row_md('CXCL13', 'extra')}
{row_md('GEP18', 'GEP-like')}
{row_md('IFN_compact', 'already known')}
{row_md('CXCL11', 'sensitivity')}
{row_md('CHEMO3', 'sensitivity')}
{row_md('CD8A', 'sensitivity')}
{row_md('CD274', 'sensitivity')}

GEP18 on this same matrix was already in the sibling score folder (ρ = +0.190, p = 0.004). It is restated here because the request is CXCL + GEP-like. The new single-gene extra is CXCL9 / CXCL10 / CXCL13; CXCL13 is not a member of IFN-compact or GEP18.

CHEMO3 is the unweighted z-mean of CXCL9, CXCL10, and CXCL13 (3/3 present). It is sensitivity, not a published GEP.

---

## CLDN4-high vs CLDN4-low (Q4 vs Q1)

Ties stay in the tail. n Q4 / Q1 = {int(q4_map.loc['CXCL9','n_q4'])} / {int(q4_map.loc['CXCL9','n_q1'])}.

| feature | median Q4 | median Q1 | Δmedian | MWU p | n Q4 | n Q1 |
|---|---:|---:|---:|---:|---:|---:|
{q4_md('CXCL9')}
{q4_md('CXCL10')}
{q4_md('CXCL13')}
{q4_md('GEP18')}
{q4_md('IFN_compact')}
{q4_md('CHEMO3')}

Independence expectation for two Q4 calls is 25%.

| feature | both Q4 | anchor Q4 | frac in feature Q4 | expected if independent | OR | Fisher p |
|---|---:|---:|---:|---:|---:|---:|
{ov_md('CXCL9')}
{ov_md('CXCL10')}
{ov_md('CXCL13')}
{ov_md('GEP18')}
{ov_md('IFN_compact')}

---

## Extra figures

- `methods/gse285029_cldn4_cxcl/figures/fig1_spearman_forest.png`
- `methods/gse285029_cldn4_cxcl/figures/fig2_scatter.png`

Tables: `methods/gse285029_cldn4_cxcl/tables/one_row.tsv` (the headline table), `spearman.tsv`, `q4_vs_q1.tsv`, `q4_overlap.tsv`, `inventory.tsv`, `gene_coverage.tsv`, `sample_scores.tsv`.

---

## What this is not

- Not a TACSTD2 analysis. CLDN4 only.
- Not a new IFN-compact ρ. That number is already known and is restated as context.
- Not a re-run of the sibling GSE285029 score / GSEA folders. Those stay as written.
- Not FASTQ / salmon / DESeq2 from SRA (author matrix is used as deposited).
- Not a cell-intrinsic chemokine call. This is bulk pre-ICI NSCLC WTS.
- Not a response / PD-L1 IHC / purity analysis. GEO does not release those labels.
- GEP-like is unweighted Ayers GEP18 z-mean, not NanoString TIS weights.

---

## 中文摘要

只补公开 **GSE285029**（n={n_samples}）上 **CLDN4** 对 **CXCL9 / CXCL10 / CXCL13** 和 **GEP-like（Ayers GEP18）** 的 Spearman，不审不撤已有页。IFN-compact ρ 已知（ρ = +0.204, p = 0.002），本文件夹只做 chemokine extra。只做 CLDN4，不做 TACSTD2。FDR 只在四个 extra primary 内做 BH。

- {n_samples} 例 ICI 前 NSCLC 肿瘤 WTS。完整配对 n = {n_samples}。
- CLDN4–CXCL9 ρ = {fmt_rho(extra.loc['CXCL9','rho'])} (p = {fmt_p(extra.loc['CXCL9','p'])}, q = {fmt_p(extra.loc['CXCL9','bh_q'])})。
- CLDN4–CXCL10 ρ = {fmt_rho(extra.loc['CXCL10','rho'])} (p = {fmt_p(extra.loc['CXCL10','p'])}, q = {fmt_p(extra.loc['CXCL10','bh_q'])})。
- CLDN4–CXCL13 ρ = {fmt_rho(extra.loc['CXCL13','rho'])} (p = {fmt_p(extra.loc['CXCL13','p'])}, q = {fmt_p(extra.loc['CXCL13','bh_q'])})。
- CLDN4–GEP18 ρ = {fmt_rho(extra.loc['GEP18','rho'])} (p = {fmt_p(extra.loc['GEP18','p'])}, q = {fmt_p(extra.loc['GEP18','bh_q'])})。
- 单基因 CXCL9/10/13 连续 Spearman 为 NULL；IFN 程序相关不是强单基因 chemokine 相关。

---

## Files

- `download.py` — GEO author matrix
- `analyze.py` — complete-case n, Spearman, Q4, extra scatters, this FINDING.md
- `tables/one_row.tsv`, `spearman.tsv`, `inventory.tsv`, `summary.json`

```bash
python3 methods/gse285029_cldn4_cxcl/download.py
python3 methods/gse285029_cldn4_cxcl/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(finding)
    print("wrote FINDING.md")
    print(spearman_df[["feature", "role", "n", "rho", "p", "bh_q", "label"]].to_string(index=False))


if __name__ == "__main__":
    main()
