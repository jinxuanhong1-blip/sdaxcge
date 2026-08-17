#!/usr/bin/env python3
"""CLDN4-only tumor-cell-intrinsic patient-pseudobulk DE.

Triple: GSE131907 + GSE205335 + GSE189357 (no GSE148071).
The GSE123902+GSE189357 pair DE is already filled (PR #469);
this extra is the remaining 189357-axis triple that differs in PR #459.

Not infiltrate. Not CellChat. No dual-high.

Method: patient/sample UMI-sum → TMM → log2(CPM+1) → OLS.
Families: IFN, MHC-I/APM, TJ (CLDN4 held out).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"

COMBO = "GSE131907+GSE205335+GSE189357"
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}

# PR #459 given cuts (T/NK ρ not re-audited)
GIVEN_PAIR = {
    "combo": "GSE123902+GSE189357",
    "score": "pct",
    "n": 22,
    "rho": -0.638,
    "note": "PR #459 pair that differs; pair DE already filled in PR #469",
}
GIVEN_TRIPLE = {
    "combo": COMBO,
    "score": "pct",
    "n": 52,
    "rho": -0.497,
    "p": 0.00035,
    "note": "PR #459 triple %pos; T/NK ρ not re-audited; no GSE148071",
}

FOCAL_BOX = [
    "CLDN4", "STAT1", "IRF1", "HLA-A", "HLA-B", "B2M", "TAP1", "TAP2",
    "NLRC5", "OCLN", "TJP1", "CLDN1", "CLDN7", "ISG15", "MX1", "PSMB8",
]
JUNCTION_FOCAL = [
    "CLDN1", "CLDN7", "F11R", "PARD3", "OCLN", "TJP1", "CDH1",
]
FAM_COLORS = {
    "IFN": "#d62728",
    "MHC-I/APM": "#1f77b4",
    "TJ": "#2ca02c",
}


def load_families() -> dict[str, set[str]]:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    tj = (
        set(sets["KEGG_TIGHT_JUNCTION"])
        | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
        | set(sets.get("GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY", []))
        | set(JUNCTION_FOCAL)
        | set(g for g in a8.get("focal_genes", []) if g not in sets.get("KRT_EPITHELIAL", []))
    )
    tj.discard("CLDN4")  # held out: split gene
    return {
        "IFN": set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
        | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]),
        "MHC-I/APM": set(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]),
        "TJ": tj,
    }


def family_of(gene: str, fam: dict[str, set[str]]) -> str:
    hits = [k for k, vs in fam.items() if gene in vs]
    return "|".join(hits) if hits else "other"


def tmm_norm_factors(counts: pd.DataFrame) -> pd.Series:
    lib = counts.sum(axis=0).astype(float).replace(0, np.nan)
    rel = counts.div(lib, axis=1)
    f75 = rel.quantile(0.75, axis=0)
    ref = (f75 - f75.mean()).abs().idxmin()
    ref_c = counts[ref].astype(float)
    ref_lib = float(lib[ref])
    factors = {}
    for col in counts.columns:
        obs = counts[col].astype(float)
        obs_lib = float(lib[col])
        keep = (obs > 0) & (ref_c > 0)
        if keep.sum() < 50:
            factors[col] = 1.0
            continue
        m = np.log2((obs[keep] / obs_lib) / (ref_c[keep] / ref_lib))
        a = 0.5 * np.log2((obs[keep] / obs_lib) * (ref_c[keep] / ref_lib))
        w = (obs_lib - obs[keep]) / (obs_lib * obs[keep]) + (ref_lib - ref_c[keep]) / (
            ref_lib * ref_c[keep]
        )
        ok = np.isfinite(m) & np.isfinite(a) & np.isfinite(w) & (w > 0)
        m, a, w = m[ok], a[ok], w[ok]
        if len(m) < 50:
            factors[col] = 1.0
            continue
        lo_m, hi_m = np.quantile(m, [0.30, 0.70])
        lo_a, hi_a = np.quantile(a, [0.05, 0.95])
        trim = (m >= lo_m) & (m <= hi_m) & (a >= lo_a) & (a <= hi_a)
        if trim.sum() < 20:
            factors[col] = 1.0
            continue
        tmm = float(np.average(m[trim], weights=1.0 / w[trim]))
        factors[col] = 2 ** tmm
    fac = pd.Series(factors)
    return fac / fac.mean()


def log_cpm(counts: pd.DataFrame, factors: pd.Series) -> pd.DataFrame:
    lib = counts.sum(axis=0).astype(float) * factors.reindex(counts.columns).astype(float)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)


def filter_genes(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 3) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def assign_quartiles(values: pd.Series) -> pd.Series:
    s = values.astype(float)
    ranks = s.rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return pd.Series(qs.astype(str), index=s.index)


def triple_meta() -> pd.DataFrame:
    """PR #459 locked units: 21 + 22 + 9 = 52. Quartiles within-cohort %pos."""
    a = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    a = a.loc[a["origin"].isin(TUMOR_ORIGINS) & (a["n_malignant"] >= 20)].copy()
    a["patient"] = a["sample"].astype(str)
    a["cohort"] = "GSE131907"
    a["unit"] = "sample"
    a["malig_def"] = "author_malig"
    a["cldn4_pct"] = a["mal_CLDN4_pct"].astype(float)  # 0–100
    a["cldn4_mean"] = a["mal_CLDN4_mean"].astype(float)
    a["n_malignant"] = a["n_malignant"].astype(int)
    a["n_tnk"] = a["n_tnk"].astype(int)
    a["quartile"] = assign_quartiles(a.set_index("patient")["cldn4_pct"]).values
    a["in_malig_matrix"] = True

    b = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    b["patient"] = b["patient"].astype(str)
    b["cohort"] = "GSE205335"
    b["unit"] = "patient"
    b["malig_def"] = "author_malig"
    b["cldn4_pct"] = b["mal_CLDN4_pct_pos"].astype(float)  # 0–100
    b["cldn4_mean"] = b["mal_CLDN4_mean"].astype(float)
    b["n_malignant"] = b["n_malignant"].astype(int)
    b["n_tnk"] = b["n_tnk"].astype(int)
    b["quartile"] = assign_quartiles(b.set_index("patient")["cldn4_pct"]).values
    # P4001 is in the locked n=22 T/NK vector (27 malignant cells) but out of the
    # author-malignant UMI-sum matrix (n_mal≥30, same as PR #456).
    b["in_malig_matrix"] = b["n_malignant"] >= 30

    c = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    c = c.loc[c["eligible"].astype(str).str.lower() == "true"].copy()
    c["patient"] = c["patient"].astype(str)
    c["cohort"] = "GSE189357"
    c["unit"] = "patient"
    c["malig_def"] = "marker_malig"
    c["cldn4_pct"] = c["mal_CLDN4_pct"].astype(float) * 100.0  # fraction → percent
    c["cldn4_mean"] = c["mal_CLDN4_mean"].astype(float)
    c["n_malignant"] = c["n_malignant"].astype(int)
    c["n_tnk"] = c["n_tnk"].astype(int)
    c["quartile"] = assign_quartiles(c.set_index("patient")["cldn4_pct"]).values
    c["in_malig_matrix"] = True

    cols = [
        "patient", "cohort", "unit", "malig_def", "cldn4_pct", "cldn4_mean",
        "n_malignant", "n_tnk", "quartile", "in_malig_matrix",
    ]
    return pd.concat([a[cols], b[cols], c[cols]], ignore_index=True)


def read_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper()
    df.columns = df.columns.astype(str)
    return df.groupby(df.index).sum()


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def ols_de(logcpm: pd.DataFrame, design: pd.DataFrame, coef: str) -> pd.DataFrame:
    X = design.reindex(logcpm.columns).astype(float)
    if X.isna().any().any():
        keep = ~X.isna().any(axis=1)
        X = X.loc[keep]
        Y = logcpm.loc[:, keep].astype(float).to_numpy()
    else:
        Y = logcpm.astype(float).to_numpy()
    Xv = X.to_numpy()
    n, p = Xv.shape
    df_res = n - p
    if df_res < 1:
        return pd.DataFrame()
    xtx = Xv.T @ Xv
    try:
        xtx_inv = np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        xtx_inv = np.linalg.pinv(xtx)
    beta = xtx_inv @ Xv.T @ Y.T
    fitted = Xv @ beta
    resid = Y.T - fitted
    sse = np.sum(resid**2, axis=0)
    sigma2 = sse / df_res
    j = list(X.columns).index(coef)
    se = np.sqrt(np.maximum(sigma2 * xtx_inv[j, j], 0.0))
    est = beta[j]
    t = np.divide(est, se, out=np.zeros_like(est), where=se > 0)
    pval = 2.0 * stats.t.sf(np.abs(t), df_res)
    out = pd.DataFrame(
        {
            "gene": logcpm.index.astype(str),
            "logFC": est,
            "AveExpr": Y.mean(axis=1),
            "t": t,
            "p": pval,
            "se": se,
            "df": df_res,
        }
    )
    out["fdr"] = _bh(out["p"].values)
    return out.sort_values("p")


def _add_cohort_dummies(design: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Two dummies vs GSE131907 when ≥2 cohorts are present."""
    m = meta.set_index("patient")
    cohorts = m["cohort"]
    if cohorts.nunique() > 1:
        if (cohorts == "GSE205335").any():
            design["cohort_GSE205335"] = (cohorts.reindex(design.index) == "GSE205335").astype(float)
        if (cohorts == "GSE189357").any():
            design["cohort_GSE189357"] = (cohorts.reindex(design.index) == "GSE189357").astype(float)
    return design


def run_q4q1(counts: pd.DataFrame, meta: pd.DataFrame, cohort: str | None) -> tuple[pd.DataFrame, dict]:
    m = meta.copy()
    if cohort:
        m = m.loc[m["cohort"] == cohort]
    m = m.loc[m["quartile"].isin(["Q1", "Q4"])]
    m = m.loc[m["patient"].isin(counts.columns)].copy()
    n1 = int((m["quartile"] == "Q1").sum())
    n4 = int((m["quartile"] == "Q4").sum())
    info = {
        "n_q1": n1,
        "n_q4": n4,
        "n": n1 + n4,
        "patients_q1": ",".join(sorted(m.loc[m["quartile"] == "Q1", "patient"])),
        "patients_q4": ",".join(sorted(m.loc[m["quartile"] == "Q4", "patient"])),
        "thin": n1 < 6 or n4 < 6,
    }
    if n1 < 3 or n4 < 3:
        info["skip"] = f"thin tail n_Q1={n1} n_Q4={n4} (need ≥3 each)"
        return pd.DataFrame(), info
    cts = filter_genes(counts.loc[:, m["patient"]])
    fac = tmm_norm_factors(cts)
    lc = log_cpm(cts, fac)
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_Q4"] = (m.set_index("patient")["quartile"] == "Q4").astype(float)
    if cohort is None:
        design = _add_cohort_dummies(design, m)
    de = ols_de(lc, design, "CLDN4_Q4")
    info["n_genes"] = int(len(de))
    return de, info


def run_continuous(counts: pd.DataFrame, meta: pd.DataFrame, cohort: str | None = None) -> tuple[pd.DataFrame, dict]:
    m = meta.loc[meta["patient"].isin(counts.columns)].copy()
    if cohort:
        m = m.loc[m["cohort"] == cohort]
    info = {"n": int(len(m)), "n_q1": np.nan, "n_q4": np.nan}
    if len(m) < 8:
        info["skip"] = f"n={len(m)} < 8"
        return pd.DataFrame(), info
    cts = filter_genes(counts.loc[:, m["patient"]])
    fac = tmm_norm_factors(cts)
    lc = log_cpm(cts, fac)
    z = m.set_index("patient")["cldn4_pct"].astype(float)
    z = (z - z.mean()) / z.std(ddof=1)
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_pct_z"] = z
    if cohort is None:
        design = _add_cohort_dummies(design, m)
    de = ols_de(lc, design, "CLDN4_pct_z")
    info["n_genes"] = int(len(de))
    info["patients"] = ",".join(sorted(m["patient"]))
    return de, info


def family_score_de(
    logcpm: pd.DataFrame,
    meta: pd.DataFrame,
    fam: dict[str, set[str]],
    cohort: str | None,
    split: str,
) -> pd.DataFrame:
    """One OLS per family on the mean log2(CPM+1) of family genes."""
    m = meta.copy()
    if cohort:
        m = m.loc[m["cohort"] == cohort]
    if split == "q4q1":
        m = m.loc[m["quartile"].isin(["Q1", "Q4"])]
    m = m.loc[m["patient"].isin(logcpm.columns)].copy()
    if split == "q4q1" and ((m["quartile"] == "Q1").sum() < 3 or (m["quartile"] == "Q4").sum() < 3):
        return pd.DataFrame()
    if split == "continuous" and len(m) < 8:
        return pd.DataFrame()
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    if split == "q4q1":
        coef = "CLDN4_Q4"
        design[coef] = (m.set_index("patient")["quartile"] == "Q4").astype(float)
    else:
        coef = "CLDN4_pct_z"
        z = m.set_index("patient")["cldn4_pct"].astype(float)
        design[coef] = (z - z.mean()) / z.std(ddof=1)
    if cohort is None:
        design = _add_cohort_dummies(design, m)
    rows = []
    for name, genes in fam.items():
        present = [g for g in sorted(genes) if g in logcpm.index]
        if len(present) < 3:
            continue
        score = logcpm.loc[present, m["patient"]].astype(float).mean(axis=0).to_frame().T
        score.index = [name]
        de = ols_de(score, design, coef)
        if de.empty:
            continue
        r = de.iloc[0]
        rows.append(
            {
                "family": name,
                "split": split,
                "cohort": cohort or COMBO,
                "n_q1": int((m["quartile"] == "Q1").sum()) if split == "q4q1" else np.nan,
                "n_q4": int((m["quartile"] == "Q4").sum()) if split == "q4q1" else np.nan,
                "n": int(len(m)),
                "n_genes": int(len(present)),
                "logFC": float(r.logFC),
                "p": float(r.p),
                "t": float(r.t),
                "se": float(r.se),
                "df": float(r.df),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = _bh(out["p"].values)
    return out


def volcano(de: pd.DataFrame, fam: dict[str, set[str]], title: str, path: Path) -> None:
    if de.empty:
        return
    fig, ax = plt.subplots(figsize=(7.4, 5.5))
    x = de["logFC"].values
    y = -np.log10(np.clip(de["p"].values, 1e-300, 1))
    ax.scatter(x, y, s=8, c="#c8c8c8", linewidths=0, alpha=0.7, label="other")
    for fam_name, color in FAM_COLORS.items():
        sub = de[de["gene"].isin(fam[fam_name])]
        if sub.empty:
            continue
        ax.scatter(
            sub["logFC"],
            -np.log10(np.clip(sub["p"], 1e-300, 1)),
            s=22,
            c=color,
            linewidths=0.2,
            edgecolors="white",
            label=fam_name,
            zorder=3,
        )
    if "CLDN4" in set(de["gene"]):
        hit = de[de["gene"] == "CLDN4"]
        ax.scatter(
            hit["logFC"],
            -np.log10(np.clip(hit["p"], 1e-300, 1)),
            s=60,
            c="#111",
            marker="D",
            label="CLDN4 (held out)",
            zorder=4,
        )
    ax.axhline(-np.log10(0.05), ls="--", c="#666", lw=0.7)
    ax.axvline(0, ls="-", c="#999", lw=0.6)
    ax.set_xlabel("log2FC (CLDN4-high − low)")
    ax.set_ylabel("−log10 p (OLS)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def heatmap_family(logcpm: pd.DataFrame, meta: pd.DataFrame, genes: list[str], title: str, path: Path) -> None:
    genes = [g for g in genes if g in logcpm.index]
    if len(genes) < 4 or meta.empty:
        return
    m = meta.loc[meta["patient"].isin(logcpm.columns)].copy()
    m = m.sort_values(["cohort", "cldn4_pct"])
    mat = logcpm.loc[genes, m["patient"]].astype(float)
    z = mat.sub(mat.mean(axis=1), axis=0)
    sd = mat.std(axis=1, ddof=1).replace(0, np.nan)
    z = z.div(sd, axis=0)
    fig, ax = plt.subplots(figsize=(max(8.0, 0.26 * len(m) + 2.6), max(4.8, 0.24 * len(genes) + 1.6)))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2.2, vmax=2.2)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_xticks(range(len(m)))
    labels = [f"{r.patient}\n{r.quartile}" for r in m.itertuples()]
    ax.set_xticklabels(labels, fontsize=5.5, rotation=90)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="row z (log2 CPM+1)")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _boxplot(ax, q1, q4):
    ax.boxplot([q1, q4], tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"], widths=0.55)


def box_key_genes(logcpm: pd.DataFrame, meta: pd.DataFrame, title: str, path: Path) -> None:
    genes = [g for g in FOCAL_BOX if g in logcpm.index]
    if not genes:
        return
    m = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(logcpm.columns)].copy()
    n = len(genes)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 2.6 * nrows), squeeze=False)
    rng = np.random.default_rng(0)
    for i, gene in enumerate(genes):
        ax = axes[i // ncols][i % ncols]
        q1 = logcpm.loc[gene, m.loc[m["quartile"] == "Q1", "patient"]].astype(float)
        q4 = logcpm.loc[gene, m.loc[m["quartile"] == "Q4", "patient"]].astype(float)
        _boxplot(ax, q1, q4)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, len(q1)), q1, s=14, c="#4c78a8", zorder=3)
        ax.scatter(2 + rng.uniform(-0.08, 0.08, len(q4)), q4, s=14, c="#e45756", zorder=3)
        if len(q1) >= 3 and len(q4) >= 3:
            _, p = stats.mannwhitneyu(q4, q1, alternative="two-sided")
            ax.set_title(f"{gene}  p={p:.3g}", fontsize=9)
        else:
            ax.set_title(gene, fontsize=9)
        ax.set_ylabel("log2(CPM+1)")
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(title, y=1.01)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def box_family_scores(logcpm: pd.DataFrame, meta: pd.DataFrame, fam: dict[str, set[str]], title: str, path: Path) -> None:
    m = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(logcpm.columns)].copy()
    names = [k for k in FAM_COLORS if any(g in logcpm.index for g in fam[k])]
    if not names:
        return
    fig, axes = plt.subplots(1, len(names), figsize=(2.8 * len(names), 3.6), squeeze=False)
    rng = np.random.default_rng(1)
    for i, name in enumerate(names):
        ax = axes[0][i]
        present = [g for g in fam[name] if g in logcpm.index]
        score = logcpm.loc[present, m["patient"]].astype(float).mean(axis=0)
        q1 = score.loc[m.loc[m["quartile"] == "Q1", "patient"]]
        q4 = score.loc[m.loc[m["quartile"] == "Q4", "patient"]]
        _boxplot(ax, q1, q4)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, len(q1)), q1, s=16, c="#4c78a8", zorder=3)
        ax.scatter(2 + rng.uniform(-0.08, 0.08, len(q4)), q4, s=16, c="#e45756", zorder=3)
        if len(q1) >= 3 and len(q4) >= 3:
            _, p = stats.mannwhitneyu(q4, q1, alternative="two-sided")
            ax.set_title(f"{name}\np={p:.3g}", fontsize=9)
        else:
            ax.set_title(name, fontsize=9)
        ax.set_ylabel("mean log2(CPM+1)")
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def forest_families(de: pd.DataFrame, fam: dict[str, set[str]], title: str, path: Path) -> None:
    if de.empty:
        return
    rows = []
    for name, genes in fam.items():
        sub = de[de["gene"].isin(genes)]
        if sub.empty:
            continue
        for r in sub.nsmallest(8, "p").itertuples():
            rows.append({"family": name, "gene": r.gene, "logFC": r.logFC, "p": r.p, "fdr": r.fdr})
    if not rows:
        return
    plot = pd.DataFrame(rows).sort_values(["family", "p"])
    fig, ax = plt.subplots(figsize=(7.6, max(4.2, 0.28 * len(plot) + 1.2)))
    y = np.arange(len(plot))
    ax.axvline(0, c="#666", lw=0.7)
    ax.scatter(plot["logFC"], y, c=[FAM_COLORS[f] for f in plot["family"]], s=28, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.family}:{r.gene}" for r in plot.itertuples()], fontsize=7)
    ax.set_xlabel("log2FC (CLDN4-high − low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def forest_family_scores(fam_de: pd.DataFrame, title: str, path: Path) -> None:
    if fam_de.empty:
        return
    plot = fam_de.copy()
    plot["family"] = pd.Categorical(plot["family"], list(FAM_COLORS), ordered=True)
    plot = plot.sort_values("family")
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    y = np.arange(len(plot))
    ax.axvline(0, c="#666", lw=0.7)
    ax.errorbar(
        plot["logFC"],
        y,
        xerr=1.96 * plot["se"],
        fmt="o",
        color="#333",
        ecolor="#888",
        elinewidth=1.2,
        capsize=3,
    )
    ax.scatter(plot["logFC"], y, c=[FAM_COLORS[f] for f in plot["family"]], s=40, zorder=3)
    ax.set_yticks(y)
    labels = [f"{r.family}  n={int(r.n)}  p={r.p:.3g}" for r in plot.itertuples()]
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("family-score log2FC (CLDN4-high − low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def family_median_bar(de: pd.DataFrame, fam: dict[str, set[str]], title: str, path: Path) -> None:
    if de.empty:
        return
    names, meds, cols = [], [], []
    for name, genes in fam.items():
        sub = de[de["gene"].isin(genes)]
        if sub.empty:
            continue
        names.append(name)
        meds.append(float(sub["logFC"].median()))
        cols.append(FAM_COLORS[name])
    if not names:
        return
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.axvline(0, c="#666", lw=0.7)
    ax.barh(names, meds, color=cols, height=0.55)
    ax.set_xlabel("median log2FC (CLDN4-high − low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def strip_cldn4(meta: pd.DataFrame, path: Path) -> None:
    cohorts = ["GSE131907", "GSE205335", "GSE189357"]
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    cols = {"Q1": "#4c78a8", "Q2": "#9ecae1", "Q3": "#f4a582", "Q4": "#e45756"}
    for ax, cohort in zip(axes, cohorts):
        g = meta.loc[meta["cohort"] == cohort].sort_values("cldn4_pct")
        face = [cols[q] for q in g["quartile"]]
        edge = ["#111" if not flag else "white" for flag in g["in_malig_matrix"]]
        ax.scatter(
            range(len(g)),
            g["cldn4_pct"],
            c=face,
            s=36,
            edgecolors=edge,
            linewidths=0.6,
        )
        ax.set_title(f"{cohort}  n={len(g)}  unit={g['unit'].iloc[0]}")
        ax.set_ylabel("malignant CLDN4 %pos")
        ax.set_xlabel("units (sorted)")
        for q, c in cols.items():
            ax.scatter([], [], c=c, label=q, s=28)
        ax.legend(frameon=False, fontsize=7, ncol=4, loc="lower right")
    fig.suptitle("PR #459 triple CLDN4-only quartiles (within-cohort %pos). Black edge = out of malignant matrix (P4001).")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def n_bar(inv: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    labels = inv["contrast"].tolist()
    x = np.arange(len(labels))
    ax.bar(x - 0.18, inv["n_low"].fillna(0), 0.36, label="CLDN4-low / Q1", color="#4c78a8")
    ax.bar(x + 0.18, inv["n_high"].fillna(0), 0.36, label="CLDN4-high / Q4", color="#e45756")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("honest n (patients / samples)")
    ax.set_title("Malignant DE sample sizes — GSE131907+GSE205335+GSE189357 (no 148071)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (not np.isfinite(p))):
        return "NA"
    p = float(p)
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def _fam_summary(de: pd.DataFrame, fam: dict[str, set[str]], n_q1, n_q4) -> list[dict]:
    rows = []
    for name, genes in fam.items():
        sub = de[de["gene"].isin(genes)].copy()
        if sub.empty:
            rows.append({"family": name, "n_tested": 0, "n_p05": 0, "n_fdr05": 0, "n_up": 0, "n_down": 0})
            continue
        sig = sub[sub["p"] < 0.05]
        rows.append(
            {
                "family": name,
                "n_tested": int(len(sub)),
                "n_p05": int((sub["p"] < 0.05).sum()),
                "n_fdr05": int((sub["fdr"] < 0.05).sum()),
                "n_up": int((sig["logFC"] > 0).sum()),
                "n_down": int((sig["logFC"] < 0).sum()),
                "median_logFC": float(sub["logFC"].median()),
                "top_gene": str(sub.iloc[0]["gene"]),
                "top_logFC": float(sub.iloc[0]["logFC"]),
                "top_p": float(sub.iloc[0]["p"]),
                "top_fdr": float(sub.iloc[0]["fdr"]),
                "n_q1": n_q1,
                "n_q4": n_q4,
            }
        )
    return rows


def write_finding(
    meta: pd.DataFrame,
    inv: pd.DataFrame,
    fam_rows: pd.DataFrame,
    fam_score: pd.DataFrame,
    headline: dict,
    cldn4_row: dict,
) -> None:
    def row(contrast: str) -> pd.Series:
        hit = inv[inv["contrast"] == contrast]
        return hit.iloc[0] if len(hit) else pd.Series(dtype=object)

    mal = row("malignant_q4q1_combined")
    mal_c = row("malignant_continuous_combined")

    def fam_md(contrast: str) -> str:
        sub = fam_rows[fam_rows["contrast"] == contrast] if not fam_rows.empty else pd.DataFrame()
        if sub.empty:
            return "_No family rows (contrast not run or tail too thin)._"
        lines = [
            "| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |",
            "|---|---:|---|---:|---:|---|",
        ]
        for r in sub.itertuples():
            med = getattr(r, "median_logFC", float("nan"))
            med_s = f"{med:+.3f}" if np.isfinite(med) else "NA"
            top = getattr(r, "top_gene", "")
            if top:
                lines.append(
                    f"| {r.family} | {r.n_tested} | {r.n_p05} ({r.n_up}/{r.n_down}) | {r.n_fdr05} | "
                    f"{med_s} | {top} ({r.top_logFC:+.3f}, {_fmt_p(r.top_p)}, {_fmt_p(r.top_fdr)}) |"
                )
            else:
                lines.append(f"| {r.family} | 0 | 0 (0/0) | 0 | NA | — |")
        return "\n".join(lines)

    def score_md(split: str, cohort: str = COMBO) -> str:
        if fam_score.empty:
            return "_No family-score DE._"
        sub = fam_score[(fam_score["split"] == split) & (fam_score["cohort"] == cohort)]
        if sub.empty:
            return "_No family-score DE._"
        lines = [
            "| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |",
            "|---|---:|---|---:|---:|---:|---:|",
        ]
        for r in sub.itertuples():
            nq = "—" if split != "q4q1" else f"{int(r.n_q1)} / {int(r.n_q4)}"
            lines.append(
                f"| {r.family} | {int(r.n)} | {nq} | {int(r.n_genes)} | "
                f"{r.logFC:+.3f} | {_fmt_p(r.p)} | {_fmt_p(r.fdr)} |"
            )
        return "\n".join(lines)

    def score_md_cohorts(split: str) -> str:
        if fam_score.empty:
            return "_No per-cohort family-score DE._"
        sub = fam_score[(fam_score["split"] == split) & (fam_score["cohort"] != COMBO)]
        if sub.empty:
            return "_No per-cohort family-score DE._"
        lines = [
            "| cohort | family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |",
            "|---|---|---:|---|---:|---:|---:|---:|",
        ]
        for r in sub.itertuples():
            nq = "—" if split != "q4q1" else f"{int(r.n_q1)} / {int(r.n_q4)}"
            lines.append(
                f"| {r.cohort} | {r.family} | {int(r.n)} | {nq} | {int(r.n_genes)} | "
                f"{r.logFC:+.3f} | {_fmt_p(r.p)} | {_fmt_p(r.fdr)} |"
            )
        return "\n".join(lines)

    counts = meta["cohort"].value_counts()
    qtab = {c: meta.loc[meta["cohort"] == c, "quartile"].value_counts() for c in counts.index}
    q1 = meta[meta["quartile"] == "Q1"]
    q4 = meta[meta["quartile"] == "Q4"]
    dropped = meta.loc[~meta["in_malig_matrix"]]
    drop_txt = ", ".join(
        f"{r.patient} ({r.cohort}, {int(r.n_malignant)} malignant cells, {r.quartile})"
        for r in dropped.itertuples()
    ) or "none"

    cldn4_txt = "CLDN4 not in the combined matrix."
    if cldn4_row:
        cldn4_txt = (
            f"CLDN4 itself {cldn4_row['logFC']:+.2f} "
            f"(p={_fmt_p(cldn4_row.get('p'))}, FDR={_fmt_p(cldn4_row.get('fdr', float('nan')))}) "
            f"— direction check on the split gene (held out of the TJ family). "
            f"n_Q1={cldn4_row.get('n_q1', 'NA')} n_Q4={cldn4_row.get('n_q4', 'NA')}."
        )

    skip189 = row("malignant_q4q1_GSE189357").get("note", "")

    text = f"""# Triple GSE131907+GSE205335+GSE189357: CLDN4-only malignant IFN/MHC/TJ DE

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. **Not infiltrate.**
Tumor-cell-intrinsic program only. No GSE148071.

The pair GSE123902+GSE189357 %pos n=22 ρ=−0.638 (PR #459) already has a
malignant IFN/MHC/TJ DE (`methods/pair_123902_189357_malig_ifn_de_cldn4`,
PR #469). This extra is the remaining 189357-axis triple that differs.

Given triple cut (PR #459, **not re-audited**): GSE131907 + GSE205335 +
GSE189357 malignant CLDN4 %pos vs same-unit T/NK, n=52, Spearman ρ=−0.497
(p=0.00035). This extra does **not** re-audit that T/NK ρ.

**Thesis (already correct):** CLDN4-high malignant cells should be IFN/MHC-I
down, TJ up on the malignant cell itself.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
GSE131907 is sample-level. GSE205335 / GSE189357 are patient-level.
p-values are descriptive.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE131907 + GSE205335 + GSE189357 | Pair 123902+189357 already filled; **no GSE148071** |
| Units | PR #459 %pos n=52 (21 samples + 22 patients + 9 patients) | Mixed units (sample vs patient); do not quote a cell-level n |
| Split | Within-cohort marker/author-malignant CLDN4 **%pos Q4 vs Q1** | Quartile tails; mid quartiles unused in binary DE |
| Continuous | CLDN4 %pos z, units in the malignant matrix, cohort covariates | Linear; not a causal model |
| Malignant | GSE131907/GSE205335 author-malig UMI-sum; GSE189357 marker-malig UMI-sum | Marker gate ≠ author annotation |
| Dropped from DE | P4001 (27 malignant cells) | In the locked n=52 T/NK vector; out of malignant matrix (n_mal≥30) |
| T/NK | **not a DE compartment** | T/NK ρ taken as given from PR #459; not re-audited |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (single) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I/APM (custom), TJ (KEGG/GO + junction focal; **CLDN4 held out**) | MHC-II and chemokine panels are not the claim |

## Honest n

PR #459 locked labels (malignant CLDN4 %pos; T/NK ρ not re-scored):

- GSE131907: n={int(counts.get('GSE131907', 0))} tumor-origin samples (author malignant, n_mal≥20). Q1={int(qtab.get('GSE131907', pd.Series()).get('Q1', 0))} Q4={int(qtab.get('GSE131907', pd.Series()).get('Q4', 0))}.
- GSE205335: n={int(counts.get('GSE205335', 0))} patients (author malignant). Q1={int(qtab.get('GSE205335', pd.Series()).get('Q1', 0))} Q4={int(qtab.get('GSE205335', pd.Series()).get('Q4', 0))}.
- GSE189357: n={int(counts.get('GSE189357', 0))} patients (marker-malignant). Q1={int(qtab.get('GSE189357', pd.Series()).get('Q1', 0))} Q4={int(qtab.get('GSE189357', pd.Series()).get('Q4', 0))}.

Locked Q1: {", ".join(sorted(q1["patient"]))} (n={len(q1)}).
Locked Q4: {", ".join(sorted(q4["patient"]))} (n={len(q4)}).

Out of malignant DE: {drop_txt}.

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 combined | {int(mal.get('n_low', 0) or 0)}/{int(mal.get('n_high', 0) or 0)} | {int(mal.get('n_genes', 0) or 0)} | cohort covariates; P4001 out |
| Q4 vs Q1 GSE131907 | {int(row('malignant_q4q1_GSE131907').get('n_low', 0) or 0)}/{int(row('malignant_q4q1_GSE131907').get('n_high', 0) or 0)} | {int(row('malignant_q4q1_GSE131907').get('n_genes', 0) or 0)} | sample-level |
| Q4 vs Q1 GSE205335 | {int(row('malignant_q4q1_GSE205335').get('n_low', 0) or 0)}/{int(row('malignant_q4q1_GSE205335').get('n_high', 0) or 0)} | {int(row('malignant_q4q1_GSE205335').get('n_genes', 0) or 0)} | P4001 out of matrix |
| Q4 vs Q1 GSE189357 | {int(row('malignant_q4q1_GSE189357').get('n_low', 0) or 0)}/{int(row('malignant_q4q1_GSE189357').get('n_high', 0) or 0)} | {int(row('malignant_q4q1_GSE189357').get('n_genes', 0) or 0)} | {skip189 or "thin if Q4<3"} |
| continuous combined | n={int(mal_c.get('n', 0) or 0)} | {int(mal_c.get('n_genes', 0) or 0)} | CLDN4 %pos z; units in the matrix (51 if P4001 out) |

Do not quote n=52 for the malignant DE. n=52 is the PR #459 T/NK cut.
Malignant DE n is the units actually in the UMI-sum matrices.

## Family DE (combined Q4 vs Q1) — n / logFC / p

Family score = mean log2(TMM-CPM+1) of family genes present. Positive logFC =
higher in CLDN4-high (Q4) than CLDN4-low (Q1). Machine table:
`tables/family_de.tsv`.

Expected under the thesis: IFN down, MHC-I/APM down, TJ up.

{score_md("q4q1")}

Continuous family-score DE (units in the malignant matrix, CLDN4 %pos z):

{score_md("continuous")}

Per-cohort family scores (GSE189357 Q4 vs Q1 is thin and may be skipped):

{score_md_cohorts("q4q1")}

{cldn4_txt}

## Gene-level family members (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high.

{fam_md("q4q1_combined")}

Headline family genes (combined Q4 vs Q1, lowest p):

{headline.get("md", "_none_")}

### Continuous CLDN4 %pos (cohort covariates)

{fam_md("continuous_combined")}

## What this is not

- Not infiltrate / T/NK fraction DE and not a re-audit of PR #459 ρ.
- Not CellChat / LIANA / NicheNet.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not GSE148071 and not the already-filled GSE123902+GSE189357 pair DE.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ change.
- Genome-wide FDR on these n is expected to be thin; family scores and
  family median logFC / sign counts are the claim.

## Files

- `tables/family_de.tsv` — **headline family DE** (n / logFC / p)
- `tables/de_q4q1_combined_families.tsv` — gene-level family members
- `tables/family_summary.tsv` — IFN / MHC-I/APM / TJ counts
- `tables/de_families.tsv` — family rows, all contrasts
- `tables/de_all.tsv` — all genes, all contrasts
- `tables/n_honest.tsv`
- `tables/sample_inventory.tsv`
- `figures/` — volcano, family heatmap, key-gene boxes, family-score boxes, forests, family-median bar, CLDN4 strip, n bars

Reproduce:

```bash
python3 methods/pair_or_triple_189357_malig_ifn_de_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(text)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    fam = load_families()
    meta = triple_meta()
    meta.to_csv(TABLES / "sample_inventory.tsv", sep="\t", index=False)
    strip_cldn4(meta, FIGS / "cldn4_quartile_strip")

    paths = {
        "GSE131907": DATA / "GSE131907_malignant_counts.tsv.gz",
        "GSE205335": DATA / "GSE205335_malignant_counts.tsv.gz",
        "GSE189357": DATA / "GSE189357_malignant_counts.tsv.gz",
    }
    parts = {}
    all_de = []
    inv_rows = []
    fam_rows = []
    fam_score_rows = []

    for cohort, path in paths.items():
        if not path.exists():
            raise SystemExit(f"missing {path}")
        cts = read_counts(path)
        parts[cohort] = cts
        de, info = run_q4q1(cts, meta, cohort)
        inv_rows.append(
            {
                "contrast": f"malignant_q4q1_{cohort}",
                "compartment": "malignant",
                "cohort": cohort,
                "split": "q4q1",
                "n_low": info["n_q1"],
                "n_high": info["n_q4"],
                "n": info["n"],
                "n_genes": info.get("n_genes", 0),
                "patients_q1": info.get("patients_q1", ""),
                "patients_q4": info.get("patients_q4", ""),
                "note": info.get("skip", "thin tails" if info.get("thin") else ""),
            }
        )
        if not de.empty:
            de = de.copy()
            de["compartment"] = "malignant"
            de["cohort"] = cohort
            de["contrast"] = f"q4q1_{cohort}"
            de["n_q1"] = info["n_q1"]
            de["n_q4"] = info["n_q4"]
            de["family"] = de["gene"].map(lambda g: family_of(g, fam))
            all_de.append(de)
            for r in _fam_summary(de, fam, info["n_q1"], info["n_q4"]):
                fam_rows.append({**r, "compartment": "malignant", "contrast": f"q4q1_{cohort}"})
        de_c, info_c = run_continuous(cts, meta, cohort)
        inv_rows.append(
            {
                "contrast": f"malignant_continuous_{cohort}",
                "compartment": "malignant",
                "cohort": cohort,
                "split": "continuous",
                "n_low": np.nan,
                "n_high": np.nan,
                "n": info_c["n"],
                "n_genes": info_c.get("n_genes", 0),
                "patients_q1": "",
                "patients_q4": "",
                "note": info_c.get("skip", ""),
            }
        )
        if not de_c.empty:
            de_c = de_c.copy()
            de_c["compartment"] = "malignant"
            de_c["cohort"] = cohort
            de_c["contrast"] = f"continuous_{cohort}"
            de_c["n_q1"] = np.nan
            de_c["n_q4"] = np.nan
            de_c["family"] = de_c["gene"].map(lambda g: family_of(g, fam))
            all_de.append(de_c)
            for r in _fam_summary(de_c, fam, np.nan, np.nan):
                fam_rows.append({**r, "compartment": "malignant", "contrast": f"continuous_{cohort}"})

    genes = sorted(set.intersection(*[set(p.index) for p in parts.values()]))
    combined = pd.concat([p.reindex(genes).fillna(0) for p in parts.values()], axis=1)
    combined = combined.loc[:, ~combined.columns.duplicated()]

    de, info = run_q4q1(combined, meta, None)
    inv_rows.append(
        {
            "contrast": "malignant_q4q1_combined",
            "compartment": "malignant",
            "cohort": COMBO,
            "split": "q4q1",
            "n_low": info["n_q1"],
            "n_high": info["n_q4"],
            "n": info["n"],
            "n_genes": info.get("n_genes", 0),
            "patients_q1": info.get("patients_q1", ""),
            "patients_q4": info.get("patients_q4", ""),
            "note": "cohort covariates; P4001 out of malignant matrix",
        }
    )
    de_c, info_c = run_continuous(combined, meta, None)
    inv_rows.append(
        {
            "contrast": "malignant_continuous_combined",
            "compartment": "malignant",
            "cohort": COMBO,
            "split": "continuous",
            "n_low": np.nan,
            "n_high": np.nan,
            "n": info_c["n"],
            "n_genes": info_c.get("n_genes", 0),
            "patients_q1": "",
            "patients_q4": "",
            "note": f"CLDN4 %pos z; n={info_c['n']} in matrix (locked T/NK n=52)",
        }
    )

    cldn4_row = {}
    if not de.empty:
        de = de.copy()
        de["compartment"] = "malignant"
        de["cohort"] = COMBO
        de["contrast"] = "q4q1_combined"
        de["n_q1"] = info["n_q1"]
        de["n_q4"] = info["n_q4"]
        de["family"] = de["gene"].map(lambda g: family_of(g, fam))
        all_de.append(de)
        for r in _fam_summary(de, fam, info["n_q1"], info["n_q4"]):
            fam_rows.append({**r, "compartment": "malignant", "contrast": "q4q1_combined"})
        hit = de[de["gene"] == "CLDN4"]
        if len(hit):
            cldn4_row = {
                "logFC": float(hit.iloc[0]["logFC"]),
                "p": float(hit.iloc[0]["p"]),
                "fdr": float(hit.iloc[0]["fdr"]),
                "n_q1": info["n_q1"],
                "n_q4": info["n_q4"],
            }
        volcano(
            de,
            fam,
            f"Malignant Q4 vs Q1  n={info['n_q4']} vs {info['n_q1']}  ({COMBO}; cohort covariates)",
            FIGS / "volcano_malignant_q4q1_combined",
        )
        forest_families(
            de,
            fam,
            "Malignant family genes, Q4 vs Q1 combined (CLDN4 held out of TJ)",
            FIGS / "forest_malignant_families_q4q1",
        )
        family_median_bar(
            de,
            fam,
            "Family median logFC, Q4 vs Q1 (thesis: IFN/MHC down, TJ up)",
            FIGS / "bar_family_median_logfc",
        )

    if not de_c.empty:
        de_c = de_c.copy()
        de_c["compartment"] = "malignant"
        de_c["cohort"] = COMBO
        de_c["contrast"] = "continuous_combined"
        de_c["n_q1"] = np.nan
        de_c["n_q4"] = np.nan
        de_c["family"] = de_c["gene"].map(lambda g: family_of(g, fam))
        all_de.append(de_c)
        for r in _fam_summary(de_c, fam, np.nan, np.nan):
            fam_rows.append({**r, "compartment": "malignant", "contrast": "continuous_combined"})

    m_all = meta.loc[meta["patient"].isin(combined.columns)].copy()
    cts_all = filter_genes(combined.loc[:, m_all["patient"]])
    lc_all = log_cpm(cts_all, tmm_norm_factors(cts_all))
    for cohort in (None, "GSE131907", "GSE205335", "GSE189357"):
        for split in ("q4q1", "continuous"):
            fs = family_score_de(lc_all, meta, fam, cohort, split)
            if not fs.empty:
                fam_score_rows.append(fs)

    m_q = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(combined.columns)]
    if len(m_q) >= 6:
        cts_q = filter_genes(combined.loc[:, m_q["patient"]])
        lc = log_cpm(cts_q, tmm_norm_factors(cts_q))
        show = []
        prefer = {
            "IFN": ["STAT1", "IRF1", "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IFI6", "BST2", "GBP1"],
            "MHC-I/APM": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5", "PSMB8", "PSMB9", "ERAP1"],
            "TJ": ["CLDN1", "CLDN7", "OCLN", "TJP1", "F11R", "PARD3", "CDH1", "TJP2"],
        }
        for name in ("IFN", "MHC-I/APM", "TJ"):
            present = [g for g in prefer[name] if g in lc.index]
            show.extend(present)
        if "CLDN4" in lc.index:
            show = ["CLDN4"] + show
        heatmap_family(
            lc,
            m_q,
            show,
            "Malignant family genes, Q1/Q4 only (row z; CLDN4 shown, held out of TJ test)",
            FIGS / "heatmap_malignant_families_q4q1",
        )
        box_key_genes(
            lc,
            m_q,
            f"Malignant key genes, CLDN4 Q4 vs Q1 (n={info.get('n_q4', '?')} vs {info.get('n_q1', '?')})",
            FIGS / "box_malignant_key_genes",
        )
        box_family_scores(
            lc,
            m_q,
            fam,
            "Malignant family scores, CLDN4 Q4 vs Q1",
            FIGS / "box_malignant_family_scores",
        )

    inv = pd.DataFrame(inv_rows)
    inv.to_csv(TABLES / "n_honest.tsv", sep="\t", index=False)
    n_bar(inv[inv["split"] == "q4q1"], FIGS / "n_honest_q4q1")

    fam_score = pd.concat(fam_score_rows, ignore_index=True) if fam_score_rows else pd.DataFrame()
    if not fam_score.empty:
        fam_score.to_csv(TABLES / "family_de.tsv", sep="\t", index=False)
        head_fs = fam_score[(fam_score["cohort"] == COMBO) & (fam_score["split"] == "q4q1")]
        forest_family_scores(
            head_fs,
            "Malignant family scores, Q4 vs Q1 combined",
            FIGS / "forest_family_scores_q4q1",
        )

    headline = {"md": "_DE table empty._"}
    if all_de:
        de_all = pd.concat(all_de, ignore_index=True)
        de_all.to_csv(TABLES / "de_all.tsv", sep="\t", index=False)
        fam_de = de_all[de_all["family"] != "other"].copy()
        fam_de.to_csv(TABLES / "de_families.tsv", sep="\t", index=False)
        head = fam_de[fam_de["contrast"] == "q4q1_combined"].sort_values("p")
        head.to_csv(TABLES / "de_q4q1_combined_families.tsv", sep="\t", index=False)
        lines = [
            "| family | gene | n_Q1 | n_Q4 | logFC | p | FDR |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
        for r in head.head(12).itertuples():
            lines.append(
                f"| {r.family} | {r.gene} | {int(r.n_q1)} | {int(r.n_q4)} | "
                f"{r.logFC:+.3f} | {_fmt_p(r.p)} | {_fmt_p(r.fdr)} |"
            )
        headline = {"md": "\n".join(lines)}

    fam_df = pd.DataFrame(fam_rows)
    if not fam_df.empty:
        fam_df.to_csv(TABLES / "family_summary.tsv", sep="\t", index=False)
    write_finding(
        meta,
        inv,
        fam_df if not fam_df.empty else pd.DataFrame(),
        fam_score if not fam_score.empty else pd.DataFrame(),
        headline,
        cldn4_row,
    )
    print("wrote FINDING.md and tables/")
    print(inv.to_string(index=False))
    if not fam_score.empty:
        print(fam_score.to_string(index=False))
    if not fam_df.empty:
        print(fam_df[fam_df["contrast"] == "q4q1_combined"].to_string(index=False))


if __name__ == "__main__":
    main()
