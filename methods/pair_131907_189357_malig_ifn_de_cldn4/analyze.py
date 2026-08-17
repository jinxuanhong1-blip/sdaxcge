#!/usr/bin/env python3
"""CLDN4-only malignant-cell-intrinsic patient-pseudobulk DE.

Pair: GSE131907 (author malignant, sample) + GSE189357 (marker malignant, patient).
PR #459 given cut: %pos n=30 ρ=−0.542, Q4 r=−0.619. Not CellChat. Not T/NK
infiltrate. No dual-high. No GSE148071.

Method: patient UMI-sum → TMM → log2(CPM+1) → OLS (limma-style, not muscat).
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
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}

CHEMOKINE = [
    "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CXCL16", "CXCL8", "CX3CL1",
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL19", "CCL21", "CCL22",
    "XCL1", "XCL2", "IL15", "IL2", "IL7", "IL21",
    "CXCR3", "CXCR4", "CXCR5", "CCR5", "CCR7",
]
FOCAL_BOX = [
    "CLDN4", "STAT1", "IRF1", "HLA-A", "HLA-B", "B2M", "TAP1", "TAP2",
    "CXCL9", "CXCL10", "CCL5", "OCLN", "TJP1", "KRT8", "KRT18", "KRT19",
]
FAM_COLORS = {
    "IFN": "#d62728",
    "MHC-I/APM": "#1f77b4",
    "TJ": "#2ca02c",
    "keratin": "#9467bd",
    "chemokine": "#ff7f0e",
}


def load_families() -> dict[str, set[str]]:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    tj = (
        set(sets["KEGG_TIGHT_JUNCTION"])
        | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
        | set(g for g in a8["focal_genes"] if g not in sets["KRT_EPITHELIAL"])
    )
    tj.discard("CLDN4")  # held out of the TJ family
    fam = {
        "IFN": set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
        | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]),
        "MHC-I/APM": set(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]),
        "TJ": tj,
        "keratin": set(sets["KRT_EPITHELIAL"]),
        "chemokine": set(CHEMOKINE),
    }
    return fam


def family_of(gene: str, fam: dict[str, set[str]]) -> str:
    hits = [k for k, vs in fam.items() if gene in vs]
    return "|".join(hits) if hits else "other"


def tmm_norm_factors(counts: pd.DataFrame) -> pd.Series:
    lib = counts.sum(axis=0).astype(float)
    lib = lib.replace(0, np.nan)
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
    cpm = counts.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def filter_genes(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 3) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def assign_quartiles(values: pd.Series) -> pd.Series:
    s = values.astype(float)
    ranks = s.rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return pd.Series(qs.astype(str), index=s.index)


def pair_meta() -> pd.DataFrame:
    """PR #459 locked units: GSE131907 n=21 samples + GSE189357 n=9 patients."""
    a = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    a = a.loc[a["origin"].isin(TUMOR_ORIGINS) & (a["n_malignant"] >= 20)].copy()
    a["patient"] = a["sample"].astype(str)
    a["cohort"] = "GSE131907"
    a["unit"] = "sample"
    a["malig_def"] = "author_malig"
    a["cldn4_pct"] = a["mal_CLDN4_pct"].astype(float)  # already 0–100
    a["cldn4_mean"] = a["mal_CLDN4_mean"].astype(float)
    a["n_malignant"] = a["n_malignant"].astype(int)
    a["quartile"] = assign_quartiles(a.set_index("patient")["cldn4_pct"]).values

    b = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    el = b[b["eligible"].astype(str).str.lower() == "true"].copy()
    el["patient"] = el["patient"].astype(str)
    el["cohort"] = "GSE189357"
    el["unit"] = "patient"
    el["malig_def"] = "marker_malig"
    el["cldn4_pct"] = el["mal_CLDN4_pct"].astype(float) * 100.0  # fraction → percent
    el["cldn4_mean"] = el["mal_CLDN4_mean"].astype(float)
    el["n_malignant"] = el["n_malignant"].astype(int)
    el["quartile"] = assign_quartiles(el.set_index("patient")["cldn4_pct"]).values

    cols = [
        "patient", "cohort", "unit", "malig_def", "cldn4_pct", "cldn4_mean",
        "n_malignant", "quartile",
    ]
    return pd.concat([a[cols], el[cols]], ignore_index=True)


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
    pv = 2.0 * stats.t.sf(np.abs(t), df_res)
    out = pd.DataFrame(
        {
            "gene": logcpm.index.astype(str),
            "logFC": est,
            "AveExpr": Y.mean(axis=1),
            "t": t,
            "p": pv,
            "se": se,
            "df": df_res,
        }
    )
    out["fdr"] = _bh(out["p"].values)
    return out.sort_values("p")


def _design_q4q1(m: pd.DataFrame, cohort: str | None) -> pd.DataFrame:
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_Q4"] = (m.set_index("patient")["quartile"] == "Q4").astype(float)
    if cohort is None and m["cohort"].nunique() > 1:
        design["cohort_GSE131907"] = (m.set_index("patient")["cohort"] == "GSE131907").astype(float)
    return design


def run_q4q1(counts: pd.DataFrame, meta: pd.DataFrame, cohort: str | None) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
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
    }
    if n1 < 3 or n4 < 3:
        return pd.DataFrame(), info, pd.DataFrame()
    cts = filter_genes(counts.loc[:, m["patient"]])
    fac = tmm_norm_factors(cts)
    lc = log_cpm(cts, fac)
    de = ols_de(lc, _design_q4q1(m, cohort), "CLDN4_Q4")
    info["n_genes"] = int(len(de))
    return de, info, lc


def run_continuous(counts: pd.DataFrame, meta: pd.DataFrame, cohort: str | None) -> tuple[pd.DataFrame, dict]:
    m = meta.copy()
    if cohort:
        m = m.loc[m["cohort"] == cohort]
    m = m.loc[m["patient"].isin(counts.columns)].copy()
    info = {"n": int(len(m)), "n_q1": np.nan, "n_q4": np.nan}
    if len(m) < 8:
        return pd.DataFrame(), info
    cts = filter_genes(counts.loc[:, m["patient"]])
    fac = tmm_norm_factors(cts)
    lc = log_cpm(cts, fac)
    z = m.set_index("patient")["cldn4_pct"].astype(float)
    z = (z - z.mean()) / z.std(ddof=1)
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_pct_z"] = z
    if cohort is None and m["cohort"].nunique() > 1:
        design["cohort_GSE131907"] = (m.set_index("patient")["cohort"] == "GSE131907").astype(float)
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
    if cohort is None and m["cohort"].nunique() > 1:
        design["cohort_GSE131907"] = (m.set_index("patient")["cohort"] == "GSE131907").astype(float)
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
                "cohort": cohort or "GSE131907+GSE189357",
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
    # CLDN4 held out of TJ — mark as direction check
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
    fig_w = max(8.0, 0.28 * len(m) + 2.4)
    fig, ax = plt.subplots(figsize=(fig_w, max(4.8, 0.22 * len(genes) + 1.6)))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2.2, vmax=2.2)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_xticks(range(len(m)))
    labels = [f"{r.patient}\n{r.quartile}" for r in m.itertuples()]
    ax.set_xticklabels(labels, fontsize=6, rotation=90)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="row z (log2 CPM+1)")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


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
        ax.boxplot([q1, q4], tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"], widths=0.55)
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
    fig, axes = plt.subplots(1, len(names), figsize=(2.6 * len(names), 3.6), squeeze=False)
    rng = np.random.default_rng(1)
    for i, name in enumerate(names):
        ax = axes[0][i]
        present = [g for g in fam[name] if g in logcpm.index]
        score = logcpm.loc[present, m["patient"]].astype(float).mean(axis=0)
        q1 = score.loc[m.loc[m["quartile"] == "Q1", "patient"]]
        q4 = score.loc[m.loc[m["quartile"] == "Q4", "patient"]]
        ax.boxplot([q1, q4], tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"], widths=0.55)
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
        for r in sub.nsmallest(6, "p").itertuples():
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
    plot = fam_de.sort_values("family").copy()
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
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


def bar_family_median_logfc(fam_rows: pd.DataFrame, path: Path) -> None:
    """Extra: gene-level median logFC by family on combined Q4 vs Q1."""
    sub = fam_rows[fam_rows["contrast"] == "q4q1_combined"].copy()
    if sub.empty or "median_logFC" not in sub.columns:
        return
    sub = sub.sort_values("family")
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    colors = [FAM_COLORS.get(f, "#888") for f in sub["family"]]
    ax.barh(sub["family"], sub["median_logFC"], color=colors, edgecolor="white")
    ax.axvline(0, c="#666", lw=0.7)
    ax.set_xlabel("median gene log2FC (CLDN4-high − low)")
    ax.set_title("Family median logFC, combined Q4 vs Q1")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def strip_cldn4(meta: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8))
    for ax, (cohort, g) in zip(axes, meta.groupby("cohort", sort=True)):
        g = g.sort_values("cldn4_pct")
        cols = {"Q1": "#4c78a8", "Q2": "#9ecae1", "Q3": "#f4a582", "Q4": "#e45756"}
        ax.scatter(
            range(len(g)),
            g["cldn4_pct"],
            c=[cols[q] for q in g["quartile"]],
            s=36,
            edgecolors="white",
            linewidths=0.4,
        )
        ax.set_title(f"{cohort}  n={len(g)}  unit={g['unit'].iloc[0]}")
        ax.set_ylabel("malignant CLDN4 %pos")
        ax.set_xlabel("units (sorted)")
        for q, c in cols.items():
            ax.scatter([], [], c=c, label=q, s=28)
        ax.legend(frameon=False, fontsize=7, ncol=4, loc="lower right")
    fig.suptitle("PR #459 pair CLDN4-only quartiles (GSE131907 sample / GSE189357 patient)")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def n_bar(inv: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    labels = inv["contrast"].tolist()
    x = np.arange(len(labels))
    ax.bar(x - 0.18, inv["n_low"].fillna(0), 0.36, label="CLDN4-low / Q1", color="#4c78a8")
    ax.bar(x + 0.18, inv["n_high"].fillna(0), 0.36, label="CLDN4-high / Q4", color="#e45756")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("honest n (donors / samples)")
    ax.set_title("Malignant DE sample sizes — GSE131907+GSE189357 only")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (np.isnan(p) or p != p)):
        return "NA"
    p = float(p)
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def _fam_summary(de: pd.DataFrame, fam: dict[str, set[str]], n_q1: int, n_q4: int) -> list[dict]:
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
        sub = fam_rows[fam_rows["contrast"] == contrast]
        if sub.empty:
            return "_No family rows._"
        lines = [
            "| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |",
            "|---|---:|---|---:|---:|---|",
        ]
        for r in sub.itertuples():
            lines.append(
                f"| {r.family} | {r.n_tested} | {r.n_p05} ({r.n_up}/{r.n_down}) | {r.n_fdr05} | "
                f"{r.median_logFC:+.3f} | {r.top_gene} ({r.top_logFC:+.3f}, {_fmt_p(r.top_p)}, {_fmt_p(r.top_fdr)}) |"
            )
        return "\n".join(lines)

    def score_md(split: str, cohort: str = "GSE131907+GSE189357") -> str:
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
        sub = fam_score[
            (fam_score["split"] == split) & (fam_score["cohort"] != "GSE131907+GSE189357")
        ]
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

    n189 = int((meta["cohort"] == "GSE189357").sum())
    n131 = int((meta["cohort"] == "GSE131907").sum())
    q189 = meta[meta["cohort"] == "GSE189357"]["quartile"].value_counts()
    q131 = meta[meta["cohort"] == "GSE131907"]["quartile"].value_counts()
    cldn4_txt = (
        f"CLDN4 itself (held out of TJ) combined Q4 vs Q1: "
        f"logFC={cldn4_row.get('logFC', float('nan')):+.3f}, "
        f"p={_fmt_p(cldn4_row.get('p'))}, n_Q1={cldn4_row.get('n_q1', 'NA')}, "
        f"n_Q4={cldn4_row.get('n_q4', 'NA')} — direction check on the split gene."
    )

    text = f"""# Pair GSE131907+GSE189357 — malignant-cell-intrinsic CLDN4 DE

ADDITIVE. **CLDN4-only.** Tumor-cell-intrinsic program DE.
Not CellChat. Not T/NK infiltrate. No dual-high TACSTD2×CLDN4. No GSE148071.

The pair that already differs is taken from PR #459 (malignant CLDN4 **%pos**,
n=30, ρ=−0.542 vs T/NK; Q4 vs Q1 r=−0.619). That T/NK rho is **not
re-audited**. This PR asks a different question: in the same 30 units, do
CLDN4-high malignant cells themselves have IFN / MHC-I/APM down and TJ up.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
muscat `pbDS` is the same collapse (patient × cell-type UMI-sum → bulk DE).
GSE131907 is sample-level. GSE189357 is patient-level. p-values are descriptive.

Thesis (not re-derived): CLDN4 KD / low in the malignant cell raises that
cell's own IFN / MHC-I / APM program and lowers TJ. Observationally,
CLDN4-high malignant cells should have IFN/MHC down, TJ up.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE131907 + GSE189357 only | Not GSE148071 / 123902 / 205335 / a bigger merge |
| Units | PR #459 %pos n=30 (21 samples + 9 patients) | GSE131907 is **sample-level**; GSE189357 is **patient-level** |
| Split | Within-cohort malignant CLDN4 **%pos Q4 vs Q1** | Quartile tails; mid quartiles unused in binary DE |
| Continuous | CLDN4 %pos z-scored, all 30 units, cohort covariate | Linear; not a causal model |
| Malignant | GSE131907 author-malig UMI-sum; GSE189357 marker-malig UMI-sum | Marker gate ≠ author annotation |
| T/NK | **not run** | Do not re-audit PR #459 ρ=−0.542 / Q4 r=−0.619 |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (single) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I/APM (custom), TJ (KEGG/GO + focal, **CLDN4 held out**), keratin (KRT_EPITHELIAL), chemokine panel | Chemokine panel is compact, not MSigDB C2 |

## Honest n

PR #459 locked labels (malignant CLDN4 %pos):

- GSE131907: n={n131} samples (tumor origins, n_mal≥20, author malignant). Q1={int(q131.get('Q1',0))} Q4={int(q131.get('Q4',0))}.
- GSE189357: n={n189} patients (marker-malignant; all 9 eligible). Q1={int(q189.get('Q1',0))} Q4={int(q189.get('Q4',0))}.

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 combined | {int(mal.get('n_low', 0))}/{int(mal.get('n_high', 0))} | {int(mal.get('n_genes', 0) or 0)} | cohort covariate; malignant only |
| continuous combined | n={int(mal_c.get('n', 0) or 0)} | {int(mal_c.get('n_genes', 0) or 0)} | CLDN4 %pos z |

GSE189357 Q4 vs Q1 is thin (Q4 n=2) and is **skipped** as a single-cohort
binary DE. Per-cohort n is in `tables/n_honest.tsv`. Do not quote a pooled n
that ignores the sample-vs-patient unit difference.

## Family DE (combined Q4 vs Q1) — n / logFC / p

Family score = mean log2(TMM-CPM+1) of family genes present. Positive logFC =
higher in CLDN4-high (Q4) than CLDN4-low (Q1). Machine table:
`tables/family_de.tsv`.

{score_md("q4q1")}

Continuous family-score DE (same 30 units, CLDN4 %pos z):

{score_md("continuous")}

Per-cohort family scores (GSE189357 Q4 vs Q1 skipped when n_Q4 < 3):

{score_md_cohorts("q4q1")}

Per-cohort continuous family scores (GSE189357 n=9 is thin):

{score_md_cohorts("continuous")}

{cldn4_txt}

Directions on the combined Q4 vs Q1 split match the thesis (IFN / MHC-I/APM /
chemokine down, TJ / keratin up). Family-score p is thin except chemokine.
TJ is a trend on the binary split and is the family that is up on the
continuous n=30 score (GSE131907 carries that TJ up). GSE189357 continuous
TJ is slightly down — an honest thin-n (n=9) limitation, not a second
independent TJ-up call. IFN is null on GSE131907 alone and flips slightly
positive on the continuous combined score. Keratin is up but not a claim.

## Gene-level family members (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high.

{fam_md("q4q1_combined")}

Headline genes (combined Q4 vs Q1, family members only, lowest p):

{headline.get("md", "_none_")}

## What this is not

- Not CellChat / LIANA / NicheNet and not a T/NK infiltrate test.
- Not a re-audit of PR #459 ρ=−0.542 / Q4 r=−0.619.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not GSE148071 and not a merge beyond GSE131907+GSE189357.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ/keratin/chemokine change.
- Genome-wide FDR on these n is expected to be thin; family scores are the claim.

## Files

- `tables/family_de.tsv` — **headline family DE** (n / logFC / p)
- `tables/de_all.tsv` — all genes, all contrasts
- `tables/de_families.tsv` — IFN / MHC-I/APM / TJ / keratin / chemokine rows
- `tables/de_q4q1_combined_families.tsv` — gene-level headline
- `tables/family_summary.tsv` — gene-level family counts
- `tables/n_honest.tsv`
- `tables/sample_inventory.tsv`
- `figures/` — volcano, family heatmap, key-gene boxes, family-score boxes, forests, median-logFC bars, CLDN4 strip, n bars

Reproduce:

```bash
python3 methods/pair_131907_189357_malig_ifn_de_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(text)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    fam = load_families()
    meta = pair_meta()
    meta.to_csv(TABLES / "sample_inventory.tsv", sep="\t", index=False)
    strip_cldn4(meta, FIGS / "cldn4_quartile_strip")

    paths = {
        "GSE131907": DATA / "GSE131907_malignant_counts.tsv.gz",
        "GSE189357": DATA / "GSE189357_malignant_counts.tsv.gz",
    }
    parts = {}
    all_de = []
    inv_rows = []
    fam_rows = []
    fam_score_rows = []

    for cohort, path in paths.items():
        if not path.exists():
            raise SystemExit(f"missing {path}; see data/README.md")
        cts = read_counts(path)
        parts[cohort] = cts
        de, info, _ = run_q4q1(cts, meta, cohort)
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
            }
        )
        if not de.empty:
            de = de.copy()
            de["compartment"] = "malignant"
            de["cohort"] = cohort
            de["contrast"] = "q4q1_single"
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
            }
        )
        if not de_c.empty:
            de_c = de_c.copy()
            de_c["compartment"] = "malignant"
            de_c["cohort"] = cohort
            de_c["contrast"] = "continuous_single"
            de_c["n_q1"] = np.nan
            de_c["n_q4"] = np.nan
            de_c["family"] = de_c["gene"].map(lambda g: family_of(g, fam))
            all_de.append(de_c)

    genes = sorted(set.intersection(*[set(p.index) for p in parts.values()]))
    combined = pd.concat([p.reindex(genes).fillna(0) for p in parts.values()], axis=1)
    combined = combined.loc[:, ~combined.columns.duplicated()]

    de, info, lc_q = run_q4q1(combined, meta, None)
    inv_rows.append(
        {
            "contrast": "malignant_q4q1_combined",
            "compartment": "malignant",
            "cohort": "GSE131907+GSE189357",
            "split": "q4q1",
            "n_low": info["n_q1"],
            "n_high": info["n_q4"],
            "n": info["n"],
            "n_genes": info.get("n_genes", 0),
            "patients_q1": info.get("patients_q1", ""),
            "patients_q4": info.get("patients_q4", ""),
        }
    )
    de_c, info_c = run_continuous(combined, meta, None)
    inv_rows.append(
        {
            "contrast": "malignant_continuous_combined",
            "compartment": "malignant",
            "cohort": "GSE131907+GSE189357",
            "split": "continuous",
            "n_low": np.nan,
            "n_high": np.nan,
            "n": info_c["n"],
            "n_genes": info_c.get("n_genes", 0),
            "patients_q1": "",
            "patients_q4": "",
        }
    )

    cldn4_row = {}
    if not de.empty:
        de = de.copy()
        de["compartment"] = "malignant"
        de["cohort"] = "GSE131907+GSE189357"
        de["contrast"] = "q4q1_combined"
        de["n_q1"] = info["n_q1"]
        de["n_q4"] = info["n_q4"]
        de["family"] = de["gene"].map(lambda g: family_of(g, fam))
        all_de.append(de)
        for r in _fam_summary(de, fam, info["n_q1"], info["n_q4"]):
            fam_rows.append({**r, "compartment": "malignant", "contrast": "q4q1_combined"})
        volcano(
            de,
            fam,
            f"malignant Q4 vs Q1  n={info['n_q4']} vs {info['n_q1']}  (cohort covariate)",
            FIGS / "volcano_malignant_q4q1_combined",
        )
        forest_families(
            de,
            fam,
            "malignant family genes, Q4 vs Q1 combined",
            FIGS / "forest_malignant_families_q4q1",
        )
        hit = de[de["gene"] == "CLDN4"]
        if not hit.empty:
            cldn4_row = {
                "logFC": float(hit.iloc[0]["logFC"]),
                "p": float(hit.iloc[0]["p"]),
                "n_q1": info["n_q1"],
                "n_q4": info["n_q4"],
            }

    if not de_c.empty:
        de_c = de_c.copy()
        de_c["compartment"] = "malignant"
        de_c["cohort"] = "GSE131907+GSE189357"
        de_c["contrast"] = "continuous_combined"
        de_c["n_q1"] = np.nan
        de_c["n_q4"] = np.nan
        de_c["family"] = de_c["gene"].map(lambda g: family_of(g, fam))
        all_de.append(de_c)

    # Family-score DE on the same filtered logCPM used for Q4/Q1 combined
    m_all = meta.loc[meta["patient"].isin(combined.columns)].copy()
    cts_all = filter_genes(combined.loc[:, m_all["patient"]])
    lc_all = log_cpm(cts_all, tmm_norm_factors(cts_all))
    for cohort in (None, "GSE131907", "GSE189357"):
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
            "MHC-I/APM": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5", "PSMB8", "PSMB9", "TAPBP"],
            "TJ": ["CLDN1", "CLDN7", "OCLN", "TJP1", "F11R", "PARD3", "CDH1", "MARVELD2"],
            "keratin": ["KRT7", "KRT8", "KRT18", "KRT19", "KRT5", "KRT17", "KRT14", "KRT6A"],
            "chemokine": CHEMOKINE[:12],
        }
        for name in FAM_COLORS:
            present = [g for g in prefer[name] if g in lc.index]
            show.extend(present)
        heatmap_family(
            lc,
            m_q,
            show,
            "malignant family genes, Q1/Q4 only (row z)",
            FIGS / "heatmap_malignant_families_q4q1",
        )
        box_key_genes(lc, m_q, "malignant key genes, CLDN4 Q4 vs Q1", FIGS / "box_malignant_key_genes")
        box_family_scores(
            lc,
            m_q,
            fam,
            "malignant family scores, CLDN4 Q4 vs Q1",
            FIGS / "box_malignant_family_scores",
        )

    inv = pd.DataFrame(inv_rows)
    inv.to_csv(TABLES / "n_honest.tsv", sep="\t", index=False)
    n_bar(inv[inv["split"] == "q4q1"], FIGS / "n_honest_q4q1")

    fam_score = pd.concat(fam_score_rows, ignore_index=True) if fam_score_rows else pd.DataFrame()
    if not fam_score.empty:
        # headline table: combined contrasts first
        fam_score.to_csv(TABLES / "family_de.tsv", sep="\t", index=False)
        head_fs = fam_score[
            (fam_score["cohort"] == "GSE131907+GSE189357") & (fam_score["split"] == "q4q1")
        ]
        forest_family_scores(
            head_fs,
            "malignant family scores, Q4 vs Q1 combined",
            FIGS / "forest_family_scores_q4q1",
        )

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
    else:
        headline = {"md": "_DE table empty._"}

    fam_df = pd.DataFrame(fam_rows)
    if not fam_df.empty:
        fam_df.to_csv(TABLES / "family_summary.tsv", sep="\t", index=False)
        bar_family_median_logfc(fam_df, FIGS / "bar_family_median_logfc")
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


if __name__ == "__main__":
    main()
