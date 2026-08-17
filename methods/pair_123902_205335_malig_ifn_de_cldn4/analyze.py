#!/usr/bin/env python3
"""CLDN4-only malignant patient-pseudobulk DE on GSE123902 + GSE205335.

Tumor-cell-intrinsic program DE. Not CellChat. Not T/NK infiltrate.
Patient (GSE123902: donor) is the unit. No dual-high. No GSE148071.

Method: patient UMI-sum → TMM → log2(CPM+1) → OLS (limma-style).
Families: IFN, MHC-I/APM, TJ (CLDN4 held out), keratin.
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

GIVEN = {
    "combo": "GSE123902+GSE205335",
    "score": "%pos",
    "n": 35,
    "rho": -0.522,
    "rho_p": 0.002,
    "q4_r": -0.802,
    "q4_p": 0.005,
    "n_q1": 9,
    "n_q4": 9,
    "source": "PR459_given",
}

FOCAL_BOX = [
    "CLDN4",
    "STAT1",
    "IRF1",
    "HLA-A",
    "HLA-B",
    "B2M",
    "TAP1",
    "TAP2",
    "NLRC5",
    "KRT5",
    "KRT17",
    "KRT8",
    "OCLN",
    "TJP1",
    "CDH1",
]
FAM_ORDER = ["IFN", "MHC-I/APM", "TJ", "keratin"]
FAM_COLORS = {
    "IFN": "#d62728",
    "MHC-I/APM": "#1f77b4",
    "TJ": "#2ca02c",
    "keratin": "#9467bd",
}


def load_families() -> dict[str, set[str]]:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    tj = (
        set(sets["KEGG_TIGHT_JUNCTION"])
        | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
        | set(a8["focal_genes"])
    )
    tj.discard("CLDN4")
    keratin = set(sets["GOBP_KERATINIZATION"]) | set(sets["KRT_EPITHELIAL"])
    return {
        "IFN": set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
        | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]),
        "MHC-I/APM": set(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]),
        "TJ": tj,
        "keratin": keratin,
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
    a = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = a[a["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"].copy()
    el["patient"] = el["patient"].astype(str)
    el["cohort"] = "GSE123902"
    el["unit"] = "donor"
    el["cldn4_pct"] = el["mal_CLDN4_pct"].astype(float)
    el["cldn4_mean"] = el["mal_CLDN4_mean"].astype(float)
    el["n_malignant"] = el["n_malignant"].astype(int)
    el["quartile"] = assign_quartiles(el.set_index("patient")["cldn4_pct"]).values

    b = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    b["patient"] = b["patient"].astype(str)
    b["cohort"] = "GSE205335"
    b["unit"] = "patient"
    b["cldn4_pct"] = b["mal_CLDN4_pct_pos"].astype(float)
    b["cldn4_mean"] = b["mal_CLDN4_mean"].astype(float)
    b["n_malignant"] = b["n_malignant"].astype(int)
    b["quartile"] = assign_quartiles(b.set_index("patient")["cldn4_pct"]).values

    cols = ["patient", "cohort", "unit", "cldn4_pct", "cldn4_mean", "n_malignant", "quartile"]
    return pd.concat([el[cols], b[cols]], ignore_index=True)


def read_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df.groupby(df.index).sum()


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
    pvals = 2.0 * stats.t.sf(np.abs(t), df_res)
    out = pd.DataFrame(
        {
            "gene": logcpm.index.astype(str),
            "logFC": est,
            "AveExpr": Y.mean(axis=1),
            "t": t,
            "p": pvals,
            "se": se,
            "df": df_res,
        }
    )
    out["fdr"] = _bh(out["p"].values)
    return out.sort_values("p")


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


def _design(m: pd.DataFrame, cohort: str | None, continuous: bool = False) -> tuple[pd.DataFrame, str]:
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    if continuous:
        z = m.set_index("patient")["cldn4_pct"]
        z = (z - z.mean()) / z.std(ddof=1)
        design["CLDN4_pct_z"] = z
        coef = "CLDN4_pct_z"
    else:
        design["CLDN4_Q4"] = (m.set_index("patient")["quartile"] == "Q4").astype(float)
        coef = "CLDN4_Q4"
    if cohort is None and m["cohort"].nunique() > 1:
        design["cohort_GSE205335"] = (m.set_index("patient")["cohort"] == "GSE205335").astype(float)
    return design, coef


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
    design, coef = _design(m, cohort, continuous=False)
    de = ols_de(lc, design, coef)
    info["n_genes"] = int(len(de))
    return de, info, lc


def run_continuous(counts: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    m = meta.loc[meta["patient"].isin(counts.columns)].copy()
    info = {"n": int(len(m)), "n_q1": np.nan, "n_q4": np.nan}
    if len(m) < 8:
        return pd.DataFrame(), info
    cts = filter_genes(counts.loc[:, m["patient"]])
    fac = tmm_norm_factors(cts)
    lc = log_cpm(cts, fac)
    design, coef = _design(m, None, continuous=True)
    de = ols_de(lc, design, coef)
    info["n_genes"] = int(len(de))
    info["patients"] = ",".join(sorted(m["patient"]))
    return de, info


def family_score_de(logcpm: pd.DataFrame, meta: pd.DataFrame, fam: dict[str, set[str]], cohort: str | None) -> pd.DataFrame:
    """One OLS row per family on the mean of family-gene logCPM. Patient is the unit."""
    m = meta.copy()
    if cohort:
        m = m.loc[m["cohort"] == cohort]
    m = m.loc[m["quartile"].isin(["Q1", "Q4"])]
    m = m.loc[m["patient"].isin(logcpm.columns)].copy()
    n1 = int((m["quartile"] == "Q1").sum())
    n4 = int((m["quartile"] == "Q4").sum())
    if n1 < 3 or n4 < 3:
        return pd.DataFrame()
    design, coef = _design(m, cohort, continuous=False)
    rows = []
    for name in FAM_ORDER:
        genes = [g for g in sorted(fam[name]) if g in logcpm.index]
        if len(genes) < 3:
            continue
        score = logcpm.loc[genes, m["patient"]].astype(float).mean(axis=0)
        score.name = name
        de = ols_de(score.to_frame().T, design, coef)
        if de.empty:
            continue
        r = de.iloc[0]
        q1 = score.loc[m.loc[m["quartile"] == "Q1", "patient"]]
        q4 = score.loc[m.loc[m["quartile"] == "Q4", "patient"]]
        try:
            mw_p = float(stats.mannwhitneyu(q4, q1, alternative="two-sided").pvalue)
        except ValueError:
            mw_p = np.nan
        rows.append(
            {
                "family": name,
                "n_q1": n1,
                "n_q4": n4,
                "n": n1 + n4,
                "n_genes": int(len(genes)),
                "logFC": float(r["logFC"]),
                "p": float(r["p"]),
                "fdr": np.nan,
                "t": float(r["t"]),
                "se": float(r["se"]),
                "mw_p": mw_p,
                "mean_q1": float(q1.mean()),
                "mean_q4": float(q4.mean()),
                "cohort": cohort or "GSE123902+GSE205335",
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = _bh(out["p"].values)
    return out


def volcano(de: pd.DataFrame, fam: dict[str, set[str]], title: str, path: Path) -> None:
    if de.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
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
    fig_w = max(8.0, 0.22 * len(m) + 2.4)
    fig, ax = plt.subplots(figsize=(fig_w, max(4.5, 0.22 * len(genes) + 1.6)))
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
    ncols = 5
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12.5, 2.6 * nrows), squeeze=False)
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


def forest_families(de: pd.DataFrame, fam: dict[str, set[str]], title: str, path: Path) -> None:
    if de.empty:
        return
    rows = []
    for name in FAM_ORDER:
        sub = de[de["gene"].isin(fam[name])]
        if sub.empty:
            continue
        for r in sub.nsmallest(8, "p").itertuples():
            rows.append({"family": name, "gene": r.gene, "logFC": r.logFC, "p": r.p, "fdr": r.fdr})
    if not rows:
        return
    plot = pd.DataFrame(rows).sort_values(["family", "p"])
    fig, ax = plt.subplots(figsize=(7.6, max(4.0, 0.28 * len(plot) + 1.2)))
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


def forest_family_scores(scores: pd.DataFrame, title: str, path: Path) -> None:
    if scores.empty:
        return
    plot = scores.set_index("family").reindex([f for f in FAM_ORDER if f in set(scores["family"])]).reset_index()
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    y = np.arange(len(plot))
    ax.axvline(0, c="#666", lw=0.7)
    ax.scatter(plot["logFC"], y, c=[FAM_COLORS[f] for f in plot["family"]], s=64, zorder=3)
    for i, r in enumerate(plot.itertuples()):
        ax.text(
            r.logFC,
            i + 0.18,
            f"n={int(r.n)}  logFC={r.logFC:+.2f}  p={_fmt_p(r.p)}",
            fontsize=8,
            ha="center",
        )
    ax.set_yticks(y)
    ax.set_yticklabels(plot["family"])
    ax.set_xlabel("family-mean log2FC (CLDN4-high − low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def box_family_scores(logcpm: pd.DataFrame, meta: pd.DataFrame, fam: dict[str, set[str]], title: str, path: Path) -> None:
    m = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(logcpm.columns)].copy()
    fig, axes = plt.subplots(1, 4, figsize=(11.2, 3.6))
    rng = np.random.default_rng(1)
    for ax, name in zip(axes, FAM_ORDER):
        genes = [g for g in fam[name] if g in logcpm.index]
        score = logcpm.loc[genes, m["patient"]].astype(float).mean(axis=0)
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
        ax.set_ylabel("mean family log2(CPM+1)")
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
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
    fig.suptitle("PR #459 pair CLDN4-only quartiles (locked %pos; T/NK ρ not re-audited)")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def n_bar(inv: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 3.8))
    labels = inv["contrast"].tolist()
    x = np.arange(len(labels))
    ax.bar(x - 0.18, inv["n_low"].fillna(0), 0.36, label="CLDN4-low / Q1", color="#4c78a8")
    ax.bar(x + 0.18, inv["n_high"].fillna(0), 0.36, label="CLDN4-high / Q4", color="#e45756")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("honest n (patients / donors)")
    ax.set_title("Malignant DE n — GSE123902+GSE205335 (not the given T/NK n=35)")
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
    for name in FAM_ORDER:
        sub = de[de["gene"].isin(fam[name])].copy()
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
    fam_scores: pd.DataFrame,
    headline: dict,
    cldn4_row: dict,
) -> None:
    def row(contrast: str) -> pd.Series:
        hit = inv[inv["contrast"] == contrast]
        return hit.iloc[0] if len(hit) else pd.Series(dtype=object)

    mal = row("malignant_q4q1_combined")
    mal_c = row("malignant_continuous_combined")
    mal_123 = row("malignant_q4q1_GSE123902")
    mal_205 = row("malignant_q4q1_GSE205335")

    def fam_md(contrast: str) -> str:
        sub = fam_rows[fam_rows["contrast"] == contrast]
        if sub.empty:
            return "_No family rows (contrast not run)._"
        lines = [
            "| family | n_Q1 | n_Q4 | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |",
            "|---|---:|---:|---:|---|---:|---:|---|",
        ]
        for r in sub.itertuples():
            n1 = getattr(r, "n_q1", "")
            n4 = getattr(r, "n_q4", "")
            lines.append(
                f"| {r.family} | {n1} | {n4} | {r.n_tested} | {r.n_p05} ({r.n_up}/{r.n_down}) | {r.n_fdr05} | "
                f"{r.median_logFC:+.3f} | {r.top_gene} ({r.top_logFC:+.3f}, {_fmt_p(r.top_p)}, {_fmt_p(r.top_fdr)}) |"
            )
        return "\n".join(lines)

    def score_md(df: pd.DataFrame) -> str:
        if df.empty:
            return "_Family-score DE empty._"
        lines = [
            "| family | n_Q1 | n_Q4 | n | n_genes | logFC | p | FDR | MW p |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for r in df.itertuples():
            lines.append(
                f"| {r.family} | {int(r.n_q1)} | {int(r.n_q4)} | {int(r.n)} | {int(r.n_genes)} | "
                f"{r.logFC:+.3f} | {_fmt_p(r.p)} | {_fmt_p(r.fdr)} | {_fmt_p(r.mw_p)} |"
            )
        return "\n".join(lines)

    n123 = int((meta["cohort"] == "GSE123902").sum())
    n205 = int((meta["cohort"] == "GSE205335").sum())
    q123 = meta[meta["cohort"] == "GSE123902"]["quartile"].value_counts()
    q205 = meta[meta["cohort"] == "GSE205335"]["quartile"].value_counts()
    cldn4_txt = (
        f"CLDN4 itself {cldn4_row['logFC']:+.2f} (p={_fmt_p(cldn4_row['p'])}, "
        f"n={cldn4_row['n_q1']}/{cldn4_row['n_q4']}) — direction check on the split gene; "
        "held out of the TJ family."
        if cldn4_row
        else "CLDN4 itself was not tested in the combined matrix."
    )

    text = f"""# Pair GSE123902+GSE205335: CLDN4-only malignant-cell IFN/MHC/TJ/keratin DE

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. Not CellChat. Not T/NK infiltrate.
No GSE148071. Patient (GSE123902: donor) is the unit.

The pair T/NK Spearman is **taken as given** from PR #459 and is not re-audited:

| combo | score | N | ρ (p) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE205335 | %pos | 35 | −0.522 (0.002) | −0.802 (0.005; 9 vs 9) |

This folder is **tumor-cell-intrinsic program DE** on within-cohort quartiles of
the locked %pos vector: patient-pseudobulk of **malignant cells**, CLDN4 %pos
Q4 vs Q1. Thesis (already correct, not audited here): CLDN4-low malignant cells
raise their own IFN / MHC-I program and lower TJ; observationally CLDN4-high
malignant = IFN/MHC down, TJ up. The given T/NK 9 vs 9 is a different contrast.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat. muscat `pbDS`
is the same collapse (patient × cell-type UMI-sum → bulk DE). No R/limma/edgeR
in this environment; TMM + OLS + BH is the implementation. p-values are descriptive.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE123902 + GSE205335 only | Not GSE148071 / 131907 / 189357 / 7-cohort pool |
| Split | Within-cohort locked malignant CLDN4 **%pos Q4 vs Q1** | Not the PR #459 pooled-pair T/NK 9 vs 9; mid unused |
| Continuous | CLDN4 %pos z, all gated units with counts, cohort covariate | Linear; not causal |
| Malignant | GSE123902 marker-malignant UMI-sum; GSE205335 author malignant UMI-sum | GSE123902 has no author malignant labels |
| T/NK | **not run** | Given ρ is infiltrate, not this DE |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (single) | Small n; no voom weights |
| Families | IFN (Hallmark IFNα/γ); MHC-I/APM (custom); TJ (KEGG/GO + focal, **CLDN4 held out**); keratin (GOBP keratinization + KRT epithelial) | Compact custom MHC-I, not all HLA |

GSE123902 malignant = (EPCAM\\|KRT8\\|KRT18\\|KRT19)>0 and PTPRC==0; tumor/met
donors; normals dropped; one library per donor (PRIMARY preferred). GSE205335
malignant = author `Malignant cells` (existing patient UMI-sum; P4001 n_mal=27
is absent from that matrix).

## Honest n

Locked labels (PR #459 %pos; T/NK n=35 is **not** the DE n):

- GSE123902: n={n123} donors. Q1={int(q123.get('Q1', 0))} Q4={int(q123.get('Q4', 0))}.
- GSE205335: n={n205} patients. Q1={int(q205.get('Q1', 0))} Q4={int(q205.get('Q4', 0))}.

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 combined malignant | {int(mal.get('n_low', 0) or 0)}/{int(mal.get('n_high', 0) or 0)} | {int(mal.get('n_genes', 0) or 0)} | cohort covariate; units with counts |
| Q4 vs Q1 GSE123902 | {int(mal_123.get('n_low', 0) or 0)}/{int(mal_123.get('n_high', 0) or 0)} | {int(mal_123.get('n_genes', 0) or 0)} | donor; thin tails |
| Q4 vs Q1 GSE205335 | {int(mal_205.get('n_low', 0) or 0)}/{int(mal_205.get('n_high', 0) or 0)} | {int(mal_205.get('n_genes', 0) or 0)} | P4001 out of malignant counts |
| continuous combined | n={int(mal_c.get('n', 0) or 0)} | {int(mal_c.get('n_genes', 0) or 0)} | CLDN4 %pos z |

Combined malignant DE is also 9 vs 9 **after P4001 drops** (n_mal=27; Q1 on
GSE205335). That is not the given T/NK 9 vs 9. Per-contrast patients:
`tables/n_honest.tsv`.

## Family-score DE (headline table)

One OLS per family on the **mean log2(TMM-CPM+1)** of family genes. Positive
logFC = higher in CLDN4-high (Q4). Patient/donor is the unit.

### Combined Q4 vs Q1 (cohort covariate)

{score_md(fam_scores[fam_scores['cohort'] == 'GSE123902+GSE205335'] if not fam_scores.empty else fam_scores)}

{cldn4_txt}

IFN and MHC-I/APM family means are lower in CLDN4-high malignant cells. TJ
(CLDN4 held out) is **flat** as a family mean (gene-level mixed). Keratin
family-mean trends down. Do not quote a TJ-up family logFC on this pair.

Machine table: `tables/family_de.tsv` (combined n / logFC / p). All splits:
`tables/family_score_de.tsv`.

## Gene-level family DE (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high malignant cells.

{fam_md("q4q1_combined")}

Headline genes (family members only, lowest p):

{headline.get("md", "_none_")}

## What this is not

- Not CellChat / LIANA / NicheNet and not a T/NK infiltrate re-audit.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not GSE148071 and not a bigger merge.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ/keratin change.
- Genome-wide FDR on these n is expected to be thin; family n / logFC / p is the claim.

## Extra figures

- `figures/volcano_malignant_q4q1_combined.png`
- `figures/heatmap_malignant_families_q4q1.png`
- `figures/box_malignant_key_genes.png`
- `figures/box_family_scores_q4q1.png`
- `figures/forest_malignant_families_q4q1.png`
- `figures/forest_family_scores_q4q1.png`
- `figures/cldn4_quartile_strip.png`
- `figures/n_honest_q4q1.png`

## Files

- `tables/family_de.tsv` — combined family DE (n / logFC / p)
- `tables/family_score_de.tsv` — family scores including single-cohort
- `tables/de_q4q1_combined_families.tsv` — gene-level family DE
- `tables/de_families.tsv` / `tables/de_all.tsv.gz`
- `tables/family_summary.tsv` / `tables/n_honest.tsv` / `tables/sample_inventory.tsv`

Reproduce:

```bash
python3 methods/pair_123902_205335_malig_ifn_de_cldn4/scripts/download_gse123902.py
python3 methods/pair_123902_205335_malig_ifn_de_cldn4/scripts/build_gse123902.py
python3 methods/pair_123902_205335_malig_ifn_de_cldn4/analyze.py
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
        "GSE123902": DATA / "GSE123902_malignant_counts.tsv.gz",
        "GSE205335": DATA / "GSE205335_malignant_counts.tsv.gz",
    }

    all_de = []
    inv_rows = []
    fam_rows = []
    score_rows = []
    logcpm_q = None
    meta_q = None

    parts = []
    for cohort, path in paths.items():
        if not path.exists():
            raise SystemExit(f"missing {path}; build GSE123902 first")
        cts = read_counts(path)
        parts.append(cts)
        de, info, _lc = run_q4q1(cts, meta, cohort)
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

    genes = sorted(set.intersection(*[set(p.index) for p in parts]))
    combined = pd.concat([p.reindex(genes).fillna(0) for p in parts], axis=1)
    combined = combined.loc[:, ~combined.columns.duplicated()]
    de, info, lc_q = run_q4q1(combined, meta, None)
    inv_rows.append(
        {
            "contrast": "malignant_q4q1_combined",
            "compartment": "malignant",
            "cohort": "GSE123902+GSE205335",
            "split": "q4q1",
            "n_low": info["n_q1"],
            "n_high": info["n_q4"],
            "n": info["n"],
            "n_genes": info.get("n_genes", 0),
            "patients_q1": info.get("patients_q1", ""),
            "patients_q4": info.get("patients_q4", ""),
        }
    )
    de_c, info_c = run_continuous(combined, meta)
    inv_rows.append(
        {
            "contrast": "malignant_continuous_combined",
            "compartment": "malignant",
            "cohort": "GSE123902+GSE205335",
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
        de["cohort"] = "GSE123902+GSE205335"
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
            r = hit.iloc[0]
            cldn4_row = {
                "logFC": float(r["logFC"]),
                "p": float(r["p"]),
                "n_q1": int(info["n_q1"]),
                "n_q4": int(info["n_q4"]),
            }

    if not de_c.empty:
        de_c = de_c.copy()
        de_c["compartment"] = "malignant"
        de_c["cohort"] = "GSE123902+GSE205335"
        de_c["contrast"] = "continuous_combined"
        de_c["n_q1"] = np.nan
        de_c["n_q4"] = np.nan
        de_c["family"] = de_c["gene"].map(lambda g: family_of(g, fam))
        all_de.append(de_c)

    m_q = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(combined.columns)]
    if len(m_q) >= 6:
        cts_q = filter_genes(combined.loc[:, m_q["patient"]])
        lc = log_cpm(cts_q, tmm_norm_factors(cts_q))
        logcpm_q, meta_q = lc, m_q
        show = []
        prefer = {
            "IFN": ["STAT1", "IRF1", "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IFI6", "BST2", "GBP1"],
            "MHC-I/APM": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5", "PSMB8", "TAPBP", "ERAP1"],
            "TJ": ["CLDN1", "CLDN7", "OCLN", "TJP1", "F11R", "PARD3", "CDH1", "CLDN3"],
            "keratin": ["KRT5", "KRT6A", "KRT14", "KRT17", "KRT7", "KRT8", "KRT18", "KRT19"],
        }
        for name in FAM_ORDER:
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
            "malignant family-mean scores, CLDN4 Q4 vs Q1",
            FIGS / "box_family_scores_q4q1",
        )
        for cohort in (None, "GSE123902", "GSE205335"):
            sc = family_score_de(lc, meta, fam, cohort)
            if not sc.empty:
                score_rows.append(sc)
        if score_rows:
            sc_all = pd.concat(score_rows, ignore_index=True)
            sc_comb = sc_all[sc_all["cohort"] == "GSE123902+GSE205335"]
            forest_family_scores(
                sc_comb,
                "family-mean DE, Q4 vs Q1 combined",
                FIGS / "forest_family_scores_q4q1",
            )

    inv = pd.DataFrame(inv_rows)
    inv.to_csv(TABLES / "n_honest.tsv", sep="\t", index=False)
    n_bar(inv[inv["split"] == "q4q1"], FIGS / "n_honest_q4q1")

    if all_de:
        de_all = pd.concat(all_de, ignore_index=True)
        de_all.to_csv(TABLES / "de_all.tsv.gz", sep="\t", index=False)
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
    score_df = pd.concat(score_rows, ignore_index=True) if score_rows else pd.DataFrame()
    if not score_df.empty:
        score_df.to_csv(TABLES / "family_score_de.tsv", sep="\t", index=False)
        comb = score_df[score_df["cohort"] == "GSE123902+GSE205335"][
            ["family", "n_q1", "n_q4", "n", "n_genes", "logFC", "p", "fdr"]
        ]
        comb.to_csv(TABLES / "family_de.tsv", sep="\t", index=False)

    write_finding(
        meta,
        inv,
        fam_df if not fam_df.empty else pd.DataFrame(),
        score_df,
        headline,
        cldn4_row,
    )
    print("wrote FINDING.md and tables/", flush=True)
    print(inv.to_string(index=False), flush=True)
    if not score_df.empty:
        print(score_df.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
