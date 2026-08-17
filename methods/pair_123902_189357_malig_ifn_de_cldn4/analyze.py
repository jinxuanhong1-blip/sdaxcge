#!/usr/bin/env python3
"""CLDN4-only tumor-cell-intrinsic patient-pseudobulk DE.

Pair: GSE123902 + GSE189357 (PR #459 given cut that differs).
Not CellChat. No dual-high. Do not re-audit the T/NK ρ.

Method: patient/donor UMI-sum of marker-malignant cells → TMM →
log2(CPM+1) → OLS. Families: IFN, MHC-I/APM, TJ (CLDN4 held out).
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

# PR #459 given (not re-audited)
GIVEN = {
    "combo": "GSE123902+GSE189357",
    "score": "pct",
    "n": 22,
    "rho": -0.638,
    "p": 0.003,
    "n_q1": 7,
    "n_q4": 5,
    "note": "PR #459 given; T/NK ρ not re-audited; Q4 vs Q1 tails thin",
}

FOCAL_BOX = [
    "CLDN4", "STAT1", "IRF1", "HLA-A", "HLA-B", "B2M", "TAP1", "TAP2",
    "NLRC5", "OCLN", "TJP1", "CLDN1", "CLDN7", "ISG15", "MX1", "PSMB8",
]
JUNCTION_FOCAL = [
    "CLDN1", "CLDN7", "F11R", "PARD3", "OCLN", "TJP1", "CDH1",
]


def load_families() -> dict[str, set[str]]:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    tj = (
        set(sets["KEGG_TIGHT_JUNCTION"])
        | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
        | set(sets.get("GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY", []))
        | set(JUNCTION_FOCAL)
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


def pair_meta() -> pd.DataFrame:
    a = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = a[a["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"] == True].copy()
    el["patient"] = el["patient"].astype(str)
    el["cohort"] = "GSE123902"
    el["unit"] = "donor"
    el["cldn4_pct"] = el["mal_CLDN4_pct"].astype(float)
    el["cldn4_mean"] = el["mal_CLDN4_mean"].astype(float)
    el["n_malignant"] = el["n_malignant"].astype(int)
    el["n_tnk"] = el["n_tnk"].astype(int)
    el["quartile"] = assign_quartiles(el.set_index("patient")["cldn4_pct"]).values

    b = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    el2 = b[b["eligible"] == True].copy()
    el2["patient"] = el2["patient"].astype(str)
    el2["cohort"] = "GSE189357"
    el2["unit"] = "patient"
    el2["cldn4_pct"] = el2["mal_CLDN4_pct"].astype(float)
    el2["cldn4_mean"] = el2["mal_CLDN4_mean"].astype(float)
    el2["n_malignant"] = el2["n_malignant"].astype(int)
    el2["n_tnk"] = el2["n_tnk"].astype(int)
    el2["quartile"] = assign_quartiles(el2.set_index("patient")["cldn4_pct"]).values

    cols = [
        "patient", "cohort", "unit", "cldn4_pct", "cldn4_mean",
        "n_malignant", "n_tnk", "quartile",
    ]
    return pd.concat([el[cols], el2[cols]], ignore_index=True)


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
    if cohort is None and m["cohort"].nunique() > 1:
        design["cohort_GSE189357"] = (m.set_index("patient")["cohort"] == "GSE189357").astype(float)
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
    z = m.set_index("patient")["cldn4_pct"]
    z = (z - z.mean()) / z.std(ddof=1)
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_pct_z"] = z
    if cohort is None and m["cohort"].nunique() > 1:
        design["cohort_GSE189357"] = (m.set_index("patient")["cohort"] == "GSE189357").astype(float)
    de = ols_de(lc, design, "CLDN4_pct_z")
    info["n_genes"] = int(len(de))
    info["patients"] = ",".join(sorted(m["patient"]))
    return de, info


def volcano(de: pd.DataFrame, fam: dict[str, set[str]], title: str, path: Path) -> None:
    if de.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    x = de["logFC"].values
    y = -np.log10(np.clip(de["p"].values, 1e-300, 1))
    ax.scatter(x, y, s=8, c="#c8c8c8", linewidths=0, alpha=0.7, label="other")
    colors = {"IFN": "#d62728", "MHC-I/APM": "#1f77b4", "TJ": "#2ca02c"}
    for fam_name, color in colors.items():
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
    fig, ax = plt.subplots(figsize=(max(8.0, 0.28 * len(m) + 2.6), max(4.8, 0.24 * len(genes) + 1.6)))
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
    fig, ax = plt.subplots(figsize=(7.6, max(4.0, 0.28 * len(plot) + 1.2)))
    y = np.arange(len(plot))
    colors = {"IFN": "#d62728", "MHC-I/APM": "#1f77b4", "TJ": "#2ca02c"}
    ax.axvline(0, c="#666", lw=0.7)
    ax.scatter(plot["logFC"], y, c=[colors[f] for f in plot["family"]], s=28, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.family}:{r.gene}" for r in plot.itertuples()], fontsize=7)
    ax.set_xlabel("log2FC (CLDN4-high − low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def family_median_bar(de: pd.DataFrame, fam: dict[str, set[str]], title: str, path: Path) -> None:
    if de.empty:
        return
    names, meds, cols = [], [], []
    colors = {"IFN": "#d62728", "MHC-I/APM": "#1f77b4", "TJ": "#2ca02c"}
    for name, genes in fam.items():
        sub = de[de["gene"].isin(genes)]
        if sub.empty:
            continue
        names.append(name)
        meds.append(float(sub["logFC"].median()))
        cols.append(colors[name])
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
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8))
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
    fig.suptitle("PR #459 pair CLDN4-only quartiles (within-cohort %pos) — tails thin")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def n_bar(inv: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    labels = inv["contrast"].tolist()
    x = np.arange(len(labels))
    ax.bar(x - 0.18, inv["n_low"].fillna(0), 0.36, label="CLDN4-low / Q1", color="#4c78a8")
    ax.bar(x + 0.18, inv["n_high"].fillna(0), 0.36, label="CLDN4-high / Q4", color="#e45756")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("honest n (patients/donors)")
    ax.set_title("DE sample sizes — GSE123902+GSE189357 (thin Q4)")
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


def write_finding(meta: pd.DataFrame, inv: pd.DataFrame, fam_rows: pd.DataFrame, headline: dict, cldn4_row: dict) -> None:
    def row(contrast: str) -> pd.Series:
        hit = inv[inv["contrast"] == contrast]
        return hit.iloc[0] if len(hit) else pd.Series(dtype=object)

    mal = row("malignant_q4q1_combined")
    mal_c = row("malignant_continuous_combined")
    s123 = row("malignant_q4q1_GSE123902")
    s189 = row("malignant_q4q1_GSE189357")

    def fam_md(contrast: str) -> str:
        sub = fam_rows[fam_rows["contrast"] == contrast]
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

    n123 = int((meta["cohort"] == "GSE123902").sum())
    n189 = int((meta["cohort"] == "GSE189357").sum())
    q123 = meta[meta["cohort"] == "GSE123902"]["quartile"].value_counts()
    q189 = meta[meta["cohort"] == "GSE189357"]["quartile"].value_counts()
    q1 = meta[meta["quartile"] == "Q1"]
    q4 = meta[meta["quartile"] == "Q4"]

    cldn4_txt = "CLDN4 not in the combined matrix."
    if cldn4_row:
        cldn4_txt = (
            f"CLDN4 itself {cldn4_row['logFC']:+.2f} "
            f"(p={_fmt_p(cldn4_row['p'])}, FDR={_fmt_p(cldn4_row['fdr'])}) "
            "— direction check on the split gene (held out of the TJ family)."
        )

    skip189 = s189.get("note", "")
    text = f"""# Pair GSE123902+GSE189357: CLDN4-only malignant IFN/MHC/TJ DE

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. **Not CellChat.**
Tumor-cell-intrinsic program only (marker-malignant / epithelial cells).

Given cut (PR #459, **not re-audited**): GSE123902 + GSE189357 malignant
CLDN4 %pos vs same-unit T/NK, n=22, Spearman ρ=−0.638 (p=0.003).
This extra does **not** re-audit that T/NK ρ.

**Thesis (already correct):** CLDN4 KD / low raises the malignant cell's own
IFN / MHC-I; CLDN4-high should be IFN/MHC down, TJ up.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
Patient/donor is the unit (GSE123902 donor-level; GSE189357 patient-level).
p-values are descriptive. Genome-wide FDR is thin on these n.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE123902 + GSE189357 only | The PR #459 pair that differs; not a bigger merge |
| Split | Within-cohort marker-malignant CLDN4 **%pos Q4 vs Q1** | Quartile tails, not median; mid quartiles unused in binary DE |
| Tails | Combined Q1 n=7 / Q4 n=5 | **Thin.** GSE189357 Q4 n=2 — single-cohort binary DE skipped |
| Continuous | CLDN4 %pos z, all 22 units, cohort covariate | Linear; not a causal model |
| Malignant | Marker-malignant UMI-sum: (EPCAM\\|KRT8\\|KRT18\\|KRT19)>0 and PTPRC==0 | Not author CNV-malignant; same gate as PR #459 |
| T/NK | not a DE compartment here | T/NK ρ taken as given from PR #459; not re-audited |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (GSE123902 only) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I/APM (custom), TJ (KEGG/GO + junction focal; **CLDN4 held out**) | MHC-II and chemokine panels are not the claim |

## Honest n

PR #459 pair labels (marker-malignant %pos; T/NK ρ not re-scored):

- GSE123902: n={n123} tumor/met donors (normals dropped). Q1={int(q123.get('Q1',0))} Q4={int(q123.get('Q4',0))}.
- GSE189357: n={n189} patients. Q1={int(q189.get('Q1',0))} Q4={int(q189.get('Q4',0))}.

Combined Q1: {", ".join(sorted(q1["patient"]))} (n={len(q1)}).
Combined Q4: {", ".join(sorted(q4["patient"]))} (n={len(q4)}).

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 combined | {int(mal.get('n_low', 0) or 0)}/{int(mal.get('n_high', 0) or 0)} | {int(mal.get('n_genes', 0) or 0)} | cohort covariate; **thin tails** |
| Q4 vs Q1 GSE123902 | {int(s123.get('n_low', 0) or 0)}/{int(s123.get('n_high', 0) or 0)} | {int(s123.get('n_genes', 0) or 0)} | donor-level; thin |
| Q4 vs Q1 GSE189357 | {int(s189.get('n_low', 0) or 0)}/{int(s189.get('n_high', 0) or 0)} | {int(s189.get('n_genes', 0) or 0)} | {skip189 or "skipped if Q4<3"} |
| continuous combined | n={int(mal_c.get('n', 0) or 0)} | {int(mal_c.get('n_genes', 0) or 0)} | CLDN4 %pos z; all gated units |

Q1 includes LX699 (46 malignant cells) and LX701 (90). Those units stay in
because they are in the given n=22 pair. Do not quote a cell-level n.

## IFN / MHC-I/APM / TJ (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high (Q4) than CLDN4-low (Q1).
Expected under the thesis: IFN down, MHC-I/APM down, TJ up.

{fam_md("q4q1_combined")}

{cldn4_txt}

Headline family genes (combined Q4 vs Q1, lowest p):

{headline.get("md", "_none_")}

Continuous (n=22) family summary is in `tables/family_summary.tsv` (contrast
`continuous_combined`). It is the better-powered direction check when the
Q4 tail is only 5 donors/patients.

## What this is not

- Not CellChat / LIANA / NicheNet.
- Not a re-audit of the PR #459 T/NK ρ.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not a merge beyond GSE123902+GSE189357.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ change.
- Q4 n=5 (and GSE189357 Q4 n=2) is a thin tail; genome-wide FDR on 7 vs 5
  is expected to be empty. Family median logFC and sign counts are the claim.

## Files

- `tables/de_q4q1_combined_families.tsv` — headline family DE table
- `tables/family_summary.tsv` — IFN / MHC-I/APM / TJ counts
- `tables/de_families.tsv` — family rows, all contrasts
- `tables/de_all.tsv` — all genes, all contrasts
- `tables/n_honest.tsv`
- `tables/sample_inventory.tsv`
- `figures/` — volcano, family heatmap, key-gene boxes, forest, family-median bar, CLDN4 strip, n bars

Reproduce:

```bash
python3 methods/pair_123902_189357_malig_ifn_de_cldn4/download.py
python3 methods/pair_123902_189357_malig_ifn_de_cldn4/build_malignant_pseudobulk.py
python3 methods/pair_123902_189357_malig_ifn_de_cldn4/analyze.py
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
        "GSE189357": DATA / "GSE189357_malignant_counts.tsv.gz",
    }
    parts = []
    all_de = []
    inv_rows = []
    fam_rows = []

    for cohort, path in paths.items():
        if not path.exists():
            raise SystemExit(f"missing {path}; run build_malignant_pseudobulk.py")
        cts = read_counts(path)
        parts.append(cts)
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

    genes = sorted(set.intersection(*[set(p.index) for p in parts]))
    combined = pd.concat([p.reindex(genes).fillna(0) for p in parts], axis=1)
    combined = combined.loc[:, ~combined.columns.duplicated()]

    de, info = run_q4q1(combined, meta, None)
    inv_rows.append(
        {
            "contrast": "malignant_q4q1_combined",
            "compartment": "malignant",
            "cohort": "GSE123902+GSE189357",
            "split": "q4q1",
            "n_low": info["n_q1"],
            "n_high": info["n_q4"],
            "n": info["n"],
            "n_genes": info.get("n_genes", 0),
            "patients_q1": info.get("patients_q1", ""),
            "patients_q4": info.get("patients_q4", ""),
            "note": "thin tails 7/5; cohort covariate",
        }
    )
    de_c, info_c = run_continuous(combined, meta, None)
    inv_rows.append(
        {
            "contrast": "malignant_continuous_combined",
            "compartment": "malignant",
            "cohort": "GSE123902+GSE189357",
            "split": "continuous",
            "n_low": np.nan,
            "n_high": np.nan,
            "n": info_c["n"],
            "n_genes": info_c.get("n_genes", 0),
            "patients_q1": "",
            "patients_q4": "",
            "note": "CLDN4 %pos z; n=22",
        }
    )

    cldn4_row = {}
    if not de.empty:
        de = de.copy()
        de["compartment"] = "malignant"
        de["cohort"] = "GSE123902+GSE189357"
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
            }
        volcano(
            de,
            fam,
            f"Malignant Q4 vs Q1  n={info['n_q4']} vs {info['n_q1']}  (thin tails; cohort covariate)",
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
        de_c["cohort"] = "GSE123902+GSE189357"
        de_c["contrast"] = "continuous_combined"
        de_c["n_q1"] = np.nan
        de_c["n_q4"] = np.nan
        de_c["family"] = de_c["gene"].map(lambda g: family_of(g, fam))
        all_de.append(de_c)
        for r in _fam_summary(de_c, fam, np.nan, np.nan):
            fam_rows.append({**r, "compartment": "malignant", "contrast": "continuous_combined"})

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
            "Malignant key genes, CLDN4 Q4 vs Q1 (n=7 vs 5; thin)",
            FIGS / "box_malignant_key_genes",
        )

    inv = pd.DataFrame(inv_rows)
    inv.to_csv(TABLES / "n_honest.tsv", sep="\t", index=False)
    n_bar(inv[inv["split"] == "q4q1"], FIGS / "n_honest_q4q1")

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
    write_finding(meta, inv, fam_df if not fam_df.empty else pd.DataFrame(), headline, cldn4_row)
    print("wrote FINDING.md and tables/")
    print(inv.to_string(index=False))
    if not fam_df.empty:
        print(fam_df[fam_df["contrast"] == "q4q1_combined"].to_string(index=False))


if __name__ == "__main__":
    main()
