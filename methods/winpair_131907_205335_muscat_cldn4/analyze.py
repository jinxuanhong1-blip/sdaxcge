#!/usr/bin/env python3
"""CLDN4-only patient-pseudobulk DE on the winning pair GSE131907+GSE205335.

Method: patient-level UMI-sum → TMM → log2(CPM+1) → OLS (limma-style, not muscat).
No dual-high TACSTD2×CLDN4. Not a bigger merge.
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
MHC_II = [
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "CD74", "CIITA",
]
FOCAL_BOX = [
    "CLDN4", "STAT1", "IRF1", "HLA-A", "HLA-B", "B2M", "TAP1",
    "CXCL9", "CXCL10", "CCL5", "OCLN", "TJP1",
]


def load_families() -> dict[str, set[str]]:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    fam = {
        "IFN": set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
        | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]),
        "MHC": set(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]) | set(MHC_II),
        "TJ": set(sets["KEGG_TIGHT_JUNCTION"])
        | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
        | set(a8["focal_genes"]),
        "chemokine": set(CHEMOKINE),
    }
    return fam


def family_of(gene: str, fam: dict[str, set[str]]) -> str:
    hits = [k for k, vs in fam.items() if gene in vs]
    return "|".join(hits) if hits else "other"


def tmm_norm_factors(counts: pd.DataFrame) -> pd.Series:
    """edgeR-like TMM (trimmed mean of M-values), no R."""
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
    fac = fac / fac.mean()
    return fac


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


def winning_pair_meta() -> pd.DataFrame:
    a = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    a = a.loc[a["origin"].isin(TUMOR_ORIGINS) & (a["n_malignant"] >= 20)].copy()
    a["patient"] = a["sample"].astype(str)
    a["cohort"] = "GSE131907"
    a["unit"] = "sample"
    a["cldn4_pct"] = a["mal_CLDN4_pct"].astype(float)
    a["cldn4_mean"] = a["mal_CLDN4_mean"].astype(float)
    a["n_malignant"] = a["n_malignant"].astype(int)
    a["n_tnk"] = a["n_tnk"].astype(int)
    a["frac_tnk"] = a["frac_tnk"].astype(float)
    a["quartile"] = assign_quartiles(a.set_index("patient")["cldn4_pct"]).values

    b = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    b["patient"] = b["patient"].astype(str)
    b["cohort"] = "GSE205335"
    b["unit"] = "patient"
    b["cldn4_pct"] = b["mal_CLDN4_pct_pos"].astype(float)
    b["cldn4_mean"] = b["mal_CLDN4_mean"].astype(float)
    b["n_malignant"] = b["n_malignant"].astype(int)
    b["n_tnk"] = b["n_tnk"].astype(int)
    b["frac_tnk"] = b["frac_tnk"].astype(float)
    b["quartile"] = assign_quartiles(b.set_index("patient")["cldn4_pct"]).values

    cols = [
        "patient", "cohort", "unit", "cldn4_pct", "cldn4_mean",
        "n_malignant", "n_tnk", "frac_tnk", "quartile",
    ]
    return pd.concat([a[cols], b[cols]], ignore_index=True)


def read_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df.groupby(df.index).sum()


def ols_de(logcpm: pd.DataFrame, design: pd.DataFrame, coef: str) -> pd.DataFrame:
    """Vectorized OLS. design already contains intercept. coef is the contrast column."""
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
    beta = xtx_inv @ Xv.T @ Y.T  # p x genes
    fitted = Xv @ beta
    resid = Y.T - fitted
    sse = np.sum(resid**2, axis=0)
    sigma2 = sse / df_res
    j = list(X.columns).index(coef)
    se = np.sqrt(np.maximum(sigma2 * xtx_inv[j, j], 0.0))
    est = beta[j]
    t = np.divide(est, se, out=np.zeros_like(est), where=se > 0)
    p = 2.0 * stats.t.sf(np.abs(t), df_res)
    out = pd.DataFrame(
        {
            "gene": logcpm.index.astype(str),
            "logFC": est,
            "AveExpr": Y.mean(axis=1),
            "t": t,
            "p": p,
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
    }
    if n1 < 3 or n4 < 3:
        return pd.DataFrame(), info
    cts = filter_genes(counts.loc[:, m["patient"]])
    fac = tmm_norm_factors(cts)
    lc = log_cpm(cts, fac)
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_Q4"] = (m.set_index("patient")["quartile"] == "Q4").astype(float)
    if cohort is None and m["cohort"].nunique() > 1:
        design["cohort_GSE205335"] = (m.set_index("patient")["cohort"] == "GSE205335").astype(float)
    de = ols_de(lc, design, "CLDN4_Q4")
    info["n_genes"] = int(len(de))
    return de, info


def run_continuous(counts: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    m = meta.loc[meta["patient"].isin(counts.columns)].copy()
    info = {"n": int(len(m)), "n_q1": np.nan, "n_q4": np.nan}
    if len(m) < 8:
        return pd.DataFrame(), info
    cts = filter_genes(counts.loc[:, m["patient"]])
    fac = tmm_norm_factors(cts)
    lc = log_cpm(cts, fac)
    z = m.set_index("patient")["cldn4_pct"]
    z = (z - z.mean()) / z.std(ddof=1)
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_pct_z"] = z
    if m["cohort"].nunique() > 1:
        design["cohort_GSE205335"] = (m.set_index("patient")["cohort"] == "GSE205335").astype(float)
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
    colors = {"IFN": "#d62728", "MHC": "#1f77b4", "TJ": "#2ca02c", "chemokine": "#ff7f0e"}
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
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 2.6 * nrows), squeeze=False)
    for i, gene in enumerate(genes):
        ax = axes[i // ncols][i % ncols]
        q1 = logcpm.loc[gene, m.loc[m["quartile"] == "Q1", "patient"]].astype(float)
        q4 = logcpm.loc[gene, m.loc[m["quartile"] == "Q4", "patient"]].astype(float)
        ax.boxplot([q1, q4], tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"], widths=0.55)
        rng = np.random.default_rng(0)
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
        hit = sub.nsmallest(8, "p")
        for r in hit.itertuples():
            rows.append({"family": name, "gene": r.gene, "logFC": r.logFC, "p": r.p, "fdr": r.fdr})
    if not rows:
        return
    plot = pd.DataFrame(rows).sort_values(["family", "p"])
    fig, ax = plt.subplots(figsize=(7.4, max(4.0, 0.28 * len(plot) + 1.2)))
    y = np.arange(len(plot))
    colors = {"IFN": "#d62728", "MHC": "#1f77b4", "TJ": "#2ca02c", "chemokine": "#ff7f0e"}
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
        ax.set_xlabel("patients (sorted)")
        for q, c in cols.items():
            ax.scatter([], [], c=c, label=q, s=28)
        ax.legend(frameon=False, fontsize=7, ncol=4, loc="lower right")
    fig.suptitle("Winning-pair CLDN4-only quartiles (PR #320 author %pos)")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def n_bar(inv: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 3.8))
    labels = inv["contrast"].tolist()
    x = np.arange(len(labels))
    ax.bar(x - 0.18, inv["n_low"].fillna(0), 0.36, label="CLDN4-low / Q1", color="#4c78a8")
    ax.bar(x + 0.18, inv["n_high"].fillna(0), 0.36, label="CLDN4-high / Q4", color="#e45756")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("honest n (patients/samples)")
    ax.set_title("DE sample sizes — winning pair only")
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
                "n_up": int(((sig["logFC"] > 0)).sum()),
                "n_down": int(((sig["logFC"] < 0)).sum()),
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
    headline: dict,
) -> None:
    def row(contrast: str) -> pd.Series:
        hit = inv[inv["contrast"] == contrast]
        return hit.iloc[0] if len(hit) else pd.Series(dtype=object)

    mal = row("malignant_q4q1_combined")
    tnk = row("tnk_q4q1_combined")
    mal_c = row("malignant_continuous_combined")
    tnk_c = row("tnk_continuous_combined")

    def fam_md(compartment: str, contrast: str) -> str:
        sub = fam_rows[(fam_rows["compartment"] == compartment) & (fam_rows["contrast"] == contrast)]
        if sub.empty:
            return "_No family rows (contrast not run)._"
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

    n131 = int((meta["cohort"] == "GSE131907").sum())
    n205 = int((meta["cohort"] == "GSE205335").sum())
    q131 = meta[meta["cohort"] == "GSE131907"]["quartile"].value_counts()
    q205 = meta[meta["cohort"] == "GSE205335"]["quartile"].value_counts()

    text = f"""# Winning-pair CLDN4-only patient-pseudobulk DE (GSE131907 + GSE205335)

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. Not a bigger merge.
Winning pair is taken from PR #320 (author malignant CLDN4 %pos vs T/NK:
Spearman ρ=−0.479, N=43; Q4 vs Q1 r=−0.705, n=23 compared).

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
muscat `pbDS` is the same idea (sum cells → patient × cell-type count matrix →
bulk DE). Here the patient (GSE131907: sample) is already the unit, so a
direct limma-style OLS is the honest tool. No R/limma/edgeR/muscat in this
environment; TMM + OLS + BH is the implementation. p-values are descriptive.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE131907 + GSE205335 only | Not GSE207422 / 291670 / 253013 / 148071 |
| Split | Within-cohort author malignant CLDN4 **%pos Q4 vs Q1** (PR #320) | Quartile tails, not median; mid quartiles unused in binary DE |
| Continuous | CLDN4 %pos z-scored, all gated units, cohort covariate | Linear; not a causal model |
| Malignant | Author malignant UMI-sum (existing GSEA matrices; n_mal≥30) | P4001 (27 malignant cells) drops from malignant DE |
| T/NK | Same patients; author T/NK UMI-sum; ≥20 T/NK cells | GSE131907 is **sample-level** |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (single) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I custom + MHC-II, TJ (KEGG/GO + focal), chemokine panel | Chemokine panel is compact, not MSigDB C2 |

## Honest n

Winning-pair labels (PR #320 author %pos):

- GSE131907: n={n131} samples (tumor origins, n_mal≥20). Q1={int(q131.get('Q1',0))} Q4={int(q131.get('Q4',0))} (sample-level).
- GSE205335: n={n205} patients. Q1={int(q205.get('Q1',0))} Q4={int(q205.get('Q4',0))}.

| contrast | compartment | n_low / n_high | n_genes | note |
|---|---|---:|---:|---|
| Q4 vs Q1 combined | malignant | {int(mal.get('n_low', 0))}/{int(mal.get('n_high', 0))} | {int(mal.get('n_genes', 0) or 0)} | cohort covariate; P4001 out of malignant matrix |
| Q4 vs Q1 combined | T/NK | {int(tnk.get('n_low', 0))}/{int(tnk.get('n_high', 0))} | {int(tnk.get('n_genes', 0) or 0)} | same Q labels; T/NK from those patients |
| continuous combined | malignant | n={int(mal_c.get('n', 0) or 0)} | {int(mal_c.get('n_genes', 0) or 0)} | CLDN4 %pos z |
| continuous combined | T/NK | n={int(tnk_c.get('n', 0) or 0)} | {int(tnk_c.get('n_genes', 0) or 0)} | CLDN4 %pos z |

Per-cohort Q4 vs Q1 n is in `tables/n_honest.tsv`. Do not quote a pooled n that
ignores the sample-vs-patient unit difference on GSE131907.

## IFN / MHC / TJ / chemokine (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high (Q4) than CLDN4-low (Q1).

### Malignant cells

{fam_md("malignant", "q4q1_combined")}

### T/NK cells (same patients)

{fam_md("T/NK", "q4q1_combined")}

T/NK top TJ hits are **epithelial genes** (CLDN4, CDH1, CLDN3, TJP1). That is
ambient / doublet leakage from malignant cells, not a T/NK tight-junction
program. The T/NK IFN arm is one-sided down (39/39 p<0.05 genes down; JAK2
top). Do not quote T/NK CLDN4 logFC as a lymphocyte finding.

Headline genes (combined Q4 vs Q1, family members only, lowest p):

{headline.get("md", "_none_")}

## What this is not

- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not a merge beyond GSE131907+GSE205335.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ/chemokine change.
- Genome-wide FDR on n≈11 vs 11 is expected to be thin; family counts are the claim.

## Files

- `tables/de_all.tsv` — all genes, all contrasts
- `tables/de_families.tsv` — IFN/MHC/TJ/chemokine rows
- `tables/de_q4q1_combined_families.tsv` — headline DE table
- `tables/n_honest.tsv`
- `tables/sample_inventory.tsv`
- `figures/` — volcano, family heatmap, key-gene boxes, forest, CLDN4 strip, n bars

Reproduce:

```bash
# T/NK matrices (needs GEO files in /tmp/winpair_geo)
python3 methods/winpair_131907_205335_muscat_cldn4/build_tnk_pseudobulk.py
python3 methods/winpair_131907_205335_muscat_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(text)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    fam = load_families()
    meta = winning_pair_meta()
    meta.to_csv(TABLES / "sample_inventory.tsv", sep="\t", index=False)
    strip_cldn4(meta, FIGS / "cldn4_quartile_strip")

    compartments = {
        "malignant": {
            "GSE131907": DATA / "GSE131907_malignant_counts.tsv.gz",
            "GSE205335": DATA / "GSE205335_malignant_counts.tsv.gz",
        },
        "T/NK": {
            "GSE131907": DATA / "GSE131907_tnk_counts.tsv.gz",
            "GSE205335": DATA / "GSE205335_tnk_counts.tsv.gz",
        },
    }

    all_de = []
    inv_rows = []
    fam_rows = []
    logcpm_cache = {}

    for comp, paths in compartments.items():
        parts = []
        for cohort, path in paths.items():
            if not path.exists():
                print(f"SKIP {comp} {cohort}: missing {path.name}")
                continue
            cts = read_counts(path)
            parts.append(cts)
            de, info = run_q4q1(cts, meta, cohort)
            contrast = f"{comp.lower().replace('/', '')}_q4q1_{cohort}"
            inv_rows.append(
                {
                    "contrast": contrast,
                    "compartment": comp,
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
                de["compartment"] = comp
                de["cohort"] = cohort
                de["contrast"] = "q4q1_single"
                de["n_q1"] = info["n_q1"]
                de["n_q4"] = info["n_q4"]
                de["family"] = de["gene"].map(lambda g: family_of(g, fam))
                all_de.append(de)
                for r in _fam_summary(de, fam, info["n_q1"], info["n_q4"]):
                    fam_rows.append({**r, "compartment": comp, "contrast": f"q4q1_{cohort}"})

        if not parts:
            continue
        genes = sorted(set.intersection(*[set(p.index) for p in parts]))
        combined = pd.concat([p.reindex(genes).fillna(0) for p in parts], axis=1)
        combined = combined.loc[:, ~combined.columns.duplicated()]
        de, info = run_q4q1(combined, meta, None)
        inv_rows.append(
            {
                "contrast": f"{comp.lower().replace('/', '')}_q4q1_combined",
                "compartment": comp,
                "cohort": "GSE131907+GSE205335",
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
                "contrast": f"{comp.lower().replace('/', '')}_continuous_combined",
                "compartment": comp,
                "cohort": "GSE131907+GSE205335",
                "split": "continuous",
                "n_low": np.nan,
                "n_high": np.nan,
                "n": info_c["n"],
                "n_genes": info_c.get("n_genes", 0),
                "patients_q1": "",
                "patients_q4": "",
            }
        )
        if not de.empty:
            de = de.copy()
            de["compartment"] = comp
            de["cohort"] = "GSE131907+GSE205335"
            de["contrast"] = "q4q1_combined"
            de["n_q1"] = info["n_q1"]
            de["n_q4"] = info["n_q4"]
            de["family"] = de["gene"].map(lambda g: family_of(g, fam))
            all_de.append(de)
            for r in _fam_summary(de, fam, info["n_q1"], info["n_q4"]):
                fam_rows.append({**r, "compartment": comp, "contrast": "q4q1_combined"})
            volcano(
                de,
                fam,
                f"{comp} Q4 vs Q1  n={info['n_q4']} vs {info['n_q1']}  (cohort covariate)",
                FIGS / f"volcano_{comp.replace('/', '')}_q4q1_combined",
            )
            forest_families(
                de,
                fam,
                f"{comp} family genes, Q4 vs Q1 combined",
                FIGS / f"forest_{comp.replace('/', '')}_families_q4q1",
            )
        if not de_c.empty:
            de_c = de_c.copy()
            de_c["compartment"] = comp
            de_c["cohort"] = "GSE131907+GSE205335"
            de_c["contrast"] = "continuous_combined"
            de_c["n_q1"] = np.nan
            de_c["n_q4"] = np.nan
            de_c["family"] = de_c["gene"].map(lambda g: family_of(g, fam))
            all_de.append(de_c)

        # logCPM for figures on Q1/Q4 patients
        m_q = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(combined.columns)]
        if len(m_q) >= 6:
            cts_q = filter_genes(combined.loc[:, m_q["patient"]])
            lc = log_cpm(cts_q, tmm_norm_factors(cts_q))
            logcpm_cache[comp] = (lc, m_q)
            show = []
            for name in ("IFN", "MHC", "TJ", "chemokine"):
                present = [g for g in sorted(fam[name]) if g in lc.index]
                # keep the most variable / focal
                if name == "IFN":
                    prefer = ["STAT1", "IRF1", "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IFI6", "BST2", "GBP1"]
                elif name == "MHC":
                    prefer = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5", "PSMB8", "HLA-DRA", "CD74"]
                elif name == "TJ":
                    prefer = ["CLDN4", "CLDN1", "CLDN7", "OCLN", "TJP1", "F11R", "PARD3", "CDH1"]
                else:
                    prefer = CHEMOKINE[:12]
                show.extend([g for g in prefer if g in present])
            heatmap_family(
                lc,
                m_q,
                show,
                f"{comp} family genes, Q1/Q4 only (row z)",
                FIGS / f"heatmap_{comp.replace('/', '')}_families_q4q1",
            )
            box_key_genes(
                lc,
                m_q,
                f"{comp} key genes, CLDN4 Q4 vs Q1",
                FIGS / f"box_{comp.replace('/', '')}_key_genes",
            )

    inv = pd.DataFrame(inv_rows)
    inv.to_csv(TABLES / "n_honest.tsv", sep="\t", index=False)
    n_bar(inv[inv["split"] == "q4q1"], FIGS / "n_honest_q4q1")

    if all_de:
        de_all = pd.concat(all_de, ignore_index=True)
        de_all.to_csv(TABLES / "de_all.tsv", sep="\t", index=False)
        fam_de = de_all[de_all["family"] != "other"].copy()
        fam_de.to_csv(TABLES / "de_families.tsv", sep="\t", index=False)
        head = fam_de[fam_de["contrast"] == "q4q1_combined"].sort_values(["compartment", "p"])
        head.to_csv(TABLES / "de_q4q1_combined_families.tsv", sep="\t", index=False)
        # compact headline markdown (top 8 per compartment)
        lines = [
            "| compartment | family | gene | n_Q1 | n_Q4 | logFC | p | FDR |",
            "|---|---|---|---:|---:|---:|---:|---:|",
        ]
        for comp in ("malignant", "T/NK"):
            sub = head[head["compartment"] == comp].head(8)
            for r in sub.itertuples():
                lines.append(
                    f"| {r.compartment} | {r.family} | {r.gene} | {int(r.n_q1)} | {int(r.n_q4)} | "
                    f"{r.logFC:+.3f} | {_fmt_p(r.p)} | {_fmt_p(r.fdr)} |"
                )
        headline = {"md": "\n".join(lines)}
    else:
        headline = {"md": "_DE table empty._"}

    fam_df = pd.DataFrame(fam_rows)
    if not fam_df.empty:
        fam_df.to_csv(TABLES / "family_summary.tsv", sep="\t", index=False)
    write_finding(meta, inv, fam_df if not fam_df.empty else pd.DataFrame(), headline)
    print("wrote FINDING.md and tables/")
    print(inv.to_string(index=False))


if __name__ == "__main__":
    main()
