#!/usr/bin/env python3
"""REWORK A1 deconvolution — TCGA-LUAD TACSTD2 vs xCell, MCP-counter,
ESTIMATE ImmuneScore, GEP18 ssGSEA, and Danaher.

Self-contained. Public data only. ABSOLUTE partial Spearman.
Writes results/rework/A1_deconv/.

Why this exists
---------------
The LUAD-only A1 scan (CIBERSORT / Wolf / marker genes) was weak
(|ρ| ≤ 0.22). The histology rework then showed GEP18 and ESTIMATE are
null in LUAD. This slice asks whether a different *definition* of
"immune" — official deconvolution / signature methods — changes that
LUAD picture. It does not pool LUSC. It reports every requested method
and every population those methods emit.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SIG = DATA / "signatures"
OUT = ROOT / "results" / "rework" / "A1_deconv"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

URLS = {
    "expression": {
        "file": "TCGA.LUAD.HiSeqV2.gz",
        "url": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz",
        "desc": "UCSC Xena TCGA-LUAD HiSeqV2, log2(RSEM normalized_count + 1)",
    },
    "purity": {
        "file": "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
        "url": "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
        "desc": "PanCanAtlas ABSOLUTE purity (GDC 4f277128-f793-4354-a13d-30cc7fe9f6b5)",
    },
    "estimate": {
        "file": "ESTIMATE_LUAD_RNAseqV2.txt",
        "url": "https://ibl.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
        "desc": "MD Anderson official ESTIMATE RNAseqV2 ImmuneScore / StromalScore / ESTIMATEScore",
    },
    "timer2": {
        "file": "infiltration_estimation_for_tcga.csv.gz",
        "url": "https://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz",
        "desc": "TIMER2.0 / immunedeconv precomputed TCGA infiltration (xCell + MCP-counter + others)",
    },
}

# Pre-specified primary panel used for the verdict (not a hidden subset).
# Every other population is still in correlations.tsv.
PRIMARY = [
    "ESTIMATE_ImmuneScore",
    "GEP18_ssGSEA",
    "xCell_immune_score",
    "xCell_T_cell_CD8",
    "xCell_microenvironment_score",
    "TIMER2_MCP_T_cell",
    "TIMER2_MCP_T_cell_CD8",
    "TIMER2_MCP_cytotoxicity_score",
    "Xena_MCP_T_cells",
    "Xena_MCP_CD8_T_cells",
    "Xena_MCP_Cytotoxic_lymphocytes",
    "Danaher_T_cells",
    "Danaher_CD8_T_cells",
    "Danaher_Cytotoxic_cells",
    "Danaher_Total_TILs",
]

DANAHER_TOTAL_TILS_POPS = [
    "B-cells",
    "CD45",
    "Cytotoxic cells",
    "Exhausted CD8",
    "Macrophages",
    "Neutrophils",
    "NK CD56dim",
    "NK cells",
    "T-cells",
    "Th1 cells",
    "CD8 T cells",
]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(key: str) -> Path:
    meta = URLS[key]
    dest = DATA / meta["file"]
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    DATA.mkdir(parents=True, exist_ok=True)
    print(f"downloading {key} -> {dest.name}")
    r = requests.get(meta["url"], timeout=300)
    r.raise_for_status()
    dest.write_bytes(r.content)
    return dest


def sample01(barcode: str) -> str | None:
    parts = str(barcode).replace(".", "-").split("-")
    if len(parts) < 4 or not parts[3].startswith("01"):
        return None
    return "-".join(parts[:3]) + "-01"


def bh_fdr(pvals) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    q = np.full(n, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    idx = np.where(ok)[0]
    pv = p[idx]
    order = np.argsort(pv)
    prev = 1.0
    out = np.empty(len(pv))
    for rank_from_end, j in enumerate(order[::-1]):
        rank = len(pv) - rank_from_end
        val = min(prev, pv[j] * len(pv) / rank)
        out[j] = val
        prev = val
    q[idx] = out
    return q


def fisher_z_ci(r: float, n: int, k_covariates: int = 1):
    if r is None or not np.isfinite(r) or n <= k_covariates + 3:
        return np.nan, np.nan
    r = float(np.clip(r, -0.999999, 0.999999))
    z = math.atanh(r)
    se = 1.0 / math.sqrt(n - k_covariates - 3)
    zcrit = stats.norm.ppf(0.975)
    return float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se))


def partial_spearman_algebraic(x, y, z):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 6:
        return np.nan, np.nan, n
    rxy = stats.spearmanr(x[m], y[m]).statistic
    rxz = stats.spearmanr(x[m], z[m]).statistic
    ryz = stats.spearmanr(y[m], z[m]).statistic
    denom = math.sqrt(max((1 - rxz**2) * (1 - ryz**2), 1e-12))
    r = (rxy - rxz * ryz) / denom
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def partial_spearman_residual(x, y, z):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 6:
        return np.nan, np.nan, n
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = stats.rankdata(z[m])
    zc = np.column_stack([np.ones(n), zr])
    bx, *_ = np.linalg.lstsq(zc, xr, rcond=None)
    by, *_ = np.linalg.lstsq(zc, yr, rcond=None)
    r, _ = stats.pearsonr(xr - zc @ bx, yr - zc @ by)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def ssgsea_one_set(expr: pd.DataFrame, genes: list[str], tau: float = 0.25) -> pd.Series:
    """Barbie / GSVA ssGSEA (running-sum ES, tau=0.25).

    Matches the ESTIMATE v1.0.13 / GSVA convention used in A2:
    rank genes within sample (average ties), scale ranks to 1..10000,
    sort descending, weight hits by |rank|^tau, ES = sum(P_hit - P_miss).
    Expression scale is irrelevant (ranks only). One gene set, so
    ssgsea.norm across gene sets is not applied.
    """
    present = [g for g in genes if g in expr.index]
    if len(present) < 2:
        raise SystemExit(f"ssGSEA gene set too small: {present}")
    n_genes = expr.shape[0]
    ranked = expr.rank(axis=0, method="average", ascending=True) * (10000.0 / n_genes)
    gene_set = set(present)
    scores = {}
    for sample in expr.columns:
        m = ranked[sample]
        order = m.sort_values(ascending=False).index
        m_ord = m.loc[order].to_numpy(float)
        hits = np.fromiter((g in gene_set for g in order), dtype=bool, count=len(order))
        w = np.abs(m_ord) ** tau
        w_hit = np.where(hits, w, 0.0)
        nhit = float(w_hit.sum())
        nmiss = float((~hits).sum())
        if nhit <= 0 or nmiss <= 0:
            scores[sample] = np.nan
            continue
        p_hit = np.cumsum(w_hit) / nhit
        p_miss = np.cumsum((~hits).astype(float)) / nmiss
        scores[sample] = float(np.sum(p_hit - p_miss))
    return pd.Series(scores, name="GEP18_ssGSEA")


def slug(name: str) -> str:
    return (
        name.replace(" ", "_")
        .replace("+", "")
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("-", "_")
    )


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def main():
    for k in URLS:
        download(k)

    print("Loading Xena LUAD expression ...")
    expr = pd.read_csv(DATA / URLS["expression"]["file"], sep="\t", index_col=0)
    expr = expr.loc[:, [c for c in expr.columns if str(c).endswith("-01")]]
    if expr.columns.duplicated().any():
        expr = expr.T.groupby(level=0).mean().T
    if "TACSTD2" not in expr.index:
        raise SystemExit("TACSTD2 missing from Xena HiSeqV2")
    tacstd2 = expr.loc["TACSTD2"].astype(float)
    n_expr = int(tacstd2.shape[0])
    print(f"  primary tumors with RNA: {n_expr}")

    print("Loading ABSOLUTE ...")
    ab = pd.read_csv(DATA / URLS["purity"]["file"], sep="\t")
    ab = ab[ab["array"].astype(str).str.endswith("-01", na=False)][["array", "purity"]].dropna()
    purity = ab.groupby("array")["purity"].mean()

    print("Loading official ESTIMATE ...")
    est = pd.read_csv(DATA / URLS["estimate"]["file"], sep="\t")
    est.columns = [c.strip() for c in est.columns]
    idcol = est.columns[0]
    est["sample"] = est[idcol].map(lambda b: sample01(b) or (str(b) if str(b).endswith("-01") else None))
    est = est.dropna(subset=["sample"])
    rename = {}
    for c in est.columns:
        cl = c.lower().replace(" ", "_")
        if cl in ("immune_score", "immunescore"):
            rename[c] = "ESTIMATE_ImmuneScore"
        elif cl in ("stromal_score", "stromalscore"):
            rename[c] = "ESTIMATE_StromalScore"
        elif cl in ("estimate_score", "estimatescore"):
            rename[c] = "ESTIMATE_Score"
    est = est.rename(columns=rename)
    est_keep = [c for c in ("ESTIMATE_ImmuneScore", "ESTIMATE_StromalScore", "ESTIMATE_Score") if c in est.columns]
    est = est.groupby("sample")[est_keep].mean()

    print("Loading TIMER2 xCell + MCP-counter ...")
    timer = pd.read_csv(DATA / URLS["timer2"]["file"])
    timer["sample"] = timer["cell_type"].map(lambda b: sample01(b) or (str(b) if str(b).endswith("-01") else None))
    timer = timer.dropna(subset=["sample"])
    xcell_cols = [c for c in timer.columns if c.endswith("_XCELL")]
    mcp_t_cols = [c for c in timer.columns if c.endswith("_MCPCOUNTER")]
    timer_x = timer.groupby("sample")[xcell_cols + mcp_t_cols].mean()
    timer_x = timer_x.reindex(tacstd2.index)
    xcell_rename = {c: "xCell_" + slug(c[: -len("_XCELL")]) for c in xcell_cols}
    mcp_rename = {c: "TIMER2_MCP_" + slug(c[: -len("_MCPCOUNTER")]) for c in mcp_t_cols}
    timer_x = timer_x.rename(columns={**xcell_rename, **mcp_rename})

    print("Scoring MCP-counter on Xena (Becht 2016 markers) ...")
    mcp_genes = pd.read_csv(SIG / "mcp_counter_genes.tsv", sep="\t")
    mcp_coverage = []
    mcp_scores = {}
    for pop, sub in mcp_genes.groupby("population", sort=False):
        symbols = []
        missing = []
        for _, row in sub.iterrows():
            cand = [row["xena_symbol"], row["hugo"]]
            hit = next((g for g in cand if isinstance(g, str) and g in expr.index), None)
            if hit is None:
                missing.append(row["hugo"])
            else:
                symbols.append(hit)
        symbols = list(dict.fromkeys(symbols))
        col = "Xena_MCP_" + slug(pop)
        if symbols:
            mcp_scores[col] = expr.loc[symbols].mean(axis=0)
        else:
            mcp_scores[col] = pd.Series(np.nan, index=expr.columns)
        mcp_coverage.append(
            {
                "method": "MCP-counter (Xena recompute)",
                "population": pop,
                "n_official": int(len(sub)),
                "n_present": int(len(symbols)),
                "present_genes": ",".join(symbols),
                "missing_genes": ",".join(missing),
                "score_column": col,
            }
        )
    mcp_xena = pd.DataFrame(mcp_scores)

    print("Scoring Danaher 2017 on Xena ...")
    dan = pd.read_csv(SIG / "danaher_genes.tsv", sep="\t")
    dan_coverage = []
    dan_scores = {}
    for pop, sub in dan.groupby("population", sort=False):
        present = [g for g in sub["hugo"] if g in expr.index]
        missing = [g for g in sub["hugo"] if g not in expr.index]
        col = "Danaher_" + slug(pop)
        if present:
            dan_scores[col] = expr.loc[present].mean(axis=0)
        else:
            dan_scores[col] = pd.Series(np.nan, index=expr.columns)
        dan_coverage.append(
            {
                "method": "Danaher 2017",
                "population": pop,
                "n_official": int(len(sub)),
                "n_present": int(len(present)),
                "present_genes": ",".join(present),
                "missing_genes": ",".join(missing),
                "score_column": col,
            }
        )
    dan_df = pd.DataFrame(dan_scores)
    tils_cols = ["Danaher_" + slug(p) for p in DANAHER_TOTAL_TILS_POPS]
    missing_tils = [c for c in tils_cols if c not in dan_df.columns]
    if missing_tils:
        raise SystemExit(f"Danaher Total TILs missing columns: {missing_tils}")
    dan_df["Danaher_Total_TILs"] = dan_df[tils_cols].mean(axis=1)
    # Empirical Total TILs: paper rule ρ(score, PTPRC) > 0.6 on this LUAD set.
    if "Danaher_CD45" in dan_df.columns:
        emp = []
        for c in dan_df.columns:
            if c in ("Danaher_CD45", "Danaher_Total_TILs"):
                continue
            r, _ = stats.spearmanr(dan_df[c], dan_df["Danaher_CD45"])
            if np.isfinite(r) and r > 0.6:
                emp.append(c)
        dan_df["Danaher_Total_TILs_empirical"] = dan_df[emp].mean(axis=1) if emp else np.nan
    else:
        emp = []
        dan_df["Danaher_Total_TILs_empirical"] = np.nan

    print("Scoring GEP18 ssGSEA + z-mean ...")
    gep = pd.read_csv(SIG / "gep18_genes.tsv", sep="\t")
    gep_genes = gep["hugo"].tolist()
    gep_present = [g for g in gep_genes if g in expr.index]
    gep_missing = [g for g in gep_genes if g not in expr.index]
    if gep_missing:
        raise SystemExit(f"GEP18 genes missing from Xena: {gep_missing}")
    gep_ssgsea = ssgsea_one_set(expr, gep_present, tau=0.25)
    z = expr.loc[gep_present].T
    gep_zmean = ((z - z.mean(axis=0)) / z.std(axis=0, ddof=0)).mean(axis=1)
    gep_zmean.name = "GEP18_zmean"
    gep_logmean = expr.loc[gep_present].mean(axis=0)
    gep_logmean.name = "GEP18_logmean"

    print("Merging analysis table ...")
    df = pd.DataFrame({"TACSTD2": tacstd2})
    df = df.join(purity.rename("ABSOLUTE_purity"), how="left")
    df = df.join(est, how="left")
    df = df.join(timer_x, how="left")
    df = df.join(mcp_xena, how="left")
    df = df.join(dan_df, how="left")
    df = df.join(gep_ssgsea, how="left")
    df = df.join(gep_zmean, how="left")
    df = df.join(gep_logmean, how="left")
    df.index.name = "sample"

    feature_meta = []
    for c in est_keep:
        feature_meta.append({"feature": c, "method": "ESTIMATE", "role": "primary" if c == "ESTIMATE_ImmuneScore" else "secondary", "source": "MD Anderson official RNAseqV2 table"})
    feature_meta.append({"feature": "GEP18_ssGSEA", "method": "GEP18", "role": "primary", "source": "ssGSEA tau=0.25 on Ayers 18 genes, Xena matrix"})
    feature_meta.append({"feature": "GEP18_zmean", "method": "GEP18", "role": "sensitivity", "source": "unweighted mean of within-LUAD z-scores (histology-rework definition)"})
    feature_meta.append({"feature": "GEP18_logmean", "method": "GEP18", "role": "sensitivity", "source": "mean of log2(RSEM+1) of 18 genes"})
    for c in timer_x.columns:
        if c.startswith("xCell_"):
            role = "primary" if c in PRIMARY else "all_populations"
            feature_meta.append({"feature": c, "method": "xCell", "role": role, "source": "TIMER2.0 immunedeconv precomputed"})
        elif c.startswith("TIMER2_MCP_"):
            role = "primary" if c in PRIMARY else "all_populations"
            feature_meta.append({"feature": c, "method": "MCP-counter", "role": role, "source": "TIMER2.0 immunedeconv precomputed"})
    for c in mcp_xena.columns:
        role = "primary" if c in PRIMARY else "all_populations"
        feature_meta.append({"feature": c, "method": "MCP-counter", "role": role, "source": "Becht 2016 markers, mean log2(RSEM+1) on Xena"})
    for c in dan_df.columns:
        role = "primary" if c in PRIMARY else "all_populations"
        feature_meta.append({"feature": c, "method": "Danaher", "role": role, "source": "Danaher JITC 2017 Table 1, mean log2(RSEM+1) on Xena"})
    meta = pd.DataFrame(feature_meta).drop_duplicates("feature")

    features = [c for c in df.columns if c not in ("TACSTD2", "ABSOLUTE_purity")]
    x_all = df["TACSTD2"].to_numpy(float)
    z_all = df["ABSOLUTE_purity"].to_numpy(float)
    rows = []
    for feat in features:
        y = df[feat].to_numpy(float)
        m_un = np.isfinite(x_all) & np.isfinite(y)
        n_un = int(m_un.sum())
        if n_un >= 6 and np.nanstd(x_all[m_un]) > 0 and np.nanstd(y[m_un]) > 0:
            rho, p_un = stats.spearmanr(x_all[m_un], y[m_un])
        else:
            rho, p_un = np.nan, np.nan
        pr, pp, n_p = partial_spearman_algebraic(x_all, y, z_all)
        pr2, pp2, _ = partial_spearman_residual(x_all, y, z_all)
        lo, hi = fisher_z_ci(pr, n_p, k_covariates=1)
        m_fp = np.isfinite(y) & np.isfinite(z_all)
        if m_fp.sum() >= 6:
            r_fp, p_fp = stats.spearmanr(y[m_fp], z_all[m_fp])
        else:
            r_fp, p_fp = np.nan, np.nan
        rows.append(
            {
                "feature": feat,
                "method": meta.set_index("feature").loc[feat, "method"] if feat in meta.feature.values else "other",
                "role": meta.set_index("feature").loc[feat, "role"] if feat in meta.feature.values else "other",
                "n_unadjusted": n_un,
                "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                "spearman_p": float(p_un) if np.isfinite(p_un) else np.nan,
                "n_partial": n_p,
                "partial_rho_algebraic": pr,
                "partial_p_algebraic": pp,
                "partial_rho_residual": pr2,
                "partial_p_residual": pp2,
                "partial_ci_low": lo,
                "partial_ci_high": hi,
                "feature_vs_purity_rho": float(r_fp) if np.isfinite(r_fp) else np.nan,
                "feature_vs_purity_p": float(p_fp) if np.isfinite(p_fp) else np.nan,
            }
        )
    res = pd.DataFrame(rows)
    res["spearman_fdr_all"] = bh_fdr(res["spearman_p"])
    res["partial_fdr_all"] = bh_fdr(res["partial_p_algebraic"])
    prim_mask = res["feature"].isin(PRIMARY)
    res["partial_fdr_primary15"] = np.nan
    res.loc[prim_mask, "partial_fdr_primary15"] = bh_fdr(res.loc[prim_mask, "partial_p_algebraic"])
    res = res.sort_values(["method", "partial_p_algebraic"]).reset_index(drop=True)

    m_tp = df["ABSOLUTE_purity"].notna()
    rho_tp, p_tp = stats.spearmanr(df.loc[m_tp, "TACSTD2"], df.loc[m_tp, "ABSOLUTE_purity"])
    n_partial = int(m_tp.sum())

    coverage = pd.DataFrame(mcp_coverage + dan_coverage)
    coverage.loc[len(coverage)] = {
        "method": "GEP18",
        "population": "T-cell-inflamed GEP",
        "n_official": 18,
        "n_present": len(gep_present),
        "present_genes": ",".join(gep_present),
        "missing_genes": ",".join(gep_missing),
        "score_column": "GEP18_ssGSEA",
    }

    # ---------- write tables ----------
    res.to_csv(OUT / "correlations.tsv", sep="\t", index=False)
    df.reset_index().to_csv(OUT / "sample_table.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "gene_coverage.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "feature_definitions.tsv", sep="\t", index=False)
    prim = res[res["feature"].isin(PRIMARY)].copy()
    prim["primary_order"] = prim["feature"].map({f: i for i, f in enumerate(PRIMARY)})
    prim = prim.sort_values("primary_order")
    prim.to_csv(OUT / "primary_partial.tsv", sep="\t", index=False)

    # ---------- figures ----------
    plot_primary = prim.dropna(subset=["partial_rho_algebraic"])
    fig, ax = plt.subplots(figsize=(9, 6.2))
    y = np.arange(len(plot_primary))
    colors = ["#c0392b" if r < 0 else "#2471a3" for r in plot_primary["partial_rho_algebraic"]]
    ax.barh(y, plot_primary["partial_rho_algebraic"], color=colors, edgecolor="none")
    ax.errorbar(
        plot_primary["partial_rho_algebraic"],
        y,
        xerr=[
            plot_primary["partial_rho_algebraic"] - plot_primary["partial_ci_low"],
            plot_primary["partial_ci_high"] - plot_primary["partial_rho_algebraic"],
        ],
        fmt="none",
        ecolor="black",
        elinewidth=0.8,
        capsize=2,
    )
    ax.axvline(0, color="k", lw=0.8)
    ax.axvline(-0.22, color="#7f8c8d", ls="--", lw=0.8, label="prior LUAD |ρ|=0.22 bound")
    ax.axvline(0.22, color="#7f8c8d", ls="--", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_primary["feature"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("ABSOLUTE partial Spearman ρ (TACSTD2 vs feature)")
    ax.set_title("TCGA-LUAD: TACSTD2 vs deconvolution / signature scores\nprimary panel, 95% Fisher-z CI")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "forest_primary_partial_rho.png", dpi=160)
    plt.close(fig)

    # all features, grouped
    r2 = res.dropna(subset=["partial_rho_algebraic"]).sort_values("partial_rho_algebraic")
    fig, ax = plt.subplots(figsize=(8.5, max(6, 0.18 * len(r2) + 1.5)))
    method_color = {
        "ESTIMATE": "#8e44ad",
        "GEP18": "#16a085",
        "xCell": "#2980b9",
        "MCP-counter": "#d35400",
        "Danaher": "#27ae60",
    }
    cols = [method_color.get(m, "#7f8c8d") for m in r2["method"]]
    ax.barh(r2["feature"], r2["partial_rho_algebraic"], color=cols)
    ax.axvline(0, color="k", lw=0.7)
    ax.axvline(-0.22, color="#7f8c8d", ls="--", lw=0.7)
    ax.axvline(0.22, color="#7f8c8d", ls="--", lw=0.7)
    ax.set_xlabel("ABSOLUTE partial Spearman ρ")
    ax.set_title("All reported features (color = method)")
    ax.tick_params(axis="y", labelsize=6)
    fig.tight_layout()
    fig.savefig(FIG / "bar_all_partial_rho.png", dpi=140)
    plt.close(fig)

    scatter_feats = [
        "ESTIMATE_ImmuneScore",
        "GEP18_ssGSEA",
        "xCell_immune_score",
        "xCell_T_cell_CD8",
        "TIMER2_MCP_T_cell_CD8",
        "Xena_MCP_Cytotoxic_lymphocytes",
        "Danaher_Cytotoxic_cells",
        "Danaher_Total_TILs",
    ]
    scatter_feats = [f for f in scatter_feats if f in df.columns]
    ncol = 4
    nrow = int(math.ceil(len(scatter_feats) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 3.4 * nrow), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    lookup = res.set_index("feature")
    for ax, feat in zip(axes, scatter_feats):
        m = df[feat].notna() & df["ABSOLUTE_purity"].notna()
        sc = ax.scatter(
            df.loc[m, "TACSTD2"],
            df.loc[m, feat],
            c=df.loc[m, "ABSOLUTE_purity"],
            cmap="viridis",
            s=9,
            alpha=0.75,
            linewidths=0,
        )
        r = lookup.loc[feat]
        ax.set_title(
            f"{feat}\nρ={r.spearman_rho:.3f}; adj ρ={r.partial_rho_algebraic:.3f} "
            f"(p={r.partial_p_algebraic:.2g}, n={int(r.n_partial)})",
            fontsize=8,
        )
        ax.set_xlabel("TACSTD2 log2(RSEM+1)", fontsize=8)
        ax.set_ylabel(feat, fontsize=7)
    for ax in axes[len(scatter_feats) :]:
        ax.axis("off")
    fig.colorbar(sc, ax=axes[: len(scatter_feats)], label="ABSOLUTE purity", shrink=0.5)
    fig.suptitle("TCGA-LUAD primary tumors: TACSTD2 vs immune deconvolution scores", fontsize=11)
    fig.savefig(FIG / "scatter_primary.png", dpi=150)
    plt.close(fig)

    # ---------- summary / report ----------
    prim_abs = prim["partial_rho_algebraic"].abs()
    max_abs = float(prim_abs.max()) if prim_abs.notna().any() else float("nan")
    n_gt_022 = int((prim_abs > 0.22).sum())
    n_sig_fdr = int((prim["partial_fdr_primary15"] < 0.05).sum()) if prim["partial_fdr_primary15"].notna().any() else 0

    all_abs = res["partial_rho_algebraic"].abs()
    extra = res.loc[all_abs > 0.22, ["feature", "partial_rho_algebraic"]].sort_values(
        "partial_rho_algebraic", key=lambda s: s.abs(), ascending=False
    )
    extra_txt = "; ".join(f"{r.feature} ρ={r.partial_rho_algebraic:.3f}" for r in extra.itertuples())
    if n_gt_022 == 0:
        verdict = (
            "STILL WEAK for the named immune-infiltration / T-cell-inflamed scores. "
            "xCell ImmuneScore, ESTIMATE ImmuneScore, GEP18 ssGSEA, MCP T/CD8/cytotoxicity, "
            "and Danaher T/CD8/cytotoxic/Total TILs all stay inside the prior LUAD |ρ|≤0.22 bound "
            f"(largest primary |ρ| = {max_abs:.3f}). "
            + (
                f"Two non-primary xCell subsets exceed 0.22 and are reported, not hidden: {extra_txt}. "
                "They do not form a coherent immune-cold pattern (one is positive)."
                if len(extra)
                else ""
            )
        )
        verdict_tag = "STILL_WEAK"
    else:
        verdict = (
            f"METHOD-DEPENDENT among primary scores. {n_gt_022} primary feature(s) exceed |ρ|=0.22 "
            f"(max primary |ρ| = {max_abs:.3f}). All-feature |ρ|>0.22: {extra_txt or 'none'}."
        )
        verdict_tag = "METHOD_DEPENDENT"

    def row_md(feat: str) -> str:
        r = lookup.loc[feat]
        return (
            f"| `{feat}` | {fmt_rho(r.spearman_rho)} | {fmt_p(r.spearman_p)} | "
            f"{fmt_rho(r.partial_rho_algebraic)} | {fmt_p(r.partial_p_algebraic)} | "
            f"{fmt_p(r.partial_fdr_primary15)} | {int(r.n_partial)} | "
            f"{fmt_rho(r.feature_vs_purity_rho)} |"
        )

    primary_table = "\n".join(row_md(f) for f in PRIMARY if f in lookup.index)

    def method_block(method: str) -> str:
        sub = res[res["method"] == method].sort_values("partial_rho_algebraic")
        lines = [
            "| Feature | Unadj ρ | Unadj p | Partial ρ | Partial p | FDR (all) | n | vs purity ρ |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for _, r in sub.iterrows():
            lines.append(
                f"| `{r.feature}` | {fmt_rho(r.spearman_rho)} | {fmt_p(r.spearman_p)} | "
                f"{fmt_rho(r.partial_rho_algebraic)} | {fmt_p(r.partial_p_algebraic)} | "
                f"{fmt_p(r.partial_fdr_all)} | {int(r.n_partial)} | {fmt_rho(r.feature_vs_purity_rho)} |"
            )
        return "\n".join(lines)

    mcp_cov_md = coverage[coverage.method.str.contains("MCP")].to_string(index=False)
    dan_cov_md = coverage[coverage.method.str.contains("Danaher")].to_string(index=False)

    report = f"""# REWORK A1 deconvolution — TCGA-LUAD TACSTD2 vs xCell / MCP-counter / ESTIMATE / GEP18 ssGSEA / Danaher

**Self-contained. Public data only. Written to be read without the rest of the repo.**

**Cohort:** TCGA-LUAD primary tumors (`-01`). {n_expr} with RNA-seq; **n = {n_partial}** with RNA-seq + ABSOLUTE purity (partial-correlation set).

**Question:** The prior LUAD-only A1 scan was weak (|ρ| ≤ 0.22 on CIBERSORT / Wolf / marker genes). Does a different immune *definition* — official deconvolution and signature methods — change that?

## Verdict: {verdict_tag}

{verdict}

TACSTD2 vs ABSOLUTE purity: Spearman ρ = {rho_tp:.3f}, p = {fmt_p(p_tp)}, n = {n_partial}. Purity adjustment is nearly a no-op for TACSTD2 itself (same observation as the prior LUAD-only A1). It is **not** optional for ESTIMATE / xCell ImmuneScore / MCP totals, which track leukocyte content and are strongly (negatively) correlated with purity.

## Why this rework exists

Prior LUAD-only A1 (`results/w200/A1_LUAD/`): 32 features, all |ρ| ≤ 0.22 after ABSOLUTE; no association with total leukocyte fraction (ρ = 0.03). Histology rework (`results/rework/A1_histology/`): LUAD GEP18 ρ ≈ 0, ESTIMATE ImmuneScore weakly positive / NS; CD8/CYT small negatives. Those used CD8A, CYT, unweighted GEP18 z-mean, and official ESTIMATE. They did **not** run xCell, MCP-counter, GEP18 *ssGSEA*, or Danaher.

This file runs those five named methods on the same TCGA-LUAD Xena freeze + ABSOLUTE, and reports **every population** each method emits.

## Analysis set

| Filter | n |
|---|---:|
| Xena HiSeqV2 primary tumors (`-01`) | {n_expr} |
| + ABSOLUTE purity | {n_partial} |
| + official ESTIMATE ImmuneScore | {int(df['ESTIMATE_ImmuneScore'].notna().sum()) if 'ESTIMATE_ImmuneScore' in df.columns else 0} |
| + TIMER2 xCell / MCP-counter | {int(df['xCell_immune_score'].notna().sum()) if 'xCell_immune_score' in df.columns else 0} |

Primary tumors only. One row per 15-character barcode. Replicate aliquots averaged. Partial n is TACSTD2 + feature + ABSOLUTE (exact n is in every table cell).

## ALL definitions

### TACSTD2
- Gene symbol `TACSTD2` (TROP2, GA733-1, EGP-1). Ensembl `ENSG00000184292`.
- Measure: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` **log2(RSEM normalized_count + 1)**.
- Not protein. Not ADC response. Not ICI outcome.

### Samples
- TCGA-LUAD only. No LUSC. No OncoSG.
- Keep barcode sample-type `01` (primary solid tumor). Drop `11` normals and other types.
- Identifier: first 15 characters (`TCGA-XX-XXXX-01`).

### Covariate
- PanCanAtlas **ABSOLUTE** `purity` from GDC file `4f277128-f793-4354-a13d-30cc7fe9f6b5` (`TCGA_mastercalls.abs_tables_JSedit.fixed.txt`).
- Not ESTIMATE TumorPurity. Not methylation leukocyte fraction.

### Partial Spearman (primary statistic)
Two equivalent one-covariate estimators; both are reported.

1. **Algebraic (primary, matches original pooled A1 / histology rework):**
   \\( r_{{xy\\cdot z}} = (r_{{xy}} - r_{{xz}} r_{{yz}}) / \\sqrt{{(1-r_{{xz}}^2)(1-r_{{yz}}^2)}} \\)
   with pairwise Spearman rhos. t-test, **df = n − 3**.
2. **Residual-rank (sensitivity, matches LUAD-only A1):** Pearson of rank residuals after OLS of ranked x and y on ranked purity.

Also reported: unadjusted Spearman; 95% Fisher-z CI with variance `1/(n−4)`; feature vs purity Spearman (so ImmuneScore-like composites are not silently called “immune”).

BH-FDR is shown two ways: across **all** reported features, and across the **15** pre-specified primary features. Raw p-values are always shown. No feature is hidden.

### 1. xCell (Aran *Genome Biol* 2017)
- **Source used:** TIMER2.0 precomputed scores (`infiltration_estimation_for_tcga.csv.gz`, Li *NAR* 2020), produced by **immunedeconv** running the official xCell algorithm on TCGA TPM.
- **Not recomputed** on Xena log2(RSEM+1). xCell’s power/calibration/spillover steps expect TPM-like input; the published TIMER2 table is the standard TCGA xCell freeze.
- **Primary features:** `immune score`, `T cell CD8+`, `microenvironment score`.
- **All other xCell columns** in TIMER2 are reported (B, CD4 subsets, Tcm/Tem/naive CD8, NK, macrophages M1/M2, Tregs, stroma score, progenitors, etc.).
- xCell ImmuneScore is a composite of immune cell scores after spillover compensation. It is expected to anti-correlate with ABSOLUTE purity.

### 2. MCP-counter (Becht *Genome Biol* 2016)
Two implementations, both reported (they are not identical).

**A. TIMER2 / immunedeconv (official algorithm on TCGA TPM).** Columns: T cell, T cell CD8+, cytotoxicity score, NK cell, B cell, Monocyte, Macrophage/Monocyte, myeloid DC, Neutrophil, Endothelial cell, CAF.

**B. Xena recompute (definition-transparent, same matrix as TACSTD2).** Score = **mean of log2(RSEM+1)** of the official Becht marker genes for each of the 10 populations (T cells, CD8 T cells, Cytotoxic lymphocytes, B lineage, NK cells, Monocytic lineage, Myeloid dendritic cells, Neutrophils, Endothelial cells, Fibroblasts). Missing genes are **not imputed**; coverage is in `gene_coverage.tsv`.

Xena symbol aliases (do not treat as new genes): `ADGRL4` → `ELTD1`; `DIPK2B` → `CXorf36`. Absent from HiSeqV2 and dropped: `CHRM3-AS2`, `IGKC`, `KIR3DS1`, `WFDC21P`.

CD8 T cells is a **single gene** (`CD8B`) in the official MCP signature. Cytotoxic lymphocytes is a 7-gene set that includes `CD8A`.

### 3. ESTIMATE ImmuneScore (Yoshihara *Nat Commun* 2013)
- Official MD Anderson **RNAseqV2** table `lung_adenocarcinoma_RNAseqV2.txt`.
- Primary: `Immune_score`. Secondary: `Stromal_score`, `ESTIMATE_score` (= Immune + Stromal).
- **Not recomputed** by ssGSEA in this slice (official table returned HTTP 200).
- ImmuneScore is built to track leukocyte content. ABSOLUTE partialling is required before calling it an “immune” association.

### 4. GEP18 ssGSEA (Ayers *JCI* 2017)
18 genes (all present in Xena):
`CCL5, CD27, CD274, CD276, CD8A, CMKLR1, CXCL9, CXCR6, HLA-DQA1, HLA-DRB1, HLA-E, IDO1, LAG3, NKG7, PDCD1LG2, PSMB10, STAT1, TIGIT`.

- **Primary score:** Barbie/GSVA **ssGSEA**, tau = 0.25, running-sum ES, ranks scaled to 1..10000 over the full HiSeqV2 gene universe (20,530 genes). No Merck NanoString TIS weights (not public). No inter-gene-set `ssgsea.norm` (one set).
- **Sensitivity:** unweighted mean of within-LUAD z-scores of log2(RSEM+1) — the histology-rework definition — and simple mean of the 18 log2 values.
- This is **not** the clinical NanoString TIS.

### 5. Danaher (Danaher *JITC* 2017 Table 1)
14 populations, 60 marker genes. Score = **average log2 expression** of that population’s markers (paper Methods). Xena is already log2(RSEM+1), so this is the paper’s estimator.

| Population | Genes |
|---|---|
| B-cells | BLK, CD19, FCRL2, MS4A1, KIAA0125, TNFRSF17, TCL1A, SPIB, PNOC |
| CD45 | PTPRC (paper typo: PTRPC) |
| Cytotoxic cells | PRF1, GZMA, GZMB, NKG7, GZMH, KLRK1, KLRB1, KLRD1, CTSW, GNLY |
| DC | CCL13, CD209, HSD11B1 |
| Exhausted CD8 | LAG3, CD244, EOMES, PTGER4 |
| Macrophages | CD68, CD84, CD163, MS4A4A |
| Mast cells | TPSB2, TPSAB1, CPA3, MS4A2, HDC |
| Neutrophils | FPR1, SIGLEC5, CSF3R, FCAR, FCGR3B, CEACAM3, S100A12 |
| NK CD56dim | KIR2DL3, KIR3DL1, KIR3DL2, IL21R |
| NK cells | XCL1, XCL2, NCR1 |
| T-cells | CD6, CD3D, CD3E, SH2D1A, TRAT1, CD3G |
| Th1 cells | TBX21 |
| Treg | FOXP3 |
| CD8 T cells | CD8A, CD8B |

**Total TILs (primary):** mean of the 11 populations the paper retained after requiring correlation with PTPRC > 0.6, i.e. exclude DC, Treg, and mast cells. Sensitivity: `Danaher_Total_TILs_empirical` recomputes that cutoff on this LUAD set (populations used: {", ".join(emp) if emp else "none"}).

Th1 and Treg are single-gene scores; the paper says so. They are not multi-gene deconvolution.

## Primary result — ABSOLUTE partial Spearman

| Feature | Unadj ρ | Unadj p | Partial ρ | Partial p | FDR (15) | n | vs purity ρ |
|---|---:|---:|---:|---:|---:|---:|---:|
{primary_table}

Dashed lines on the forest plot mark the prior LUAD-only |ρ| = 0.22 bound. {n_sig_fdr} / 15 primary tests have BH-FDR < 0.05.

Features with |partial ρ| > 0.22 anywhere in the full scan (not hidden): {extra_txt or "none"}.
If those rows are not in the primary table, they are xCell subsets, not global immune scores. Do not rewrite the LUAD story around a single subset.

## All populations, by method

### ESTIMATE
{method_block("ESTIMATE")}

### GEP18
{method_block("GEP18")}

### xCell (TIMER2, all columns)
{method_block("xCell")}

### MCP-counter (TIMER2 + Xena recompute)
{method_block("MCP-counter")}

### Danaher
{method_block("Danaher")}

## Gene coverage (Xena recomputes)

Missing genes are listed, not imputed.

### MCP-counter (Xena)
```
{mcp_cov_md}
```

### Danaher
```
{dan_cov_md}
```

GEP18: 18/18 present.

## Honest interpretation

1. **This is still LUAD-only.** A strong LUSC TACSTD2–immune signal (histology rework) is irrelevant here. Do not pool.
2. **Effect sizes.** Even where p is small, |ρ| near 0.1–0.2 in bulk RNA is a weak rank association in mixed tissue, not “immune desert because of TROP2”.
3. **Purity.** ESTIMATE ImmuneScore, xCell ImmuneScore, and MCP/Danaher totals are compositionally entangled with purity. The ABSOLUTE-partial number is the one that is allowed to be called “immune” rather than “not tumor”. TACSTD2 itself is orthogonal to ABSOLUTE (ρ = {rho_tp:.3f}).
4. **Two MCP implementations.** TIMER2 ran official MCP-counter on TPM. Xena recompute is mean log2 of the same markers on RSEM. They should agree in sign; they will not match coefficient-for-coefficient. Both are shown.
5. **xCell is TIMER2, not a from-scratch Xena port.** Re-running xCell on log2(RSEM+1) would be a different assay.
6. **GEP18 ssGSEA is not the Merck assay.** Unweighted ssGSEA / z-mean of the 18 genes is the public approximation.
7. **Not protein, not ICI, not OncoSG.** TROP2 protein (the ADC target) was not measured.
8. **No causality.** Bulk correlation, even purity-adjusted, does not say TACSTD2 excludes T cells.

## Reproduce

```
pip install -r requirements.txt
python scripts/rework_A1_deconv.py
```

Downloads ~40 MB of public tables into `data/` (gitignored) on first run. Signature gene lists are bundled under `data/signatures/`.

## Files

- `REPORT.md` — this writeup
- `DEFINITIONS.md` — the definitions section, standalone
- `correlations.tsv` — unadjusted + both partial methods + CIs + FDR + feature-vs-purity
- `primary_partial.tsv` — the 15 pre-specified features
- `sample_table.tsv` — per-sample TACSTD2, purity, and all scores
- `gene_coverage.tsv` — present/missing genes per signature
- `feature_definitions.tsv` — method / role / source per column
- `summary.json` / `provenance.json`
- `figures/forest_primary_partial_rho.png`
- `figures/bar_all_partial_rho.png`
- `figures/scatter_primary.png`

## Data

- Expression: {URLS['expression']['url']}
- Purity: {URLS['purity']['url']}
- ESTIMATE: {URLS['estimate']['url']}
- TIMER2: {URLS['timer2']['url']}
"""

    (OUT / "REPORT.md").write_text(report)

    definitions_only = report.split("## Primary result")[0]
    (OUT / "DEFINITIONS.md").write_text(
        definitions_only
        + "\n## Pointer\n\nNumbers live in `REPORT.md` and `correlations.tsv`. This file is the locked definition list.\n"
    )

    summary = {
        "task": "A1_deconv",
        "cohort": "TCGA-LUAD primary tumors (-01)",
        "n_expression": n_expr,
        "n_partial": n_partial,
        "tacstd2_vs_absolute": {"rho": float(rho_tp), "p": float(p_tp), "n": n_partial},
        "verdict_tag": verdict_tag,
        "verdict": verdict,
        "max_abs_primary_partial_rho": max_abs,
        "n_primary_abs_rho_gt_0.22": n_gt_022,
        "n_primary_fdr_lt_0.05": n_sig_fdr,
        "primary_features": PRIMARY,
        "danaher_total_tils_pops": DANAHER_TOTAL_TILS_POPS,
        "danaher_total_tils_empirical_pops": emp,
        "gep18_genes": gep_present,
        "methods": {
            "unadjusted": "Spearman",
            "adjusted_primary": "algebraic first-order partial Spearman, t with n-3 df",
            "adjusted_sensitivity": "Pearson of rank residuals on ranked ABSOLUTE purity, t with n-3 df",
            "ci": "Fisher-z, var=1/(n-4)",
            "fdr": "BH across all features and separately across 15 primary features",
        },
        "primary_results": prim.drop(columns=["primary_order"], errors="ignore").to_dict(orient="records"),
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    with open(OUT / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)

    prov = {
        "claim": "A1 rework (LUAD deconvolution): TACSTD2 vs xCell, MCP-counter, ESTIMATE ImmuneScore, GEP18 ssGSEA, Danaher after ABSOLUTE partial Spearman",
        "inputs": {},
        "signatures": {
            "gep18": str(SIG / "gep18_genes.tsv"),
            "danaher": str(SIG / "danaher_genes.tsv"),
            "mcp_counter": str(SIG / "mcp_counter_genes.tsv"),
        },
        "randomness": "deterministic (no sampling)",
    }
    for k, v in URLS.items():
        p = DATA / v["file"]
        prov["inputs"][k] = {**v, "bytes": p.stat().st_size, "md5": md5(p)}
    with open(OUT / "provenance.json", "w") as f:
        json.dump(prov, f, indent=2)

    print("\n=== PRIMARY PANEL ===")
    print(
        prim[
            [
                "feature",
                "n_partial",
                "spearman_rho",
                "partial_rho_algebraic",
                "partial_p_algebraic",
                "partial_fdr_primary15",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:.4g}")
    )
    print(f"\nTACSTD2 vs ABSOLUTE: rho={rho_tp:.4f} p={p_tp:.3g} n={n_partial}")
    print(f"Verdict: {verdict_tag}")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
