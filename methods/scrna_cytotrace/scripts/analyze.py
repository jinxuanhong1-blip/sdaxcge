#!/usr/bin/env python3
"""CytoTRACE-like stemness + Tirosh cycling on malignant cells.

Public GSE207422 and GSE241934 only. Sample is the inferential unit.
Cell-level Spearman is stored and labeled exploratory (pseudoreplication).

CytoTRACE2 R package is not run. Stemness = residual n_genes after OLS on
log1p(total UMI), rank-scaled to [0, 1] within malignant cells of each dataset
(Gulati et al. Science 2020 gene-count idea).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from pathlib import Path as _P
import sys as _sys
_sys.path.insert(0, str(_P(__file__).resolve().parent))
from gene_sets import (
    ALVEOLAR_DIFF,
    G2M_GENES,
    KERATIN_ALL,
    KERATIN_BASAL,
    KERATIN_SIMPLE,
    LINEAGE,
    NORMAL_LUNG,
    S_GENES,
    STEM_MARKERS,
)

DATA = Path("/tmp/scrna_cytotrace")
HERE = Path(__file__).resolve().parent.parent
RES = HERE / "results"
FIG = RES / "figures"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

MIN_MAL_PRIMARY = 20
MIN_MAL_SENS = 10
MIN_STRATUM = 8  # min malignant cells in cycling/non-cycling (or TACSTD2 high/low) per sample


def present(genes, columns) -> list[str]:
    return [g for g in genes if g in columns]


def log1p_cp10k(counts: pd.DataFrame, total: pd.Series) -> pd.DataFrame:
    tot = np.asarray(total, float)
    tot[tot <= 0] = np.nan
    return pd.DataFrame(
        np.log1p(np.asarray(counts, float) / tot[:, None] * 1e4),
        index=counts.index,
        columns=counts.columns,
    )


def module_mean(ln: pd.DataFrame, genes: list[str]) -> pd.Series:
    cols = present(genes, ln.columns)
    if not cols:
        return pd.Series(np.nan, index=ln.index)
    return ln[cols].mean(axis=1)


def cytotrace_like(n_genes: np.ndarray, total_umi: np.ndarray) -> np.ndarray:
    """Residual n_genes ~ log1p(UMI), rank-scaled to [0, 1]. Higher = more stem-like."""
    y = np.asarray(n_genes, float)
    x = np.log1p(np.asarray(total_umi, float))
    m = np.isfinite(x) & np.isfinite(y)
    out = np.full(y.shape, np.nan, dtype=float)
    if m.sum() < 10 or np.nanstd(x[m]) < 1e-8:
        if m.sum() >= 2:
            out[m] = stats.rankdata(y[m]) / m.sum()
        return out
    xm, ym = x[m], y[m]
    b = np.cov(xm, ym, ddof=1)[0, 1] / np.var(xm, ddof=1)
    a = ym.mean() - b * xm.mean()
    resid = ym - (a + b * xm)
    out[m] = stats.rankdata(resid) / m.sum()
    return out


def mwu(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, len(a), len(b)
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(u), float(p), len(a), len(b)


def spear(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def wilcoxon_paired(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if len(a) < 3:
        return np.nan, np.nan, int(len(a))
    # Wilcoxon signed-rank; if all diffs 0, p is undefined
    d = a - b
    if np.allclose(d, 0):
        return 0.0, 1.0, int(len(a))
    try:
        w, p = stats.wilcoxon(a, b, alternative="two-sided", zero_method="wilcox")
        return float(w), float(p), int(len(a))
    except ValueError:
        return np.nan, np.nan, int(len(a))


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"


def fmt_r(r) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def add_stat(rows, dataset, analysis, comparison, metric, n, stat_name, value, p, note=""):
    rows.append(dict(
        dataset=dataset, analysis=analysis, comparison=comparison, metric=metric,
        n=n, stat=stat_name, value=value, p_value=p, note=note,
    ))


def score_malignant(cells: pd.DataFrame, count_genes: list[str]) -> pd.DataFrame:
    ln = log1p_cp10k(cells[count_genes], cells["total_umi"])
    out = cells.copy()
    for gene in ("TACSTD2", "CLDN4"):
        if gene in ln.columns:
            out[f"{gene}_log1p"] = ln[gene]
        else:
            out[f"{gene}_log1p"] = np.nan
    out["score_S"] = module_mean(ln, S_GENES)
    out["score_G2M"] = module_mean(ln, G2M_GENES)
    out["cycle_score"] = (out["score_S"] + out["score_G2M"]) / 2.0
    out["keratin"] = module_mean(ln, KERATIN_ALL)
    out["keratin_simple"] = module_mean(ln, KERATIN_SIMPLE)
    out["keratin_basal"] = module_mean(ln, KERATIN_BASAL)
    out["alveolar_diff"] = module_mean(ln, ALVEOLAR_DIFF)
    out["stem_markers"] = module_mean(ln, STEM_MARKERS)
    return out


def assign_cycle_strata(mal: pd.DataFrame) -> pd.DataFrame:
    mal = mal.copy()
    q25 = float(mal["cycle_score"].quantile(0.25))
    q75 = float(mal["cycle_score"].quantile(0.75))
    mal["cycling"] = mal["cycle_score"] >= q75
    mal["noncycling"] = mal["cycle_score"] <= q25
    mal["cycle_q25"] = q25
    mal["cycle_q75"] = q75
    return mal


def per_sample_table(mal: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    rows = []
    for sid, g in mal.groupby(sample_col, dropna=False):
        n = len(g)
        cyc = g.loc[g["cycling"]]
        nc = g.loc[g["noncycling"]]
        tac = g["TACSTD2_log1p"]
        hi = g.loc[tac >= tac.quantile(0.75)] if n >= MIN_STRATUM else g.iloc[0:0]
        lo = g.loc[tac <= tac.quantile(0.25)] if n >= MIN_STRATUM else g.iloc[0:0]
        rec = {
            sample_col: sid,
            "n_malignant": n,
            "mean_TACSTD2": float(g["TACSTD2_log1p"].mean()),
            "mean_CLDN4": float(g["CLDN4_log1p"].mean()),
            "mean_cytotrace_like": float(g["cytotrace_like"].mean()),
            "mean_cycle": float(g["cycle_score"].mean()),
            "mean_keratin": float(g["keratin"].mean()),
            "mean_keratin_basal": float(g["keratin_basal"].mean()),
            "mean_alveolar": float(g["alveolar_diff"].mean()),
            "mean_stem_markers": float(g["stem_markers"].mean()),
            "frac_cycling": float(g["cycling"].mean()),
            "n_cycling": int(g["cycling"].sum()),
            "n_noncycling": int(g["noncycling"].sum()),
            "TACSTD2_cycling": float(cyc["TACSTD2_log1p"].mean()) if len(cyc) >= MIN_STRATUM else np.nan,
            "TACSTD2_noncycling": float(nc["TACSTD2_log1p"].mean()) if len(nc) >= MIN_STRATUM else np.nan,
            "CLDN4_cycling": float(cyc["CLDN4_log1p"].mean()) if len(cyc) >= MIN_STRATUM else np.nan,
            "CLDN4_noncycling": float(nc["CLDN4_log1p"].mean()) if len(nc) >= MIN_STRATUM else np.nan,
            "cyto_cycling": float(cyc["cytotrace_like"].mean()) if len(cyc) >= MIN_STRATUM else np.nan,
            "cyto_noncycling": float(nc["cytotrace_like"].mean()) if len(nc) >= MIN_STRATUM else np.nan,
            "keratin_TACSTD2_high": float(hi["keratin"].mean()) if len(hi) >= MIN_STRATUM else np.nan,
            "keratin_TACSTD2_low": float(lo["keratin"].mean()) if len(lo) >= MIN_STRATUM else np.nan,
            "cyto_TACSTD2_high": float(hi["cytotrace_like"].mean()) if len(hi) >= MIN_STRATUM else np.nan,
            "cyto_TACSTD2_low": float(lo["cytotrace_like"].mean()) if len(lo) >= MIN_STRATUM else np.nan,
            "alveolar_TACSTD2_high": float(hi["alveolar_diff"].mean()) if len(hi) >= MIN_STRATUM else np.nan,
            "alveolar_TACSTD2_low": float(lo["alveolar_diff"].mean()) if len(lo) >= MIN_STRATUM else np.nan,
        }
        if "mpr_group" in g.columns:
            rec["mpr_group"] = g["mpr_group"].iloc[0]
        if "cohort" in g.columns:
            rec["cohort"] = g["cohort"].iloc[0]
        if "timing" in g.columns:
            rec["timing"] = g["timing"].iloc[0]
        if "Pathology" in g.columns:
            rec["Pathology"] = g["Pathology"].iloc[0]
        if "Histology" in g.columns:
            rec["Histology"] = g["Histology"].iloc[0]
        rows.append(rec)
    return pd.DataFrame(rows)


def sample_tests(rows, dataset, samp: pd.DataFrame, min_mal: int, note: str):
    ok = samp.loc[samp["n_malignant"] >= min_mal].copy()
    pairs = [
        ("mean_TACSTD2", "mean_cytotrace_like", "TACSTD2 vs CytoTRACE-like"),
        ("mean_CLDN4", "mean_cytotrace_like", "CLDN4 vs CytoTRACE-like"),
        ("mean_TACSTD2", "mean_cycle", "TACSTD2 vs cycling score"),
        ("mean_CLDN4", "mean_cycle", "CLDN4 vs cycling score"),
        ("mean_TACSTD2", "mean_keratin", "TACSTD2 vs keratin"),
        ("mean_CLDN4", "mean_keratin", "CLDN4 vs keratin"),
        ("mean_TACSTD2", "mean_keratin_basal", "TACSTD2 vs basal keratin"),
        ("mean_TACSTD2", "mean_alveolar", "TACSTD2 vs alveolar/diff"),
        ("mean_TACSTD2", "frac_cycling", "TACSTD2 vs fraction cycling"),
        ("mean_CLDN4", "frac_cycling", "CLDN4 vs fraction cycling"),
        ("mean_cytotrace_like", "mean_cycle", "CytoTRACE-like vs cycling score"),
        ("mean_cytotrace_like", "mean_keratin", "CytoTRACE-like vs keratin"),
    ]
    for x, y, label in pairs:
        r, p, n = spear(ok[x], ok[y])
        add_stat(rows, dataset, "sample_spearman", label, f"{x}~{y}", n,
                 "spearman_rho", r, p, note)
        print(f"  {label}: n={n} ρ={fmt_r(r)} p={fmt_p(p)}", flush=True)

    # combinatorial: TACSTD2 cycling vs non-cycling (paired samples)
    w, p, n = wilcoxon_paired(ok["TACSTD2_cycling"], ok["TACSTD2_noncycling"])
    add_stat(rows, dataset, "combinatorial_paired",
             "TACSTD2 cycling vs non-cycling (paired sample means)",
             "TACSTD2_log1p", n, "wilcoxon_w", w, p,
             note + "; cycling=top quartile cycle_score; non-cycling=bottom quartile")
    print(f"  TACSTD2 cycling vs non-cycling paired n={n} W={w} p={fmt_p(p)}", flush=True)

    w, p, n = wilcoxon_paired(ok["CLDN4_cycling"], ok["CLDN4_noncycling"])
    add_stat(rows, dataset, "combinatorial_paired",
             "CLDN4 cycling vs non-cycling (paired sample means)",
             "CLDN4_log1p", n, "wilcoxon_w", w, p, note)

    # extra: TACSTD2-high more keratin / less stem-like?
    w, p, n = wilcoxon_paired(ok["keratin_TACSTD2_high"], ok["keratin_TACSTD2_low"])
    add_stat(rows, dataset, "tacstd2_high_diff",
             "keratin TACSTD2-high vs low (paired sample means)",
             "keratin", n, "wilcoxon_w", w, p,
             note + "; high/low = within-sample TACSTD2 quartiles")
    print(f"  keratin TACSTD2-high vs low paired n={n} W={w} p={fmt_p(p)}", flush=True)

    w, p, n = wilcoxon_paired(ok["cyto_TACSTD2_high"], ok["cyto_TACSTD2_low"])
    add_stat(rows, dataset, "tacstd2_high_diff",
             "CytoTRACE-like TACSTD2-high vs low (paired sample means)",
             "cytotrace_like", n, "wilcoxon_w", w, p,
             note + "; high/low = within-sample TACSTD2 quartiles")

    w, p, n = wilcoxon_paired(ok["alveolar_TACSTD2_high"], ok["alveolar_TACSTD2_low"])
    add_stat(rows, dataset, "tacstd2_high_diff",
             "alveolar/diff TACSTD2-high vs low (paired sample means)",
             "alveolar_diff", n, "wilcoxon_w", w, p, note)

    if "mpr_group" not in ok.columns:
        return ok
    mpr = ok.loc[ok.mpr_group == "MPR"]
    nmpr = ok.loc[ok.mpr_group == "NMPR"]
    for metric, label in (
        ("mean_TACSTD2", "malignant TACSTD2"),
        ("mean_CLDN4", "malignant CLDN4"),
        ("mean_cytotrace_like", "malignant CytoTRACE-like"),
        ("mean_cycle", "malignant cycling score"),
        ("frac_cycling", "malignant fraction cycling"),
        ("mean_keratin", "malignant keratin"),
    ):
        u, p, na, nb = mwu(nmpr[metric], mpr[metric])
        add_stat(rows, dataset, "mpr_split", f"NMPR vs MPR {label}",
                 metric, f"{na}+{nb}", "mannwhitney_u", u, p,
                 f"{note}; MPR median={mpr[metric].median():.4f} (n={na}); "
                 f"NMPR median={nmpr[metric].median():.4f} (n={nb}); pCR counted as MPR")
        print(f"  NMPR vs MPR {label}: n={na}+{nb} U={u} p={fmt_p(p)} "
              f"med NMPR={nmpr[metric].median():.4f} MPR={mpr[metric].median():.4f}", flush=True)

    # Spearman within MPR / NMPR
    for grp, sub in (("MPR", mpr), ("NMPR", nmpr)):
        r, p, n = spear(sub["mean_TACSTD2"], sub["mean_cytotrace_like"])
        add_stat(rows, dataset, "sample_spearman_by_mpr",
                 f"TACSTD2 vs CytoTRACE-like within {grp}",
                 "mean_TACSTD2~mean_cytotrace_like", n, "spearman_rho", r, p, note)
        r, p, n = spear(sub["mean_TACSTD2"], sub["mean_cycle"])
        add_stat(rows, dataset, "sample_spearman_by_mpr",
                 f"TACSTD2 vs cycling score within {grp}",
                 "mean_TACSTD2~mean_cycle", n, "spearman_rho", r, p, note)
        r, p, n = spear(sub["mean_TACSTD2"], sub["mean_keratin"])
        add_stat(rows, dataset, "sample_spearman_by_mpr",
                 f"TACSTD2 vs keratin within {grp}",
                 "mean_TACSTD2~mean_keratin", n, "spearman_rho", r, p, note)
        r, p, n = spear(sub["mean_TACSTD2"], sub["frac_cycling"])
        add_stat(rows, dataset, "sample_spearman_by_mpr",
                 f"TACSTD2 vs fraction cycling within {grp}",
                 "mean_TACSTD2~frac_cycling", n, "spearman_rho", r, p, note)
    return ok


def cell_exploratory(rows, dataset, mal: pd.DataFrame, note: str):
    pairs = [
        ("TACSTD2_log1p", "cytotrace_like", "TACSTD2 vs CytoTRACE-like"),
        ("CLDN4_log1p", "cytotrace_like", "CLDN4 vs CytoTRACE-like"),
        ("TACSTD2_log1p", "cycle_score", "TACSTD2 vs cycling score"),
        ("CLDN4_log1p", "cycle_score", "CLDN4 vs cycling score"),
        ("TACSTD2_log1p", "keratin", "TACSTD2 vs keratin"),
        ("CLDN4_log1p", "keratin", "CLDN4 vs keratin"),
        ("TACSTD2_log1p", "keratin_basal", "TACSTD2 vs basal keratin"),
        ("TACSTD2_log1p", "alveolar_diff", "TACSTD2 vs alveolar/diff"),
        ("TACSTD2_log1p", "CLDN4_log1p", "TACSTD2 vs CLDN4"),
    ]
    for x, y, label in pairs:
        r, p, n = spear(mal[x], mal[y])
        add_stat(rows, dataset, "cell_spearman_exploratory", label, f"{x}~{y}",
                 n, "spearman_rho", r, p, note + "; cells as replicates — exploratory")
        print(f"  [cell] {label}: n={n} ρ={fmt_r(r)} p={fmt_p(p)}", flush=True)

    cyc = mal.loc[mal["cycling"]]
    nc = mal.loc[mal["noncycling"]]
    r, p, n = spear(cyc["TACSTD2_log1p"], cyc["cytotrace_like"])
    add_stat(rows, dataset, "cell_spearman_exploratory",
             "TACSTD2 vs CytoTRACE-like within cycling cells",
             "TACSTD2~cytotrace_like", n, "spearman_rho", r, p,
             note + "; cycling=top quartile; exploratory")
    r, p, n = spear(nc["TACSTD2_log1p"], nc["cytotrace_like"])
    add_stat(rows, dataset, "cell_spearman_exploratory",
             "TACSTD2 vs CytoTRACE-like within non-cycling cells",
             "TACSTD2~cytotrace_like", n, "spearman_rho", r, p,
             note + "; non-cycling=bottom quartile; exploratory")
    u, p, na, nb = mwu(cyc["TACSTD2_log1p"], nc["TACSTD2_log1p"])
    add_stat(rows, dataset, "cell_mwu_exploratory",
             "TACSTD2 cycling vs non-cycling cells",
             "TACSTD2_log1p", f"{na}+{nb}", "mannwhitney_u", u, p,
             note + "; exploratory")
    u, p, na, nb = mwu(cyc["keratin"], nc["keratin"])
    add_stat(rows, dataset, "cell_mwu_exploratory",
             "keratin cycling vs non-cycling cells",
             "keratin", f"{na}+{nb}", "mannwhitney_u", u, p,
             note + "; exploratory")


def scatter_mpr(ax, samp, x, y, title):
    colors = {"MPR": "#2ca02c", "NMPR": "#d62728"}
    for grp, sub in samp.groupby("mpr_group"):
        ax.scatter(sub[x], sub[y], c=colors.get(grp, "#7f7f7f"),
                   label=f"{grp} n={len(sub)}", s=36, edgecolors="k", linewidths=0.3)
    r, p, n = spear(samp[x], samp[y])
    ax.set_title(f"{title}\nn={n} ρ={fmt_r(r)} p={fmt_p(p)}", fontsize=9)
    ax.legend(fontsize=7, loc="best")


def box_mpr(ax, samp, metric, title):
    mpr = samp.loc[samp.mpr_group == "MPR", metric].dropna()
    nmpr = samp.loc[samp.mpr_group == "NMPR", metric].dropna()
    ax.boxplot([mpr, nmpr], tick_labels=[f"MPR n={len(mpr)}", f"NMPR n={len(nmpr)}"])
    u, p, na, nb = mwu(nmpr, mpr)
    ax.set_title(f"{title}\nMWU p={fmt_p(p)}", fontsize=9)


# ---------------------------------------------------------------------------
# GSE207422
# ---------------------------------------------------------------------------
def run_gse207422(rows) -> pd.DataFrame:
    print("=== GSE207422 ===", flush=True)
    cells = pd.read_csv(DATA / "gse207422_panel_cells.tsv.gz", sep="\t")
    cells["sample"] = cells["barcode"].str.rsplit("_", n=1).str[0]
    cells = cells.loc[cells.total_umi >= 200].copy()
    count_genes = [c for c in cells.columns if c not in ("barcode", "total_umi", "n_genes", "sample")]
    ln = log1p_cp10k(cells[count_genes], cells["total_umi"])
    scores = pd.DataFrame(
        {lin: module_mean(ln, gs) for lin, gs in LINEAGE.items()},
        index=cells.index,
    )
    cells["lineage"] = scores.idxmax(axis=1)
    cells["normal_lung"] = module_mean(ln, NORMAL_LUNG)
    epi = cells["lineage"] == "Epithelial"
    nl_cut = float(cells.loc[epi, "normal_lung"].quantile(0.75)) if epi.any() else 0.0
    cells["malignant_like"] = epi & (cells["normal_lung"] <= nl_cut)
    print(f"  cells={len(cells)} epi={int(epi.sum())} mal={int(cells.malignant_like.sum())} nl_cut={nl_cut:.4f}",
          flush=True)

    scored = score_malignant(cells, count_genes)
    mal = scored.loc[scored.malignant_like].copy()
    mal["cytotrace_like"] = cytotrace_like(mal["n_genes"].to_numpy(), mal["total_umi"].to_numpy())
    mal = assign_cycle_strata(mal)

    meta = pd.read_excel(DATA / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")].copy()
    meta = meta.rename(columns={"Sample": "sample"})
    meta["mpr_group"] = meta["Pathologic Response"].map({"MPR": "MPR", "pCR": "MPR", "NMPR": "NMPR"})
    meta["timing"] = meta["Resource"].map({
        "Pre-treatment biopsy": "pre",
        "Post-treatment surgery": "post",
    })
    mal = mal.merge(meta[["sample", "mpr_group", "timing", "Pathology", "Pathologic Response", "RECIST"]],
                    on="sample", how="left")

    cell_exploratory(rows, "GSE207422", mal,
                     "malignant-like = epithelial AND normal-lung ≤ epi p75; UMI≥200")

    samp = per_sample_table(mal, "sample")
    extra_cols = [c for c in ("mpr_group", "timing", "Pathology", "Patient") if c not in samp.columns]
    if extra_cols:
        samp = samp.merge(meta[["sample"] + extra_cols], on="sample", how="left")
    samp.to_csv(RES / "gse207422_per_sample.tsv", sep="\t", index=False)

    post = samp.loc[samp.timing == "post"].copy()
    print(f"  post samples={len(post)}", flush=True)
    sample_tests(rows, "GSE207422", post, MIN_MAL_PRIMARY,
                 "GSE207422 post-tx; malignant-like; min 20 cells; unit=sample")
    sample_tests(rows, "GSE207422_min10", post, MIN_MAL_SENS,
                 "GSE207422 post-tx; malignant-like; min 10 cells; sensitivity")

    # all-epithelial sensitivity on post-tx
    epi_cells = scored.loc[scored.lineage == "Epithelial"].copy()
    epi_cells["cytotrace_like"] = cytotrace_like(epi_cells["n_genes"].to_numpy(),
                                                 epi_cells["total_umi"].to_numpy())
    epi_cells = assign_cycle_strata(epi_cells)
    epi_cells = epi_cells.merge(meta[["sample", "mpr_group", "timing"]], on="sample", how="left")
    epi_samp = per_sample_table(epi_cells, "sample")
    extra_cols = [c for c in ("mpr_group", "timing") if c not in epi_samp.columns]
    if extra_cols:
        epi_samp = epi_samp.merge(meta[["sample"] + extra_cols], on="sample", how="left")
    epi_post = epi_samp.loc[epi_samp.timing == "post"]
    sample_tests(rows, "GSE207422_allEpi", epi_post, MIN_MAL_SENS,
                 "GSE207422 post-tx; ALL epithelial (no normal-lung gate); min 10; sensitivity")

    ok = post.loc[post.n_malignant >= MIN_MAL_SENS]
    fig, axes = plt.subplots(2, 3, figsize=(11, 7))
    scatter_mpr(axes[0, 0], ok, "mean_TACSTD2", "mean_cytotrace_like", "TACSTD2 vs CytoTRACE-like")
    scatter_mpr(axes[0, 1], ok, "mean_TACSTD2", "mean_cycle", "TACSTD2 vs cycling")
    scatter_mpr(axes[0, 2], ok, "mean_TACSTD2", "mean_keratin", "TACSTD2 vs keratin")
    box_mpr(axes[1, 0], ok, "mean_cytotrace_like", "CytoTRACE-like by MPR")
    box_mpr(axes[1, 1], ok, "frac_cycling", "Fraction cycling by MPR")
    box_mpr(axes[1, 2], ok, "mean_keratin", "Keratin by MPR")
    for ax in axes[0]:
        ax.set_xlabel("malignant TACSTD2 log1p CP10k")
    axes[0, 0].set_ylabel("CytoTRACE-like")
    axes[0, 1].set_ylabel("cycling score")
    axes[0, 2].set_ylabel("keratin")
    fig.suptitle("GSE207422 post-tx malignant-like (sample unit)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "gse207422_stemness_mpr.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    paired = ok.dropna(subset=["TACSTD2_cycling", "TACSTD2_noncycling"])
    for _, r in paired.iterrows():
        col = "#2ca02c" if r.get("mpr_group") == "MPR" else "#d62728"
        ax.plot([0, 1], [r["TACSTD2_noncycling"], r["TACSTD2_cycling"]],
                color=col, alpha=0.55, lw=1)
        ax.scatter([0, 1], [r["TACSTD2_noncycling"], r["TACSTD2_cycling"]],
                   color=col, s=22)
    w, p, n = wilcoxon_paired(paired["TACSTD2_cycling"], paired["TACSTD2_noncycling"])
    ax.set_xticks([0, 1], ["non-cycling\n(bottom Q)", "cycling\n(top Q)"])
    ax.set_ylabel("malignant TACSTD2 log1p CP10k")
    ax.set_title(f"GSE207422 combinatorial TACSTD2\npaired n={n} Wilcoxon p={fmt_p(p)}")
    fig.tight_layout()
    fig.savefig(FIG / "gse207422_combinatorial_cycling.png", dpi=140)
    plt.close(fig)
    return samp


# ---------------------------------------------------------------------------
# GSE241934
# ---------------------------------------------------------------------------
def mpr_of(x):
    x = str(x)
    if x in {"non-MPR", "NMPR", "nonMPR"}:
        return "NMPR"
    if x in {"MPR", "pCR", "MPR (pCR)"}:
        return "MPR"
    return np.nan


def load_gse241934_cohort(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    major_col = "major.cell.type" if "major.cell.type" in df.columns else "major_cell_type"
    df["is_epi"] = df[major_col].astype(str).eq("Epi")
    resp_col = "Pathological Response" if "Pathological Response" in df.columns else None
    if resp_col:
        df["mpr_group"] = df[resp_col].map(mpr_of)
    else:
        df["mpr_group"] = np.nan
    if "sampleID" not in df.columns:
        raise KeyError("sampleID missing")
    return df


def run_gse241934(rows) -> pd.DataFrame:
    print("=== GSE241934 ===", flush=True)
    iit = load_gse241934_cohort(DATA / "gse241934_iit_panel_cells.tsv.gz")
    real = load_gse241934_cohort(DATA / "gse241934_real_panel_cells.tsv.gz")
    cells = pd.concat([iit, real], ignore_index=True)
    cells = cells.loc[cells.total_umi >= 200].copy()
    print(f"  cells={len(cells)} epi={int(cells.is_epi.sum())} cohorts={cells.cohort.value_counts().to_dict()}",
          flush=True)
    print(f"  meta cols sample={list(cells.columns)[:20]}", flush=True)

    skip = {
        "barcode", "total_umi", "n_genes", "cohort", "sampleID", "cellID",
        "is_epi", "mpr_group",
    }
    count_genes = [c for c in cells.columns if c in set(
        ["TACSTD2", "CLDN4"] + S_GENES + G2M_GENES + KERATIN_ALL + ALVEOLAR_DIFF + STEM_MARKERS
        + NORMAL_LUNG + [g for gs in LINEAGE.values() for g in gs]
    ) and c not in skip]
    # keep only numeric gene columns that exist
    count_genes = [c for c in count_genes if c in cells.columns and pd.api.types.is_numeric_dtype(cells[c])]
    scored = score_malignant(cells, count_genes)
    mal = scored.loc[scored.is_epi].copy()
    # CytoTRACE-like within each cohort (depth / capture differs)
    mal["cytotrace_like"] = np.nan
    for coh, idx in mal.groupby("cohort").groups.items():
        sub = mal.loc[idx]
        mal.loc[idx, "cytotrace_like"] = cytotrace_like(
            sub["n_genes"].to_numpy(), sub["total_umi"].to_numpy()
        )
    # cycle strata within each cohort
    mal["cycling"] = False
    mal["noncycling"] = False
    for coh, idx in mal.groupby("cohort").groups.items():
        sub = assign_cycle_strata(mal.loc[idx])
        mal.loc[idx, "cycling"] = sub["cycling"].to_numpy()
        mal.loc[idx, "noncycling"] = sub["noncycling"].to_numpy()

    cell_exploratory(rows, "GSE241934", mal,
                     "author major.cell.type==Epi; UMI≥200; CytoTRACE-like within cohort")

    samp = per_sample_table(mal, "sampleID")
    extra = mal.groupby("sampleID").agg(
        mpr_group=("mpr_group", "first"),
        cohort=("cohort", "first"),
    ).reset_index()
    need = [c for c in ("mpr_group", "cohort") if c not in samp.columns]
    if need:
        samp = samp.merge(extra[["sampleID"] + need], on="sampleID", how="left")
    samp.to_csv(RES / "gse241934_per_sample.tsv", sep="\t", index=False)
    print(samp[["sampleID", "cohort", "mpr_group", "n_malignant", "mean_TACSTD2",
                "mean_cytotrace_like", "frac_cycling"]].to_string(index=False), flush=True)

    print("  -- combined IIT+REAL --", flush=True)
    sample_tests(rows, "GSE241934", samp, MIN_MAL_PRIMARY,
                 "GSE241934 IIT+REAL; author Epi; min 20 cells; unit=sample")
    sample_tests(rows, "GSE241934_min10", samp, MIN_MAL_SENS,
                 "GSE241934 IIT+REAL; author Epi; min 10 cells; sensitivity")

    for coh, label in (("IIT_EGFRmut", "GSE241934_IIT"), ("REAL_WT", "GSE241934_REAL")):
        print(f"  -- {coh} --", flush=True)
        sample_tests(rows, label, samp.loc[samp.cohort == coh], MIN_MAL_SENS,
                     f"GSE241934 {coh} only; author Epi; min 10; unit=sample")

    ok = samp.loc[samp.n_malignant >= MIN_MAL_SENS]
    fig, axes = plt.subplots(2, 3, figsize=(11, 7))
    scatter_mpr(axes[0, 0], ok, "mean_TACSTD2", "mean_cytotrace_like", "TACSTD2 vs CytoTRACE-like")
    scatter_mpr(axes[0, 1], ok, "mean_TACSTD2", "mean_cycle", "TACSTD2 vs cycling")
    scatter_mpr(axes[0, 2], ok, "mean_TACSTD2", "mean_keratin", "TACSTD2 vs keratin")
    box_mpr(axes[1, 0], ok, "mean_cytotrace_like", "CytoTRACE-like by MPR")
    box_mpr(axes[1, 1], ok, "frac_cycling", "Fraction cycling by MPR")
    box_mpr(axes[1, 2], ok, "mean_keratin", "Keratin by MPR")
    for ax in axes[0]:
        ax.set_xlabel("epithelial TACSTD2 log1p CP10k")
    axes[0, 0].set_ylabel("CytoTRACE-like")
    axes[0, 1].set_ylabel("cycling score")
    axes[0, 2].set_ylabel("keratin")
    fig.suptitle("GSE241934 author Epi (sample unit; IIT+REAL)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "gse241934_stemness_mpr.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    paired = ok.dropna(subset=["TACSTD2_cycling", "TACSTD2_noncycling"])
    for _, r in paired.iterrows():
        col = "#2ca02c" if r.get("mpr_group") == "MPR" else "#d62728"
        ax.plot([0, 1], [r["TACSTD2_noncycling"], r["TACSTD2_cycling"]],
                color=col, alpha=0.45, lw=1)
        ax.scatter([0, 1], [r["TACSTD2_noncycling"], r["TACSTD2_cycling"]],
                   color=col, s=18)
    w, p, n = wilcoxon_paired(paired["TACSTD2_cycling"], paired["TACSTD2_noncycling"])
    ax.set_xticks([0, 1], ["non-cycling\n(bottom Q)", "cycling\n(top Q)"])
    ax.set_ylabel("epithelial TACSTD2 log1p CP10k")
    ax.set_title(f"GSE241934 combinatorial TACSTD2\npaired n={n} Wilcoxon p={fmt_p(p)}")
    fig.tight_layout()
    fig.savefig(FIG / "gse241934_combinatorial_cycling.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    paired = ok.dropna(subset=["keratin_TACSTD2_high", "keratin_TACSTD2_low"])
    for _, r in paired.iterrows():
        col = "#2ca02c" if r.get("mpr_group") == "MPR" else "#d62728"
        ax.plot([0, 1], [r["keratin_TACSTD2_low"], r["keratin_TACSTD2_high"]],
                color=col, alpha=0.45, lw=1)
        ax.scatter([0, 1], [r["keratin_TACSTD2_low"], r["keratin_TACSTD2_high"]],
                   color=col, s=18)
    w, p, n = wilcoxon_paired(paired["keratin_TACSTD2_high"], paired["keratin_TACSTD2_low"])
    ax.set_xticks([0, 1], ["TACSTD2-low\n(bottom Q)", "TACSTD2-high\n(top Q)"])
    ax.set_ylabel("keratin score")
    ax.set_title(f"GSE241934 TACSTD2-high vs keratin\npaired n={n} Wilcoxon p={fmt_p(p)}")
    fig.tight_layout()
    fig.savefig(FIG / "gse241934_tacstd2_keratin.png", dpi=140)
    plt.close(fig)
    return samp


def main() -> None:
    rows: list[dict] = []
    s207 = run_gse207422(rows)
    s241 = run_gse241934(rows)
    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(RES / "stats.tsv", sep="\t", index=False)

    def pick(ds, analysis, comparison):
        hit = stats_df[(stats_df.dataset == ds) & (stats_df.analysis == analysis)
                       & (stats_df.comparison == comparison)]
        if hit.empty:
            return None
        r = hit.iloc[0]
        return {"n": r["n"], "stat": r["stat"], "value": r["value"], "p": r["p_value"]}

    summary = {
        "note": "CytoTRACE2 R package was not run. Stemness is residual n_genes "
                "(Gulati 2020 idea) rank-scaled within malignant cells. "
                "Primary unit is the sample. Cell-level rows are exploratory.",
        "GSE207422_post_min20": {
            "TACSTD2_vs_cyto": pick("GSE207422", "sample_spearman", "TACSTD2 vs CytoTRACE-like"),
            "TACSTD2_vs_cycle": pick("GSE207422", "sample_spearman", "TACSTD2 vs cycling score"),
            "TACSTD2_vs_keratin": pick("GSE207422", "sample_spearman", "TACSTD2 vs keratin"),
            "NMPR_vs_MPR_cyto": pick("GSE207422", "mpr_split", "NMPR vs MPR malignant CytoTRACE-like"),
            "paired_TACSTD2_cycling": pick("GSE207422", "combinatorial_paired",
                                           "TACSTD2 cycling vs non-cycling (paired sample means)"),
        },
        "GSE241934_min20": {
            "TACSTD2_vs_cyto": pick("GSE241934", "sample_spearman", "TACSTD2 vs CytoTRACE-like"),
            "TACSTD2_vs_cycle": pick("GSE241934", "sample_spearman", "TACSTD2 vs cycling score"),
            "TACSTD2_vs_keratin": pick("GSE241934", "sample_spearman", "TACSTD2 vs keratin"),
            "NMPR_vs_MPR_cyto": pick("GSE241934", "mpr_split", "NMPR vs MPR malignant CytoTRACE-like"),
            "paired_TACSTD2_cycling": pick("GSE241934", "combinatorial_paired",
                                           "TACSTD2 cycling vs non-cycling (paired sample means)"),
            "paired_keratin_TACSTD2hi": pick("GSE241934", "tacstd2_high_diff",
                                             "keratin TACSTD2-high vs low (paired sample means)"),
        },
        "n_samples_GSE207422_post": int((s207.timing == "post").sum()) if "timing" in s207.columns else None,
        "n_samples_GSE241934": int(len(s241)),
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    print(f"wrote {RES / 'stats.tsv'} rows={len(stats_df)}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
