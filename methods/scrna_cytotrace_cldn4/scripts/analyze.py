#!/usr/bin/env python3
"""CLDN4 vs CytoTRACE-like potency in GSE131907 epithelial cells.

CytoTRACE2 is not run. Potency = residual n_genes after OLS on log1p(UMI),
rank-scaled to [0, 1] within the epithelial universe of that contrast
(Gulati et al. Science 2020 gene-count idea).

Primary unit = sample. Cell-level Spearman is exploratory (pseudoreplication).
CLDN4 only — TACSTD2 is not tested.
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

DATA = Path("/tmp/gse131907")
HERE = Path(__file__).resolve().parent.parent
RES = HERE / "results"
FIG = RES / "figures"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

MIN_EPI = 20
MIN_STRATUM = 8
MIN_UMI = 200
TUMOR_ORIGINS = ("tLung", "tL/B", "mLN", "mBrain", "PE")


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


def add_stat(rows, analysis, comparison, n, stat_name, value, p, note=""):
    rows.append(dict(
        dataset="GSE131907",
        analysis=analysis,
        comparison=comparison,
        n=n,
        stat=stat_name,
        value=value,
        p_value=p,
        note=note,
    ))


def score_universe(cells: pd.DataFrame) -> pd.DataFrame:
    out = cells.copy()
    tot = np.asarray(out["total_umi"], float)
    tot[tot <= 0] = np.nan
    cldn4 = np.asarray(out["CLDN4"], float)
    out["CLDN4_log1p_cp10k"] = np.log1p(cldn4 / tot * 1e4)
    out["cytotrace_like"] = cytotrace_like(out["n_genes"].to_numpy(), out["total_umi"].to_numpy())
    return out


def per_sample(epi: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sid, origin), g in epi.groupby(["Sample", "Sample_Origin"], dropna=False):
        n = len(g)
        cldn = g["CLDN4_log1p_cp10k"]
        hi = g.loc[cldn >= cldn.quantile(0.75)] if n >= MIN_STRATUM else g.iloc[0:0]
        lo = g.loc[cldn <= cldn.quantile(0.25)] if n >= MIN_STRATUM else g.iloc[0:0]
        rows.append({
            "Sample": sid,
            "Sample_Origin": origin,
            "n_epithelial": n,
            "mean_CLDN4": float(cldn.mean()),
            "median_CLDN4": float(cldn.median()),
            "pct_CLDN4_pos": float((g["CLDN4"] > 0).mean() * 100.0),
            "mean_cytotrace_like": float(g["cytotrace_like"].mean()),
            "median_cytotrace_like": float(g["cytotrace_like"].median()),
            "mean_n_genes": float(g["n_genes"].mean()),
            "mean_total_umi": float(g["total_umi"].mean()),
            "mean_EPCAM_log1p_cp10k": float(np.log1p(g["EPCAM"] / g["total_umi"] * 1e4).mean()) if "EPCAM" in g else np.nan,
            "n_CLDN4_high": int(len(hi)),
            "n_CLDN4_low": int(len(lo)),
            "mean_cyto_CLDN4_high": float(hi["cytotrace_like"].mean()) if len(hi) else np.nan,
            "mean_cyto_CLDN4_low": float(lo["cytotrace_like"].mean()) if len(lo) else np.nan,
            "eligible": int(n >= MIN_EPI),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ann = pd.read_csv(DATA / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t", dtype=str)
    expr = pd.read_csv(DATA / "gse131907_cldn4_cells.tsv.gz", sep="\t")
    cells = ann.merge(expr, on="Index", how="inner")
    if len(cells) != len(ann):
        raise SystemExit(f"annotation/matrix mismatch: ann={len(ann)} merged={len(cells)}")
    if "CLDN4" not in cells.columns:
        raise SystemExit("CLDN4 missing from streamed matrix")

    n_before_qc = len(cells)
    cells = cells[cells["total_umi"] >= MIN_UMI].copy()
    epi_all = cells[cells["Cell_type"] == "Epithelial cells"].copy()

    coverage = {
        "dataset": "GSE131907",
        "pmid": "32385277",
        "n_cells_matrix": int(n_before_qc),
        "n_cells_umi_ge_200": int(len(cells)),
        "n_epithelial_umi_ge_200": int(len(epi_all)),
        "n_samples_total": int(cells["Sample"].nunique()),
        "cytotrace2_run": False,
        "potency_definition": "residual n_genes ~ log1p(total_UMI), rank-scaled within the epithelial universe of each contrast",
        "test_gene": "CLDN4",
        "unit": "sample",
        "min_epithelial_per_sample": MIN_EPI,
        "min_umi": MIN_UMI,
        "ici_or_mpr_labels": False,
    }

    universes = {
        "tLung_epithelial": epi_all[epi_all["Sample_Origin"] == "tLung"],
        "tumor_epithelial": epi_all[epi_all["Sample_Origin"].isin(TUMOR_ORIGINS)],
        "nLung_epithelial": epi_all[epi_all["Sample_Origin"] == "nLung"],
        "all_epithelial": epi_all,
    }

    stat_rows = []
    sample_tables = {}
    cell_exploratory = {}

    for name, uni in universes.items():
        uni = score_universe(uni.copy())
        samp = per_sample(uni)
        sample_tables[name] = samp
        elig = samp[samp["eligible"] == 1]
        r, p, n = spear(elig["mean_CLDN4"], elig["mean_cytotrace_like"])
        add_stat(
            stat_rows, name, "sample-mean CLDN4 vs CytoTRACE-like",
            n, "spearman_rho", r, p,
            note=f"min {MIN_EPI} epithelial cells; potency ranked within {name}",
        )
        paired = elig[(elig["n_CLDN4_high"] >= MIN_STRATUM) & (elig["n_CLDN4_low"] >= MIN_STRATUM)]
        w, pw, nw = wilcoxon_paired(paired["mean_cyto_CLDN4_high"], paired["mean_cyto_CLDN4_low"])
        delta = float((paired["mean_cyto_CLDN4_high"] - paired["mean_cyto_CLDN4_low"]).median()) if len(paired) else np.nan
        add_stat(
            stat_rows, name, "paired CytoTRACE-like CLDN4-high vs CLDN4-low",
            nw, "wilcoxon_W", w, pw,
            note=f"within-sample quartiles; median Δ(high-low)={delta:.4f}" if np.isfinite(delta) else "too few paired samples",
        )
        cr, cp, cn = spear(uni["CLDN4_log1p_cp10k"], uni["cytotrace_like"])
        add_stat(
            stat_rows, name, "cell-level CLDN4 vs CytoTRACE-like (exploratory)",
            cn, "spearman_rho", cr, cp,
            note="pseudoreplication; not confirmatory",
        )
        cell_exploratory[name] = {
            "n_cells": int(len(uni)),
            "n_eligible_samples": int(n),
            "n_paired_quartile_samples": int(nw),
            "sample_spearman_rho": r,
            "sample_spearman_p": p,
            "paired_wilcoxon_W": w,
            "paired_wilcoxon_p": pw,
            "paired_median_delta_high_minus_low": delta,
            "cell_spearman_rho": cr,
            "cell_spearman_p": cp,
        }
        samp.to_csv(RES / f"per_sample_{name}.tsv", sep="\t", index=False)

    # Origin-stratified tumor-site sensitivity (not mixed into primary)
    tumor_scored = score_universe(universes["tumor_epithelial"].copy())
    tumor_samp = per_sample(tumor_scored)
    for origin, g in tumor_samp[tumor_samp["eligible"] == 1].groupby("Sample_Origin"):
        r, p, n = spear(g["mean_CLDN4"], g["mean_cytotrace_like"])
        add_stat(
            stat_rows, f"tumor_epithelial:{origin}",
            "sample-mean CLDN4 vs CytoTRACE-like",
            n, "spearman_rho", r, p,
            note="origin-stratified; potency ranked within all tumor epithelium",
        )

    stats_df = pd.DataFrame(stat_rows)
    stats_df.to_csv(RES / "stats.tsv", sep="\t", index=False)

    primary = next(s for s in stat_rows if s["analysis"] == "tLung_epithelial" and s["comparison"].startswith("sample-mean"))
    tumor = next(s for s in stat_rows if s["analysis"] == "tumor_epithelial" and s["comparison"].startswith("sample-mean"))
    nlung = next(s for s in stat_rows if s["analysis"] == "nLung_epithelial" and s["comparison"].startswith("sample-mean"))
    paired_t = next(s for s in stat_rows if s["analysis"] == "tLung_epithelial" and s["comparison"].startswith("paired"))
    paired_u = next(s for s in stat_rows if s["analysis"] == "tumor_epithelial" and s["comparison"].startswith("paired"))

    summary = {
        "coverage": coverage,
        "primary": {
            "universe": "tLung author epithelial cells",
            "n_samples": primary["n"],
            "spearman_rho": primary["value"],
            "spearman_p": primary["p_value"],
            "paired_n": paired_t["n"],
            "paired_p": paired_t["p_value"],
            "paired_median_delta": cell_exploratory["tLung_epithelial"]["paired_median_delta_high_minus_low"],
        },
        "sensitivity_tumor_sites": {
            "universe": "tLung + tL/B + mLN + mBrain + PE epithelial",
            "n_samples": tumor["n"],
            "spearman_rho": tumor["value"],
            "spearman_p": tumor["p_value"],
            "paired_n": paired_u["n"],
            "paired_p": paired_u["p_value"],
            "paired_median_delta": cell_exploratory["tumor_epithelial"]["paired_median_delta_high_minus_low"],
        },
        "nLung_contrast": {
            "n_samples": nlung["n"],
            "spearman_rho": nlung["value"],
            "spearman_p": nlung["p_value"],
        },
        "cell_exploratory": cell_exploratory,
        "verdict_notes": [
            "CytoTRACE2 was not run.",
            "Primary n is tLung samples with ≥20 epithelial cells, not cells.",
            "GSE131907 is treatment-naive; no ICI/MPR test.",
            "CLDN4 only.",
        ],
    }
    with (RES / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    # Figures
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6))
    panels = [
        (axes[0], sample_tables["tLung_epithelial"], "tLung epithelium", primary),
        (axes[1], sample_tables["tumor_epithelial"], "tumor-site epithelium", tumor),
    ]
    for ax, samp, title, st in panels:
        elig = samp[samp["eligible"] == 1]
        for origin, g in elig.groupby("Sample_Origin"):
            ax.scatter(g["mean_CLDN4"], g["mean_cytotrace_like"], s=36, alpha=0.85, label=origin)
        ax.set_xlabel("sample-mean CLDN4 (log1p CP10k)")
        ax.set_ylabel("sample-mean CytoTRACE-like")
        ax.set_title(f"{title}\nn={st['n']}  ρ={fmt_r(st['value'])}  p={fmt_p(st['p_value'])}")
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "gse131907_cldn4_vs_potency.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))
    for ax, name, title in (
        (axes[0], "tLung_epithelial", "tLung"),
        (axes[1], "tumor_epithelial", "tumor sites"),
    ):
        samp = sample_tables[name]
        paired = samp[(samp["eligible"] == 1) & (samp["n_CLDN4_high"] >= MIN_STRATUM) & (samp["n_CLDN4_low"] >= MIN_STRATUM)]
        x = np.array([0, 1])
        for _, row in paired.iterrows():
            ax.plot(x, [row["mean_cyto_CLDN4_low"], row["mean_cyto_CLDN4_high"]], color="0.7", lw=0.8)
        ax.scatter(np.zeros(len(paired)), paired["mean_cyto_CLDN4_low"], s=28, label="CLDN4 low Q1")
        ax.scatter(np.ones(len(paired)), paired["mean_cyto_CLDN4_high"], s=28, label="CLDN4 high Q4")
        st = next(s for s in stat_rows if s["analysis"] == name and s["comparison"].startswith("paired"))
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["CLDN4 Q1", "CLDN4 Q4"])
        ax.set_ylabel("mean CytoTRACE-like")
        ax.set_title(f"{title} paired n={st['n']}  p={fmt_p(st['p_value'])}")
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "gse131907_cldn4_quartile_potency.png", dpi=160)
    plt.close(fig)

    print(json.dumps(summary, indent=2))
    print(stats_df.to_string(index=False))


if __name__ == "__main__":
    main()
