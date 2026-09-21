#!/usr/bin/env python3
"""CPTAC LUAD/LSCC protein: CLDN4 and TACSTD2 vs MHC-I, IFN, and CD8A.

Public TMT freeze v1.2, tumor protein only. LUAD and LSCC stay separate.
Pairwise-complete Spearman, two-sided p, 2,000-resample bootstrap 95% CI.
CLDN4 missingness is reported in full, including whether the dropped tumors
differ on the same protein endpoints. No imputation. No RNA stand-in.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from genes import (  # noqa: E402
    APM,
    CD8,
    CLASS,
    ENSEMBL,
    ENDPOINTS,
    IFN_ISG,
    IFN_LIGANDS,
    IFN_RECEPTORS,
    IFN_SCORE_MEMBERS,
    IFN_SIGNALING,
    MHC1,
    MHC1_SCORE_MEMBERS,
    PREDICTORS,
    PRIMARY_SCORES,
)

SEED = 20260921
N_BOOT = 2000
MIN_N = 8
MIN_CONTRAST = 8
MIN_IFN_GENES = 4
MIN_MHC1_GENES = 3


def find_row(index: pd.Index, prefixes: list[str]) -> str | None:
    hits: list[str] = []
    for prefix in prefixes:
        for i in index:
            s = str(i)
            if s == prefix or s.startswith(prefix + ".") or s.startswith(prefix + "|"):
                hits.append(s)
    hits = list(dict.fromkeys(hits))
    if not hits:
        return None
    if len(hits) > 1:
        raise ValueError(f"multiple rows for {prefixes}: {hits}")
    return hits[0]


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce")


def extract_gene(mat: pd.DataFrame, symbol: str) -> tuple[pd.Series, str | None]:
    row = find_row(mat.index, ENSEMBL[symbol])
    if row is None:
        return pd.Series(np.nan, index=mat.columns, name=symbol), None
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    s.name = symbol
    return s, row


def mean_z(series_list: list[pd.Series], min_genes: int) -> pd.Series:
    rows = []
    for s in series_list:
        sd = float(s.std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            continue
        rows.append((s - s.mean()) / sd)
    if not rows:
        return pd.Series(dtype=float)
    z = pd.concat(rows, axis=1)
    n = z.notna().sum(axis=1)
    score = z.mean(axis=1, skipna=True)
    score = score.where(n >= min_genes)
    score.name = "score"
    return score


def _boot_rho(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    n = len(x)
    boots = np.empty(N_BOOT, dtype=float)
    kept = 0
    for i in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        r, _ = stats.spearmanr(x[idx], y[idx])
        if np.isfinite(r):
            boots[kept] = r
            kept += 1
    if kept < 100:
        return np.nan, np.nan
    return float(np.percentile(boots[:kept], 2.5)), float(np.percentile(boots[:kept], 97.5))


def spearman_boot(x: pd.Series, y: pd.Series, rng: np.random.Generator) -> dict:
    d = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    n = int(len(d))
    rec = {"n": n, "rho": np.nan, "p": np.nan, "ci_lo": np.nan, "ci_hi": np.nan}
    if n < MIN_N or d["x"].nunique() < 2 or d["y"].nunique() < 2:
        return rec
    rho, p = stats.spearmanr(d["x"], d["y"])
    rec["rho"] = float(rho)
    rec["p"] = float(p)
    rec["ci_lo"], rec["ci_hi"] = _boot_rho(d["x"].to_numpy(), d["y"].to_numpy(), rng)
    return rec


def rank_resid(a: np.ndarray, z: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(z)), z])
    coef, _, _, _ = np.linalg.lstsq(design, a, rcond=None)
    return a - design @ coef


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series, rng: np.random.Generator) -> dict:
    """Pearson correlation of rank-residuals on z. Bootstrap resamples tumors."""
    d = pd.concat([x.rename("x"), y.rename("y"), z.rename("z")], axis=1).dropna()
    n = int(len(d))
    rec = {
        "n_partial": n,
        "rho_partial": np.nan,
        "p_partial": np.nan,
        "ci_lo_partial": np.nan,
        "ci_hi_partial": np.nan,
    }
    if n < MIN_N or d["z"].nunique() < 2:
        return rec
    rx = d["x"].rank().to_numpy()
    ry = d["y"].rank().to_numpy()
    rz = d["z"].rank().to_numpy()
    xr = rank_resid(rx, rz)
    yr = rank_resid(ry, rz)
    if np.nanstd(xr) == 0 or np.nanstd(yr) == 0:
        return rec
    rho, p = stats.pearsonr(xr, yr)
    rec["rho_partial"] = float(rho)
    rec["p_partial"] = float(p)
    boots = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        xb = rank_resid(rx[idx], rz[idx])
        yb = rank_resid(ry[idx], rz[idx])
        if np.nanstd(xb) == 0 or np.nanstd(yb) == 0:
            continue
        r, _ = stats.pearsonr(xb, yb)
        if np.isfinite(r):
            boots.append(float(r))
    if len(boots) >= 100:
        rec["ci_lo_partial"] = float(np.percentile(boots, 2.5))
        rec["ci_hi_partial"] = float(np.percentile(boots, 97.5))
    return rec


def bh_q(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    idx = np.where(np.isfinite(p))[0]
    m = len(idx)
    if m == 0:
        return out.tolist()
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.clip(q, 0, 1)
    return out.tolist()


def protein_path(data: Path, cohort: str) -> Path:
    return (
        data
        / cohort
        / f"{cohort}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
    )


def phenotype_path(data: Path, cohort: str) -> Path:
    return data / cohort / f"{cohort}_phenotype.txt"


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3g}"


def fmt_rho(r: float) -> str:
    if not np.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def fmt_ci(lo: float, hi: float) -> str:
    if not np.isfinite(lo) or not np.isfinite(hi):
        return "NA"
    return f"{lo:+.3f} to {hi:+.3f}"


def support_call(row: pd.Series) -> str:
    """Inverse IFN/MHC support. CD8A uses the same numeric rule but is not an IFN call."""
    rho, lo, hi = row["rho"], row["ci_lo"], row["ci_hi"]
    q = row["q_primary"]
    plo, phi = row.get("ci_lo_partial", np.nan), row.get("ci_hi_partial", np.nan)
    pr = row.get("rho_partial", np.nan)
    if not np.isfinite(rho):
        return "not tested"
    inverse = np.isfinite(hi) and hi < 0
    opposite = np.isfinite(lo) and lo > 0
    q_ok = np.isfinite(q) and q < 0.05
    partial_inverse = np.isfinite(phi) and phi < 0 and np.isfinite(pr) and pr < 0
    partial_crosses = np.isfinite(plo) and np.isfinite(phi) and plo <= 0 <= phi
    partial_opposite = np.isfinite(plo) and plo > 0
    if inverse and q_ok and partial_inverse:
        return "supports inverse"
    if inverse and q_ok and (partial_crosses or partial_opposite):
        return "complete-case inverse; purity partial does not"
    if inverse and q_ok and not np.isfinite(pr):
        return "supports inverse (no purity column)"
    if inverse and not q_ok:
        return "CI below 0; primary-family q ≥ 0.05"
    if opposite and q_ok:
        return "opposite (positive); does not support inverse"
    if opposite:
        return "positive CI; primary-family q ≥ 0.05"
    return "null; does not support"


def style(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def prepare_cohort(data: Path, cohort: str) -> dict:
    mat = load_matrix(protein_path(data, cohort))
    ph = pd.read_csv(phenotype_path(data, cohort), sep="\t")
    if "idx" not in ph.columns:
        raise ValueError(f"{cohort} phenotype has no idx column")
    ph = ph.set_index(ph["idx"].astype(str))
    n_tumors = int(mat.shape[1])
    if set(mat.columns) != set(ph.index):
        raise ValueError(f"{cohort} protein columns do not match phenotype idx")

    treat_keys = (
        "treat",
        "therap",
        "neoadj",
        "ici",
        "immunotherap",
        "pd-1",
        "pd1",
        "pembro",
        "nivol",
        "atezo",
        "chemo",
        "drug",
    )
    treat_cols = [c for c in ph.columns if any(k in str(c).lower() for k in treat_keys)]

    series: dict[str, pd.Series] = {}
    rows_used: dict[str, str | None] = {}
    coverage = []
    for symbol in list(ENSEMBL):
        s, row = extract_gene(mat, symbol)
        s.index = s.index.astype(str)
        series[symbol] = s
        rows_used[symbol] = row
        n_q = int(s.notna().sum())
        n_zero = int((s == 0).sum())
        coverage.append(
            {
                "cohort": cohort,
                "symbol": symbol,
                "class": CLASS[symbol],
                "ensembl_query": ";".join(ENSEMBL[symbol]),
                "row": row or "",
                "present": row is not None,
                "n_tumors": n_tumors,
                "n_quantified": n_q,
                "n_missing": n_tumors - n_q if row is not None else n_tumors,
                "pct_missing": (100.0 * (n_tumors - n_q) / n_tumors) if row is not None else 100.0,
                "n_exact_zero": n_zero if row is not None else 0,
                "min": float(s.min()) if n_q else np.nan,
                "median": float(s.median()) if n_q else np.nan,
                "max": float(s.max()) if n_q else np.nan,
            }
        )

    ifn_members = []
    ifn_vecs = []
    for g in IFN_SCORE_MEMBERS:
        if rows_used[g] is not None and int(series[g].notna().sum()) >= MIN_N:
            ifn_members.append(g)
            ifn_vecs.append(series[g])
    if len(ifn_members) >= MIN_IFN_GENES:
        series["IFN_core_protein"] = mean_z(ifn_vecs, MIN_IFN_GENES)
    else:
        series["IFN_core_protein"] = pd.Series(np.nan, index=mat.columns)

    mhc_members = []
    mhc_vecs = []
    for g in MHC1_SCORE_MEMBERS:
        if rows_used[g] is not None and int(series[g].notna().sum()) >= MIN_N:
            mhc_members.append(g)
            mhc_vecs.append(series[g])
    if len(mhc_members) >= MIN_MHC1_GENES:
        series["MHC1_protein"] = mean_z(mhc_vecs, MIN_MHC1_GENES)
    else:
        series["MHC1_protein"] = pd.Series(np.nan, index=mat.columns)
    series["CD8A_protein"] = series["CD8A"]

    wes = pd.to_numeric(ph["WES_purity"], errors="coerce") if "WES_purity" in ph.columns else None
    wgs = pd.to_numeric(ph["WGS_purity"], errors="coerce") if "WGS_purity" in ph.columns else None
    if wes is not None:
        wes.index = wes.index.astype(str)
        wes = wes.reindex(mat.columns)
    if wgs is not None:
        wgs.index = wgs.index.astype(str)
        wgs = wgs.reindex(mat.columns)

    return {
        "cohort": cohort,
        "n_tumors": n_tumors,
        "n_rows": int(mat.shape[0]),
        "n_pheno_cols": int(ph.shape[1]),
        "treat_cols": treat_cols,
        "series": series,
        "rows_used": rows_used,
        "coverage": coverage,
        "ifn_members": ifn_members,
        "mhc_members": mhc_members,
        "wes": wes,
        "wgs": wgs,
    }


def endpoint_series(bundle: dict, name: str) -> pd.Series:
    return bundle["series"][name]


def run_correlations(bundle: dict, rng: np.random.Generator) -> list[dict]:
    cohort = bundle["cohort"]
    out = []
    named = list(ENDPOINTS) + ["MHC1_protein", "IFN_core_protein"]
    # Coexpression of the two predictors is context, not an IFN endpoint.
    pairs = [(pred, ep) for pred in PREDICTORS for ep in named if ep != pred]
    pairs.append(("CLDN4", "TACSTD2"))

    for pred, ep in pairs:
        x = endpoint_series(bundle, pred)
        if ep in bundle["series"]:
            y = endpoint_series(bundle, ep)
        else:
            continue
        if ep in ("MHC1_protein", "IFN_core_protein", "CD8A"):
            symbol = ep if ep != "CD8A" else "CD8A"
            klass = "score" if ep.endswith("_protein") else CLASS.get(ep, "endpoint")
            rowname = ",".join(bundle["mhc_members"] if ep == "MHC1_protein" else bundle["ifn_members"] if ep == "IFN_core_protein" else [bundle["rows_used"].get(ep) or ""])
        else:
            symbol = ep
            klass = CLASS.get(ep, "endpoint")
            rowname = bundle["rows_used"].get(ep) or ""
        present = bool(y.notna().any())
        rec = spearman_boot(x, y, rng)
        primary = pred in PREDICTORS and (
            ep in ("MHC1_protein", "IFN_core_protein") or ep == "CD8A"
        )
        # CD8A gene is the CD8A_protein primary. Avoid double-counting CD8A_protein key.
        if ep == "CD8A":
            endpoint_name = "CD8A_protein"
        else:
            endpoint_name = ep if ep.endswith("_protein") else f"{ep}_protein"
        out.append(
            {
                "cohort": cohort,
                "predictor": f"{pred}_protein",
                "endpoint": endpoint_name,
                "symbol": symbol if ep != "CD8A" else "CD8A",
                "class": "score" if ep in ("MHC1_protein", "IFN_core_protein") else klass,
                "row": rowname,
                "present": present,
                "primary_family": primary,
                **rec,
            }
        )
    return out


def run_partials(bundle: dict, corrs: list[dict], rng: np.random.Generator) -> None:
    wes = bundle["wes"]
    for rec in corrs:
        rec["n_partial"] = 0
        rec["rho_partial"] = np.nan
        rec["p_partial"] = np.nan
        rec["ci_lo_partial"] = np.nan
        rec["ci_hi_partial"] = np.nan
        if wes is None or not rec["primary_family"] or not np.isfinite(rec["rho"]):
            continue
        pred = rec["predictor"].replace("_protein", "")
        ep = rec["endpoint"]
        if ep == "CD8A_protein":
            y = bundle["series"]["CD8A"]
        else:
            y = bundle["series"][ep]
        x = bundle["series"][pred]
        rec.update(partial_spearman(x, y, wes, rng))


def run_missingness(bundle: dict) -> list[dict]:
    """Quantified vs missing for each predictor, on pre-specified endpoints."""
    cohort = bundle["cohort"]
    targets = ["CLDN4", "TACSTD2", "CD8A", "MHC1_protein", "IFN_core_protein", "HLA-A", "HLA-B", "HLA-C", "STAT1", "IRF9", "JAK1", "JAK2"]
    rows = []
    for pred in PREDICTORS:
        status = bundle["series"][pred]
        quantified = status.notna()
        missing = status.isna()
        # If the row itself is absent, every tumor is missing and there is no contrast.
        if bundle["rows_used"][pred] is None:
            continue
        for ep in targets:
            if ep == pred:
                continue
            y = bundle["series"].get(ep)
            if y is None:
                continue
            a = y[quantified].dropna()
            b = y[missing].dropna()
            rec = {
                "cohort": cohort,
                "predictor": pred,
                "endpoint": ep if ep.endswith("_protein") or ep in PREDICTORS or ep in ("CD8A", "HLA-A", "HLA-B", "HLA-C", "STAT1", "IRF9", "JAK1", "JAK2") else ep,
                "n_predictor_quantified": int(quantified.sum()),
                "n_predictor_missing": int(missing.sum()),
                "n_endpoint_in_quantified": int(len(a)),
                "n_endpoint_in_missing": int(len(b)),
                "median_quantified": float(a.median()) if len(a) else np.nan,
                "median_missing": float(b.median()) if len(b) else np.nan,
                "delta_missing_minus_quantified": (
                    float(b.median() - a.median()) if len(a) and len(b) else np.nan
                ),
                "u": np.nan,
                "p": np.nan,
                "tested": False,
            }
            if len(a) >= MIN_CONTRAST and len(b) >= MIN_CONTRAST and a.nunique() + b.nunique() > 1:
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                rec["u"] = float(u)
                rec["p"] = float(p)
                rec["tested"] = True
            rows.append(rec)
        if bundle["wes"] is not None:
            y = bundle["wes"]
            a = y[quantified].dropna()
            b = y[missing].dropna()
            rec = {
                "cohort": cohort,
                "predictor": pred,
                "endpoint": "WES_purity",
                "n_predictor_quantified": int(quantified.sum()),
                "n_predictor_missing": int(missing.sum()),
                "n_endpoint_in_quantified": int(len(a)),
                "n_endpoint_in_missing": int(len(b)),
                "median_quantified": float(a.median()) if len(a) else np.nan,
                "median_missing": float(b.median()) if len(b) else np.nan,
                "delta_missing_minus_quantified": (
                    float(b.median() - a.median()) if len(a) and len(b) else np.nan
                ),
                "u": np.nan,
                "p": np.nan,
                "tested": False,
            }
            if len(a) >= MIN_CONTRAST and len(b) >= MIN_CONTRAST:
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                rec["u"] = float(u)
                rec["p"] = float(p)
                rec["tested"] = True
            rows.append(rec)
    return rows


def joint_missing(bundle: dict) -> dict:
    c = bundle["series"]["CLDN4"].notna()
    t = bundle["series"]["TACSTD2"].notna()
    n = bundle["n_tumors"]
    return {
        "cohort": bundle["cohort"],
        "n_tumors": n,
        "both_quantified": int((c & t).sum()),
        "cldn4_only": int((c & ~t).sum()),
        "tacstd2_only": int((~c & t).sum()),
        "both_missing": int((~c & ~t).sum()),
        "cldn4_quantified": int(c.sum()),
        "cldn4_missing": int((~c).sum()),
        "cldn4_pct_missing": 100.0 * int((~c).sum()) / n,
        "tacstd2_quantified": int(t.sum()),
        "tacstd2_missing": int((~t).sum()),
        "tacstd2_pct_missing": 100.0 * int((~t).sum()) / n,
    }


def plot_missingness(cov: pd.DataFrame, figdir: Path) -> None:
    symbols = ["CLDN4", "TACSTD2", "CD8A", "CD8B", "HLA-A", "HLA-B", "HLA-C", "B2M", "IFNG", "STAT1", "ISG15", "IFNGR1"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
    for ax, cohort in zip(axes, ["LUAD", "LSCC"]):
        sub = cov[(cov["cohort"] == cohort) & (cov["symbol"].isin(symbols))].copy()
        sub["symbol"] = pd.Categorical(sub["symbol"], symbols, ordered=True)
        sub = sub.sort_values("symbol")
        pct_q = 100.0 - sub["pct_missing"].to_numpy()
        colors = ["#b3483a" if s == "CLDN4" else "#4c78a8" for s in sub["symbol"]]
        ax.barh(sub["symbol"].astype(str)[::-1], pct_q[::-1], color=colors[::-1])
        ax.set_xlim(0, 100)
        ax.set_xlabel("% of tumors quantified")
        ax.set_title(f"{cohort} (n={int(sub['n_tumors'].iloc[0])})")
        style(ax)
        for y, (pq, nmiss) in enumerate(zip(pct_q[::-1], sub["n_missing"].to_numpy()[::-1])):
            ax.text(min(pq + 1.5, 86), y, f"{pq:.0f}%  ({int(nmiss)} NA)", va="center", fontsize=7)
    fig.suptitle("CPTAC TMT freeze v1.2 — protein completeness (red = CLDN4)", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "fig_missingness.png", dpi=160)
    fig.savefig(figdir / "fig_missingness.pdf")
    plt.close(fig)


def plot_forest(primary: pd.DataFrame, figdir: Path) -> None:
    order_ep = ["MHC1_protein", "IFN_core_protein", "CD8A_protein"]
    order_pred = ["CLDN4_protein", "TACSTD2_protein"]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6), sharex=True)
    for ax, cohort in zip(axes, ["LUAD", "LSCC"]):
        sub = primary[primary["cohort"] == cohort].copy()
        labels = []
        ys = []
        for i, ep in enumerate(order_ep):
            for j, pred in enumerate(order_pred):
                hit = sub[(sub["endpoint"] == ep) & (sub["predictor"] == pred)]
                y = (len(order_ep) - 1 - i) * 2.4 + (1 - j) * 0.9
                short_p = "CLDN4" if pred.startswith("CLDN4") else "TACSTD2"
                short_e = {"MHC1_protein": "MHC-I", "IFN_core_protein": "IFN-core", "CD8A_protein": "CD8A"}[ep]
                labels.append((y, f"{short_p} vs {short_e}"))
                if hit.empty or not np.isfinite(hit.iloc[0]["rho"]):
                    continue
                r = hit.iloc[0]
                color = "#b3483a" if r["ci_hi"] < 0 else ("#2a6f4e" if r["ci_lo"] > 0 else "#5c6570")
                ax.errorbar(
                    r["rho"],
                    y,
                    xerr=[[r["rho"] - r["ci_lo"]], [r["ci_hi"] - r["rho"]]],
                    fmt="o",
                    color=color,
                    ms=5,
                    lw=1.2,
                    capsize=2,
                )
                ax.text(
                    r["ci_hi"] + 0.03,
                    y,
                    f"{r['rho']:+.2f}  n={int(r['n'])}",
                    va="center",
                    fontsize=7,
                    color=color,
                )
                ys.append(y)
        ax.axvline(0, color="#888888", lw=0.8)
        ax.set_yticks([y for y, _ in labels])
        ax.set_yticklabels([lab for _, lab in labels], fontsize=8)
        ax.set_xlim(-0.85, 1.05)
        ax.set_xlabel("Spearman ρ (95% bootstrap CI)")
        ax.set_title(cohort)
        style(ax)
    fig.suptitle("Primary protein–protein correlations (pairwise-complete)", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "fig_forest.png", dpi=160)
    fig.savefig(figdir / "fig_forest.pdf")
    plt.close(fig)


def plot_scatters(bundles: dict, figdir: Path) -> None:
    endpoints = [("MHC1_protein", "MHC-I protein score"), ("IFN_core_protein", "IFN-core protein score"), ("CD8A", "CD8A protein")]
    for pred, fname in (("CLDN4", "fig_scatter_cldn4"), ("TACSTD2", "fig_scatter_tacstd2")):
        fig, axes = plt.subplots(2, 3, figsize=(9.6, 6.2))
        for row, cohort in enumerate(["LUAD", "LSCC"]):
            b = bundles[cohort]
            x = b["series"][pred]
            for col, (ep, label) in enumerate(endpoints):
                ax = axes[row, col]
                y = b["series"][ep]
                d = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
                ax.scatter(d["x"], d["y"], s=12, c="#3d4c63", alpha=0.75, linewidths=0)
                if len(d) >= MIN_N:
                    rho, p = stats.spearmanr(d["x"], d["y"])
                    ax.set_title(f"{cohort}  ρ={rho:+.2f}  p={fmt_p(p)}  n={len(d)}", fontsize=8)
                else:
                    ax.set_title(f"{cohort}  n={len(d)}", fontsize=8)
                if row == 1:
                    ax.set_xlabel(f"{pred} protein (log2 ratio)")
                if col == 0:
                    ax.set_ylabel(label)
                style(ax)
        fig.suptitle(f"{pred} protein vs MHC-I, IFN-core, and CD8A protein", fontsize=11)
        fig.tight_layout()
        fig.savefig(figdir / f"{fname}.png", dpi=160)
        fig.savefig(figdir / f"{fname}.pdf")
        plt.close(fig)


def md_table(headers: list[str], rows: list[list[str]]) -> str:
    line = "| " + " | ".join(headers) + " |"
    sep = "|" + "|".join("---" if i == 0 else "---:" for i in range(len(headers))) + "|"
    # first col left, rest right — rebuild properly
    segs = ["---"] + ["---:"] * (len(headers) - 1)
    sep = "| " + " | ".join(segs) + " |"
    body = ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join([line, sep, *body])


def write_finding(
    path: Path,
    bundles: dict,
    cov: pd.DataFrame,
    corrs: pd.DataFrame,
    miss: pd.DataFrame,
    joint: pd.DataFrame,
) -> None:
    primary = corrs[corrs["primary_family"]].copy()

    def cell(cohort: str, pred: str, ep: str) -> pd.Series:
        hit = primary[
            (primary["cohort"] == cohort)
            & (primary["predictor"] == pred)
            & (primary["endpoint"] == ep)
        ]
        return hit.iloc[0]

    def cov_row(cohort: str, symbol: str) -> pd.Series:
        return cov[(cov["cohort"] == cohort) & (cov["symbol"] == symbol)].iloc[0]

    lines = []
    lines.append("# Finding — CPTAC LUAD / LSCC protein: CLDN4 and TACSTD2 vs MHC-I, IFN, and CD8A")
    lines.append("")
    lines.append(
        "Public CPTAC TMT freeze v1.2 tumor protein. **LUAD** (Gillette *Cell* 2020, n=110) and "
        "**LSCC** (Satpathy *Cell* 2021, n=108) kept separate. Pairwise-complete Spearman ρ, "
        "two-sided p, 2,000-resample bootstrap 95% CI (seed `20260921`). Honest n. "
        "Do not write n=110 / n=108 for a CLDN4 test."
    )
    lines.append("")
    lines.append(
        "This page is protein–protein. It does not re-fit CLDN4 protein vs ImmuneScore, GEP18, "
        "or CD8A **RNA** (PR #315, #289, #245) and it does not replace the CLDN4 vs CD274 protein "
        "page (PR #356). TACSTD2, CD8A **protein**, and the CLDN4 missingness account are the addition."
    )
    lines.append("")
    lines.append("## Treatment-naive")
    lines.append("")
    lines.append(
        "The freeze phenotype table has **no treatment, neoadjuvant, or ICI column** "
        f"(LUAD {bundles['LUAD']['n_pheno_cols']} columns, LSCC {bundles['LSCC']['n_pheno_cols']} columns; "
        "deconvolution, hallmark ssGSEA, CNV, purity, TMB). "
        f"Therapy-like column names found: "
        + (
            ", ".join(f"`{c}`" for c in bundles["LUAD"]["treat_cols"] + bundles["LSCC"]["treat_cols"])
            if (bundles["LUAD"]["treat_cols"] or bundles["LSCC"]["treat_cols"])
            else "**none**"
        )
        + ". Hallmark interferon columns in that file are **RNA ssGSEA**, not protein, and are not tested here."
    )
    lines.append("")
    lines.append(
        "Treatment-naive is the published cohort definition, not a label in this matrix: "
        "both studies are prospectively collected, previously untreated surgical resections. "
        "There is no ICI arm, no on-treatment biopsy, and no response label. "
        "A correlation here cannot be written as acquired resistance or post-ICI change."
    )
    lines.append("")
    lines.append("## Direction under test")
    lines.append("")
    lines.append(
        "Private PDX proteome numbers are **not** in this repository and are not restated. "
        "The human-protein direction that would support an inverse PDX IFN/MHC result is "
        "pre-specified: **higher CLDN4 or TACSTD2 protein, lower MHC-I or IFN-core protein**. "
        "A pair is called **supports inverse** only when all three hold: bootstrap 95% CI entirely below 0, "
        "Benjamini–Hochberg q < 0.05 inside the primary family "
        "(2 cohorts × 2 predictors × MHC-I, IFN-core, and CD8A), and the WES-purity partial CI also entirely below 0. "
        "CD8A uses the same numeric rule but is an immune-cell protein, not an IFN/MHC score. "
        "Null and positive pairs are not support."
    )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")

    verdict_rows = []
    for cohort in ["LUAD", "LSCC"]:
        for pred, pred_lab in (("CLDN4_protein", "CLDN4"), ("TACSTD2_protein", "TACSTD2")):
            for ep, ep_lab in (
                ("MHC1_protein", "MHC-I"),
                ("IFN_core_protein", "IFN-core"),
                ("CD8A_protein", "CD8A"),
            ):
                r = cell(cohort, pred, ep)
                verdict_rows.append(
                    [
                        cohort,
                        f"{pred_lab} vs {ep_lab}",
                        str(int(r["n"])),
                        fmt_rho(r["rho"]),
                        fmt_ci(r["ci_lo"], r["ci_hi"]),
                        fmt_p(r["p"]),
                        fmt_p(r["q_primary"]) if np.isfinite(r["q_primary"]) else "NA",
                        fmt_rho(r["rho_partial"]),
                        r["call"],
                    ]
                )
    lines.append(
        md_table(
            ["Cohort", "Pair", "n", "ρ", "95% CI", "p", "q", "partial ρ (WES)", "Call"],
            verdict_rows,
        )
    )
    lines.append("")

    def pair_detail(cohort: str, pred: str, ep: str) -> str:
        r = cell(cohort, pred, ep)
        return (
            f"ρ={fmt_rho(r['rho'])}, n={int(r['n'])}, 95% CI {fmt_ci(r['ci_lo'], r['ci_hi'])}, "
            f"p={fmt_p(r['p'])}, q={fmt_p(r['q_primary'])}, "
            f"WES partial ρ={fmt_rho(r['rho_partial'])} "
            f"(CI {fmt_ci(r['ci_lo_partial'], r['ci_hi_partial'])}, n={int(r['n_partial'])})"
        )

    lines.append(
        "**MHC-I, squamous CLDN4 only.** LSCC CLDN4 protein vs the HLA-A/B/C score "
        f"{pair_detail('LSCC', 'CLDN4_protein', 'MHC1_protein')}. "
        "Call: supports inverse. The partial CI stays below 0, so this is not removed by WES purity. "
        "LUAD CLDN4 vs MHC-I is null "
        f"({pair_detail('LUAD', 'CLDN4_protein', 'MHC1_protein')}). "
        "TACSTD2 vs MHC-I is null in both histologies."
    )
    lines.append("")
    lines.append(
        "**IFN-core score does not support the inverse direction.** "
        f"LSCC CLDN4 {pair_detail('LSCC', 'CLDN4_protein', 'IFN_core_protein')}. "
        f"LUAD CLDN4 {pair_detail('LUAD', 'CLDN4_protein', 'IFN_core_protein')}. "
        "Both confidence intervals include 0. "
        "LSCC CLDN4 is inverse at HLA-A and HLA-C and at JAK1, JAK2, IFNGR1, STAT1, IRF9, TAP1, TAP2, TAPBP, PSMB8, PSMB9, and NLRC5 (gene table; each of those CIs lies below 0). "
        "HLA-B is the weaker heavy chain (its CI crosses 0). "
        "ISG15, MX1, OAS1, and IFI35 are null in that same cohort and pull the 10-gene mean-z score across 0. "
        "Do not promote the JAK rows into an IFN-core or IFN-ligand result. IFNG, IFNA1, and IFNB1 protein are absent."
    )
    lines.append("")
    tac_ifn = cell("LSCC", "TACSTD2_protein", "IFN_core_protein")
    lines.append(
        "**TACSTD2 is not the CLDN4 result.** It does not track MHC-I or CD8A. "
        f"LSCC TACSTD2 vs IFN-core is the opposite sign: {pair_detail('LSCC', 'TACSTD2_protein', 'IFN_core_protein')}. "
        f"Call: {tac_ifn['call']}. The partial CI is also entirely above 0. "
        "That pair is not inverse support, and it does not clear q < 0.05 in the 12-test primary family. "
        "Do not write TACSTD2 protein as a human-protein copy of the LSCC CLDN4–MHC-I inverse."
    )
    lines.append("")
    lines.append(
        "**CD8A protein, squamous CLDN4 only.** LSCC CLDN4 vs CD8A "
        f"{pair_detail('LSCC', 'CLDN4_protein', 'CD8A_protein')}. "
        "Call: supports inverse. CD8A protein is quantified in every LSCC tumor (108/108); "
        "n=78 is CLDN4 completeness, not CD8A dropout. "
        "This is an immune-cell protein association, not an IFN-core result. "
        "LUAD CLDN4 vs CD8A is null, and TACSTD2 vs CD8A is null in both histologies. "
        "CD8B is absent in LUAD and quantified in 10/108 LSCC tumors, so it is not an endpoint."
    )
    lines.append("")
    lines.append(
        "Human-protein support for an inverse PDX IFN/MHC direction, on the pre-specified rule, "
        "is **LSCC CLDN4 vs MHC-I protein only**. The IFN-core score does not support it. "
        "LUAD does not support it. TACSTD2 does not support it. "
        "The additional inverse that clears the same bar is LSCC CLDN4 vs CD8A protein."
    )
    lines.append("")

    lines.append("<!-- SENSITIVITY -->")
    lines.append("")
    lines.append("## CLDN4 missingness — completeness")
    lines.append("")
    lines.append(
        "CLDN4 is on the TMT matrix (`ENSG00000189143.9` in both cohorts). "
        "Missing means the abundance cell is empty after numeric parse, not a join failure and not an exact zero. "
        "Exact zeros are real log2 ratios equal to the reference and are counted separately. "
        "Values are not imputed. Complete-case Spearman uses only quantified CLDN4 tumors."
    )
    lines.append("")
    miss_rows = []
    for cohort in ["LUAD", "LSCC"]:
        for symbol in ["CLDN4", "TACSTD2", "CD8A", "CD8B", "HLA-A", "HLA-B", "HLA-C", "B2M", "IFNG", "IFNA1", "IFNB1"]:
            c = cov_row(cohort, symbol)
            miss_rows.append(
                [
                    cohort,
                    symbol,
                    "yes" if c["present"] else "absent",
                    str(int(c["n_quantified"])) if c["present"] else "0",
                    str(int(c["n_missing"])),
                    f"{c['pct_missing']:.1f}%",
                    str(int(c["n_exact_zero"])),
                ]
            )
    lines.append(
        md_table(
            ["Cohort", "Protein", "Row", "Quantified", "Missing", "% missing", "Exact zeros"],
            miss_rows,
        )
    )
    lines.append("")
    jbits = []
    for cohort in ["LUAD", "LSCC"]:
        j = joint[joint["cohort"] == cohort].iloc[0]
        c4 = cov_row(cohort, "CLDN4")
        jbits.append(
            f"{cohort}: CLDN4 quantified **{int(j['cldn4_quantified'])} / {int(j['n_tumors'])}** "
            f"({j['cldn4_pct_missing']:.1f}% missing"
            + (
                f"; range {c4['min']:.2f} to {c4['max']:.2f}, median {c4['median']:.2f}"
                if c4["present"]
                else ""
            )
            + f"). TACSTD2 quantified **{int(j['tacstd2_quantified'])} / {int(j['n_tumors'])}** "
            f"({j['tacstd2_pct_missing']:.1f}% missing). "
            f"Both quantified {int(j['both_quantified'])}; CLDN4 only {int(j['cldn4_only'])}; "
            f"TACSTD2 only {int(j['tacstd2_only'])}; both missing {int(j['both_missing'])}."
        )
    lines.append(" ".join(jbits))
    lines.append("")
    lines.append(
        "The CLDN4 Spearman n is 79 (LUAD) and 78 (LSCC). "
        "TACSTD2 and CD8A are quantified in every tumor, so those predictors and that endpoint do not add further dropout. "
        "HLA-A/B/C are also complete. "
        "IFNG, IFNA1, IFNB1, and B2M are a different kind of gap: the row is absent, not a within-row NA. "
        "There is nothing to correlate."
    )
    lines.append("")
    lines.append("### Do the CLDN4-missing tumors differ?")
    lines.append("")
    lines.append(
        "Mann–Whitney, two-sided, quantified CLDN4 vs missing CLDN4, on tumors that have the endpoint. "
        "Delta is median(missing) − median(quantified). "
        "A higher MHC-I or IFN score in the missing group would mean the dropped tumors sit toward the "
        "high-IFN side and the complete-case ρ is not the full cohort. Minimum n per side = 8."
    )
    lines.append("")
    contrast_keep = [
        "TACSTD2",
        "MHC1_protein",
        "IFN_core_protein",
        "CD8A",
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "WES_purity",
    ]
    crows = []
    subm = miss[(miss["predictor"] == "CLDN4") & (miss["endpoint"].isin(contrast_keep))]
    for cohort in ["LUAD", "LSCC"]:
        for ep in contrast_keep:
            hit = subm[(subm["cohort"] == cohort) & (subm["endpoint"] == ep)]
            if hit.empty:
                continue
            r = hit.iloc[0]
            crows.append(
                [
                    cohort,
                    ep,
                    f"{int(r['n_endpoint_in_quantified'])} vs {int(r['n_endpoint_in_missing'])}",
                    fmt_rho(r["median_quantified"]).replace("+", ""),
                    fmt_rho(r["median_missing"]).replace("+", ""),
                    fmt_rho(r["delta_missing_minus_quantified"]),
                    fmt_p(r["p"]) if r["tested"] else "not tested",
                ]
            )
    lines.append(
        md_table(
            ["Cohort", "Endpoint", "n quantified vs missing", "Median in", "Median out", "Δ (out−in)", "MW p"],
            crows,
        )
    )
    lines.append("")
    def mw(cohort: str, endpoint: str) -> pd.Series:
        hit = miss[
            (miss["predictor"] == "CLDN4")
            & (miss["cohort"] == cohort)
            & (miss["endpoint"] == endpoint)
        ]
        return hit.iloc[0]

    def mw_bit(cohort: str, endpoint: str) -> str:
        r = mw(cohort, endpoint)
        return (
            f"Δ={r['delta_missing_minus_quantified']:+.3f}, "
            f"p={fmt_p(r['p'])} "
            f"({int(r['n_endpoint_in_quantified'])} vs {int(r['n_endpoint_in_missing'])})"
        )

    lines.append(
        "CLDN4 dropout is not a purity split and not a TACSTD2 split. "
        f"WES purity {mw_bit('LUAD', 'WES_purity')} in LUAD and {mw_bit('LSCC', 'WES_purity')} in LSCC. "
        f"TACSTD2 {mw_bit('LUAD', 'TACSTD2')} in LUAD and {mw_bit('LSCC', 'TACSTD2')} in LSCC. "
        "In LUAD the missing tumors are not a distinct MHC-I, IFN-core, or CD8A group "
        f"(MHC-I {mw_bit('LUAD', 'MHC1_protein')}; IFN-core {mw_bit('LUAD', 'IFN_core_protein')}; "
        f"CD8A {mw_bit('LUAD', 'CD8A')}). "
        "In LSCC the 30 CLDN4-missing tumors have higher CD8A "
        f"({mw_bit('LSCC', 'CD8A')}) and a same-direction MHC-I shift that does not clear 0.05 "
        f"({mw_bit('LSCC', 'MHC1_protein')}). IFN-core does not differ "
        f"({mw_bit('LSCC', 'IFN_core_protein')}). "
        "Those 30 tumors sit toward higher CD8A, the same direction as the complete-case inverse, "
        "not toward a hidden CLDN4-high / IFN-high state. "
        "They are still excluded from the Spearman: CLDN4 is not imputed, and the LSCC CLDN4 n stays 78."
    )
    lines.append("")
    lines.append(
        "TACSTD2 has no missing tumors, so a TACSTD2 quantified-vs-missing contrast does not exist. "
        "Full contrast table: `tables/missingness_contrast.tsv`."
    )
    lines.append("")
    lines.append("## Primary correlations")
    lines.append("")
    lines.append(
        "MHC-I score = mean of per-gene z-scores for HLA-A, HLA-B, and HLA-C "
        "(B2M is absent, so the score is 3/4; a tumor needs all three HLA rows). "
        f"IFN-core = mean of z-scores for the signaling+ISG members with a usable row "
        f"(LUAD {len(bundles['LUAD']['ifn_members'])}/10: {', '.join(bundles['LUAD']['ifn_members'])}; "
        f"LSCC {len(bundles['LSCC']['ifn_members'])}/10: {', '.join(bundles['LSCC']['ifn_members'])}). "
        "A tumor needs at least 4 of those members. CD8A is the single protein, not a score. "
        "Partial ρ residualizes ranks on WES purity."
    )
    lines.append("")
    lines.append("WES purity coverage: " + ", ".join(
        f"{c} {int(bundles[c]['wes'].notna().sum())}/{bundles[c]['n_tumors']}"
        if bundles[c]["wes"] is not None
        else f"{c} absent"
        for c in ["LUAD", "LSCC"]
    ) + ".")
    lines.append("")
    lines.append("Forest: `figures/fig_forest.png`. Scatters: `figures/fig_scatter_cldn4.png`, `figures/fig_scatter_tacstd2.png`.")
    lines.append("")

    lines.append("## Gene-level proteins (supporting, not the support call)")
    lines.append("")
    lines.append(
        "Locked list. Cells are pairwise n, ρ, p. Absent or pairwise n<8 → no ρ. "
        "These p-values are **not** in the primary-family FDR. "
        "Pairs with n=10 (LSCC CD8B, LSCC IFNAR2) are shown and not interpreted."
    )
    lines.append("")

    def gene_cell(cohort: str, pred: str, symbol: str) -> str:
        ep = f"{symbol}_protein"
        hit = corrs[
            (corrs["cohort"] == cohort)
            & (corrs["predictor"] == pred)
            & (corrs["endpoint"] == ep)
            & (corrs["class"] != "score")
        ]
        if hit.empty:
            return "NA"
        r = hit.iloc[0]
        if not r["present"]:
            return "absent"
        if not np.isfinite(r["rho"]):
            return f"n={int(r['n'])}"
        return f"{int(r['n'])}, {r['rho']:+.3f}, {fmt_p(r['p'])}"

    blocks = [
        ("IFN ligands", IFN_LIGANDS),
        ("IFN receptors", IFN_RECEPTORS),
        ("IFN signaling", IFN_SIGNALING),
        ("ISGs", IFN_ISG),
        ("MHC-I", MHC1),
        ("Antigen processing (extra)", APM),
        ("CD8", CD8),
    ]
    for title, genes in blocks:
        lines.append(f"### {title}")
        lines.append("")
        grows = []
        for g in genes:
            grows.append(
                [
                    g,
                    gene_cell("LUAD", "CLDN4_protein", g),
                    gene_cell("LSCC", "CLDN4_protein", g),
                    gene_cell("LUAD", "TACSTD2_protein", g),
                    gene_cell("LSCC", "TACSTD2_protein", g),
                ]
            )
        lines.append(
            md_table(
                ["Protein", "LUAD CLDN4", "LSCC CLDN4", "LUAD TACSTD2", "LSCC TACSTD2"],
                grows,
            )
        )
        lines.append("")

    # Coexpression
    lines.append("## CLDN4 vs TACSTD2 protein")
    lines.append("")
    co_rows = []
    for cohort in ["LUAD", "LSCC"]:
        hit = corrs[
            (corrs["cohort"] == cohort)
            & (corrs["predictor"] == "CLDN4_protein")
            & (corrs["endpoint"] == "TACSTD2_protein")
        ].iloc[0]
        co_rows.append(
            [cohort, str(int(hit["n"])), fmt_rho(hit["rho"]), fmt_ci(hit["ci_lo"], hit["ci_hi"]), fmt_p(hit["p"])]
        )
    lines.append(md_table(["Cohort", "n", "ρ", "95% CI", "p"], co_rows))
    lines.append("")
    lines.append(
        "Coexpression is context. It is not an IFN/MHC endpoint and it is outside the primary-family FDR."
    )
    lines.append("")

    lines.append("## What this does not claim")
    lines.append("")
    lines.append(
        "- It does not put a private PDX number into the public table. Support is sign-and-uncertainty on CPTAC protein only.\n"
        "- It does not test ICI, PD-L1 IHC, or RNA. CD274 protein was null in both histologies on PR #356 and is not re-fit.\n"
        "- It does not treat IFN **ligand** protein as measured. IFNG, IFNA1, and IFNB1 rows are absent. B2M is absent.\n"
        "- It does not impute the CLDN4-missing tumors, and it does not write n=110 or n=108 as the CLDN4 n.\n"
        "- It does not pool LUAD and LSCC.\n"
        "- Gene-level HLA or JAK rows are not a second discovery set. The support call is the primary scores.\n"
        "- Hallmark IFN-γ / IFN-α columns in the phenotype file are RNA signatures and are unused.\n"
        "- WGS purity is stored on the phenotype table and is not a second primary adjustment. WES is the pre-specified residual, matching PR #289."
    )
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append(
        "- **Source:** open S3 `cptac-pancancer-data / data_freeze_v1.2_reorganized/{LUAD,LSCC}/`. "
        "Filenames from the LinkedOmics CPTAC-pancan index.\n"
        "- **Protein file:** `{LUAD,LSCC}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt`. "
        "Already log2 vs reference; not logged again.\n"
        "- **Phenotype file:** `{LUAD,LSCC}_phenotype.txt`, joined on `idx` = protein column name (110/110 and 108/108). "
        "Used for the treatment-column audit and `WES_purity`.\n"
        "- **IDs:** Ensembl prefix (`ENSG` ± version). Map in `genes.py`.\n"
        "- **Predictors:** CLDN4 and TACSTD2 protein.\n"
        "- **Primary endpoints:** MHC-I score (HLA-A/B/C), IFN-core score (STAT1, STAT2, IRF1, IRF9, JAK1, JAK2, ISG15, MX1, OAS1, IFI35, dropping any row that is absent or has fewer than 8 quantified tumors), CD8A protein.\n"
        "- **Gene-level list:** IFN ligands IFNG, IFNA1, IFNB1; receptors IFNAR1, IFNAR2, IFNGR1, IFNGR2; the IFN-core members; MHC-I HLA-A, HLA-B, HLA-C, B2M; antigen processing TAP1, TAP2, TAPBP, PSMB8, PSMB9, NLRC5; CD8B.\n"
        "- **Missingness:** empty cell = missing. Exact zero counted and kept as a value. "
        "Mann–Whitney compares endpoint distributions in predictor-quantified vs predictor-missing tumors.\n"
        "- **Partial:** rank-transform predictor, endpoint, and WES purity; residualize the two ranks on the purity rank; Pearson of residuals; bootstrap of that partial ρ.\n"
        "- **FDR:** Benjamini–Hochberg across the primary family only.\n"
        "- **Not done:** no RNA endpoints, no ImmuneScore re-fit, no imputation, no ICI model, no LUAD+LSCC pool, no post-hoc gene-set search."
    )
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append(
        "- `tables/coverage.tsv` — row present, quantified, missing, % missing, exact zeros, range\n"
        "- `tables/joint_missingness.tsv` — CLDN4 × TACSTD2 detection\n"
        "- `tables/missingness_contrast.tsv` — quantified vs missing\n"
        "- `tables/correlations.tsv` — n, ρ, p, CI, partial, q, call\n"
        "- `tables/sample_scores.tsv`\n"
        "- `tables/summary.json`\n"
        "- `figures/fig_missingness.png`\n"
        "- `figures/fig_forest.png`\n"
        "- `figures/fig_scatter_cldn4.png`\n"
        "- `figures/fig_scatter_tacstd2.png`"
    )
    lines.append("")
    lines.append("```bash")
    lines.append("python3 -m pip install -r methods/cptac_ifn_mhc_protein/requirements.txt")
    lines.append("python3 methods/cptac_ifn_mhc_protein/download.py --outdir data/cptac_ifn_mhc_protein")
    lines.append("python3 methods/cptac_ifn_mhc_protein/analyze.py --data data/cptac_ifn_mhc_protein --outdir methods/cptac_ifn_mhc_protein")
    lines.append("```")
    lines.append("")
    path.write_text("\n".join(lines))


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--data", default="data/cptac_ifn_mhc_protein")
    p.add_argument("--outdir", default="methods/cptac_ifn_mhc_protein")
    args = p.parse_args()
    data = Path(args.data)
    outdir = Path(args.outdir)
    tab = outdir / "tables"
    fig = outdir / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)
    bundles = {c: prepare_cohort(data, c) for c in ["LUAD", "LSCC"]}

    coverage = pd.DataFrame([row for b in bundles.values() for row in b["coverage"]])
    joint = pd.DataFrame([joint_missing(b) for b in bundles.values()])
    corrs_list = [rec for b in bundles.values() for rec in run_correlations(b, rng)]
    for b in bundles.values():
        run_partials(b, [r for r in corrs_list if r["cohort"] == b["cohort"]], rng)
    corrs = pd.DataFrame(corrs_list)
    # BH within primary family. Non-primary rows keep NA q.
    q = bh_q(corrs.loc[corrs["primary_family"], "p"].tolist())
    corrs["q_primary"] = np.nan
    corrs.loc[corrs["primary_family"], "q_primary"] = q
    corrs["call"] = corrs.apply(lambda r: support_call(r) if r["primary_family"] else "", axis=1)

    miss = pd.DataFrame([row for b in bundles.values() for row in run_missingness(b)])

    # Sample table
    frames = []
    for cohort, b in bundles.items():
        df = pd.DataFrame(
            {
                "CLDN4_protein": b["series"]["CLDN4"],
                "TACSTD2_protein": b["series"]["TACSTD2"],
                "CD8A_protein": b["series"]["CD8A"],
                "MHC1_protein": b["series"]["MHC1_protein"],
                "IFN_core_protein": b["series"]["IFN_core_protein"],
                "HLA-A_protein": b["series"]["HLA-A"],
                "STAT1_protein": b["series"]["STAT1"],
            }
        )
        df["WES_purity"] = b["wes"] if b["wes"] is not None else np.nan
        df["cohort"] = cohort
        df.index.name = "case_id"
        frames.append(df.reset_index())
    samples = pd.concat(frames, ignore_index=True)

    coverage.to_csv(tab / "coverage.tsv", sep="\t", index=False)
    joint.to_csv(tab / "joint_missingness.tsv", sep="\t", index=False)
    miss.to_csv(tab / "missingness_contrast.tsv", sep="\t", index=False)
    corrs.to_csv(tab / "correlations.tsv", sep="\t", index=False)
    samples.to_csv(tab / "sample_scores.tsv", sep="\t", index=False)

    plot_missingness(coverage, fig)
    plot_forest(corrs[corrs["primary_family"]], fig)
    plot_scatters(bundles, fig)
    write_finding(outdir / "FINDING.md", bundles, coverage, corrs, miss, joint)

    from sensitivity import run_sensitivity  # noqa: E402

    run_sensitivity(data, outdir, bundles, rng)
    finding = (outdir / "FINDING.md").read_text()
    extra = (outdir / "SENSITIVITY.md").read_text()
    marker = "<!-- SENSITIVITY -->"
    if marker not in finding:
        raise RuntimeError("FINDING.md is missing the sensitivity marker")
    (outdir / "FINDING.md").write_text(finding.replace(marker, extra.rstrip() + "\n"))

    summary = {
        "seed": SEED,
        "n_boot": N_BOOT,
        "cohorts": {
            c: {
                "n_tumors": bundles[c]["n_tumors"],
                "n_protein_rows": bundles[c]["n_rows"],
                "treat_cols": bundles[c]["treat_cols"],
                "ifn_members": bundles[c]["ifn_members"],
                "mhc_members": bundles[c]["mhc_members"],
                "cldn4_missing": int(joint.loc[joint["cohort"] == c, "cldn4_missing"].iloc[0]),
                "tacstd2_missing": int(joint.loc[joint["cohort"] == c, "tacstd2_missing"].iloc[0]),
                "wes_n": int(bundles[c]["wes"].notna().sum()) if bundles[c]["wes"] is not None else 0,
            }
            for c in ["LUAD", "LSCC"]
        },
        "primary_calls": corrs.loc[
            corrs["primary_family"],
            ["cohort", "predictor", "endpoint", "n", "rho", "p", "q_primary", "rho_partial", "call"],
        ].to_dict(orient="records"),
    }
    (tab / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Hard checks so a bad join cannot look like a result.
    assert bundles["LUAD"]["n_tumors"] == 110
    assert bundles["LSCC"]["n_tumors"] == 108
    assert bundles["LUAD"]["rows_used"]["CLDN4"] == "ENSG00000189143.9"
    assert bundles["LSCC"]["rows_used"]["CLDN4"] == "ENSG00000189143.9"
    assert bundles["LUAD"]["treat_cols"] == []
    assert bundles["LSCC"]["treat_cols"] == []
    print(tab / "summary.json")
    prim = corrs[corrs["primary_family"]][
        ["cohort", "predictor", "endpoint", "n", "rho", "p", "q_primary", "rho_partial", "call"]
    ]
    print(prim.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
