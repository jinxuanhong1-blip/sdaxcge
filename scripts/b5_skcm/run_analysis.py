#!/usr/bin/env python3
"""Run the pre-specified B5_SKCM analysis: CLDN4 vs ICI response in melanoma.

Implements results/w200/B5_SKCM/ANALYSIS_PLAN.md. Writes tables, figures, and a
machine-readable summary used by the results report. Does not decide the narrative:
every pre-specified test is computed and written regardless of significance.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.b5_skcm.common import (  # noqa: E402
    AYERS_IFNG_6,
    EPITHELIAL_CONTROLS,
    GENE_OF_INTEREST,
    KERATINOCYTE_SCORE_GENES,
    MELANOCYTE_GENES,
    OUT,
    POSITIVE_CONTROLS,
    PROC,
    load_hugo,
    load_mgh,
    load_riaz,
    load_tcga_skcm,
    log2t,
)

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# --------------------------------------------------------------------------------------------
# Sample selection (plan §2)
# --------------------------------------------------------------------------------------------


def _patient_mean(expr: pd.DataFrame, clin: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Average log2 expression across replicate biopsies of the same patient; keep first clin row."""
    patients = clin["patient"].astype(str)
    if patients.nunique() == len(clin):
        return expr, clin
    grouped = expr.T.groupby(patients.values).mean().T
    keep = clin.reset_index().groupby("patient", as_index=False).first().set_index("patient")
    keep.index.name = "sample"
    # reindex expression columns to the patient ids used as the new sample index
    grouped = grouped.reindex(columns=keep.index)
    return grouped, keep


def select_primary(expr: pd.DataFrame, clin: pd.DataFrame, *, first_episode_only: bool = True,
                   response_col: str = "response") -> tuple[pd.DataFrame, pd.DataFrame]:
    mask = clin["pretreatment"].fillna(False) & clin[response_col].isin(["R", "NR"])
    if first_episode_only and "episode_order" in clin.columns:
        mask = mask & (clin["episode_order"].fillna(1) == 1)
    sub_c = clin.loc[mask].copy()
    sub_e = expr.loc[:, sub_c.index]
    return _patient_mean(sub_e, sub_c)


def zscore(s: pd.Series) -> pd.Series:
    sd = s.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean()) / sd


def gene_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns)
    z = expr.loc[present].apply(zscore, axis=1)
    return z.mean(axis=0)


# --------------------------------------------------------------------------------------------
# Association tests
# --------------------------------------------------------------------------------------------


def mwu_auc(x_r: np.ndarray, x_nr: np.ndarray) -> dict:
    x_r = np.asarray(x_r, dtype=float)
    x_nr = np.asarray(x_nr, dtype=float)
    n_r, n_nr = len(x_r), len(x_nr)
    if n_r < 2 or n_nr < 2:
        return {"n_r": n_r, "n_nr": n_nr, "mwu_U": np.nan, "mwu_p": np.nan,
                "auc": np.nan, "rank_biserial": np.nan, "median_r": np.nan,
                "median_nr": np.nan, "median_diff": np.nan}
    res = stats.mannwhitneyu(x_r, x_nr, alternative="two-sided", method="auto")
    auc = float(res.statistic) / (n_r * n_nr)  # P(R > NR) with ties 0.5
    return {
        "n_r": n_r,
        "n_nr": n_nr,
        "mwu_U": float(res.statistic),
        "mwu_p": float(res.pvalue),
        "auc": auc,
        "rank_biserial": 2 * auc - 1,
        "median_r": float(np.median(x_r)),
        "median_nr": float(np.median(x_nr)),
        "median_diff": float(np.median(x_r) - np.median(x_nr)),
    }


def logit_or_per_sd(z: np.ndarray, y: np.ndarray) -> dict:
    """Unadjusted logistic regression of binary y on a z-scored predictor. Returns OR per 1 SD."""
    z = np.asarray(z, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(z) & np.isfinite(y)
    z, y = z[ok], y[ok]
    out = {"n": int(len(y)), "n_r": int(y.sum()), "n_nr": int((1 - y).sum()),
           "logor": np.nan, "se": np.nan, "or": np.nan, "or_lo": np.nan, "or_hi": np.nan,
           "logit_p": np.nan, "converged": False, "note": ""}
    if len(y) < 8 or y.sum() < 2 or (1 - y).sum() < 2 or np.nanstd(z) == 0:
        out["note"] = "too_few_or_no_variance"
        return out
    X = np.column_stack([np.ones(len(z)), z])
    beta = np.zeros(2)
    converged = False
    for _ in range(40):
        eta = X @ beta
        eta = np.clip(eta, -20, 20)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = p * (1.0 - p)
        if w.sum() < 1e-8:
            break
        # IRLS
        z_wls = eta + (y - p) / np.clip(w, 1e-8, None)
        XtW = X.T * w
        try:
            beta_new = np.linalg.solve(XtW @ X, XtW @ z_wls)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(beta_new - beta)) < 1e-8:
            beta = beta_new
            converged = True
            break
        beta = beta_new
    eta = np.clip(X @ beta, -20, 20)
    p = 1.0 / (1.0 + np.exp(-eta))
    w = p * (1.0 - p)
    try:
        cov = np.linalg.inv(X.T @ (w[:, None] * X))
        se = float(np.sqrt(cov[1, 1]))
    except np.linalg.LinAlgError:
        out["note"] = "hessian_singular"
        return out
    logor = float(beta[1])
    # Wald
    zstat = logor / se if se > 0 else np.nan
    pval = float(2 * stats.norm.sf(abs(zstat))) if np.isfinite(zstat) else np.nan
    out.update({
        "logor": logor,
        "se": se,
        "or": float(np.exp(logor)),
        "or_lo": float(np.exp(logor - 1.96 * se)),
        "or_hi": float(np.exp(logor + 1.96 * se)),
        "logit_p": pval,
        "converged": bool(converged),
    })
    return out


def logit_or_adjusted(z: np.ndarray, y: np.ndarray, cov: pd.DataFrame) -> dict:
    """Logistic y ~ z + covariates. Categorical columns are dummy-coded; drop-first."""
    z = pd.Series(z, index=cov.index, dtype=float)
    y = pd.Series(y, index=cov.index, dtype=float)
    X_parts = [pd.DataFrame({"z": z})]
    for c in cov.columns:
        s = cov[c]
        if s.dtype == object or str(s.dtype) == "category" or s.nunique(dropna=True) <= 8:
            dummies = pd.get_dummies(s.astype("string"), prefix=c, drop_first=True, dummy_na=False)
            X_parts.append(dummies)
        else:
            X_parts.append(pd.DataFrame({c: pd.to_numeric(s, errors="coerce")}))
    Xdf = pd.concat(X_parts, axis=1).apply(pd.to_numeric, errors="coerce")
    ok = Xdf.notna().all(axis=1) & y.notna()
    Xdf, y = Xdf.loc[ok], y.loc[ok]
    out = {"n": int(len(y)), "n_r": int(y.sum()), "logor": np.nan, "se": np.nan,
           "or": np.nan, "or_lo": np.nan, "or_hi": np.nan, "logit_p": np.nan,
           "converged": False, "covariates": ",".join(cov.columns)}
    if len(y) < 10 or y.sum() < 2 or (1 - y).sum() < 2:
        return out
    X = np.column_stack([np.ones(len(Xdf)), Xdf.to_numpy(dtype=float)])
    beta = np.zeros(X.shape[1])
    converged = False
    for _ in range(50):
        eta = np.clip(X @ beta, -20, 20)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = p * (1.0 - p)
        if w.sum() < 1e-8:
            break
        z_wls = eta + (y.to_numpy() - p) / np.clip(w, 1e-8, None)
        XtW = X.T * w
        try:
            beta_new = np.linalg.solve(XtW @ X, XtW @ z_wls)
        except np.linalg.LinAlgError:
            return out
        if np.max(np.abs(beta_new - beta)) < 1e-7:
            beta = beta_new
            converged = True
            break
        beta = beta_new
    eta = np.clip(X @ beta, -20, 20)
    p = 1.0 / (1.0 + np.exp(-eta))
    w = p * (1.0 - p)
    try:
        covm = np.linalg.inv(X.T @ (w[:, None] * X))
        se = float(np.sqrt(covm[1, 1]))
    except np.linalg.LinAlgError:
        return out
    logor = float(beta[1])
    zstat = logor / se if se > 0 else np.nan
    out.update({
        "logor": logor, "se": se, "or": float(np.exp(logor)),
        "or_lo": float(np.exp(logor - 1.96 * se)),
        "or_hi": float(np.exp(logor + 1.96 * se)),
        "logit_p": float(2 * stats.norm.sf(abs(zstat))) if np.isfinite(zstat) else np.nan,
        "converged": bool(converged),
    })
    return out


def ivw_meta(rows: list[dict]) -> dict:
    """Inverse-variance meta-analysis of log-ORs. Fixed-effect + DerSimonian–Laird RE."""
    use = [r for r in rows if np.isfinite(r.get("logor", np.nan)) and np.isfinite(r.get("se", np.nan))
           and r["se"] > 0]
    out = {"k": len(use), "logor_fe": np.nan, "se_fe": np.nan, "or_fe": np.nan,
           "or_fe_lo": np.nan, "or_fe_hi": np.nan, "p_fe": np.nan,
           "logor_re": np.nan, "se_re": np.nan, "or_re": np.nan, "or_re_lo": np.nan,
           "or_re_hi": np.nan, "p_re": np.nan, "Q": np.nan, "Q_p": np.nan, "I2": np.nan,
           "tau2": np.nan, "direction_pos": 0, "direction_neg": 0}
    if len(use) < 1:
        return out
    b = np.array([r["logor"] for r in use])
    se = np.array([r["se"] for r in use])
    w = 1.0 / se ** 2
    logor_fe = float(np.sum(w * b) / np.sum(w))
    se_fe = float(np.sqrt(1.0 / np.sum(w)))
    Q = float(np.sum(w * (b - logor_fe) ** 2))
    df = len(use) - 1
    Q_p = float(stats.chi2.sf(Q, df)) if df > 0 else np.nan
    I2 = float(max(0.0, (Q - df) / Q * 100.0)) if Q > 0 and df > 0 else 0.0
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w) if df > 0 else np.nan
    tau2 = float(max(0.0, (Q - df) / c)) if df > 0 and np.isfinite(c) and c > 0 else 0.0
    w_re = 1.0 / (se ** 2 + tau2)
    logor_re = float(np.sum(w_re * b) / np.sum(w_re))
    se_re = float(np.sqrt(1.0 / np.sum(w_re)))

    def pack(logor, se_):
        z = logor / se_
        return {
            "logor": logor, "se": se_, "or": float(np.exp(logor)),
            "or_lo": float(np.exp(logor - 1.96 * se_)),
            "or_hi": float(np.exp(logor + 1.96 * se_)),
            "p": float(2 * stats.norm.sf(abs(z))),
        }

    fe, re = pack(logor_fe, se_fe), pack(logor_re, se_re)
    out.update({
        "k": len(use),
        "logor_fe": fe["logor"], "se_fe": fe["se"], "or_fe": fe["or"],
        "or_fe_lo": fe["or_lo"], "or_fe_hi": fe["or_hi"], "p_fe": fe["p"],
        "logor_re": re["logor"], "se_re": re["se"], "or_re": re["or"],
        "or_re_lo": re["or_lo"], "or_re_hi": re["or_hi"], "p_re": re["p"],
        "Q": Q, "Q_p": Q_p, "I2": I2, "tau2": tau2,
        "direction_pos": int((b > 0).sum()),
        "direction_neg": int((b < 0).sum()),
    })
    return out


def test_gene(expr: pd.DataFrame, clin: pd.DataFrame, gene: str) -> dict:
    if gene not in expr.index:
        return {"gene": gene, "present": False}
    x = expr.loc[gene]
    y = (clin["response"] == "R").astype(int)
    z = zscore(x)
    mwu = mwu_auc(x[y == 1].to_numpy(), x[y == 0].to_numpy())
    logit = logit_or_per_sd(z.to_numpy(), y.to_numpy())
    out = {"gene": gene, "present": True, **mwu, **{f"logit_{k}" if not k.startswith("logit")
                                                    and k not in mwu else k: v
                                                    for k, v in logit.items()}}
    # flatten without clobbering
    out = {"gene": gene, "present": True}
    out.update(mwu)
    out.update({
        "logor": logit["logor"], "se": logit["se"], "or": logit["or"],
        "or_lo": logit["or_lo"], "or_hi": logit["or_hi"],
        "logit_p": logit["logit_p"], "logit_converged": logit["converged"],
        "logit_note": logit.get("note", ""),
    })
    return out


# --------------------------------------------------------------------------------------------
# Survival (Hugo only)
# --------------------------------------------------------------------------------------------


def cox_and_logrank(time, event, z, high) -> dict:
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test

    df = pd.DataFrame({"T": time, "E": event, "z": z, "high": high}).dropna()
    out = {"n": int(len(df)), "n_event": int(df["E"].sum()),
           "hr_cont": np.nan, "hr_cont_lo": np.nan, "hr_cont_hi": np.nan, "hr_cont_p": np.nan,
           "hr_split": np.nan, "hr_split_lo": np.nan, "hr_split_hi": np.nan, "hr_split_p": np.nan,
           "logrank_p": np.nan}
    if len(df) < 8 or df["E"].sum() < 3:
        return out
    cph = CoxPHFitter()
    try:
        cph.fit(df[["T", "E", "z"]], duration_col="T", event_col="E")
        s = cph.summary.loc["z"]
        out.update({
            "hr_cont": float(s["exp(coef)"]),
            "hr_cont_lo": float(s["exp(coef) lower 95%"]),
            "hr_cont_hi": float(s["exp(coef) upper 95%"]),
            "hr_cont_p": float(s["p"]),
        })
    except Exception as e:
        out["cox_cont_error"] = str(e)
    try:
        cph.fit(df[["T", "E", "high"]], duration_col="T", event_col="E")
        s = cph.summary.loc["high"]
        out.update({
            "hr_split": float(s["exp(coef)"]),
            "hr_split_lo": float(s["exp(coef) lower 95%"]),
            "hr_split_hi": float(s["exp(coef) upper 95%"]),
            "hr_split_p": float(s["p"]),
        })
        lr = logrank_test(df.loc[df["high"] == 1, "T"], df.loc[df["high"] == 0, "T"],
                          df.loc[df["high"] == 1, "E"], df.loc[df["high"] == 0, "E"])
        out["logrank_p"] = float(lr.p_value)
    except Exception as e:
        out["cox_split_error"] = str(e)
    return out


# --------------------------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------------------------


def _save(fig, name: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / "figures" / name, dpi=160, bbox_inches="tight")
    plt.close(fig)


def fig_boxplots(per_cohort: dict, gene: str = GENE_OF_INTEREST) -> None:
    names = list(per_cohort)
    fig, axes = plt.subplots(1, len(names), figsize=(3.4 * len(names), 4.2), sharey=False)
    if len(names) == 1:
        axes = [axes]
    for ax, name in zip(axes, names):
        expr, clin = per_cohort[name]
        x = expr.loc[gene]
        groups = [("NR", x[clin["response"] == "NR"]), ("R", x[clin["response"] == "R"])]
        data = [g[1].to_numpy() for g in groups]
        bp = ax.boxplot(data, labels=["NR", "R"], widths=0.55, patch_artist=True,
                        medianprops=dict(color="black", lw=1.4))
        bp["boxes"][0].set_facecolor("#cfcfcf")
        bp["boxes"][1].set_facecolor("#6baed6")
        rng = np.random.default_rng(0)
        for i, d in enumerate(data, start=1):
            ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(d)), d, s=14, color="#333", zorder=3)
        ax.set_title(f"{name}\nn={len(x)} ({(clin['response']=='R').sum()} R)")
        ax.set_ylabel(f"log2({gene} + 0.1)" if ax is axes[0] else "")
    fig.suptitle(f"{gene} in pre-treatment biopsies, by ICI response", y=1.02)
    _save(fig, "cldn4_boxplot_by_response.png")


def fig_forest(primary_rows: pd.DataFrame, meta: dict) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    labels = list(primary_rows["cohort"]) + ["Fixed-effect meta", "Random-effects meta"]
    ors = list(primary_rows["or"]) + [meta["or_fe"], meta["or_re"]]
    lo = list(primary_rows["or_lo"]) + [meta["or_fe_lo"], meta["or_re_lo"]]
    hi = list(primary_rows["or_hi"]) + [meta["or_fe_hi"], meta["or_re_hi"]]
    y = np.arange(len(labels))[::-1]
    ax.axvline(1.0, color="#888", lw=1, ls="--")
    for i, (o, a, b) in enumerate(zip(ors, lo, hi)):
        color = "#222" if i < len(primary_rows) else "#b2182b"
        if not np.isfinite(o):
            continue
        ax.plot([a, b], [y[i], y[i]], color=color, lw=1.6)
        ax.plot(o, y[i], "o", color=color, ms=6)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio for response per 1 SD log2 CLDN4 (95% CI)")
    ax.set_title("Primary analysis: CLDN4 vs ICI response")
    _save(fig, "cldn4_forest.png")


def fig_genomewide(gw: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    p = gw["p_fe"].clip(lower=1e-12)
    ax.scatter(gw["logor_fe"], -np.log10(p), s=6, c="#9e9e9e", alpha=0.45, linewidths=0)
    hit = gw[gw["gene"] == GENE_OF_INTEREST]
    if len(hit):
        ax.scatter(hit["logor_fe"], -np.log10(hit["p_fe"].clip(lower=1e-12)),
                   s=50, c="#b2182b", zorder=3, label="CLDN4")
        ax.annotate("CLDN4", (hit["logor_fe"].iloc[0], -np.log10(hit["p_fe"].iloc[0])),
                    textcoords="offset points", xytext=(6, 6), color="#b2182b")
    for g in POSITIVE_CONTROLS:
        row = gw[gw["gene"] == g]
        if len(row):
            ax.scatter(row["logor_fe"], -np.log10(row["p_fe"].clip(lower=1e-12)),
                       s=22, c="#2166ac", zorder=2)
    ax.axhline(-np.log10(0.05), color="#888", ls="--", lw=0.8)
    ax.set_xlabel("Fixed-effect meta log-odds ratio (higher = more expression in responders)")
    ax.set_ylabel("-log10 p (fixed-effect meta)")
    ax.set_title("Genome-wide meta-analysis; CLDN4 in red, ICI positive controls in blue")
    _save(fig, "genomewide_volcano.png")


def fig_epithelial_corr(corr_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    mat = corr_df.set_index("gene")
    im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=30, ha="right")
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index)
    fig.colorbar(im, ax=ax, shrink=0.8, label="Spearman rho vs CLDN4")
    ax.set_title("CLDN4 vs epithelial / melanocyte genes (pre-treatment, primary set)")
    _save(fig, "cldn4_epithelial_correlation.png")


def fig_tcga(x: pd.Series) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.hist(x.dropna(), bins=40, color="#6baed6", edgecolor="white")
    ax.axvline(x.median(), color="#b2182b", ls="--", label=f"median={x.median():.2f}")
    ax.set_xlabel("log2(CLDN4 + 0.1)  TCGA-SKCM tumors")
    ax.set_ylabel("Tumors")
    ax.legend()
    ax.set_title("CLDN4 is low in most TCGA melanomas")
    _save(fig, "tcga_skcm_cldn4_hist.png")


def fig_km(time, event, high, title: str) -> None:
    from lifelines import KaplanMeierFitter
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    kmf = KaplanMeierFitter()
    df = pd.DataFrame({"T": time, "E": event, "high": high}).dropna()
    for lab, val, color in [("CLDN4 low", 0, "#4d4d4d"), ("CLDN4 high", 1, "#b2182b")]:
        sub = df[df["high"] == val]
        if len(sub) == 0:
            continue
        kmf.fit(sub["T"], sub["E"], label=f"{lab} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, color=color, ci_show=True)
    ax.set_xlabel("Overall survival (days)")
    ax.set_ylabel("Survival probability")
    ax.set_title(title)
    _save(fig, "hugo_os_km.png")


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    PROC.mkdir(parents=True, exist_ok=True)

    print("Loading cohorts...", flush=True)
    loaders = {
        "Hugo_GSE78220": load_hugo,
        "Riaz_GSE91061": load_riaz,
        "MGH_GSE115821": load_mgh,
    }
    raw = {}
    load_notes = {}
    for name, fn in loaders.items():
        expr, clin, notes = fn()
        raw[name] = (expr, clin)
        load_notes[name] = notes
        print(f"  {name}: {expr.shape[0]} genes x {expr.shape[1]} samples; notes={notes}", flush=True)

    # ---- primary analysis set ----
    primary = {}
    inventory_rows = []
    for name, (expr, clin) in raw.items():
        e, c = select_primary(expr, clin, first_episode_only=True, response_col="response")
        primary[name] = (e, c)
        inventory_rows.append({
            "cohort": name,
            "n_samples_downloaded": int(clin.shape[0]),
            "n_pretreatment": int(clin["pretreatment"].fillna(False).sum()),
            "n_primary_patients": int(c.shape[0]),
            "n_R": int((c["response"] == "R").sum()),
            "n_NR": int((c["response"] == "NR").sum()),
            "unit": load_notes[name]["unit"],
            "CLDN4_present": GENE_OF_INTEREST in e.index,
        })
    inventory = pd.DataFrame(inventory_rows)
    inventory.to_csv(OUT / "tables" / "cohort_inventory.tsv", sep="\t", index=False)

    # ---- detectability (plan §3) ----
    det_rows = []
    for name, (e, c) in primary.items():
        if GENE_OF_INTEREST not in e.index:
            det_rows.append({"cohort": name, "present": False})
            continue
        x = e.loc[GENE_OF_INTEREST]
        # back-transform to linear units for "detectable" = original > 0
        lin = np.power(2.0, x) - 0.1
        det_rows.append({
            "cohort": name,
            "present": True,
            "n": int(len(x)),
            "n_detectable_linear_gt0": int((lin > 0).sum()),
            "frac_detectable": float((lin > 0).mean()),
            "median_log2": float(x.median()),
            "q25_log2": float(x.quantile(0.25)),
            "q75_log2": float(x.quantile(0.75)),
            "min_log2": float(x.min()),
            "max_log2": float(x.max()),
            "frac_below_log2_0": float((x < 0).mean()),  # < 1 linear unit
            "low_information_flag": bool((lin > 0).mean() < 0.5 or x.median() < 0),
        })
    detect = pd.DataFrame(det_rows)
    detect.to_csv(OUT / "tables" / "cldn4_detectability.tsv", sep="\t", index=False)
    print("Detectability:\n", detect.to_string(index=False), flush=True)

    # ---- primary per-cohort tests ----
    primary_rows = []
    for name, (e, c) in primary.items():
        r = test_gene(e, c, GENE_OF_INTEREST)
        r["cohort"] = name
        primary_rows.append(r)
    primary_df = pd.DataFrame(primary_rows)
    primary_df.to_csv(OUT / "tables" / "cldn4_primary_per_cohort.tsv", sep="\t", index=False)
    meta = ivw_meta(primary_rows)
    pd.DataFrame([meta]).to_csv(OUT / "tables" / "cldn4_primary_meta.tsv", sep="\t", index=False)
    print("Primary per-cohort:\n", primary_df.to_string(index=False), flush=True)
    print("Primary meta:", meta, flush=True)

    # ---- keratinocyte / melanocyte scores and correlations ----
    corr_rows = []
    score_corrs = []
    for name, (e, c) in primary.items():
        if GENE_OF_INTEREST not in e.index:
            continue
        x = e.loc[GENE_OF_INTEREST]
        krt = gene_score(e, KERATINOCYTE_SCORE_GENES)
        mel = gene_score(e, MELANOCYTE_GENES)
        if krt.notna().sum() >= 5:
            rho, p = stats.spearmanr(x, krt, nan_policy="omit")
            score_corrs.append({"cohort": name, "score": "keratinocyte", "rho": float(rho), "p": float(p)})
        if mel.notna().sum() >= 5:
            rho, p = stats.spearmanr(x, mel, nan_policy="omit")
            score_corrs.append({"cohort": name, "score": "melanocyte", "rho": float(rho), "p": float(p)})
        for g in EPITHELIAL_CONTROLS + MELANOCYTE_GENES:
            if g not in e.index:
                corr_rows.append({"cohort": name, "gene": g, "rho": np.nan, "p": np.nan})
                continue
            rho, p = stats.spearmanr(x, e.loc[g], nan_policy="omit")
            corr_rows.append({"cohort": name, "gene": g, "rho": float(rho), "p": float(p)})
    corr_long = pd.DataFrame(corr_rows)
    corr_long.to_csv(OUT / "tables" / "cldn4_gene_correlations.tsv", sep="\t", index=False)
    corr_wide = corr_long.pivot(index="gene", columns="cohort", values="rho").reset_index()
    corr_wide.to_csv(OUT / "tables" / "cldn4_gene_correlations_wide.tsv", sep="\t", index=False)
    pd.DataFrame(score_corrs).to_csv(OUT / "tables" / "cldn4_score_correlations.tsv", sep="\t", index=False)

    # ---- adjusted models (plan §5) ----
    adj_rows = []
    for name, (e, c) in primary.items():
        if GENE_OF_INTEREST not in e.index:
            continue
        z = zscore(e.loc[GENE_OF_INTEREST])
        y = (c["response"] == "R").astype(int)
        krt = gene_score(e, KERATINOCYTE_SCORE_GENES)
        cov = pd.DataFrame({"krt": krt}, index=c.index)
        r = logit_or_adjusted(z.to_numpy(), y.to_numpy(), cov)
        r.update({"cohort": name, "model": "zCLDN4 + keratinocyte_score"})
        adj_rows.append(r)
        extra_cols = []
        if name == "Hugo_GSE78220":
            extra_cols = [col for col in ["library", "site"] if col in c.columns]
        if name == "MGH_GSE115821":
            extra_cols = [col for col in ["batch", "therapy"] if col in c.columns]
        if extra_cols:
            cov2 = cov.copy()
            for col in extra_cols:
                cov2[col] = c[col]
            r2 = logit_or_adjusted(z.to_numpy(), y.to_numpy(), cov2)
            r2.update({"cohort": name, "model": "zCLDN4 + keratinocyte_score + " + "+".join(extra_cols)})
            adj_rows.append(r2)
    adj_df = pd.DataFrame(adj_rows)
    adj_df.to_csv(OUT / "tables" / "cldn4_adjusted_logit.tsv", sep="\t", index=False)
    meta_adj = ivw_meta([r for r in adj_rows if r["model"] == "zCLDN4 + keratinocyte_score"])
    pd.DataFrame([meta_adj]).to_csv(OUT / "tables" / "cldn4_adjusted_krt_meta.tsv", sep="\t", index=False)

    # ---- positive / negative controls ----
    ctrl_rows = []
    for name, (e, c) in primary.items():
        for g in POSITIVE_CONTROLS + EPITHELIAL_CONTROLS + ["CLDN4"]:
            r = test_gene(e, c, g)
            r["cohort"] = name
            r["class"] = ("positive_control" if g in POSITIVE_CONTROLS
                          else "epithelial_control" if g in EPITHELIAL_CONTROLS
                          else "gene_of_interest")
            ctrl_rows.append(r)
        # Ayers 6-gene score
        sc = gene_score(e, AYERS_IFNG_6)
        y = (c["response"] == "R").astype(int)
        mwu = mwu_auc(sc[y == 1].to_numpy(), sc[y == 0].to_numpy())
        logit = logit_or_per_sd(zscore(sc).to_numpy(), y.to_numpy())
        ctrl_rows.append({
            "gene": "Ayers_IFNG_6", "present": True, "cohort": name, "class": "positive_control_score",
            **mwu, "logor": logit["logor"], "se": logit["se"], "or": logit["or"],
            "or_lo": logit["or_lo"], "or_hi": logit["or_hi"], "logit_p": logit["logit_p"],
        })
    ctrl_df = pd.DataFrame(ctrl_rows)
    ctrl_df.to_csv(OUT / "tables" / "control_genes_per_cohort.tsv", sep="\t", index=False)
    # meta per control gene
    ctrl_meta = []
    for g, sub in ctrl_df.groupby("gene"):
        m = ivw_meta(sub.to_dict("records"))
        m["gene"] = g
        m["class"] = sub["class"].iloc[0]
        ctrl_meta.append(m)
    ctrl_meta_df = pd.DataFrame(ctrl_meta).sort_values("p_fe")
    ctrl_meta_df.to_csv(OUT / "tables" / "control_genes_meta.tsv", sep="\t", index=False)
    print("Control meta (top):\n", ctrl_meta_df.head(15).to_string(index=False), flush=True)

    # ---- genome-wide sweep (plan §5) ----
    print("Genome-wide sweep...", flush=True)
    common_genes = set.intersection(*[set(e.index) for e, _ in primary.values()])
    print(f"  genes in all 3 cohorts: {len(common_genes)}", flush=True)
    gw_rows = []
    for i, g in enumerate(sorted(common_genes)):
        cohort_stats = []
        for name, (e, c) in primary.items():
            cohort_stats.append(test_gene(e, c, g))
        m = ivw_meta(cohort_stats)
        m["gene"] = g
        # require finite FE p
        gw_rows.append(m)
        if (i + 1) % 2000 == 0:
            print(f"  ... {i+1}/{len(common_genes)}", flush=True)
    gw = pd.DataFrame(gw_rows)
    gw = gw[np.isfinite(gw["p_fe"])].copy()
    gw["p_fe_rank"] = gw["p_fe"].rank(method="min")
    gw["p_fe_percentile"] = 100.0 * gw["p_fe_rank"] / len(gw)
    gw["bh_fdr"] = multipletests(gw["p_fe"].to_numpy(), method="fdr_bh")[1]
    gw = gw.sort_values(["p_fe", "gene"])
    gw.to_csv(OUT / "tables" / "genomewide_meta.tsv", sep="\t", index=False)
    cldn4_gw = gw[gw["gene"] == GENE_OF_INTEREST]
    print("CLDN4 genome-wide rank:\n", cldn4_gw.to_string(index=False), flush=True)
    n_sig_nom = int((gw["p_fe"] < 0.05).sum())
    n_sig_fdr = int((gw["bh_fdr"] < 0.05).sum())

    # ---- survival (Hugo, exploratory) ----
    e_h, c_h = primary["Hugo_GSE78220"]
    surv = {}
    if GENE_OF_INTEREST in e_h.index and {"os_days", "os_event"}.issubset(c_h.columns):
        z = zscore(e_h.loc[GENE_OF_INTEREST])
        high = (e_h.loc[GENE_OF_INTEREST] >= e_h.loc[GENE_OF_INTEREST].median()).astype(int)
        surv = cox_and_logrank(c_h["os_days"], c_h["os_event"], z, high)
        pd.DataFrame([surv]).to_csv(OUT / "tables" / "hugo_os_cldn4.tsv", sep="\t", index=False)
        fig_km(c_h["os_days"], c_h["os_event"], high,
               "Hugo GSE78220 OS by median-split CLDN4 (exploratory)")
        print("Hugo OS:", surv, flush=True)

    # ---- sensitivities ----
    sens = []

    # S1: Riaz SD as NR
    e_r, c_r_full = raw["Riaz_GSE91061"]
    e1, c1 = select_primary(e_r, c_r_full, first_episode_only=True, response_col="response_sd_as_nr")
    # temporarily put response
    c1 = c1.copy()
    c1["response"] = c1["response_sd_as_nr"]
    s1_riaz = test_gene(e1, c1, GENE_OF_INTEREST)
    s1_rows = []
    for name, (e, c) in primary.items():
        if name == "Riaz_GSE91061":
            s1_rows.append(s1_riaz)
        else:
            s1_rows.append(test_gene(e, c, GENE_OF_INTEREST))
    s1_meta = ivw_meta(s1_rows)
    s1_meta["analysis"] = "S1_Riaz_SD_as_NR"
    s1_riaz_out = {"analysis": "S1_Riaz_SD_as_NR_cohort", **s1_riaz,
                   "n_R": int((c1["response"] == "R").sum()),
                   "n_NR": int((c1["response"] == "NR").sum())}
    sens.append({**s1_meta})
    pd.DataFrame([s1_riaz_out]).to_csv(OUT / "tables" / "sensitivity_S1_riaz_cohort.tsv", sep="\t", index=False)

    # S2: median-dichotomized CLDN4, logistic high vs low
    s2_rows = []
    for name, (e, c) in primary.items():
        if GENE_OF_INTEREST not in e.index:
            continue
        high = (e.loc[GENE_OF_INTEREST] >= e.loc[GENE_OF_INTEREST].median()).astype(float)
        y = (c["response"] == "R").astype(float)
        # 2x2 fisher + logit on the binary predictor (not z-scored; OR is high vs low)
        tab = pd.crosstab(high, y)
        if tab.shape == (2, 2):
            or_f, p_f = stats.fisher_exact(tab.to_numpy())
        else:
            or_f, p_f = np.nan, np.nan
        logit = logit_or_per_sd(zscore(high).to_numpy(), y.to_numpy())
        s2_rows.append({"cohort": name, "fisher_or": float(or_f) if np.isfinite(or_f) else np.nan,
                        "fisher_p": float(p_f) if np.isfinite(p_f) else np.nan, **logit})
    s2_meta = ivw_meta(s2_rows)
    s2_meta["analysis"] = "S2_median_split"
    pd.DataFrame(s2_rows).to_csv(OUT / "tables" / "sensitivity_S2_per_cohort.tsv", sep="\t", index=False)
    sens.append(s2_meta)

    # S3: Riaz from raw counts, UQ-normalized log-CPM
    raw_counts = pd.read_csv(
        Path(__file__).resolve().parents[2] / "data" / "raw" /
        "GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz", index_col=0
    )
    raw_counts.index = raw_counts.index.astype(str)
    from scripts.b5_skcm.common import entrez_to_symbol, collapse_duplicate_symbols
    sym = entrez_to_symbol()
    mapped = raw_counts.index.map(sym)
    raw_counts = raw_counts[~pd.isna(mapped)]
    raw_counts.index = pd.Index([s for s in mapped if not pd.isna(s)])
    raw_counts = collapse_duplicate_symbols(raw_counts)
    # upper-quartile normalization
    uq = raw_counts.apply(lambda s: np.quantile(s[s > 0], 0.75) if (s > 0).any() else np.nan)
    cpm = raw_counts.div(uq, axis=1) * uq.median()
    # library-size CPM-like: already UQ scaled counts; convert to per-million of UQ-scaled lib
    lib = cpm.sum(axis=0)
    lcpm = log2t(cpm.div(lib, axis=1) * 1e6)
    e3, c3 = select_primary(lcpm, raw["Riaz_GSE91061"][1])
    s3 = test_gene(e3, c3, GENE_OF_INTEREST)
    s3_rows = []
    for name, (e, c) in primary.items():
        if name == "Riaz_GSE91061":
            s3_rows.append(s3)
        else:
            s3_rows.append(test_gene(e, c, GENE_OF_INTEREST))
    s3_meta = ivw_meta(s3_rows)
    s3_meta["analysis"] = "S3_Riaz_UQ_logCPM"
    pd.DataFrame([{"analysis": "S3_Riaz_UQ_logCPM_cohort", **s3}]).to_csv(
        OUT / "tables" / "sensitivity_S3_riaz_cohort.tsv", sep="\t", index=False)
    sens.append(s3_meta)

    # S4: MGH all baseline episodes (not just first)
    e_m, c_m = raw["MGH_GSE115821"]
    e4, c4 = select_primary(e_m, c_m, first_episode_only=False)
    s4 = test_gene(e4, c4, GENE_OF_INTEREST)
    s4_rows = []
    for name, (e, c) in primary.items():
        if name == "MGH_GSE115821":
            s4_rows.append(s4)
        else:
            s4_rows.append(test_gene(e, c, GENE_OF_INTEREST))
    s4_meta = ivw_meta(s4_rows)
    s4_meta["analysis"] = "S4_MGH_all_baseline_episodes"
    pd.DataFrame([{"analysis": "S4_MGH_all_episodes_cohort", **s4,
                   "n_patients": int(c4.shape[0])}]).to_csv(
        OUT / "tables" / "sensitivity_S4_mgh_cohort.tsv", sep="\t", index=False)
    sens.append(s4_meta)

    # S5: one-stage IPD logistic with cohort fixed effect
    frames = []
    for name, (e, c) in primary.items():
        if GENE_OF_INTEREST not in e.index:
            continue
        frames.append(pd.DataFrame({
            "z": zscore(e.loc[GENE_OF_INTEREST]),
            "y": (c["response"] == "R").astype(int),
            "cohort": name,
        }))
    ipd = pd.concat(frames, axis=0)
    dummies = pd.get_dummies(ipd["cohort"], drop_first=True)
    cov = dummies.copy()
    s5 = logit_or_adjusted(ipd["z"].to_numpy(), ipd["y"].to_numpy(), cov)
    s5["analysis"] = "S5_IPD_cohort_FE"
    s5["n_total"] = int(len(ipd))
    s5["n_R"] = int(ipd["y"].sum())
    sens.append({
        "analysis": "S5_IPD_cohort_FE",
        "or_fe": s5["or"], "or_fe_lo": s5["or_lo"], "or_fe_hi": s5["or_hi"],
        "p_fe": s5["logit_p"], "k": 1, "n": s5["n"],
    })
    pd.DataFrame([s5]).to_csv(OUT / "tables" / "sensitivity_S5_ipd.tsv", sep="\t", index=False)
    pd.DataFrame(sens).to_csv(OUT / "tables" / "sensitivity_meta_summary.tsv", sep="\t", index=False)

    # ---- TCGA context ----
    print("Loading TCGA-SKCM...", flush=True)
    tcga_e, tcga_c, tcga_notes = load_tcga_skcm()
    tumor = tcga_c["is_tumor"].fillna(False)
    tcga_x = tcga_e.loc[GENE_OF_INTEREST, tumor[tumor].index] if GENE_OF_INTEREST in tcga_e.index else pd.Series(dtype=float)
    lin_t = np.power(2.0, tcga_x) - 0.1
    tcga_det = {
        "n_tumor": int(tumor.sum()),
        "n_with_CLDN4": int(len(tcga_x)),
        "median_log2": float(tcga_x.median()) if len(tcga_x) else np.nan,
        "q25_log2": float(tcga_x.quantile(0.25)) if len(tcga_x) else np.nan,
        "q75_log2": float(tcga_x.quantile(0.75)) if len(tcga_x) else np.nan,
        "frac_detectable": float((lin_t > 0).mean()) if len(tcga_x) else np.nan,
        "frac_below_log2_0": float((tcga_x < 0).mean()) if len(tcga_x) else np.nan,
        "unit": tcga_notes["unit"],
    }
    # correlations in TCGA
    tcga_corr = []
    if len(tcga_x):
        for g in EPITHELIAL_CONTROLS + MELANOCYTE_GENES + KERATINOCYTE_SCORE_GENES:
            if g not in tcga_e.index:
                continue
            rho, p = stats.spearmanr(tcga_x, tcga_e.loc[g, tcga_x.index], nan_policy="omit")
            tcga_corr.append({"gene": g, "rho": float(rho), "p": float(p)})
        krt = gene_score(tcga_e.loc[:, tcga_x.index], KERATINOCYTE_SCORE_GENES)
        rho, p = stats.spearmanr(tcga_x, krt, nan_policy="omit")
        tcga_corr.append({"gene": "keratinocyte_score", "rho": float(rho), "p": float(p)})
        mel = gene_score(tcga_e.loc[:, tcga_x.index], MELANOCYTE_GENES)
        rho, p = stats.spearmanr(tcga_x, mel, nan_policy="omit")
        tcga_corr.append({"gene": "melanocyte_score", "rho": float(rho), "p": float(p)})
    pd.DataFrame([tcga_det]).to_csv(OUT / "tables" / "tcga_skcm_cldn4.tsv", sep="\t", index=False)
    pd.DataFrame(tcga_corr).to_csv(OUT / "tables" / "tcga_skcm_cldn4_correlations.tsv", sep="\t", index=False)

    # ---- figures ----
    fig_boxplots(primary)
    fig_forest(primary_df, meta)
    fig_genomewide(gw)
    fig_epithelial_corr(corr_wide)
    if len(tcga_x):
        fig_tcga(tcga_x)

    # ---- decision against the pre-specified positive-result rule ----
    cldn4_rank = int(cldn4_gw["p_fe_rank"].iloc[0]) if len(cldn4_gw) else None
    cldn4_pct = float(cldn4_gw["p_fe_percentile"].iloc[0]) if len(cldn4_gw) else None
    cldn4_fdr = float(cldn4_gw["bh_fdr"].iloc[0]) if len(cldn4_gw) else None
    crit = {
        "meta_p_lt_0.05": bool(meta["p_fe"] < 0.05) if np.isfinite(meta["p_fe"]) else False,
        "consistent_direction_ge_2_of_3": bool(max(meta["direction_pos"], meta["direction_neg"]) >= 2),
        "top_5pct_genomewide": bool(cldn4_pct is not None and cldn4_pct <= 5.0),
        "survives_krt_adjustment": bool(np.isfinite(meta_adj["p_fe"]) and meta_adj["p_fe"] < 0.05
                                        and np.sign(meta_adj["logor_fe"]) == np.sign(meta["logor_fe"])),
    }
    crit["all_four"] = all(crit.values())
    # positive-control sanity: is the assay capable of seeing a known ICI signal?
    ayers = ctrl_meta_df[ctrl_meta_df["gene"] == "Ayers_IFNG_6"]
    cd8a = ctrl_meta_df[ctrl_meta_df["gene"] == "CD8A"]
    pos_ok = False
    if len(ayers):
        pos_ok = bool(ayers["p_fe"].iloc[0] < 0.05 and ayers["or_fe"].iloc[0] > 1)
    if len(cd8a) and cd8a["p_fe"].iloc[0] < 0.05 and cd8a["or_fe"].iloc[0] > 1:
        pos_ok = True
    # if neither CD8A nor Ayers is nominally + and in the expected direction, underpowered
    underpowered = not pos_ok

    summary = {
        "gene": GENE_OF_INTEREST,
        "n_patients_primary": int(sum(c.shape[0] for _, c in primary.values())),
        "n_R": int(sum((c["response"] == "R").sum() for _, c in primary.values())),
        "n_NR": int(sum((c["response"] == "NR").sum() for _, c in primary.values())),
        "primary_meta": meta,
        "adjusted_krt_meta": meta_adj,
        "criteria": crit,
        "cldn4_genomewide_rank": cldn4_rank,
        "cldn4_genomewide_percentile": cldn4_pct,
        "cldn4_genomewide_fdr": cldn4_fdr,
        "n_genes_genomewide": int(len(gw)),
        "n_genes_nominal_p_lt_0.05": n_sig_nom,
        "n_genes_bh_fdr_lt_0.05": n_sig_fdr,
        "positive_control_detectable": pos_ok,
        "underpowered_flag": underpowered,
        "verdict": (
            "POSITIVE — meets all four pre-specified criteria"
            if crit["all_four"]
            else "NULL / INCONCLUSIVE — does not meet the pre-specified positive-result rule"
        ),
        "detectability": detect.to_dict("records"),
        "per_cohort": primary_df.to_dict("records"),
        "hugo_os": surv,
        "tcga": tcga_det,
        "load_notes": load_notes,
        "inventory": inventory.to_dict("records"),
        "s5": s5,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print("\nVERDICT:", summary["verdict"], flush=True)
    print("Criteria:", crit, flush=True)
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
