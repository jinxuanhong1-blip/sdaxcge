#!/usr/bin/env python3
"""Association of TACSTD2/CLDN4 with ICI RESPONSE after purity adjustment.

Core question (this continuation): in open lung ICI bulk cohorts that actually
measure both genes, is expression associated with clinical response once tumor
purity / epithelial content is accounted for?

Methods (all real, no fabricated numbers):
  * OLS residualize log-expression on an epithelial (tumor-content) score and,
    separately, on a leukocyte (infiltrate) score.
  * Mann-Whitney U on residuals (NR vs R), two-sided + one-sided NR>R.
  * Spearman of expression vs binary NR (1=NR, 0=R); partial Spearman controlling
    epithelial or leukocyte score.
  * Logistic regression: NR ~ gene  and  NR ~ gene + epithelial (Firth not used;
    if MLE fails / separates, that is recorded, not invented).
  * GSE135222: Cox PFS ~ gene  and  PFS ~ gene + epithelial; plus a secondary
    binary DCB endpoint defined as PFS time >= 180 days (published 6-month DCB
    convention used on this exact GEO series; GEO itself deposits PFS, not RECIST).

Durvalumab GSE253564 / GSE248378: TACSTD2+CLDN4 present, but GEO has no
responder/MPR label and no public per-sample table was retrievable — excluded
from the RESPONSE tests (kept only for CD8/NK correlation elsewhere).
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from lifelines import CoxPHFitter

sys.path.insert(0, os.path.dirname(__file__))
from purity_corrected import (
    load_gse126044, load_gse166449, load_gse207422, load_gse135222,
)
import signatures as S

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = "/tmp/ici_bulk_data"
RES = "/workspace/results/fable_ici_bulk"
os.makedirs(f"{RES}/tables", exist_ok=True)
os.makedirs(f"{RES}/figures", exist_ok=True)

mwu_rows = []
corr_rows = []
logit_rows = []
cox_rows = []


def ols_residual(y, x):
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    m = ~(np.isnan(y) | np.isnan(x))
    resid = np.full_like(y, np.nan)
    if m.sum() < 3:
        return resid
    X = np.column_stack([np.ones(m.sum()), x[m]])
    b, *_ = np.linalg.lstsq(X, y[m], rcond=None)
    resid[m] = y[m] - X @ b
    return resid


def fit_logit(y, X_df):
    """Return dict; never invent coefficients if MLE fails."""
    X = sm.add_constant(X_df, has_constant="add")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = sm.Logit(y, X).fit(disp=False, maxiter=200)
        if not m.mle_retvals.get("converged", False):
            return {"status": "not_converged"}
        out = {"status": "ok", "n": int(len(y)), "n_NR": int(y.sum()),
               "n_R": int((1 - y).sum())}
        for name in X_df.columns:
            out[f"OR_{name}"] = float(np.exp(m.params[name]))
            ci = m.conf_int().loc[name]
            out[f"OR_{name}_CI_low"] = float(np.exp(ci.iloc[0]))
            out[f"OR_{name}_CI_high"] = float(np.exp(ci.iloc[1]))
            out[f"p_{name}"] = float(m.pvalues[name])
        return out
    except Exception as e:  # PerfectSeparation, LinAlgError, etc.
        return {"status": f"failed:{type(e).__name__}"}


def load_gse135222_dcb():
    """Secondary binary endpoint: DCB = PFS time >= 180 days.

    This 6-month DCB/NDB cut is the convention used by multiple published
    reanalyses of GSE135222 (Jung et al.). GEO deposits PFS event+time only;
    no patient was censored before 180 days in this series (verified below).
    """
    gse, expr, labels = load_gse135222()
    surv = labels["surv"]
    resp = {}
    n_early_censor = 0
    for c, (ev, t) in surv.items():
        if t >= 180:
            resp[c] = "responder"  # DCB
        else:
            if ev == 0:
                n_early_censor += 1
                resp[c] = None
            else:
                resp[c] = "non-responder"  # NDB
    if n_early_censor:
        print(f"WARNING GSE135222: {n_early_censor} censored before 180d excluded")
    return "GSE135222_DCB", expr, {
        "type": "binary", "resp": resp,
        "responder": "responder", "nonresponder": "non-responder",
        "endpoint": "DCB=PFS>=180d (published 6-month convention on this series)",
    }


def analyze_binary(gse, expr, labels):
    epi, _ = S.signature_score(expr, S.EPITHELIAL)
    leuk, _ = S.signature_score(expr, S.LEUKOCYTE)
    resp = labels["resp"]
    r_lab, n_lab = labels["responder"], labels["nonresponder"]
    cols = [c for c in expr.columns if resp.get(c) in (r_lab, n_lab)]
    y = np.array([1 if resp[c] == n_lab else 0 for c in cols], float)  # 1=NR
    endpoint = labels.get("endpoint", "binary responder/non-responder")

    for gene in S.TARGET_GENES:
        if gene not in expr.index:
            continue
        raw = expr.loc[gene, cols].astype(float).values
        e = epi.reindex(cols).values
        l = leuk.reindex(cols).values
        res_epi = ols_residual(raw, e)
        res_leuk = ols_residual(raw, l)

        for measure, vals in [("raw", raw),
                              ("resid_epithelial", res_epi),
                              ("resid_leukocyte", res_leuk)]:
            nr = vals[y == 1]
            r = vals[y == 0]
            mw = S.mannwhitney_directional(nr, r)
            mw.update({"dataset": gse, "gene": gene, "measure": measure,
                       "endpoint": endpoint, "hypothesis": "NR > R"})
            mwu_rows.append(mw)

        # Spearman vs binary NR
        for measure, vals, cov in [
            ("raw", raw, None),
            ("partial_epi", raw, [e]),
            ("partial_leuk", raw, [l]),
            ("resid_epithelial", res_epi, None),
            ("resid_leukocyte", res_leuk, None),
        ]:
            if cov is None:
                rho, p, n = S.spearman(vals, y)
            else:
                rho, p, n = S.partial_spearman(vals, y, cov)
            corr_rows.append({
                "dataset": gse, "gene": gene, "measure": measure,
                "n": n, "n_NR": int(y.sum()), "n_R": int((1 - y).sum()),
                "spearman_vs_NR": rho, "p": p, "endpoint": endpoint,
                "note": "positive rho => higher expr in NR",
            })

        # logistic
        dfX = pd.DataFrame({"gene": raw, "epi": e, "leuk": l}, index=cols)
        uni = fit_logit(y, dfX[["gene"]])
        uni.update({"dataset": gse, "gene": gene, "model": "NR ~ gene",
                    "endpoint": endpoint})
        logit_rows.append(uni)
        multi = fit_logit(y, dfX[["gene", "epi"]])
        multi.update({"dataset": gse, "gene": gene,
                      "model": "NR ~ gene + epithelial", "endpoint": endpoint})
        logit_rows.append(multi)

        _box(gse, gene, raw, res_epi, y)


def analyze_survival(gse, expr, labels):
    epi, _ = S.signature_score(expr, S.EPITHELIAL)
    leuk, _ = S.signature_score(expr, S.LEUKOCYTE)
    surv = labels["surv"]
    cols = [c for c in expr.columns if c in surv]
    time = np.array([surv[c][1] for c in cols], float)
    event = np.array([surv[c][0] for c in cols], float)
    e = epi.reindex(cols).values
    l = leuk.reindex(cols).values
    for gene in S.TARGET_GENES:
        raw = expr.loc[gene, cols].astype(float).values
        res_epi = ols_residual(raw, e)
        df = pd.DataFrame({"time": time, "event": event, "gene": raw,
                           "epi": e, "leuk": l, "resid_epi": res_epi})
        for cols_fit, model in [
            (["time", "event", "gene"], "PFS ~ gene"),
            (["time", "event", "gene", "epi"], "PFS ~ gene + epithelial"),
            (["time", "event", "gene", "leuk"], "PFS ~ gene + leukocyte"),
            (["time", "event", "resid_epi"], "PFS ~ resid(gene|epithelial)"),
        ]:
            cph = CoxPHFitter()
            try:
                cph.fit(df[cols_fit], "time", "event")
                key = "gene" if "gene" in cols_fit else "resid_epi"
                cox_rows.append({
                    "dataset": gse, "gene": gene, "model": model,
                    "n": len(df), "n_events": int(event.sum()),
                    "HR": float(np.exp(cph.params_[key])),
                    "HR_CI_low": float(np.exp(cph.confidence_intervals_.loc[key].iloc[0])),
                    "HR_CI_high": float(np.exp(cph.confidence_intervals_.loc[key].iloc[1])),
                    "p": float(cph.summary.loc[key, "p"]),
                    "status": "ok",
                    "note": "HR>1 => higher expr worse PFS",
                })
            except Exception as ex:
                cox_rows.append({"dataset": gse, "gene": gene, "model": model,
                                 "status": f"failed:{type(ex).__name__}"})


def _box(gse, gene, raw, resid, y):
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.6))
    for ax, vals, title in [
        (axes[0], raw, "raw log-expr"),
        (axes[1], resid, "OLS residual | epithelial"),
    ]:
        data = [vals[y == 0], vals[y == 1]]
        bp = ax.boxplot(data, tick_labels=[f"R\n(n={int((y==0).sum())})",
                                           f"NR\n(n={int((y==1).sum())})"],
                        widths=0.55, showfliers=False, patch_artist=True)
        for p, c in zip(bp["boxes"], ["#7dcea0", "#e6b0aa"]):
            p.set_facecolor(c)
        rng = np.random.default_rng(1)
        for i, d in enumerate(data, start=1):
            ax.scatter(rng.normal(i, 0.05, len(d)), d, s=16, color="#222",
                       alpha=0.75, zorder=3)
        mw = S.mannwhitney_directional(data[1], data[0])
        ax.set_title(f"{title}\nAUC(NR>R)={mw['auc_high_gt_low']:.2f} "
                     f"p1={mw['p_one_sided_greater']:.2f}")
        ax.set_ylabel(gene)
    fig.suptitle(f"{gse}: {gene} vs response", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{RES}/figures/{gse}_{gene}_response_purity.png", dpi=140)
    plt.close(fig)


def forest_auc(mwu):
    """Forest of AUC(NR>R) for raw vs residual-epithelial, core genes."""
    sub = mwu[mwu["measure"].isin(["raw", "resid_epithelial"])].copy()
    sub = sub.sort_values(["gene", "dataset", "measure"])
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    y = 0
    yticks, ylabels = [], []
    for gene in ["TACSTD2", "CLDN4"]:
        for ds in sub["dataset"].unique():
            rows = sub[(sub["gene"] == gene) & (sub["dataset"] == ds)]
            if rows.empty:
                continue
            for _, r in rows.iterrows():
                color = "#c0392b" if r["measure"] == "raw" else "#2471a3"
                ax.plot(r["auc_high_gt_low"], y, "o", color=color, ms=7)
                ax.plot([0.5, r["auc_high_gt_low"]], [y, y], color=color, lw=1.4)
                yticks.append(y)
                ylabels.append(f"{ds} {gene} [{r['measure']}]")
                y -= 1
            y -= 0.3
        y -= 0.4
    ax.axvline(0.5, color="#888", ls="--", lw=1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("AUC (NR > R)   >0.5 supports higher expr in non-responders")
    ax.set_xlim(0.15, 0.9)
    ax.set_title("Response association: raw vs epithelial-residualized")
    fig.tight_layout()
    fig.savefig(f"{RES}/figures/response_auc_forest.png", dpi=140)
    plt.close(fig)


def main():
    binaries = [load_gse126044(), load_gse166449(), load_gse207422(),
                load_gse135222_dcb()]
    for gse, expr, labels in binaries:
        print(f"binary {gse}: n labeled =",
              sum(1 for v in labels["resp"].values() if v))
        analyze_binary(gse, expr, labels)

    gse, expr, labels = load_gse135222()
    analyze_survival(gse, expr, labels)

    mwu = pd.DataFrame(mwu_rows)
    corr = pd.DataFrame(corr_rows)
    logit = pd.DataFrame(logit_rows)
    cox = pd.DataFrame(cox_rows)
    mwu.to_csv(f"{RES}/tables/response_after_purity_mwu.csv", index=False)
    corr.to_csv(f"{RES}/tables/response_after_purity_spearman.csv", index=False)
    logit.to_csv(f"{RES}/tables/response_after_purity_logistic.csv", index=False)
    cox.to_csv(f"{RES}/tables/response_after_purity_cox.csv", index=False)
    forest_auc(mwu)

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)
    print("\n=== MWU NR>R (raw / resid|epi / resid|leuk) ===")
    print(mwu[["dataset", "gene", "measure", "n_high", "n_low",
               "auc_high_gt_low", "p_one_sided_greater", "p_two_sided"]].to_string(index=False))
    print("\n=== Spearman vs NR (positive = higher in NR) ===")
    print(corr.to_string(index=False))
    print("\n=== Logistic (OR per unit gene; status recorded if MLE fails) ===")
    cols = [c for c in ["dataset", "gene", "model", "status", "n", "n_NR", "n_R",
                        "OR_gene", "OR_gene_CI_low", "OR_gene_CI_high", "p_gene"]
            if c in logit.columns]
    print(logit[cols].to_string(index=False))
    print("\n=== Cox PFS (GSE135222) ===")
    print(cox.to_string(index=False))


if __name__ == "__main__":
    main()
