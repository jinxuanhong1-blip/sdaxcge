#!/usr/bin/env python3
"""Expression <-> ICI response axis (GEO GSE135222, Jung 2019).

27 advanced NSCLC patients treated with anti-PD-1/PD-L1, RNA-seq (TPM) + PFS.
Tests whether TACSTD2 (TROP2) / CLDN4 (Claudin-4) expression relates to:
  * durable clinical benefit proxy (PFS >= 6 months)   -> Mann-Whitney
  * progression-free survival (continuous)             -> univariable Cox HR
  * median-split PFS                                    -> Kaplan-Meier + log-rank
"""
from __future__ import annotations

import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import fable_common as fc
import fable_stats as fs

warnings.filterwarnings("ignore")

GENES = fc.TARGET_GENES + ["CD274", "CD8A"]


def load():
    expr = pd.read_csv(os.path.join(fc.PROC, "GSE135222_expression.csv"))
    clin = pd.read_csv(os.path.join(fc.PROC, "GSE135222_clinical.csv"))
    df = expr.merge(clin, on="sample", how="inner")
    for g in GENES:
        df[f"log2_{g}"] = np.log2(df[g].clip(lower=0) + 1)
    df["pfs_months"] = df["pfs_time_days"] / 30.44
    return df


def response_tests(df):
    rows = []
    for tgt in fc.TARGET_GENES:
        col = f"log2_{tgt}"
        yes = df.loc[df.DCB == "YES", col]
        no = df.loc[df.DCB == "NO", col]
        r = fs.group_compare(yes, no)  # "mut"=DCB-yes, "wt"=DCB-no
        r = {"target": tgt, "n_DCB": r["n_mut"], "n_noDCB": r["n_wt"],
             "median_DCB": r["median_mut"], "median_noDCB": r["median_wt"],
             "log2_diff": r["log2fc_median"], "cliffs_delta": r["cliffs_delta"],
             "mwu_p": r["mwu_p"]}
        rows.append(r)
    return pd.DataFrame(rows)


def cox_tests(df):
    from lifelines import CoxPHFitter
    rows = []
    for tgt in GENES:
        col = f"log2_{tgt}"
        d = df[[col, "pfs_months", "pfs_event"]].dropna()
        cph = CoxPHFitter()
        try:
            cph.fit(d, duration_col="pfs_months", event_col="pfs_event")
            hr = float(np.exp(cph.params_[col]))
            ci = np.exp(cph.confidence_intervals_.loc[col].values)
            p = float(cph.summary.loc[col, "p"])
        except Exception:
            hr, ci, p = np.nan, (np.nan, np.nan), np.nan
        rows.append({"target": tgt, "cox_hr_per_log2": hr,
                     "hr_lo": float(ci[0]), "hr_hi": float(ci[1]), "cox_p": p,
                     "n": int(d.shape[0])})
    return pd.DataFrame(rows)


def median_split_km(df):
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    rows = []
    fig, axes = plt.subplots(1, len(fc.TARGET_GENES), figsize=(12, 5))
    for ax, tgt in zip(axes, fc.TARGET_GENES):
        col = f"log2_{tgt}"
        med = df[col].median()
        df["_grp"] = np.where(df[col] > med, "high", "low")
        hi = df[df._grp == "high"]
        lo = df[df._grp == "low"]
        lr = logrank_test(hi.pfs_months, lo.pfs_months,
                          event_observed_A=hi.pfs_event, event_observed_B=lo.pfs_event)
        kmf = KaplanMeierFitter()
        for g, name, c in [(hi, "high", "#e41a1c"), (lo, "low", "#377eb8")]:
            kmf.fit(g.pfs_months, g.pfs_event, label=f"{tgt} {name} (n={g.shape[0]})")
            kmf.plot_survival_function(ax=ax, color=c, ci_show=False)
        ax.set_title(f"{tgt} median split (logrank p={lr.p_value:.3f})")
        ax.set_xlabel("Months")
        ax.set_ylabel("PFS probability")
        rows.append({"target": tgt, "logrank_p": float(lr.p_value),
                     "median_pfs_high": float(hi.pfs_months.median()),
                     "median_pfs_low": float(lo.pfs_months.median())})
    fig.tight_layout()
    out = os.path.join(fc.FIG, "geo_pfs_by_target_median_split.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return pd.DataFrame(rows), out


def corr_tests(df):
    rows = []
    for tgt in fc.TARGET_GENES:
        for ctx in ["CD8A", "CD274"]:
            rho, p = stats.spearmanr(df[f"log2_{tgt}"], df[f"log2_{ctx}"])
            rows.append({"target": tgt, "context": ctx,
                         "spearman_rho": float(rho), "p": float(p)})
    return pd.DataFrame(rows)


def main():
    df = load()
    print(f"GSE135222: n={df.shape[0]} patients, events={int(df.pfs_event.sum())}, "
          f"DCB(YES/NO)={dict(df.DCB.value_counts())}")
    resp = response_tests(df)
    cox = cox_tests(df)
    km, fig = median_split_km(df)
    corr = corr_tests(df)

    resp.to_csv(os.path.join(fc.TAB, "geo_expression_vs_dcb.csv"), index=False)
    cox.to_csv(os.path.join(fc.TAB, "geo_expression_vs_pfs_cox.csv"), index=False)
    km.to_csv(os.path.join(fc.TAB, "geo_expression_median_split_km.csv"), index=False)
    corr.to_csv(os.path.join(fc.TAB, "geo_target_immune_correlation.csv"), index=False)

    print("\n=== Expression by DCB proxy (PFS>=6mo), Mann-Whitney ===")
    print(resp.round(4).to_string(index=False))
    print("\n=== Univariable Cox: expression -> PFS (HR per log2 unit) ===")
    print(cox.round(4).to_string(index=False))
    print("\n=== Median-split KM log-rank ===")
    print(km.round(4).to_string(index=False))
    print("\n=== Target vs immune context (Spearman) ===")
    print(corr.round(4).to_string(index=False))
    print(f"\nfigure: {os.path.relpath(fig, fc.REPO)}")


if __name__ == "__main__":
    np.random.seed(0)
    main()
