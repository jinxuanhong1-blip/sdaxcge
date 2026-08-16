#!/usr/bin/env python3
"""Genotype <-> ICI response axis (public NSCLC immunotherapy cohorts).

Cohorts (cBioPortal): Rizvi2015, Hellmann2018, Rizvi2018.
Endpoints:
  * Durable clinical benefit (DCB, harmonized YES/NO) -> Fisher exact per genotype
  * Progression-free survival -> Kaplan-Meier + log-rank + univariable Cox (HR)
Pooled DCB association via Mantel-Haenszel odds ratio.
"""
from __future__ import annotations

import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fable_common as fc
import fable_stats as fs

warnings.filterwarnings("ignore")

DCB_MAP = {
    "Durable clinical benefit beyond 6 months": "YES",
    "No durable benefit": "NO",
    "Not reached 6 months follow-up": np.nan,
    "Durable Clinical Benefit": "YES",
    "No Durable Benefit": "NO",
    "YES": "YES", "NO": "NO", "NE": np.nan,
}

ICI = ["Rizvi2015", "Hellmann2018", "Rizvi2018"]
GENO_TESTS = ["STK11", "KEAP1", "KRAS", "STK11_or_KEAP1"]


def load(label):
    clin = pd.read_csv(os.path.join(fc.PROC, f"{label}_clinical.csv"))
    geno = pd.read_csv(os.path.join(fc.PROC, f"{label}_genotype.csv"))
    df = clin.merge(geno, on="sampleId", how="inner")
    df["DCB"] = df["DURABLE_CLINICAL_BENEFIT"].map(DCB_MAP)
    df["STK11_or_KEAP1_mut"] = ((df.STK11_mut == 1) | (df.KEAP1_mut == 1)).astype(int)
    if "PFS_STATUS" in df.columns:
        df["pfs_event"] = df["PFS_STATUS"].astype(str).str.split(":").str[0]
        df["pfs_event"] = pd.to_numeric(df["pfs_event"], errors="coerce")
    if "PFS_MONTHS" in df.columns:
        df["pfs_months"] = pd.to_numeric(df["PFS_MONTHS"], errors="coerce")
    return df


def dcb_tests(df, label):
    rows = []
    sub = df[df.DCB.isin(["YES", "NO"])]
    for gt in GENO_TESTS:
        flag = f"{gt}_mut"
        r = fs.fisher_dcb(sub[flag] == 1, sub.DCB == "YES")
        r.update({"cohort": label, "genotype": gt})
        rows.append(r)
    return pd.DataFrame(rows)


def kras_restricted_dcb(df, label):
    """Within KRAS-mutant tumors, effect of STK11 / KEAP1 co-mutation on DCB
    (the Skoulidis 2018 hypothesis: STK11/LKB1 drives ICI resistance in KRAS-mut)."""
    rows = []
    sub = df[(df.KRAS_mut == 1) & df.DCB.isin(["YES", "NO"])]
    for co in ["STK11", "KEAP1"]:
        r = fs.fisher_dcb(sub[f"{co}_mut"] == 1, sub.DCB == "YES")
        r.update({"cohort": label, "background": "KRAS-mut", "co_mutation": co})
        rows.append(r)
    return pd.DataFrame(rows)


def cox_km(df, label):
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test
    out = []
    if "pfs_event" not in df.columns or "pfs_months" not in df.columns:
        return pd.DataFrame(out)
    d = df.dropna(subset=["pfs_event", "pfs_months"])
    d = d[d.pfs_months >= 0]
    for gt in GENO_TESTS:
        flag = f"{gt}_mut"
        if d[flag].nunique() < 2:
            continue
        g1 = d[d[flag] == 1]
        g0 = d[d[flag] == 0]
        lr = logrank_test(g1.pfs_months, g0.pfs_months,
                          event_observed_A=g1.pfs_event, event_observed_B=g0.pfs_event)
        cph = CoxPHFitter()
        try:
            cph.fit(d[[flag, "pfs_months", "pfs_event"]],
                    duration_col="pfs_months", event_col="pfs_event")
            hr = float(np.exp(cph.params_[flag]))
            ci = np.exp(cph.confidence_intervals_.loc[flag].values)
            cph_p = float(cph.summary.loc[flag, "p"])
        except Exception:
            hr, ci, cph_p = np.nan, (np.nan, np.nan), np.nan
        out.append({
            "cohort": label, "genotype": gt,
            "n_mut": int((d[flag] == 1).sum()), "n_wt": int((d[flag] == 0).sum()),
            "median_pfs_mut": float(g1.pfs_months.median()),
            "median_pfs_wt": float(g0.pfs_months.median()),
            "logrank_p": float(lr.p_value), "cox_hr": hr,
            "cox_hr_lo": float(ci[0]), "cox_hr_hi": float(ci[1]), "cox_p": cph_p,
        })
    return pd.DataFrame(out)


def km_figure(dfs):
    from lifelines import KaplanMeierFitter
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, label in zip(axes, ["Hellmann2018", "Rizvi2018"]):
        d = dfs[label].dropna(subset=["pfs_event", "pfs_months"])
        kmf = KaplanMeierFitter()
        for val, name, color in [(0, "STK11/KEAP1 WT", "#377eb8"),
                                 (1, "STK11/KEAP1 MUT", "#e41a1c")]:
            m = d.STK11_or_KEAP1_mut == val
            if m.sum() == 0:
                continue
            kmf.fit(d.loc[m, "pfs_months"], d.loc[m, "pfs_event"],
                    label=f"{name} (n={int(m.sum())})")
            kmf.plot_survival_function(ax=ax, color=color, ci_show=True)
        ax.set_title(f"{label}: PFS by STK11/KEAP1 status")
        ax.set_xlabel("Months")
        ax.set_ylabel("PFS probability")
    fig.tight_layout()
    out = os.path.join(fc.FIG, "ici_pfs_by_stk11_keap1.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def dcb_barplot(dcb_all):
    piv = dcb_all[dcb_all.genotype == "STK11_or_KEAP1"]
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(piv))
    w = 0.38
    ax.bar(x - w / 2, piv.dcb_rate_wt * 100, w, label="STK11/KEAP1 WT", color="#377eb8")
    ax.bar(x + w / 2, piv.dcb_rate_mut * 100, w, label="STK11/KEAP1 MUT", color="#e41a1c")
    ax.set_xticks(x)
    ax.set_xticklabels(piv.cohort)
    ax.set_ylabel("Durable clinical benefit rate (%)")
    ax.set_title("ICI durable clinical benefit by STK11/KEAP1 status")
    for i, row in enumerate(piv.itertuples()):
        ax.text(i + w / 2, row.dcb_rate_mut * 100 + 1, f"p={row.fisher_p:.3f}",
                ha="center", fontsize=8)
    ax.legend()
    fig.tight_layout()
    out = os.path.join(fc.FIG, "ici_dcb_by_stk11_keap1.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main():
    dfs = {lab: load(lab) for lab in ICI}
    dcb_all = pd.concat([dcb_tests(dfs[l], l) for l in ICI], ignore_index=True)
    cox_all = pd.concat([cox_km(dfs[l], l) for l in ICI], ignore_index=True)
    kras_all = pd.concat([kras_restricted_dcb(dfs[l], l) for l in ICI], ignore_index=True)

    # pooled Mantel-Haenszel OR for each genotype across cohorts
    pooled = []
    for gt in GENO_TESTS:
        tabs = dcb_all[dcb_all.genotype == gt]["table"].tolist()
        pooled.append({"genotype": gt, "mh_odds_ratio": fs.mantel_haenszel_or(tabs),
                       "n_cohorts": len(tabs)})
    pooled = pd.DataFrame(pooled)

    # pooled KRAS-restricted co-mutation OR
    kras_pooled = []
    for co in ["STK11", "KEAP1"]:
        tabs = kras_all[kras_all.co_mutation == co]["table"].tolist()
        kras_pooled.append({"co_mutation": co,
                            "mh_odds_ratio": fs.mantel_haenszel_or(tabs)})
    kras_pooled = pd.DataFrame(kras_pooled)

    dcb_all.to_csv(os.path.join(fc.TAB, "ici_genotype_vs_dcb.csv"), index=False)
    cox_all.to_csv(os.path.join(fc.TAB, "ici_genotype_vs_pfs.csv"), index=False)
    pooled.to_csv(os.path.join(fc.TAB, "ici_genotype_dcb_pooled_MH.csv"), index=False)
    kras_all.to_csv(os.path.join(fc.TAB, "ici_kras_restricted_dcb.csv"), index=False)

    fig1 = km_figure(dfs)
    fig2 = dcb_barplot(dcb_all)

    print("=== DCB by genotype (Fisher exact) ===")
    print(dcb_all[["cohort", "genotype", "n_mut", "n_wt", "dcb_rate_mut",
                   "dcb_rate_wt", "odds_ratio", "fisher_p"]].round(4).to_string(index=False))
    print("\n=== PFS by genotype (log-rank + Cox HR) ===")
    print(cox_all[["cohort", "genotype", "n_mut", "n_wt", "median_pfs_mut",
                   "median_pfs_wt", "logrank_p", "cox_hr", "cox_hr_lo",
                   "cox_hr_hi", "cox_p"]].round(4).to_string(index=False))
    print("\n=== Within KRAS-mutant: co-mutation vs DCB (Fisher) ===")
    print(kras_all[["cohort", "co_mutation", "n_mut", "n_wt", "dcb_rate_mut",
                    "dcb_rate_wt", "odds_ratio", "fisher_p"]].round(4).to_string(index=False))
    print("\n=== Pooled Mantel-Haenszel OR (DCB, mut vs wt) ===")
    print(pooled.round(4).to_string(index=False))
    print("--- KRAS-restricted pooled MH OR ---")
    print(kras_pooled.round(4).to_string(index=False))
    print(f"\nfigures: {os.path.relpath(fig1, fc.REPO)}, {os.path.relpath(fig2, fc.REPO)}")


if __name__ == "__main__":
    main()
