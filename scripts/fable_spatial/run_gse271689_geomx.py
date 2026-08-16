"""GSE271689 (DSP-GeoMx WTA, NSCLC first-line ICI) — TACSTD2/CLDN4 analysis.

Uses the per-patient compartment matrices + clinical outcomes published as
source data of Nature Genetics s41588-025-02351-7 (Fig. 6b/6d/6f/6h), which
correspond to the GEO series GSE271689 (Yale discovery + Greek validation).

Analyses:
 1) tumor-compartment (PanCK AOI) TACSTD2/CLDN4 vs stromal-compartment
    immune activation (T-effector score) — Spearman per cohort;
 2) overall survival: Cox PH per gene (per-SD of transformed expression),
    univariable and adjusted for the stromal immune score; KM median split
    log-rank.
"""

import os
import sys

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from common import TEFF_MARKERS, TARGETS, ensure_results, spearman

DATA = "/workspace/data/fable_spatial/gse271689"
OUT = ensure_results()

CLIN_6B = [
    "No", "Core.Type", "Biopsy_ITx", "Smoking_Status", "PD_L1_TPS_IC_other_biopsy",
    "Mutation_Status", "Stage_at_biopsy", "Stage_at_ITx", "Site_of_biopsy", "Histotype",
    "Metastatic_Site_at_ITx", "Line_ITx", "Agent_ITx", "Concurrent_treatment",
    "OR_1st_scan", "BOR", "response", "Start_date_ITx", "End_date_ITx", "Cycles_ITx",
    "Progression_date", "PFS_Days", "PFS_Index", "PFS_2Years_months", "PFS_2Years_Index",
    "PFS_5Years_months", "PFS_5Years_Index", "DOD_or_Last_FU", "OS_Days", "OS_Index",
    "OS_5Years_months", "OS_5Years_Index", "OS_2Years_months", "OS_2Years_Index",
    "irAEs", "Previous_surgery", "Previous_RT", "Previous_adjuvant_chemotherapy",
    "Previous_line", "DCB6", "LTB12", "LTB24", "group", "binary_score",
    "Unnamed: 0", "Spot_ID", "X",
]


def rint(v):
    """Rank-based inverse normal transform (for the batch-corrected Greek data)."""
    r = stats.rankdata(v)
    return stats.norm.ppf((r - 0.5) / len(r))


def teff_score(df, genes, transform):
    avail = [g for g in genes if g in df.columns]
    Z = np.column_stack([transform(df[g].astype(float).to_numpy()) for g in avail])
    Z = (Z - Z.mean(0)) / np.where(Z.std(0) == 0, 1, Z.std(0))
    return Z.mean(1), avail


def cox_row(time, event, x, covar=None, label=""):
    d = pd.DataFrame({"time": time, "event": event, "x": (x - np.mean(x)) / np.std(x)})
    if covar is not None:
        d["covar"] = (covar - np.mean(covar)) / np.std(covar)
    d = d.dropna()
    d = d[d.time > 0]
    cph = CoxPHFitter()
    cph.fit(d, duration_col="time", event_col="event")
    s = cph.summary.loc["x"]
    return dict(model=label, n=len(d), events=int(d.event.sum()),
                hr_per_sd=round(float(s["exp(coef)"]), 3),
                ci_low=round(float(s["exp(coef) lower 95%"]), 3),
                ci_high=round(float(s["exp(coef) upper 95%"]), 3),
                p=float(s["p"]))


def km_median_split(time, event, x):
    med = np.median(x)
    hi = x > med
    res = logrank_test(time[hi], time[~hi], event_observed_A=event[hi], event_observed_B=event[~hi])
    med_hi = kmq(time[hi], event[hi])
    med_lo = kmq(time[~hi], event[~hi])
    return dict(n_high=int(hi.sum()), n_low=int((~hi).sum()),
                median_surv_high=med_hi, median_surv_low=med_lo,
                logrank_p=float(res.p_value))


def kmq(t, e):
    from lifelines import KaplanMeierFitter
    km = KaplanMeierFitter().fit(t, e)
    m = km.median_survival_time_
    return float(m) if np.isfinite(m) else np.nan


def main():
    xl = pd.ExcelFile(os.path.join(DATA, "moesm10.xlsx"))
    yale_t = pd.read_excel(xl, "source_data_Figure_6b")
    yale_s = pd.read_excel(xl, "source_data_Figure_6f")
    greek_t = pd.read_excel(xl, "source_data_Figure_6d")
    greek_s = pd.read_excel(xl, "source_data_Figure_6h")

    corr_rows, cox_rows, km_rows = [], [], []

    # ------------------------------------------------------------- Yale ----
    yale_t = yale_t.set_index("Spot_ID")
    yale_s = yale_s.set_index("Spot_ID")
    common_pt = yale_t.index.intersection(yale_s.index)
    yt, ys = yale_t.loc[common_pt], yale_s.loc[common_pt]
    log2p = lambda v: np.log2(np.asarray(v, float) + 1)

    imm_y, genes_y = teff_score(ys, TEFF_MARKERS, log2p)
    os_days = yt["OS_Days"].astype(float).to_numpy()
    os_evt = yt["OS_Index"].astype(int).to_numpy()
    # univariable models use all tumor-compartment patients (n=37)
    os_days_all = yale_t["OS_Days"].astype(float).to_numpy()
    os_evt_all = yale_t["OS_Index"].astype(int).to_numpy()

    for gene in TARGETS:
        expr = log2p(yt[gene])
        expr_all = log2p(yale_t[gene])
        rho, p = spearman(expr, imm_y)
        corr_rows.append(dict(cohort="Yale (WTA, n=%d)" % len(yt), gene=gene,
                              vs="stromal T-effector score", spearman_rho=round(rho, 3),
                              p=p, n=len(yt)))
        cox_rows.append(dict(cohort="Yale", gene=gene, endpoint="OS",
                             **cox_row(os_days_all, os_evt_all, expr_all,
                                       label="univariable")))
        cox_rows.append(dict(cohort="Yale", gene=gene, endpoint="OS",
                             **cox_row(os_days, os_evt, expr, covar=imm_y,
                                       label="adjusted for stromal T-eff score")))
        km_rows.append(dict(cohort="Yale", gene=gene, endpoint="OS",
                            **km_median_split(os_days_all, os_evt_all, expr_all)))
    # immune score itself
    cox_rows.append(dict(cohort="Yale", gene="stromal_Teff_score", endpoint="OS",
                         **cox_row(os_days, os_evt, imm_y, label="univariable")))

    # ------------------------------------------------------------ Greek ----
    greek_t = greek_t.set_index("ROILabel")
    greek_s = greek_s.set_index("ROILabel")
    common_roi = greek_t.index.intersection(greek_s.index)
    gt, gs = greek_t.loc[common_roi], greek_s.loc[common_roi]

    imm_g, genes_g = teff_score(gs, TEFF_MARKERS, rint)
    os1 = gt["OS_1"].astype(float).to_numpy()
    death = (gt["Death"].astype(str).str.strip().str.lower() == "yes").astype(int).to_numpy()

    for gene in TARGETS:
        expr = rint(gt[gene].astype(float).to_numpy())
        rho, p = spearman(expr, imm_g)
        corr_rows.append(dict(cohort="Greek (WTA, n=%d)" % len(gt), gene=gene,
                              vs="stromal T-effector score", spearman_rho=round(rho, 3),
                              p=p, n=len(gt)))
        cox_rows.append(dict(cohort="Greek", gene=gene, endpoint="OS",
                             **cox_row(os1, death, expr, label="univariable")))
        cox_rows.append(dict(cohort="Greek", gene=gene, endpoint="OS",
                             **cox_row(os1, death, expr, covar=imm_g,
                                       label="adjusted for stromal T-eff score")))
        km_rows.append(dict(cohort="Greek", gene=gene, endpoint="OS",
                            **km_median_split(os1, death, expr)))
    cox_rows.append(dict(cohort="Greek", gene="stromal_Teff_score", endpoint="OS",
                         **cox_row(os1, death, imm_g, label="univariable")))

    pd.DataFrame(corr_rows).to_csv(os.path.join(OUT, "gse271689_tumor_vs_stromal_immune_corr.csv"), index=False)
    pd.DataFrame(cox_rows).to_csv(os.path.join(OUT, "gse271689_os_cox.csv"), index=False)
    pd.DataFrame(km_rows).to_csv(os.path.join(OUT, "gse271689_os_km_median_split.csv"), index=False)

    print("T-effector genes used  Yale:", genes_y)
    print("T-effector genes used Greek:", genes_g)
    print("\n== correlations ==")
    print(pd.DataFrame(corr_rows).to_string(index=False))
    print("\n== Cox OS ==")
    print(pd.DataFrame(cox_rows).to_string(index=False))
    print("\n== KM median split ==")
    print(pd.DataFrame(km_rows).to_string(index=False))


if __name__ == "__main__":
    main()
