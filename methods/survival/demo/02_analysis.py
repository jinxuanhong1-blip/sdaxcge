#!/usr/bin/env python3
"""Worked ICI-biomarker survival analysis: TACSTD2 and CLDN4 on two public GEO PFS cohorts.

Everything printed here is computed from the harmonised files in data/ (produced by
01_fetch_geo.py). No number is hard-coded. Run 01_fetch_geo.py first.

Primary genes: TACSTD2 (TROP2), CLDN4.
Comparator (not a claim): CD8A.
Exploratory: mean z of TACSTD2 and CLDN4.

Discovery = GSE135222 (NSCLC, anti-PD-1/PD-L1, PFS in days).
Validation = GSE190265 (NSCLC, anti-PD-1, PFS in months).
Neither series deposits OS.

Outputs: results/demo_results.md, results/demo_results.json, results/fig_*.png
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import AalenJohansenFitter, CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
from lifelines.utils import median_survival_times, restricted_mean_survival_time
from scipy.stats import chi2, norm
from sksurv.metrics import concordance_index_ipcw, cumulative_dynamic_auc
from sksurv.util import Surv

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RESULTS = HERE / "results"
SEED = 20240816
PRIMARY = ["TACSTD2", "CLDN4"]
COMPARATOR = ["CD8A"]
GENES = PRIMARY + COMPARATOR

R: dict[str, object] = {}
LINES: list[str] = []


def say(text: str = "") -> None:
    LINES.append(text)
    print(text)


def fmt(x, digits: int = 3) -> str:
    if x is None:
        return "NA"
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return "NA"
    return "NA" if not np.isfinite(xf) else f"{xf:.{digits}f}"


def load(gse: str, ref_stats: pd.DataFrame | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    clin = pd.read_csv(DATA / f"{gse}_clinical.csv")
    expr = pd.read_csv(DATA / f"{gse}_genes_tpm.csv", index_col=0)
    expr.index = expr.index.astype(str)
    clin["sample_id"] = clin["sample_id"].astype(str)
    logexpr = np.log2(expr[GENES] + 1.0)
    stats = pd.DataFrame({"mean": logexpr.mean(), "sd": logexpr.std(ddof=1)})
    z_internal = (logexpr - stats["mean"]) / stats["sd"]
    for g in GENES:
        clin = clin.merge(
            logexpr[g].rename(f"{g}_log2"),
            left_on="sample_id",
            right_index=True,
        )
        clin = clin.merge(
            z_internal[g].rename(f"{g}_z"),
            left_on="sample_id",
            right_index=True,
        )
    clin["TJ_mean_z"] = clin[[f"{g}_z" for g in PRIMARY]].mean(axis=1)
    if ref_stats is not None:
        z_ref = (logexpr - ref_stats["mean"]) / ref_stats["sd"]
        for g in GENES:
            clin[f"{g}_z_ref"] = clin["sample_id"].map(z_ref[g])
        clin["TJ_mean_z_ref"] = clin[[f"{g}_z_ref" for g in PRIMARY]].mean(axis=1)
    else:
        for g in GENES:
            clin[f"{g}_z_ref"] = clin[f"{g}_z"]
        clin["TJ_mean_z_ref"] = clin["TJ_mean_z"]
    return clin, stats


def surv_array(df: pd.DataFrame) -> np.ndarray:
    return Surv.from_arrays(event=df["event"].astype(bool), time=df["time_months"])


def events_needed(hr: float, alpha: float = 0.05, power: float = 0.80, p: float = 0.5) -> float:
    """Schoenfeld sample-size for a two-group log-rank / Cox, equal allocation."""
    za = norm.ppf(1 - alpha / 2)
    zb = norm.ppf(power)
    return (za + zb) ** 2 / (p * (1 - p) * np.log(hr) ** 2)


def describe(df: pd.DataFrame, label: str) -> dict:
    kmf = KaplanMeierFitter().fit(df["time_months"], df["event"])
    med = kmf.median_survival_time_
    med_ci = median_survival_times(kmf.confidence_interval_)
    lo = float(med_ci.iloc[0, 0])
    hi = float(med_ci.iloc[0, 1])
    rev = KaplanMeierFitter().fit(df["time_months"], 1 - df["event"])
    milestones = {}
    for t in (3.0, 6.0, 12.0):
        at_risk = int((df["time_months"] >= t).sum())
        if t <= df["time_months"].max():
            s = float(kmf.predict(t))
            ci = kmf.confidence_interval_survival_function_
            idx = ci.index.asof(t)
            milestones[f"{t:g}mo"] = {
                "survival": s,
                "ci": [float(ci.loc[idx].iloc[0]), float(ci.loc[idx].iloc[1])],
                "n_at_risk_at_t": at_risk,
            }
    out = {
        "label": label,
        "n": int(len(df)),
        "events": int(df["event"].sum()),
        "censored": int((df["event"] == 0).sum()),
        "median_pfs_months": float(med),
        "median_pfs_ci": [lo, hi],
        "median_followup_reverse_km_months": float(rev.median_survival_time_),
        "max_observed_months": float(df["time_months"].max()),
        "milestone_pfs": milestones,
        "events_needed_hr0.50": float(events_needed(0.50)),
        "events_needed_hr0.67": float(events_needed(0.67)),
    }
    say(f"### A. {label}")
    say(
        f"- n = {out['n']}, events = {out['events']}, censored = {out['censored']} "
        f"({100 * out['censored'] / out['n']:.0f}% censored)"
    )
    say(
        f"- median PFS = {fmt(med, 2)} months (95% CI {fmt(lo, 2)}-{fmt(hi, 2)}); "
        f"max observed time = {fmt(out['max_observed_months'], 1)} months"
    )
    say(
        f"- median follow-up (reverse Kaplan-Meier) = "
        f"{fmt(out['median_followup_reverse_km_months'], 2)} months"
    )
    for key, val in milestones.items():
        say(
            f"- PFS at {key} = {100 * val['survival']:.1f}% "
            f"(95% CI {100 * val['ci'][0]:.1f}-{100 * val['ci'][1]:.1f}); "
            f"n at risk entering {key} = {val['n_at_risk_at_t']}"
        )
    say(
        f"- Schoenfeld events needed for 80% power, two-sided alpha 0.05, equal split: "
        f"{out['events_needed_hr0.50']:.0f} for HR = 0.50, "
        f"{out['events_needed_hr0.67']:.0f} for HR = 0.67. "
        f"This cohort has {out['events']} events."
    )
    say()
    return out


def cox_continuous(df: pd.DataFrame, score: str, label: str, covariates: list[str] | None = None) -> dict:
    cols = ["time_months", "event", score] + (covariates or [])
    d = df[cols].dropna().copy()
    d[score] = (d[score] - d[score].mean()) / d[score].std(ddof=1)
    cph = CoxPHFitter().fit(d, "time_months", "event")
    row = cph.summary.loc[score]
    n_var = 1 + len(covariates or [])
    epv = int(d["event"].sum()) / n_var
    ci_train = float(cph.concordance_index_)
    y = surv_array(d)
    uno = float(
        concordance_index_ipcw(y, y, cph.predict_partial_hazard(d).to_numpy(), tau=None)[0]
    )
    reduced_cols = ["time_months", "event"] + (covariates or [])
    if covariates:
        ll_reduced = CoxPHFitter().fit(d[reduced_cols], "time_months", "event").log_likelihood_
    else:
        ll_reduced = cph._ll_null_
    lr_stat = 2 * (cph.log_likelihood_ - ll_reduced)
    out = {
        "label": label,
        "score": score,
        "covariates": covariates or [],
        "n": int(len(d)),
        "events": int(d["event"].sum()),
        "events_per_variable": epv,
        "hr_per_sd": float(row["exp(coef)"]),
        "hr_ci": [float(row["exp(coef) lower 95%"]), float(row["exp(coef) upper 95%"])],
        "wald_p": float(row["p"]),
        "partial_lr_p_score": float(chi2.sf(lr_stat, 1)),
        "global_lr_p_all_covariates": float(cph.log_likelihood_ratio_test().p_value),
        "harrell_c_apparent": ci_train,
        "uno_c_ipcw": uno,
    }
    say(f"### B. Continuous score - {label}")
    say(f"- model: {score}{' + ' + ' + '.join(covariates) if covariates else ''}")
    say(f"- events per variable (EPV) = {epv:.1f}")
    say(
        f"- HR per 1 SD = {fmt(out['hr_per_sd'], 2)} "
        f"(95% CI {fmt(out['hr_ci'][0], 2)}-{fmt(out['hr_ci'][1], 2)}), "
        f"Wald p = {fmt(out['wald_p'], 4)}, "
        f"partial LR-test p for the score (1 df) = {fmt(out['partial_lr_p_score'], 4)}"
    )
    say(
        f"- Harrell C (apparent, in-sample) = {fmt(ci_train, 3)}; "
        f"Uno C (IPCW) = {fmt(uno, 3)}"
    )
    say()
    return out


def ph_diagnostics(df: pd.DataFrame, score: str, label: str) -> dict:
    d = df[["time_months", "event", score]].dropna().copy()
    d[score] = (d[score] - d[score].mean()) / d[score].std(ddof=1)
    cph = CoxPHFitter().fit(d, "time_months", "event")
    zph = proportional_hazard_test(cph, d, time_transform=["km", "rank", "identity"])
    tests = {}
    for (covariate, transform), row in zph.summary.iterrows():
        tests[f"{covariate}|{transform}"] = {
            "test_statistic": float(row["test_statistic"]),
            "p": float(row["p"]),
        }
    cut = float(np.median(d.loc[d["event"] == 1, "time_months"]))
    early = d.copy()
    early["event"] = np.where(early["time_months"] <= cut, early["event"], 0)
    early["time_months"] = np.minimum(early["time_months"], cut)
    late = d[d["time_months"] > cut].copy()
    piece = {"split_month": cut}
    for name, part in (("early", early), ("late", late)):
        if part["event"].sum() >= 3 and part[score].std(ddof=1) > 0:
            if name == "late":
                part = part.assign(entry=cut)
                m = CoxPHFitter().fit(
                    part[["entry", "time_months", "event", score]],
                    "time_months",
                    "event",
                    entry_col="entry",
                )
            else:
                m = CoxPHFitter().fit(part[["time_months", "event", score]], "time_months", "event")
            r = m.summary.loc[score]
            piece[name] = {
                "events": int(part["event"].sum()),
                "hr_per_sd": float(r["exp(coef)"]),
                "hr_ci": [float(r["exp(coef) lower 95%"]), float(r["exp(coef) upper 95%"])],
            }
        else:
            piece[name] = {"events": int(part["event"].sum()), "hr_per_sd": None, "hr_ci": None}

    tau = float(min(d.groupby(d[score] > d[score].median())["time_months"].max()))
    hi = d[d[score] > d[score].median()]
    lo = d[d[score] <= d[score].median()]
    rmst_hi = restricted_mean_survival_time(
        KaplanMeierFitter().fit(hi["time_months"], hi["event"]), t=tau
    )
    rmst_lo = restricted_mean_survival_time(
        KaplanMeierFitter().fit(lo["time_months"], lo["event"]), t=tau
    )
    rng = np.random.default_rng(SEED)
    boot = []
    for _ in range(1000):
        bh = hi.sample(len(hi), replace=True, random_state=int(rng.integers(1 << 31)))
        bl = lo.sample(len(lo), replace=True, random_state=int(rng.integers(1 << 31)))
        try:
            boot.append(
                restricted_mean_survival_time(
                    KaplanMeierFitter().fit(bh["time_months"], bh["event"]), t=tau
                )
                - restricted_mean_survival_time(
                    KaplanMeierFitter().fit(bl["time_months"], bl["event"]), t=tau
                )
            )
        except Exception:
            continue
    rmst_ci = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))] if boot else [None, None]
    out = {
        "label": label,
        "schoenfeld": tests,
        "piecewise": piece,
        "rmst": {
            "tau_months": tau,
            "rmst_high_months": float(rmst_hi),
            "rmst_low_months": float(rmst_lo),
            "difference_months": float(rmst_hi - rmst_lo),
            "difference_boot_ci": rmst_ci,
            "n_bootstrap": len(boot),
        },
    }
    say(f"### C. PH diagnostics and RMST - {label}")
    for key, val in tests.items():
        say(
            f"- Schoenfeld ({key}): chi2 = {fmt(val['test_statistic'], 2)}, p = {fmt(val['p'], 3)}"
        )
    say(f"- piecewise HR split at the median event time ({fmt(cut, 2)} months):")
    for name in ("early", "late"):
        p = piece[name]
        if p["hr_per_sd"] is None:
            say(f"  - {name}: {p['events']} events, not estimable")
        else:
            say(
                f"  - {name}: {p['events']} events, HR per SD = {fmt(p['hr_per_sd'], 2)} "
                f"(95% CI {fmt(p['hr_ci'][0], 2)}-{fmt(p['hr_ci'][1], 2)})"
            )
    say(
        f"- RMST to tau = {fmt(tau, 2)} months, median split: "
        f"high = {fmt(rmst_hi, 2)} mo, low = {fmt(rmst_lo, 2)} mo, "
        f"difference = {fmt(rmst_hi - rmst_lo, 2)} mo "
        f"(bootstrap 95% CI {fmt(rmst_ci[0], 2)} to {fmt(rmst_ci[1], 2)}, {len(boot)} resamples)"
    )
    say()
    return out


def candidate_cuts(values: np.ndarray, min_frac: float) -> np.ndarray:
    lo, hi = np.quantile(values, [min_frac, 1 - min_frac])
    uniq = np.unique(values)
    return uniq[(uniq >= lo) & (uniq < hi)]


def max_logrank(time: np.ndarray, event: np.ndarray, score: np.ndarray, min_frac: float):
    best = (-np.inf, np.nan)
    for cut in candidate_cuts(score, min_frac):
        grp = score > cut
        if grp.sum() < 2 or (~grp).sum() < 2:
            continue
        res = logrank_test(time[grp], time[~grp], event[grp], event[~grp])
        if res.test_statistic > best[0]:
            best = (float(res.test_statistic), float(cut))
    return best


def cutpoint_analysis(df: pd.DataFrame, score: str, label: str, min_frac: float = 0.25) -> dict:
    d = df[["time_months", "event", score]].dropna()
    time = d["time_months"].to_numpy()
    event = d["event"].to_numpy()
    s = d[score].to_numpy()

    med = float(np.median(s))
    grp = s > med
    lr_med = logrank_test(time[grp], time[~grp], event[grp], event[~grp])
    cph_med = CoxPHFitter().fit(
        pd.DataFrame({"time_months": time, "event": event, "high": grp.astype(int)}),
        "time_months",
        "event",
    )
    row_med = cph_med.summary.loc["high"]

    stat_opt, cut_opt = max_logrank(time, event, s, min_frac)
    grp_opt = s > cut_opt
    lr_opt = logrank_test(time[grp_opt], time[~grp_opt], event[grp_opt], event[~grp_opt])
    cph_opt = CoxPHFitter().fit(
        pd.DataFrame({"time_months": time, "event": event, "high": grp_opt.astype(int)}),
        "time_months",
        "event",
    )
    row_opt = cph_opt.summary.loc["high"]

    rng = np.random.default_rng(SEED)
    n_perm = 1000
    null = np.empty(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(len(s))
        null[i] = max_logrank(time, event, s[perm], min_frac)[0]
    p_perm = float((1 + np.sum(null >= stat_opt)) / (1 + n_perm))

    boot_cuts, boot_hr = [], []
    for _ in range(500):
        idx = rng.integers(0, len(s), len(s))
        if event[idx].sum() < 5 or np.unique(s[idx]).size < 5:
            continue
        _, c = max_logrank(time[idx], event[idx], s[idx], min_frac)
        if not np.isfinite(c):
            continue
        boot_cuts.append(float((s < c).mean()))
        g = s[idx] > c
        if g.sum() >= 2 and (~g).sum() >= 2:
            try:
                m = CoxPHFitter().fit(
                    pd.DataFrame(
                        {"time_months": time[idx], "event": event[idx], "high": g.astype(int)}
                    ),
                    "time_months",
                    "event",
                )
                boot_hr.append(float(m.summary.loc["high", "exp(coef)"]))
            except Exception:
                pass
    boot_cuts = np.array(boot_cuts) if boot_cuts else np.array([np.nan])
    pct_opt = float((s < cut_opt).mean())
    out = {
        "label": label,
        "min_group_fraction": min_frac,
        "median_cut": {
            "cut": med,
            "n_high": int(grp.sum()),
            "n_low": int((~grp).sum()),
            "logrank_chi2": float(lr_med.test_statistic),
            "logrank_p": float(lr_med.p_value),
            "hr_high_vs_low": float(row_med["exp(coef)"]),
            "hr_ci": [
                float(row_med["exp(coef) lower 95%"]),
                float(row_med["exp(coef) upper 95%"]),
            ],
        },
        "optimal_cut": {
            "cut": float(cut_opt),
            "cut_percentile": pct_opt,
            "n_high": int(grp_opt.sum()),
            "n_low": int((~grp_opt).sum()),
            "n_candidate_cuts": int(len(candidate_cuts(s, min_frac))),
            "logrank_chi2": float(lr_opt.test_statistic),
            "logrank_p_naive": float(lr_opt.p_value),
            "logrank_p_permutation": p_perm,
            "n_permutations": n_perm,
            "hr_high_vs_low": float(row_opt["exp(coef)"]),
            "hr_ci": [
                float(row_opt["exp(coef) lower 95%"]),
                float(row_opt["exp(coef) upper 95%"]),
            ],
        },
        "bootstrap_cut_stability": {
            "n_resamples_used": int(np.isfinite(boot_cuts).sum()),
            "selected_percentile_median": float(np.nanmedian(boot_cuts)),
            "selected_percentile_iqr": [
                float(np.nanpercentile(boot_cuts, 25)),
                float(np.nanpercentile(boot_cuts, 75)),
            ],
            "selected_percentile_range": [
                float(np.nanmin(boot_cuts)),
                float(np.nanmax(boot_cuts)),
            ],
            "fraction_within_10pct_of_original": float(
                np.nanmean(np.abs(boot_cuts - pct_opt) <= 0.10)
            ),
            "hr_boot_median": float(np.median(boot_hr)) if boot_hr else None,
            "hr_boot_2_5_97_5": (
                [float(np.percentile(boot_hr, 2.5)), float(np.percentile(boot_hr, 97.5))]
                if boot_hr
                else None
            ),
        },
    }
    m, o, b = out["median_cut"], out["optimal_cut"], out["bootstrap_cut_stability"]
    say(f"### D. Median cut vs data-driven cut - {label}")
    say(
        f"- median cut ({m['n_high']} high / {m['n_low']} low): HR = {fmt(m['hr_high_vs_low'], 2)} "
        f"(95% CI {fmt(m['hr_ci'][0], 2)}-{fmt(m['hr_ci'][1], 2)}), log-rank p = {fmt(m['logrank_p'], 4)}"
    )
    say(
        f"- maximally selected cut over {o['n_candidate_cuts']} candidates with >= "
        f"{100 * min_frac:.0f}% per arm: cut at the {100 * o['cut_percentile']:.0f}th percentile "
        f"({o['n_high']} high / {o['n_low']} low), HR = {fmt(o['hr_high_vs_low'], 2)} "
        f"(95% CI {fmt(o['hr_ci'][0], 2)}-{fmt(o['hr_ci'][1], 2)})"
    )
    say(
        f"- log-rank p for that cut: naive = {fmt(o['logrank_p_naive'], 4)} vs "
        f"permutation-corrected = {fmt(o['logrank_p_permutation'], 4)} "
        f"({o['n_permutations']} permutations of the score)"
    )
    say(
        f"- bootstrap ({b['n_resamples_used']} resamples): selected cut moves over percentiles "
        f"{100 * b['selected_percentile_range'][0]:.0f}-{100 * b['selected_percentile_range'][1]:.0f} "
        f"(IQR {100 * b['selected_percentile_iqr'][0]:.0f}-{100 * b['selected_percentile_iqr'][1]:.0f}); "
        f"{100 * b['fraction_within_10pct_of_original']:.0f}% of resamples land within "
        f"10 percentile points of the original cut"
    )
    if b["hr_boot_median"] is not None:
        say(
            f"- bootstrap HR for the re-selected cut: median {fmt(b['hr_boot_median'], 2)} "
            f"(2.5-97.5th percentile {fmt(b['hr_boot_2_5_97_5'][0], 2)}-{fmt(b['hr_boot_2_5_97_5'][1], 2)})"
        )
    say()
    return out


def time_dependent_auc(train: pd.DataFrame, score: str, label: str, times=(3.0, 6.0, 12.0)) -> dict:
    d = train[["time_months", "event", score]].dropna()
    y = surv_array(d)
    cph = CoxPHFitter().fit(d[["time_months", "event", score]], "time_months", "event")
    risk = cph.predict_partial_hazard(d).to_numpy()
    raw = d[score].to_numpy()
    upper = float(d.loc[d["event"] == 1, "time_months"].max())
    usable = [t for t in times if t < upper]
    auc, mean_auc = cumulative_dynamic_auc(y, y, risk, np.array(usable))
    auc_raw, _ = cumulative_dynamic_auc(y, y, raw, np.array(usable))
    out = {
        "label": label,
        "times_months": usable,
        "times_dropped": [t for t in times if t not in usable],
        "risk_score": "Cox linear predictor (correct orientation)",
        "auc": [float(a) for a in np.atleast_1d(auc)],
        "auc_if_raw_score_used_as_risk": [float(a) for a in np.atleast_1d(auc_raw)],
        "mean_auc": float(mean_auc),
        "cox_beta_sign": float(cph.summary.loc[score, "coef"]),
        "n_at_risk": {f"{t:g}": int((d["time_months"] >= t).sum()) for t in usable},
        "n_events_by": {
            f"{t:g}": int(((d["time_months"] <= t) & (d["event"] == 1)).sum()) for t in usable
        },
    }
    say(f"### E. Time-dependent (cumulative/dynamic) AUC, IPCW - {label}")
    for t, a, ar in zip(usable, np.atleast_1d(auc), np.atleast_1d(auc_raw)):
        say(
            f"- AUC({t:g} months) = {fmt(float(a), 3)} using the Cox linear predictor, vs "
            f"{fmt(float(ar), 3)} if the raw marker is passed as a risk score "
            f"(events by {t:g} mo = {out['n_events_by'][f'{t:g}']}, "
            f"still at risk = {out['n_at_risk'][f'{t:g}']})"
        )
    say(f"- integrated mean AUC over {usable} = {fmt(mean_auc, 3)}")
    say(f"- fitted log-hazard ratio = {fmt(out['cox_beta_sign'], 3)}")
    if out["times_dropped"]:
        say(
            f"- dropped time points beyond the last event ({fmt(upper, 1)} months): "
            f"{out['times_dropped']} — IPCW weights are undefined there"
        )
    say()
    return out


def dcb_vs_pfs(df: pd.DataFrame, score: str, label: str, horizon: float = 6.0) -> dict:
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from scipy import stats as sps

    d = df[["time_months", "event", score]].dropna().copy()
    z = (d[score] - d[score].mean()) / d[score].std(ddof=1)
    d["z"] = z
    d["dcb"] = np.where(
        (d["time_months"] >= horizon), 1, np.where(d["event"] == 1, 0, np.nan)
    )
    n_unclassifiable = int(d["dcb"].isna().sum())
    evaluable = d.dropna(subset=["dcb"])
    x = evaluable[["z"]].to_numpy()
    yb = evaluable["dcb"].to_numpy().astype(int)
    auc_bin = float(roc_auc_score(yb, x.ravel()))
    mw = sps.mannwhitneyu(
        evaluable.loc[yb == 1, "z"], evaluable.loc[yb == 0, "z"], alternative="two-sided"
    )
    cph_full = CoxPHFitter().fit(d[["time_months", "event", "z"]], "time_months", "event")
    landmark = d[d["time_months"] >= horizon].copy()
    landmark["time_months"] = landmark["time_months"] - horizon
    landmark = landmark[landmark["time_months"] > 0]
    lm: dict[str, object] = {"n": int(len(landmark)), "events": int(landmark["event"].sum())}
    if lm["events"] >= 3 and len(landmark) >= 6:
        m = CoxPHFitter().fit(landmark[["time_months", "event", "z"]], "time_months", "event")
        lm["hr_per_sd"] = float(m.summary.loc["z", "exp(coef)"])
        lm["hr_ci"] = [
            float(m.summary.loc["z", "exp(coef) lower 95%"]),
            float(m.summary.loc["z", "exp(coef) upper 95%"]),
        ]
        lm["p"] = float(m.summary.loc["z", "p"])
    out = {
        "label": label,
        "horizon_months": horizon,
        "n_total": int(len(d)),
        "n_dcb": int((d["dcb"] == 1).sum()),
        "n_no_dcb": int((d["dcb"] == 0).sum()),
        "n_unclassifiable_censored_before_horizon": n_unclassifiable,
        "binary_auc_apparent": auc_bin,
        "binary_mannwhitney_p": float(mw.pvalue),
        "cox_all_patients": {
            "n": int(len(d)),
            "events": int(d["event"].sum()),
            "hr_per_sd": float(cph_full.summary.loc["z", "exp(coef)"]),
            "p": float(cph_full.summary.loc["z", "p"]),
        },
        "landmark_at_horizon": lm,
    }
    say(f"### F. DCB >= {horizon:g} months (binary) vs PFS (time-to-event) - {label}")
    say(
        f"- PFS-based DCB proxy: {out['n_dcb']} with DCB, {out['n_no_dcb']} without, "
        f"{n_unclassifiable} not classifiable because censoring occurred before "
        f"{horizon:g} months"
    )
    say(
        f"- binary analysis on evaluable patients: apparent AUC = {fmt(auc_bin, 3)}, "
        f"Mann-Whitney p = {fmt(out['binary_mannwhitney_p'], 4)}"
    )
    say(
        f"- Cox on all {out['cox_all_patients']['n']} patients "
        f"({out['cox_all_patients']['events']} events): HR per SD = "
        f"{fmt(out['cox_all_patients']['hr_per_sd'], 2)}, p = {fmt(out['cox_all_patients']['p'], 4)}"
    )
    if n_unclassifiable == 0:
        say(
            f"- every censored patient was followed past {horizon:g} months, so nobody was "
            "excluded and the guarantee-time trap does not bite here; the exclusion count "
            "is exactly the number you must report"
        )
    if "hr_per_sd" in lm:
        say(
            f"- landmark analysis restarting the clock at {horizon:g} months "
            f"(n = {lm['n']}, {lm['events']} subsequent events): HR per SD = "
            f"{fmt(lm['hr_per_sd'], 2)} (95% CI {fmt(lm['hr_ci'][0], 2)}-{fmt(lm['hr_ci'][1], 2)}), "
            f"p = {fmt(lm['p'], 3)}"
        )
    else:
        say(
            f"- landmark analysis at {horizon:g} months is not estimable: only {lm['n']} patients "
            f"reach the landmark with {lm['events']} subsequent events"
        )
    say()
    return out


def transfer(disc: pd.DataFrame, val: pd.DataFrame, score_disc: str, score_val: str, cut: float, label: str) -> dict:
    def hr_for(df: pd.DataFrame, high: np.ndarray) -> dict:
        if high.sum() < 2 or (~high).sum() < 2:
            return {"n_high": int(high.sum()), "n_low": int((~high).sum()), "hr": None}
        m = CoxPHFitter().fit(
            pd.DataFrame(
                {
                    "time_months": df["time_months"].to_numpy(),
                    "event": df["event"].to_numpy(),
                    "high": high.astype(int),
                }
            ),
            "time_months",
            "event",
        )
        r = m.summary.loc["high"]
        lr = logrank_test(
            df["time_months"][high],
            df["time_months"][~high],
            df["event"][high],
            df["event"][~high],
        )
        return {
            "n_high": int(high.sum()),
            "n_low": int((~high).sum()),
            "hr": float(r["exp(coef)"]),
            "hr_ci": [float(r["exp(coef) lower 95%"]), float(r["exp(coef) upper 95%"])],
            "logrank_p": float(lr.p_value),
        }

    val_ref = val[score_val].to_numpy()
    out = {
        "label": label,
        "discovery_cut_value": float(cut),
        "validation_fixed_cut": hr_for(val, val_ref > cut),
        "validation_own_median": hr_for(val, val_ref > np.median(val_ref)),
        "score_ref_distribution": {
            "discovery_mean": float(disc[score_disc].mean()),
            "discovery_sd": float(disc[score_disc].std(ddof=1)),
            "validation_mean": float(np.mean(val_ref)),
            "validation_sd": float(np.std(val_ref, ddof=1)),
            "fraction_above_discovery_cut": float(np.mean(val_ref > cut)),
        },
    }
    say(f"### G. Moving the cut to an independent cohort - {label}")
    dist = out["score_ref_distribution"]
    say(
        f"- score rescaled with the discovery cohort's gene-wise mean/SD: validation mean = "
        f"{fmt(dist['validation_mean'], 2)}, SD = {fmt(dist['validation_sd'], 2)} "
        f"(discovery is 0 / 1 by construction)"
    )
    say(
        f"- applying the discovery cut ({fmt(cut, 3)}) verbatim puts "
        f"{100 * dist['fraction_above_discovery_cut']:.0f}% of the validation cohort in the high group"
    )
    for key in ("validation_fixed_cut", "validation_own_median"):
        v = out[key]
        if v["hr"] is None:
            say(f"- {key}: {v['n_high']} high / {v['n_low']} low — not estimable")
        else:
            say(
                f"- {key}: {v['n_high']} high / {v['n_low']} low, HR = {fmt(v['hr'], 2)} "
                f"(95% CI {fmt(v['hr_ci'][0], 2)}-{fmt(v['hr_ci'][1], 2)}), "
                f"log-rank p = {fmt(v['logrank_p'], 3)}"
            )
    say()
    return out


def competing_risks_demo(n: int = 400) -> dict:
    rng = np.random.default_rng(SEED)
    t1 = rng.exponential(12.0, n)
    t2 = rng.exponential(18.0, n)
    cens = rng.exponential(30.0, n)
    time = np.minimum(np.minimum(t1, t2), cens)
    cause = np.where(cens <= np.minimum(t1, t2), 0, np.where(t1 <= t2, 1, 2))
    horizon = 12.0
    aj = AalenJohansenFitter(calculate_variance=False, seed=SEED)
    aj.fit(time, cause, event_of_interest=1)
    cif = aj.cumulative_density_
    cif_at = float(cif.loc[: cif.index.asof(horizon)].iloc[-1, 0])
    km_naive = KaplanMeierFitter().fit(time, (cause == 1).astype(int))
    one_minus_km = float(1 - km_naive.predict(horizon))
    out = {
        "note": "simulated data with a known data-generating process, not GEO",
        "n": n,
        "horizon_months": horizon,
        "n_cause1_progression": int((cause == 1).sum()),
        "n_cause2_competing_death": int((cause == 2).sum()),
        "n_censored": int((cause == 0).sum()),
        "aalen_johansen_cif_cause1": cif_at,
        "one_minus_km_censoring_competing_events": one_minus_km,
        "absolute_overestimation": one_minus_km - cif_at,
    }
    say("### H. Competing risks (simulated, known truth)")
    say(
        f"- n = {n}: {out['n_cause1_progression']} progressions, "
        f"{out['n_cause2_competing_death']} competing deaths, {out['n_censored']} censored"
    )
    say(
        f"- risk of progression by {horizon:g} months: Aalen-Johansen CIF = "
        f"{100 * cif_at:.1f}% vs 1 − Kaplan-Meier (competing deaths censored) = "
        f"{100 * one_minus_km:.1f}%, i.e. an absolute overestimate of "
        f"{100 * out['absolute_overestimation']:.1f} percentage points"
    )
    say(
        "- neither GSE135222 nor GSE190265 deposits a cause-of-event field, so this "
        "contrast cannot be computed on the GEO files"
    )
    say()
    return out


def figures(disc: pd.DataFrame) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, gene in zip(axes, PRIMARY):
        s = disc[f"{gene}_z"].to_numpy()
        t = disc["time_months"].to_numpy()
        e = disc["event"].to_numpy()
        high = s > np.median(s)
        for mask, lab in ((high, "high"), (~high, "low")):
            KaplanMeierFitter().fit(t[mask], e[mask], label=f"{lab} (n={mask.sum()})").plot_survival_function(
                ax=ax, ci_show=True
            )
        lr = logrank_test(t[high], t[~high], e[high], e[~high])
        ax.set_title(f"GSE135222 PFS, {gene} median split\nlog-rank p = {lr.p_value:.3f}")
        ax.set_xlabel("months")
        ax.set_ylabel("progression-free")
        ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(RESULTS / "fig_km_tacstd2_cldn4_median.png", dpi=110)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    t = disc["time_months"].to_numpy()
    e = disc["event"].to_numpy()
    for ax, gene in zip(axes, PRIMARY):
        s = disc[f"{gene}_z"].to_numpy()
        stats = []
        for cut in candidate_cuts(s, 0.25):
            g = s > cut
            stats.append(
                (float((s < cut).mean() * 100), logrank_test(t[g], t[~g], e[g], e[~g]).test_statistic)
            )
        if not stats:
            continue
        xs, ys = zip(*stats)
        ax.plot(xs, ys, marker="o", lw=1)
        ax.axhline(3.841, ls="--", c="grey", label="chi2 = 3.84 (nominal p = 0.05)")
        ax.set_xlabel("cut position (percentile)")
        ax.set_ylabel("log-rank chi-square")
        ax.set_title(f"Every candidate cut is a separate test\n{gene}, GSE135222")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "fig_cutpoint_profile.png", dpi=110)
    plt.close(fig)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    disc, disc_stats = load("gse135222")
    val, _ = load("gse190265", ref_stats=disc_stats)

    say("# Worked example: TACSTD2 / CLDN4 vs PFS on two public GEO ICI cohorts")
    say()
    say(
        "Generated by `02_analysis.py`. Discovery = GSE135222 (Jung et al., NSCLC, "
        "anti-PD-1/PD-L1, PFS in days). Validation = GSE190265 (Limagne / Ghiringhelli, "
        "NSCLC, anti-PD-1, PFS in months). Neither series deposits OS. Primary genes are "
        "TACSTD2 and CLDN4 as log2(TPM+1); CD8A is a known-direction comparator, not a "
        "TACSTD2/CLDN4 claim. The mean-z of TACSTD2+CLDN4 is exploratory."
    )
    say()

    R["provenance"] = json.loads((DATA / "provenance.json").read_text())
    R["describe_discovery"] = describe(disc, "GSE135222 (discovery)")
    R["describe_validation"] = describe(val, "GSE190265 (validation)")

    R["cox"] = {}
    R["ph"] = {}
    R["cuts"] = {}
    R["tdauc"] = {}
    R["dcb"] = {}
    R["transfer"] = {}

    scores = [(g, f"{g}_z") for g in GENES] + [("TJ_mean", "TJ_mean_z")]
    for name, col in scores:
        R["cox"][f"disc_{name}"] = cox_continuous(disc, col, f"GSE135222 {name}")
        R["ph"][f"disc_{name}"] = ph_diagnostics(disc, col, f"GSE135222 {name}")
        R["cuts"][f"disc_{name}"] = cutpoint_analysis(disc, col, f"GSE135222 {name}")
        R["tdauc"][f"disc_{name}"] = time_dependent_auc(disc, col, f"GSE135222 {name}")
        R["dcb"][f"disc_{name}"] = dcb_vs_pfs(disc, col, f"GSE135222 {name}")
        R["cox"][f"val_{name}"] = cox_continuous(val, f"{col}_ref" if name != "TJ_mean" else "TJ_mean_z_ref", f"GSE190265 {name} (discovery-scaled)")
        R["dcb"][f"val_{name}"] = dcb_vs_pfs(val, f"{col}_ref" if name != "TJ_mean" else "TJ_mean_z_ref", f"GSE190265 {name}")
        R["transfer"][name] = transfer(
            disc,
            val,
            col,
            f"{col}_ref" if name != "TJ_mean" else "TJ_mean_z_ref",
            R["cuts"][f"disc_{name}"]["optimal_cut"]["cut"],
            f"{name}: GSE135222 → GSE190265",
        )

    if "age" in disc.columns and "sex" in disc.columns:
        R["cox"]["disc_TACSTD2_adj"] = cox_continuous(
            disc.assign(male=(disc["sex"] == "male").astype(int)),
            "TACSTD2_z",
            "GSE135222 TACSTD2 adjusted for age and sex",
            covariates=["age", "male"],
        )
        R["cox"]["disc_CLDN4_adj"] = cox_continuous(
            disc.assign(male=(disc["sex"] == "male").astype(int)),
            "CLDN4_z",
            "GSE135222 CLDN4 adjusted for age and sex",
            covariates=["age", "male"],
        )

    R["competing_risks_simulation"] = competing_risks_demo()
    figures(disc)
    say("Figures: `fig_km_tacstd2_cldn4_median.png`, `fig_cutpoint_profile.png`.")
    say()
    say(
        f"Reminder: with {R['describe_discovery']['events']} and "
        f"{R['describe_validation']['events']} events these cohorts can illustrate how the "
        "methods behave but cannot establish TACSTD2 or CLDN4 as ICI biomarkers. "
        "Read `../playbook.md` before quoting any number above."
    )

    (RESULTS / "demo_results.md").write_text("\n".join(LINES) + "\n")
    (RESULTS / "demo_results.json").write_text(json.dumps(R, indent=2, default=float) + "\n")
    print(f"\nwrote {RESULTS / 'demo_results.md'}")


if __name__ == "__main__":
    main()
