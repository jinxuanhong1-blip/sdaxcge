#!/usr/bin/env python3
"""B5_BLCA: does CLDN4-high predict worse ICI response in urothelial carcinoma?

Executes the prespecified plan in results/w200/B5_BLCA/analysis_plan.md against
three public urothelial ICI cohorts (IMvigor210, BACI/GSE176307, Snyder 2017).

The user-asserted effect size is fixed as a constant before any computation and is
never used to select a cutoff, cohort, filter or covariate set:

    USER_OR = 0.42

Does not tune any specification to match a pre-specified user number.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats.contingency import odds_ratio

# --------------------------------------------------------------------------- #
# Prespecified constants (fixed before computation; see analysis_plan.md)
# --------------------------------------------------------------------------- #
USER_OR = 0.42
EXACT_DECIMALS = 2
NEARBY_ABS = 0.05
DID_WE_TUNE = False

PRIMARY_GENE = "CLDN4"
POOLED_COHORTS = ["IMvigor210", "BACI", "Snyder"]   # IMvigor210_PredictIO excluded
CONCORDANCE_COHORT = "IMvigor210_PredictIO"
EXPLORATORY_GENES = ["CLDN1", "CLDN2", "CLDN3", "CLDN7", "CLDN18", "TACSTD2",
                     "EPCAM", "CDH1", "OCLN", "TJP1", "CD274"]
SEED = 20260816


# --------------------------------------------------------------------------- #
# 2x2 machinery
# --------------------------------------------------------------------------- #
def table_2x2(high, resp):
    """Rows = [CLDN4-high, CLDN4-low], cols = [responder, non-responder]."""
    high = np.asarray(high, dtype=bool)
    resp = np.asarray(resp, dtype=bool)
    a = int(np.sum(high & resp))        # high, responder
    b = int(np.sum(high & ~resp))       # high, non-responder
    c = int(np.sum(~high & resp))       # low, responder
    d = int(np.sum(~high & ~resp))      # low, non-responder
    return [[a, b], [c, d]]


def fisher_block(tab):
    """Conditional-MLE odds ratio of response for high vs low, with exact CI."""
    (a, b), (c, d) = tab
    _, p = stats.fisher_exact(tab, alternative="two-sided")
    res = odds_ratio(tab, kind="conditional")
    ci = res.confidence_interval(confidence_level=0.95)
    # Woolf log-OR and SE for inverse-variance pooling; Haldane-Anscombe 0.5 only
    # when a cell is zero.
    if min(a, b, c, d) == 0:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        corrected = True
    else:
        aa, bb, cc, dd = float(a), float(b), float(c), float(d)
        corrected = False
    log_or = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return {
        "n": a + b + c + d,
        "high_responder": a, "high_nonresponder": b,
        "low_responder": c, "low_nonresponder": d,
        "orr_high_pct": 100.0 * a / (a + b) if (a + b) else float("nan"),
        "orr_low_pct": 100.0 * c / (c + d) if (c + d) else float("nan"),
        "or_conditional_mle": float(res.statistic),
        "ci_low": float(ci.low), "ci_high": float(ci.high),
        "fisher_p": float(p),
        "log_or_woolf": log_or, "se_log_or": se,
        "haldane_corrected": corrected,
    }


def split_high(values, rule="median"):
    """Response-blind dichotomisation of the exposure. Ties go to the low group."""
    v = np.asarray(values, dtype=float)
    if rule == "median":
        thr = float(np.median(v))
        return v > thr, {"rule": "median", "threshold": thr}, np.ones(len(v), dtype=bool)
    if rule == "tertile_T3_vs_T1":
        lo, hi = np.quantile(v, [1 / 3, 2 / 3])
        keep = (v <= lo) | (v > hi)
        return v > hi, {"rule": rule, "lower": float(lo), "upper": float(hi)}, keep
    if rule == "quartile_Q4_vs_Q1":
        lo, hi = np.quantile(v, [0.25, 0.75])
        keep = (v <= lo) | (v > hi)
        return v > hi, {"rule": rule, "lower": float(lo), "upper": float(hi)}, keep
    if rule == "upper_quartile_vs_rest":
        hi = float(np.quantile(v, 0.75))
        return v > hi, {"rule": rule, "threshold": hi}, np.ones(len(v), dtype=bool)
    if rule == "split_60_40":
        thr = float(np.quantile(v, 0.60))
        return v > thr, {"rule": rule, "threshold": thr}, np.ones(len(v), dtype=bool)
    raise ValueError(rule)


# --------------------------------------------------------------------------- #
# Meta-analysis
# --------------------------------------------------------------------------- #
def meta_analyze(effects):
    """Pool log ORs. Random effects (DerSimonian-Laird) is the prespecified primary;
    Hartung-Knapp limits are reported because k is small."""
    y = np.array([e["log_or_woolf"] for e in effects], dtype=float)
    se = np.array([e["se_log_or"] for e in effects], dtype=float)
    k = len(y)
    w = 1.0 / se ** 2

    fe = float(np.sum(w * y) / np.sum(w))
    se_fe = float(math.sqrt(1.0 / np.sum(w)))

    Q = float(np.sum(w * (y - fe) ** 2))
    df = k - 1
    C = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    I2 = max(0.0, (Q - df) / Q * 100.0) if Q > 0 else 0.0
    Q_p = float(stats.chi2.sf(Q, df)) if df > 0 else float("nan")

    w_re = 1.0 / (se ** 2 + tau2)
    re = float(np.sum(w_re * y) / np.sum(w_re))
    se_re = float(math.sqrt(1.0 / np.sum(w_re)))

    # Hartung-Knapp: rescale the variance and use a t reference with k-1 df.
    if k > 1:
        q_hk = float(np.sum(w_re * (y - re) ** 2) / (k - 1) / np.sum(w_re))
        se_hk = math.sqrt(max(q_hk, 0.0))
        tcrit = float(stats.t.ppf(0.975, k - 1))
        hk_lo, hk_hi = re - tcrit * se_hk, re + tcrit * se_hk
        hk_p = float(2 * stats.t.sf(abs(re / se_hk), k - 1)) if se_hk > 0 else float("nan")
    else:
        se_hk = float("nan"); hk_lo = hk_hi = hk_p = float("nan")

    def ci(est, s):
        return est - 1.96 * s, est + 1.96 * s

    fe_lo, fe_hi = ci(fe, se_fe)
    re_lo, re_hi = ci(re, se_re)

    # Mantel-Haenszel with Robins-Breslow-Greenland variance
    num = den = 0.0
    P = Qs = R = S = 0.0
    for e in effects:
        a, b = e["high_responder"], e["high_nonresponder"]
        c, d = e["low_responder"], e["low_nonresponder"]
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n
        den += b * c / n
        P += (a + d) / n
        Qs += (b + c) / n
        R += a * d / n
        S += b * c / n
    mh = num / den if den > 0 else float("nan")
    if R > 0 and S > 0:
        var_mh = 0.0
        sp = sq = sr = 0.0
        for e in effects:
            a, b = e["high_responder"], e["high_nonresponder"]
            c, d = e["low_responder"], e["low_nonresponder"]
            n = a + b + c + d
            if n == 0:
                continue
            p_i, q_i = (a + d) / n, (b + c) / n
            r_i, s_i = a * d / n, b * c / n
            sp += p_i * r_i
            sq += p_i * s_i + q_i * r_i
            sr += q_i * s_i
        var_mh = sp / (2 * R ** 2) + sq / (2 * R * S) + sr / (2 * S ** 2)
        se_mh = math.sqrt(var_mh)
        mh_lo, mh_hi = math.exp(math.log(mh) - 1.96 * se_mh), math.exp(math.log(mh) + 1.96 * se_mh)
    else:
        se_mh = float("nan"); mh_lo = mh_hi = float("nan")

    return {
        "k": k,
        "fixed_iv": {"or": math.exp(fe), "ci_low": math.exp(fe_lo), "ci_high": math.exp(fe_hi),
                     "log_or": fe, "se": se_fe,
                     "p": float(2 * stats.norm.sf(abs(fe / se_fe)))},
        "random_dl": {"or": math.exp(re), "ci_low": math.exp(re_lo), "ci_high": math.exp(re_hi),
                      "log_or": re, "se": se_re,
                      "p": float(2 * stats.norm.sf(abs(re / se_re)))},
        "random_dl_hartung_knapp": {"or": math.exp(re), "ci_low": math.exp(hk_lo),
                                    "ci_high": math.exp(hk_hi), "se": se_hk, "p": hk_p},
        "mantel_haenszel": {"or": mh, "ci_low": mh_lo, "ci_high": mh_hi, "se_log": se_mh},
        "heterogeneity": {"Q": Q, "df": df, "p": Q_p, "I2_percent": I2, "tau2": tau2},
    }


# --------------------------------------------------------------------------- #
def analysis_set(df, cohort, gene=PRIMARY_GENE):
    sub = df[(df["cohort"] == cohort) & df["responder"].notna() & df[gene].notna()].copy()
    sub["responder"] = sub["responder"].astype(int)
    return sub


def continuous_logistic(sub, gene=PRIMARY_GENE):
    """OR per 1 log2 unit and per 1 within-cohort SD, via logistic regression."""
    import statsmodels.api as sm
    x = sub[gene].astype(float).values
    y = sub["responder"].astype(int).values
    sd = float(np.std(x, ddof=1))
    X = sm.add_constant(x)
    try:
        res = sm.Logit(y, X).fit(disp=0)
        beta, se = float(res.params[1]), float(res.bse[1])
        p = float(res.pvalues[1])
    except Exception as exc:  # separation / non-convergence
        return {"error": str(exc), "n": int(len(y)), "sd": sd}
    return {
        "n": int(len(y)), "sd_of_exposure": sd,
        "or_per_log2_unit": math.exp(beta),
        "ci_low_per_log2_unit": math.exp(beta - 1.96 * se),
        "ci_high_per_log2_unit": math.exp(beta + 1.96 * se),
        "or_per_sd": math.exp(beta * sd),
        "ci_low_per_sd": math.exp((beta - 1.96 * se) * sd),
        "ci_high_per_sd": math.exp((beta + 1.96 * se) * sd),
        "p": p,
    }


def mannwhitney(sub, gene=PRIMARY_GENE):
    r = sub.loc[sub["responder"] == 1, gene].astype(float).values
    n = sub.loc[sub["responder"] == 0, gene].astype(float).values
    if len(r) < 2 or len(n) < 2:
        return {"error": "too few in one group", "n_responder": len(r), "n_nonresponder": len(n)}
    u, p = stats.mannwhitneyu(r, n, alternative="two-sided")
    rb = 2 * u / (len(r) * len(n)) - 1  # rank-biserial correlation
    return {
        "n_responder": int(len(r)), "n_nonresponder": int(len(n)),
        "median_responder": float(np.median(r)), "median_nonresponder": float(np.median(n)),
        "u": float(u), "p": float(p), "rank_biserial_r": float(rb),
    }


def test_against_claim(pooled_log_or, se_log_or, claimed_or=USER_OR):
    """Formally test H0: true OR equals the claimed value.

    A confidence interval that excludes the claimed value is a stronger and more
    useful statement than 'we failed to reject 1', so it is reported explicitly.
    """
    z = (pooled_log_or - math.log(claimed_or)) / se_log_or
    p = float(2 * stats.norm.sf(abs(z)))
    return {
        "null_hypothesis": f"true pooled OR = {claimed_or}",
        "z": float(z),
        "p": p,
        "rejected_at_0.05": bool(p < 0.05),
        "interpretation": ("the urothelial data are statistically incompatible with an OR of "
                           f"{claimed_or}" if p < 0.05 else
                           f"the data cannot exclude an OR of {claimed_or}"),
    }


def power_for_claimed_effect(effects, claimed_or=USER_OR, n_sim=20000, seed=SEED):
    """Simulation power: if the claimed effect were real, how often would this
    analysis have detected it?

    This separates 'no evidence of an effect' from 'evidence against an effect'.
    Group sizes and the CLDN4-low response rate are taken from the observed data;
    the CLDN4-high rate is set to whatever the claimed OR implies.
    """
    rng = np.random.default_rng(seed)
    design = []
    for e in effects:
        n_hi = e["high_responder"] + e["high_nonresponder"]
        n_lo = e["low_responder"] + e["low_nonresponder"]
        p_lo = e["low_responder"] / n_lo if n_lo else 0.0
        p_lo = min(max(p_lo, 1e-6), 1 - 1e-6)
        odds_lo = p_lo / (1 - p_lo)
        odds_hi = claimed_or * odds_lo
        p_hi = odds_hi / (1 + odds_hi)
        design.append((n_hi, p_hi, n_lo, p_lo))

    detected = 0
    detected_ci = 0
    for _ in range(n_sim):
        ys, ses = [], []
        for n_hi, p_hi, n_lo, p_lo in design:
            a = rng.binomial(n_hi, p_hi)
            c = rng.binomial(n_lo, p_lo)
            b, d = n_hi - a, n_lo - c
            if min(a, b, c, d) == 0:
                a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
            ys.append(math.log((a * d) / (b * c)))
            ses.append(math.sqrt(1 / a + 1 / b + 1 / c + 1 / d))
        y = np.array(ys); se = np.array(ses)
        w = 1 / se ** 2
        fe = float(np.sum(w * y) / np.sum(w))
        Q = float(np.sum(w * (y - fe) ** 2))
        dfree = len(y) - 1
        C = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
        tau2 = max(0.0, (Q - dfree) / C) if C > 0 else 0.0
        w_re = 1 / (se ** 2 + tau2)
        re = float(np.sum(w_re * y) / np.sum(w_re))
        se_re = math.sqrt(1 / np.sum(w_re))
        p = 2 * stats.norm.sf(abs(re / se_re))
        if p < 0.05:
            detected += 1
            if re < 0:
                detected_ci += 1
    return {
        "claimed_or_simulated": claimed_or,
        "n_simulations": n_sim,
        "power_two_sided_alpha_0.05": detected / n_sim,
        "power_correct_direction": detected_ci / n_sim,
        "design_note": ("group sizes and CLDN4-low response rates taken from the observed "
                        "cohorts; CLDN4-high rate implied by the claimed OR"),
    }


def verdict_block(pooled_or, ci_low, ci_high):
    rounded = round(pooled_or, EXACT_DECIMALS)
    diff = abs(pooled_or - USER_OR)
    matches = rounded == USER_OR
    nearby = diff <= NEARBY_ABS
    if matches:
        label = "MATCH_AT_2DP"
    elif nearby:
        label = "NEARBY_NOT_EXACT"
    else:
        label = "DOES_NOT_MATCH"
    compatible = bool(ci_low <= USER_OR <= ci_high)
    direction = ("same direction (OR<1, CLDN4-high responds less often)" if pooled_or < 1
                 else "opposite direction (OR>1, CLDN4-high responds MORE often)")
    return {
        "stat_compared": "random-effects (DerSimonian-Laird) pooled OR of objective response, CLDN4-high vs CLDN4-low",
        "user_claimed_or": USER_OR,
        "observed_or": pooled_or,
        "rounded_2dp": rounded,
        f"matches_user_{USER_OR}_at_2dp": matches,
        "abs_diff": diff,
        f"within_{NEARBY_ABS}": nearby,
        "label": label,
        "claimed_value_inside_observed_95ci": compatible,
        "direction_vs_claim": direction,
        "did_we_tune_to_0.42": DID_WE_TUNE,
        "statement": "",  # filled in main()
    }


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", default="results/w200/B5_BLCA")
    ap.add_argument("--out-dir", default="results/w200/B5_BLCA")
    args = ap.parse_args()
    np.random.seed(SEED)
    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(os.path.join(args.in_dir, "per_sample_expression.csv"))

    # ===================== PRIMARY ========================================== #
    per_cohort, effects = [], []
    for cohort in POOLED_COHORTS:
        sub = analysis_set(df, cohort)
        high, meta, keep = split_high(sub[PRIMARY_GENE].values, "median")
        tab = table_2x2(high, sub["responder"].values == 1)
        blk = fisher_block(tab)
        blk.update({"cohort": cohort, "cutoff_rule": "median",
                    "threshold": meta.get("threshold"),
                    "expression_scale": sub["expression_scale"].iloc[0]})
        per_cohort.append(blk)
        effects.append(blk)
        print(f"[primary] {cohort}: n={blk['n']} ORR high={blk['orr_high_pct']:.1f}% "
              f"low={blk['orr_low_pct']:.1f}% OR={blk['or_conditional_mle']:.3f} "
              f"[{blk['ci_low']:.3f},{blk['ci_high']:.3f}] p={blk['fisher_p']:.4f}")

    pooled = meta_analyze(effects)
    print(f"[meta] random-effects OR={pooled['random_dl']['or']:.3f} "
          f"[{pooled['random_dl']['ci_low']:.3f},{pooled['random_dl']['ci_high']:.3f}] "
          f"I2={pooled['heterogeneity']['I2_percent']:.1f}% "
          f"tau2={pooled['heterogeneity']['tau2']:.4f}")

    verdict = verdict_block(pooled["random_dl"]["or"],
                            pooled["random_dl"]["ci_low"],
                            pooled["random_dl"]["ci_high"])

    claim_test = test_against_claim(pooled["random_dl"]["log_or"], pooled["random_dl"]["se"])
    power = power_for_claimed_effect(effects)
    verdict["test_against_claimed_value"] = claim_test
    verdict["power_to_detect_claimed_effect"] = power
    print(f"[claim] H0: OR=={USER_OR} -> z={claim_test['z']:.2f} p={claim_test['p']:.4g} "
          f"rejected={claim_test['rejected_at_0.05']}")
    print(f"[power] if OR were {USER_OR}, detection probability = "
          f"{power['power_two_sided_alpha_0.05']:.3f}")

    # ===================== CONTROLS ========================================= #
    controls, controls_continuous = [], []
    for label, col, expectation in (
        ("positive_control_cd8_teff", "cd8_teff_z", "OR>1 expected (published to track response in IMvigor210)"),
        ("negative_control_housekeeping", "housekeeping_z", "OR~1 expected (null)"),
    ):
        ctrl_effects = []
        for cohort in POOLED_COHORTS:
            sub = df[(df["cohort"] == cohort) & df["responder"].notna() & df[col].notna()].copy()
            if len(sub) < 10:
                continue
            sub["responder"] = sub["responder"].astype(int)
            high, meta, _ = split_high(sub[col].values, "median")
            blk = fisher_block(table_2x2(high, sub["responder"].values == 1))
            blk.update({"control": label, "score": col, "cohort": cohort,
                        "expectation": expectation})
            controls.append(blk)
            ctrl_effects.append(blk)
            # A median split of a continuous score is lossy, so also test it
            # continuously: this is the fairer check of whether the pipeline can
            # detect a signal that is known to be present.
            cc = continuous_logistic(sub, gene=col)
            cc.update({"control": label, "score": col, "cohort": cohort})
            controls_continuous.append(cc)
            print(f"[control] {label} {cohort}: median-split OR={blk['or_conditional_mle']:.3f} "
                  f"[{blk['ci_low']:.3f},{blk['ci_high']:.3f}] p={blk['fisher_p']:.4g}"
                  + (f" | per-SD OR={cc['or_per_sd']:.3f} p={cc['p']:.4g}" if "or_per_sd" in cc else ""))
        if ctrl_effects:
            cm = meta_analyze(ctrl_effects)
            controls.append({
                "control": label, "score": col, "cohort": "POOLED_random_effects",
                "expectation": expectation, "n": sum(e["n"] for e in ctrl_effects),
                "or_conditional_mle": cm["random_dl"]["or"],
                "ci_low": cm["random_dl"]["ci_low"], "ci_high": cm["random_dl"]["ci_high"],
                "fisher_p": cm["random_dl"]["p"],
            })
            print(f"[control] {label} POOLED: OR={cm['random_dl']['or']:.3f} "
                  f"[{cm['random_dl']['ci_low']:.3f},{cm['random_dl']['ci_high']:.3f}] "
                  f"p={cm['random_dl']['p']:.4g}")

    # ===================== S1 cutoff robustness ============================== #
    sens_cutoffs = []
    for rule in ["median", "tertile_T3_vs_T1", "quartile_Q4_vs_Q1",
                 "upper_quartile_vs_rest", "split_60_40"]:
        rule_effects = []
        for cohort in POOLED_COHORTS:
            sub = analysis_set(df, cohort)
            high, meta, keep = split_high(sub[PRIMARY_GENE].values, rule)
            h = high[keep]
            r = (sub["responder"].values == 1)[keep]
            blk = fisher_block(table_2x2(h, r))
            blk.update({"cohort": cohort, "cutoff_rule": rule})
            sens_cutoffs.append(blk)
            rule_effects.append(blk)
        pm = meta_analyze(rule_effects)
        sens_cutoffs.append({
            "cohort": "POOLED_random_effects", "cutoff_rule": rule,
            "n": sum(e["n"] for e in rule_effects),
            "or_conditional_mle": pm["random_dl"]["or"],
            "ci_low": pm["random_dl"]["ci_low"], "ci_high": pm["random_dl"]["ci_high"],
            "fisher_p": pm["random_dl"]["p"],
            "I2_percent": pm["heterogeneity"]["I2_percent"],
        })
        print(f"[S1] {rule}: pooled OR={pm['random_dl']['or']:.3f} "
              f"[{pm['random_dl']['ci_low']:.3f},{pm['random_dl']['ci_high']:.3f}]")

    # ===================== S2 continuous, S3 Mann-Whitney =================== #
    cont_rows, mw_rows = [], []
    for cohort in POOLED_COHORTS:
        sub = analysis_set(df, cohort)
        c = continuous_logistic(sub); c["cohort"] = cohort
        cont_rows.append(c)
        m = mannwhitney(sub); m["cohort"] = cohort
        mw_rows.append(m)
        if "or_per_sd" in c:
            print(f"[S2] {cohort}: OR per SD={c['or_per_sd']:.3f} "
                  f"[{c['ci_low_per_sd']:.3f},{c['ci_high_per_sd']:.3f}] p={c['p']:.4g}")

    # ===================== S4 ITT (NE as non-responder) ===================== #
    itt_effects = []
    for cohort in POOLED_COHORTS:
        sub = df[(df["cohort"] == cohort) & df[PRIMARY_GENE].notna()].copy()
        resp = sub["responder"].fillna(0).astype(int)
        high, _, _ = split_high(sub[PRIMARY_GENE].values, "median")
        blk = fisher_block(table_2x2(high, resp.values == 1))
        blk.update({"cohort": cohort, "variant": "S4_ITT_NE_as_nonresponder"})
        itt_effects.append(blk)
    itt_pooled = meta_analyze(itt_effects)

    # ===================== S5 bladder-site only (IMvigor210) ================ #
    sub = analysis_set(df, "IMvigor210")
    bl = sub[sub["tissue"].astype(str).str.lower() == "bladder"]
    high, _, _ = split_high(bl[PRIMARY_GENE].values, "median")
    s5 = fisher_block(table_2x2(high, bl["responder"].values == 1))
    s5.update({"cohort": "IMvigor210", "variant": "S5_bladder_biopsy_only"})

    # ===================== S6 DESeq normalisation =========================== #
    sub6 = df[(df["cohort"] == "IMvigor210") & df["responder"].notna()
              & df["CLDN4_deseq_norm"].notna()].copy()
    high, _, _ = split_high(sub6["CLDN4_deseq_norm"].values, "median")
    s6 = fisher_block(table_2x2(high, sub6["responder"].values == 1))
    s6.update({"cohort": "IMvigor210", "variant": "S6_deseq_size_factor_normalised"})

    # ===================== S7 pipeline concordance ========================== #
    # PredictIO labels the same patients as "P" + the IMvigor210 ANONPT_ID.
    def norm_pid(s):
        return s.astype(str).str.replace(r"^P", "", regex=True)

    own = df[df["cohort"] == "IMvigor210"][["patient_id", PRIMARY_GENE, "responder"]].copy()
    harm = df[df["cohort"] == CONCORDANCE_COHORT][["patient_id", PRIMARY_GENE]].copy()
    own["pid"] = norm_pid(own["patient_id"])
    harm["pid"] = norm_pid(harm["patient_id"])
    merged = own.merge(harm, on="pid", suffixes=("_own", "_predictio")).dropna(
        subset=[f"{PRIMARY_GENE}_own", f"{PRIMARY_GENE}_predictio"])
    if len(merged) > 3:
        rho, rho_p = stats.spearmanr(merged[f"{PRIMARY_GENE}_own"],
                                     merged[f"{PRIMARY_GENE}_predictio"])
        pear = float(np.corrcoef(merged[f"{PRIMARY_GENE}_own"],
                                 merged[f"{PRIMARY_GENE}_predictio"])[0, 1])
    else:
        rho, rho_p, pear = float("nan"), float("nan"), float("nan")
    sub7 = df[(df["cohort"] == CONCORDANCE_COHORT) & df["responder"].notna()
              & df[PRIMARY_GENE].notna()].copy()
    high, _, _ = split_high(sub7[PRIMARY_GENE].values, "median")
    s7 = fisher_block(table_2x2(high, sub7["responder"].values == 1))
    s7.update({"cohort": CONCORDANCE_COHORT, "variant": "S7_predictio_harmonised_same_trial"})
    concordance = {
        "n_matched_patients": int(len(merged)),
        "spearman_rho": float(rho), "spearman_p": float(rho_p),
        "pearson_r": pear,
        "primary_test_on_harmonised_values": s7,
        "note": "Same trial, independent processing pipeline. Not pooled; concordance check only.",
    }
    print(f"[S7] CLDN4 own vs PredictIO: rho={rho:.3f} on n={len(merged)}")

    # ===================== S8 multivariable ================================= #
    import statsmodels.formula.api as smf
    m8 = df[(df["cohort"] == "IMvigor210") & df["responder"].notna()
            & df[PRIMARY_GENE].notna()].copy()
    m8["tmb"] = pd.to_numeric(m8["tmb_per_mb"], errors="coerce")
    m8["log_tmb"] = np.log2(m8["tmb"] + 1)
    m8["bladder"] = (m8["tissue"].astype(str).str.lower() == "bladder").astype(int)
    m8["ecog_n"] = pd.to_numeric(m8["ecog"], errors="coerce")
    m8["platinum"] = m8["received_platinum"].astype(str)
    m8["phenotype"] = m8["immune_phenotype"].astype(str)
    m8["cldn4_high"] = (m8[PRIMARY_GENE] > m8[PRIMARY_GENE].median()).astype(int)
    m8["resp"] = m8["responder"].astype(int)
    need = ["resp", "cldn4_high", PRIMARY_GENE, "log_tmb", "phenotype", "ecog_n",
            "platinum", "bladder"]
    m8v = m8.dropna(subset=need)
    m8v = m8v[m8v["phenotype"].isin(["desert", "excluded", "inflamed"])]
    multivar = {"n_analysed": int(len(m8v)),
                "n_dropped_missing_covariate": int(len(m8) - len(m8v)),
                "adjustment_set": "log2(TMB+1), immune phenotype, ECOG, received platinum, bladder vs non-bladder biopsy"}
    for label, exposure in (("dichotomous_cldn4_high", "cldn4_high"),
                            ("continuous_cldn4", PRIMARY_GENE)):
        try:
            fit = smf.logit(f"resp ~ {exposure} + log_tmb + C(phenotype) + ecog_n "
                            f"+ C(platinum) + bladder", data=m8v).fit(disp=0)
            b, se = float(fit.params[exposure]), float(fit.bse[exposure])
            multivar[label] = {
                "adjusted_or": math.exp(b),
                "ci_low": math.exp(b - 1.96 * se), "ci_high": math.exp(b + 1.96 * se),
                "p": float(fit.pvalues[exposure]),
            }
            print(f"[S8] {label}: adjusted OR={math.exp(b):.3f} p={fit.pvalues[exposure]:.4g}")
        except Exception as exc:
            multivar[label] = {"error": str(exc)}

    # ===================== Exploratory (BH within family) =================== #
    expl = []
    sub_i = df[(df["cohort"] == "IMvigor210") & df["responder"].notna()].copy()
    sub_i["responder"] = sub_i["responder"].astype(int)
    for gene in EXPLORATORY_GENES:
        s = sub_i[sub_i[gene].notna()]
        if len(s) < 20:
            continue
        high, _, _ = split_high(s[gene].values, "median")
        blk = fisher_block(table_2x2(high, s["responder"].values == 1))
        blk.update({"gene": gene, "cohort": "IMvigor210"})
        expl.append(blk)
    if expl:
        ps = np.array([e["fisher_p"] for e in expl])
        order = np.argsort(ps)
        m = len(ps)
        q = np.empty(m)
        prev = 1.0
        for rank, idx in enumerate(order[::-1]):
            i = m - rank
            prev = min(prev, ps[idx] * m / i)
            q[idx] = prev
        for e, qq in zip(expl, q):
            e["bh_q_value"] = float(qq)
        expl_meta = {"family": "claudin/junction/checkpoint genes, median split, IMvigor210",
                     "m_tests": m, "method": "Benjamini-Hochberg",
                     "threshold_q": 0.05,
                     "n_passing": int(np.sum(q < 0.05)),
                     "note": "CLDN4 is the primary hypothesis and is deliberately NOT in this family."}
    else:
        expl_meta = {}

    # Exploratory: CLDN4 vs immune context
    ctx = df[(df["cohort"] == "IMvigor210") & df[PRIMARY_GENE].notna()].copy()
    rho_cd8, p_cd8 = stats.spearmanr(ctx[PRIMARY_GENE], ctx["cd8_teff_z"])
    pheno_stats = {}
    groups = []
    for ph in ["desert", "excluded", "inflamed"]:
        v = ctx.loc[ctx["immune_phenotype"] == ph, PRIMARY_GENE].dropna().values
        if len(v):
            pheno_stats[ph] = {"n": int(len(v)), "median_cldn4": float(np.median(v))}
            groups.append(v)
    kw = stats.kruskal(*groups) if len(groups) >= 2 else None
    immune_context = {
        "spearman_cldn4_vs_cd8_teff": {"rho": float(rho_cd8), "p": float(p_cd8), "n": int(len(ctx))},
        "cldn4_by_immune_phenotype": pheno_stats,
        "kruskal_wallis": {"H": float(kw.statistic), "p": float(kw.pvalue)} if kw else None,
        "label": "EXPLORATORY - describes the immune-exclusion thesis, not a test of the B5 claim",
    }

    # ===================== write tables ===================================== #
    def write_csv(name, rows, cols=None):
        if not rows:
            return
        cols = cols or sorted({k for r in rows for k in r})
        p = os.path.join(args.out_dir, name)
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"[write] {p}")

    primary_cols = ["cohort", "cutoff_rule", "threshold", "expression_scale", "n",
                    "high_responder", "high_nonresponder", "low_responder",
                    "low_nonresponder", "orr_high_pct", "orr_low_pct",
                    "or_conditional_mle", "ci_low", "ci_high", "fisher_p",
                    "log_or_woolf", "se_log_or", "haldane_corrected"]
    write_csv("cohort_effects.csv", per_cohort, primary_cols)
    write_csv("sensitivity_cutoffs.csv", sens_cutoffs,
              ["cohort", "cutoff_rule", "n", "high_responder", "high_nonresponder",
               "low_responder", "low_nonresponder", "orr_high_pct", "orr_low_pct",
               "or_conditional_mle", "ci_low", "ci_high", "fisher_p", "I2_percent"])
    write_csv("controls.csv", controls,
              ["control", "score", "cohort", "expectation", "n", "orr_high_pct",
               "orr_low_pct", "or_conditional_mle", "ci_low", "ci_high", "fisher_p"])
    write_csv("controls_continuous.csv", controls_continuous,
              ["control", "score", "cohort", "n", "sd_of_exposure", "or_per_sd",
               "ci_low_per_sd", "ci_high_per_sd", "p", "error"])
    write_csv("continuous_logistic.csv", cont_rows,
              ["cohort", "n", "sd_of_exposure", "or_per_log2_unit",
               "ci_low_per_log2_unit", "ci_high_per_log2_unit", "or_per_sd",
               "ci_low_per_sd", "ci_high_per_sd", "p", "error"])
    write_csv("mannwhitney.csv", mw_rows,
              ["cohort", "n_responder", "n_nonresponder", "median_responder",
               "median_nonresponder", "u", "p", "rank_biserial_r", "error"])
    write_csv("exploratory_genes.csv", expl,
              ["gene", "cohort", "n", "orr_high_pct", "orr_low_pct",
               "or_conditional_mle", "ci_low", "ci_high", "fisher_p", "bh_q_value"])

    meta_rows = []
    for model, key in (("fixed_iv", "fixed_iv"), ("random_DL", "random_dl"),
                       ("random_DL_hartung_knapp", "random_dl_hartung_knapp"),
                       ("mantel_haenszel", "mantel_haenszel")):
        blk = pooled[key]
        meta_rows.append({
            "model": model, "k": pooled["k"], "or": blk.get("or"),
            "ci_low": blk.get("ci_low"), "ci_high": blk.get("ci_high"),
            "p": blk.get("p"), "Q": pooled["heterogeneity"]["Q"],
            "Q_df": pooled["heterogeneity"]["df"], "Q_p": pooled["heterogeneity"]["p"],
            "I2_percent": pooled["heterogeneity"]["I2_percent"],
            "tau2": pooled["heterogeneity"]["tau2"],
            "user_claimed_or": USER_OR,
        })
    write_csv("meta_summary.csv", meta_rows,
              ["model", "k", "or", "ci_low", "ci_high", "p", "Q", "Q_df", "Q_p",
               "I2_percent", "tau2", "user_claimed_or"])

    # ===================== figures ========================================== #
    figdir = args.out_dir
    # Fig 1: forest
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    labels = [e["cohort"] for e in effects] + ["Random effects (DL)"]
    ors = [e["or_conditional_mle"] for e in effects] + [pooled["random_dl"]["or"]]
    los = [e["ci_low"] for e in effects] + [pooled["random_dl"]["ci_low"]]
    his = [e["ci_high"] for e in effects] + [pooled["random_dl"]["ci_high"]]
    ns = [e["n"] for e in effects] + [sum(e["n"] for e in effects)]
    ypos = list(range(len(labels)))[::-1]
    for i, (y, o, lo, hi) in enumerate(zip(ypos, ors, los, his)):
        is_pooled = i == len(labels) - 1
        hi_plot = min(hi, 12)
        ax.plot([lo, hi_plot], [y, y], color="black", lw=1.6, zorder=2)
        ax.plot([o], [y], marker="D" if is_pooled else "o",
                ms=11 if is_pooled else 8,
                color="#c0392b" if is_pooled else "#2c6fbb", zorder=3)
        if hi > 12:
            ax.annotate("", xy=(12.4, y), xytext=(11.4, y),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.4))
    ax.axvline(1.0, color="grey", ls="--", lw=1)
    ax.axvline(USER_OR, color="#e67e22", ls=":", lw=1.8,
               label=f"claimed OR = {USER_OR}")
    ax.set_xscale("log")
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{l}  (n={n})" for l, n in zip(labels, ns)])
    ax.set_xlabel("Odds ratio of objective response, CLDN4-high vs CLDN4-low (log scale)")
    ax.set_title("B5_BLCA: CLDN4-high and ICI objective response in urothelial carcinoma\n"
                 f"median split; random-effects OR = {pooled['random_dl']['or']:.2f} "
                 f"({pooled['random_dl']['ci_low']:.2f}\u2013{pooled['random_dl']['ci_high']:.2f}), "
                 f"I\u00b2 = {pooled['heterogeneity']['I2_percent']:.0f}%", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "fig_forest_meta.png"), dpi=160)
    plt.close(fig)

    # Fig 2: ORR bars
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    x = np.arange(len(effects))
    ax.bar(x - 0.19, [e["orr_low_pct"] for e in effects], 0.38,
           label="CLDN4-low", color="#7fb3d5", edgecolor="black", lw=0.6)
    ax.bar(x + 0.19, [e["orr_high_pct"] for e in effects], 0.38,
           label="CLDN4-high", color="#e59866", edgecolor="black", lw=0.6)
    for i, e in enumerate(effects):
        ax.text(i - 0.19, e["orr_low_pct"] + 1.2,
                f"{e['low_responder']}/{e['low_responder']+e['low_nonresponder']}",
                ha="center", fontsize=8)
        ax.text(i + 0.19, e["orr_high_pct"] + 1.2,
                f"{e['high_responder']}/{e['high_responder']+e['high_nonresponder']}",
                ha="center", fontsize=8)
        ax.text(i, max(e["orr_low_pct"], e["orr_high_pct"]) + 6.0,
                f"p = {e['fisher_p']:.2f}", ha="center", fontsize=8, style="italic")
    ax.set_xticks(x)
    ax.set_xticklabels([e["cohort"] for e in effects])
    ax.set_ylabel("Objective response rate (%)")
    ax.set_ylim(0, max(max(e["orr_high_pct"], e["orr_low_pct"]) for e in effects) * 1.34)
    ax.set_title("Objective response rate by CLDN4 median split", fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "fig_orr_by_cldn4.png"), dpi=160)
    plt.close(fig)

    # Fig 3: CLDN4 distribution by response
    fig, axes = plt.subplots(1, len(POOLED_COHORTS), figsize=(11, 4), sharey=False)
    for axx, cohort, mw in zip(np.atleast_1d(axes), POOLED_COHORTS, mw_rows):
        sub = analysis_set(df, cohort)
        data = [sub.loc[sub["responder"] == 0, PRIMARY_GENE].values,
                sub.loc[sub["responder"] == 1, PRIMARY_GENE].values]
        bp = axx.boxplot(data, tick_labels=["non-resp", "resp"], widths=0.55,
                         patch_artist=True, showfliers=False)
        for patch, col in zip(bp["boxes"], ["#aeb6bf", "#82c99a"]):
            patch.set_facecolor(col)
        for j, d in enumerate(data):
            jitter = np.random.normal(0, 0.055, len(d))
            axx.plot(np.full(len(d), j + 1) + jitter, d, "k.", ms=3.2, alpha=0.55)
        axx.set_title(f"{cohort}\nMann-Whitney p = {mw.get('p', float('nan')):.2f}", fontsize=10)
        axx.set_ylabel(f"CLDN4  {sub['expression_scale'].iloc[0]}")
    fig.suptitle("CLDN4 expression in ICI responders vs non-responders (urothelial carcinoma)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "fig_cldn4_by_response.png"), dpi=160)
    plt.close(fig)

    # Fig 4: cutoff robustness (pooled)
    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    pr = [r for r in sens_cutoffs if r["cohort"] == "POOLED_random_effects"]
    y = np.arange(len(pr))[::-1]
    for yy, r in zip(y, pr):
        ax.plot([r["ci_low"], min(r["ci_high"], 12)], [yy, yy], color="black", lw=1.5)
        ax.plot([r["or_conditional_mle"]], [yy], "D", color="#c0392b", ms=9)
    ax.axvline(1.0, color="grey", ls="--", lw=1)
    ax.axvline(USER_OR, color="#e67e22", ls=":", lw=1.8, label=f"claimed OR = {USER_OR}")
    ax.set_yticks(y)
    ax.set_yticklabels([r["cutoff_rule"] for r in pr])
    ax.set_xscale("log")
    ax.set_xlabel("Pooled random-effects OR of response, CLDN4-high vs low")
    ax.set_title("S1: sensitivity of the pooled estimate to the CLDN4 cutoff", fontsize=11)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "fig_cutoff_robustness.png"), dpi=160)
    plt.close(fig)

    # Fig 5: controls
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    rows5 = controls + [dict(e, control="CLDN4 (primary)") for e in effects]
    y = np.arange(len(rows5))[::-1]
    colmap = {"positive_control_cd8_teff": "#27ae60",
              "negative_control_housekeeping": "#7f8c8d",
              "CLDN4 (primary)": "#2c6fbb"}
    for yy, r in zip(y, rows5):
        ax.plot([r["ci_low"], min(r["ci_high"], 30)], [yy, yy], color="black", lw=1.4)
        ax.plot([r["or_conditional_mle"]], [yy], "o", ms=8,
                color=colmap.get(r.get("control"), "#2c6fbb"))
    ax.axvline(1.0, color="grey", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.get('control','')} | {r['cohort']}" for r in rows5], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("OR of objective response, score-high vs score-low")
    ax.set_title("Pipeline calibration: controls vs the CLDN4 primary", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "fig_controls.png"), dpi=160)
    plt.close(fig)
    print("[write] 5 figures")

    # ===================== summary.json ===================================== #
    n_total = sum(e["n"] for e in effects)
    n_resp = sum(e["high_responder"] + e["low_responder"] for e in effects)
    ndir = sum(1 for e in effects if e["or_conditional_mle"] < 1)
    sig = [e["cohort"] for e in effects if e["fisher_p"] < 0.05]

    re = pooled["random_dl"]
    verdict["statement"] = (
        f"Across {pooled['k']} public urothelial ICI cohorts (n={n_total} response-evaluable "
        f"patients, {n_resp} objective responders), the random-effects pooled odds ratio of "
        f"objective response for CLDN4-high versus CLDN4-low tumours is "
        f"{re['or']:.2f} (95% CI {re['ci_low']:.2f}-{re['ci_high']:.2f}, p={re['p']:.2f}). "
        f"The user-claimed value of {USER_OR} is "
        f"{'inside' if verdict['claimed_value_inside_observed_95ci'] else 'outside'} "
        f"this confidence interval, and the observed estimate does not match {USER_OR} at "
        f"2 decimal places (label: {verdict['label']}). "
        f"{ndir} of {pooled['k']} cohorts point in the claimed direction (OR<1); "
        f"{len(sig)} cohort-level test(s) reached p<0.05"
        + (f" ({', '.join(sig)})" if sig else "") + "."
    )

    if verdict["claimed_value_inside_observed_95ci"]:
        conclusion = "suggestive"
    elif claim_test["rejected_at_0.05"] and power["power_two_sided_alpha_0.05"] >= 0.8:
        conclusion = "inconsistent"
    elif claim_test["rejected_at_0.05"]:
        conclusion = "inconsistent"
    else:
        conclusion = "not testable"

    summary = {
        "task": "B5_BLCA",
        "conclusion": conclusion,
        "conclusion_vocabulary": "supported | suggestive | inconsistent | not testable",
        "question": "In public urothelial/bladder ICI cohorts, do CLDN4-high tumours have lower ORR?",
        "claim_under_test": "B5: 11-cohort ICI meta-analysis, CLDN4-high OR=0.42 (urothelial arm only tested here)",
        "user_claimed_or": USER_OR,
        "user_claimed_stat": "odds ratio of objective response for CLDN4-high vs CLDN4-low",
        "claim_provenance": ("No published CLDN4 x ICI meta-analysis could be located "
                            "(PubMed '\"claudin 4\" AND (\"checkpoint inhibitor\" OR \"anti-PD-1\")' "
                            "returns 0 records). OR=0.42 is treated as an unsourced assertion under "
                            "test, not a published result being reproduced. We do not name 11 cohorts "
                            "we cannot enumerate; only 3 eligible public urothelial ICI cohorts with "
                            "both RNA and RECIST were found."),
        "prespecified": True,
        "analysis_plan": "results/w200/B5_BLCA/analysis_plan.md (committed before results)",
        "seed": SEED,
        "primary_endpoint": "ORR: responder = CR or PR; non-responder = SD or PD; non-evaluable excluded",
        "primary_exposure": "CLDN4 > cohort median (response-blind), ties to low group",
        "cohorts_pooled": POOLED_COHORTS,
        "n_response_evaluable_total": n_total,
        "n_objective_responders_total": n_resp,
        "per_cohort": per_cohort,
        "meta": pooled,
        "honest_verdict": verdict,
        "controls": {
            "rows": controls,
            "continuous": controls_continuous,
            "interpretation": ("The positive control must show OR>1 for the pipeline to be "
                               "credible; the negative control must be near 1. Read the CLDN4 "
                               "result only in light of these."),
        },
        "sensitivity": {
            "S1_cutoffs": sens_cutoffs,
            "S2_continuous_logistic": cont_rows,
            "S3_mannwhitney": mw_rows,
            "S4_itt_ne_as_nonresponder": {"per_cohort": itt_effects, "meta": itt_pooled},
            "S5_bladder_biopsy_only": s5,
            "S6_deseq_normalisation": s6,
            "S7_pipeline_concordance": concordance,
            "S8_multivariable": multivar,
        },
        "exploratory": {
            "gene_family": expl,
            "gene_family_meta": expl_meta,
            "immune_context": immune_context,
        },
        "methods": {
            "expression_units": {"IMvigor210": "log2(TPM+1) from counts/gene-length",
                                 "BACI": "log2(TPM+1) from salmon gene TPM",
                                 "Snyder": "log2(TPM+0.001) as deposited by PredictIO"},
            "per_cohort_test": "two-sided Fisher exact; conditional-MLE OR with exact 95% CI",
            "pooling": "log OR (Woolf, Haldane-Anscombe 0.5 only if a zero cell) pooled by "
                       "DerSimonian-Laird random effects (primary), with Hartung-Knapp limits, "
                       "fixed-effect inverse variance and Mantel-Haenszel also reported",
            "multiplicity": "primary = 1 prespecified pooled test; exploratory gene family "
                            "BH-adjusted with m reported; sensitivity analyses are robustness "
                            "checks and are not alpha-spending confirmations",
            "did_we_tune_to_0.42": DID_WE_TUNE,
            "double_counting_control": "IMvigor210 contributes one estimate; the PredictIO "
                                       "harmonisation of the same trial is used only for S7",
        },
        "limitations": [
            "Three cohorts, not eleven; the urothelial arm alone cannot confirm or refute a pooled 11-cohort number.",
            "Snyder contributes only 21 evaluable patients, so its interval is very wide.",
            "BACI is a real-world cohort with mixed ICI agents and investigator-assessed response.",
            "IMvigor210 was obtained from a third-party mirror because the official distribution now 404s; accepted only after it reproduced published RECIST and immune-phenotype distributions exactly.",
            "Expression scales differ across cohorts; the median split is within-cohort, which makes the contrast scale-free but does not harmonise the underlying platforms.",
            "Observational association only; no causal or predictive-utility claim is implied.",
        ],
    }
    spath = os.path.join(args.out_dir, "summary.json")
    with open(spath, "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    print(f"[write] {spath}")
    print("\n=== VERDICT ===")
    print(verdict["statement"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
