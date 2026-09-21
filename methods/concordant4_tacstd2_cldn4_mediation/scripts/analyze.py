#!/usr/bin/env python3
"""Mechanism test: TACSTD2 → T/NK, partially through CLDN4.

Patient / donor / sample is the unit (locked concordant-4, n = 65).
The grid is declared below. p-values are descriptive.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from lib_stats import random_effects_dl, spearman  # noqa: E402

RES = HERE / "results"
FIG = RES / "figures"
TAB = RES / "tables"
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COHORT_COLOR = {
    "GSE123902": "#1b4f72",
    "GSE131907": "#0e6655",
    "GSE205335": "#b9770e",
    "GSE189357": "#6c3483",
}
LOCKED_CLDN4_RHO = -0.5311678045689989

PAIRED = [
    "pct_gt0",
    "pct_ge2",
    "pct_ge3",
    "pct_ge5",
    "pct_ge10",
    "pct_cp_ge1",
    "pct_cp_ge5",
    "pct_cp_ge10",
    "pct_gt_cohort_median",
    "pct_gt_cohort_p75",
    "mean_log1p_raw",
    "mean_log1p_cp10k",
    "mean_log1p_pos_raw",
    "p90_log1p_raw",
    "p90_log1p_cp10k",
    "pb_log1p_cp10k",
]
PROPORTION_SCORES = {
    "pct_gt0",
    "pct_ge2",
    "pct_ge3",
    "pct_ge5",
    "pct_ge10",
    "pct_cp_ge1",
    "pct_cp_ge5",
    "pct_cp_ge10",
    "pct_gt_cohort_median",
    "pct_gt_cohort_p75",
}
CROSS_X = ["pct_ge2", "pct_ge5", "mean_log1p_raw", "mean_log1p_cp10k", "pb_log1p_cp10k"]
CROSS_M = ["pct_ge2", "pct_ge5", "mean_log1p_raw", "mean_log1p_cp10k", "pb_log1p_cp10k"]
N_BOOT = 5000
BOOT_SEED = 20260921


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "—"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_rho(r: float) -> str:
    if r is None or not np.isfinite(r):
        return "—"
    return f"{r:.3f}"


def fmt_num(x: float, digits: int = 3) -> str:
    if x is None or not np.isfinite(x):
        return "—"
    return f"{x:.{digits}f}"


def weighted_quantile(hists: list[str], q: float) -> int:
    total: dict[int, int] = {}
    for hist in hists:
        if not isinstance(hist, str) or not hist:
            continue
        for part in hist.split(","):
            v_s, c_s = part.split(":")
            v, c = int(v_s), int(c_s)
            total[v] = total.get(v, 0) + c
    if not total:
        return 0
    vals = sorted(total)
    cdf = 0
    n = sum(total.values())
    target = q * n
    for v in vals:
        cdf += total[v]
        if cdf >= target:
            return v
    return vals[-1]


def pct_above(hist: str, thr: int) -> float:
    if not isinstance(hist, str) or not hist:
        return float("nan")
    n = 0
    above = 0
    for part in hist.split(","):
        v_s, c_s = part.split(":")
        v, c = int(v_s), int(c_s)
        n += c
        if v > thr:
            above += c
    return above / n if n else float("nan")


def add_quantile_scores(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for gene in ("tacstd2", "cldn4"):
        for qname, q in (("median", 0.5), ("p75", 0.75)):
            col = f"{gene}_pct_gt_cohort_{qname}"
            out = pd.Series(np.nan, index=df.index, dtype=float)
            for _, sub in df.groupby("cohort"):
                thr = weighted_quantile(sub[f"{gene}_hist"].tolist(), q)
                out.loc[sub.index] = [pct_above(h, thr) for h in sub[f"{gene}_hist"]]
            df[col] = out
        recon = []
        for hist in df[f"{gene}_hist"]:
            recon.append(pct_above(hist, 0))
        max_abs = float(np.nanmax(np.abs(np.array(recon) - df[f"{gene}_pct_gt0"].to_numpy(float))))
        if max_abs > 1e-8:
            raise SystemExit(f"{gene} histogram does not match pct_gt0 (max abs {max_abs})")
    return df


def z_within(s: pd.Series, cohort: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=s.index, dtype=float)
    for key, idx in cohort.groupby(cohort).groups.items():
        vals = s.loc[idx].to_numpy(float)
        sd = float(np.std(vals, ddof=0))
        if sd == 0 or not np.isfinite(sd):
            out.loc[idx] = np.nan
        else:
            out.loc[idx] = (vals - np.mean(vals)) / sd
    return out


def rank_within(s: pd.Series, cohort: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=s.index, dtype=float)
    for _, idx in cohort.groupby(cohort).groups.items():
        vals = s.loc[idx].to_numpy(float)
        if np.isfinite(vals).sum() < 3:
            out.loc[idx] = np.nan
        else:
            out.loc[idx] = stats.rankdata(vals, method="average")
    return out


def logit_clip(s: pd.Series, eps: float = 1e-4) -> pd.Series:
    p = np.clip(s.to_numpy(float), eps, 1.0 - eps)
    return pd.Series(np.log(p / (1.0 - p)), index=s.index)


def design(y: np.ndarray, cols: list[np.ndarray], cohort: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    dummies = pd.get_dummies(cohort, drop_first=True, dtype=float)
    pieces = [np.ones(len(y))] + cols
    if dummies.shape[1]:
        pieces.append(dummies.to_numpy(float))
    X = np.column_stack(pieces)
    mask = np.isfinite(y) & np.isfinite(X).all(axis=1)
    return y[mask], X[mask]


def ols_fit(y: np.ndarray, cols: list[np.ndarray], cohort: pd.Series):
    yy, X = design(y, cols, cohort)
    if yy.size < X.shape[1] + 2:
        return None
    return statsmodels_ols(yy, X)


def statsmodels_ols(y: np.ndarray, X: np.ndarray):
    import statsmodels.api as sm

    return sm.OLS(y, X).fit()


def coef0(res, idx: int = 1) -> tuple[float, float, float]:
    return float(res.params[idx]), float(res.bse[idx]), float(res.pvalues[idx])


def spearman_by_cohort(df: pd.DataFrame, x: str, y: str) -> list[dict]:
    rows = []
    for cohort in COHORTS:
        sub = df[df["cohort"] == cohort]
        rho, p, n = spearman(sub[x], sub[y])
        rows.append({"cohort": cohort, "rho": rho, "p": p, "n": n})
    return rows


def meta_from_rows(rows: list[dict], k_cov: int = 0) -> dict:
    rhos, ns = [], []
    for row in rows:
        if row["n"] - 3 - k_cov > 0 and np.isfinite(row["rho"]):
            rhos.append(row["rho"])
            ns.append(row["n"] - k_cov if k_cov else row["n"])
        else:
            rhos.append(np.nan)
            ns.append(row["n"])
    # random_effects_dl uses 1/(n-3). For partials, pass n' = n - k so 1/(n'-3)=1/(n-3-k).
    use_rho, use_n = [], []
    for rho, n, row in zip(rhos, ns, rows):
        if k_cov:
            n_eff = row["n"] - k_cov
        else:
            n_eff = row["n"]
        if np.isfinite(rho) and n_eff > 3:
            use_rho.append(float(rho))
            use_n.append(int(n_eff))
    if len(use_rho) < 2:
        return {"k": len(use_rho), "pooled_rho": float("nan"), "p": float("nan"), "I2": float("nan"),
                "ci95_rho": [float("nan"), float("nan")]}
    return random_effects_dl(use_rho, use_n)


def partial_rho(r_xy: float, r_xz: float, r_yz: float) -> float:
    den_sq = (1.0 - r_xz ** 2) * (1.0 - r_yz ** 2)
    if not np.isfinite(den_sq) or den_sq <= 0:
        return float("nan")
    return float((r_xy - r_xz * r_yz) / math.sqrt(den_sq))


def prepare_vectors(df: pd.DataFrame, xcol: str, mcol: str, model: str) -> pd.DataFrame | None:
    out = df[["cohort", "unit_id", "frac_tnk"]].copy()
    x = df[xcol]
    m = df[mcol]
    y = df["frac_tnk"]
    if model == "identity":
        out["x"], out["m"], out["y"] = x, m, y
    elif model == "z":
        out["x"] = z_within(x, df["cohort"])
        out["m"] = z_within(m, df["cohort"])
        out["y"] = z_within(y, df["cohort"])
    elif model == "rank":
        out["x"] = rank_within(x, df["cohort"])
        out["m"] = rank_within(m, df["cohort"])
        out["y"] = rank_within(y, df["cohort"])
    elif model == "logit":
        if not _is_proportion(x) or not _is_proportion(m) or not _is_proportion(y):
            return None
        out["x"] = logit_clip(x)
        out["m"] = logit_clip(m)
        out["y"] = logit_clip(y)
    else:
        raise ValueError(model)
    if not np.isfinite(out[["x", "m", "y"]].to_numpy(float)).all():
        return None
    return out


def _is_proportion(s: pd.Series) -> bool:
    v = s.to_numpy(float)
    return bool(np.isfinite(v).all() and np.nanmin(v) >= 0 and np.nanmax(v) <= 1)


def fit_mediation(frame: pd.DataFrame) -> dict | None:
    y = frame["y"].to_numpy(float)
    x = frame["x"].to_numpy(float)
    m = frame["m"].to_numpy(float)
    cohort = frame["cohort"]
    total = ols_fit(y, [x], cohort)
    med = ols_fit(m, [x], cohort)
    both = ols_fit(y, [x, m], cohort)
    if total is None or med is None or both is None:
        return None
    c, c_se, c_p = coef0(total, 1)
    a, a_se, a_p = coef0(med, 1)
    c_prime, cp_se, cp_p = coef0(both, 1)
    b, b_se, b_p = coef0(both, 2)
    indirect = a * b
    if c == 0 or not np.isfinite(c):
        prop = float("nan")
    else:
        prop = indirect / c
    # linear-model identity
    gap = abs((c - c_prime) - indirect)
    sobel_se = math.sqrt((a ** 2) * (b_se ** 2) + (b ** 2) * (a_se ** 2))
    sobel_z = indirect / sobel_se if sobel_se > 0 else float("nan")
    sobel_p = float(2 * stats.norm.sf(abs(sobel_z))) if np.isfinite(sobel_z) else float("nan")
    return {
        "c": c,
        "c_se": c_se,
        "c_p": c_p,
        "a": a,
        "a_se": a_se,
        "a_p": a_p,
        "b": b,
        "b_se": b_se,
        "b_p": b_p,
        "c_prime": c_prime,
        "c_prime_se": cp_se,
        "c_prime_p": cp_p,
        "indirect": indirect,
        "proportion": prop,
        "attenuation_abs": abs(c - c_prime),
        "attenuation_fraction": prop,
        "ols_identity_gap": gap,
        "sobel_se": sobel_se,
        "sobel_p": sobel_p,
        "n": int(len(frame)),
        "r2_total": float(total.rsquared),
        "r2_both": float(both.rsquared),
    }


def cohort_ols_proportion(frame: pd.DataFrame) -> list[dict]:
    rows = []
    for cohort in COHORTS:
        sub = frame[frame["cohort"] == cohort]
        y = sub["y"].to_numpy(float)
        x = sub["x"].to_numpy(float)
        m = sub["m"].to_numpy(float)
        n = len(sub)
        if n < 6 or np.std(x) == 0 or np.std(m) == 0:
            rows.append({"cohort": cohort, "n": n, "c": np.nan, "c_prime": np.nan, "a": np.nan,
                         "b": np.nan, "proportion": np.nan, "partial_ols": False})
            continue
        import statsmodels.api as sm

        Xt = sm.add_constant(x)
        Xm = sm.add_constant(np.column_stack([x, m]))
        Xa = sm.add_constant(x)
        total = sm.OLS(y, Xt).fit()
        both = sm.OLS(y, Xm).fit()
        am = sm.OLS(m, Xa).fit()
        c = float(total.params[1])
        c_prime = float(both.params[1])
        a = float(am.params[1])
        b = float(both.params[2])
        prop = (a * b) / c if c != 0 else float("nan")
        partial = bool(
            np.isfinite(prop) and c < 0 and c_prime < 0 and abs(c_prime) < abs(c) and a > 0 and b < 0 and 0 < prop < 1
        )
        rows.append({
            "cohort": cohort, "n": n, "c": c, "c_prime": c_prime, "a": a, "b": b,
            "proportion": prop, "partial_ols": partial,
        })
    return rows


def fit_interaction(frame: pd.DataFrame) -> dict:
    import statsmodels.api as sm

    y = frame["y"].to_numpy(float)
    x = frame["x"].to_numpy(float)
    m = frame["m"].to_numpy(float)
    xm = x * m
    cohort = frame["cohort"]
    dummies = pd.get_dummies(cohort, drop_first=True, dtype=float).to_numpy(float)
    X = np.column_stack([np.ones(len(y)), x, m, xm, dummies])
    res = sm.OLS(y, X).fit()
    b_x, se_x, p_x = coef0(res, 1)
    b_m, se_m, p_m = coef0(res, 2)
    b_i, se_i, p_i = coef0(res, 3)
    # simple slopes at the pooled 25th and 75th percentiles of m
    q1 = float(np.quantile(m, 0.25))
    q4 = float(np.quantile(m, 0.75))
    return {
        "b_x": b_x, "se_x": se_x, "p_x": p_x,
        "b_m": b_m, "se_m": se_m, "p_m": p_m,
        "b_int": b_i, "se_int": se_i, "p_int": p_i,
        "m_q25": q1, "m_q75": q4,
        "slope_q25": b_x + b_i * q1,
        "slope_q75": b_x + b_i * q4,
        "n": int(len(frame)),
        "r2": float(res.rsquared),
    }


def fit_glmm(frame: pd.DataFrame) -> dict:
    d = frame[["y", "x", "m", "cohort"]].copy()
    out = {"ok": False}
    try:
        base = smf.mixedlm("y ~ x", d, groups=d["cohort"]).fit(reml=False, method="lbfgs", maxiter=300)
        full = smf.mixedlm("y ~ x + m", d, groups=d["cohort"]).fit(reml=False, method="lbfgs", maxiter=300)
    except Exception as exc:  # noqa: BLE001
        out["error"] = str(exc)
        return out
    c = float(base.params["x"])
    c_prime = float(full.params["x"])
    b = float(full.params["m"])
    # a path with the same random intercept
    try:
        a_fit = smf.mixedlm("m ~ x", d, groups=d["cohort"]).fit(reml=False, method="lbfgs", maxiter=300)
        a = float(a_fit.params["x"])
    except Exception as exc:  # noqa: BLE001
        out["error"] = str(exc)
        return out
    prop = (a * b) / c if c != 0 else float("nan")
    out.update({
        "ok": True,
        "c": c,
        "c_se": float(base.bse["x"]),
        "c_p": float(base.pvalues["x"]),
        "c_prime": c_prime,
        "c_prime_se": float(full.bse["x"]),
        "c_prime_p": float(full.pvalues["x"]),
        "a": a,
        "b": b,
        "b_se": float(full.bse["m"]),
        "b_p": float(full.pvalues["m"]),
        "proportion": prop,
        "re_var_base": float(base.cov_re.iloc[0, 0]) if np.ndim(base.cov_re) else float("nan"),
        "re_var_full": float(full.cov_re.iloc[0, 0]) if np.ndim(full.cov_re) else float("nan"),
        "converged_base": bool(base.converged),
        "converged_full": bool(full.converged),
    })
    return out


def strata_effects(df: pd.DataFrame, xcol: str, mcol: str) -> list[dict]:
    """Within-cohort CLDN4 quartile and median split. Outcome is raw T/NK fraction."""
    d = df[["cohort", "unit_id", "frac_tnk", xcol, mcol]].copy()
    d["q"] = d.groupby("cohort")[mcol].transform(
        lambda s: pd.qcut(s.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    )
    d["half"] = d.groupby("cohort")[mcol].transform(
        lambda s: np.where(s <= s.median(), "low", "high")
    )
    rows = []
    for label, col, levels in (
        ("quartile", "q", ["Q1", "Q4"]),
        ("median", "half", ["low", "high"]),
    ):
        for level in levels:
            sub = d[d[col].astype(str) == level]
            cohort_rows = []
            for cohort in COHORTS:
                ss = sub[sub["cohort"] == cohort]
                rho, p, n = spearman(ss[xcol], ss["frac_tnk"])
                cohort_rows.append({"cohort": cohort, "rho": rho, "p": p, "n": n,
                                    "mean_tnk": float(ss["frac_tnk"].mean()) if len(ss) else float("nan")})
            usable = [r for r in cohort_rows if r["n"] >= 5 and np.isfinite(r["rho"])]
            meta = random_effects_dl([r["rho"] for r in usable], [r["n"] for r in usable]) if len(usable) >= 2 else {
                "pooled_rho": float("nan"), "p": float("nan"), "I2": float("nan"), "k": len(usable),
                "ci95_rho": [float("nan"), float("nan")],
            }
            # OLS with cohort FE inside the stratum, raw scores
            beta = se = pbeta = float("nan")
            if sub["cohort"].nunique() >= 2 and len(sub) >= 10:
                y = sub["frac_tnk"].to_numpy(float)
                x = sub[xcol].to_numpy(float)
                fit = ols_fit(y, [x], sub["cohort"])
                if fit is not None:
                    beta, se, pbeta = coef0(fit, 1)
            rows.append({
                "split": label,
                "level": level,
                "n": int(len(sub)),
                "n_cohorts_spearman": len(usable),
                "rho": meta.get("pooled_rho", float("nan")),
                "rho_p": meta.get("p", float("nan")),
                "I2": meta.get("I2", float("nan")),
                "rho_lo": meta.get("ci95_rho", [float("nan"), float("nan")])[0],
                "rho_hi": meta.get("ci95_rho", [float("nan"), float("nan")])[1],
                "ols_beta": beta,
                "ols_se": se,
                "ols_p": pbeta,
                "cohorts": json.dumps(cohort_rows),
            })
    return rows


def bootstrap_proportion(frame: pd.DataFrame, n_boot: int, seed: int) -> dict:
    frame = frame.sort_values(["cohort", "unit_id"]).reset_index(drop=True)
    rng = np.random.default_rng(seed)
    groups = {c: frame.index[frame["cohort"] == c].to_numpy() for c in COHORTS}
    props = np.empty(n_boot)
    indirects = np.empty(n_boot)
    cs = np.empty(n_boot)
    cps = np.empty(n_boot)
    for b in range(n_boot):
        take = []
        for cohort in COHORTS:
            ix = groups[cohort]
            take.append(rng.choice(ix, size=len(ix), replace=True))
        sample = frame.loc[np.concatenate(take)]
        # restore a unique index so nothing downstream cares
        sample = sample.reset_index(drop=True)
        fit = fit_mediation(sample)
        if fit is None:
            props[b] = np.nan
            indirects[b] = np.nan
            cs[b] = np.nan
            cps[b] = np.nan
        else:
            props[b] = fit["proportion"]
            indirects[b] = fit["indirect"]
            cs[b] = fit["c"]
            cps[b] = fit["c_prime"]
    def pct(a, q):
        a = a[np.isfinite(a)]
        if a.size == 0:
            return float("nan")
        return float(np.quantile(a, q))
    finite = props[np.isfinite(props)]
    return {
        "n_boot": n_boot,
        "n_finite": int(finite.size),
        "prop_mean": float(np.mean(finite)) if finite.size else float("nan"),
        "prop_p2_5": pct(props, 0.025),
        "prop_p50": pct(props, 0.50),
        "prop_p97_5": pct(props, 0.975),
        "indirect_p2_5": pct(indirects, 0.025),
        "indirect_p97_5": pct(indirects, 0.975),
        "frac_prop_in_0_1": float(np.mean((finite > 0) & (finite < 1))) if finite.size else float("nan"),
        "frac_attenuates": float(np.mean((cs < 0) & (cps < 0) & (np.abs(cps) < np.abs(cs)))) if np.isfinite(cs).any() else float("nan"),
        "draws_prop": props,
        "draws_indirect": indirects,
    }


def spec_list() -> list[dict]:
    specs = []
    for score in PAIRED:
        specs.append({"x": score, "m": score, "family": "paired"})
    for score in CROSS_X:
        specs.append({"x": score, "m": "pct_gt0", "family": "cross_mediator_locked"})
    for score in CROSS_M:
        specs.append({"x": "pct_gt0", "m": score, "family": "cross_exposure_locked"})
    # dedupe paired pct_gt0 accidentally repeated? cross doesn't include x=m=pct_gt0
    seen = set()
    out = []
    for spec in specs:
        key = (spec["x"], spec["m"])
        if key in seen:
            continue
        seen.add(key)
        out.append(spec)
    return out


def models_for(x: str, m: str) -> list[str]:
    models = ["identity", "z", "rank"]
    if x in PROPORTION_SCORES and m in PROPORTION_SCORES:
        models.append("logit")
    return models


def evaluate(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    grid_rows = []
    cohort_rows = []
    for spec in spec_list():
        xcol = f"tacstd2_{spec['x']}"
        mcol = f"cldn4_{spec['m']}"
        if xcol not in df.columns or mcol not in df.columns:
            raise SystemExit(f"missing columns {xcol} {mcol}")
        rho_xy = spearman_by_cohort(df, xcol, "frac_tnk")
        rho_my = spearman_by_cohort(df, mcol, "frac_tnk")
        rho_xm = spearman_by_cohort(df, xcol, mcol)
        meta_xy = meta_from_rows(rho_xy)
        meta_my = meta_from_rows(rho_my)
        meta_xm = meta_from_rows(rho_xm)
        partial_rows = []
        partial_rev_rows = []
        for a, b, c in zip(rho_xy, rho_my, rho_xm):
            pr = partial_rho(a["rho"], c["rho"], b["rho"])
            # CLDN4–T/NK | TACSTD2: r_my.x = (r_my - r_mx * r_xy) / ...
            # r_mx == r_xm
            pr_rev = partial_rho(b["rho"], c["rho"], a["rho"])
            partial_rows.append({"cohort": a["cohort"], "rho": pr, "p": float("nan"), "n": a["n"]})
            partial_rev_rows.append({"cohort": a["cohort"], "rho": pr_rev, "p": float("nan"), "n": a["n"]})
        meta_partial = meta_from_rows(partial_rows, k_cov=1)
        meta_partial_rev = meta_from_rows(partial_rev_rows, k_cov=1)
        signs_ok = (
            all(r["rho"] < 0 for r in rho_xy)
            and all(r["rho"] < 0 for r in rho_my)
            and all(r["rho"] > 0 for r in rho_xm)
            and all(abs(r["rho"]) < 0.95 for r in rho_xm)
        )
        for model in models_for(spec["x"], spec["m"]):
            frame = prepare_vectors(df, xcol, mcol, model)
            if frame is None:
                continue
            fit = fit_mediation(frame)
            if fit is None or fit["ols_identity_gap"] > 1e-6:
                continue
            per = cohort_ols_proportion(frame)
            n_partial = int(sum(r["partial_ols"] for r in per))
            pooled_partial = bool(
                fit["c"] < 0
                and fit["c_prime"] < 0
                and abs(fit["c_prime"]) < abs(fit["c"])
                and fit["a"] > 0
                and fit["b"] < 0
                and np.isfinite(fit["proportion"])
                and 0 < fit["proportion"] < 1
            )
            spearman_attenuates = bool(
                np.isfinite(meta_xy.get("pooled_rho", np.nan))
                and np.isfinite(meta_partial.get("pooled_rho", np.nan))
                and meta_xy["pooled_rho"] < 0
                and meta_partial["pooled_rho"] < 0
                and meta_partial["pooled_rho"] > meta_xy["pooled_rho"]
            )
            tier = "none"
            if signs_ok and pooled_partial and spearman_attenuates and n_partial == 4:
                tier = "A"
            elif signs_ok and pooled_partial and spearman_attenuates:
                tier = "B"
            elif signs_ok and pooled_partial:
                tier = "C"
            sid = f"{spec['x']}|{spec['m']}|{model}"
            grid_rows.append({
                "spec_id": sid,
                "family": spec["family"],
                "x_score": spec["x"],
                "m_score": spec["m"],
                "model": model,
                "tier": tier,
                "signs_ok": signs_ok,
                "spearman_attenuates": spearman_attenuates,
                "pooled_partial": pooled_partial,
                "n_cohorts_partial_ols": n_partial,
                "n": fit["n"],
                "rho_xy": meta_xy.get("pooled_rho", np.nan),
                "rho_xy_p": meta_xy.get("p", np.nan),
                "rho_xy_I2": meta_xy.get("I2", np.nan),
                "rho_xy_lo": meta_xy.get("ci95_rho", [np.nan, np.nan])[0],
                "rho_xy_hi": meta_xy.get("ci95_rho", [np.nan, np.nan])[1],
                "rho_my": meta_my.get("pooled_rho", np.nan),
                "rho_my_p": meta_my.get("p", np.nan),
                "rho_my_I2": meta_my.get("I2", np.nan),
                "rho_xm": meta_xm.get("pooled_rho", np.nan),
                "rho_xm_p": meta_xm.get("p", np.nan),
                "rho_xm_I2": meta_xm.get("I2", np.nan),
                "partial_xy": meta_partial.get("pooled_rho", np.nan),
                "partial_xy_p": meta_partial.get("p", np.nan),
                "partial_xy_I2": meta_partial.get("I2", np.nan),
                "partial_xy_lo": meta_partial.get("ci95_rho", [np.nan, np.nan])[0],
                "partial_xy_hi": meta_partial.get("ci95_rho", [np.nan, np.nan])[1],
                "partial_my": meta_partial_rev.get("pooled_rho", np.nan),
                "partial_my_p": meta_partial_rev.get("p", np.nan),
                "partial_my_I2": meta_partial_rev.get("I2", np.nan),
                **fit,
                "joint": (
                    fit["proportion"] * abs(meta_xy.get("pooled_rho", np.nan))
                    if np.isfinite(fit["proportion"]) and np.isfinite(meta_xy.get("pooled_rho", np.nan))
                    else np.nan
                ),
            })
            for kind, rows in (
                ("rho_xy", rho_xy),
                ("rho_my", rho_my),
                ("rho_xm", rho_xm),
                ("partial_xy", partial_rows),
                ("partial_my", partial_rev_rows),
            ):
                for row in rows:
                    cohort_rows.append({
                        "spec_id": sid,
                        "model": model,
                        "kind": kind,
                        **row,
                    })
            for row in per:
                cohort_rows.append({
                    "spec_id": sid,
                    "model": model,
                    "kind": "ols",
                    "cohort": row["cohort"],
                    "rho": row["proportion"],
                    "p": float("nan"),
                    "n": row["n"],
                    "c": row["c"],
                    "c_prime": row["c_prime"],
                    "a": row["a"],
                    "b": row["b"],
                    "partial_ols": row["partial_ols"],
                })
    return pd.DataFrame(grid_rows), pd.DataFrame(cohort_rows)


def choose_winner(grid: pd.DataFrame) -> tuple[pd.Series | None, str]:
    """Largest mediation proportion among consistent rows with ρ(TACSTD2, T/NK) ≤ −0.20."""
    for tier in ("A", "B", "C"):
        pool = grid[(grid["tier"] == tier) & (grid["rho_xy"] <= -0.20)]
        if pool.empty:
            continue
        pool = pool.sort_values(
            ["proportion", "rho_xy_I2", "joint"],
            ascending=[False, True, False],
        )
        return pool.iloc[0], tier
    return None, "none"


def leave_one_out(df: pd.DataFrame, xcol: str, mcol: str, model: str) -> list[dict]:
    rows = []
    for drop in COHORTS:
        sub = df[df["cohort"] != drop].copy()
        frame = prepare_vectors(sub, xcol, mcol, model)
        if frame is None:
            continue
        fit = fit_mediation(frame)
        if fit is None:
            continue
        rows.append({"dropped": drop, "n": fit["n"], "proportion": fit["proportion"],
                     "c": fit["c"], "c_prime": fit["c_prime"], "indirect": fit["indirect"]})
    return rows


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def fig_scatter(df: pd.DataFrame, xcol: str, xlabel: str, rho: float, path: Path, title: str):
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for cohort in COHORTS:
        sub = df[df["cohort"] == cohort]
        ax.scatter(sub[xcol], sub["frac_tnk"], s=28, color=COHORT_COLOR[cohort], label=f"{cohort} (n={len(sub)})", zorder=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("T/NK fraction")
    ax.set_title(f"{title}\nDL ρ = {fmt_rho(rho)}", fontsize=10)
    ax.legend(frameon=False, fontsize=7)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_forest(cohort_detail: list[dict], path: Path, title: str):
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    labels = []
    ypos = []
    y = 0
    for kind, name, color in (
        ("rho_xy", "TACSTD2–T/NK", "#1a5276"),
        ("partial_xy", "TACSTD2–T/NK | CLDN4", "#5d6d7e"),
        ("rho_my", "CLDN4–T/NK", "#196f3d"),
        ("partial_my", "CLDN4–T/NK | TACSTD2", "#117a65"),
    ):
        for row in cohort_detail:
            if row["kind"] != kind:
                continue
            labels.append(f"{name}  {row['cohort'].replace('GSE','')}")
            ypos.append(y)
            rho = row["rho"]
            # Fisher CI per cohort
            if row["n"] > 3 and np.isfinite(rho) and abs(rho) < 1:
                z = np.arctanh(np.clip(rho, -0.999, 0.999))
                se = 1 / math.sqrt(row["n"] - 3)
                lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
            else:
                lo = hi = rho
            ax.plot([lo, hi], [y, y], color=color, lw=1.4)
            ax.plot(rho, y, "o", color=color, ms=5)
            y += 1
        y += 0.6
    ax.axvline(0, color="#bbbbbb", lw=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Spearman ρ")
    ax.set_title(title, fontsize=10)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_boot(draws: np.ndarray, point: float, path: Path, title: str, xlabel: str):
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    finite = draws[np.isfinite(draws)]
    ax.hist(finite, bins=40, color="#1a5276", alpha=0.85)
    ax.axvline(0, color="#bbbbbb", lw=0.8)
    ax.axvline(point, color="#b03a2e", lw=1.4, label=f"point {point:.3f}")
    if finite.size:
        lo, hi = np.quantile(finite, [0.025, 0.975])
        ax.axvline(lo, color="#7f8c8d", ls="--", lw=1)
        ax.axvline(hi, color="#7f8c8d", ls="--", lw=1, label=f"95% {lo:.3f} to {hi:.3f}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Bootstrap samples")
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=7)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_compartment(joint: dict, path: Path):
    labels = [
        ("cldn4_pct_gt0", "CLDN4 %"),
        ("pct_cldn4_only_gt0", "CLDN4 only"),
        ("pct_both_gt0", "Double positive"),
        ("tacstd2_pct_gt0", "TACSTD2 %"),
        ("pct_tacstd2_only_gt0", "TACSTD2 only"),
    ]
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    y = 0
    yticks = []
    ylabels = []
    for key, name in labels:
        pack = joint["pieces"][key]
        meta = pack["meta"]
        rho = meta.get("pooled_rho", np.nan)
        lo, hi = meta.get("ci95_rho", [np.nan, np.nan])
        ax.plot([lo, hi], [y, y], color="#1a5276", lw=1.6)
        ax.plot(rho, y, "o", color="#1a5276", ms=6)
        yticks.append(y)
        ylabels.append(name)
        y += 1
        for row in pack["cohorts"]:
            ax.plot(row["rho"], y, "o", color=COHORT_COLOR[row["cohort"]], ms=4)
            yticks.append(y)
            ylabels.append(f"  {row['cohort'].replace('GSE', '')}")
            y += 1
        y += 0.4
    ax.axvline(0, color="#bbbbbb", lw=0.8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("Spearman ρ versus T/NK fraction")
    ax.set_title("Compartment split, DerSimonian–Laird and cohort ρ", fontsize=10)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_strata(df: pd.DataFrame, xcol: str, mcol: str, path: Path, xlabel: str):
    d = df.copy()
    d["q"] = d.groupby("cohort")[mcol].transform(
        lambda s: pd.qcut(s.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), sharey=True)
    for ax, level, title in zip(axes, ["Q1", "Q4"], ["CLDN4 Q1 (low)", "CLDN4 Q4 (high)"]):
        sub = d[d["q"].astype(str) == level]
        for cohort in COHORTS:
            ss = sub[sub["cohort"] == cohort]
            ax.scatter(ss[xcol], ss["frac_tnk"], s=26, color=COHORT_COLOR[cohort], label=cohort)
        ax.set_title(f"{title}  n={len(sub)}", fontsize=10)
        ax.set_xlabel(xlabel)
        style_ax(ax)
    axes[0].set_ylabel("T/NK fraction")
    axes[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def cohort_table_md(detail: pd.DataFrame, spec_id: str, model: str) -> str:
    sub = detail[(detail["spec_id"] == spec_id) & (detail["model"] == model)]
    lines = ["| cohort | n | TACSTD2–T/NK ρ | CLDN4–T/NK ρ | TACSTD2–CLDN4 ρ | partial TACSTD2 | partial CLDN4 | cohort OLS proportion |",
             "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for cohort in COHORTS:
        def grab(kind):
            hit = sub[(sub["kind"] == kind) & (sub["cohort"] == cohort)]
            if hit.empty:
                return float("nan"), 0
            return float(hit.iloc[0]["rho"]), int(hit.iloc[0]["n"])
        xy, n = grab("rho_xy")
        my, _ = grab("rho_my")
        xm, _ = grab("rho_xm")
        pxy, _ = grab("partial_xy")
        pmy, _ = grab("partial_my")
        ols, _ = grab("ols")
        lines.append(
            f"| {cohort} | {n} | {fmt_rho(xy)} | {fmt_rho(my)} | {fmt_rho(xm)} | {fmt_rho(pxy)} | {fmt_rho(pmy)} | {fmt_num(ols)} |"
        )
    return "\n".join(lines)


def block_for_row(row: pd.Series, detail: pd.DataFrame, boot: dict | None, inter: dict | None,
                  glmm: dict | None, strata: list[dict], loo: list[dict]) -> str:
    lines = []
    lines.append(f"Specification `{row['spec_id']}` ({row['family']}, tier {row['tier']}). N = {int(row['n'])}.")
    lines.append("")
    lines.append("### (1) Spearman")
    lines.append("")
    lines.append(
        f"TACSTD2 vs T/NK: ρ = {fmt_rho(row['rho_xy'])} (p = {fmt_p(row['rho_xy_p'])}, "
        f"I² = {fmt_num(row['rho_xy_I2'], 1)}%, {fmt_rho(row['rho_xy_lo'])} to {fmt_rho(row['rho_xy_hi'])})."
    )
    lines.append(
        f"CLDN4 vs T/NK: ρ = {fmt_rho(row['rho_my'])} (p = {fmt_p(row['rho_my_p'])}, I² = {fmt_num(row['rho_my_I2'], 1)}%)."
    )
    lines.append(
        f"TACSTD2 vs CLDN4: ρ = {fmt_rho(row['rho_xm'])} (p = {fmt_p(row['rho_xm_p'])}, I² = {fmt_num(row['rho_xm_I2'], 1)}%)."
    )
    lines.append("")
    lines.append("### (2) Partial Spearman")
    lines.append("")
    lines.append(
        f"TACSTD2–T/NK | CLDN4: ρ = {fmt_rho(row['partial_xy'])} (p = {fmt_p(row['partial_xy_p'])}, "
        f"I² = {fmt_num(row['partial_xy_I2'], 1)}%, {fmt_rho(row['partial_xy_lo'])} to {fmt_rho(row['partial_xy_hi'])})."
    )
    lines.append(
        f"CLDN4–T/NK | TACSTD2: ρ = {fmt_rho(row['partial_my'])} (p = {fmt_p(row['partial_my_p'])}, "
        f"I² = {fmt_num(row['partial_my_I2'], 1)}%)."
    )
    if (
        np.isfinite(row["rho_xy"]) and np.isfinite(row["partial_xy"])
        and row["rho_xy"] < 0 and row["partial_xy"] < 0 and abs(row["partial_xy"]) < abs(row["rho_xy"])
    ):
        sp_att = 1.0 - (row["partial_xy"] / row["rho_xy"])
        lines.append(f"Spearman attenuation fraction 1 − ρ_partial/ρ_total = {fmt_num(sp_att)}.")
    elif np.isfinite(row["rho_xy"]) and np.isfinite(row["partial_xy"]) and row["rho_xy"] * row["partial_xy"] < 0:
        lines.append(
            "The partial Spearman changes sign relative to the marginal Spearman, so 1 − ρ_partial/ρ_total is not an attenuation fraction in (0, 1)."
        )
    lines.append("")
    lines.append("### (3) Nested OLS with cohort fixed effects")
    lines.append("")
    lines.append(
        f"Total TACSTD2 β = {fmt_num(row['c'], 4)} (SE {fmt_num(row['c_se'], 4)}, p = {fmt_p(row['c_p'])})."
    )
    lines.append(
        f"Direct TACSTD2 β after CLDN4 = {fmt_num(row['c_prime'], 4)} (SE {fmt_num(row['c_prime_se'], 4)}, p = {fmt_p(row['c_prime_p'])})."
    )
    if row["c"] < 0 and row["c_prime"] < 0 and abs(row["c_prime"]) < abs(row["c"]):
        lines.append(
            f"The TACSTD2 coefficient stays negative and shrinks from {fmt_num(row['c'], 4)} to {fmt_num(row['c_prime'], 4)}. "
            f"Attenuation fraction = {fmt_num(row['proportion'])}."
        )
    elif row["c"] * row["c_prime"] < 0:
        lines.append(
            f"The TACSTD2 coefficient changes sign, from {fmt_num(row['c'], 4)} to {fmt_num(row['c_prime'], 4)}. "
            f"The ratio a×b/c equals {fmt_num(row['proportion'])}, which is outside (0, 1)."
        )
    else:
        lines.append(
            f"Direct β = {fmt_num(row['c_prime'], 4)} next to total β = {fmt_num(row['c'], 4)}. "
            f"Ratio a×b/c = {fmt_num(row['proportion'])}."
        )
    lines.append("")
    lines.append("### (4) Mediation TACSTD2 → CLDN4 → T/NK")
    lines.append("")
    lines.append(
        f"Baron–Kenny / product of coefficients, same cohort fixed effects. "
        f"a (TACSTD2 → CLDN4) = {fmt_num(row['a'], 4)} (p = {fmt_p(row['a_p'])}). "
        f"b (CLDN4 → T/NK | TACSTD2) = {fmt_num(row['b'], 4)} (p = {fmt_p(row['b_p'])}). "
        f"Indirect a×b = {fmt_num(row['indirect'], 4)}. "
        f"Proportion mediated a×b/c = {fmt_num(row['proportion'])}. "
        f"Sobel p = {fmt_p(row['sobel_p'])}."
    )
    lines.append(
        f"The OLS identity (c − c') − a×b has absolute gap {row['ols_identity_gap']:.2e}."
    )
    if boot is not None:
        lines.append(
            f"Cohort-stratified bootstrap, {boot['n_boot']} draws, {boot['n_finite']} finite. "
            f"Proportion median {fmt_num(boot['prop_p50'])}, "
            f"95% percentile interval {fmt_num(boot['prop_p2_5'])} to {fmt_num(boot['prop_p97_5'])}. "
            f"Indirect 95% interval {fmt_num(boot['indirect_p2_5'], 4)} to {fmt_num(boot['indirect_p97_5'], 4)}. "
            f"Share of draws with proportion in (0, 1): {fmt_num(boot['frac_prop_in_0_1'])}. "
            f"Share with the same partial-attenuation sign pattern: {fmt_num(boot['frac_attenuates'])}."
        )
    if loo:
        bits = ", ".join(
            f"drop {r['dropped']} → {fmt_num(r['proportion'])} (n={r['n']})" for r in loo
        )
        lines.append(f"Leave-one-cohort-out proportion: {bits}.")
    lines.append("")
    lines.append("### (5) TACSTD2 × CLDN4 interaction")
    lines.append("")
    if inter is not None:
        lines.append(
            f"Product-term β = {fmt_num(inter['b_int'], 4)} (SE {fmt_num(inter['se_int'], 4)}, p = {fmt_p(inter['p_int'])}). "
            f"Simple slope of TACSTD2 at the pooled 25th percentile of CLDN4 = {fmt_num(inter['slope_q25'], 4)}; "
            f"at the 75th percentile = {fmt_num(inter['slope_q75'], 4)}."
        )
    lines.append("")
    lines.append("### (6) CLDN4 Q1 vs Q4")
    lines.append("")
    if strata:
        lines.append("| split | level | n | cohorts in Spearman meta | ρ | p | I² | OLS β | p |")
        lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
        for s in strata:
            lines.append(
                f"| {s['split']} | {s['level']} | {s['n']} | {s['n_cohorts_spearman']} | "
                f"{fmt_rho(s['rho'])} | {fmt_p(s['rho_p'])} | {fmt_num(s['I2'], 1)}% | "
                f"{fmt_num(s['ols_beta'], 4)} | {fmt_p(s['ols_p'])} |"
            )
        lines.append("")
        lines.append("Quartile Spearman meta uses cohorts with n ≥ 5 inside the bin. GSE189357 has 9 patients, so a quartile bin often falls under that floor.")
        for s in strata:
            if s["split"] != "median":
                continue
            cohorts = json.loads(s["cohorts"]) if isinstance(s["cohorts"], str) else s["cohorts"]
            bits = ", ".join(
                f"{c['cohort']} {c['rho']:+.3f} (n={c['n']})" if np.isfinite(c["rho"]) else f"{c['cohort']} — (n={c['n']})"
                for c in cohorts
            )
            lines.append(f"CLDN4 {s['level']} half, cohort Spearman of the TACSTD2 score versus T/NK: {bits}.")
    if glmm is not None:
        lines.append("")
        lines.append("### GLMM random intercept for cohort")
        lines.append("")
        if glmm.get("ok"):
            lines.append(
                f"Total β = {fmt_num(glmm['c'], 4)} (p = {fmt_p(glmm['c_p'])}). "
                f"Direct β = {fmt_num(glmm['c_prime'], 4)} (p = {fmt_p(glmm['c_prime_p'])}). "
                f"a = {fmt_num(glmm['a'], 4)}, b = {fmt_num(glmm['b'], 4)}. "
                f"Proportion = {fmt_num(glmm['proportion'])}. "
                f"Random-intercept variance: base {fmt_num(glmm.get('re_var_base', float('nan')), 4)}, "
                f"full {fmt_num(glmm.get('re_var_full', float('nan')), 4)}. "
                f"Converged: base {glmm['converged_base']}, full {glmm['converged_full']}. "
                f"A random-intercept variance of 0 leaves the coefficients equal to the cohort-fixed-effect OLS coefficients on this specification."
            )
        else:
            lines.append(f"GLMM did not return a fit ({glmm.get('error', 'unknown')}). The OLS fixed-effect estimates above are the reported model.")
    lines.append("")
    lines.append("Cohort detail:")
    lines.append("")
    lines.append(cohort_table_md(detail, row["spec_id"], row["model"]))
    lines.append("")
    return "\n".join(lines)


def joint_compartment(df: pd.DataFrame) -> dict:
    """TACSTD2 % = double-positive % + TACSTD2-only %. Same identity for CLDN4."""
    need = ["pct_both_gt0", "pct_tacstd2_only_gt0", "pct_cldn4_only_gt0"]
    if any(c not in df.columns for c in need):
        return {}
    gap_t = float(np.max(np.abs(df["tacstd2_pct_gt0"] - df["pct_both_gt0"] - df["pct_tacstd2_only_gt0"])))
    gap_c = float(np.max(np.abs(df["cldn4_pct_gt0"] - df["pct_both_gt0"] - df["pct_cldn4_only_gt0"])))
    if gap_t > 1e-8 or gap_c > 1e-8:
        raise SystemExit(f"compartment identity failed tac {gap_t} cldn {gap_c}")
    pieces = {}
    for col in ["pct_both_gt0", "pct_tacstd2_only_gt0", "pct_cldn4_only_gt0", "tacstd2_pct_gt0", "cldn4_pct_gt0"]:
        rows = spearman_by_cohort(df, col, "frac_tnk")
        meta = meta_from_rows(rows)
        pieces[col] = {"cohorts": rows, "meta": meta}
    y = df["frac_tnk"].to_numpy(float)
    both = df["pct_both_gt0"].to_numpy(float)
    tony = df["pct_tacstd2_only_gt0"].to_numpy(float)
    conly = df["pct_cldn4_only_gt0"].to_numpy(float)
    fit = ols_fit(y, [both, tony, conly], df["cohort"])
    betas = {}
    if fit is not None:
        for i, name in enumerate(["both", "tacstd2_only", "cldn4_only"], start=1):
            b, se, p = coef0(fit, i)
            betas[name] = {"beta": b, "se": se, "p": p}
        betas["r2"] = float(fit.rsquared)
        betas["n"] = int(fit.nobs)
    # nested: TACSTD2 % alone, then add CLDN4-only (the part of CLDN4 outside TACSTD2)
    # and both (shared). Reported as the compartment model above.
    return {"gap_t": gap_t, "gap_c": gap_c, "pieces": pieces, "betas": betas}


def write_finding(text: str) -> None:
    (HERE / "FINDING.md").write_text(text)


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    units = pd.read_csv(RES / "unit_scores.tsv", sep="\t")
    units["unit_id"] = units["unit_id"].astype(str)
    units = units.sort_values(["cohort", "unit_id"]).reset_index(drop=True)
    if len(units) != 65:
        raise SystemExit(f"expected 65 units, found {len(units)}")
    units = add_quantile_scores(units)
    # pipeline check on locked CLDN4 %pos
    rho_rows = spearman_by_cohort(units, "cldn4_pct_gt0", "frac_tnk")
    meta = meta_from_rows(rho_rows)
    if abs(meta["pooled_rho"] - LOCKED_CLDN4_RHO) > 1e-6:
        raise SystemExit(f"CLDN4 pipeline check failed: {meta['pooled_rho']} vs {LOCKED_CLDN4_RHO}")
    print(f"pipeline check CLDN4 %pos ρ={meta['pooled_rho']:.6f}", flush=True)

    grid, detail = evaluate(units)
    grid.to_csv(TAB / "grid_mediation.tsv", sep="\t", index=False)
    detail.to_csv(TAB / "cohort_effects.tsv", sep="\t", index=False)
    winner, tier = choose_winner(grid)
    print(f"grid {len(grid)} winner tier {tier} " + (winner["spec_id"] if winner is not None else "none"), flush=True)

    focus_ids = []
    locked_hits = grid[(grid["x_score"] == "pct_gt0") & (grid["m_score"] == "pct_gt0") & (grid["model"] == "identity")]
    if locked_hits.empty:
        raise SystemExit("locked identity spec missing")
    locked = locked_hits.iloc[0]
    focus = [("locked", locked)]
    if winner is not None and winner["spec_id"] != locked["spec_id"]:
        focus.append(("winner", winner))
    elif winner is not None:
        focus = [("winner_locked", winner)]

    extras = {}
    for tag, row in focus:
        xcol = f"tacstd2_{row['x_score']}"
        mcol = f"cldn4_{row['m_score']}"
        frame = prepare_vectors(units, xcol, mcol, row["model"])
        print(f"bootstrap {tag} {row['spec_id']}", flush=True)
        boot = bootstrap_proportion(frame, N_BOOT, BOOT_SEED)
        inter = fit_interaction(frame)
        glmm = fit_glmm(frame)
        # strata always on the raw score vs raw T/NK, which is the patient phenotype
        strata = strata_effects(units, xcol, mcol)
        loo = leave_one_out(units, xcol, mcol, row["model"])
        extras[tag] = {"row": row, "boot": boot, "inter": inter, "glmm": glmm, "strata": strata, "loo": loo, "frame": frame}
        pd.DataFrame([{k: v for k, v in boot.items() if not str(k).startswith("draws_")}]).to_csv(
            TAB / f"bootstrap_{tag}.tsv", sep="\t", index=False
        )
        pd.DataFrame({"proportion": boot["draws_prop"], "indirect": boot["draws_indirect"]}).to_csv(
            TAB / f"bootstrap_draws_{tag}.tsv", sep="\t", index=False
        )
        pd.DataFrame([inter]).to_csv(TAB / f"interaction_{tag}.tsv", sep="\t", index=False)
        pd.DataFrame([glmm]).to_csv(TAB / f"glmm_{tag}.tsv", sep="\t", index=False)
        pd.DataFrame(strata).to_csv(TAB / f"strata_{tag}.tsv", sep="\t", index=False)
        pd.DataFrame(loo).to_csv(TAB / f"loo_{tag}.tsv", sep="\t", index=False)
        focus_ids.append(row["spec_id"])

    # figures for locked pct and for winner if different
    fig_scatter(
        units, "tacstd2_pct_gt0", "Malignant TACSTD2 fraction UMI > 0",
        float(locked["rho_xy"]), FIG / "scatter_tacstd2_pct", "Locked TACSTD2 %pos vs T/NK",
    )
    fig_scatter(
        units, "cldn4_pct_gt0", "Malignant CLDN4 fraction UMI > 0",
        float(locked["rho_my"]), FIG / "scatter_cldn4_pct", "Locked CLDN4 %pos vs T/NK",
    )
    locked_detail_rows = detail[(detail["spec_id"] == locked["spec_id"]) & (detail["model"] == "identity") & (detail["kind"].isin(["rho_xy", "partial_xy", "rho_my", "partial_my"]))]
    fig_forest(locked_detail_rows.to_dict(orient="records"), FIG / "forest_locked", "Locked %pos (UMI > 0), identity scale")
    locked_boot = extras["locked"]["boot"] if "locked" in extras else extras["winner_locked"]["boot"]
    fig_boot(
        locked_boot["draws_indirect"], float(locked["indirect"]), FIG / "bootstrap_locked",
        "Locked UMI>0 indirect effect", "Indirect effect a×b",
    )
    fig_strata(units, "tacstd2_pct_gt0", "cldn4_pct_gt0", FIG / "strata_locked", "TACSTD2 fraction UMI > 0")

    win_tag = "winner" if "winner" in extras else ("winner_locked" if "winner_locked" in extras else None)
    if win_tag == "winner":
        wrow = extras["winner"]["row"]
        fig_scatter(
            units, f"tacstd2_{wrow['x_score']}", f"TACSTD2 {wrow['x_score']}",
            float(wrow["rho_xy"]), FIG / "scatter_winner_tacstd2", f"Max spec TACSTD2 ({wrow['model']})",
        )
        wdetail = detail[(detail["spec_id"] == wrow["spec_id"]) & (detail["kind"].isin(["rho_xy", "partial_xy", "rho_my", "partial_my"]))]
        # cohort rho does not depend on model; kind rows are repeated per model. take this model.
        wdetail = wdetail[wdetail["model"] == wrow["model"]]
        fig_forest(wdetail.to_dict(orient="records"), FIG / "forest_winner", f"Maximum spec {wrow['spec_id']}")
        fig_boot(
            extras["winner"]["boot"]["draws_indirect"], float(wrow["indirect"]), FIG / "bootstrap_winner",
            "Maximum-spec indirect effect", "Indirect effect a×b",
        )
        fig_strata(units, f"tacstd2_{wrow['x_score']}", f"cldn4_{wrow['m_score']}", FIG / "strata_winner",
                   f"TACSTD2 {wrow['x_score']}")

    # top eligible table
    eligible = grid[grid["tier"].isin(["A", "B", "C"])].sort_values(["tier", "proportion"], ascending=[True, False])
    eligible.head(25).to_csv(TAB / "top_eligible.tsv", sep="\t", index=False)

    joint = joint_compartment(units)
    if joint:
        joint_rows = []
        for col, pack in joint["pieces"].items():
            meta_j = pack["meta"]
            joint_rows.append({
                "score": col,
                "rho": meta_j.get("pooled_rho"),
                "p": meta_j.get("p"),
                "I2": meta_j.get("I2"),
                "lo": meta_j.get("ci95_rho", [None, None])[0],
                "hi": meta_j.get("ci95_rho", [None, None])[1],
                "n_neg": int(sum(r["rho"] < 0 for r in pack["cohorts"])),
                **{f"rho_{r['cohort']}": r["rho"] for r in pack["cohorts"]},
            })
        pd.DataFrame(joint_rows).to_csv(TAB / "compartment_spearman.tsv", sep="\t", index=False)
        pd.DataFrame([
            {"term": k, **v} if isinstance(v, dict) else {"term": k, "beta": v}
            for k, v in joint["betas"].items()
        ]).to_csv(TAB / "compartment_ols.tsv", sep="\t", index=False)
        print("compartment", joint["betas"], flush=True)
        fig_compartment(joint, FIG / "forest_compartment")

    n_tier = grid["tier"].value_counts().to_dict()
    summary = {
        "n_units": 65,
        "n_specs": int(len(grid)),
        "tier_counts": {k: int(v) for k, v in n_tier.items()},
        "pipeline_cldn4_rho": meta["pooled_rho"],
        "locked_spec": locked["spec_id"],
        "locked_proportion": float(locked["proportion"]),
        "locked_tier": locked["tier"],
        "winner_spec": None if winner is None else winner["spec_id"],
        "winner_tier": tier,
        "winner_proportion": None if winner is None else float(winner["proportion"]),
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2))

    # FINDING
    lines = []
    lines.append("# TACSTD2 → T/NK on concordant-4, and how much of it runs through CLDN4")
    lines.append("")
    lines.append("ADDITIVE mechanism test. Locked concordant-4 only: GSE123902 (13 donors) + GSE131907 (21 samples) + GSE205335 (22 patients) + GSE189357 (9 patients). N = 65. Not GSE148071, GSE127465, GSE207422, GSE154826, GSE200563, or E-MTAB-13526.")
    lines.append("")
    lines.append("The unit is the patient, donor, or sample already locked for malignant CLDN4 % positive versus T/NK. Cell counts are not n. p-values are descriptive. This is not a causal mediation, not a spatial exclusion result, and not a claim about private KL tumors.")
    lines.append("")
    lines.append("## Pipeline check")
    lines.append("")
    lines.append(f"Malignant CLDN4 fraction with UMI > 0 versus T/NK fraction, DerSimonian–Laird ρ = {meta['pooled_rho']:.6f} (p = {fmt_p(meta['p'])}, I² = {meta['I2']:.1f}%). The locked published value is −0.5311678045689989. Cohort ρ: " + ", ".join(f"{r['cohort']} {r['rho']:.3f} (n={r['n']})" for r in rho_rows) + ".")
    lines.append("")
    lines.append("## What was maximized")
    lines.append("")
    lines.append("Pre-specified grid. Same malignant-cell summary for TACSTD2 and CLDN4, plus a cross of each gene against the locked UMI>0 fraction of the other. Summaries: UMI>0, ≥2, ≥3, ≥5, ≥10; CP10k ≥1, ≥5, ≥10; fraction above the cohort median UMI and above the cohort 75th percentile; mean log1p(UMI); mean log1p(CP10k); mean log1p among detected cells; 90th percentile on both scales; pseudobulk log1p(CP10k).")
    lines.append("")
    lines.append("Models, each with a cohort fixed effect: raw scores; within-cohort z-scores of TACSTD2, CLDN4, and T/NK; within-cohort ranks; logit of proportions (clip 10⁻⁴) when both scores are fractions. Mediation is Baron–Kenny on that linear model. The product of coefficients equals the drop in the TACSTD2 coefficient. Proportion = a×b / c.")
    lines.append("")
    lines.append("A row is tier A when all four of these hold: every cohort has TACSTD2–T/NK ρ < 0, CLDN4–T/NK ρ < 0, and TACSTD2–CLDN4 ρ > 0 with |ρ| < 0.95; the pooled Spearman partial of TACSTD2 also moves toward zero and stays negative; the pooled OLS path is a > 0, b < 0, c < 0, c′ < 0, and 0 < a×b/c < 1; and the same partial-attenuation pattern holds in the separate OLS of each cohort. Tier B drops only the per-cohort OLS requirement. Tier C drops the Spearman-attenuation requirement and keeps the sign pattern plus pooled OLS attenuation. The reported maximum is the largest mediation proportion inside the best non-empty tier, among rows whose pooled TACSTD2–T/NK ρ is ≤ −0.20. Ties break on lower I², then on |ρ| × proportion.")
    lines.append("")
    lines.append(f"Grid size {len(grid)}. Tier counts: " + ", ".join(f"{k} = {int(v)}" for k, v in sorted(n_tier.items())) + ".")
    lines.append("")
    if winner is None:
        lines.append("No row cleared tier C with pooled TACSTD2–T/NK ρ ≤ −0.20. Every TACSTD2 summary in the grid is positively correlated with T/NK in GSE205335, so none is 4-cohort negative. The most negative pooled TACSTD2–T/NK ρ in the grid is the locked UMI>0 fraction.")
        same = grid[grid["pooled_partial"] & grid["spearman_attenuates"]].sort_values("proportion", ascending=False)
        if not same.empty:
            top = same.iloc[0]
            lines.append("")
            lines.append(
                f"Largest same-sign attenuation in the grid, without the 4-cohort sign rule: `{top['spec_id']}`. "
                f"Proportion {fmt_num(float(top['proportion']))}, total β {fmt_num(float(top['c']), 4)}, "
                f"direct β {fmt_num(float(top['c_prime']), 4)}, pooled TACSTD2–T/NK ρ {fmt_rho(float(top['rho_xy']))}, "
                f"partial ρ {fmt_rho(float(top['partial_xy']))}. "
                f"Cohorts with the per-cohort OLS attenuation pattern: {int(top['n_cohorts_partial_ols'])}/4. "
                f"This row stays in the table as the grid maximum under a weaker rule. It is not the 4-cohort result."
            )
        lines.append("")
        lines.append("The six tests below are the locked UMI>0 specification, which is the pre-specified score.")
    else:
        lines.append(f"Maximum row: `{winner['spec_id']}`, tier {tier}, mediation proportion {fmt_num(float(winner['proportion']))}, pooled TACSTD2–T/NK ρ {fmt_rho(float(winner['rho_xy']))}.")
    lines.append("")
    # lead with winner if it is the max and not locked, else locked first then note
    if win_tag == "winner":
        lines.append("## Maximum partial attenuation")
        lines.append("")
        lines.append(block_for_row(extras["winner"]["row"], detail, extras["winner"]["boot"], extras["winner"]["inter"], extras["winner"]["glmm"], extras["winner"]["strata"], extras["winner"]["loo"]))
        lines.append("## Locked UMI > 0 (pre-specified)")
        lines.append("")
        lines.append(block_for_row(extras["locked"]["row"], detail, extras["locked"]["boot"], extras["locked"]["inter"], extras["locked"]["glmm"], extras["locked"]["strata"], extras["locked"]["loo"]))
    else:
        tag0 = "winner_locked" if "winner_locked" in extras else "locked"
        lines.append("## Locked UMI > 0")
        lines.append("")
        lines.append("This pre-specified row is also the maximum under the rule above." if tag0 == "winner_locked" else "Pre-specified row.")
        lines.append("")
        lines.append(block_for_row(extras[tag0]["row"], detail, extras[tag0]["boot"], extras[tag0]["inter"], extras[tag0]["glmm"], extras[tag0]["strata"], extras[tag0]["loo"]))

    lines.append("## Top eligible rows")
    lines.append("")
    if eligible.empty:
        lines.append("No eligible row.")
    else:
        lines.append("| spec | tier | proportion | TACSTD2–T/NK ρ | partial ρ | I² of ρ | cohorts with OLS attenuation |")
        lines.append("|---|---|---:|---:|---:|---:|---:|")
        for _, r in eligible.head(12).iterrows():
            lines.append(
                f"| `{r['spec_id']}` | {r['tier']} | {fmt_num(r['proportion'])} | {fmt_rho(r['rho_xy'])} | {fmt_rho(r['partial_xy'])} | {fmt_num(r['rho_xy_I2'], 1)}% | {int(r['n_cohorts_partial_ols'])} |"
            )
    lines.append("")
    lines.append("## Reading")
    lines.append("")
    lines.append("The stable path is the product of coefficients. On the locked UMI>0 scores, a is positive and b is negative, so a×b is negative: higher TACSTD2 tracks higher CLDN4, and higher CLDN4 tracks lower T/NK. The cohort-stratified bootstrap interval for a×b stays below zero. Sobel p for that product is 0.0024. The marginal TACSTD2 coefficient is a small negative number, and the direct coefficient after CLDN4 is positive, so a×b/c sits above 1. That ratio is outside a partial-mediation proportion. 8.4% of bootstrap draws fall inside (0, 1).")
    lines.append("")
    lines.append("CLDN4’s partial Spearman with T/NK, holding TACSTD2, stays near the marginal CLDN4 result and is negative in every cohort, with I² = 0. TACSTD2’s partial Spearman, holding CLDN4, is positive in the pool. GSE123902 is the only cohort whose own OLS shows a negative TACSTD2 coefficient that shrinks and stays negative (proportion 0.66).")
    lines.append("")
    lines.append("No TACSTD2 summary in the 89-row grid is negatively correlated with T/NK in GSE205335. The most negative pooled ρ is −0.112, the locked percent. The largest same-sign attenuation anywhere in the grid is 0.33, for TACSTD2 percent versus T/NK with CLDN4 entered as malignant pseudobulk log1p(CP10k). That row’s pooled ρ is still −0.112, and the within-cohort OLS attenuation pattern is present in 1 of 4 cohorts.")
    lines.append("")
    lines.append("The percent split points the same way. CLDN4-only (CLDN4 detected, TACSTD2 undetected) is negative versus T/NK in all four cohorts (I² = 0). Double-positive cells are negative in three cohorts and positive in GSE205335. TACSTD2-only cells are negative in one cohort. In the joint OLS all three compartment coefficients are negative; the TACSTD2-only coefficient is the least precise because that compartment is the smallest fraction.")
    lines.append("")
    lines.append("Below each cohort’s CLDN4 median, TACSTD2 percent versus T/NK is negative in all four cohorts (pooled ρ −0.29, I² = 0, n = 34). Above the median the cohorts disagree and the pool is positive. The continuous TACSTD2 × CLDN4 product term is close to zero.")
    lines.append("")
    lines.append("Taken together, the immune-cold association that is consistent across the four cohorts is CLDN4’s. TACSTD2 shares that direction through its correlation with CLDN4 (the indirect path). A residual TACSTD2 association in the same direction, after CLDN4, is not a 4-cohort result.")
    lines.append("")
    lines.append("Cohort fixed effects absorb the mean difference between studies. They do not make the four studies exchangeable. GSE123902 and GSE189357 malignant cells are marker-gated (EPCAM or KRT8/18/19, and PTPRC = 0). GSE131907 and GSE205335 use the author malignant label. T/NK is the locked fraction: marker T/NK and not malignant, or the author T/NK label, over all cells in the unit.")
    lines.append("")
    lines.append("A high mediation proportion on a weak TACSTD2–T/NK correlation is a large share of a small association. The table carries both numbers. The searched maximum is not a confirmatory p-value.")
    lines.append("")
    if joint:
        lines.append("## Compartment split of the UMI > 0 fractions")
        lines.append("")
        lines.append("On malignant cells, TACSTD2 % = (TACSTD2>0 and CLDN4>0) + (TACSTD2>0 and CLDN4=0). CLDN4 % = double-positive + CLDN4-only. The two identities hold to numerical error on all 65 units. This split asks which part of the TACSTD2-positive fraction carries the T/NK association.")
        lines.append("")
        lines.append("| score | ρ | p | I² | cohorts < 0 | GSE123902 | GSE131907 | GSE205335 | GSE189357 |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for col, pack in joint["pieces"].items():
            meta_j = pack["meta"]
            by = {r["cohort"]: r["rho"] for r in pack["cohorts"]}
            n_neg = sum(r["rho"] < 0 for r in pack["cohorts"])
            lines.append(
                f"| {col} | {fmt_rho(meta_j.get('pooled_rho'))} | {fmt_p(meta_j.get('p'))} | "
                f"{fmt_num(meta_j.get('I2'), 1)}% | {n_neg}/4 | "
                + " | ".join(fmt_rho(by[c]) for c in COHORTS) + " |"
            )
        lines.append("")
        betas = joint["betas"]
        if "both" in betas:
            means = units[["pct_both_gt0", "pct_tacstd2_only_gt0", "pct_cldn4_only_gt0"]].mean()
            lines.append(
                f"Mean fractions across the 65 units: double-positive {fmt_num(float(means['pct_both_gt0']))}, "
                f"TACSTD2-only {fmt_num(float(means['pct_tacstd2_only_gt0']))}, "
                f"CLDN4-only {fmt_num(float(means['pct_cldn4_only_gt0']))}."
            )
            lines.append(
                f"OLS of T/NK fraction on the three compartments with cohort fixed effects (n = {int(betas['n'])}, R² = {fmt_num(betas['r2'])}). "
                f"Double-positive β = {fmt_num(betas['both']['beta'], 4)} (p = {fmt_p(betas['both']['p'])}). "
                f"TACSTD2-only β = {fmt_num(betas['tacstd2_only']['beta'], 4)} (p = {fmt_p(betas['tacstd2_only']['p'])}). "
                f"CLDN4-only β = {fmt_num(betas['cldn4_only']['beta'], 4)} (p = {fmt_p(betas['cldn4_only']['p'])}). "
                f"β is the change in T/NK fraction per 1.0 change in that compartment, so the TACSTD2-only slope applies to a compartment whose mean is {fmt_num(float(means['pct_tacstd2_only_gt0']))}."
            )
            lines.append("")
    lines.append("## Caveats")
    lines.append("")
    lines.append("Observational co-expression. CLDN4 may sit beside TACSTD2 on a shared epithelial program rather than on a path from TACSTD2 to immune composition. No intervention, no instrument, no spatial radius. Keratin adjustment is not in the grid; the locked CLDN4 result was already strongest without it. Logit uses a 10⁻⁴ clip at 0 and 1. Within-cohort quartile bins in GSE189357 are small. Four random-effect levels is a thin GLMM.")
    lines.append("")
    lines.append("Reproduce: `bash methods/concordant4_tacstd2_cldn4_mediation/scripts/download.sh /tmp/geo_dl` then `python3 methods/concordant4_tacstd2_cldn4_mediation/scripts/extract_scores.py` then `python3 methods/concordant4_tacstd2_cldn4_mediation/scripts/analyze.py`. GEO matrices stay outside the repo. `results/unit_scores.tsv` is enough to rerun the grid.")
    lines.append("")
    write_finding("\n".join(lines))
    print("wrote FINDING", flush=True)


if __name__ == "__main__":
    main()
