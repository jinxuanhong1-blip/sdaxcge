"""Gaussian SEM likelihoods for three observed variables.

The three graphs have the same number of free parameters. AIC and BIC
therefore rank them by likelihood. Means are removed before fitting
(centered Gaussian). Variances use the maximum-likelihood divisor n.

Columns of X are always (exposure-or-first, second, outcome) as documented
by the caller. This module names them T, C, and I:

  T  TACSTD2
  C  CLDN4
  I  immune
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def center(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[1] != 3:
        raise ValueError(f"X must be n x 3, got {X.shape}")
    if X.shape[0] < 8:
        raise ValueError(f"need at least 8 rows, got {X.shape[0]}")
    if not np.isfinite(X).all():
        raise ValueError("X contains non-finite values")
    return X - X.mean(axis=0)


def _var(x: np.ndarray) -> float:
    return float(np.dot(x, x) / x.size)


def _ll_var(v: float, n: int) -> float:
    if v <= 0:
        return -np.inf
    return -0.5 * n * (np.log(2.0 * np.pi) + np.log(v) + 1.0)


def _reg(y: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, float, float]:
    """OLS through the origin. Returns beta, MLE residual variance, Gaussian ll."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    n = y.size
    ve = float(np.dot(resid, resid) / n)
    return beta, ve, _ll_var(ve, n)


def mvn_ll(Xc: np.ndarray, sigma: np.ndarray) -> float:
    """Centered multivariate-normal log-likelihood at covariance sigma."""
    n, p = Xc.shape
    sign, logdet = np.linalg.slogdet(sigma)
    if sign <= 0 or not np.isfinite(logdet):
        return -np.inf
    sample = (Xc.T @ Xc) / n
    inv = np.linalg.inv(sigma)
    quad = float(np.sum(inv * sample))
    return -0.5 * n * (p * np.log(2.0 * np.pi) + logdet + quad)


def _pack(name: str, ll: float, k: int, n: int, paths: dict, sigma: np.ndarray, Xc: np.ndarray) -> dict:
    sigma = np.asarray(sigma, dtype=float)
    sigma = 0.5 * (sigma + sigma.T)
    aic = -2.0 * ll + 2.0 * k
    bic = -2.0 * ll + k * np.log(n)
    aicc = aic + (2.0 * k * (k + 1.0) / (n - k - 1.0))
    mvn = mvn_ll(Xc, sigma)
    return {
        "model": name,
        "ll": float(ll),
        "k": int(k),
        "n": int(n),
        "aic": float(aic),
        "bic": float(bic),
        "aicc": float(aicc),
        "ll_mvn": float(mvn),
        "paths": paths,
        "sigma": sigma,
    }


def fit_graphs(X: np.ndarray) -> dict[str, dict]:
    """Fit M1, M2, M3, and the saturated linear SEM.

    M1  TACSTD2 → CLDN4 → immune, no direct TACSTD2 → immune
    M2  CLDN4 → TACSTD2 → immune, no direct CLDN4 → immune
    M3  TACSTD2 → immune ← CLDN4, with TACSTD2 independent of CLDN4
    SAT saturated chain with the direct path (6 covariance parameters)
    """
    Xc = center(X)
    T, C, I = Xc[:, 0], Xc[:, 1], Xc[:, 2]
    n = int(T.size)
    vt = _var(T)
    vc = _var(C)
    ll_t = _ll_var(vt, n)
    ll_c_marg = _ll_var(vc, n)

    # M1: C on T, I on C.
    a1, ve_c1, ll_c1 = _reg(C, T[:, None])
    b1, ve_i1, ll_i1 = _reg(I, C[:, None])
    a1, b1 = float(a1[0]), float(b1[0])
    var_c1 = a1 * a1 * vt + ve_c1
    sigma1 = np.array(
        [
            [vt, a1 * vt, a1 * b1 * vt],
            [a1 * vt, var_c1, b1 * var_c1],
            [a1 * b1 * vt, b1 * var_c1, b1 * b1 * var_c1 + ve_i1],
        ]
    )
    m1 = _pack(
        "M1_TACSTD2_to_CLDN4_to_immune",
        ll_t + ll_c1 + ll_i1,
        5,
        n,
        {"a_CLDN4_on_TACSTD2": a1, "b_immune_on_CLDN4": b1, "var_T": vt, "var_C_resid": ve_c1, "var_I_resid": ve_i1},
        sigma1,
        Xc,
    )

    # M2: T on C, I on T.
    a2, ve_t2, ll_t2 = _reg(T, C[:, None])
    b2, ve_i2, ll_i2 = _reg(I, T[:, None])
    a2, b2 = float(a2[0]), float(b2[0])
    var_t2 = a2 * a2 * vc + ve_t2
    sigma2 = np.array(
        [
            [var_t2, a2 * vc, b2 * var_t2],
            [a2 * vc, vc, a2 * b2 * vc],
            [b2 * var_t2, a2 * b2 * vc, b2 * b2 * var_t2 + ve_i2],
        ]
    )
    m2 = _pack(
        "M2_CLDN4_to_TACSTD2_to_immune",
        ll_c_marg + ll_t2 + ll_i2,
        5,
        n,
        {"a_TACSTD2_on_CLDN4": a2, "b_immune_on_TACSTD2": b2, "var_C": vc, "var_T_resid": ve_t2, "var_I_resid": ve_i2},
        sigma2,
        Xc,
    )

    # M3: T ⊥ C, I on both. Conditional MLE is OLS on the observed predictors.
    bt, ve_i3, ll_i3 = _reg(I, np.column_stack([T, C]))
    bt_t, bt_c = float(bt[0]), float(bt[1])
    sigma3 = np.array(
        [
            [vt, 0.0, bt_t * vt],
            [0.0, vc, bt_c * vc],
            [bt_t * vt, bt_c * vc, bt_t * bt_t * vt + bt_c * bt_c * vc + ve_i3],
        ]
    )
    m3 = _pack(
        "M3_independent",
        ll_t + ll_c_marg + ll_i3,
        5,
        n,
        {
            "b_immune_on_TACSTD2": bt_t,
            "b_immune_on_CLDN4": bt_c,
            "var_T": vt,
            "var_C": vc,
            "var_I_resid": ve_i3,
        },
        sigma3,
        Xc,
    )

    # Saturated: T → C and I on T and C. Same likelihood as C → T and I on both.
    a_s, ve_c_s, ll_c_s = _reg(C, T[:, None])
    b_s, ve_i_s, ll_i_s = _reg(I, np.column_stack([T, C]))
    a_s = float(a_s[0])
    c_prime, b_med = float(b_s[0]), float(b_s[1])
    # Rebuild covariance from the recursive factorization so it matches S.
    # C = a T + e_c, I = c' T + b C + e_i
    # cov(T,C) = a vt
    # cov(T,I) = c' vt + b a vt
    # cov(C,I) = a c' vt + b var(C)
    var_c_s = a_s * a_s * vt + ve_c_s
    cov_ti = (c_prime + b_med * a_s) * vt
    cov_ci = a_s * c_prime * vt + b_med * var_c_s
    var_i = _var(I)
    sigma_s = np.array(
        [
            [vt, a_s * vt, cov_ti],
            [a_s * vt, var_c_s, cov_ci],
            [cov_ti, cov_ci, var_i],
        ]
    )
    sat = _pack(
        "saturated",
        ll_t + ll_c_s + ll_i_s,
        6,
        n,
        {
            "a_CLDN4_on_TACSTD2": a_s,
            "direct_immune_on_TACSTD2": c_prime,
            "b_immune_on_CLDN4": b_med,
            "var_T": vt,
            "var_C_resid": ve_c_s,
            "var_I_resid": ve_i_s,
        },
        sigma_s,
        Xc,
    )
    out = {"M1": m1, "M2": m2, "M3": m3, "saturated": sat}
    for key, model in out.items():
        if key == "saturated":
            model["lr_vs_saturated"] = 0.0
            model["df_vs_saturated"] = 0
            model["p_vs_saturated"] = 1.0
        else:
            lr = 2.0 * (sat["ll"] - model["ll"])
            lr = float(max(0.0, lr))
            model["lr_vs_saturated"] = lr
            model["df_vs_saturated"] = 1
            model["p_vs_saturated"] = float(stats.chi2.sf(lr, 1))
        model["delta_bic_vs_best_restricted"] = None  # filled by the caller if needed
    _mark_delta_bic(out)
    return out


def _mark_delta_bic(models: dict[str, dict]) -> None:
    best = min(models[k]["bic"] for k in ("M1", "M2", "M3"))
    for key in models:
        models[key]["delta_bic_vs_best_restricted"] = float(models[key]["bic"] - best)


def srmr(Xc: np.ndarray, sigma: np.ndarray) -> float:
    """Standardized root mean square residual of the unique covariance entries."""
    n = Xc.shape[0]
    sample = (Xc.T @ Xc) / n
    d = np.sqrt(np.diag(sample))
    d_hat = np.sqrt(np.diag(sigma))
    # Correlation residuals, including the diagonal (should be ~0 when variances are free).
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = sample / np.outer(d, d)
        corr_hat = sigma / np.outer(d_hat, d_hat)
    resid = corr - corr_hat
    iu = np.triu_indices(3)
    return float(np.sqrt(np.mean(resid[iu] ** 2)))


def mediation(X: np.ndarray, ordering: str) -> dict:
    """Product-of-coefficients mediation with the direct path included.

    ordering 'T_to_C': mediator CLDN4, exposure TACSTD2
    ordering 'C_to_T': mediator TACSTD2, exposure CLDN4
    Proportion = indirect / total. It is not bounded to [0, 1] when the
    direct and indirect paths have opposite signs.
    """
    Xc = center(X)
    T, C, I = Xc[:, 0], Xc[:, 1], Xc[:, 2]
    if ordering == "T_to_C":
        exposure, mediator = T, C
        label = "TACSTD2_to_CLDN4_to_immune"
    elif ordering == "C_to_T":
        exposure, mediator = C, T
        label = "CLDN4_to_TACSTD2_to_immune"
    else:
        raise ValueError(ordering)
    a, ve_m, _ = _reg(mediator, exposure[:, None])
    both, ve_i, _ = _reg(I, np.column_stack([exposure, mediator]))
    total_beta, _, _ = _reg(I, exposure[:, None])
    a = float(a[0])
    direct = float(both[0])
    b = float(both[1])
    indirect = a * b
    total = float(total_beta[0])
    # Algebraic identity for centered OLS: direct + indirect == total.
    pct = indirect / total if abs(total) > 1e-12 else np.nan
    return {
        "ordering": label,
        "a": a,
        "b": b,
        "direct": direct,
        "indirect": indirect,
        "total": total,
        "proportion": float(pct) if np.isfinite(pct) else np.nan,
        "percent": float(100.0 * pct) if np.isfinite(pct) else np.nan,
        "var_mediator_resid": ve_m,
        "var_outcome_resid": ve_i,
    }


def two_stage_residual(X: np.ndarray) -> dict:
    """Two residual regressions.

    Conditional-independence check (the graph discriminator):
      residualize immune on the hypothesized mediator, then correlate the
      residual with the other gene.

    Effect decomposition (residual inclusion):
      stage 1 residualizes the mediator on the exposure;
      stage 2 regresses immune on the exposure and that residual.
      The exposure coefficient equals the total effect. The residual
      coefficient equals the mediator coefficient in the outcome regression.
    """
    Xc = center(X)
    T, C, I = Xc[:, 0], Xc[:, 1], Xc[:, 2]
    n = T.size

    def partial_r(y, x, control):
        _, _, _ = _reg(y, control[:, None])  # touch for symmetry; recompute explicitly
        by, *_ = np.linalg.lstsq(np.column_stack([np.ones(n), control]), y, rcond=None)
        bx, *_ = np.linalg.lstsq(np.column_stack([np.ones(n), control]), x, rcond=None)
        ry = y - np.column_stack([np.ones(n), control]) @ by
        rx = x - np.column_stack([np.ones(n), control]) @ bx
        if np.std(rx) < 1e-12 or np.std(ry) < 1e-12:
            return np.nan, np.nan, np.nan
        r, p = stats.pearsonr(rx, ry)
        return float(r), float(p), float(np.dot(ry, ry) / n)

    r_t_given_c, p_t_given_c, _ = partial_r(I, T, C)
    r_c_given_t, p_c_given_t, _ = partial_r(I, C, T)
    r_tc, p_tc = stats.pearsonr(T, C)

    def inclusion(exposure, mediator):
        a, *_ = np.linalg.lstsq(exposure[:, None], mediator, rcond=None)
        resid = mediator - exposure * float(a[0])
        beta, *_ = np.linalg.lstsq(np.column_stack([exposure, resid]), I, rcond=None)
        return float(a[0]), float(beta[0]), float(beta[1])

    a_tc, total_t, b_c = inclusion(T, C)
    a_ct, total_c, b_t = inclusion(C, T)
    return {
        "partial_r_TACSTD2_immune_given_CLDN4": r_t_given_c,
        "partial_p_TACSTD2_immune_given_CLDN4": p_t_given_c,
        "partial_r_CLDN4_immune_given_TACSTD2": r_c_given_t,
        "partial_p_CLDN4_immune_given_TACSTD2": p_c_given_t,
        "pearson_TACSTD2_CLDN4": float(r_tc),
        "pearson_p_TACSTD2_CLDN4": float(p_tc),
        "inclusion_T_then_C_a": a_tc,
        "inclusion_T_then_C_total": total_t,
        "inclusion_T_then_C_b": b_c,
        "inclusion_C_then_T_a": a_ct,
        "inclusion_C_then_T_total": total_c,
        "inclusion_C_then_T_b": b_t,
    }


def information_criteria_table(models: dict[str, dict]) -> list[dict]:
    rows = []
    for key in ("M1", "M2", "M3", "saturated"):
        m = models[key]
        rows.append(
            {
                "model": m["model"],
                "ll": m["ll"],
                "k": m["k"],
                "n": m["n"],
                "aic": m["aic"],
                "bic": m["bic"],
                "aicc": m["aicc"],
                "delta_bic_vs_best_restricted": m["delta_bic_vs_best_restricted"],
                "lr_vs_saturated": m["lr_vs_saturated"],
                "df_vs_saturated": m["df_vs_saturated"],
                "p_vs_saturated": m["p_vs_saturated"],
                "ll_minus_mvn": m["ll"] - m["ll_mvn"],
            }
        )
    return rows
