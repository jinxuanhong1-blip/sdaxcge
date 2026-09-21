"""Patient-level compositional models for concordant-4 cell fractions.

ILR uses a pre-specified sequential binary partition. Dirichlet-multinomial
uses a reference-part softmax mean and a shared precision. Zeros in the ILR
path get a +1/2 cell pseudocount; the Dirichlet-multinomial uses raw counts.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp, softmax


def proportions_with_pseudocount(counts: np.ndarray, pseudo: float = 0.5) -> np.ndarray:
    x = np.asarray(counts, dtype=float) + pseudo
    if np.any(x <= 0):
        raise ValueError("pseudocount must keep every part positive")
    return x / x.sum(axis=1, keepdims=True)


def ilr_from_sbp(proportions: np.ndarray, sbp: np.ndarray) -> np.ndarray:
    """ILR coordinates for a sequential binary partition.

    sbp rows are balances. Entries are +1 (numerator), -1 (denominator), 0 (absent).
    """
    p = np.asarray(proportions, dtype=float)
    codes = np.asarray(sbp, dtype=float)
    if p.ndim != 2 or codes.ndim != 2 or codes.shape[1] != p.shape[1]:
        raise ValueError("proportion and SBP shapes do not match")
    out = np.zeros((p.shape[0], codes.shape[0]), dtype=float)
    logp = np.log(p)
    for b, row in enumerate(codes):
        pos = row > 0
        neg = row < 0
        r = float(pos.sum())
        s = float(neg.sum())
        if r < 1 or s < 1:
            raise ValueError("each balance needs a numerator and a denominator")
        scale = math.sqrt(r * s / (r + s))
        gpos = logp[:, pos].mean(axis=1)
        gneg = logp[:, neg].mean(axis=1)
        out[:, b] = scale * (gpos - gneg)
    return out


def sbp_orthonorm_check(sbp: np.ndarray, tol: float = 1e-8) -> bool:
    """Return True when the SBP contrast matrix is orthonormal."""
    codes = np.asarray(sbp, dtype=float)
    rows = []
    for row in codes:
        pos = row > 0
        neg = row < 0
        r = float(pos.sum())
        s = float(neg.sum())
        scale = math.sqrt(r * s / (r + s))
        v = np.zeros_like(row)
        v[pos] = scale / r
        v[neg] = -scale / s
        rows.append(v)
    psi = np.vstack(rows)
    gram = psi @ psi.T
    return bool(np.allclose(gram, np.eye(len(rows)), atol=tol))


# 4-part order: T, NK, Myeloid, Rest
SBP_TNK_MYE_REST = np.array(
    [
        [1, 1, 1, -1],  # T/NK/myeloid vs rest of the sample
        [1, 1, -1, 0],  # T/NK vs myeloid
        [1, -1, 0, 0],  # T vs NK
    ],
    dtype=float,
)
BALANCE_NAMES_4 = (
    "immune_focus_vs_rest",
    "lymphoid_vs_myeloid",
    "T_vs_NK",
)

# 3-part order: T, NK, Myeloid (renormalized; mix among the three)
SBP_WITHIN_IMMUNE = np.array(
    [
        [1, 1, -1],
        [1, -1, 0],
    ],
    dtype=float,
)
BALANCE_NAMES_3 = (
    "lymphoid_vs_myeloid",
    "T_vs_NK",
)


def within_cohort_z(values: np.ndarray, cohorts: np.ndarray) -> np.ndarray:
    z = np.zeros(len(values), dtype=float)
    for coh in np.unique(cohorts):
        m = cohorts == coh
        x = np.asarray(values[m], dtype=float)
        sd = float(x.std(ddof=1)) if m.sum() > 1 else 1.0
        if not np.isfinite(sd) or sd == 0:
            sd = 1.0
        z[m] = (x - x.mean()) / sd
    return z


def design_matrix(cohorts: np.ndarray, covariate: np.ndarray, ref: str) -> tuple[np.ndarray, list[str]]:
    levels = [ref] + [c for c in sorted(set(cohorts.tolist())) if c != ref]
    names = ["intercept"] + [f"cohort_{c}" for c in levels[1:]] + ["cldn4"]
    x = np.zeros((len(cohorts), len(names)), dtype=float)
    x[:, 0] = 1.0
    for j, coh in enumerate(levels[1:], start=1):
        x[:, j] = (cohorts == coh).astype(float)
    x[:, -1] = np.asarray(covariate, dtype=float)
    return x, names


def ols_hc1(y: np.ndarray, x: np.ndarray) -> dict:
    """OLS with HC1 robust standard errors. y is 1d."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    n, p = x.shape
    xtx = x.T @ x
    try:
        xtx_inv = np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        xtx_inv = np.linalg.pinv(xtx)
    beta = xtx_inv @ (x.T @ y)
    resid = y - x @ beta
    meat = x.T @ ((resid**2)[:, None] * x)
    # HC1
    scale = n / max(n - p, 1)
    cov = scale * (xtx_inv @ meat @ xtx_inv)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    dof = max(n - p, 1)
    return {"beta": beta, "se": se, "resid": resid, "dof": dof, "n": n, "p": p}


def student_p(beta: float, se: float, dof: int) -> float:
    if not np.isfinite(se) or se <= 0:
        return float("nan")
    from scipy.stats import t

    return float(2 * t.sf(abs(beta / se), dof))


def spearman_rho(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    from scipy.stats import spearmanr

    if len(x) < 4:
        return float("nan"), float("nan")
    rho, p = spearmanr(x, y)
    return float(rho), float(p)


def dl_fisher_z(rhos: np.ndarray, ns: np.ndarray) -> dict:
    """DerSimonian–Laird pool of Spearman rhos on the Fisher-z scale."""
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    ok = np.isfinite(rhos) & np.isfinite(ns) & (ns > 3)
    rhos, ns = rhos[ok], ns[ok]
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    var = 1.0 / (ns - 3.0)
    w = 1.0 / var
    zbar = np.sum(w * z) / np.sum(w)
    q = float(np.sum(w * (z - zbar) ** 2))
    k = int(len(rhos))
    dfree = k - 1
    cdenom = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (q - dfree) / cdenom) if dfree > 0 and cdenom > 0 else 0.0
    wstar = 1.0 / (var + tau2)
    zre = float(np.sum(wstar * z) / np.sum(wstar))
    se = math.sqrt(1.0 / float(np.sum(wstar)))
    from scipy.stats import norm

    p = float(2 * norm.sf(abs(zre / se)))
    i2 = max(0.0, (q - dfree) / q) if q > 0 else 0.0
    ci = np.tanh(zre + np.array([-1.0, 1.0]) * 1.96 * se)
    return {
        "rho": float(np.tanh(zre)),
        "p": p,
        "I2": float(i2),
        "ci_lo": float(ci[0]),
        "ci_hi": float(ci[1]),
        "k": k,
        "N": int(ns.sum()),
        "tau2": float(tau2),
    }


def permute_within_cohort(covariate: np.ndarray, cohorts: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = np.array(covariate, copy=True)
    for coh in np.unique(cohorts):
        m = np.flatnonzero(cohorts == coh)
        out[m] = rng.permutation(out[m])
    return out


@dataclass
class DMFit:
    beta: np.ndarray  # (K-1, p) reference-coded; last part is the reference
    log_phi: float
    success: bool
    nll: float
    part_names: list[str]
    coef_names: list[str]


def _dm_nll(theta: np.ndarray, counts: np.ndarray, x: np.ndarray, k: int) -> float:
    p = x.shape[1]
    beta = theta[: (k - 1) * p].reshape(k - 1, p)
    log_phi = float(theta[-1])
    # keep precision in a range that the gamma functions can handle
    log_phi = float(np.clip(log_phi, math.log(0.2), math.log(5000.0)))
    phi = math.exp(log_phi)
    eta = x @ beta.T  # (n, K-1)
    # reference part has eta 0
    eta_full = np.concatenate([eta, np.zeros((eta.shape[0], 1))], axis=1)
    mu = softmax(eta_full, axis=1)
    alpha = phi * mu
    n = counts.sum(axis=1)
    # Dirichlet-multinomial log-likelihood
    ll = (
        gammaln(n + 1.0)
        - gammaln(counts + 1.0).sum(axis=1)
        + gammaln(phi)
        - gammaln(n + phi)
        + (gammaln(counts + alpha) - gammaln(alpha)).sum(axis=1)
    )
    if not np.isfinite(ll).all():
        return 1e12
    return float(-ll.sum())


def fit_dirichlet_multinomial(
    counts: np.ndarray,
    x: np.ndarray,
    part_names: list[str],
    coef_names: list[str],
    n_starts: int = 3,
    seed: int = 1,
) -> DMFit:
    counts = np.asarray(counts, dtype=float)
    x = np.asarray(x, dtype=float)
    k = counts.shape[1]
    p = x.shape[1]
    if len(part_names) != k:
        raise ValueError("part_names length")
    rng = np.random.default_rng(seed)
    # moment start: CLR of mean composition onto X, reference = last part
    props = proportions_with_pseudocount(counts, 0.5)
    logp = np.log(props)
    clr = logp - logp.mean(axis=1, keepdims=True)
    # reference-coded eta ≈ clr_k - clr_ref
    eta0 = clr[:, :-1] - clr[:, -1:]
    beta0 = np.linalg.lstsq(x, eta0, rcond=None)[0].T  # (K-1, p)
    n = counts.sum(axis=1)
    # rough precision from the average part
    mean_mu = props.mean(axis=0)
    var = props.var(axis=0).mean()
    m = max(float(mean_mu.mean() * (1 - mean_mu.mean())), 1e-6)
    phi0 = max(1.0, (m / max(var, 1e-8) - 1.0))
    starts = [np.concatenate([beta0.ravel(), [math.log(phi0)]])]
    starts.append(np.concatenate([np.zeros((k - 1) * p), [math.log(20.0)]]))
    for _ in range(max(0, n_starts - 2)):
        noise = rng.normal(scale=0.05, size=(k - 1) * p)
        starts.append(np.concatenate([beta0.ravel() + noise, [math.log(phi0)]]))

    best = None
    for theta in starts:
        opt = minimize(
            _dm_nll,
            theta,
            args=(counts, x, k),
            method="L-BFGS-B",
            options={"maxiter": 400, "ftol": 1e-10},
        )
        if best is None or opt.fun < best.fun:
            best = opt
    theta = best.x
    beta = theta[: (k - 1) * p].reshape(k - 1, p)
    return DMFit(
        beta=beta,
        log_phi=float(theta[-1]),
        success=bool(best.success),
        nll=float(best.fun),
        part_names=list(part_names),
        coef_names=list(coef_names),
    )


def dm_mean_proportions(fit: DMFit, x: np.ndarray) -> np.ndarray:
    eta = x @ fit.beta.T
    eta_full = np.concatenate([eta, np.zeros((eta.shape[0], 1))], axis=1)
    return softmax(eta_full, axis=1)


def finite_difference_proportions(fit: DMFit, x: np.ndarray, coef_index: int, delta: float = 1.0) -> np.ndarray:
    """Mean predicted proportion at covariate +delta/2 minus -delta/2, other columns fixed."""
    hi = np.array(x, copy=True)
    lo = np.array(x, copy=True)
    hi[:, coef_index] = hi[:, coef_index] + delta / 2.0
    lo[:, coef_index] = lo[:, coef_index] - delta / 2.0
    return dm_mean_proportions(fit, hi).mean(axis=0) - dm_mean_proportions(fit, lo).mean(axis=0)
