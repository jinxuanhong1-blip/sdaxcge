"""Dirichlet-multinomial and ALR composition models.

scCODA (Büttner, Ostner et al. 2021) is a Bayesian DM with a spike-and-slab
on a reference-cell-type ALR. The tensorflow/tfp pin is brittle in this
environment, so the *default* engine here is a frequentist DM-GLM with the
same estimand: log(π_k / π_ref) = β0_k + β1_k x.

A spike-and-slab is *not* faked. If `sccoda` imports, it is run as an
optional Bayesian companion and written to a separate table. Otherwise the
run is DM-GLM + ALR-LM + fraction Wilcoxon/Spearman, and scCODA is marked
`not_available`.

Small-n honesty: with n≈9–15 the DM Hessian is often ill-conditioned.
Ridge (default 1.0 on slopes) and a shared precision are used; failures
are returned as `not_estimable` rather than p=0.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, special, stats


def _softmax(eta: np.ndarray) -> np.ndarray:
    z = eta - np.max(eta, axis=-1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=-1, keepdims=True)


def dm_loglik(counts: np.ndarray, alpha: np.ndarray) -> float:
    """Dirichlet-multinomial log-likelihood (multinomial constant omitted)."""
    counts = np.asarray(counts, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    n = counts.sum(axis=1)
    a0 = alpha.sum(axis=1)
    ll = special.gammaln(a0) - special.gammaln(a0 + n)
    ll += special.gammaln(alpha + counts).sum(axis=1) - special.gammaln(alpha).sum(axis=1)
    return float(np.sum(ll))


@dataclass
class DMResult:
    celltypes: list[str]
    reference: str
    covariate: str
    intercept: np.ndarray
    slope: np.ndarray
    precision: float
    se_slope: np.ndarray
    z_slope: np.ndarray
    p_slope: np.ndarray
    loglik: float
    converged: bool
    n_samples: int
    n_celltypes: int
    note: str = ""


def fit_dm_glm(
    counts: np.ndarray,
    x: np.ndarray,
    celltypes: list[str],
    reference: str | None = None,
    ridge: float = 1.0,
    maxiter: int = 400,
) -> DMResult:
    """One-covariate DM-GLM in ALR coordinates.

    counts : (n, K) non-negative integers
    x      : (n,) covariate (0/1 or continuous, already centered if desired)
    reference cell type has slope/intercept fixed at 0 (compositional identifiability).
    """
    counts = np.asarray(counts, dtype=float)
    x = np.asarray(x, dtype=float)
    n, k = counts.shape
    names = list(celltypes)
    if reference is None:
        # automatic: least variable relative abundance among types with mean > 2%
        props = counts / np.clip(counts.sum(1, keepdims=True), 1, None)
        mean_p = props.mean(0)
        sd_p = props.std(0, ddof=1) if n > 1 else np.ones(k)
        eligible = np.where(mean_p >= 0.02)[0]
        if eligible.size == 0:
            eligible = np.arange(k)
        ref_i = int(eligible[np.argmin(sd_p[eligible])])
    else:
        if reference not in names:
            raise ValueError(f"reference {reference} not in {names}")
        ref_i = names.index(reference)
    others = [i for i in range(k) if i != ref_i]
    p = len(others)  # intercepts + slopes + log_precision

    def unpack(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
        b0 = np.zeros(k)
        b1 = np.zeros(k)
        b0[others] = theta[:p]
        b1[others] = theta[p : 2 * p]
        log_phi = theta[-1]
        phi = float(np.exp(np.clip(log_phi, -4, 8)))
        return b0, b1, phi

    def nll(theta: np.ndarray, freeze_slope: int | None = None) -> float:
        b0, b1, phi = unpack(theta)
        if freeze_slope is not None:
            b1[freeze_slope] = 0.0
        eta = b0[None, :] + np.outer(x, b1)
        pi = _softmax(eta)
        alpha = np.clip(phi * pi, 1e-6, None)
        ll = dm_loglik(counts, alpha)
        pen = 0.5 * ridge * (np.sum(b1**2) + 0.1 * np.sum(b0**2))
        return -ll + pen

    theta0 = np.zeros(2 * p + 1)
    # start precision near the mean library size / 4
    theta0[-1] = np.log(max(counts.sum(1).mean() / 4.0, 2.0))
    opt = optimize.minimize(
        nll,
        theta0,
        method="L-BFGS-B",
        options={"maxiter": maxiter, "ftol": 1e-8},
    )
    b0, b1, phi = unpack(opt.x)
    se = np.full(k, np.nan)
    z = np.full(k, np.nan)
    pv = np.full(k, np.nan)
    note = "p_from_LRT"
    if not opt.success:
        note = f"not_converged:{opt.message};p_from_LRT"
    # Likelihood-ratio test per non-reference slope. Do not use L-BFGS hess_inv
    # (it produced 1e-40-scale p-values on n=9 in an earlier draft).
    ll_full = -nll(opt.x) + 0.5 * ridge * (np.sum(b1**2) + 0.1 * np.sum(b0**2))
    for j, ct_i in enumerate(others):
        theta_r = opt.x.copy()
        theta_r[p + j] = 0.0
        opt_r = optimize.minimize(
            lambda th, idx=ct_i: nll(th, freeze_slope=idx),
            theta_r,
            method="L-BFGS-B",
            options={"maxiter": maxiter, "ftol": 1e-8},
        )
        b0_r, b1_r, _ = unpack(opt_r.x)
        b1_r[ct_i] = 0.0
        ll_r = -nll(opt_r.x, freeze_slope=ct_i) + 0.5 * ridge * (np.sum(b1_r**2) + 0.1 * np.sum(b0_r**2))
        lrt = max(0.0, 2.0 * (ll_full - ll_r))
        pv[ct_i] = float(stats.chi2.sf(lrt, 1))
        se[ct_i] = np.nan
        z[ct_i] = np.sign(b1[ct_i]) * np.sqrt(lrt) if lrt > 0 else 0.0
    return DMResult(
        celltypes=names,
        reference=names[ref_i],
        covariate="x",
        intercept=b0,
        slope=b1,
        precision=phi,
        se_slope=se,
        z_slope=z,
        p_slope=pv,
        loglik=-float(opt.fun) if np.isfinite(opt.fun) else np.nan,
        converged=bool(opt.success),
        n_samples=int(n),
        n_celltypes=int(k),
        note=note,
    )


def alr_lm(
    counts: np.ndarray,
    x: np.ndarray,
    celltypes: list[str],
    reference: str,
    hc: bool = True,
) -> list[dict]:
    """Per-type ALR linear model: log((y_k+0.5)/(y_ref+0.5)) ~ x.

    Uses OLS + HC3 if statsmodels is available, else scipy linregress.
    """
    counts = np.asarray(counts, dtype=float)
    x = np.asarray(x, dtype=float)
    names = list(celltypes)
    ref_i = names.index(reference)
    yref = counts[:, ref_i] + 0.5
    rows = []
    try:
        import statsmodels.api as sm
    except Exception:
        sm = None
    for i, name in enumerate(names):
        if i == ref_i:
            rows.append(
                {
                    "celltype": name,
                    "method": "alr_lm",
                    "effect": 0.0,
                    "se": 0.0,
                    "p_value": 1.0,
                    "n": int(len(x)),
                    "note": "reference",
                }
            )
            continue
        y = np.log((counts[:, i] + 0.5) / yref)
        mask = np.isfinite(y) & np.isfinite(x)
        if mask.sum() < 5:
            rows.append(
                {
                    "celltype": name,
                    "method": "alr_lm",
                    "effect": np.nan,
                    "se": np.nan,
                    "p_value": np.nan,
                    "n": int(mask.sum()),
                    "note": "too_few_samples",
                }
            )
            continue
        xx, yy = x[mask], y[mask]
        if sm is not None:
            X = sm.add_constant(xx)
            fit = sm.OLS(yy, X).fit(cov_type="HC3") if hc else sm.OLS(yy, X).fit()
            rows.append(
                {
                    "celltype": name,
                    "method": "alr_lm",
                    "effect": float(fit.params[1]),
                    "se": float(fit.bse[1]),
                    "p_value": float(fit.pvalues[1]),
                    "n": int(mask.sum()),
                    "note": "HC3" if hc else "OLS",
                }
            )
        else:
            lr = stats.linregress(xx, yy)
            rows.append(
                {
                    "celltype": name,
                    "method": "alr_lm",
                    "effect": float(lr.slope),
                    "se": float(lr.stderr),
                    "p_value": float(lr.pvalue),
                    "n": int(mask.sum()),
                    "note": "linregress",
                }
            )
    return rows


def fraction_tests(counts: np.ndarray, x: np.ndarray, celltypes: list[str], binary: bool) -> list[dict]:
    """Naive (non-compositional) tests on fractions. Reported as comparators only."""
    props = counts / np.clip(counts.sum(1, keepdims=True), 1, None)
    rows = []
    for i, name in enumerate(celltypes):
        y = props[:, i]
        mask = np.isfinite(y) & np.isfinite(x)
        xx, yy = x[mask], y[mask]
        if mask.sum() < 4:
            rows.append(
                {
                    "celltype": name,
                    "method": "fraction_naive",
                    "effect": np.nan,
                    "se": np.nan,
                    "p_value": np.nan,
                    "n": int(mask.sum()),
                    "note": "too_few_samples",
                }
            )
            continue
        if binary and set(np.unique(xx)).issubset({0.0, 1.0}) and (xx == 0).sum() >= 2 and (xx == 1).sum() >= 2:
            a, b = yy[xx == 1], yy[xx == 0]
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append(
                {
                    "celltype": name,
                    "method": "fraction_mwu",
                    "effect": float(np.median(a) - np.median(b)),
                    "se": np.nan,
                    "p_value": float(p),
                    "n": int(mask.sum()),
                    "note": f"U={u:.1f}; median_high-low",
                }
            )
        else:
            rho, p = stats.spearmanr(xx, yy)
            rows.append(
                {
                    "celltype": name,
                    "method": "fraction_spearman",
                    "effect": float(rho),
                    "se": np.nan,
                    "p_value": float(p),
                    "n": int(mask.sum()),
                    "note": "spearman_on_fraction",
                }
            )
    return rows


def permutation_p_alr(
    counts: np.ndarray,
    x: np.ndarray,
    ref_i: int,
    target_i: int,
    nperm: int = 499,
    seed: int = 1,
) -> float:
    """Patient-level permutation p for one ALR slope. Does not scale with n_cells."""
    y = np.log((counts[:, target_i] + 0.5) / (counts[:, ref_i] + 0.5))
    mask = np.isfinite(y) & np.isfinite(x)
    xx, yy = x[mask], y[mask]
    if mask.sum() < 5:
        return np.nan
    obs = stats.linregress(xx, yy).slope
    rng = np.random.default_rng(seed)
    # exact enumeration when the binary split is small
    uniq = np.unique(xx)
    if set(uniq).issubset({0.0, 1.0}):
        from itertools import combinations

        n = len(xx)
        k = int((xx == 1).sum())
        n_comb = int(special.comb(n, k))
        if 1 < n_comb <= 2000:
            idx = np.arange(n)
            extreme = 0
            for comb in combinations(idx, k):
                xp = np.zeros(n)
                xp[list(comb)] = 1.0
                sl = stats.linregress(xp, yy).slope
                if abs(sl) >= abs(obs) - 1e-15:
                    extreme += 1
            return float(extreme / n_comb)
    extreme = 0
    for _ in range(nperm):
        xp = rng.permutation(xx)
        sl = stats.linregress(xp, yy).slope
        if abs(sl) >= abs(obs) - 1e-15:
            extreme += 1
    return float((1 + extreme) / (1 + nperm))


def try_sccoda(counts, covariates, celltypes, formula: str, reference: str) -> dict:
    """Optional scCODA HMC. Returns a status dict; never invents credible effects."""
    try:
        import anndata as ad
        from sccoda.util import comp_ana as mod
    except Exception as exc:  # noqa: BLE001
        return {"status": "not_available", "reason": str(exc)}
    try:
        adata = ad.AnnData(X=np.asarray(counts, dtype=float))
        adata.obs = covariates.copy()
        adata.var_names = list(celltypes)
        model = mod.CompositionalAnalysis(adata, formula=formula, reference_cell_type=reference)
        res = model.sample_hmc()
        return {"status": "ok", "summary": str(res.summary())}
    except Exception as exc:  # noqa: BLE001
        return {"status": "failed", "reason": str(exc)}
