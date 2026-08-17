"""scCODA-style composition tests (frequentist fallback).

scCODA (Büttner / Ostner 2021) is a Bayesian Dirichlet-multinomial with a
spike-and-slab on reference-cell-type ALR. The tensorflow / rpy2 pin is not
available here, so the engine is:

- ALR linear model + patient-level permutation p (primary)
- CLR / fraction Wilcoxon as descriptive companions
- Dirichlet-multinomial two-group LRT as a diagnostic (library-size inflated)

A spike-and-slab is not faked. Patient is the unit. Cells are library size.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy import optimize, special, stats


def clr_matrix(counts: np.ndarray, pc: float = 0.5) -> np.ndarray:
    x = np.asarray(counts, dtype=float) + pc
    x = x / x.sum(axis=1, keepdims=True)
    g = np.exp(np.mean(np.log(x), axis=1, keepdims=True))
    return np.log(x / g)


def fractions(counts: np.ndarray) -> np.ndarray:
    x = np.asarray(counts, dtype=float)
    n = x.sum(axis=1, keepdims=True)
    n[n == 0] = np.nan
    return x / n


def alr_coordinate(counts: np.ndarray, target_i: int, ref_i: int, pc: float = 0.5) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    return np.log((counts[:, target_i] + pc) / (counts[:, ref_i] + pc))


def pick_reference(names: list[str], counts: np.ndarray, prefer: tuple[str, ...] = ("Stromal", "Myeloid", "Other", "other", "stromal", "myeloid")) -> str:
    props = counts / np.clip(counts.sum(1, keepdims=True), 1, None)
    mean_p = props.mean(0)
    sd_p = props.std(0, ddof=1) if counts.shape[0] > 1 else np.ones(counts.shape[1])
    for name in prefer:
        if name in names and mean_p[names.index(name)] >= 0.02:
            return name
    eligible = np.where(mean_p >= 0.02)[0]
    if eligible.size == 0:
        eligible = np.arange(len(names))
    # do not pick a primary test compartment as reference if avoidable
    blocked = {i for i, n in enumerate(names) if n in {"TNK", "T_NK", "B", "TLS", "T", "NK"}}
    cand = [i for i in eligible if i not in blocked] or list(eligible)
    return names[int(min(cand, key=lambda i: sd_p[i]))]


def permutation_p_alr(
    counts: np.ndarray,
    x: np.ndarray,
    ref_i: int,
    target_i: int,
    nperm: int = 499,
    seed: int = 1,
) -> tuple[float, str]:
    """Patient-level permutation p for one ALR slope. Does not scale with n_cells."""
    y = alr_coordinate(counts, target_i, ref_i)
    mask = np.isfinite(y) & np.isfinite(x)
    xx, yy = x[mask], y[mask]
    if mask.sum() < 4:
        return float("nan"), "too_few"
    obs = stats.linregress(xx, yy).slope
    uniq = np.unique(xx)
    if set(uniq).issubset({0.0, 1.0}):
        n = len(xx)
        k = int((xx == 1).sum())
        n_comb = int(special.comb(n, k))
        if 1 < n_comb <= 3000:
            extreme = 0
            idx = np.arange(n)
            for comb in combinations(idx, k):
                xp = np.zeros(n)
                xp[list(comb)] = 1.0
                sl = stats.linregress(xp, yy).slope
                if abs(sl) >= abs(obs) - 1e-15:
                    extreme += 1
            return float(extreme / n_comb), f"exact_enumeration nC={n_comb}"
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(nperm):
        sl = stats.linregress(rng.permutation(xx), yy).slope
        if abs(sl) >= abs(obs) - 1e-15:
            extreme += 1
    return float((1 + extreme) / (1 + nperm)), f"monte_carlo nperm={nperm}"


def alr_slope(counts: np.ndarray, x: np.ndarray, ref_i: int, target_i: int) -> dict:
    y = alr_coordinate(counts, target_i, ref_i)
    mask = np.isfinite(y) & np.isfinite(x)
    xx, yy = x[mask], y[mask]
    out = {"effect": float("nan"), "se": float("nan"), "p_ols": float("nan"), "n": int(mask.sum())}
    if mask.sum() < 4:
        return out
    lr = stats.linregress(xx, yy)
    out.update({"effect": float(lr.slope), "se": float(lr.stderr), "p_ols": float(lr.pvalue)})
    return out


def exact_mwu_p(a: np.ndarray, b: np.ndarray) -> tuple[float, str]:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return float("nan"), "too_few"
    if n1 + n2 > 16:
        return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue), "asymptotic_mwu"
    pooled = np.concatenate([a, b])
    obs = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    null = n1 * n2 / 2.0
    extreme = 0
    total = 0
    for idx in combinations(range(n1 + n2), n1):
        mask = np.zeros(n1 + n2, dtype=bool)
        mask[list(idx)] = True
        u = stats.mannwhitneyu(pooled[mask], pooled[~mask], alternative="two-sided").statistic
        total += 1
        if abs(u - null) + 1e-12 >= abs(obs - null):
            extreme += 1
    return extreme / total, "exact_enumeration"


def fraction_mwu(high: np.ndarray, low: np.ndarray) -> dict:
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    high = high[np.isfinite(high)]
    low = low[np.isfinite(low)]
    p, method = exact_mwu_p(high, low)
    return {
        "n_high": int(len(high)),
        "n_low": int(len(low)),
        "median_high": float(np.median(high)) if len(high) else float("nan"),
        "median_low": float(np.median(low)) if len(low) else float("nan"),
        "delta_median_high_minus_low": float(np.median(high) - np.median(low)) if len(high) and len(low) else float("nan"),
        "p": p,
        "p_method": method,
    }


def spearman(x, y) -> dict:
    x = np.asarray(list(x), dtype=float)
    y = np.asarray(list(y), dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    out = {"n": int(len(x)), "rho": float("nan"), "p": float("nan")}
    if len(x) < 4:
        return out
    rho, p = stats.spearmanr(x, y)
    out["rho"] = float(rho)
    out["p"] = float(p)
    return out


def dm_loglik(counts: np.ndarray, alpha: np.ndarray) -> float:
    counts = np.asarray(counts, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    n = counts.sum(axis=1)
    a0 = alpha.sum(axis=1) if alpha.ndim == 2 else float(alpha.sum())
    if alpha.ndim == 1:
        ll = special.gammaln(a0) - special.gammaln(a0 + n)
        ll += special.gammaln(counts + alpha).sum(axis=1) - special.gammaln(alpha).sum()
        return float(np.sum(ll))
    ll = special.gammaln(a0) - special.gammaln(a0 + n)
    ll += special.gammaln(alpha + counts).sum(axis=1) - special.gammaln(alpha).sum(axis=1)
    return float(np.sum(ll))


def fit_dm_shared(counts: np.ndarray) -> tuple[np.ndarray, float]:
    counts = np.asarray(counts, dtype=float)
    p = counts / np.maximum(counts.sum(axis=1, keepdims=True), 1.0)
    mu = np.clip(p.mean(axis=0), 1e-6, 1.0)
    start = np.log(mu * max(float(counts.sum() / max(len(counts), 1) / 10.0), 1.0))

    def nll(log_a):
        a = np.exp(log_a)
        if np.any(~np.isfinite(a)) or np.any(a <= 0):
            return 1e12
        return -dm_loglik(counts, a)

    res = optimize.minimize(nll, start, method="L-BFGS-B")
    alpha = np.exp(res.x)
    return alpha, dm_loglik(counts, alpha)


def dm_two_group(counts: np.ndarray, group: np.ndarray, n_perm: int = 499, seed: int = 1) -> dict:
    """Shared-alpha vs group-specific-alpha LRT. Diagnostic only."""
    counts = np.asarray(counts, dtype=float)
    group = np.asarray(group)
    labs = np.unique(group)
    if len(labs) != 2:
        return {"note": "need exactly two groups", "p_perm": float("nan")}
    g0, g1 = labs
    a0 = group == g0
    a1 = group == g1
    if a0.sum() < 2 or a1.sum() < 2:
        return {"note": "need >=2 per group", "p_perm": float("nan"), "n0": int(a0.sum()), "n1": int(a1.sum())}
    _, ll0 = fit_dm_shared(counts)
    _, ll_a = fit_dm_shared(counts[a0])
    _, ll_b = fit_dm_shared(counts[a1])
    stat = 2.0 * (ll_a + ll_b - ll0)
    p_chi = float(stats.chi2.sf(max(stat, 0.0), counts.shape[1]))
    rng = np.random.default_rng(seed)
    extreme = 0
    n_ok = 0
    labels = group.copy()
    for _ in range(n_perm):
        rng.shuffle(labels)
        try:
            _, p_ll0 = fit_dm_shared(counts)
            _, p_lla = fit_dm_shared(counts[labels == g0])
            _, p_llb = fit_dm_shared(counts[labels == g1])
            pstat = 2.0 * (p_lla + p_llb - p_ll0)
        except Exception:
            continue
        n_ok += 1
        if pstat + 1e-9 >= stat:
            extreme += 1
    p_perm = (extreme + 1) / (n_ok + 1) if n_ok else float("nan")
    return {
        "n0": int(a0.sum()),
        "n1": int(a1.sum()),
        "lrt_stat": float(stat),
        "p_chi2": p_chi,
        "p_perm": float(p_perm),
        "n_perm": int(n_ok),
        "note": "diagnostic_dm_lrt",
    }


def try_sccoda() -> dict:
    try:
        import anndata  # noqa: F401
        from sccoda.util import comp_ana  # noqa: F401
    except Exception as exc:  # noqa: BLE001
        return {"status": "not_available", "reason": str(exc)}
    return {"status": "import_ok_not_run"}


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(n, dtype=float)
    tmp[order] = q
    out[ok] = tmp
    return out
