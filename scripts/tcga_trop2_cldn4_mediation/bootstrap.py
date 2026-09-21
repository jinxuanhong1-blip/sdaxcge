"""Patient bootstrap of partial correlations and mediation proportions.

One resample ranks every column once, then reuses those ranks for every
covariate set. Seeds stay local to the call so cohorts can run in parallel.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from stats import _mediation_ranked, _partial_pearson, coef_given

# Primary analysis columns. The driver builds this matrix and nothing else.
COLUMNS = [
    "TACSTD2",
    "CLDN4",
    "CLDN7",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "CD8_score",
    "CD3_score",
]

KERATIN = ["KRT8", "KRT18", "KRT19"]
RHO_MODELS = {
    "keratin": KERATIN,
    "keratin_CLDN4": KERATIN + ["CLDN4"],
    "keratin_CLDN7": KERATIN + ["CLDN7"],
    "keratin_EPCAM": KERATIN + ["EPCAM"],
    "keratin_CLDN7_EPCAM": KERATIN + ["CLDN7", "EPCAM"],
    "keratin_CLDN7_EPCAM_CLDN4": KERATIN + ["CLDN7", "EPCAM", "CLDN4"],
}
OUTCOMES = ["CD8_score", "CD3_score"]
MEDIATORS = ["CLDN4", "CLDN7", "EPCAM"]
PREDICTOR = "TACSTD2"
MODEL_K = {name: len(cols) for name, cols in RHO_MODELS.items()}

# LUSC sensitivity: squamous keratins replace KRT8/18/19. Not in the pool.
SQUAMOUS = ["KRT5", "KRT6A", "KRT14"]
SQUAMOUS_MODELS = {
    "squamous": SQUAMOUS,
    "squamous_CLDN4": SQUAMOUS + ["CLDN4"],
    "squamous_CLDN7": SQUAMOUS + ["CLDN7"],
    "squamous_EPCAM": SQUAMOUS + ["EPCAM"],
}


def _rank_sample(sample: np.ndarray) -> np.ndarray:
    ranked = np.empty(sample.shape, dtype=float)
    for j in range(sample.shape[1]):
        ranked[:, j] = stats.rankdata(sample[:, j])
    return ranked


def _metrics_ranked(ranked, index, models, outcomes, mediators, predictor, base_covariates, incremental_pair):
    x = ranked[:, index[predictor]]
    rhos = {}
    for outcome in outcomes:
        y = ranked[:, index[outcome]]
        rhos[outcome] = {}
        for model, cov_names in models.items():
            Z = ranked[:, [index[name] for name in cov_names]]
            rhos[outcome][model] = _partial_pearson(x, y, Z)
    pms = {}
    Z_base = ranked[:, [index[name] for name in base_covariates]]
    for outcome in outcomes:
        y = ranked[:, index[outcome]]
        pms[outcome] = {}
        for mediator in mediators:
            m = ranked[:, index[mediator]]
            pms[outcome][mediator] = _mediation_ranked(x, m, y, Z_base)["proportion"]
    incremental = {}
    if incremental_pair is not None:
        ctrl_names, full_names = incremental_pair
        for outcome in outcomes:
            y = ranked[:, index[outcome]]
            c = coef_given(y, x, Z_base)
            c_ctrl = coef_given(y, x, ranked[:, [index[name] for name in ctrl_names]])
            c_full = coef_given(y, x, ranked[:, [index[name] for name in full_names]])
            incremental[outcome] = (c_ctrl - c_full) / c if abs(c) > 1e-12 else np.nan
    return rhos, pms, incremental


def bootstrap_design(
    data: dict[str, np.ndarray],
    models: dict[str, list[str]],
    outcomes: list[str],
    mediators: list[str],
    predictor: str,
    base_covariates: list[str],
    n_boot: int,
    seed: int,
    incremental_pair: tuple[list[str], list[str]] | None = None,
) -> dict:
    """Bootstrap partial rho and rank-OLS mediation for one patient table.

    `data` values are aligned 1d arrays with no missing entries. The returned
    arrays have length n_boot. Point estimates are not mixed into the draws.
    """
    if n_boot < 1:
        raise ValueError("n_boot must be positive")
    names = [predictor, *outcomes, *mediators, *base_covariates]
    for covs in models.values():
        names.extend(covs)
    if incremental_pair is not None:
        names.extend(incremental_pair[0])
        names.extend(incremental_pair[1])
    ordered = list(dict.fromkeys(names))
    missing = [name for name in ordered if name not in data]
    if missing:
        raise KeyError(f"bootstrap data missing {missing}")
    matrix = np.column_stack([np.asarray(data[name], dtype=float) for name in ordered])
    if not np.isfinite(matrix).all():
        raise ValueError("bootstrap data contains non-finite values")
    index = {name: i for i, name in enumerate(ordered)}
    n = int(matrix.shape[0])
    rng = np.random.default_rng(seed)
    rho = {outcome: {model: np.empty(n_boot) for model in models} for outcome in outcomes}
    pm = {outcome: {mediator: np.empty(n_boot) for mediator in mediators} for outcome in outcomes}
    incremental = {outcome: np.empty(n_boot) for outcome in outcomes} if incremental_pair else {}
    indices = rng.integers(0, n, size=(n_boot, n))
    for b in range(n_boot):
        ranked = _rank_sample(matrix[indices[b]])
        rhos, pms, inc = _metrics_ranked(
            ranked, index, models, outcomes, mediators, predictor, base_covariates, incremental_pair
        )
        for outcome in outcomes:
            for model in models:
                rho[outcome][model][b] = rhos[outcome][model]
            for mediator in mediators:
                pm[outcome][mediator][b] = pms[outcome][mediator]
            if incremental_pair is not None:
                incremental[outcome][b] = inc[outcome]
    full = _metrics_ranked(
        _rank_sample(matrix), index, models, outcomes, mediators, predictor, base_covariates, incremental_pair
    )
    return {
        "n": n,
        "n_boot": int(n_boot),
        "seed": int(seed),
        "rho": rho,
        "pm": pm,
        "incremental": incremental,
        "full_rho": full[0],
        "full_pm": full[1],
        "full_incremental": full[2],
    }


def bootstrap_primary(data: dict[str, np.ndarray], n_boot: int, seed: int) -> dict:
    """Primary keratin design: CLDN4 mediator, CLDN7 and EPCAM controls."""
    return bootstrap_design(
        data=data,
        models=RHO_MODELS,
        outcomes=OUTCOMES,
        mediators=MEDIATORS,
        predictor=PREDICTOR,
        base_covariates=KERATIN,
        n_boot=n_boot,
        seed=seed,
        incremental_pair=(KERATIN + ["CLDN7", "EPCAM"], KERATIN + ["CLDN7", "EPCAM", "CLDN4"]),
    )
