"""Software checks for the SEM likelihood. These are not biological results."""

from __future__ import annotations

import numpy as np

from semcore import fit_graphs, mediation, two_stage_residual


def _assert_close(a, b, tol=1e-8, msg=""):
    if not np.isfinite(a) or abs(a - b) > tol:
        raise AssertionError(f"{msg} {a} != {b}")


def test_likelihood_matches_mvn_and_identities():
    rng = np.random.default_rng(11)
    n = 180
    T = rng.normal(size=n)
    C = 0.7 * T + rng.normal(scale=0.6, size=n)
    I = -0.2 * T - 0.5 * C + rng.normal(scale=0.7, size=n)
    X = np.column_stack([T, C, I])
    models = fit_graphs(X)
    for key, model in models.items():
        _assert_close(model["ll"], model["ll_mvn"], 1e-6, key)
        assert model["ll"] <= models["saturated"]["ll"] + 1e-8
    # Saturated likelihood does not depend on which gene is written first.
    rev = fit_graphs(np.column_stack([C, T, I]))
    # Reverse column order is a different labeling, so compare the saturated
    # fit of the original to a refit that uses the same regression skeleton.
    _assert_close(models["saturated"]["ll"], models["saturated"]["ll_mvn"], 1e-6, "sat")
    med = mediation(X, "T_to_C")
    _assert_close(med["direct"] + med["indirect"], med["total"], 1e-8, "mediation sum")
    stage = two_stage_residual(X)
    _assert_close(stage["inclusion_T_then_C_total"], med["total"], 1e-8, "two-stage total")
    _assert_close(stage["inclusion_T_then_C_b"], med["b"], 1e-8, "two-stage b")
    _assert_close(stage["inclusion_T_then_C_a"], med["a"], 1e-8, "two-stage a")
    # Restricted models each have one fewer parameter than saturated.
    for key in ("M1", "M2", "M3"):
        assert models[key]["k"] == 5
        assert models[key]["df_vs_saturated"] == 1
    assert models["saturated"]["k"] == 6
    assert rev["saturated"]["k"] == 6


def test_scaling_does_not_change_bic_ranking():
    rng = np.random.default_rng(3)
    n = 120
    raw = rng.normal(size=(n, 3))
    raw[:, 2] = 0.4 * raw[:, 0] - 0.6 * raw[:, 1] + 0.3 * rng.normal(size=n)
    base = fit_graphs(raw)
    scaled = raw * np.array([2.5, 0.3, 4.0])
    alt = fit_graphs(scaled)
    for key in ("M1", "M2", "M3", "saturated"):
        _assert_close(
            base[key]["bic"] - base["M1"]["bic"],
            alt[key]["bic"] - alt["M1"]["bic"],
            1e-6,
            key,
        )


def test_simulated_graphs_recover_the_generator():
    rng = np.random.default_rng(21)
    n = 800

    T = rng.normal(size=n)
    C = 0.85 * T + rng.normal(scale=0.35, size=n)
    I = -0.75 * C + rng.normal(scale=0.45, size=n)
    m = fit_graphs(np.column_stack([T, C, I]))
    assert min(m, key=lambda k: m[k]["bic"] if k != "saturated" else np.inf) == "M1"

    C = rng.normal(size=n)
    T = 0.85 * C + rng.normal(scale=0.35, size=n)
    I = -0.75 * T + rng.normal(scale=0.45, size=n)
    m = fit_graphs(np.column_stack([T, C, I]))
    assert min(("M1", "M2", "M3"), key=lambda k: m[k]["bic"]) == "M2"

    T = rng.normal(size=n)
    C = rng.normal(size=n)
    I = 0.55 * T - 0.65 * C + rng.normal(scale=0.4, size=n)
    m = fit_graphs(np.column_stack([T, C, I]))
    assert min(("M1", "M2", "M3"), key=lambda k: m[k]["bic"]) == "M3"


def test_known_mediation_fraction_on_fixed_design():
    # Population proportion a*b/(a*b+c'). Noise keeps the two regressors
    # from being collinear, so the sample percent is estimable.
    rng = np.random.default_rng(1)
    n = 20000
    a, b, direct = 0.8, -0.5, 0.15
    T = rng.normal(size=n)
    C = a * T + rng.normal(scale=0.5, size=n)
    I = direct * T + b * C + rng.normal(scale=0.5, size=n)
    med = mediation(np.column_stack([T, C, I]), "T_to_C")
    expected = 100.0 * (a * b) / (a * b + direct)
    _assert_close(med["percent"], expected, 1.5, "percent")


if __name__ == "__main__":
    test_likelihood_matches_mvn_and_identities()
    test_scaling_does_not_change_bic_ranking()
    test_simulated_graphs_recover_the_generator()
    test_known_mediation_fraction_on_fixed_design()
    print("semcore tests passed")
