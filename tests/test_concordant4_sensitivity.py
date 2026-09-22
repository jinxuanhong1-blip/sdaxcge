"""Locked concordant-4 ρ is reproduced and is not replaced by a larger joint panel."""

import pandas as pd

from scripts.concordant4_outcome_sensitivity import (
    LOCKED_RHO,
    cliffs_delta,
    prepare,
    run_grid,
    UNITS,
)


def test_cliffs_delta_sign():
    assert cliffs_delta([1, 2], [3, 4]) < 0
    assert cliffs_delta([3, 4], [1, 2]) > 0


def test_locked_panel_matches_published_rho():
    units = prepare(pd.read_csv(UNITS, sep="\t"))
    grid = run_grid(units)
    locked = grid.loc[grid["panel_id"] == "frac_tnk|pct_gt0|none|none|all"].iloc[0]
    assert abs(float(locked["pooled_rho"]) - LOCKED_RHO) < 5e-4
    assert float(locked["I2"]) == 0.0
    assert int(locked["n_total"]) == 65
    assert float(locked["cliff"]) < -0.7
    frac = grid[(grid["outcome"] == "frac_tnk") & (grid["consistent"]) & (grid["I2"] <= 1e-8)]
    assert float(frac["pooled_rho"].min()) == float(locked["pooled_rho"])
