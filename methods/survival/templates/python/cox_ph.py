#!/usr/bin/env python3
"""Template: Cox PH for a continuous ICI biomarker (lifelines).

Primary analysis is the continuous marker (HR per 1 SD). Dichotomisation is a
sensitivity analysis, not the primary model. See playbook.md sections 3–5.

Required columns: time, event (1 = event, 0 = censored), marker.
Optional: a small prespecified covariate list. Refuse the model when events
per variable fall below 10 (Peduzzi / Vittinghoff rule of thumb).
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="CSV with time, event, marker")
    ap.add_argument("--time-col", default="time")
    ap.add_argument("--event-col", default="event")
    ap.add_argument("--marker-col", default="marker")
    ap.add_argument("--covariates", default="", help="comma-separated extra columns")
    ap.add_argument("--min-epv", type=float, default=10.0)
    args = ap.parse_args()

    covs = [c.strip() for c in args.covariates.split(",") if c.strip()]
    cols = [args.time_col, args.event_col, args.marker_col] + covs
    df = pd.read_csv(args.input)[cols].dropna()
    n_events = int(df[args.event_col].sum())
    n_var = 1 + len(covs)
    epv = n_events / n_var
    if epv < args.min_epv:
        sys.exit(
            f"refuse: events per variable = {epv:.1f} "
            f"({n_events} events / {n_var} variables) < {args.min_epv}. "
            "Fit a univariable model or collect more events."
        )

    df = df.copy()
    df[args.marker_col] = (
        df[args.marker_col] - df[args.marker_col].mean()
    ) / df[args.marker_col].std(ddof=1)

    cph = CoxPHFitter()
    cph.fit(df, duration_col=args.time_col, event_col=args.event_col)
    print(cph.summary.to_string())
    print(f"\nHarrell C (apparent) = {cph.concordance_index_:.3f}")
    print(f"events per variable = {epv:.1f}")

    zph = proportional_hazard_test(cph, df, time_transform=["km", "rank"])
    print("\nSchoenfeld PH tests:")
    print(zph.summary.to_string())
    if (zph.summary["p"] < 0.05).any():
        print(
            "PH is questionable. Report RMST or a time-varying coefficient; "
            "do not quote a single HR as if it were constant."
        )


if __name__ == "__main__":
    main()
