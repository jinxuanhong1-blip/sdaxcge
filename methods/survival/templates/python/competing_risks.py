#!/usr/bin/env python3
"""Template: Aalen–Johansen CIF vs 1−KM that censors the competing event (lifelines).

Use when the event of interest (progression, irAE, TTR) can be precluded by
death or treatment switch. Requires a cause code: 0 = censored, 1 = event of
interest, 2+ = competing events. Most GEO ICI series do not deposit this field.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd
from lifelines import AalenJohansenFitter, KaplanMeierFitter


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--time-col", default="time")
    ap.add_argument("--cause-col", default="cause", help="0=censored, 1=interest, 2+=competing")
    ap.add_argument("--horizon", type=float, required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.input)[[args.time_col, args.cause_col]].dropna()
    causes = set(df[args.cause_col].astype(int))
    if causes <= {0, 1}:
        sys.exit(
            "no competing event is coded. Do not run a competing-risks model "
            "just because the playbook mentions one."
        )
    time = df[args.time_col]
    cause = df[args.cause_col].astype(int)
    aj = AalenJohansenFitter(calculate_variance=True)
    aj.fit(time, cause, event_of_interest=1)
    cif = float(aj.predict(args.horizon))
    km = KaplanMeierFitter().fit(time, (cause == 1).astype(int))
    naive = float(1 - km.predict(args.horizon))
    print(f"Aalen-Johansen CIF(cause=1, t={args.horizon:g}) = {cif:.3f}")
    print(f"1-KM treating competing events as censoring     = {naive:.3f}")
    print(f"absolute overestimate                           = {naive - cif:.3f}")


if __name__ == "__main__":
    main()
