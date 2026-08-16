#!/usr/bin/env python3
"""Template: time-dependent (cumulative/dynamic) AUC with IPCW (scikit-survival).

Pass a RISK score (higher = earlier event). For a protective marker, use the
Cox linear predictor, not the raw expression. Apparent AUC on the same data
used to fit the Cox model is optimistic; nest this inside CV or apply it to
a locked external cohort.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from sksurv.metrics import cumulative_dynamic_auc
from sksurv.util import Surv


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--time-col", default="time")
    ap.add_argument("--event-col", default="event")
    ap.add_argument("--marker-col", default="marker")
    ap.add_argument("--times", default="3,6,12", help="evaluation times in the same unit as time")
    args = ap.parse_args()

    df = pd.read_csv(args.input)[[args.time_col, args.event_col, args.marker_col]].dropna()
    y = Surv.from_arrays(event=df[args.event_col].astype(bool), time=df[args.time_col])
    cph = CoxPHFitter().fit(df, args.time_col, args.event_col)
    risk = cph.predict_partial_hazard(df).to_numpy()
    last_event = float(df.loc[df[args.event_col] == 1, args.time_col].max())
    times = [float(t) for t in args.times.split(",") if float(t) < last_event]
    if not times:
        sys.exit(f"no evaluation time is earlier than the last event ({last_event})")
    auc, mean_auc = cumulative_dynamic_auc(y, y, risk, np.array(times))
    print(f"Cox log-HR for {args.marker_col} = {cph.summary.loc[args.marker_col, 'coef']:.3f}")
    for t, a in zip(times, np.atleast_1d(auc)):
        n_evt = int(((df[args.time_col] <= t) & (df[args.event_col] == 1)).sum())
        n_risk = int((df[args.time_col] >= t).sum())
        print(f"AUC({t:g}) = {a:.3f}  events_by={n_evt}  still_at_risk={n_risk}")
    print(f"mean AUC over {times} = {mean_auc:.3f} (apparent; not externally validated)")


if __name__ == "__main__":
    main()
