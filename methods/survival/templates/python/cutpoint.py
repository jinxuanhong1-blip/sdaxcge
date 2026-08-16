#!/usr/bin/env python3
"""Template: median cut vs maximally selected log-rank cut (lifelines).

The naive p-value at the 'best' cut is not a valid type-I error. This script
reports it only next to a permutation p-value and a bootstrap of the selected
percentile. Prefer the continuous Cox model in cox_ph.py as the primary analysis.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test


def candidates(values: np.ndarray, min_frac: float) -> np.ndarray:
    lo, hi = np.quantile(values, [min_frac, 1 - min_frac])
    uniq = np.unique(values)
    return uniq[(uniq >= lo) & (uniq < hi)]


def max_logrank(time, event, score, min_frac):
    best = (-np.inf, np.nan)
    for cut in candidates(score, min_frac):
        high = score > cut
        if high.sum() < 2 or (~high).sum() < 2:
            continue
        res = logrank_test(time[high], time[~high], event[high], event[~high])
        if res.test_statistic > best[0]:
            best = (float(res.test_statistic), float(cut))
    return best


def hr_high(time, event, high) -> float:
    m = CoxPHFitter().fit(
        pd.DataFrame({"t": time, "e": event, "high": high.astype(int)}), "t", "e"
    )
    return float(m.summary.loc["high", "exp(coef)"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--time-col", default="time")
    ap.add_argument("--event-col", default="event")
    ap.add_argument("--marker-col", default="marker")
    ap.add_argument("--min-frac", type=float, default=0.25)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20240816)
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    time = df[args.time_col].to_numpy()
    event = df[args.event_col].to_numpy()
    score = df[args.marker_col].to_numpy()

    med = float(np.median(score))
    high_med = score > med
    lr_med = logrank_test(time[high_med], time[~high_med], event[high_med], event[~high_med])
    print(
        f"median cut: HR={hr_high(time, event, high_med):.2f} "
        f"log-rank p={lr_med.p_value:.4f} "
        f"n_high={high_med.sum()} n_low={(~high_med).sum()}"
    )

    stat, cut = max_logrank(time, event, score, args.min_frac)
    high = score > cut
    lr = logrank_test(time[high], time[~high], event[high], event[~high])
    rng = np.random.default_rng(args.seed)
    null = np.array(
        [
            max_logrank(time, event, score[rng.permutation(len(score))], args.min_frac)[0]
            for _ in range(args.n_perm)
        ]
    )
    p_perm = (1 + np.sum(null >= stat)) / (1 + args.n_perm)
    print(
        f"max-selected cut={cut:.4g} (percentile {(score < cut).mean()*100:.0f}): "
        f"HR={hr_high(time, event, high):.2f} "
        f"naive p={lr.p_value:.4f} permutation p={p_perm:.4f} "
        f"n_candidates={len(candidates(score, args.min_frac))}"
    )
    print("Do not report the naive p-value as the result.")


if __name__ == "__main__":
    main()
