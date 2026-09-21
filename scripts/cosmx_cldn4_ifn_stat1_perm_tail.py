#!/usr/bin/env python3
"""One-sided spatial permutation p for a positive neighbor ρ (descriptive).

The anti-correlation tail is already in summary.json. This records the
other tail and the null SD so the positive neighbor ρ can be compared
with the within-FOV shuffle null.
"""

from __future__ import annotations

import json
from pathlib import Path

import importlib.util
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "cosmx_ifn", ROOT / "scripts" / "cosmx_cldn4_ifn_stat1_spatial.py"
)
base = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(base)

N = 199
RNG = np.random.default_rng(20260922)


def main() -> None:
    meta = base.load_matrix()
    sections, _ = base.prepare_sections(meta)
    del meta
    keys = ["neigh", "partial"]
    patients, observed = [], {k: [] for k in keys}
    for sec in sections:
        g = sec["graphs"][base.PRIMARY_RADIUS_UM]
        nm = base.neighbor_means(sec["core"], g["indices"], g["indptr"])
        met = base.metrics_for_score(sec, "core", base.PRIMARY_RADIUS_UM, nm)
        patients.append(sec["patient"])
        for k in keys:
            observed[k].append(met[k])
    obs = {k: base.patient_mean_from_sections(patients, observed[k]) for k in keys}
    nulls = {k: np.empty(N) for k in keys}
    for b in range(N):
        pats, vals = [], {k: [] for k in keys}
        for sec in sections:
            g = sec["graphs"][base.PRIMARY_RADIUS_UM]
            shuffled = base.shuffle_within_fov(sec["core"], sec["fov"], RNG)
            nm = base.neighbor_means(shuffled, g["indices"], g["indptr"])
            elig = g["elig"]
            cldn = sec["cldn4"][elig]
            own = sec["core"][elig]
            vals["neigh"].append(base.spearman(cldn, nm))
            vals["partial"].append(base.partial_neighbor_rho(cldn, own, nm, None))
            pats.append(sec["patient"])
        for k in keys:
            nulls[k][b] = base.patient_mean_from_sections(pats, vals[k])
        if (b + 1) % 50 == 0:
            print(b + 1, flush=True)
    out = {"n_perm": N}
    for k in keys:
        p_greater = (1 + int(np.sum(nulls[k] >= obs[k]))) / (N + 1)
        p_less = (1 + int(np.sum(nulls[k] <= obs[k]))) / (N + 1)
        out[k] = {
            "observed_patient_mean_rho": obs[k],
            "null_mean": float(nulls[k].mean()),
            "null_sd": float(nulls[k].std()),
            "perm_p_onesided_greater": float(p_greater),
            "perm_p_onesided_less": float(p_less),
        }
        print(k, out[k], flush=True)
    path = ROOT / "results" / "cosmx_cldn4_ifn_stat1" / "tables" / "perm_positive_tail.json"
    path.write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
