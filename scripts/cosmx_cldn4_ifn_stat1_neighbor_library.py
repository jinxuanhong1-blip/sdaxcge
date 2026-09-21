#!/usr/bin/env python3
"""Does the positive CLDN4–neighbor IFN ρ survive library size?

Same tumor gate and 50 µm tumor-neighbor graph as
cosmx_cldn4_ifn_stat1_spatial.py. Residualizes CLDN4 and neighbor IFN
on the index cell's log n_counts and on the mean log n_counts of the
tumor neighbors.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "cosmx_ifn", ROOT / "scripts" / "cosmx_cldn4_ifn_stat1_spatial.py"
)
base = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(base)


def main() -> None:
    meta = base.load_matrix()
    sections, _audit = base.prepare_sections(meta)
    del meta
    rows = []
    for sec in sections:
        graph = sec["graphs"][base.PRIMARY_RADIUS_UM]
        elig = graph["elig"]
        neigh_ifn = base.neighbor_means(sec["core"], graph["indices"], graph["indptr"])
        neigh_logn = base.neighbor_means(sec["logn"], graph["indices"], graph["indptr"])
        cldn = sec["cldn4"][elig]
        own = sec["core"][elig]
        logn = sec["logn"][elig]
        own_lib = logn.reshape(-1, 1)
        both_lib = np.column_stack([logn, neigh_logn])
        rows.append(
            {
                "sample": sec["sample"],
                "patient": sec["patient"],
                "neigh": base.spearman(cldn, neigh_ifn),
                "neigh_ownlib": base.residual_spearman(cldn, neigh_ifn, own_lib),
                "neigh_bothlib": base.residual_spearman(cldn, neigh_ifn, both_lib),
                "partial": base.partial_neighbor_rho(cldn, own, neigh_ifn, None),
                "partial_ownlib": base.partial_neighbor_rho(cldn, own, neigh_ifn, own_lib),
                "partial_bothlib": base.partial_neighbor_rho(cldn, own, neigh_ifn, both_lib),
                "cldn_vs_neigh_logn": base.spearman(cldn, neigh_logn),
            }
        )
        print(sec["sample"], rows[-1]["neigh_bothlib"], flush=True)
    df = pd.DataFrame(rows)
    out = ROOT / "results" / "cosmx_cldn4_ifn_stat1" / "tables" / "neighbor_library.tsv"
    df.to_csv(out, sep="\t", index=False)
    pat = df.groupby("patient").mean(numeric_only=True)
    print(pat.mean().round(3).to_string())


if __name__ == "__main__":
    main()
