#!/usr/bin/env python3
"""Lock a few honest numbers so the writeup cannot drift."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from lib import meta_fisher_z, spearman_rho_p

ROOT = Path(__file__).resolve().parents[1]


def approx(a, b, tol=1e-6):
    if abs(float(a) - float(b)) > tol:
        raise SystemExit(f"mismatch: {a} vs {b}")


def main() -> None:
    gse = pd.read_csv(ROOT / "harvested" / "gse253013_per_patient.tsv", sep="\t")
    sub = gse[(gse["tissue"] == "Tumor") & (gse["eligible_malig"])]
    rho, p, n = spearman_rho_p(sub["TACSTD2_mean_log1p"], sub["tnk_fraction"])
    approx(n, 9)
    approx(rho, -0.7166666666666667, 1e-9)
    approx(p, 0.02981803569584528, 1e-9)

    effects = pd.read_csv(ROOT / "results" / "cohort_effects.tsv", sep="\t")
    luad = effects[(effects["gene"] == "TACSTD2") & (effects["primary_meta"] == 1) & (effects["histology"] == "LUAD")]
    meta = meta_fisher_z(luad.to_dict("records"), min_n=6)
    approx(meta["k"], 5)
    approx(meta["n_total"], 72)
    approx(meta["rho_re"], -0.208467, 1e-5)
    print("selftest ok", {"k": meta["k"], "n": meta["n_total"], "rho": meta["rho_re"]})


if __name__ == "__main__":
    main()
