#!/usr/bin/env python3
"""Fit scCODA-style ALR / DM tests: T/NK and B vs malignant CLDN4-high vs low.

Primary family (pre-specified): one TNK and one B test per primary slice.
Recovery = ALR effect < 0 and BH q < 0.10 on patient-permutation p.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lib.compositional import (  # noqa: E402
    alr_slope,
    bh_fdr,
    dm_two_group,
    fraction_mwu,
    fractions,
    permutation_p_alr,
    pick_reference,
    spearman,
    try_sccoda,
)

OUT = ROOT / "results"
MIN_N = 6
MIN_ARM = 2
PRIMARY_SLICES = {
    "GSE207422_post_A3mal",
    "GSE241934_IIT",
    "GSE241934_RWC",
    "GSE291670",
    "GSE205335",
    "GSE253013_tumor",
}
# GSE207422_drmref is sensitivity (different malignant definition)


def wide_for_slice(long: pd.DataFrame, table: pd.DataFrame, sl: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = table[table["slice"] == sl].copy()
    cl = long[long["slice"] == sl].copy()
    wide = cl.pivot_table(index="sample_id", columns="celltype", values="n", aggfunc="sum", fill_value=0)
    wide = wide.loc[:, wide.sum(axis=0) > 0]
    return wide, meta.set_index("sample_id")


def fit_slice(sl: str, wide: pd.DataFrame, meta: pd.DataFrame) -> tuple[list[dict], dict]:
    elig = meta[meta["eligible"].astype(bool)].copy()
    use = [i for i in elig.index if i in wide.index]
    wide = wide.loc[use]
    elig = elig.loc[use]
    rows = []
    dm_info = {"slice": sl, "note": "skipped"}
    if len(wide) < MIN_N:
        return rows, {"slice": sl, "note": f"n={len(wide)} < {MIN_N}"}
    x = elig["cldn4_high"].to_numpy(dtype=float)
    if (x == 1).sum() < MIN_ARM or (x == 0).sum() < MIN_ARM:
        return rows, {"slice": sl, "note": "unbalanced split"}
    counts = wide.to_numpy(dtype=float)
    names = list(wide.columns)
    ref = pick_reference(names, counts)
    ref_i = names.index(ref)
    props = fractions(counts)
    is_primary_slice = sl in PRIMARY_SLICES

    targets = []
    if "TNK" in names:
        targets.append(("TNK", names.index("TNK")))
    if "B" in names:
        targets.append(("B", names.index("B")))

    for compartment, ti in targets:
        slp = alr_slope(counts, x, ref_i, ti)
        p_perm, p_method = permutation_p_alr(counts, x, ref_i, ti)
        hi = props[x == 1, ti]
        lo = props[x == 0, ti]
        mwu = fraction_mwu(hi, lo)
        sp = spearman(elig["cldn4_score"], props[:, ti])
        rows.append(
            {
                "slice": sl,
                "cohort": elig["cohort"].iloc[0],
                "annotation": elig["annotation"].iloc[0],
                "compartment": compartment,
                "reference": ref,
                "n": int(len(wide)),
                "n_high": int((x == 1).sum()),
                "n_low": int((x == 0).sum()),
                "alr_effect": slp["effect"],
                "alr_se": slp["se"],
                "alr_p_ols": slp["p_ols"],
                "alr_p_perm": p_perm,
                "alr_p_method": p_method,
                "frac_median_high": mwu["median_high"],
                "frac_median_low": mwu["median_low"],
                "frac_delta_high_minus_low": mwu["delta_median_high_minus_low"],
                "frac_mwu_p": mwu["p"],
                "frac_mwu_method": mwu["p_method"],
                "spearman_rho": sp["rho"],
                "spearman_p": sp["p"],
                "spearman_n": sp["n"],
                "primary_family": bool(is_primary_slice),
                "direction_tnk_or_b_down": bool(np.isfinite(slp["effect"]) and slp["effect"] < 0),
            }
        )

    dm_info = dm_two_group(counts, x, n_perm=199, seed=1)
    dm_info["slice"] = sl
    dm_info["n"] = int(len(wide))
    dm_info["reference"] = ref
    dm_info["celltypes"] = names
    return rows, dm_info


def main() -> None:
    table = pd.read_csv(OUT / "composition_table.tsv", sep="\t")
    long = pd.read_csv(OUT / "composition_long.tsv", sep="\t")
    tests = []
    dms = []
    for sl in table["slice"].drop_duplicates():
        wide, meta = wide_for_slice(long, table, sl)
        rows, dm = fit_slice(sl, wide, meta)
        tests.extend(rows)
        dms.append(dm)

    # pre-specified ICI pool: cohort-centered CLDN4, primary slices only, 4-part counts
    pool_slices = ["GSE207422_post_A3mal", "GSE241934_IIT", "GSE241934_RWC", "GSE291670", "GSE205335", "GSE253013_tumor"]
    pool_rows = []
    for sl in pool_slices:
        wide, meta = wide_for_slice(long, table, sl)
        elig = meta[meta["eligible"].astype(bool)]
        use = [i for i in elig.index if i in wide.index]
        if not use:
            continue
        tmp = elig.loc[use].copy()
        tmp["slice"] = sl
        # 4-part: Epithelial, TNK, B, Other-or-rest
        w = wide.loc[use].copy()
        rest_cols = [c for c in w.columns if c not in {"Epithelial", "TNK", "B"}]
        tmp["Epithelial"] = w["Epithelial"] if "Epithelial" in w.columns else 0
        tmp["TNK"] = w["TNK"]
        tmp["B"] = w["B"] if "B" in w.columns else 0
        tmp["Other"] = w[rest_cols].sum(axis=1) if rest_cols else 0
        tmp["cldn4_cc"] = tmp["cldn4_score"] - tmp["cldn4_score"].median()
        pool_rows.append(tmp.reset_index())
    if pool_rows:
        pool = pd.concat(pool_rows, ignore_index=True)
        # within-pool median of cohort-centered CLDN4
        cut = float(pool["cldn4_cc"].median())
        pool["cldn4_high"] = (pool["cldn4_cc"] >= cut).astype(float)
        w4 = pool[["Epithelial", "TNK", "B", "Other"]]
        counts = w4.to_numpy(dtype=float)
        names = list(w4.columns)
        x = pool["cldn4_high"].to_numpy(dtype=float)
        ref = pick_reference(names, counts)
        ref_i = names.index(ref)
        props = fractions(counts)
        for compartment in ["TNK", "B"]:
            ti = names.index(compartment)
            slp = alr_slope(counts, x, ref_i, ti)
            p_perm, p_method = permutation_p_alr(counts, x, ref_i, ti)
            mwu = fraction_mwu(props[x == 1, ti], props[x == 0, ti])
            sp = spearman(pool["cldn4_cc"], props[:, ti])
            tests.append(
                {
                    "slice": "ICI_pool_cohort_centered",
                    "cohort": "pool",
                    "annotation": "mixed_4part",
                    "compartment": compartment,
                    "reference": ref,
                    "n": int(len(pool)),
                    "n_high": int((x == 1).sum()),
                    "n_low": int((x == 0).sum()),
                    "alr_effect": slp["effect"],
                    "alr_se": slp["se"],
                    "alr_p_ols": slp["p_ols"],
                    "alr_p_perm": p_perm,
                    "alr_p_method": p_method,
                    "frac_median_high": mwu["median_high"],
                    "frac_median_low": mwu["median_low"],
                    "frac_delta_high_minus_low": mwu["delta_median_high_minus_low"],
                    "frac_mwu_p": mwu["p"],
                    "frac_mwu_method": mwu["p_method"],
                    "spearman_rho": sp["rho"],
                    "spearman_p": sp["p"],
                    "spearman_n": sp["n"],
                    "primary_family": False,
                    "direction_tnk_or_b_down": bool(np.isfinite(slp["effect"]) and slp["effect"] < 0),
                }
            )
        dms.append(dm_two_group(counts, x, n_perm=199, seed=1) | {"slice": "ICI_pool_cohort_centered", "n": int(len(pool))})

    testdf = pd.DataFrame(tests)
    prim = testdf["primary_family"].astype(bool)
    testdf["q_bh_primary"] = np.nan
    testdf.loc[prim, "q_bh_primary"] = bh_fdr(testdf.loc[prim, "alr_p_perm"].to_numpy())
    n_prim = int(prim.sum())
    testdf["bonferroni_primary"] = np.nan
    testdf.loc[prim, "bonferroni_primary"] = np.clip(testdf.loc[prim, "alr_p_perm"] * n_prim, 0, 1)
    testdf["recovery"] = "not_primary"
    for i, row in testdf.iterrows():
        if not row["primary_family"]:
            testdf.at[i, "recovery"] = "not_primary"
            continue
        if not np.isfinite(row.get("alr_effect", np.nan)) or not np.isfinite(row.get("q_bh_primary", np.nan)):
            testdf.at[i, "recovery"] = "not_estimable"
        elif row["q_bh_primary"] < 0.10 and row["alr_effect"] < 0:
            testdf.at[i, "recovery"] = "recovered_down"
        elif row["q_bh_primary"] < 0.10 and row["alr_effect"] > 0:
            testdf.at[i, "recovery"] = "opposite_up"
        else:
            testdf.at[i, "recovery"] = "not_recovered"

    testdf.to_csv(OUT / "tests_full.tsv", sep="\t", index=False)
    testdf[testdf["primary_family"]].to_csv(OUT / "tests_primary.tsv", sep="\t", index=False)
    pd.DataFrame(dms).to_csv(OUT / "dm_diagnostic.tsv", sep="\t", index=False)

    recovered = testdf[testdf["recovery"] == "recovered_down"]
    summary = {
        "sccoda": try_sccoda(),
        "n_primary_tests": n_prim,
        "n_recovered_tnk_or_b_down": int(len(recovered)),
        "recovered": recovered[["slice", "compartment", "n", "alr_effect", "alr_p_perm", "q_bh_primary"]].to_dict(orient="records"),
        "primary_pvalue": "patient-level permutation of the ALR slope (does not scale with n_cells)",
        "honest_n_note": "n is patients. Cells are compositional library size, not the replicate unit.",
        "exposure": "within-slice median split of malignant / author-epithelial CLDN4. TACSTD2 is never a gate.",
        "q_cut": 0.10,
        "excluded_non_ici": ["GSE131907 (treatment-naive atlas)", "GSE267108 / GSE274595 leftover (not ICI; malignant CLDN4 empty)"],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(testdf[["slice", "compartment", "n", "n_high", "n_low", "alr_effect", "alr_p_perm", "q_bh_primary", "recovery"]].to_string(index=False))


if __name__ == "__main__":
    main()
