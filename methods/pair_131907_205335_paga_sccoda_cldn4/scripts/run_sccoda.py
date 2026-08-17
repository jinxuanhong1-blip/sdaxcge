#!/usr/bin/env python3
"""scCODA-style composition: CLDN4 Q4 vs Q1 on the winning pair.

Table-only. No Harmony, no UMI, no GSE148071, no 7-cohort pool.
Patient (GSE205335) / sample (GSE131907) is the unit.

Primary family (6 tests): T/NK and B × {GSE131907, GSE205335, merge}
on within-cohort author-malignant CLDN4 %pos Q4 vs Q1 — the same cut
that already differs in PR #320 (n=23, r=−0.705). Median split is
sensitivity only.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
    spearman,
    try_sccoda,
)

TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
MIN_MALIGNANT = 20


def assign_quartiles(values: pd.Series) -> pd.Series:
    """Match PR #320: rank then qcut into Q1–Q4."""
    ranks = values.rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return pd.Series(index=values.index, dtype="object")
    return qs.astype(str)


def load_units(given: Path) -> pd.DataFrame:
    a = pd.read_csv(given / "GSE131907_samples.tsv", sep="\t")
    a = a.loc[
        a["origin"].isin(TUMOR_ORIGINS) & (pd.to_numeric(a["n_malignant"], errors="coerce") >= MIN_MALIGNANT)
    ].copy()
    a["cohort"] = "GSE131907"
    a["unit_id"] = "GSE131907:" + a["sample"].astype(str)
    a["unit_type"] = "sample"
    a["cldn4_pct"] = pd.to_numeric(a["mal_CLDN4_pct"], errors="coerce")
    a["cldn4_mean"] = pd.to_numeric(a["mal_CLDN4_mean"], errors="coerce")
    a["n_TNK"] = pd.to_numeric(a["n_tnk"], errors="coerce")
    a["n_B"] = pd.to_numeric(a["n_b"], errors="coerce")
    a["n_cells"] = pd.to_numeric(a["n_cells"], errors="coerce")
    a["frac_TNK"] = pd.to_numeric(a["frac_tnk"], errors="coerce")
    a["frac_B"] = pd.to_numeric(a["frac_b"], errors="coerce")
    a["origin_or_tissue"] = a["origin"].astype(str)
    a["display"] = a["sample"].astype(str)

    b = pd.read_csv(given / "GSE205335_patients.tsv", sep="\t")
    b["cohort"] = "GSE205335"
    b["unit_id"] = "GSE205335:" + b["patient"].astype(str)
    b["unit_type"] = "patient"
    b["cldn4_pct"] = pd.to_numeric(b["mal_CLDN4_pct_pos"], errors="coerce")
    b["cldn4_mean"] = pd.to_numeric(b["mal_CLDN4_mean"], errors="coerce")
    b["n_TNK"] = pd.to_numeric(b["n_tnk"], errors="coerce")
    b["n_B"] = pd.to_numeric(b["n_b_plasma"], errors="coerce")
    b["n_cells"] = pd.to_numeric(b["n_cells"], errors="coerce")
    b["frac_TNK"] = pd.to_numeric(b["frac_tnk"], errors="coerce")
    b["frac_B"] = pd.to_numeric(b["frac_b_plasma"], errors="coerce")
    b["origin_or_tissue"] = b["tissue"].astype(str)
    b["display"] = b["patient"].astype(str)
    b["n_malignant"] = pd.to_numeric(b["n_malignant"], errors="coerce")

    keep = [
        "cohort",
        "unit_id",
        "unit_type",
        "display",
        "origin_or_tissue",
        "n_cells",
        "n_malignant",
        "n_TNK",
        "n_B",
        "frac_TNK",
        "frac_B",
        "cldn4_pct",
        "cldn4_mean",
    ]
    frame = pd.concat([a[keep], b[keep]], ignore_index=True)
    frame["n_Other"] = (frame["n_cells"] - frame["n_TNK"] - frame["n_B"]).clip(lower=0)
    frame["frac_Other"] = frame["n_Other"] / frame["n_cells"]
    frame["cldn4_score"] = frame["cldn4_pct"]
    frame["quartile"] = "NA"
    frame["cldn4_q4q1"] = np.nan
    frame["cldn4_median_high"] = np.nan
    for cohort, idx in frame.groupby("cohort").groups.items():
        sub = frame.loc[list(idx), "cldn4_score"]
        qs = assign_quartiles(sub)
        frame.loc[list(idx), "quartile"] = qs
        frame.loc[list(idx), "cldn4_q4q1"] = qs.map({"Q4": 1, "Q1": 0})
        cut = float(sub.median())
        frame.loc[list(idx), "cldn4_median_cut"] = cut
        frame.loc[list(idx), "cldn4_median_high"] = (sub >= cut).astype(int)
    return frame


def _fit_block(comp: pd.DataFrame, x: np.ndarray, slice_name: str, split: str) -> list[dict]:
    counts = comp[["n_TNK", "n_B", "n_Other"]].to_numpy(dtype=float)
    names = ["TNK", "B", "Other"]
    ref_i = names.index("Other")
    props = fractions(counts)
    rows = []
    for compartment, ti in (("TNK", 0), ("B", 1)):
        slp = alr_slope(counts, x, ref_i, ti)
        p_perm, p_method = permutation_p_alr(counts, x, ref_i, ti, nperm=499, seed=13)
        hi = props[x == 1, ti]
        lo = props[x == 0, ti]
        mwu = fraction_mwu(hi, lo)
        sp = spearman(comp["cldn4_score"], props[:, ti])
        rows.append(
            {
                "slice": slice_name,
                "split": split,
                "cohort": ",".join(sorted(comp["cohort"].unique())),
                "compartment": compartment,
                "reference": "Other",
                "n": int(len(comp)),
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
                "primary_family": split == "q4q1",
                "direction_down": bool(np.isfinite(slp["effect"]) and slp["effect"] < 0),
            }
        )
    return rows


def run_tests(comp: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    rows: list[dict] = []
    dm_rows: list[dict] = []
    slices = {
        "GSE131907": comp[comp["cohort"] == "GSE131907"],
        "GSE205335": comp[comp["cohort"] == "GSE205335"],
        "merge_within_cohort": comp,
    }
    for sl, frame in slices.items():
        q = frame[frame["cldn4_q4q1"].isin([0, 1])].copy()
        if len(q) >= 6:
            x = q["cldn4_q4q1"].to_numpy(dtype=float)
            rows.extend(_fit_block(q, x, sl, "q4q1"))
            dm = dm_two_group(q[["n_TNK", "n_B", "n_Other"]].to_numpy(dtype=float), x, n_perm=199, seed=13)
            dm.update({"slice": sl, "split": "q4q1"})
            dm_rows.append(dm)
        med = frame[frame["cldn4_median_high"].isin([0, 1])].copy()
        if len(med) >= 6:
            x = med["cldn4_median_high"].to_numpy(dtype=float)
            rows.extend(_fit_block(med, x, sl, "median"))
            dm = dm_two_group(med[["n_TNK", "n_B", "n_Other"]].to_numpy(dtype=float), x, n_perm=199, seed=13)
            dm.update({"slice": sl, "split": "median"})
            dm_rows.append(dm)
    tests = pd.DataFrame(rows)
    primary = tests["primary_family"].to_numpy(dtype=bool)
    tests["q_bh_primary"] = np.nan
    if primary.any():
        tests.loc[primary, "q_bh_primary"] = bh_fdr(tests.loc[primary, "alr_p_perm"].to_numpy())
    tests["recover_q10"] = (
        tests["primary_family"]
        & tests["direction_down"]
        & (tests["q_bh_primary"] < 0.10)
    )
    return tests, dm_rows


def honest_n(comp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort, sub in comp.groupby("cohort"):
        q = sub[sub["cldn4_q4q1"].isin([0, 1])]
        rows.append(
            {
                "cohort": cohort,
                "unit_type": sub["unit_type"].iloc[0],
                "n_eligible": int(len(sub)),
                "n_q1": int((sub["quartile"] == "Q1").sum()),
                "n_q4": int((sub["quartile"] == "Q4").sum()),
                "n_q4q1_compared": int(len(q)),
                "n_median_high": int((sub["cldn4_median_high"] == 1).sum()),
                "n_median_low": int((sub["cldn4_median_high"] == 0).sum()),
                "cldn4": "author-malignant %pos",
                "note": (
                    "GSE131907 sample is the unit (not patient)"
                    if cohort == "GSE131907"
                    else "GSE205335 B = B+plasma; all 22 patients"
                ),
            }
        )
    q = comp[comp["cldn4_q4q1"].isin([0, 1])]
    rows.append(
        {
            "cohort": "merge_within_cohort",
            "unit_type": "sample+patient",
            "n_eligible": int(len(comp)),
            "n_q1": int((comp["quartile"] == "Q1").sum()),
            "n_q4": int((comp["quartile"] == "Q4").sum()),
            "n_q4q1_compared": int(len(q)),
            "n_median_high": int((comp["cldn4_median_high"] == 1).sum()),
            "n_median_low": int((comp["cldn4_median_high"] == 0).sum()),
            "cldn4": "within-cohort author-malignant %pos",
            "note": "PR #320 Q4 vs Q1 n=23; not a 7-cohort pool",
        }
    )
    return pd.DataFrame(rows)


def fig_fractions(comp: pd.DataFrame, out: Path) -> None:
    q = comp[comp["cldn4_q4q1"].isin([0, 1])].copy()
    q["arm"] = q["cldn4_q4q1"].map({0: "Q1", 1: "Q4"})
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), constrained_layout=True)
    for ax, col, title in (
        (axes[0], "frac_TNK", "T/NK fraction"),
        (axes[1], "frac_B", "B fraction"),
    ):
        for i, cohort in enumerate(["GSE131907", "GSE205335"]):
            sub = q[q["cohort"] == cohort]
            for j, arm in enumerate(["Q1", "Q4"]):
                vals = sub.loc[sub["arm"] == arm, col].to_numpy(dtype=float)
                x = i + (j - 0.5) * 0.28
                ax.scatter(
                    np.full(len(vals), x) + np.linspace(-0.04, 0.04, len(vals)),
                    vals,
                    s=28,
                    alpha=0.85,
                    color="#1f77b4" if arm == "Q1" else "#d62728",
                    zorder=3,
                )
                if len(vals):
                    ax.hlines(np.median(vals), x - 0.08, x + 0.08, color="black", lw=1.4)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["GSE131907\n(sample)", "GSE205335\n(patient)"])
        ax.set_ylabel(title)
        ax.set_title(title + " · CLDN4 Q4 vs Q1")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#1f77b4", label="Q1", markersize=8),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#d62728", label="Q4", markersize=8),
    ]
    axes[1].legend(handles=handles, frameon=False, loc="upper right")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def fig_alr_forest(tests: pd.DataFrame, out: Path) -> None:
    prim = tests[tests["split"] == "q4q1"].copy()
    if prim.empty:
        return
    prim = prim.sort_values(["compartment", "slice"])
    fig, ax = plt.subplots(figsize=(7.4, 3.8), constrained_layout=True)
    y = np.arange(len(prim))
    colors = ["#d62728" if e < 0 else "#1f77b4" for e in prim["alr_effect"]]
    ax.axvline(0, color="0.6", lw=0.8)
    ax.errorbar(
        prim["alr_effect"],
        y,
        xerr=1.96 * prim["alr_se"].to_numpy(dtype=float),
        fmt="none",
        ecolor="0.45",
        elinewidth=1.1,
        capsize=2,
    )
    ax.scatter(prim["alr_effect"], y, c=colors, s=36, zorder=3)
    labels = [f"{r.slice} · {r.compartment}  n={r.n_low}+{r.n_high}" for r in prim.itertuples()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("ALR effect (Q4 − Q1); <0 = compartment down in CLDN4-high")
    ax.set_title("scCODA-style ALR · Q4 vs Q1 (primary)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def fig_honest_n(n_tab: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.2), constrained_layout=True)
    show = n_tab[n_tab["cohort"] != "merge_within_cohort"]
    x = np.arange(len(show))
    ax.bar(x - 0.18, show["n_q1"], width=0.36, color="#1f77b4", label="Q1")
    ax.bar(x + 0.18, show["n_q4"], width=0.36, color="#d62728", label="Q4")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}\n{t}" for c, t in zip(show["cohort"], show["unit_type"])])
    ax.set_ylabel("units")
    ax.set_title("Honest n · within-cohort CLDN4 quartiles")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--given", type=Path, default=ROOT / "data" / "given")
    p.add_argument("--outdir", type=Path, default=ROOT / "results")
    args = p.parse_args()
    out = args.outdir
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    comp = load_units(args.given)
    tests, dm_rows = run_tests(comp)
    n_tab = honest_n(comp)
    engine = try_sccoda()

    comp.to_csv(out / "tables" / "composition_table.tsv", sep="\t", index=False)
    tests.to_csv(out / "tables" / "sccoda_tests.tsv", sep="\t", index=False)
    tests[tests["primary_family"]].to_csv(out / "tables" / "sccoda_tests_primary.tsv", sep="\t", index=False)
    n_tab.to_csv(out / "tables" / "honest_n.tsv", sep="\t", index=False)
    pd.DataFrame(dm_rows).to_csv(out / "tables" / "dm_diagnostic.tsv", sep="\t", index=False)
    (out / "tables" / "sccoda_engine.json").write_text(json.dumps(engine, indent=2))

    fig_fractions(comp, out / "figures" / "fig_sccoda_fractions.png")
    fig_alr_forest(tests, out / "figures" / "fig_sccoda_alr_forest.png")
    fig_honest_n(n_tab, out / "figures" / "fig_honest_n.png")

    summary = {
        "n_units": int(len(comp)),
        "n_q4q1": int(comp["cldn4_q4q1"].isin([0, 1]).sum()),
        "n_q1": int((comp["quartile"] == "Q1").sum()),
        "n_q4": int((comp["quartile"] == "Q4").sum()),
        "engine": engine,
        "n_primary_recover_q10": int(tests["recover_q10"].fillna(False).sum()) if len(tests) else 0,
        "gse148071_added": False,
        "seven_cohort_pool": False,
        "dual_high": False,
    }
    (out / "tables" / "sccoda_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    print(tests[tests["split"] == "q4q1"][
        ["slice", "compartment", "n", "n_high", "n_low", "alr_effect", "alr_p_perm", "q_bh_primary"]
    ].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
