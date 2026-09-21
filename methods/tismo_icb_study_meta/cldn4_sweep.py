#!/usr/bin/env python3
"""Cldn4-only specification sweep on the locked TISMO pairing.

Tacstd2 49/64 is not touched. Every Cldn4 row is a real filter of the same 64
slices. The sweep looks for a thesis-aligned increase (ICB mean above baseline).
It does not add studies, relabel LLC, or drop a negative line and call that
the result.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import analyze

OUT = analyze.OUT
FIG = analyze.FIG
TABLES = analyze.TABLES

# Extra modality on top of checkpoint blockade. Fc-silent anti-PD1 alone stays ICB-only.
EXTRA_AGENT = re.compile(
    r"birinapant|regorafenib|radiation|antiTGFb|antiGARP|TGFb|IL17|exercise|high\.fat|BRAF|swainsonine",
    re.I,
)
# Engineered or mutagenized genotype relative to the parental line in the same export.
PERTURBED = re.compile(r"KO|knockdown|Setdb1|Brca|Aire|Man2a1|MDK|Apobec|_UV(?:_|$)", re.I)

FLOORS = (
    "none",
    "max_arm_ge_0.25",
    "max_arm_ge_0.5",
    "max_arm_ge_1",
    "baseline_ge_0.5",
    "baseline_ge_1",
)
REGIMENS = ("all", "icb_only", "icb_plus_other")
GENOTYPES = ("all", "parental", "perturbed")


def annotate_slices(slices: pd.DataFrame) -> pd.DataFrame:
    out = slices.copy()
    out["icb_only"] = ~out["stem"].str.contains(EXTRA_AGENT)
    out["perturbed"] = out["stem"].str.contains(PERTURBED)
    out["max_arm"] = np.maximum(out["mean_baseline"], out["mean_icb"])
    return out


def apply_floor(df: pd.DataFrame, floor: str) -> pd.DataFrame:
    if floor == "none":
        return df
    if floor == "max_arm_ge_0.25":
        return df[df["max_arm"] >= 0.25]
    if floor == "max_arm_ge_0.5":
        return df[df["max_arm"] >= 0.5]
    if floor == "max_arm_ge_1":
        return df[df["max_arm"] >= 1.0]
    if floor == "baseline_ge_0.5":
        return df[df["mean_baseline"] >= 0.5]
    if floor == "baseline_ge_1":
        return df[df["mean_baseline"] >= 1.0]
    raise ValueError(floor)


def apply_regimen(df: pd.DataFrame, regimen: str) -> pd.DataFrame:
    if regimen == "all":
        return df
    if regimen == "icb_only":
        return df[df["icb_only"]]
    if regimen == "icb_plus_other":
        return df[~df["icb_only"]]
    raise ValueError(regimen)


def apply_genotype(df: pd.DataFrame, genotype: str) -> pd.DataFrame:
    if genotype == "all":
        return df
    if genotype == "parental":
        return df[~df["perturbed"]]
    if genotype == "perturbed":
        return df[df["perturbed"]]
    raise ValueError(genotype)


def _nan_fit() -> dict:
    keys = (
        "k", "mu", "se", "ci_low", "ci_high", "p_wald", "tau2", "i2",
        "mkh_low", "mkh_high", "p_mkh",
    )
    return {k: math.nan for k in keys}


def reml_on_slices(df: pd.DataFrame, variance: str) -> dict:
    """Study-level REML on whatever slices remain. variance is primary or technical."""
    if df["gse_id"].nunique() < 3 or len(df) < 4:
        return _nan_fit()
    try:
        borrowed = analyze.borrowed_arm_variance(df)
        studies = analyze.study_estimates(df, borrowed)
    except (RuntimeError, ValueError):
        return _nan_fit()
    y = studies["y"].to_numpy(dtype=float)
    v_tech = studies["v"].to_numpy(dtype=float)
    if variance == "technical":
        v = v_tech
    else:
        try:
            sigma = analyze.pooled_slice_variance(df)
        except RuntimeError:
            sigma = float(np.median(v_tech))
        v = np.maximum(v_tech, sigma / studies["n_slices"].to_numpy(dtype=float))
    if np.any(~np.isfinite(v)) or np.any(v <= 0):
        return _nan_fit()
    fit = analyze.fit_meta(y, v, method="REML")
    return {k: fit[k] for k in _nan_fit()}


def score_subset(df: pd.DataFrame) -> dict:
    empty = {
        "n_slices": int(len(df)),
        "n_studies": int(df["gse_id"].nunique()) if len(df) else 0,
        "n_lines": int(df["line"].nunique()) if len(df) else 0,
        "n_up": 0,
        "n_down": 0,
        "n_tie": 0,
        "mean_delta": math.nan,
        "median_delta": math.nan,
        "slice_wilcoxon_p": math.nan,
        "slice_wilcoxon_p_greater": math.nan,
        "slice_binomial_p": math.nan,
        "slice_binomial_p_greater": math.nan,
        "study_n_up": 0,
        "study_n": 0,
        "study_wilcoxon_p": math.nan,
        "study_binomial_p": math.nan,
        "one_vote_mean": math.nan,
        "one_vote_p": math.nan,
        "reml_mu": math.nan,
        "reml_ci_low": math.nan,
        "reml_ci_high": math.nan,
        "reml_p_wald": math.nan,
        "reml_p_mkh": math.nan,
        "reml_i2": math.nan,
        "technical_reml_mu": math.nan,
        "technical_reml_p_wald": math.nan,
    }
    if len(df) == 0:
        return empty
    deltas = df["delta_mean"].to_numpy(dtype=float)
    signs = analyze.sign_test(deltas)
    wil = analyze.wilcoxon_signed(deltas)
    empty.update(
        {
            "n_up": signs["n_up"],
            "n_down": signs["n_down"],
            "n_tie": signs["n_tie"],
            "mean_delta": float(np.mean(deltas)),
            "median_delta": float(np.median(deltas)),
            "slice_wilcoxon_p": wil["p_two_sided"],
            "slice_wilcoxon_p_greater": wil["p_greater"],
            "slice_binomial_p": signs["binomial_p_two_sided"],
            "slice_binomial_p_greater": signs["binomial_p_greater"],
        }
    )
    if df["gse_id"].nunique() >= 1:
        study_means = df.groupby("gse_id")["delta_mean"].mean()
        sy = study_means.to_numpy(dtype=float)
        ss = analyze.sign_test(sy)
        sw = analyze.wilcoxon_signed(sy)
        empty["study_n_up"] = ss["n_up"]
        empty["study_n"] = ss["n"]
        empty["study_wilcoxon_p"] = sw["p_two_sided"]
        empty["study_binomial_p"] = ss["binomial_p_two_sided"]
        empty["one_vote_mean"] = float(np.mean(sy))
        if len(sy) >= 3 and float(np.std(sy, ddof=1)) > 0:
            se = float(np.std(sy, ddof=1) / math.sqrt(len(sy)))
            empty["one_vote_p"] = float(2 * analyze.stats.t.sf(abs(np.mean(sy) / se), len(sy) - 1))
    primary = reml_on_slices(df, "primary")
    technical = reml_on_slices(df, "technical")
    empty.update(
        {
            "reml_mu": primary["mu"],
            "reml_ci_low": primary["ci_low"],
            "reml_ci_high": primary["ci_high"],
            "reml_p_wald": primary["p_wald"],
            "reml_p_mkh": primary["p_mkh"],
            "reml_i2": primary["i2"],
            "technical_reml_mu": technical["mu"],
            "technical_reml_p_wald": technical["p_wald"],
        }
    )
    return empty


def build_grid(slices: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for floor in FLOORS:
        for regimen in REGIMENS:
            for genotype in GENOTYPES:
                df = apply_genotype(apply_regimen(apply_floor(slices, floor), regimen), genotype)
                rec = {
                    "family": "biology",
                    "floor": floor,
                    "regimen": regimen,
                    "genotype": genotype,
                    "subset": "all_lines",
                    "spec": f"{floor}|{regimen}|{genotype}",
                }
                rec.update(score_subset(df))
                rows.append(rec)
    for cancer in sorted(slices["cancer_group"].unique()):
        base = slices[slices["cancer_group"] == cancer]
        for floor in ("none", "max_arm_ge_0.5", "baseline_ge_1"):
            df = apply_floor(base, floor)
            rec = {
                "family": "histology",
                "floor": floor,
                "regimen": "all",
                "genotype": "all",
                "subset": cancer,
                "spec": f"histology:{cancer}|{floor}",
            }
            rec.update(score_subset(df))
            rows.append(rec)
    for line in sorted(slices["line"].unique()):
        df = slices[slices["line"] != line]
        rec = {
            "family": "influence_drop_line",
            "floor": "none",
            "regimen": "all",
            "genotype": "all",
            "subset": f"drop_{line}",
            "spec": f"drop_line:{line}",
        }
        rec.update(score_subset(df))
        rows.append(rec)
    return pd.DataFrame(rows)


def positive_slice_p(row: pd.Series) -> float:
    """Two-sided slice Wilcoxon p, kept only when the mean delta is positive."""
    if not np.isfinite(row["slice_wilcoxon_p"]) or not np.isfinite(row["mean_delta"]):
        return math.inf
    if row["mean_delta"] <= 0 or row["n_slices"] < 10:
        return math.inf
    return float(row["slice_wilcoxon_p"])


def sweep_summary(grid: pd.DataFrame) -> dict:
    biology = grid[grid["family"] == "biology"].copy()
    biology["rank_p"] = biology.apply(positive_slice_p, axis=1)
    best = biology.sort_values(["rank_p", "n_slices"], ascending=[True, False]).iloc[0]
    expressed = biology[biology["floor"] == "baseline_ge_1"]
    influence = grid[grid["family"] == "influence_drop_line"].copy()
    influence["rank_p"] = influence.apply(positive_slice_p, axis=1)
    best_drop = influence.sort_values("rank_p").iloc[0]
    n_biology = int(len(biology))
    n_lt_05 = int((biology["slice_wilcoxon_p"] < 0.05).sum())
    n_greater_lt_05 = int(
        ((biology["slice_wilcoxon_p_greater"] < 0.05) & (biology["mean_delta"] > 0) & (biology["n_slices"] >= 10)).sum()
    )
    return {
        "n_biology": n_biology,
        "n_biology_slice_p_lt_0.05": n_lt_05,
        "n_biology_greater_p_lt_0.05_n_ge_10": n_greater_lt_05,
        "bonferroni_0.05": 0.05 / n_biology,
        "best_biology": best.to_dict(),
        "best_drop_line": best_drop.to_dict(),
        "expressed_baseline_ge_1_all": biology[
            (biology["floor"] == "baseline_ge_1") & (biology["regimen"] == "all") & (biology["genotype"] == "all")
        ].iloc[0].to_dict(),
        "locked_all": biology[
            (biology["floor"] == "none") & (biology["regimen"] == "all") & (biology["genotype"] == "all")
        ].iloc[0].to_dict(),
    }


def write_sweep_figure(grid: pd.DataFrame, path: Path) -> None:
    biology = grid[grid["family"] == "biology"].copy()
    biology = biology[biology["n_slices"] >= 10].sort_values("mean_delta")
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    ax.axvline(0, color="0.4", lw=0.8)
    ys = np.arange(len(biology))
    colors = np.where(biology["slice_wilcoxon_p"] < 0.05, "#E45756", "#4C78A8")
    ax.scatter(biology["mean_delta"], ys, c=colors, s=28, zorder=3)
    for y, (_, row) in zip(ys, biology.iterrows()):
        ax.plot([row["median_delta"], row["mean_delta"]], [y, y], color="0.75", lw=0.6, zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(biology["spec"], fontsize=6)
    ax.set_xlabel("Cldn4 mean slice Δ (dot) and median (tick)")
    ax.set_title("Cldn4 biology sweep (n≥10 slices)\nblue = two-sided Wilcoxon p≥0.05")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _fmt_row(row: dict) -> str:
    p = row["slice_wilcoxon_p"]
    return (
        f"{row['spec']}: {int(row['n_up'])}/{int(row['n_slices'])} up, "
        f"mean Δ {row['mean_delta']:+.3f}, median {row['median_delta']:+.3f}, "
        f"Wilcoxon p = {analyze.fmt_p(p)}, one-sided greater p = {analyze.fmt_p(row['slice_wilcoxon_p_greater'])}"
    )


def render_section(grid: pd.DataFrame, summary: dict) -> str:
    locked = summary["locked_all"]
    best = summary["best_biology"]
    drop = summary["best_drop_line"]
    expressed = summary["expressed_baseline_ge_1_all"]
    lines = []
    a = lines.append
    a("")
    a("## Cldn4 specification sweep")
    a("")
    a("Tacstd2 stays at 49/64. This section changes only Cldn4 weights, slice filters, expression-floor rules, and genotype subsets. No slice was added. LLC is still LLC. The full grid is `results/tismo_icb_study_meta/tables/cldn4_sweep.tsv`.")
    a("")
    a("Floor rules use the TISMO log2(TPM+1) arm means. `max_arm` is the larger of the baseline and ICB means. `baseline_ge_1` keeps slices whose naive mean is at least 1 (TPM about 1). ICB-only drops slices whose label also names another modality (radiation, birinapant, regorafenib, TGF-β or GARP blockade, IL-17 blockade, exercise, high-fat diet, BRAF inhibitor, swainsonine). Parental drops KO, knockdown, Setdb1, Brca, Aire, Man2a1, MDK, Apobec, and UV labels. Perturbed is the complement.")
    a("")
    a(
        f"Biology grid: **{summary['n_biology']}** specifications (6 floors × 3 regimen rules × 3 genotype rules). "
        f"Two-sided slice Wilcoxon p < 0.05 in **{summary['n_biology_slice_p_lt_0.05']}** of them. "
        f"Bonferroni threshold for this grid: {summary['bonferroni_0.05']:.4g}. "
        f"Among specs with at least 10 slices and a positive mean, one-sided greater p < 0.05 in **{summary['n_biology_greater_p_lt_0.05_n_ge_10']}**."
    )
    a("")
    a(f"Locked cell (no filter): {_fmt_row(locked)}.")
    a("")
    a(f"Smallest two-sided slice p in that biology grid, requiring a positive mean and at least 10 slices: {_fmt_row(best)}.")
    a("")
    hits = grid[(grid["family"] == "biology") & (grid["slice_wilcoxon_p"] < 0.05)].copy()
    hit_regimens = sorted(hits["regimen"].unique())
    a(
        f"Rows with two-sided slice p < 0.05: {len(hits)}. Every one of them has regimen `{', '.join(hit_regimens)}`. "
        "They are not six different findings. Parental and all genotypes are the same slices in this subset, and the max-arm floors of 0.25 and 0.5 keep the same 8 slices. "
        "The two unique sets are the 12 ICB-plus-other slices (8 up, p = 0.027) and those 12 after dropping 4 slices whose higher arm mean is below 0.25 (7/8 up, p = 0.016). "
        "Those 12 labels are radiation, birinapant, swainsonine, regorafenib, GARP:TGF-β (two Fc designs), TGF-β, high-fat diet, TGF-β trap, exercise, IL-17 blockade, and a BRAF inhibitor, each on top of a checkpoint antibody. "
        "REML on the 8-slice set is μ = +0.426 (Wald p = 2.6×10⁻⁴; modified Knapp–Hartung p = 0.019). "
        f"Neither slice p clears the biology-grid Bonferroni line ({summary['bonferroni_0.05']:.4g})."
    )
    a("")
    icb = grid[
        (grid["family"] == "biology")
        & (grid["floor"] == "none")
        & (grid["regimen"] == "icb_only")
        & (grid["genotype"] == "all")
    ].iloc[0]
    a(
        f"The thesis-aligned regimen is checkpoint blockade without a second modality. That row is {_fmt_row(icb.to_dict())}. "
        "Study-level tests on those slices stay above 0.05 as well. Adding a floor or a parental-only cut does not push the ICB-only two-sided p under 0.05 in any cell with at least 10 slices."
    )
    a("")
    a(
        f"Where Cldn4 is already on at baseline (baseline mean ≥ 1, all regimens, all genotypes): {_fmt_row(expressed)}. "
        "The increase does not get stronger in the slices where a tight-junction reading is possible. It gets weaker. "
        "High-baseline models include 4T1 and KPB25L, and those lines are not a consistent rise."
    )
    a("")
    a("### Weights on the unfiltered 64")
    a("")
    a("| Weight | Estimate | p |")
    a("|---|---|---:|")
    a(
        f"| Equal slice, Wilcoxon two-sided | {int(locked['n_up'])}/{int(locked['n_slices'])} up, mean {locked['mean_delta']:+.3f} | {analyze.fmt_p(locked['slice_wilcoxon_p'])} |"
    )
    a(
        f"| Equal slice, Wilcoxon greater | same counts | {analyze.fmt_p(locked['slice_wilcoxon_p_greater'])} |"
    )
    a(
        f"| Equal study, sign test | {int(locked['study_n_up'])}/{int(locked['study_n'])} studies up | {analyze.fmt_p(locked['study_binomial_p'])} |"
    )
    a(
        f"| Equal study, Wilcoxon | mean of study means {locked['one_vote_mean']:+.3f} | {analyze.fmt_p(locked['study_wilcoxon_p'])} |"
    )
    a(
        f"| Equal study, Student t on study means | {locked['one_vote_mean']:+.3f} | {analyze.fmt_p(locked['one_vote_p'])} |"
    )
    a(
        f"| REML, primary variance | μ {locked['reml_mu']:+.3f} ({locked['reml_ci_low']:+.3f} to {locked['reml_ci_high']:+.3f}) | Wald {analyze.fmt_p(locked['reml_p_wald'])} |"
    )
    a(
        f"| REML, technical SE | μ {locked['technical_reml_mu']:+.3f} | {analyze.fmt_p(locked['technical_reml_p_wald'])} |"
    )
    a("")
    a(
        "The one-sided greater p on the locked 64 is 0.038. That is the same signed-rank statistic as the two-sided lock (p = 0.077), read in the thesis direction. It is not a new sample. The two-sided lock stays the reference, because Tacstd2’s two-sided p is 5.84×10⁻⁵ and Cldn4’s is not in that range. "
        "The Student interval on 22 study means (p = 0.043) is the equal-study mean, already reported above, and it is pulled by a right skew (median study Δ near 0)."
    )
    a("")
    a("### What was not adopted")
    a("")
    a(
        f"Leaving out one cell line at a time, the smallest two-sided slice p is {_fmt_row(drop)}. "
        "That specification deletes a negative line. KPB25L is a mammary TISMO model inside the 64-slice lock. Removing it is influence, not a genotype subset and not an expression-floor rule. "
        "Dropping 4T1, the other high-Cldn4 mammary line, also moves the p-value down because 4T1’s deltas are mostly negative. Those rows are in the table under `influence_drop_line` and are not the result."
    )
    a("")
    a("Histology splits are in the same table (`family = histology`). They are small and are not a substitute for the 64-slice lock.")
    a("")
    a("**Sweep result.** The locked Cldn4 sentence stays 34/64, Wilcoxon p = 0.077. Checkpoint-only slices do not show a two-sided increase. The grid’s smallest p-values are checkpoint blockade plus a second named modality, on 8 or 12 slices, and they do not survive a Bonferroni cut of the 54-cell grid. Restricting to slices with baseline Cldn4 already on (mean ≥ 1) makes the paired change weaker (8/22, p = 0.61). Tacstd2 remains 49/64.")
    a("")
    a("Figure: `results/tismo_icb_study_meta/figures/cldn4_sweep.png`.")
    a("")
    return "\n".join(lines)


def run_sweep() -> dict:
    slices = annotate_slices(analyze.load_locked_slices()["slices"]["Cldn4"])
    # The sweep must see the same 64 stems as the lock.
    if len(slices) != 64:
        raise SystemExit(f"Cldn4 sweep expected 64 locked slices, found {len(slices)}")
    flags = slices[["stem", "line", "gse_id", "cancer_group", "is_llc", "mean_baseline", "mean_icb", "max_arm", "delta_mean", "icb_only", "perturbed"]].copy()
    TABLES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    flags.to_csv(TABLES / "cldn4_slice_flags.tsv", sep="\t", index=False)
    grid = build_grid(slices)
    grid.to_csv(TABLES / "cldn4_sweep.tsv", sep="\t", index=False)
    summary = sweep_summary(grid)
    write_sweep_figure(grid, FIG / "cldn4_sweep.png")
    section = render_section(grid, summary)
    for path in (analyze.ROOT / "RESULTS.md", OUT / "RESULTS.md"):
        text = path.read_text()
        marker = "\n## Cldn4 specification sweep\n"
        if marker in text:
            text = text.split(marker)[0].rstrip() + "\n"
        path.write_text(text.rstrip() + "\n" + section)
    print(
        f"Cldn4 sweep: biology n={summary['n_biology']}, "
        f"two-sided p<0.05: {summary['n_biology_slice_p_lt_0.05']}, "
        f"best {summary['best_biology']['spec']} p={summary['best_biology']['slice_wilcoxon_p']:.4g}"
    )
    return summary


if __name__ == "__main__":
    run_sweep()
