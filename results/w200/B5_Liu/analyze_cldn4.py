#!/usr/bin/env python3
"""Reproduce the CLDN4-versus-response analysis in Liu et al. (2019)."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact, mannwhitneyu, norm


TPM_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41591-019-0654-5/MediaObjects/"
    "41591_2019_654_MOESM3_ESM.txt"
)
CLINICAL_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41591-019-0654-5/MediaObjects/"
    "41591_2019_654_MOESM4_ESM.xlsx"
)
EXPECTED_SHA256 = {
    "Liu_RNA_TPM_matrix.txt":
        "80d3a6e6ec7a2c89101ed9f1c9078dda32c2f7152fb20d0e24edcadf37871553",
    "Liu_Supplementary_Tables.xlsx":
        "5bb2d19d5a73869d508d67323e84adccf00f948ec36a034e68efc1bb47b9e899",
}
SEED = 2019
N_BOOT = 20_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    if not destination.exists():
        print(f"Downloading {destination.name} ...")
        urllib.request.urlretrieve(url, destination)
    observed = sha256(destination)
    if observed != EXPECTED_SHA256[destination.name]:
        raise RuntimeError(
            f"Checksum mismatch for {destination}: expected "
            f"{EXPECTED_SHA256[destination.name]}, got {observed}"
        )


def compare(
    data: pd.DataFrame,
    responder_mask: pd.Series,
    comparator_mask: pd.Series,
    label: str,
) -> dict[str, float | int | str]:
    responder = data.loc[responder_mask, "log2_TPM_plus_1"].to_numpy()
    comparator = data.loc[comparator_mask, "log2_TPM_plus_1"].to_numpy()
    test = mannwhitneyu(responder, comparator, alternative="two-sided")
    auc = float(test.statistic / (len(responder) * len(comparator)))
    return {
        "analysis": label,
        "n_responder": len(responder),
        "n_comparator": len(comparator),
        "responder_median_log2_TPM_plus_1": float(np.median(responder)),
        "comparator_median_log2_TPM_plus_1": float(np.median(comparator)),
        "median_difference": float(np.median(responder) - np.median(comparator)),
        "rank_biserial": 2 * auc - 1,
        "auc_higher_CLDN4_predicts_response": auc,
        "mann_whitney_U": float(test.statistic),
        "p_value_two_sided": float(test.pvalue),
    }


def woolf_odds_ratio_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """Haldane–Anscombe corrected log-odds interval."""
    aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    log_or = np.log((aa * dd) / (bb * cc))
    se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return (
        float(np.exp(log_or)),
        float(np.exp(log_or - 1.959963984540054 * se)),
        float(np.exp(log_or + 1.959963984540054 * se)),
    )


def two_by_two(
    high: pd.Series,
    responder: pd.Series,
    analysis: str,
    predictor: str,
    outcome: str,
) -> dict[str, float | int | str]:
    a = int((high & responder).sum())
    b = int((high & ~responder).sum())
    c = int((~high & responder).sum())
    d = int((~high & ~responder).sum())
    table = [[a, b], [c, d]]
    fisher = fisher_exact(table, alternative="two-sided")
    _corrected_or, ci_low, ci_high = woolf_odds_ratio_ci(a, b, c, d)
    return {
        "analysis": analysis,
        "predictor": predictor,
        "outcome": outcome,
        "n": a + b + c + d,
        "n_high_responder": a,
        "n_high_nonresponder": b,
        "n_low_responder": c,
        "n_low_nonresponder": d,
        "odds_ratio": float(fisher.statistic),
        "or_ci_low": ci_low,
        "or_ci_high": ci_high,
        "p": float(fisher.pvalue),
        "method": "Fisher exact OR; Haldane–Anscombe 95% CI",
    }


def logistic_or_per_unit(
    x: pd.Series,
    y: pd.Series,
    analysis: str,
    predictor: str,
    outcome: str,
) -> dict[str, float | int | str]:
    design = np.column_stack([np.ones(len(x)), np.asarray(x, dtype=float)])
    response = np.asarray(y, dtype=float)
    beta = np.zeros(2)
    hessian = np.eye(2)
    for _ in range(50):
        linear = np.clip(design @ beta, -30, 30)
        fitted = 1.0 / (1.0 + np.exp(-linear))
        weight = fitted * (1.0 - fitted)
        gradient = design.T @ (response - fitted)
        hessian = design.T @ (weight[:, None] * design)
        try:
            step = np.linalg.solve(hessian, gradient)
        except np.linalg.LinAlgError:
            break
        beta = beta + step
        if np.max(np.abs(step)) < 1e-12:
            break
    covariance = np.linalg.inv(hessian)
    se = float(np.sqrt(covariance[1, 1]))
    z_stat = beta[1] / se
    return {
        "analysis": analysis,
        "predictor": predictor,
        "outcome": outcome,
        "n": int(len(response)),
        "n_high_responder": "",
        "n_high_nonresponder": "",
        "n_low_responder": "",
        "n_low_nonresponder": "",
        "odds_ratio": float(np.exp(beta[1])),
        "or_ci_low": float(np.exp(beta[1] - 1.959963984540054 * se)),
        "or_ci_high": float(np.exp(beta[1] + 1.959963984540054 * se)),
        "p": float(2 * (1 - norm.cdf(abs(z_stat)))),
        "method": "univariate logistic OR per +1 log2(TPM+1)",
    }


def holm_adjust(values: pd.Series) -> pd.Series:
    order = np.argsort(values.to_numpy())
    adjusted = np.empty(len(values), dtype=float)
    running_max = 0.0
    for rank, position in enumerate(order):
        candidate = min(1.0, (len(values) - rank) * float(values.iloc[position]))
        running_max = max(running_max, candidate)
        adjusted[position] = running_max
    return pd.Series(adjusted, index=values.index)


def bootstrap_primary(responder: np.ndarray, progressor: np.ndarray) -> dict[str, list[float]]:
    rng = np.random.default_rng(SEED)
    effects = np.empty((N_BOOT, 3))
    for i in range(N_BOOT):
        resample_r = rng.choice(responder, len(responder), replace=True)
        resample_p = rng.choice(progressor, len(progressor), replace=True)
        u_stat = mannwhitneyu(resample_r, resample_p, alternative="two-sided").statistic
        effects[i] = (
            2 * u_stat / (len(responder) * len(progressor)) - 1,
            np.median(resample_r) - np.median(resample_p),
            np.mean(resample_r) - np.mean(resample_p),
        )
    intervals = np.quantile(effects, [0.025, 0.975], axis=0)
    return {
        "rank_biserial_95pct_percentile_CI": intervals[:, 0].tolist(),
        "median_difference_95pct_percentile_CI": intervals[:, 1].tolist(),
        "mean_difference_95pct_percentile_CI": intervals[:, 2].tolist(),
    }


def make_forest(table: pd.DataFrame, output_dir: Path) -> None:
    plot_table = table[table["method"].str.startswith("Fisher")].copy().iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    y = np.arange(len(plot_table))
    ax.axvline(1.0, color="#888888", linewidth=1, linestyle="--")
    ax.errorbar(
        plot_table["odds_ratio"],
        y,
        xerr=[
            plot_table["odds_ratio"] - plot_table["or_ci_low"],
            plot_table["or_ci_high"] - plot_table["odds_ratio"],
        ],
        fmt="o",
        color="#2B8CBE",
        ecolor="#2B8CBE",
        capsize=3,
    )
    labels = [
        f"{row.analysis}\nOR={row.odds_ratio:.2f}, n={int(row.n)}, p={row.p:.3f}"
        for row in plot_table.itertuples()
    ]
    ax.set_yticks(y, labels)
    ax.set_xlabel("Odds ratio for response if CLDN4-high")
    ax.set_xscale("log")
    ax.set_xlim(0.15, 4.0)
    ax.set_title("Liu 2019 CLDN4 versus anti-PD-1 response")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_dir / "or_n_p.png", dpi=220)
    fig.savefig(output_dir / "or_n_p.svg")
    plt.close(fig)


def make_plot(
    primary: pd.DataFrame,
    p_value: float,
    or_n_p: dict[str, float | int | str],
    output_dir: Path,
) -> None:
    groups = ["Responder (CR/PR)", "Progressor (PD)"]
    colors = ["#2B8CBE", "#D95F0E"]
    rng = np.random.default_rng(SEED)
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    for position, (group, color) in enumerate(zip(groups, colors), start=1):
        values = primary.loc[primary["primary_group"] == group, "log2_TPM_plus_1"]
        jitter = rng.uniform(-0.12, 0.12, len(values))
        ax.scatter(
            position + jitter,
            values,
            s=30,
            alpha=0.72,
            color=color,
            edgecolor="white",
            linewidth=0.35,
            zorder=3,
        )
        quartiles = np.quantile(values, [0.25, 0.5, 0.75])
        ax.vlines(position, quartiles[0], quartiles[2], color="black", linewidth=5)
        ax.scatter(position, quartiles[1], marker="_", s=240, color="white", zorder=4)
    ax.set_xticks([1, 2], [f"{g}\n(n={sum(primary.primary_group == g)})" for g in groups])
    ax.set_ylabel("CLDN4 expression, log$_2$(TPM + 1)")
    ax.set_title("Liu 2019 pretreatment melanoma RNA-seq")
    ax.text(
        0.98,
        0.97,
        (
            f"Median-split OR = {or_n_p['odds_ratio']:.2f} "
            f"(n={or_n_p['n']}, p={or_n_p['p']:.3f})\n"
            f"Two-sided Mann–Whitney p = {p_value:.3f}"
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.7, zorder=0)
    fig.tight_layout()
    fig.savefig(output_dir / "cldn4_vs_response.png", dpi=220)
    fig.savefig(output_dir / "cldn4_vs_response.svg")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("/tmp/liu2019_cldn4_source"),
        help="Location for downloaded public supplements.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Location for generated results.",
    )
    args = parser.parse_args()
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tpm_path = args.cache_dir / "Liu_RNA_TPM_matrix.txt"
    clinical_path = args.cache_dir / "Liu_Supplementary_Tables.xlsx"
    download(TPM_URL, tpm_path)
    download(CLINICAL_URL, clinical_path)

    clinical = pd.read_excel(
        clinical_path,
        sheet_name="Supplemental Table 1",
        header=2,
        index_col=0,
    ).iloc[:144]
    tpm = pd.read_csv(
        tpm_path,
        sep="\t",
        index_col=0,
        usecols=["Unnamed: 0", "CLDN4"],
    )
    data = tpm.join(
        clinical[["BR", "Primary_Type", "daysBiopsyAfterIpiStart"]],
        how="inner",
    )
    data.index.name = "sample_id"
    data = data.rename(columns={"BR": "best_response", "CLDN4": "CLDN4_TPM"})
    data["log2_TPM_plus_1"] = np.log2(data["CLDN4_TPM"] + 1)
    data["primary_group"] = data["best_response"].map(
        {"CR": "Responder (CR/PR)", "PR": "Responder (CR/PR)", "PD": "Progressor (PD)"}
    )
    data.to_csv(args.output_dir / "cldn4_patient_data.csv")

    is_r = data["best_response"].isin(["CR", "PR"])
    is_pd = data["best_response"].eq("PD")
    primary = compare(data, is_r, is_pd, "Primary: CR/PR vs PD")
    responder = data.loc[is_r, "log2_TPM_plus_1"].to_numpy()
    progressor = data.loc[is_pd, "log2_TPM_plus_1"].to_numpy()
    primary.update(bootstrap_primary(responder, progressor))
    primary["responder_median_TPM"] = float(data.loc[is_r, "CLDN4_TPM"].median())
    primary["progressor_median_TPM"] = float(data.loc[is_pd, "CLDN4_TPM"].median())
    detection_table = [
        [int((data.loc[is_r, "CLDN4_TPM"] > 0).sum()), int((data.loc[is_r, "CLDN4_TPM"] == 0).sum())],
        [int((data.loc[is_pd, "CLDN4_TPM"] > 0).sum()), int((data.loc[is_pd, "CLDN4_TPM"] == 0).sum())],
    ]
    detection_test = fisher_exact(detection_table, alternative="two-sided")
    primary.update(
        {
            "responder_detected_n": detection_table[0][0],
            "progressor_detected_n": detection_table[1][0],
            "detection_odds_ratio": float(detection_test.statistic),
            "detection_fisher_p_two_sided": float(detection_test.pvalue),
            "bootstrap_seed": SEED,
            "bootstrap_replicates": N_BOOT,
        }
    )

    primary_subset = data.loc[is_r | is_pd].copy()
    primary_high = primary_subset["CLDN4_TPM"] > primary_subset["CLDN4_TPM"].median()
    primary_response = primary_subset["best_response"].isin(["CR", "PR"])
    orr_high = data["CLDN4_TPM"] > data["CLDN4_TPM"].median()
    orr_response = data["best_response"].isin(["CR", "PR"])
    days_after_ipi = pd.to_numeric(data["daysBiopsyAfterIpiStart"], errors="coerce")
    or_rows = [
        two_by_two(
            primary_high,
            primary_response,
            "Primary median-split high vs low",
            "CLDN4 > cohort median TPM",
            "CR/PR vs PD",
        ),
        two_by_two(
            primary_subset["CLDN4_TPM"] > 0,
            primary_response,
            "Primary detected vs undetected",
            "CLDN4 TPM > 0",
            "CR/PR vs PD",
        ),
        two_by_two(
            primary_subset["CLDN4_TPM"] > primary_subset["CLDN4_TPM"].quantile(0.75),
            primary_response,
            "Primary top-quartile vs rest",
            "CLDN4 > 75th percentile TPM",
            "CR/PR vs PD",
        ),
        logistic_or_per_unit(
            primary_subset["log2_TPM_plus_1"],
            primary_response,
            "Primary logistic per +1 log2(TPM+1)",
            "continuous log2(TPM+1)",
            "CR/PR vs PD",
        ),
        two_by_two(
            orr_high,
            orr_response,
            "All RNA median-split high vs low",
            "CLDN4 > cohort median TPM",
            "CR/PR vs PD/SD/MR",
        ),
        two_by_two(
            data["CLDN4_TPM"] > 0,
            orr_response,
            "All RNA detected vs undetected",
            "CLDN4 TPM > 0",
            "CR/PR vs PD/SD/MR",
        ),
        two_by_two(
            data.loc[(is_r | is_pd) & data["daysBiopsyAfterIpiStart"].eq("na"), "CLDN4_TPM"]
            > data.loc[(is_r | is_pd) & data["daysBiopsyAfterIpiStart"].eq("na"), "CLDN4_TPM"].median(),
            data.loc[(is_r | is_pd) & data["daysBiopsyAfterIpiStart"].eq("na"), "best_response"].isin(["CR", "PR"]),
            "Ipilimumab-naive median-split",
            "CLDN4 > subgroup median TPM",
            "CR/PR vs PD",
        ),
        two_by_two(
            data.loc[(is_r | is_pd) & (days_after_ipi > 0), "CLDN4_TPM"]
            > data.loc[(is_r | is_pd) & (days_after_ipi > 0), "CLDN4_TPM"].median(),
            data.loc[(is_r | is_pd) & (days_after_ipi > 0), "best_response"].isin(["CR", "PR"]),
            "Post-ipilimumab median-split",
            "CLDN4 > subgroup median TPM",
            "CR/PR vs PD",
        ),
        two_by_two(
            data.loc[(is_r | is_pd) & data["Primary_Type"].isin(["skin", "occult"]), "CLDN4_TPM"]
            > data.loc[(is_r | is_pd) & data["Primary_Type"].isin(["skin", "occult"]), "CLDN4_TPM"].median(),
            data.loc[(is_r | is_pd) & data["Primary_Type"].isin(["skin", "occult"]), "best_response"].isin(["CR", "PR"]),
            "Skin/occult median-split",
            "CLDN4 > subgroup median TPM",
            "CR/PR vs PD",
        ),
        two_by_two(
            data.loc[(is_r | is_pd) & data["Primary_Type"].eq("skin"), "CLDN4_TPM"]
            > data.loc[(is_r | is_pd) & data["Primary_Type"].eq("skin"), "CLDN4_TPM"].median(),
            data.loc[(is_r | is_pd) & data["Primary_Type"].eq("skin"), "best_response"].isin(["CR", "PR"]),
            "Skin-only median-split",
            "CLDN4 > subgroup median TPM",
            "CR/PR vs PD",
        ),
    ]
    or_table = pd.DataFrame(or_rows)
    or_table["holm_p_across_listed_ORs"] = holm_adjust(or_table["p"])
    or_table.to_csv(args.output_dir / "or_n_p.csv", index=False)
    markdown_lines = [
        "| analysis | n | OR | 95% CI | p | 2x2 high-R / high-NR / low-R / low-NR |",
        "|---|---:|---:|---|---:|---|",
    ]
    for row in or_table.itertuples():
        cells = (
            f"{row.n_high_responder}/{row.n_high_nonresponder}/"
            f"{row.n_low_responder}/{row.n_low_nonresponder}"
            if row.n_high_responder != ""
            else "continuous"
        )
        markdown_lines.append(
            f"| {row.analysis} | {int(row.n)} | {row.odds_ratio:.3f} | "
            f"{row.or_ci_low:.3f}–{row.or_ci_high:.3f} | {row.p:.3f} | {cells} |"
        )
    (args.output_dir / "or_n_p.md").write_text("\n".join(markdown_lines) + "\n")
    headline = or_rows[0]
    primary.update(
        {
            "n": headline["n"],
            "odds_ratio": headline["odds_ratio"],
            "or_ci_low": headline["or_ci_low"],
            "or_ci_high": headline["or_ci_high"],
            "or_p": headline["p"],
            "or_definition": (
                "response odds if CLDN4 > primary-cohort median TPM versus below/equal"
            ),
            "or_2x2_high_responder": headline["n_high_responder"],
            "or_2x2_high_progressor": headline["n_high_nonresponder"],
            "or_2x2_low_responder": headline["n_low_responder"],
            "or_2x2_low_progressor": headline["n_low_nonresponder"],
        }
    )
    with (args.output_dir / "primary_result.json").open("w") as handle:
        json.dump(primary, handle, indent=2)

    checks = [
        compare(data, is_r, data["best_response"].isin(["PD", "SD", "MR"]),
                "Objective response: CR/PR vs PD/SD/MR"),
        compare(data, is_r & (days_after_ipi > 0), is_pd & (days_after_ipi > 0),
                "Post-ipilimumab biopsy: CR/PR vs PD"),
        compare(data, is_r & data["daysBiopsyAfterIpiStart"].eq("na"),
                is_pd & data["daysBiopsyAfterIpiStart"].eq("na"),
                "Ipilimumab-naive: CR/PR vs PD"),
        compare(data, is_r & data["Primary_Type"].isin(["skin", "occult"]),
                is_pd & data["Primary_Type"].isin(["skin", "occult"]),
                "Skin/occult primary: CR/PR vs PD"),
        compare(data, is_r & data["Primary_Type"].eq("skin"),
                is_pd & data["Primary_Type"].eq("skin"),
                "Skin primary only: CR/PR vs PD"),
    ]
    sensitivity = pd.DataFrame(checks)
    sensitivity["holm_p_across_listed_sensitivities"] = holm_adjust(
        sensitivity["p_value_two_sided"]
    )
    sensitivity.to_csv(args.output_dir / "sensitivity_results.csv", index=False)

    published = pd.read_excel(
        clinical_path,
        sheet_name="Supplemental Table 4",
        header=2,
        index_col=0,
    )
    published_p = float(published.loc["CLDN4", "Overall R vs PD - MWW p-val"])
    if not np.isclose(primary["p_value_two_sided"], published_p):
        raise RuntimeError("Primary p-value does not reproduce Supplementary Table 4.")

    source_manifest = {
        "study": "Liu, Schilling, et al., Nature Medicine (2019)",
        "doi": "10.1038/s41591-019-0654-5",
        "sources": [
            {"url": TPM_URL, "sha256": sha256(tpm_path), "role": "RNA-seq TPM matrix"},
            {
                "url": CLINICAL_URL,
                "sha256": sha256(clinical_path),
                "role": "response and clinical annotations",
            },
        ],
        "software": {
            "python": ".".join(map(str, __import__("sys").version_info[:3])),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": __import__("matplotlib").__version__,
        },
        "validation": {
            "published_supplementary_table_4_CLDN4_p": published_p,
            "reproduced_primary_p": primary["p_value_two_sided"],
        },
    }
    with (args.output_dir / "source_manifest.json").open("w") as handle:
        json.dump(source_manifest, handle, indent=2)

    make_plot(
        data.loc[is_r | is_pd].copy(),
        primary["p_value_two_sided"],
        headline,
        args.output_dir,
    )
    make_forest(or_table, args.output_dir)
    print(json.dumps(primary, indent=2))


if __name__ == "__main__":
    main()
