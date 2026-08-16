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
from scipy.stats import fisher_exact, mannwhitneyu


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


def make_plot(primary: pd.DataFrame, p_value: float, output_dir: Path) -> None:
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
        f"Two-sided Mann–Whitney p = {p_value:.3f}",
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
    with (args.output_dir / "primary_result.json").open("w") as handle:
        json.dump(primary, handle, indent=2)

    days_after_ipi = pd.to_numeric(data["daysBiopsyAfterIpiStart"], errors="coerce")
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

    make_plot(data.loc[is_r | is_pd].copy(), primary["p_value_two_sided"], args.output_dir)
    print(json.dumps(primary, indent=2))


if __name__ == "__main__":
    main()
