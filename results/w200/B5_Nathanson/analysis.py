#!/usr/bin/env python3
"""CLDN4 versus clinical benefit in the public Nathanson 2017 melanoma ICI cohort."""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

REVISION = "2ce35b97c15a54d27090ac7e6b2e73f140e087b3"
BASE_URL = (
    "https://raw.githubusercontent.com/hammerlab/melanoma-reanalysis/"
    f"{REVISION}/files"
)
SOURCES = {
    "cohort.csv": {
        "url": f"{BASE_URL}/cohort.csv",
        "sha256": "33a3978ad365771a50e04dc1ac8308404627b63e978f7ae321f23637e0dd8eeb",
    },
    "cufflinks.tar.gz": {
        "url": f"{BASE_URL}/cufflinks.tar.gz",
        "sha256": "67f5fd34c6643845892eb16a50e5fa7079c165817c29963f5255a5464ce2d2bd",
    },
}
CANONICAL_CLDN4 = "ENSG00000189143"
SEED = 2017
BOOTSTRAP_DRAWS = 20_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch_sources(cache: Path) -> dict[str, Path]:
    cache.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, source in SOURCES.items():
        path = cache / name
        if not path.exists() or sha256(path) != source["sha256"]:
            urllib.request.urlretrieve(source["url"], path)
        actual = sha256(path)
        if actual != source["sha256"]:
            raise ValueError(f"Checksum mismatch for {name}: {actual}")
        paths[name] = path
    return paths


def load_cldn4_rows(archive: Path) -> pd.DataFrame:
    columns = ["tracking_id", "gene_short_name", "FPKM", "sample", "MSK_ID"]
    with tarfile.open(archive, "r:gz") as tar:
        handle = tar.extractfile(tar.getmember("cufflinks.csv"))
        if handle is None:
            raise FileNotFoundError("cufflinks.csv is missing from archive")
        chunks = [
            chunk.loc[chunk["gene_short_name"].eq("CLDN4")]
            for chunk in pd.read_csv(handle, usecols=columns, chunksize=250_000)
        ]
    rows = pd.concat(chunks, ignore_index=True)
    if rows.empty:
        raise ValueError("CLDN4 is absent from the expression data")
    return rows


def normalize_sample_id(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lstrip("0")


def collapse_cldn4(rows: pd.DataFrame, method: str) -> pd.DataFrame:
    if method == "canonical":
        chosen = rows.loc[rows["tracking_id"].eq(CANONICAL_CLDN4)].copy()
        if chosen.empty:
            raise ValueError(f"{CANONICAL_CLDN4} is missing")
        chosen = chosen.rename(columns={"tracking_id": "cldn4_id"})
        return chosen[["sample", "MSK_ID", "cldn4_id", "FPKM"]]
    if method == "symbol_median":
        collapsed = rows.groupby("sample", as_index=False).agg(
            FPKM=("FPKM", "median"),
            MSK_ID=("MSK_ID", "first"),
            n_gene_ids=("tracking_id", "nunique"),
        )
        collapsed["cldn4_id"] = "symbol_median"
        return collapsed
    raise ValueError(f"Unknown collapse method: {method}")


def bootstrap_auc_ci(
    responders: np.ndarray, nonresponders: np.ndarray
) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    estimates = np.empty(BOOTSTRAP_DRAWS)
    n_r = len(responders)
    n_nr = len(nonresponders)
    for i in range(BOOTSTRAP_DRAWS):
        a = rng.choice(responders, n_r, replace=True)
        b = rng.choice(nonresponders, n_nr, replace=True)
        estimates[i] = mannwhitneyu(a, b, alternative="two-sided").statistic / (
            n_r * n_nr
        )
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def summarize(samples: pd.DataFrame, label: str) -> dict:
    benefit = samples.loc[samples["Benefit"], "log2_FPKM_plus_1"].to_numpy()
    no_benefit = samples.loc[~samples["Benefit"], "log2_FPKM_plus_1"].to_numpy()
    test = mannwhitneyu(benefit, no_benefit, alternative="two-sided", method="exact")
    n_product = len(benefit) * len(no_benefit)
    auc = float(test.statistic / n_product)
    auc_lo, auc_hi = bootstrap_auc_ci(benefit, no_benefit)
    detected_b = int(np.sum(samples.loc[samples["Benefit"], "FPKM"] > 0))
    detected_nb = int(np.sum(samples.loc[~samples["Benefit"], "FPKM"] > 0))
    return {
        "analysis": label,
        "n_benefit": int(len(benefit)),
        "n_no_benefit": int(len(no_benefit)),
        "median_FPKM_benefit": float(samples.loc[samples["Benefit"], "FPKM"].median()),
        "median_FPKM_no_benefit": float(
            samples.loc[~samples["Benefit"], "FPKM"].median()
        ),
        "median_log2_benefit": float(np.median(benefit)),
        "median_log2_no_benefit": float(np.median(no_benefit)),
        "mann_whitney_U": float(test.statistic),
        "exact_two_sided_p": float(test.pvalue),
        "auc_higher_predicts_benefit": auc,
        "auc_bootstrap_95ci_low": float(auc_lo),
        "auc_bootstrap_95ci_high": float(auc_hi),
        "rank_biserial": float(2 * auc - 1),
        "detected_benefit": detected_b,
        "detected_no_benefit": detected_nb,
    }


def leave_one_out(samples: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for study_id in samples["Study ID"]:
        subset = samples.loc[samples["Study ID"].ne(study_id)]
        if subset["Benefit"].nunique() < 2:
            continue
        stats = summarize(subset, f"loo_without_{study_id}")
        stats["dropped_study_id"] = study_id
        rows.append(stats)
    return pd.DataFrame(rows)


def plot_panels(samples: pd.DataFrame, output: Path) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2))
    rng = np.random.default_rng(SEED)
    panels = [
        (
            "Pretreatment (primary)",
            samples.loc[samples["Biopsy pre or post ipi"].eq("pre")],
        ),
        (
            "Post-treatment (exploratory)",
            samples.loc[samples["Biopsy pre or post ipi"].eq("post")],
        ),
    ]
    colors = ["#6b7280", "#2563eb"]
    for ax, (title, frame) in zip(axes, panels, strict=True):
        groups = [
            frame.loc[~frame["Benefit"], "log2_FPKM_plus_1"].to_numpy(),
            frame.loc[frame["Benefit"], "log2_FPKM_plus_1"].to_numpy(),
        ]
        labels = [
            f"No benefit\n(n={len(groups[0])})",
            f"Benefit\n(n={len(groups[1])})",
        ]
        box = ax.boxplot(
            groups,
            positions=[1, 2],
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, color in zip(box["boxes"], colors, strict=True):
            patch.set_facecolor(color)
            patch.set_alpha(0.35)
        for pos, values, color in zip([1, 2], groups, colors, strict=True):
            jitter = rng.uniform(-0.12, 0.12, len(values))
            ax.scatter(
                pos + jitter,
                values,
                s=28,
                color=color,
                alpha=0.85,
                edgecolor="white",
                linewidth=0.4,
                zorder=3,
            )
        ax.set_xticks([1, 2], labels)
        ax.set_title(title)
        ax.set_ylabel("CLDN4 log2(FPKM + 1)")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    fig.savefig(output.with_suffix(".svg"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, default=Path(".cache"))
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    paths = fetch_sources(args.cache_dir)
    manifest = {
        "dataset_revision": REVISION,
        "citation": (
            "Nathanson T et al. Cancer Immunol Res. 2017;5:84-91. "
            "doi:10.1158/2326-6066.CIR-16-0019"
        ),
        "sources": SOURCES,
        "primary_gene_id": CANONICAL_CLDN4,
        "bootstrap_seed": SEED,
        "bootstrap_replicates": BOOTSTRAP_DRAWS,
    }
    (args.output_dir / "download_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    rows = load_cldn4_rows(paths["cufflinks.tar.gz"])
    extra_id_max = float(
        rows.loc[rows["tracking_id"].ne(CANONICAL_CLDN4), "FPKM"].max()
    )
    cohort = pd.read_csv(paths["cohort.csv"], dtype={"sample": str})
    expression = collapse_cldn4(rows, "canonical")
    expression["sample"] = normalize_sample_id(expression["sample"])
    cohort["sample"] = normalize_sample_id(cohort["sample"])
    samples = expression.merge(cohort, on="sample", how="left", validate="one_to_one")
    if samples["Benefit"].isna().any() or len(samples) != 24:
        raise ValueError("Expression samples failed to match the 24-sample RNA cohort")
    samples["Benefit"] = samples["Benefit"].astype(bool)
    samples["log2_FPKM_plus_1"] = np.log2(samples["FPKM"] + 1)

    keep = [
        "Study ID",
        "sample",
        "MSK_ID",
        "Biopsy pre or post ipi",
        "Benefit",
        "OS",
        "Alive",
        "cldn4_id",
        "FPKM",
        "log2_FPKM_plus_1",
    ]
    samples[keep].sort_values(
        ["Biopsy pre or post ipi", "Benefit", "sample"]
    ).to_csv(args.output_dir / "cldn4_samples.tsv", sep="\t", index=False)

    pre = samples.loc[samples["Biopsy pre or post ipi"].eq("pre")].copy()
    post = samples.loc[samples["Biopsy pre or post ipi"].eq("post")].copy()
    if len(pre) != 9 or len(post) != 15:
        raise ValueError("Unexpected pretreatment/post-treatment RNA counts")

    median_collapsed = collapse_cldn4(rows, "symbol_median")
    median_collapsed["sample"] = normalize_sample_id(median_collapsed["sample"])
    median_samples = median_collapsed.merge(
        cohort, on="sample", how="left", validate="one_to_one"
    )
    median_samples["Benefit"] = median_samples["Benefit"].astype(bool)
    median_samples["log2_FPKM_plus_1"] = np.log2(median_samples["FPKM"] + 1)
    median_pre = median_samples.loc[
        median_samples["Biopsy pre or post ipi"].eq("pre")
    ]

    statistics = [
        summarize(pre, "pretreatment_primary"),
        summarize(post, "posttreatment_exploratory"),
        summarize(median_pre, "pretreatment_symbol_median_sensitivity"),
    ]
    pd.DataFrame(statistics).to_csv(
        args.output_dir / "statistics.tsv", sep="\t", index=False
    )
    loo = leave_one_out(pre)
    loo.to_csv(args.output_dir / "leave_one_out.tsv", sep="\t", index=False)

    pre_alive_by_benefit = {
        "benefit_alive": int(pre.loc[pre["Benefit"], "Alive"].sum()),
        "benefit_n": int(pre["Benefit"].sum()),
        "no_benefit_alive": int(pre.loc[~pre["Benefit"], "Alive"].sum()),
        "no_benefit_n": int((~pre["Benefit"]).sum()),
    }
    primary = statistics[0]
    honest = {
        "label": "OPEN_MATRIX_ANALYZED_UNDERPOWERED",
        "primary_result": (
            "Pretreatment CLDN4 is directionally higher with deposited clinical "
            f"benefit (n=4 vs 5; exact p={primary['exact_two_sided_p']:.4f}; "
            f"AUC={primary['auc_higher_predicts_benefit']:.2f}, bootstrap 95% CI "
            f"{primary['auc_bootstrap_95ci_low']:.2f}-"
            f"{primary['auc_bootstrap_95ci_high']:.2f}), but the nine-sample "
            "result is not conclusive and is leave-one-out fragile."
        ),
        "ucla_substitution_used": False,
        "leftover_note": (
            "This is the leftover public melanoma ICI RNA matrix after the Riaz "
            "and Liu B5 analyses. No UCLA/Hugo matrix was substituted."
        ),
        "os_not_independent_of_benefit": (
            pre_alive_by_benefit["benefit_alive"]
            == pre_alive_by_benefit["benefit_n"]
            and pre_alive_by_benefit["no_benefit_alive"] == 0
        ),
        "second_cldn4_id_max_fpkm": extra_id_max,
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(
            {
                **manifest,
                "honest_verdict": honest,
                "vital_status_alignment": pre_alive_by_benefit,
                "statistics": statistics,
                "leave_one_out_p_range": [
                    float(loo["exact_two_sided_p"].min()),
                    float(loo["exact_two_sided_p"].max()),
                ],
            },
            indent=2,
        )
        + "\n"
    )
    plot_panels(samples, args.output_dir / "cldn4_vs_benefit.png")


if __name__ == "__main__":
    main()
