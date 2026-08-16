#!/usr/bin/env python3
"""Reproduce the CLDN4 check in the Nathanson anti-CTLA-4 melanoma cohort."""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

matplotlib.use("Agg")
import matplotlib.pyplot as plt


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


def load_cldn4(archive: Path) -> pd.DataFrame:
    columns = ["tracking_id", "gene_short_name", "FPKM", "sample", "MSK_ID"]
    with tarfile.open(archive, "r:gz") as tar:
        member = tar.getmember("cufflinks.csv")
        handle = tar.extractfile(member)
        if handle is None:
            raise FileNotFoundError("cufflinks.csv is missing from archive")
        chunks = []
        for chunk in pd.read_csv(handle, usecols=columns, chunksize=250_000):
            chunks.append(chunk.loc[chunk["gene_short_name"].eq("CLDN4")])
    cldn4 = pd.concat(chunks, ignore_index=True)
    if cldn4.empty:
        raise ValueError("CLDN4 is absent from the expression data")
    # Follow the paper: collapse duplicate official gene symbols by median FPKM.
    return cldn4.groupby("sample", as_index=False).agg(
        FPKM=("FPKM", "median"),
        MSK_ID=("MSK_ID", "first"),
        n_gene_ids=("tracking_id", "nunique"),
    )


def analyze(samples: pd.DataFrame, label: str) -> dict:
    benefit = samples.loc[samples["Benefit"], "log2_FPKM_plus_1"]
    no_benefit = samples.loc[~samples["Benefit"], "log2_FPKM_plus_1"]
    test = mannwhitneyu(benefit, no_benefit, alternative="two-sided", method="exact")
    n_product = len(benefit) * len(no_benefit)
    return {
        "analysis": label,
        "n_benefit": len(benefit),
        "n_no_benefit": len(no_benefit),
        "median_FPKM_benefit": float(samples.loc[samples["Benefit"], "FPKM"].median()),
        "median_FPKM_no_benefit": float(
            samples.loc[~samples["Benefit"], "FPKM"].median()
        ),
        "mann_whitney_U": float(test.statistic),
        "exact_two_sided_p": float(test.pvalue),
        "auc_higher_predicts_benefit": float(test.statistic / n_product),
        "rank_biserial": float(2 * test.statistic / n_product - 1),
    }


def plot_pre_treatment(samples: pd.DataFrame, output: Path) -> None:
    pre = samples.loc[samples["Biopsy pre or post ipi"].eq("pre")].copy()
    rng = np.random.default_rng(20260816)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for position, (benefit, color) in enumerate(
        [(False, "#9ca3af"), (True, "#2563eb")]
    ):
        values = pre.loc[pre["Benefit"].eq(benefit), "log2_FPKM_plus_1"].to_numpy()
        x = position + rng.uniform(-0.07, 0.07, size=len(values))
        ax.scatter(x, values, s=46, color=color, edgecolor="white", linewidth=0.7)
        median = np.median(values)
        ax.plot([position - 0.18, position + 0.18], [median, median], color="black")
    ax.set_xticks([0, 1], ["No benefit\n(n=5)", "Benefit\n(n=4)"])
    ax.set_ylabel("CLDN4 log2(FPKM + 1)")
    ax.set_title("Nathanson 2017, pretreatment biopsies")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, default=Path(".cache"))
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    paths = fetch_sources(args.cache_dir)
    with (args.output_dir / "download_manifest.json").open("w") as handle:
        json.dump(
            {"dataset_revision": REVISION, "sources": SOURCES},
            handle,
            indent=2,
        )
        handle.write("\n")
    expression = load_cldn4(paths["cufflinks.tar.gz"])
    cohort = pd.read_csv(paths["cohort.csv"], dtype={"sample": str})
    expression["sample"] = expression["sample"].astype(str).str.lstrip("0")
    cohort["sample"] = cohort["sample"].str.lstrip("0")
    samples = expression.merge(cohort, on="sample", how="left", validate="one_to_one")
    if samples["Benefit"].isna().any():
        raise ValueError("Expression samples failed to match cohort metadata")
    samples["Benefit"] = samples["Benefit"].astype(bool)
    samples["log2_FPKM_plus_1"] = np.log2(samples["FPKM"] + 1)

    keep = [
        "Study ID",
        "sample",
        "MSK_ID",
        "Biopsy pre or post ipi",
        "Benefit",
        "FPKM",
        "log2_FPKM_plus_1",
        "n_gene_ids",
    ]
    samples[keep].sort_values(["Biopsy pre or post ipi", "Benefit", "sample"]).to_csv(
        args.output_dir / "cldn4_samples.tsv", sep="\t", index=False
    )

    statistics = [
        analyze(
            samples.loc[samples["Biopsy pre or post ipi"].eq("pre")],
            "pretreatment_primary",
        ),
        analyze(
            samples.loc[samples["Biopsy pre or post ipi"].eq("post")],
            "posttreatment_descriptive",
        ),
    ]
    pd.DataFrame(statistics).to_csv(
        args.output_dir / "statistics.tsv", sep="\t", index=False
    )
    with (args.output_dir / "summary.json").open("w") as handle:
        json.dump(
            {
                "dataset_revision": REVISION,
                "sources": SOURCES,
                "honest_verdict": {
                    "label": "OPEN_MATRIX_ANALYZED_UNDERPOWERED",
                    "primary_result": (
                        "CLDN4 is directionally higher with pretreatment clinical "
                        "benefit, but the exact two-sided result is not conclusive "
                        "(n=9, p=0.0635)."
                    ),
                    "ucla_substitution_used": False,
                },
                "statistics": statistics,
            },
            handle,
            indent=2,
        )
        handle.write("\n")
    plot_pre_treatment(samples, args.output_dir / "cldn4_pretreatment.png")


if __name__ == "__main__":
    main()
