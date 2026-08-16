#!/usr/bin/env python3
"""Targeted patient-level extraction from the GSE241934 sparse scRNA matrix."""

from __future__ import annotations

import gzip
from pathlib import Path
import subprocess
import tempfile
import urllib.request

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "ici"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl"
FILES = {
    "matrix": "GSE241934_Real_Matrix.mtx.gz",
    "metadata": "GSE241934_Real_Meta.txt.gz",
    "features": "GSE241934_RWC_features.tsv.gz",
    "barcodes": "GSE241934_RWC_barcodes.tsv.gz",
}
MIN_EPITHELIAL_CELLS = 10


def download(name: str, cache: Path) -> Path:
    path = cache / FILES[name]
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".part")
        urllib.request.urlretrieve(f"{BASE}/{FILES[name]}", temporary)
        temporary.replace(path)
    return path


def target_rows(features_path: Path) -> dict[int, str]:
    features = pd.read_csv(features_path, sep="\t", header=None)
    output = {}
    for gene in ("TACSTD2", "CLDN4"):
        matches = np.flatnonzero(features[1].eq(gene).to_numpy())
        if len(matches) != 1:
            raise RuntimeError(f"Expected one {gene} feature, found {len(matches)}")
        output[int(matches[0] + 1)] = gene  # MatrixMarket is one-based.
    return output


def extract_coordinates(matrix: Path, rows: dict[int, str], destination: Path) -> None:
    if destination.exists():
        return
    expression = " || ".join(f"$1=={row}" for row in rows)
    with destination.open("wb") as output:
        unzip = subprocess.Popen(["gzip", "-dc", matrix], stdout=subprocess.PIPE)
        assert unzip.stdout is not None
        select = subprocess.run(
            ["awk", expression], stdin=unzip.stdout, stdout=output, check=True
        )
        unzip.stdout.close()
        return_code = unzip.wait()
        if return_code != 0 or select.returncode != 0:
            destination.unlink(missing_ok=True)
            raise RuntimeError("Sparse target-row extraction failed")


def bh_adjust(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.empty(len(values))
    adjusted[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return pd.Series(np.minimum(adjusted, 1), index=p_values.index)


def main() -> None:
    cache = Path(tempfile.gettempdir()) / "ici_gse241934"
    paths = {name: download(name, cache) for name in FILES}
    metadata = pd.read_csv(paths["metadata"], sep="\t", low_memory=False)
    barcodes = pd.read_csv(paths["barcodes"], header=None)[0].astype(str)
    if not np.array_equal(metadata["cellID"].astype(str), barcodes):
        raise RuntimeError("Barcode order does not match cell metadata")

    rows = target_rows(paths["features"])
    coordinates_path = cache / "TACSTD2_CLDN4_coordinates.tsv"
    extract_coordinates(paths["matrix"], rows, coordinates_path)
    coordinates = pd.read_csv(
        coordinates_path,
        sep=r"\s+",
        header=None,
        names=["row", "column", "value"],
    )
    values = {gene: np.zeros(len(metadata), dtype=np.int32) for gene in rows.values()}
    for row, gene in rows.items():
        selected = coordinates[coordinates["row"].eq(row)]
        values[gene][selected["column"].to_numpy(int) - 1] = selected[
            "value"
        ].to_numpy(int)

    patient_rows = []
    for patient, cells in metadata.groupby("sampleID", sort=True):
        epithelial = cells["major_cell_type"].eq("Epi")
        indices = cells.index[epithelial]
        non_epithelial = ~epithelial
        record = {
            "patient": patient,
            "pathologic_response": cells["Pathological Response"].iloc[0],
            "mpr_or_pcr": cells["Pathological Response"].iloc[0] in {"MPR", "pCR"},
            "histology": cells["Histology"].iloc[0],
            "pd1_agent": cells["PD1"].iloc[0],
            "n_cells": len(cells),
            "n_epithelial": int(epithelial.sum()),
            "tnk_fraction_all_cells": cells["major_cell_type"].isin(["T", "NK"]).mean(),
            "tnk_fraction_non_epithelial": cells.loc[
                non_epithelial, "major_cell_type"
            ].isin(["T", "NK"]).mean(),
        }
        for gene in ("TACSTD2", "CLDN4"):
            gene_values = values[gene][indices]
            record[f"{gene}_epithelial_mean_log1p"] = (
                float(np.log1p(gene_values).mean()) if len(gene_values) else np.nan
            )
            record[f"{gene}_epithelial_detection_fraction"] = (
                float((gene_values > 0).mean()) if len(gene_values) else np.nan
            )
        patient_rows.append(record)
    scores = pd.DataFrame(patient_rows)
    scores["included_min10_epithelial"] = scores["n_epithelial"] >= MIN_EPITHELIAL_CELLS
    scores.to_csv(OUT / "gse241934_patient_scores.tsv", sep="\t", index=False)
    analysis = scores[scores["included_min10_epithelial"]].copy()

    statistics = []
    for gene in ("TACSTD2", "CLDN4"):
        column = f"{gene}_epithelial_mean_log1p"
        nmpr = analysis.loc[~analysis["mpr_or_pcr"], column]
        benefit = analysis.loc[analysis["mpr_or_pcr"], column]
        _, p_value = mannwhitneyu(nmpr, benefit, alternative="two-sided")
        statistics.append(
            {
                "analysis": "epithelial score; NMPR vs MPR/pCR",
                "feature": gene,
                "n": len(analysis),
                "n_nmpr": len(nmpr),
                "n_mpr_pcr": len(benefit),
                "effect": nmpr.median() - benefit.median(),
                "effect_definition": "median(NMPR)-median(MPR/pCR)",
                "rho": np.nan,
                "p": p_value,
            }
        )
        for outcome in ("tnk_fraction_all_cells", "tnk_fraction_non_epithelial"):
            rho, p_value = spearmanr(analysis[column], analysis[outcome])
            statistics.append(
                {
                    "analysis": f"epithelial {gene} vs {outcome}",
                    "feature": gene,
                    "n": len(analysis),
                    "n_nmpr": np.nan,
                    "n_mpr_pcr": np.nan,
                    "effect": np.nan,
                    "effect_definition": "Spearman rho",
                    "rho": rho,
                    "p": p_value,
                }
            )
    rho, p_value = spearmanr(
        analysis["TACSTD2_epithelial_mean_log1p"],
        analysis["CLDN4_epithelial_mean_log1p"],
    )
    statistics.append(
        {
            "analysis": "epithelial TACSTD2 vs CLDN4",
            "feature": "TACSTD2~CLDN4",
            "n": len(analysis),
            "n_nmpr": np.nan,
            "n_mpr_pcr": np.nan,
            "effect": np.nan,
            "effect_definition": "Spearman rho",
            "rho": rho,
            "p": p_value,
        }
    )
    statistics = pd.DataFrame(statistics)
    statistics["q_bh"] = bh_adjust(statistics["p"])
    statistics.to_csv(OUT / "gse241934_statistics.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    rng = np.random.default_rng(20260816)
    for ax, gene in zip(axes, ("TACSTD2", "CLDN4")):
        column = f"{gene}_epithelial_mean_log1p"
        groups = [
            analysis.loc[~analysis["mpr_or_pcr"], column],
            analysis.loc[analysis["mpr_or_pcr"], column],
        ]
        ax.boxplot(
            groups,
            tick_labels=[f"NMPR (n={len(groups[0])})", f"MPR/pCR (n={len(groups[1])})"],
        )
        for position, group in enumerate(groups, 1):
            ax.scatter(
                position + rng.uniform(-0.06, 0.06, len(group)),
                group,
                color="#276FBF",
            )
        result = statistics[
            statistics["analysis"].eq("epithelial score; NMPR vs MPR/pCR")
            & statistics["feature"].eq(gene)
        ].iloc[0]
        ax.set(
            title=f"GSE241934: {gene}\nNMPR−MPR={result.effect:.3f}, p={result.p:.3g}",
            ylabel="Epithelial-cell mean log1p UMI",
        )
    fig.tight_layout()
    fig.savefig(OUT / "gse241934_response.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
