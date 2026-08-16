#!/usr/bin/env python3
"""Reproduce the gene-set overlaps and TACSTD2-only LuCA Census summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import cellxgene_census
import numpy as np
import pandas as pd


CENSUS_VERSION = "2025-11-08"
DATASET_ID = "1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1"
RELEVANT_CELL_TYPES = [
    "pulmonary alveolar type 1 cell",
    "club cell",
    "malignant cell",
    "epithelial cell of lung",
    "pulmonary alveolar type 2 cell",
    "neutrophil",
    "macrophage",
    "plasma cell",
    "CD8-positive, alpha-beta T cell",
    "CD4-positive, alpha-beta T cell",
    "B cell",
    "alveolar macrophage",
]


def summarize(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    return (
        frame.groupby(columns, observed=True)
        .agg(
            n_cells=("expression", "size"),
            detected_cells=("detected", "sum"),
            detected_fraction=("detected", "mean"),
            mean_census_normalized_expression=("expression", "mean"),
            mean_when_detected=(
                "expression",
                lambda values: (
                    values[values > 0].mean() if (values > 0).any() else 0.0
                ),
            ),
        )
        .reset_index()
    )


def query_tacstd2() -> pd.DataFrame:
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            measurement_name="RNA",
            X_name="normalized",
            obs_value_filter=f"dataset_id == '{DATASET_ID}'",
            obs_column_names=["cell_type", "assay"],
            var_value_filter="feature_name == 'TACSTD2'",
            var_column_names=["feature_name", "feature_id"],
        )

    if adata.shape != (1_283_972, 1):
        raise RuntimeError(f"Unexpected query shape: {adata.shape}")

    frame = adata.obs[["cell_type", "assay"]].copy()
    frame["expression"] = np.asarray(adata.X.toarray()).ravel()
    frame["detected"] = frame["expression"] > 0
    return frame


def write_expression_summaries(frame: pd.DataFrame, output_dir: Path) -> None:
    by_type = summarize(frame, ["cell_type"])
    by_type = (
        by_type[by_type["cell_type"].isin(RELEVANT_CELL_TYPES)]
        .assign(
            cell_type=lambda x: pd.Categorical(
                x["cell_type"], categories=RELEVANT_CELL_TYPES, ordered=True
            )
        )
        .sort_values("cell_type")
    )
    by_type.to_csv(output_dir / "tacstd2_relevant_cell_types.tsv", sep="\t", index=False)

    by_assay = summarize(
        frame[frame["cell_type"].isin(["malignant cell", "neutrophil"])],
        ["cell_type", "assay"],
    ).sort_values(
        ["cell_type", "n_cells"], ascending=[True, False]
    )
    by_assay.to_csv(
        output_dir / "tacstd2_malignant_vs_neutrophil_by_assay.tsv",
        sep="\t",
        index=False,
    )


def write_overlaps(signature_path: Path, output_dir: Path) -> None:
    table = pd.read_csv(signature_path, sep="\t")
    sets = {
        row.signature: set(row.genes.split(","))
        for row in table.itertuples(index=False)
    }
    sets["TACSTD2"] = {"TACSTD2"}
    comparisons = [
        ("TACSTD2", "Salcher_TRN"),
        ("TACSTD2", "Salcher_TAN"),
        ("TACSTD2", "Salcher_NAN"),
        ("TACSTD2", "Salcher_major_neutrophils"),
        ("TACSTD2", "Leader_LCAM_hi"),
        ("TACSTD2", "Leader_LCAM_lo"),
        ("Salcher_TAN", "Salcher_NAN"),
        ("Salcher_TRN", "Salcher_major_neutrophils"),
        ("Salcher_TRN", "Leader_LCAM_hi"),
        ("Salcher_TRN", "Leader_LCAM_lo"),
        ("Salcher_major_neutrophils", "Leader_LCAM_hi"),
        ("Salcher_major_neutrophils", "Leader_LCAM_lo"),
        ("Leader_LCAM_hi", "Leader_LCAM_lo"),
    ]
    rows = []
    for left, right in comparisons:
        shared = sorted(sets[left] & sets[right])
        rows.append(
            {
                "set_a": left,
                "set_b": right,
                "n_a": len(sets[left]),
                "n_b": len(sets[right]),
                "n_intersection": len(shared),
                "intersection": ",".join(shared),
            }
        )
    pd.DataFrame(rows).to_csv(
        output_dir / "signature_overlap.tsv", sep="\t", index=False
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).resolve().parent
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    signature_path = Path(__file__).resolve().with_name("signatures.tsv")
    write_overlaps(signature_path, args.output_dir)
    write_expression_summaries(query_tacstd2(), args.output_dir)


if __name__ == "__main__":
    main()
