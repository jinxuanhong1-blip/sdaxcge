#!/usr/bin/env python3
"""Analyze TACSTD2 in malignant versus T/NK cells in GSE205335."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import tempfile
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import rdata
from scipy import sparse, stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt


GENES = ("TACSTD2", "EPCAM", "PTPRC")
RESPONSE_MAP = {"PR": "R", "CR": "R", "SD": "NR", "PD": "NR", "NE": "NE"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0]
                    current["title"] = titles[0]
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith(
                "!Sample_characteristics_ch1 = "
            ):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0]
            current["title"] = titles[0]
            records.append(current)

    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False)
        + "-"
        + read_end
        + "P"
    )
    return metadata.rename(
        columns={
            "tumor stage": "tumor_stage",
            "cancer subtype": "cancer_subtype",
        }
    )


def load_selected_genes(
    path: Path,
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(matrix_path)

    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix(
        (obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False
    )
    extracted: dict[str, np.ndarray] = {}
    for gene in GENES:
        positions = np.flatnonzero(genes == gene)
        if len(positions) != 1:
            raise ValueError(f"Expected exactly one {gene} feature")
        extracted[gene] = np.asarray(
            matrix.getrow(int(positions[0])).toarray()
        ).ravel()
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    return extracted, library_umi, barcodes


def summarize(group: pd.DataFrame) -> pd.Series:
    return pd.Series(
        {
            "n_cells": len(group),
            "pct_pos": 100.0 * group["tacstd2_umi"].gt(0).mean(),
            "mean_log1p_cp10k": group["tacstd2_log1p_cp10k"].mean(),
            "pseudobulk_cpm": (
                1e6 * group["tacstd2_umi"].sum() / group["total_umi"].sum()
            ),
        }
    )


def fmt_p(value: float) -> str:
    return f"{value:.2e}" if value < 0.001 else f"{value:.4f}"


def analyze(
    cells: pd.DataFrame,
    metadata: pd.DataFrame,
    outdir: Path,
    min_cells: int,
) -> dict[str, object]:
    cells = cells.merge(
        metadata[
            [
                "orig.ident",
                "gsm",
                "patient",
                "tissue",
                "recist",
                "platform",
                "cancer_subtype",
            ]
        ],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["gsm"].isna().any():
        raise ValueError("Some identity-table samples did not match GEO metadata")
    cells["response"] = cells["recist"].map(RESPONSE_MAP)
    cells["tacstd2_log1p_cp10k"] = np.log1p(
        1e4 * cells["tacstd2_umi"] / cells["total_umi"]
    )
    cells["compartment"] = np.select(
        [
            cells["lineage.sub"].eq("Malignant cells"),
            cells["lineage.total"].eq("T/NK cells"),
        ],
        ["Malignant", "T/NK"],
        default="Other",
    )
    subset = cells[cells["compartment"].isin(["Malignant", "T/NK"])]

    sample_keys = [
        "patient",
        "orig.ident",
        "gsm",
        "tissue",
        "recist",
        "response",
        "platform",
        "cancer_subtype",
        "compartment",
    ]
    per_sample = (
        subset.groupby(sample_keys, observed=True)
        .apply(summarize, include_groups=False)
        .reset_index()
    )
    per_patient = (
        subset.groupby(
            ["patient", "recist", "response", "compartment"], observed=True
        )
        .apply(summarize, include_groups=False)
        .reset_index()
    )
    per_sample.to_csv(outdir / "per_sample_tacstd2.csv", index=False)
    per_patient.to_csv(outdir / "per_patient_tacstd2.csv", index=False)

    malignant = per_patient[
        (per_patient["compartment"] == "Malignant")
        & (per_patient["n_cells"] >= min_cells)
    ].copy()
    tnk = per_patient[
        (per_patient["compartment"] == "T/NK")
        & (per_patient["n_cells"] >= min_cells)
    ].copy()
    paired = malignant.merge(tnk, on="patient", suffixes=("_mal", "_tnk"))

    all_patients = metadata[["patient", "recist"]].drop_duplicates()
    malignant_counts = per_patient[
        per_patient["compartment"] == "Malignant"
    ][["patient", "n_cells"]]
    patient_counts = all_patients.merge(
        malignant_counts, on="patient", how="left"
    ).fillna({"n_cells": 0})
    excluded = patient_counts[patient_counts["n_cells"] < min_cells]

    lines = [
        "GSE205335 TACSTD2 (TROP2) analysis — statistical results",
        "=" * 60,
        (
            f"Cells total: {len(cells)}; "
            f"malignant: {cells['compartment'].eq('Malignant').sum()}; "
            f"T/NK: {cells['compartment'].eq('T/NK').sum()}"
        ),
        (
            f"Patients with >= {min_cells} malignant cells: {len(malignant)} "
            f"(R={malignant['response'].eq('R').sum()}, "
            f"NR={malignant['response'].eq('NR').sum()}, "
            f"NE={malignant['response'].eq('NE').sum()})"
        ),
    ]
    if len(excluded):
        lines.append(
            "Patients excluded (< min malignant cells captured): "
            + ", ".join(
                f"{row.patient}({int(row.n_cells)} cells, RECIST {row.recist})"
                for row in excluded.itertuples()
            )
        )
        lines.append(
            "    NOTE: dropout is asymmetric — "
            f"{excluded['recist'].isin(['PR', 'CR']).sum()} responder vs "
            f"{excluded['recist'].isin(['SD', 'PD']).sum()} non-responder "
            "patients had (near-)zero malignant cells captured."
        )

    lines.extend(
        [
            "",
            "[1] Primary: TACSTD2 in malignant vs T/NK cells, paired within "
            "patient (all patients incl. NE)",
            f"    Patients require >= {min_cells} cells in each compartment.",
        ]
    )
    paired_results: dict[str, dict[str, float | int]] = {}
    for metric in ["mean_log1p_cp10k", "pct_pos", "pseudobulk_cpm"]:
        malignant_values = paired[f"{metric}_mal"].to_numpy()
        tnk_values = paired[f"{metric}_tnk"].to_numpy()
        result = stats.wilcoxon(
            malignant_values, tnk_values, alternative="two-sided"
        )
        paired_results[metric] = {
            "n_pairs": len(paired),
            "malignant_median": float(np.median(malignant_values)),
            "tnk_median": float(np.median(tnk_values)),
            "wilcoxon_w": float(result.statistic),
            "p_value": float(result.pvalue),
        }
        lines.append(
            f"    {metric}: n={len(paired)} pairs | malignant "
            f"median={np.median(malignant_values):.3f} | T/NK "
            f"median={np.median(tnk_values):.3f} | Wilcoxon signed-rank "
            f"W={result.statistic:.1f}, p={fmt_p(result.pvalue)}"
        )

    lines.extend(
        [
            "",
            "[2] Secondary: malignant-cell TACSTD2 in responders (PR) vs "
            "non-responders (SD/PD)",
            "    Unit = patient; NE patients excluded.",
        ]
    )
    response_results: dict[str, dict[str, float | int]] = {}
    evaluable = malignant[malignant["response"].isin(["R", "NR"])]
    for metric in ["mean_log1p_cp10k", "pct_pos", "pseudobulk_cpm"]:
        responders = evaluable.loc[evaluable["response"] == "R", metric]
        nonresponders = evaluable.loc[evaluable["response"] == "NR", metric]
        result = stats.mannwhitneyu(
            responders, nonresponders, alternative="two-sided"
        )
        response_results[metric] = {
            "n_r": len(responders),
            "n_nr": len(nonresponders),
            "r_median": float(np.median(responders)),
            "nr_median": float(np.median(nonresponders)),
            "mannwhitney_u": float(result.statistic),
            "p_value": float(result.pvalue),
        }
        lines.append(
            f"    {metric}: R n={len(responders)} "
            f"median={np.median(responders):.3f} | NR n={len(nonresponders)} "
            f"median={np.median(nonresponders):.3f} | Mann-Whitney "
            f"U={result.statistic:.1f}, p={fmt_p(result.pvalue)}"
        )

    malignant_cells = cells[cells["compartment"] == "Malignant"]
    tnk_cells = cells[cells["compartment"] == "T/NK"]
    marker_sanity = {
        "epcam_positive_malignant_pct": float(
            100 * malignant_cells["epcam_umi"].gt(0).mean()
        ),
        "epcam_positive_tnk_pct": float(
            100 * tnk_cells["epcam_umi"].gt(0).mean()
        ),
        "ptprc_positive_malignant_pct": float(
            100 * malignant_cells["ptprc_umi"].gt(0).mean()
        ),
        "ptprc_positive_tnk_pct": float(
            100 * tnk_cells["ptprc_umi"].gt(0).mean()
        ),
    }
    lines.extend(
        [
            "",
            "[sanity] EPCAM+ fraction: malignant "
            f"{marker_sanity['epcam_positive_malignant_pct']:.1f}% vs T/NK "
            f"{marker_sanity['epcam_positive_tnk_pct']:.1f}%; PTPRC(CD45)+ "
            "fraction: malignant "
            f"{marker_sanity['ptprc_positive_malignant_pct']:.1f}% vs T/NK "
            f"{marker_sanity['ptprc_positive_tnk_pct']:.1f}%",
        ]
    )
    (outdir / "stats_results.txt").write_text("\n".join(lines) + "\n")

    colors = {"R": "#2166ac", "NR": "#b2182b", "NE": "#999999"}
    order = ["R", "NR", "NE"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for panel, metric, ylabel, title in [
        (
            axes[0],
            "mean_log1p_cp10k",
            "TACSTD2 mean log1p(CP10K), malignant",
            "Malignant TACSTD2 by RECIST response",
        ),
        (
            axes[1],
            "pct_pos",
            "% TACSTD2+ malignant cells",
            "TACSTD2+ malignant-cell fraction",
        ),
    ]:
        for index, response in enumerate(order):
            values = malignant.loc[malignant["response"] == response, metric]
            jitter = np.random.default_rng(index).normal(index, 0.06, len(values))
            panel.scatter(
                jitter, values, color=colors[response], s=60, alpha=0.85
            )
            if len(values):
                panel.hlines(
                    np.median(values), index - 0.25, index + 0.25, color="black"
                )
        panel.set_xticks(range(3))
        panel.set_xticklabels(
            [
                f"{response}\n(n={malignant['response'].eq(response).sum()})"
                for response in order
            ]
        )
        panel.set_ylabel(ylabel)
        panel.set_title(title)

    ax = axes[2]
    for row in paired.itertuples():
        ax.plot(
            [0, 1],
            [row.mean_log1p_cp10k_mal, row.mean_log1p_cp10k_tnk],
            color="grey",
            alpha=0.5,
        )
    ax.scatter(
        np.zeros(len(paired)),
        paired["mean_log1p_cp10k_mal"],
        color="#762a83",
        label="Malignant",
    )
    ax.scatter(
        np.ones(len(paired)),
        paired["mean_log1p_cp10k_tnk"],
        color="#1b7837",
        label="T/NK",
    )
    ax.set_xticks([0, 1], ["Malignant", "T/NK"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("TACSTD2 mean log1p(CP10K)")
    ax.set_title(f"Paired within patient (n={len(paired)})")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "tacstd2_summary.png", dpi=200)
    plt.close(fig)

    ordered = malignant.sort_values(
        ["response", "mean_log1p_cp10k"], ascending=[True, False]
    )
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(
        range(len(ordered)),
        ordered["mean_log1p_cp10k"],
        color=[colors[value] for value in ordered["response"]],
    )
    ax.set_xticks(
        range(len(ordered)),
        [
            f"{patient}\n{recist}"
            for patient, recist in zip(
                ordered["patient"], ordered["recist"], strict=True
            )
        ],
        fontsize=8,
    )
    ax.set_ylabel("TACSTD2 mean log1p(CP10K), malignant")
    ax.set_title("Malignant TACSTD2 per patient")
    fig.tight_layout()
    fig.savefig(outdir / "tacstd2_per_patient_bars.png", dpi=200)
    plt.close(fig)

    return {
        "dataset": "GSE205335",
        "primary_contrast": "malignant_vs_tnk_paired_within_patient",
        "min_cells_per_compartment": min_cells,
        "n_total_cells": len(cells),
        "n_malignant_cells": int(cells["compartment"].eq("Malignant").sum()),
        "n_tnk_cells": int(cells["compartment"].eq("T/NK").sum()),
        "n_evaluable_patient_pairs": len(paired),
        "paired_results": paired_results,
        "secondary_response_results": response_results,
        "marker_sanity": marker_sanity,
        "raw_data": "not accessed; controlled EGA raw intentionally skipped",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--identities", required=True, type=Path)
    parser.add_argument("--soft", required=True, type=Path)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/w200/A3_GSE205335"),
    )
    parser.add_argument("--min-cells", type=int, default=20)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    identities = pd.read_csv(args.identities, sep="\t")
    if identities["barcode"].duplicated().any():
        raise ValueError("Cell-identity barcodes are not unique")
    selected, library_umi, barcodes = load_selected_genes(args.matrix)
    indexed = identities.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise ValueError(
            f"Matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra"
        )
    cells = indexed.loc[barcodes].reset_index()
    cells["total_umi"] = library_umi
    for gene in GENES:
        cells[f"{gene.lower()}_umi"] = selected[gene]

    metadata = parse_geo_soft(args.soft)
    metadata.to_csv(args.outdir / "gsm_sample_metadata.csv", index=False)
    summary = analyze(cells, metadata, args.outdir, args.min_cells)
    (args.outdir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    provenance = pd.DataFrame(
        [
            {
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in [args.identities, args.matrix, args.soft]
        ]
    )
    provenance.to_csv(args.outdir / "provenance.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
