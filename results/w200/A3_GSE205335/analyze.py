#!/usr/bin/env python3
"""Sample-level TACSTD2/T-NK analysis for GEO GSE205335.

The GEO matrix is an R Matrix::dgCMatrix serialized as an RDS.  GEO added an
outer gzip layer to a file that was already gzip-compressed; if the input ends
in .gz, this script removes one layer into a temporary file before parsing.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import tempfile
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdata
import seaborn as sns
from scipy import sparse, stats


ENDPOINTS = {
    "tnk_fraction_nonmalignant": "T/NK / non-malignant cells",
    "tnk_fraction_all": "T/NK / all cells",
    "cd4_fraction_tnk": "CD4 T / T/NK",
    "cd8_fraction_tnk": "CD8 T / T/NK",
    "nk_fraction_tnk": "NK / T/NK",
}
EXPOSURES = {
    "tacstd2_pseudobulk_cpm": "malignant TACSTD2 pseudobulk CPM",
    "tacstd2_detection_fraction": "TACSTD2+ malignant-cell fraction",
}


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
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0]
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions = []
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


def unwrap_geo_gzip(path: Path, tempdir: Path) -> Path:
    if path.suffix != ".gz":
        return path
    unwrapped = tempdir / path.stem
    with gzip.open(path, "rb") as source, unwrapped.open("wb") as destination:
        shutil.copyfileobj(source, destination, 16 * 1024 * 1024)
    return unwrapped


def load_matrix(path: Path) -> tuple[sparse.csc_matrix, np.ndarray, np.ndarray]:
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = unwrap_geo_gzip(path, Path(tmp))
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
    return matrix, genes, barcodes


def bh_adjust(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return pd.Series(result, index=p_values.index)


def bootstrap_spearman(
    x: np.ndarray, y: np.ndarray, seed: int, iterations: int = 5000
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    estimates: list[float] = []
    for _ in range(iterations):
        selected = rng.integers(0, n, n)
        if np.unique(x[selected]).size < 2 or np.unique(y[selected]).size < 2:
            continue
        estimates.append(stats.spearmanr(x[selected], y[selected]).statistic)
    if not estimates:
        return np.nan, np.nan
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def correlation_table(metrics: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    eligible = metrics[
        (metrics["n_malignant"] >= min_cells) & (metrics["n_tnk"] >= min_cells)
    ].copy()
    cohorts = {
        "all_eligible_samples": eligible,
        "adc_samples": eligible[eligible["cancer_subtype"] == "ADC"],
        "patient_means": eligible.groupby("patient", as_index=False).agg(
            {
                **{name: "mean" for name in EXPOSURES},
                **{name: "mean" for name in ENDPOINTS},
            }
        ),
    }

    rows: list[dict[str, object]] = []
    seed = 205335
    for cohort_name, cohort in cohorts.items():
        for exposure, exposure_label in EXPOSURES.items():
            for endpoint, endpoint_label in ENDPOINTS.items():
                pair = cohort[[exposure, endpoint]].dropna()
                if len(pair) < 5:
                    continue
                result = stats.spearmanr(pair[exposure], pair[endpoint])
                low, high = bootstrap_spearman(
                    pair[exposure].to_numpy(),
                    pair[endpoint].to_numpy(),
                    seed,
                )
                seed += 1
                rows.append(
                    {
                        "cohort": cohort_name,
                        "exposure": exposure,
                        "exposure_label": exposure_label,
                        "endpoint": endpoint,
                        "endpoint_label": endpoint_label,
                        "n": len(pair),
                        "spearman_rho": result.statistic,
                        "ci95_low": low,
                        "ci95_high": high,
                        "p_value": result.pvalue,
                    }
                )
    table = pd.DataFrame(rows)
    table["fdr_bh"] = table.groupby(["cohort", "exposure"])[
        "p_value"
    ].transform(bh_adjust)
    return table


def make_plot(metrics: pd.DataFrame, correlations: pd.DataFrame, out: Path) -> None:
    eligible = metrics[metrics["analysis_eligible"]].copy()
    eligible["log1p_tacstd2_cpm"] = np.log1p(
        eligible["tacstd2_pseudobulk_cpm"]
    )
    palette = dict(
        zip(
            sorted(eligible["cancer_subtype"].dropna().unique()),
            sns.color_palette("colorblind"),
        )
    )
    sns.set_theme(style="ticks", context="notebook")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    plot_endpoints = [
        "tnk_fraction_nonmalignant",
        "tnk_fraction_all",
        "cd8_fraction_tnk",
        "nk_fraction_tnk",
    ]
    for ax, endpoint in zip(axes.flat, plot_endpoints, strict=True):
        sns.scatterplot(
            data=eligible,
            x="log1p_tacstd2_cpm",
            y=endpoint,
            hue="cancer_subtype",
            palette=palette,
            s=65,
            ax=ax,
            legend=endpoint == plot_endpoints[0],
        )
        sns.regplot(
            data=eligible,
            x="log1p_tacstd2_cpm",
            y=endpoint,
            scatter=False,
            color="0.35",
            ci=None,
            ax=ax,
        )
        row = correlations[
            (correlations["cohort"] == "all_eligible_samples")
            & (correlations["exposure"] == "tacstd2_pseudobulk_cpm")
            & (correlations["endpoint"] == endpoint)
        ].iloc[0]
        ax.set(
            xlabel="log1p(malignant TACSTD2 pseudobulk CPM)",
            ylabel=ENDPOINTS[endpoint],
            title=(
                f"Spearman ρ={row.spearman_rho:.2f}, "
                f"p={row.p_value:.3g}, n={int(row.n)}"
            ),
        )
    sns.despine(fig=fig)
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--identities", required=True, type=Path)
    parser.add_argument("--soft", required=True, type=Path)
    parser.add_argument("--outdir", type=Path, default=Path(__file__).parent)
    parser.add_argument("--min-cells", type=int, default=50)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    identities = pd.read_csv(args.identities, sep="\t")
    if identities["barcode"].duplicated().any():
        raise ValueError("Cell-identity barcodes are not unique")
    matrix, genes, barcodes = load_matrix(args.matrix)
    if list(genes).count("TACSTD2") != 1:
        raise ValueError("Expected exactly one TACSTD2 feature")

    indexed = identities.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise ValueError(
            f"Matrix/identity barcode mismatch: {len(missing)} missing, "
            f"{len(extra)} extra"
        )
    cells = indexed.loc[barcodes].reset_index()

    gene_index = int(np.flatnonzero(genes == "TACSTD2")[0])
    tacstd2 = np.asarray(matrix.getrow(gene_index).toarray()).ravel()
    library_size = np.asarray(matrix.sum(axis=0)).ravel()
    if np.any(library_size <= 0):
        raise ValueError("Found cells with zero total UMI count")
    cells["tacstd2_umi"] = tacstd2
    cells["library_umi"] = library_size

    rows: list[dict[str, object]] = []
    for sample, group in cells.groupby("orig.ident", sort=True):
        malignant = group["lineage.sub"].eq("Malignant cells")
        tnk = group["lineage.total"].eq("T/NK cells")
        cd4 = group["lineage.sub"].eq("CD4+ T cells")
        cd8 = group["lineage.sub"].eq("CD8+ T cells")
        nk = group["lineage.sub"].eq("NK cells")
        malignant_library = group.loc[malignant, "library_umi"].sum()
        malignant_gene = group.loc[malignant, "tacstd2_umi"]
        n_nonmalignant = len(group) - int(malignant.sum())
        rows.append(
            {
                "orig.ident": sample,
                "n_cells": len(group),
                "n_malignant": int(malignant.sum()),
                "n_tnk": int(tnk.sum()),
                "n_cd4": int(cd4.sum()),
                "n_cd8": int(cd8.sum()),
                "n_nk": int(nk.sum()),
                "malignant_fraction_all": malignant.mean(),
                "tnk_fraction_all": tnk.mean(),
                "tnk_fraction_nonmalignant": (
                    tnk.sum() / n_nonmalignant if n_nonmalignant else np.nan
                ),
                "cd4_fraction_tnk": cd4.sum() / tnk.sum() if tnk.any() else np.nan,
                "cd8_fraction_tnk": cd8.sum() / tnk.sum() if tnk.any() else np.nan,
                "nk_fraction_tnk": nk.sum() / tnk.sum() if tnk.any() else np.nan,
                "tacstd2_malignant_umi": malignant_gene.sum(),
                "malignant_library_umi": malignant_library,
                "tacstd2_pseudobulk_cpm": (
                    1e6 * malignant_gene.sum() / malignant_library
                    if malignant_library
                    else np.nan
                ),
                "tacstd2_detection_fraction": (
                    malignant_gene.gt(0).mean() if malignant.any() else np.nan
                ),
                "tacstd2_mean_log1p_cp10k": (
                    np.log1p(
                        1e4
                        * malignant_gene
                        / group.loc[malignant, "library_umi"]
                    ).mean()
                    if malignant.any()
                    else np.nan
                ),
            }
        )

    metrics = pd.DataFrame(rows)
    geo = parse_geo_soft(args.soft)
    metrics = metrics.merge(
        geo[
            [
                "orig.ident",
                "gsm",
                "patient",
                "tissue",
                "tumor_stage",
                "cancer_subtype",
                "recist",
                "platform",
            ]
        ],
        on="orig.ident",
        how="left",
        validate="one_to_one",
    )
    if metrics["gsm"].isna().any():
        raise ValueError("Some identity-table samples did not match GEO metadata")
    metrics["analysis_eligible"] = (
        (metrics["n_malignant"] >= args.min_cells)
        & (metrics["n_tnk"] >= args.min_cells)
    )
    metrics.to_csv(args.outdir / "sample_metrics.tsv", sep="\t", index=False)

    correlations = correlation_table(metrics, args.min_cells)
    correlations.to_csv(args.outdir / "correlations.tsv", sep="\t", index=False)
    make_plot(
        metrics,
        correlations,
        args.outdir / "A3_GSE205335_TACSTD2_vs_TNK.png",
    )

    provenance = pd.DataFrame(
        [
            {
                "file": path.name,
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in [args.identities, args.matrix, args.soft]
        ]
    )
    provenance.to_csv(args.outdir / "provenance.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
