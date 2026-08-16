#!/usr/bin/env python3
"""Sample-level co-occurrence analysis for public NSCLC scRNA-seq atlases.

This script intentionally does not call the derived quantities spatial
"neighborhoods": the source atlases are dissociated.  It computes a
sample-level co-occurrence proxy from lineage-restricted expression programs.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import math
import re
import tarfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


GENE_SETS = {
    "epi_class": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1", "TACSTD2", "CLDN4"],
    "t_class": ["CD3D", "CD3E", "TRBC1", "TRBC2", "LCK", "PTPRC"],
    "b_class": ["CD79A", "CD79B", "MS4A1", "CD37", "CD74", "CD22", "CD19"],
    "myeloid_class": ["LST1", "TYROBP", "FCER1G", "CTSS", "AIF1", "LYZ"],
    "endothelial_class": ["VWF", "PECAM1", "EMCN", "KDR"],
    "fibroblast_class": ["COL1A1", "COL1A2", "DCN", "COL3A1", "COL6A1"],
    "epithelial_module": ["TACSTD2", "CLDN4"],
    "tacstd2": ["TACSTD2"],
    "cldn4": ["CLDN4"],
    "tj_barrier": ["CLDN3", "CLDN7", "OCLN", "TJP1", "F11R"],
    "cd8": ["CD8A", "CD8B", "CCL5", "NKG7", "GZMK"],
    "tls_b": ["MS4A1", "CD79A", "CD74", "CD37"],
    "tls_tfh": ["CXCL13", "CXCR5", "ICOS", "PDCD1", "TOX"],
    "tls_mregdc": ["LAMP3", "CCL19", "CCL22", "CD274"],
    "lcam_t": ["PDCD1", "CXCL13", "TOX", "HAVCR2"],
    "lcam_plasma": ["IGHG1", "IGHG3", "MZB1", "JCHAIN"],
    "lcam_mac": ["SPP1", "APOC1", "APOE", "LPL"],
}
SELECTED_GENES = sorted({g for genes in GENE_SETS.values() for g in genes})
CLASS_NAMES = ["epithelial", "T", "B", "myeloid", "endothelial", "fibroblast"]
CLASS_SETS = [
    "epi_class",
    "t_class",
    "b_class",
    "myeloid_class",
    "endothelial_class",
    "fibroblast_class",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_dense_gene_rows(path: Path) -> tuple[list[str], dict[str, np.ndarray], np.ndarray, int]:
    """Read a gene x cell tab-delimited matrix and retain only selected genes."""
    selected = {}
    with open_text(path) as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        first_line = handle.readline()
        first_gene, sep, first_values = first_line.partition("\t")
        if not sep:
            raise ValueError(f"{path}: first data row is not tab-delimited")
        n_values = first_values.count("\t") + 1
        if n_values == len(header):
            # GSE148071 has no gene-column label in the header.
            cell_ids = header
        elif n_values == len(header) - 1:
            # GSE131907 labels the leading gene column.
            cell_ids = header[1:]
        else:
            raise ValueError(
                f"{path}: first row has {n_values} values for a {len(header)}-field header"
            )
        totals = np.zeros(len(cell_ids), dtype=np.float64)
        n_genes = 0
        for line in itertools.chain([first_line], handle):
            gene, sep, values = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(values, sep="\t", dtype=np.float64)
            if arr.size != totals.size:
                raise ValueError(f"{path}: {gene} has {arr.size} values; expected {totals.size}")
            totals += arr
            n_genes += 1
            gene = gene.split(".")[0].upper()
            if gene in SELECTED_GENES:
                if gene in selected:
                    selected[gene] += arr
                else:
                    selected[gene] = arr
    return cell_ids, selected, totals, n_genes


def normalized_selected(
    selected: dict[str, np.ndarray], totals: np.ndarray
) -> dict[str, np.ndarray]:
    scale = np.divide(10_000.0, totals, out=np.zeros_like(totals), where=totals > 0)
    return {
        gene: np.log1p(selected.get(gene, np.zeros_like(totals)) * scale)
        for gene in SELECTED_GENES
    }


def module(expr: dict[str, np.ndarray], name: str) -> np.ndarray:
    present = [expr[g] for g in GENE_SETS[name] if g in expr]
    if not present:
        return np.zeros(len(next(iter(expr.values()))), dtype=float)
    return np.mean(present, axis=0)


def classify(expr: dict[str, np.ndarray]) -> np.ndarray:
    scores = np.vstack([module(expr, name) for name in CLASS_SETS])
    labels = np.array(CLASS_NAMES, dtype=object)[np.argmax(scores, axis=0)]
    labels[np.max(scores, axis=0) == 0] = "unclassified"
    return labels


def robust_z(values: pd.Series) -> pd.Series:
    median = values.median()
    mad = np.median(np.abs(values - median))
    if not np.isfinite(mad) or mad == 0:
        sd = values.std(ddof=1)
        return (values - values.mean()) / (sd if sd > 0 else 1.0)
    return (values - median) / (1.4826 * mad)


def finalize_scores(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    component_cols = [c for c in out if c.endswith("_component")]
    for col in component_cols:
        out[f"{col}_z"] = robust_z(out[col])
    out["cd8_score"] = out["cd8_component_z"]
    out["tls_score"] = out[
        ["tls_b_component_z", "tls_tfh_component_z", "tls_mregdc_component_z"]
    ].mean(axis=1)
    out["lcam_score"] = out[
        ["lcam_t_component_z", "lcam_plasma_component_z", "lcam_mac_component_z"]
    ].mean(axis=1)
    out["immune_score"] = out[["cd8_score", "tls_score", "lcam_score"]].mean(axis=1)
    return out


def sample_summaries(
    atlas: str,
    samples: np.ndarray,
    expr: dict[str, np.ndarray],
    totals: np.ndarray,
    labels: np.ndarray,
    include: np.ndarray,
) -> pd.DataFrame:
    rows = []
    epi_score = module(expr, "epithelial_module")
    trop2 = module(expr, "tacstd2")
    cldn4 = module(expr, "cldn4")
    tj_score = module(expr, "tj_barrier")
    components = {
        "cd8_component": module(expr, "cd8"),
        "tls_b_component": module(expr, "tls_b"),
        "tls_tfh_component": module(expr, "tls_tfh"),
        "tls_mregdc_component": module(expr, "tls_mregdc"),
        "lcam_t_component": module(expr, "lcam_t"),
        "lcam_plasma_component": module(expr, "lcam_plasma"),
        "lcam_mac_component": module(expr, "lcam_mac"),
    }
    for sample in sorted(np.unique(samples[include])):
        sm = include & (samples == sample)
        masks = {label: sm & (labels == label) for label in CLASS_NAMES}
        row = {
            "atlas": atlas,
            "sample": sample,
            "n_cells": int(sm.sum()),
            "n_epithelial": int(masks["epithelial"].sum()),
            "n_T": int(masks["T"].sum()),
            "n_B": int(masks["B"].sum()),
            "n_myeloid": int(masks["myeloid"].sum()),
            "epithelial_fraction": float(masks["epithelial"].sum() / max(sm.sum(), 1)),
            "immune_fraction": float(
                (masks["T"].sum() + masks["B"].sum() + masks["myeloid"].sum())
                / max(sm.sum(), 1)
            ),
            "median_umi": float(np.median(totals[sm])),
            "epithelial_module": float(np.median(epi_score[masks["epithelial"]]))
            if masks["epithelial"].any()
            else np.nan,
            "tacstd2": float(np.median(trop2[masks["epithelial"]]))
            if masks["epithelial"].any()
            else np.nan,
            "cldn4": float(np.median(cldn4[masks["epithelial"]]))
            if masks["epithelial"].any()
            else np.nan,
            "tj_module": float(np.median(tj_score[masks["epithelial"]]))
            if masks["epithelial"].any()
            else np.nan,
        }
        lineage = {
            "cd8_component": "T",
            "tls_b_component": "B",
            "tls_tfh_component": "T",
            "tls_mregdc_component": "myeloid",
            "lcam_t_component": "T",
            "lcam_plasma_component": "B",
            "lcam_mac_component": "myeloid",
        }
        for name, values in components.items():
            mask = masks[lineage[name]]
            row[name] = float(np.mean(values[mask])) if mask.any() else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    return finalize_scores(out)


def process_gse131907(source: Path) -> tuple[pd.DataFrame, dict]:
    annotation_path = source / "GSE131907_cell_annotation.txt"
    matrix_path = source / "GSE131907_raw_UMI_matrix.txt.gz"
    annotation = pd.read_csv(annotation_path, sep="\t", dtype=str).set_index("Index")
    cell_ids, selected, totals, n_genes = read_dense_gene_rows(matrix_path)
    metadata = annotation.reindex(cell_ids)
    if metadata["Sample"].isna().any():
        raise ValueError("GSE131907 matrix/annotation cell IDs do not align")
    expr = normalized_selected(selected, totals)
    labels = np.full(len(cell_ids), "unclassified", dtype=object)
    major = metadata["Cell_type"].to_numpy()
    # Only the author's malignant-cell subtype is treated as tumor epithelium.
    labels[metadata["Cell_subtype"].to_numpy() == "Malignant cells"] = "epithelial"
    labels[major == "T lymphocytes"] = "T"
    labels[major == "B lymphocytes"] = "B"
    labels[major == "Myeloid cells"] = "myeloid"
    labels[major == "Endothelial cells"] = "endothelial"
    labels[major == "Fibroblasts"] = "fibroblast"
    malignant_samples = metadata.loc[
        metadata["Cell_subtype"] == "Malignant cells", "Sample"
    ].unique()
    include = metadata["Sample"].isin(malignant_samples).to_numpy()
    result = sample_summaries(
        "GSE131907", metadata["Sample"].to_numpy(), expr, totals, labels, include
    )
    site_by_sample = metadata.groupby("Sample")["Sample_Origin"].first()
    result["site"] = result["sample"].map(site_by_sample)
    audit = {
        "matrix_cells": len(cell_ids),
        "matrix_genes": n_genes,
        "tumor_bearing_specimen_cells": int(include.sum()),
        "tumor_bearing_specimens": int(metadata.loc[include, "Sample"].nunique()),
        "classification": (
            "published major immune cell types; published Malignant cells subtype "
            "for tumor epithelium"
        ),
        "genes_found": sorted(selected),
        "genes_missing": sorted(set(SELECTED_GENES) - set(selected)),
    }
    return result, audit


def extract_gse148071(source: Path) -> list[Path]:
    paths = sorted(source.glob("GSM*_exp.txt.gz"))
    if paths:
        return paths
    archive = source / "GSE148071_RAW.tar"
    with tarfile.open(archive) as tar:
        members = [m for m in tar.getmembers() if re.match(r"GSM.*_exp\.txt\.gz$", m.name)]
        tar.extractall(source, members=members, filter="data")
    return sorted(source.glob("GSM*_exp.txt.gz"))


def process_gse148071(source: Path) -> tuple[pd.DataFrame, dict]:
    frames = []
    audits = []
    for path in extract_gse148071(source):
        cell_ids, selected, totals, n_genes = read_dense_gene_rows(path)
        expr = normalized_selected(selected, totals)
        labels = classify(expr)
        sample_match = re.search(r"_(P\d+)_", path.name)
        sample = sample_match.group(1) if sample_match else path.stem
        summary = sample_summaries(
                "GSE148071",
                np.repeat(sample, len(cell_ids)),
                expr,
                totals,
                labels,
                np.ones(len(cell_ids), dtype=bool),
            )
        summary["site"] = "advanced_biopsy"
        frames.append(summary)
        audits.append(
            {
                "file": path.name,
                "sample": sample,
                "cells": len(cell_ids),
                "genes": n_genes,
                "class_counts": {
                    key: int(np.sum(labels == key)) for key in CLASS_NAMES + ["unclassified"]
                },
                "genes_missing": sorted(set(SELECTED_GENES) - set(selected)),
            }
        )
    result = finalize_scores(pd.concat(frames, ignore_index=True))
    audit = {
        "matrix_cells": sum(x["cells"] for x in audits),
        "tumor_samples": len(audits),
        "classification": "predeclared marker argmax (no published cell labels in GEO archive)",
        "sample_audit": audits,
    }
    return result, audit


def _read_10x_selected(directory: Path, keep_barcodes: dict[str, str]) -> tuple[list[str], dict[str, np.ndarray], np.ndarray]:
    feature_path = next(directory.glob("*features.tsv"))
    barcode_path = next(directory.glob("*barcodes.tsv"))
    mtx_path = next(directory.glob("*matrix.mtx"))
    symbols = []
    with feature_path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            symbols.append(fields[1].split(".")[0].upper() if len(fields) > 1 else fields[0].upper())
    gene_rows = {i + 1: gene for i, gene in enumerate(symbols) if gene in SELECTED_GENES}
    col_to_cell = {}
    with barcode_path.open() as handle:
        for index, line in enumerate(handle, start=1):
            barcode = line.strip()
            if barcode in keep_barcodes:
                col_to_cell[index] = keep_barcodes[barcode]
    if not col_to_cell:
        return [], {}, np.array([])
    totals = {col: 0.0 for col in col_to_cell}
    selected = {gene: {col: 0.0 for col in col_to_cell} for gene in gene_rows.values()}
    with mtx_path.open() as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            break
        for line in handle:
            gene_i, cell_j, value = line.split()
            cell_j = int(cell_j)
            if cell_j not in col_to_cell:
                continue
            count = float(value)
            totals[cell_j] += count
            gene_i = int(gene_i)
            if gene_i in gene_rows:
                selected[gene_rows[gene_i]][cell_j] += count
    columns = sorted(col_to_cell)
    cell_ids = [col_to_cell[col] for col in columns]
    total_arr = np.array([totals[col] for col in columns], dtype=float)
    selected_arr = {
        gene: np.array([values[col] for col in columns], dtype=float)
        for gene, values in selected.items()
    }
    return cell_ids, selected_arr, total_arr


def process_gse154826(source: Path) -> tuple[pd.DataFrame, dict]:
    cells = pd.read_csv(source / "GSE154826_cell_metadata.csv")
    annots = pd.read_csv(source / "GSE154826_annots_list.csv")
    meta = pd.read_csv(source / "GSE154826_sample_metadata.csv")
    lineage_map = dict(zip(annots["cluster"].astype(int), annots["lineage"].astype(str)))
    tumor = meta[(meta["tissue"] == "Tumor") & (meta["Use.in.Clustering.Model."] == "Yes")]
    tumor_samples = set(tumor["sample_ID"].astype(str))
    batch_dir = source / "gse154826_batches"
    frames = []
    audits = []
    for batch in sorted({int(x) for x in tumor["amp_batch_ID"]}):
        archive = batch_dir / f"GSE154826_amp_batch_ID_{batch}.tar.gz"
        if not archive.exists():
            raise FileNotFoundError(archive)
        extract_dir = batch_dir / f"extracted_{batch}"
        extract_dir.mkdir(exist_ok=True)
        if not any(extract_dir.glob("*matrix.mtx")):
            with tarfile.open(archive) as tar:
                tar.extractall(extract_dir, filter="data")
        batch_samples = meta.loc[meta["amp_batch_ID"].astype(int) == batch, "sample_ID"].astype(str)
        keep = {}
        for sample in batch_samples:
            sample_cells = cells.loc[cells["sample_ID"].astype(str) == sample, "cell_ID"]
            prefix = f"{sample}_"
            for cell_id in sample_cells:
                barcode = cell_id[len(prefix) :] if cell_id.startswith(prefix) else cell_id.split("_", 1)[-1]
                keep[barcode] = cell_id
        cell_ids, selected, totals = _read_10x_selected(extract_dir, keep)
        if not cell_ids:
            audits.append({"batch": batch, "matched_cells": 0})
            continue
        expr = normalized_selected(selected, totals)
        marker_labels = classify(expr)
        cell_info = cells.set_index("cell_ID").reindex(cell_ids)
        labels = np.array(["unclassified"] * len(cell_ids), dtype=object)
        for i, cluster in enumerate(cell_info["cluster_ID"].to_numpy()):
            lineage = lineage_map.get(int(cluster), "")
            if lineage == "T":
                labels[i] = "T"
            elif lineage == "B&plasma":
                labels[i] = "B"
            elif lineage == "MNP":
                labels[i] = "myeloid"
            elif lineage == "epi_endo_fibro_doublet":
                labels[i] = marker_labels[i]
        samples = cell_info["sample_ID"].astype(str).to_numpy()
        include = np.array([sample in tumor_samples for sample in samples])
        if not include.any():
            continue
        summary = sample_summaries("GSE154826", samples, expr, totals, labels, include)
        site = tumor.set_index(tumor["sample_ID"].astype(str))["disease"]
        summary["site"] = summary["sample"].map(site).fillna("Tumor")
        frames.append(summary)
        audits.append(
            {
                "batch": batch,
                "matched_cells": len(cell_ids),
                "tumor_cells": int(include.sum()),
                "genes_missing": sorted(set(SELECTED_GENES) - set(selected)),
            }
        )
    if not frames:
        raise RuntimeError("GSE154826 produced no tumor sample summaries")
    result = finalize_scores(pd.concat(frames, ignore_index=True))
    result = result.groupby(["atlas", "sample"], as_index=False).first()
    result = finalize_scores(result)
    audit = {
        "status": "included_per_sample_geo_mtx",
        "matrix_cells": int(sum(x.get("matched_cells", 0) for x in audits)),
        "tumor_samples": int(result["sample"].nunique()),
        "classification": (
            "author immune lineages; marker argmax for gated epi/endo/fibro cells; "
            "clustering-used tumor samples only"
        ),
        "batch_audit": audits,
        "genes_present": bool(all("TACSTD2" not in x.get("genes_missing", []) or "CLDN4" not in x.get("genes_missing", []) for x in audits)),
    }
    return result, audit


def perm_pvalue(x: np.ndarray, y: np.ndarray, observed: float, rng, n_perm=100_000) -> float:
    xr = stats.rankdata(x)
    yr = stats.rankdata(y)
    xr = (xr - xr.mean()) / xr.std()
    yr = (yr - yr.mean()) / yr.std()
    exceed = 0
    for start in range(0, n_perm, 2_000):
        count = min(2_000, n_perm - start)
        correlations = np.empty(count)
        for i in range(count):
            correlations[i] = np.mean(xr * rng.permutation(yr))
        exceed += int(np.sum(np.abs(correlations) >= abs(observed)))
    return (exceed + 1) / (n_perm + 1)


def bootstrap_ci(x: np.ndarray, y: np.ndarray, rng, n_boot=20_000) -> tuple[float, float]:
    n = len(x)
    estimates = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        estimates[i] = stats.spearmanr(x[idx], y[idx]).statistic
    estimates = estimates[np.isfinite(estimates)]
    if estimates.size == 0:
        return np.nan, np.nan
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def residual_rank_correlation(x: np.ndarray, y: np.ndarray, covariates: np.ndarray) -> float:
    design = np.column_stack([np.ones(len(x)), *(stats.rankdata(c) for c in covariates.T)])
    xr = stats.rankdata(x)
    yr = stats.rankdata(y)
    x_resid = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    y_resid = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    return float(stats.pearsonr(x_resid, y_resid).statistic)


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    values = pvalues.to_numpy()
    order = np.argsort(values)
    adjusted = np.empty_like(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return pd.Series(np.minimum(adjusted, 1.0), index=pvalues.index)


def run_statistics(samples: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    rows = []
    min_counts = {"n_epithelial": 20, "n_T": 20, "n_B": 10, "n_myeloid": 20}
    for atlas, atlas_df in samples.groupby("atlas"):
        eligible = np.ones(len(atlas_df), dtype=bool)
        for col, threshold in min_counts.items():
            eligible &= atlas_df[col].to_numpy() >= threshold
        for score in ["cd8_score", "tls_score", "lcam_score"]:
            cols = [
                "epithelial_module",
                score,
                "epithelial_fraction",
                "n_cells",
                "median_umi",
                "site",
            ]
            frame = atlas_df.loc[eligible, cols].dropna()
            x = frame["epithelial_module"].to_numpy()
            y = frame[score].to_numpy()
            if len(frame) < 4 or np.ptp(x) == 0 or np.ptp(y) == 0:
                rows.append(
                    {
                        "atlas": atlas,
                        "outcome": score,
                        "n_samples": len(frame),
                        "spearman_rho": np.nan,
                        "bootstrap_ci_low": np.nan,
                        "bootstrap_ci_high": np.nan,
                        "permutation_p": np.nan,
                        "asymptotic_p": np.nan,
                        "partial_spearman_rho": np.nan,
                    }
                )
                continue
            rho, asymptotic_p = stats.spearmanr(x, y)
            low, high = bootstrap_ci(x, y, rng)
            site_covariates = pd.get_dummies(
                frame["site"], drop_first=True, dtype=float
            ).to_numpy()
            partial = residual_rank_correlation(
                x,
                y,
                np.column_stack(
                    [
                        frame["epithelial_fraction"].to_numpy(),
                        np.log1p(frame["n_cells"].to_numpy()),
                        np.log1p(frame["median_umi"].to_numpy()),
                        site_covariates,
                    ]
                ),
            )
            rows.append(
                {
                    "atlas": atlas,
                    "outcome": score,
                    "n_samples": len(frame),
                    "spearman_rho": rho,
                    "bootstrap_ci_low": low,
                    "bootstrap_ci_high": high,
                    "permutation_p": perm_pvalue(x, y, rho, rng),
                    "asymptotic_p": asymptotic_p,
                    "partial_spearman_rho": partial,
                }
            )
    results = pd.DataFrame(rows)
    results["bh_q_within_atlas"] = results.groupby("atlas")["permutation_p"].transform(
        bh_adjust
    )

    meta_rows = []
    for outcome, frame in results.groupby("outcome"):
        valid = frame[(frame["n_samples"] > 3) & (np.abs(frame["spearman_rho"]) < 1)]
        z = np.arctanh(valid["spearman_rho"].to_numpy())
        weights = valid["n_samples"].to_numpy() - 3
        pooled_z = np.sum(weights * z) / np.sum(weights)
        se = 1 / math.sqrt(np.sum(weights))
        q = np.sum(weights * (z - pooled_z) ** 2)
        meta_rows.append(
            {
                "outcome": outcome,
                "n_atlases": len(valid),
                "n_samples_total": int(valid["n_samples"].sum()),
                "fixed_effect_rho": math.tanh(pooled_z),
                "ci_low": math.tanh(pooled_z - 1.96 * se),
                "ci_high": math.tanh(pooled_z + 1.96 * se),
                "two_sided_p": 2 * stats.norm.sf(abs(pooled_z / se)),
                "heterogeneity_q": q,
                "heterogeneity_p": stats.chi2.sf(q, max(len(valid) - 1, 1)),
            }
        )
    return results, pd.DataFrame(meta_rows)


def _eligible(frame: pd.DataFrame) -> pd.DataFrame:
    return frame[
        (frame["n_epithelial"] >= 20)
        & (frame["n_T"] >= 20)
        & (frame["n_B"] >= 10)
        & (frame["n_myeloid"] >= 20)
    ].copy()


def _partial_rho(frame: pd.DataFrame, exposure: str, outcome: str) -> float:
    covariates = [frame["epithelial_fraction"].to_numpy()]
    if frame["site"].nunique() > 1:
        dummies = pd.get_dummies(frame["site"], drop_first=True, dtype=float)
        if not dummies.empty:
            covariates.append(dummies.to_numpy())
    return residual_rank_correlation(
        frame[exposure].to_numpy(),
        frame[outcome].to_numpy(),
        np.column_stack(covariates),
    )


def _purity_residual(values: np.ndarray, purity: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(values)), stats.rankdata(purity)])
    ranked = stats.rankdata(values)
    return ranked - design @ np.linalg.lstsq(design, ranked, rcond=None)[0]


def _match(estimate: float, pvalue: float, predicted_sign: int) -> str:
    if not np.isfinite(estimate) or not np.isfinite(pvalue) or pvalue >= 0.05:
        return "inconclusive"
    observed = 1 if estimate > 0 else -1
    return "match" if observed == predicted_sign else "mismatch"


def run_thesis_tests(samples: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for atlas, atlas_df in samples.groupby("atlas"):
        frame = _eligible(atlas_df)
        if "site" not in frame:
            frame["site"] = "tumor"
        for exposure, predicted, hypothesis in [
            ("tacstd2", -1, "TROP2-high -> immune-low after purity"),
            ("cldn4", -1, "CLDN4 TJ barrier -> immune-low after purity"),
        ]:
            needed = [exposure, "immune_score", "immune_fraction", "epithelial_fraction"]
            usable = frame.dropna(subset=needed)
            if len(usable) < 4 or np.ptp(usable[exposure]) == 0:
                rows.append(
                    {
                        "atlas": atlas,
                        "hypothesis": hypothesis,
                        "exposure": exposure,
                        "outcome": "immune_score_purity_adjusted",
                        "n": int(len(usable)),
                        "metric": "partial_spearman_rho",
                        "estimate": np.nan,
                        "p": np.nan,
                        "match": "inconclusive",
                    }
                )
                continue
            rho = _partial_rho(usable, exposure, "immune_score")
            raw_rho, raw_p = stats.spearmanr(
                usable[exposure],
                _purity_residual(
                    usable["immune_score"].to_numpy(),
                    usable["epithelial_fraction"].to_numpy(),
                ),
            )
            rows.append(
                {
                    "atlas": atlas,
                    "hypothesis": hypothesis,
                    "exposure": exposure,
                    "outcome": "immune_score_purity_adjusted",
                    "n": int(len(usable)),
                    "metric": "partial_spearman_rho",
                    "estimate": rho,
                    "p": float(raw_p),
                    "match": _match(rho, raw_p, predicted),
                }
            )
            high = usable[exposure] > usable[exposure].median()
            immune_resid = _purity_residual(
                usable["immune_score"].to_numpy(),
                usable["epithelial_fraction"].to_numpy(),
            )
            immune_low = immune_resid < np.median(immune_resid)
            table = np.array(
                [
                    [int((high & immune_low).sum()), int((high & ~immune_low).sum())],
                    [int((~high & immune_low).sum()), int((~high & ~immune_low).sum())],
                ]
            )
            odds, fisher_p = stats.fisher_exact(table, alternative="two-sided")
            rows.append(
                {
                    "atlas": atlas,
                    "hypothesis": hypothesis,
                    "exposure": exposure,
                    "outcome": "immune_low_given_high_exposure",
                    "n": int(len(usable)),
                    "metric": "OR",
                    "estimate": float(odds),
                    "p": float(fisher_p),
                    "match": _match(np.log(odds) if odds not in (0, np.inf) else np.nan, fisher_p, 1),
                }
            )
            infil_high = usable.loc[high, "immune_fraction"].to_numpy()
            infil_low = usable.loc[~high, "immune_fraction"].to_numpy()
            if infil_high.mean() > 0 and infil_low.mean() > 0:
                log2fc = float(np.log2(infil_high.mean() / infil_low.mean()))
                mw_p = float(stats.mannwhitneyu(infil_high, infil_low, alternative="two-sided").pvalue)
                rows.append(
                    {
                        "atlas": atlas,
                        "hypothesis": hypothesis,
                        "exposure": exposure,
                        "outcome": "immune_fraction_high_vs_low_exposure",
                        "n": int(len(usable)),
                        "metric": "log2FC",
                        "estimate": log2fc,
                        "p": mw_p,
                        "match": _match(log2fc, mw_p, -1),
                    }
                )
        barrier = frame.dropna(subset=["cldn4", "tj_module", "tacstd2"])
        if len(barrier) >= 4 and np.ptp(barrier["cldn4"]) > 0:
            for outcome, predicted, label in [
                ("tj_module", 1, "CLDN4 tracks tight-junction barrier"),
                ("tacstd2", 1, "CLDN4 co-occurs with TROP2 epithelial state"),
            ]:
                if np.ptp(barrier[outcome]) == 0:
                    rho, pvalue = np.nan, np.nan
                else:
                    rho, pvalue = stats.spearmanr(barrier["cldn4"], barrier[outcome])
                rows.append(
                    {
                        "atlas": atlas,
                        "hypothesis": label,
                        "exposure": "cldn4",
                        "outcome": outcome,
                        "n": int(len(barrier)),
                        "metric": "spearman_rho",
                        "estimate": float(rho) if np.isfinite(rho) else np.nan,
                        "p": float(pvalue) if np.isfinite(pvalue) else np.nan,
                        "match": _match(float(rho) if np.isfinite(rho) else np.nan, float(pvalue) if np.isfinite(pvalue) else np.nan, predicted),
                    }
                )
    return pd.DataFrame(rows)


def plot_results(samples: pd.DataFrame, statistics: pd.DataFrame, output: Path) -> None:
    atlases = [atlas for atlas in ["GSE131907", "GSE148071", "GSE154826"] if atlas in set(samples["atlas"])]
    fig, axes = plt.subplots(len(atlases), 3, figsize=(12, 3.6 * len(atlases)), constrained_layout=True)
    if len(atlases) == 1:
        axes = np.array([axes])
    colors = {"GSE131907": "#2864a6", "GSE148071": "#d17c00", "GSE154826": "#2a7f62"}
    for row, atlas in enumerate(atlases):
        atlas_df = samples[samples["atlas"] == atlas]
        eligible = (
            (atlas_df["n_epithelial"] >= 20)
            & (atlas_df["n_T"] >= 20)
            & (atlas_df["n_B"] >= 10)
            & (atlas_df["n_myeloid"] >= 20)
        )
        frame = atlas_df[eligible]
        for col, outcome in enumerate(["cd8_score", "tls_score", "lcam_score"]):
            ax = axes[row, col]
            ax.scatter(
                frame["epithelial_module"],
                frame[outcome],
                color=colors[atlas],
                alpha=0.8,
                edgecolor="white",
                linewidth=0.5,
            )
            if len(frame) >= 2:
                coef = np.polyfit(frame["epithelial_module"], frame[outcome], 1)
                grid = np.linspace(frame["epithelial_module"].min(), frame["epithelial_module"].max())
                ax.plot(grid, np.polyval(coef, grid), color="#333333", linewidth=1)
            hit = statistics[
                (statistics["atlas"] == atlas) & (statistics["outcome"] == outcome)
            ]
            if not hit.empty:
                stat = hit.iloc[0]
                ax.text(
                    0.04,
                    0.96,
                    f"ρ={stat.spearman_rho:.2f}; pperm={stat.permutation_p:.3g}; n={stat.n_samples}",
                    transform=ax.transAxes,
                    va="top",
                    fontsize=8,
                )
            ax.set_title(f"{atlas}: {outcome.replace('_score', '').upper()}")
            ax.set_xlabel("Tumor epithelial TACSTD2/CLDN4")
            ax.set_ylabel("Sample co-occurrence score (z)")
    fig.savefig(output / "module_vs_immune_scores.png", dpi=180)
    fig.savefig(output / "module_vs_immune_scores.pdf")
    plt.close(fig)


def source_manifest(source: Path) -> pd.DataFrame:
    records = []
    for path in [
        source / "GSE131907_cell_annotation.txt.gz",
        source / "GSE131907_raw_UMI_matrix.txt.gz",
        source / "GSE148071_RAW.tar",
        source / "GSE154826_cell_metadata.csv",
        source / "GSE154826_annots_list.csv",
        source / "GSE154826_sample_metadata.csv",
    ]:
        if path.exists():
            records.append(
                {
                    "file": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "used_for_statistics": path.name.startswith(("GSE131907", "GSE148071", "GSE154826")),
                }
            )
    batch_dir = source / "gse154826_batches"
    if batch_dir.exists():
        for path in sorted(batch_dir.glob("GSE154826_amp_batch_ID_*.tar.gz")):
            records.append(
                {
                    "file": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "used_for_statistics": True,
                }
            )
    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260816)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    def load_or_build(cache_csv: Path, cache_json: Path, builder):
        if cache_csv.exists() and cache_json.exists():
            frame = pd.read_csv(cache_csv)
            if {"tacstd2", "cldn4", "tj_module"}.issubset(frame.columns):
                with cache_json.open() as handle:
                    return frame, json.load(handle)
        frame, audit = builder()
        frame.to_csv(cache_csv, index=False)
        with cache_json.open("w") as handle:
            json.dump(audit, handle, indent=2)
        return frame, audit

    g131, audit131 = load_or_build(
        args.output / "_cache_v4_gse131907_scores.csv",
        args.output / "_cache_v4_gse131907_audit.json",
        lambda: process_gse131907(args.source),
    )
    g148, audit148 = load_or_build(
        args.output / "_cache_v4_gse148071_scores.csv",
        args.output / "_cache_v4_gse148071_audit.json",
        lambda: process_gse148071(args.source),
    )
    if "site" not in g148:
        g148["site"] = "advanced_biopsy"
    g154, audit154 = load_or_build(
        args.output / "_cache_v4_gse154826_scores.csv",
        args.output / "_cache_v4_gse154826_audit.json",
        lambda: process_gse154826(args.source),
    )
    g131 = finalize_scores(g131)
    g148 = finalize_scores(g148)
    g154 = finalize_scores(g154)
    samples = pd.concat([g131, g148, g154], ignore_index=True)
    statistics, meta = run_statistics(samples, args.seed)
    thesis = run_thesis_tests(samples)
    samples.to_csv(args.output / "sample_scores.csv", index=False)
    statistics.to_csv(args.output / "association_statistics.csv", index=False)
    meta.to_csv(args.output / "meta_analysis.csv", index=False)
    thesis.to_csv(args.output / "thesis_results.csv", index=False)
    source_manifest(args.source).to_csv(args.output / "source_manifest.csv", index=False)
    with (args.output / "audit.json").open("w") as handle:
        json.dump(
            {
                "seed": args.seed,
                "gene_sets": GENE_SETS,
                "minimum_lineage_cells": {
                    "epithelial": 20,
                    "T": 20,
                    "B": 10,
                    "myeloid": 20,
                },
                "GSE131907": audit131,
                "GSE148071": audit148,
                "GSE154826": audit154,
            },
            handle,
            indent=2,
        )
    plot_results(samples, statistics, args.output)


if __name__ == "__main__":
    main()
