#!/usr/bin/env python3
"""Milo-style kNN neighbourhood DA vs malignant CLDN4 on GSE205335.

Public processed GEO only (UMI dgCMatrix + author cell identity + SOFT
characteristics). miloR / edgeR are not used. Sample/patient is the
independent unit. Do not cite cell count as n.

Primary contrast: neighbourhood abundance vs patient-level mean
log1p(CP10k) CLDN4 in author-labelled malignant cells (Spearman +
SpatialFDR). Composition: are CLDN4-high neighbourhoods T/NK-poor?

MPR/NMPR is not labelled. RECIST R vs NR is a secondary Welch test
(n=6 vs 10) and is not substituted for MPR.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
import tempfile
import time
import urllib.request
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from knn_nhood import (  # noqa: E402
    EPS,
    attach_fdr,
    count_matrix,
    greedy_independent,
    make_nhoods,
    nhood_composition,
    nhood_membership,
    pca_knn,
    refine_indices,
    sample_paired_tnk_by_tacstd2,
    sample_sizes,
    select_hvg,
    self_test,
    spearman_da,
    spearman_safe,
    welch_da,
    wilcoxon_paired,
    write_json,
)

GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl"
GEO_SOFT = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz"
MATRIX_NAME = "GSE205335_Lung_IO_UMI_matrix.rds.gz"
IDENT_NAME = "GSE205335_Lung_IO_CellIdentity.txt.gz"
RESPONSE_MAP = {"PR": "R", "CR": "R", "SD": "NR", "PD": "NR", "NE": "NE"}
NSCLC = {"ADC", "SQ"}
NORMAL_TISSUE_PREFIX = "Normal"
ALWAYS = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD3D", "NKG7"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET {url}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)
    return dest


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
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def gunzip_until_rds(src: Path, dest: Path) -> Path:
    """GEO ships a double-gzipped RDS. Peel gzip until the R XDR magic."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        with dest.open("rb") as handle:
            magic = handle.read(2)
        if magic == b"X\n":
            print(f"RDS ready {dest} ({dest.stat().st_size} bytes)", flush=True)
            return dest
    current = src
    tmp_dir = dest.parent
    for i in range(4):
        with current.open("rb") as handle:
            magic = handle.read(2)
        if magic == b"\x1f\x8b":
            nxt = tmp_dir / f"{dest.name}.peel{i}"
            print(f"gunzip peel {i}: {current}", flush=True)
            with gzip.open(current, "rb") as source, nxt.open("wb") as out:
                shutil.copyfileobj(source, out, 16 * 1024 * 1024)
            if current != src and current.exists():
                current.unlink()
            current = nxt
            continue
        break
    if current != dest:
        current.replace(dest)
    print(f"RDS ready {dest} ({dest.stat().st_size} bytes)", flush=True)
    return dest


def load_dgcmatrix(rds_path: Path):
    import rdata

    t0 = time.time()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
        obj = rdata.read_rds(rds_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    print(
        f"matrix {matrix.shape} nnz={matrix.nnz} load={time.time() - t0:.1f}s",
        flush=True,
    )
    return matrix, genes, barcodes


def gene_mean_var(matrix: sparse.spmatrix) -> tuple[np.ndarray, np.ndarray]:
    n = matrix.shape[1]
    mean = np.asarray(matrix.mean(axis=1)).ravel()
    mean_sq = np.asarray(matrix.multiply(matrix).mean(axis=1)).ravel()
    var = np.maximum(mean_sq - mean * mean, 0.0) * (n / max(n - 1, 1))
    return mean, var


def plot_volcano(df: pd.DataFrame, x: str, title: str, path: Path, xlabel: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    testable = df["testable"].to_numpy()
    ax.scatter(
        df.loc[~testable, x],
        -np.log10(np.clip(df.loc[~testable, "p"], EPS, 1)),
        s=6,
        c="#bbbbbb",
        alpha=0.4,
        label="not testable",
    )
    t = df[testable]
    sig = t["SpatialFDR"] < 0.1
    ax.scatter(
        t.loc[~sig, x],
        -np.log10(np.clip(t.loc[~sig, "p"], EPS, 1)),
        s=8,
        c="#4c72b0",
        alpha=0.55,
        label="testable SpatialFDR≥0.1",
    )
    ax.scatter(
        t.loc[sig, x],
        -np.log10(np.clip(t.loc[sig, "p"], EPS, 1)),
        s=14,
        c="#c44e52",
        alpha=0.85,
        label="SpatialFDR<0.1",
    )
    ax.axhline(-np.log10(0.05), ls="--", c="0.5", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"−log10 p (patient-level test)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_scatter(x, y, path: Path, xlabel: str, ylabel: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.scatter(x, y, s=8, c="#4c72b0", alpha=0.45)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fdr_counts(df: pd.DataFrame) -> dict:
    t = df[df["testable"]]
    p = t["p"] if len(t) else pd.Series(dtype=float)
    return {
        "n_nhoods_total": int(len(df)),
        "n_testable": int(len(t)),
        "n_p_lt_0.05": int((p < 0.05).sum()) if len(t) else 0,
        "min_p": float(p.min()) if len(t) else None,
        "min_SpatialFDR": float(t["SpatialFDR"].min()) if len(t) else None,
        "min_BH_FDR": float(t["BH_FDR"].min()) if len(t) else None,
        "n_BH_FDR_lt_0.1": int((t["BH_FDR"] < 0.1).sum()) if len(t) else 0,
        "n_BH_FDR_lt_0.05": int((t["BH_FDR"] < 0.05).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.2": int((t["SpatialFDR"] < 0.2).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.1": int((t["SpatialFDR"] < 0.1).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.05": int((t["SpatialFDR"] < 0.05).sum()) if len(t) else 0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/GSE205335"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/gse205335_milo_cldn4"))
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--d", type=int, default=30)
    ap.add_argument("--prop", type=float, default=0.1)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-interface", type=int, default=3)
    ap.add_argument("--min-malignant", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    outdir = args.outdir
    tables = outdir / "tables"
    figures = outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    self_test()

    ident_path = fetch(f"{GEO_SUPPL}/{IDENT_NAME}", args.datadir / IDENT_NAME)
    soft_path = fetch(GEO_SOFT, args.datadir / "GSE205335_family.soft.gz")
    matrix_gz = fetch(f"{GEO_SUPPL}/{MATRIX_NAME}", args.datadir / MATRIX_NAME)
    rds_path = gunzip_until_rds(matrix_gz, args.datadir / "GSE205335_Lung_IO_UMI_matrix.rds")

    provenance = {
        "ident_sha256": sha256(ident_path),
        "soft_sha256": sha256(soft_path),
        "matrix_gz_sha256": sha256(matrix_gz),
        "matrix_gz_bytes": int(matrix_gz.stat().st_size),
        "rds_bytes": int(rds_path.stat().st_size),
    }

    meta = parse_geo_soft(soft_path)
    meta.to_csv(tables / "gsm_sample_metadata.csv", index=False)
    ident = pd.read_csv(ident_path, sep="\t")
    matrix, genes, barcodes = load_dgcmatrix(rds_path)

    if len(barcodes) != len(ident):
        raise ValueError(f"barcode/identity length mismatch {len(barcodes)} vs {len(ident)}")
    if not np.array_equal(barcodes, ident["barcode"].to_numpy(dtype=str)):
        # Identity table is a permutation of the matrix barcodes.
        order = pd.Index(ident["barcode"]).get_indexer(barcodes)
        if (order < 0).any():
            raise ValueError("identity barcodes do not match matrix barcodes")
        ident = ident.iloc[order].reset_index(drop=True)
        if not np.array_equal(barcodes, ident["barcode"].to_numpy(dtype=str)):
            raise ValueError("failed to align identity table to matrix barcodes")

    cells = ident.merge(
        meta[
            [
                "orig.ident",
                "gsm",
                "patient",
                "tissue",
                "recist",
                "platform",
                "cancer_subtype",
                "tumor_stage",
            ]
        ],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()
        raise ValueError(f"identity samples missing GEO metadata: {missing}")
    cells["response"] = cells["recist"].map(RESPONSE_MAP)
    cells["is_normal_tissue"] = cells["tissue"].astype(str).str.startswith(NORMAL_TISSUE_PREFIX)
    cells["is_malig"] = cells["lineage.sub"].eq("Malignant cells")
    cells["is_tnk"] = cells["lineage.total"].eq("T/NK cells")
    cells["is_epi"] = cells["lineage.total"].eq("Epithelial cells")

    gene_idx = {g: i for i, g in enumerate(genes)}
    for g in ALWAYS:
        if g not in gene_idx:
            raise ValueError(f"{g} missing from UMI matrix")

    print("gene mean/var for HVG", flush=True)
    gene_mean, gene_var = gene_mean_var(matrix)
    hvg_idx = select_hvg(gene_mean, gene_var, n_hvg=args.n_hvg)
    hvg_names = [str(genes[i]) for i in hvg_idx]
    keep_genes = sorted(set(hvg_idx.tolist()) | {gene_idx[g] for g in ALWAYS})

    lib = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    scale = np.where(lib > 0, 1e4 / lib, 0.0)
    print(f"extract {len(keep_genes)} genes (HVG+markers)", flush=True)
    sub = matrix[keep_genes, :].astype(np.float32)
    del matrix
    dense = sub.toarray()
    del sub
    log_all = np.log1p(dense * scale).astype(np.float32)
    name_of = {keep_genes[i]: i for i in range(len(keep_genes))}
    cldn4 = log_all[name_of[gene_idx["CLDN4"]]]
    tacstd2 = log_all[name_of[gene_idx["TACSTD2"]]]
    epcam = log_all[name_of[gene_idx["EPCAM"]]]
    ptprc = log_all[name_of[gene_idx["PTPRC"]]]
    hvg_rows = [name_of[i] for i in hvg_idx if i in name_of]
    X_all = log_all[hvg_rows].T.copy()
    del dense, log_all

    cells["cldn4_log1p_cp10k"] = cldn4
    cells["tacstd2_log1p_cp10k"] = tacstd2
    cells["total_umi"] = lib

    # Graph: non-normal tissues only. DA unit = patient.
    keep = ~cells["is_normal_tissue"].to_numpy()
    print(f"non-normal cells {int(keep.sum())} / {len(cells)}", flush=True)
    idx = np.flatnonzero(keep)
    X = X_all[idx]
    sample = cells.loc[keep, "patient"].to_numpy(dtype=object)
    lineage = cells.loc[keep, "lineage.total"].to_numpy(dtype=object)
    is_malig = cells.loc[keep, "is_malig"].to_numpy()
    is_tnk = cells.loc[keep, "is_tnk"].to_numpy()
    cldn4_k = cldn4[idx]
    n_k = int(keep.sum())

    print(f"HVG matrix {X.shape}; PCA+kNN k={args.k} d={args.d}", flush=True)
    pcs, knn_idx, knn_dist = pca_knn(X, sample, n_pcs=args.d, k=args.k, random_state=args.seed)
    indices = refine_indices(pcs, knn_idx, prop=args.prop, random_state=args.seed)
    nhoods = make_nhoods(knn_idx, knn_dist, indices)
    members = nhood_membership(nhoods, n_k)
    print(f"nhoods={len(members)} median_size={np.median([m.size for m in members]):.0f}", flush=True)

    sample_levels = sorted(pd.unique(sample))
    counts = count_matrix(members, sample, sample_levels)
    sizes = sample_sizes(sample, sample_levels)

    patient_meta = (
        cells.loc[keep]
        .groupby("patient", observed=True)
        .agg(
            recist=("recist", "first"),
            response=("response", "first"),
            cancer_subtype=("cancer_subtype", "first"),
            n_samples=("orig.ident", "nunique"),
            tissues=("tissue", lambda s: ",".join(sorted(set(map(str, s))))),
        )
        .reindex(sample_levels)
    )
    mal_score = []
    mal_pct = []
    tnk_frac = []
    n_mal_s = []
    n_tnk_s = []
    for s in sample_levels:
        m = sample == s
        n_mal = int((m & is_malig).sum())
        n_tnk = int((m & is_tnk).sum())
        n_mal_s.append(n_mal)
        n_tnk_s.append(n_tnk)
        mal_score.append(float(cldn4_k[m & is_malig].mean()) if n_mal >= args.min_malignant else np.nan)
        mal_pct.append(float((cldn4_k[m & is_malig] > 0).mean()) if n_mal >= args.min_malignant else np.nan)
        denom = n_mal + n_tnk
        tnk_frac.append(float(n_tnk / denom) if denom else np.nan)
    mal_score = np.asarray(mal_score, dtype=float)
    mal_pct = np.asarray(mal_pct, dtype=float)
    tnk_frac = np.asarray(tnk_frac, dtype=float)
    group = patient_meta["response"].to_numpy(dtype=object)
    subtype = patient_meta["cancer_subtype"].to_numpy(dtype=object)

    da_cldn4 = attach_fdr(spearman_da(counts, sizes, mal_score, min_samples=5), nhoods.k_distance)
    recist_mask = np.isin(group, ["R", "NR"])
    da_recist = attach_fdr(
        welch_da(counts, sizes, group, group_a="R", group_b="NR", min_per_group=2, min_samples=4),
        nhoods.k_distance,
    )
    nsclc_score = mal_score.copy()
    nsclc_score[~np.isin(subtype, list(NSCLC))] = np.nan
    da_nsclc = attach_fdr(spearman_da(counts, sizes, nsclc_score, min_samples=5), nhoods.k_distance)

    comp = nhood_composition(members, is_tnk, is_malig, cldn4_k, lineage)
    comp = comp.rename(
        columns={
            "tacstd2_all_mean": "cldn4_all_mean",
            "tacstd2_malig_mean": "cldn4_malig_mean",
        }
    )
    merged = comp.merge(da_cldn4, on="nhood", suffixes=("", "_cldn4")).merge(
        da_recist[
            ["nhood", "logFC_B_minus_A", "t", "p", "testable", "BH_FDR", "SpatialFDR", "n_samples_A", "n_samples_B"]
        ].rename(
            columns={
                "logFC_B_minus_A": "recist_logFC_NR_minus_R",
                "t": "recist_t",
                "p": "recist_p",
                "testable": "recist_testable",
                "BH_FDR": "recist_BH_FDR",
                "SpatialFDR": "recist_SpatialFDR",
                "n_samples_A": "recist_n_R",
                "n_samples_B": "recist_n_NR",
            }
        ),
        on="nhood",
    )
    merged.to_csv(tables / "nhoods.tsv", sep="\t", index=False)
    da_cldn4.to_csv(tables / "da_malignant_cldn4.tsv", sep="\t", index=False)
    da_recist.to_csv(tables / "da_recist_r_vs_nr.tsv", sep="\t", index=False)
    da_nsclc.to_csv(tables / "da_malignant_cldn4_nsclc.tsv", sep="\t", index=False)

    iface = (comp["n_malig"] >= args.min_interface) & (comp["n_tnk"] >= args.min_interface)
    all_sp = spearman_safe(comp["cldn4_malig_mean"], comp["frac_tnk"])
    iface_sp = spearman_safe(comp.loc[iface, "cldn4_malig_mean"], comp.loc[iface, "frac_tnk"])
    indep = greedy_independent(members, max_shared=0)
    indep_mask = np.zeros(len(members), dtype=bool)
    indep_mask[indep] = True
    indep_sp = spearman_safe(
        comp.loc[indep_mask & iface.to_numpy(), "cldn4_malig_mean"],
        comp.loc[indep_mask & iface.to_numpy(), "frac_tnk"],
    )
    has_mal = comp["n_malig"] >= 5
    med = float(np.nanmedian(comp.loc[has_mal, "cldn4_malig_mean"])) if has_mal.any() else np.nan
    high = has_mal.to_numpy() & (comp["cldn4_malig_mean"].to_numpy() >= med)
    low = has_mal.to_numpy() & (comp["cldn4_malig_mean"].to_numpy() < med)
    paired = sample_paired_tnk_by_tacstd2(
        members, sample, is_tnk, comp["cldn4_malig_mean"].to_numpy(), high, low
    )
    paired = paired.merge(
        pd.DataFrame(
            {
                "sample": sample_levels,
                "response": group,
                "cancer_subtype": subtype,
                "malignant_cldn4": mal_score,
                "malignant_cldn4_pctpos": mal_pct,
                "sample_frac_tnk": tnk_frac,
                "n_malignant": n_mal_s,
                "n_tnk": n_tnk_s,
                "n_cells": sizes.astype(int),
                "n_biopsies": patient_meta["n_samples"].to_numpy(),
                "tissues": patient_meta["tissues"].to_numpy(),
                "recist": patient_meta["recist"].to_numpy(),
            }
        ),
        on="sample",
    )
    paired.to_csv(tables / "patient_paired_tnk_by_nhood_cldn4.tsv", sep="\t", index=False)
    paired_test = wilcoxon_paired(paired["frac_tnk_high_nhood"], paired["frac_tnk_low_nhood"])
    patient_sp = spearman_safe(paired["malignant_cldn4"], paired["sample_frac_tnk"])
    recist_tab = paired[paired["response"].isin(["R", "NR"])]
    recist_cldn4 = {
        "n_R": int((recist_tab["response"] == "R").sum()),
        "n_NR": int((recist_tab["response"] == "NR").sum()),
        "median_R": float(recist_tab.loc[recist_tab["response"] == "R", "malignant_cldn4"].median())
        if (recist_tab["response"] == "R").any()
        else None,
        "median_NR": float(recist_tab.loc[recist_tab["response"] == "NR", "malignant_cldn4"].median())
        if (recist_tab["response"] == "NR").any()
        else None,
        "p": None,
    }
    if recist_cldn4["n_R"] >= 2 and recist_cldn4["n_NR"] >= 2:
        recist_cldn4["p"] = float(
            stats.mannwhitneyu(
                recist_tab.loc[recist_tab["response"] == "R", "malignant_cldn4"].dropna(),
                recist_tab.loc[recist_tab["response"] == "NR", "malignant_cldn4"].dropna(),
                alternative="two-sided",
            ).pvalue
        )

    da_vs_comp = spearman_safe(
        merged.loc[merged["testable"].fillna(False), "spearman_rho"],
        merged.loc[merged["testable"].fillna(False), "frac_tnk"],
    )

    sample_tab = paired[
        [
            "sample",
            "response",
            "recist",
            "cancer_subtype",
            "n_cells",
            "n_malignant",
            "n_tnk",
            "n_biopsies",
            "tissues",
            "malignant_cldn4",
            "malignant_cldn4_pctpos",
            "sample_frac_tnk",
        ]
    ].rename(columns={"sample": "patient"})
    sample_tab.to_csv(tables / "patient_scores.tsv", sep="\t", index=False)

    dropped_normal = (
        cells.loc[cells["is_normal_tissue"], ["patient", "orig.ident", "tissue", "recist", "cancer_subtype"]]
        .drop_duplicates()
        .sort_values(["patient", "orig.ident"])
    )
    dropped_normal.to_csv(tables / "excluded_normal_samples.tsv", sep="\t", index=False)

    n_with_score = int(np.isfinite(mal_score).sum())
    n_r = int((group == "R").sum())
    n_nr = int((group == "NR").sum())
    n_ne = int((group == "NE").sum())
    n_nsclc = int(np.isfinite(nsclc_score).sum())

    summary = {
        "dataset": "GSE205335",
        "citation": "Ahn / Lee et al. eLife 2024 (GEO GSE205335; processed UMI + identity)",
        "public_only": True,
        "raw_data": "not accessed; controlled EGA EGAD00001008703 skipped",
        "author_barcode_labels": True,
        "malignant_definition": "author lineage.sub == 'Malignant cells'; CNV not re-inferred",
        "tnk_definition": "author lineage.total == 'T/NK cells'",
        "predictor": "malignant CLDN4 only; TACSTD2 is companion, not a gate",
        "miloR": False,
        "da_model": "patient-level Spearman of nhood proportion vs malignant CLDN4; Welch t-test for RECIST R vs NR; not edgeR QLF",
        "spatial_fdr": "miloR graphSpatialFDR k-distance reimplementation (Dann 2022 / cydar)",
        "independent_unit": "patient (biopsies from the same patient pooled)",
        "graph": {
            "cells_all": int(len(cells)),
            "cells_non_normal": int(n_k),
            "normal_cells_excluded": int((~keep).sum()),
            "k": args.k,
            "d": min(args.d, X.shape[1], n_k - 1),
            "prop": args.prop,
            "n_hvg": len(hvg_rows),
            "n_nhoods": int(len(members)),
            "nhood_size": args.k + 1,
            "median_nhood_n_cells": float(np.median([m.size for m in members])),
            "batch": "PCA per-patient mean centering (not Harmony)",
        },
        "patients": {
            "n_geo": int(meta["patient"].nunique()),
            "n_in_graph": int(len(sample_levels)),
            "n_with_malignant_ge10": n_with_score,
            "R": n_r,
            "NR": n_nr,
            "NE": n_ne,
            "n_nsclc_with_score": n_nsclc,
            "n_normal_only_patients_excluded": int(dropped_normal["patient"].nunique()),
        },
        "lineage_counts_non_normal": {k: int(v) for k, v in pd.Series(lineage).value_counts().items()},
        "n_malignant_non_normal": int(is_malig.sum()),
        "n_tnk_non_normal": int(is_tnk.sum()),
        "marker_sanity": {
            "epcam_pos_malignant": float((epcam[idx][is_malig] > 0).mean()) if is_malig.any() else None,
            "epcam_pos_tnk": float((epcam[idx][is_tnk] > 0).mean()) if is_tnk.any() else None,
            "ptprc_pos_malignant": float((ptprc[idx][is_malig] > 0).mean()) if is_malig.any() else None,
            "ptprc_pos_tnk": float((ptprc[idx][is_tnk] > 0).mean()) if is_tnk.any() else None,
            "cldn4_pos_malignant": float((cldn4_k[is_malig] > 0).mean()) if is_malig.any() else None,
            "cldn4_pos_tnk": float((cldn4_k[is_tnk] > 0).mean()) if is_tnk.any() else None,
        },
        "da_malignant_cldn4": fdr_counts(da_cldn4),
        "da_malignant_cldn4_nsclc_only": fdr_counts(da_nsclc),
        "da_recist_r_vs_nr": fdr_counts(da_recist),
        "patient_malignant_cldn4_vs_tnk_fraction": patient_sp,
        "patient_malignant_cldn4_recist": recist_cldn4,
        "composition": {
            "note": "transcriptional kNN != spatial niche; unrestricted CLDN4 vs T/NK is partly lineage geometry",
            "all_nhoods_maligCLDN4_vs_fracTNK": all_sp,
            "interface_nhoods": {
                "definition": f">={args.min_interface} malignant AND >={args.min_interface} T/NK cells",
                "n": int(iface.sum()),
                "spearman": iface_sp,
            },
            "disjoint_interface": {
                "n_disjoint_nhoods": int(indep.size),
                "n_disjoint_interface": int((indep_mask & iface.to_numpy()).sum()),
                "spearman": indep_sp,
            },
            "sample_paired_TNK_in_CLDN4high_vs_low_nhoods": {
                "high_low_cut": "median malignant CLDN4 among nhoods with >=5 malignant cells",
                "median_cut": med,
                "n_high_nhoods": int(high.sum()),
                "n_low_nhoods": int(low.sum()),
                "wilcoxon_high_vs_low": paired_test,
                "direction_expected": "frac_tnk_high < frac_tnk_low (immune-poor CLDN4-high nhoods)",
            },
            "nhood_cldn4_DA_rho_vs_fracTNK": da_vs_comp,
        },
        "honest_n": {
            "independent_unit_DA": "patient, not cell, not neighbourhood",
            "n_patients_in_graph": int(len(sample_levels)),
            "n_patients_with_malignant_CLDN4_score": n_with_score,
            "n_RECIST_R_vs_NR": f"{n_r} vs {n_nr}",
            "n_NE": n_ne,
            "n_NSCLC_ADC_SQ_with_score": n_nsclc,
            "MPR_labeled": False,
            "nhoods_are_overlapping": True,
            "SpatialFDR_is_overlap_aware": True,
            "disjoint_nhood_subset_n": int(indep.size),
            "do_not_cite_n_cells_as_n": int(n_k),
        },
        "mpr_gate": {
            "mpr_labeled": False,
            "endpoint_available": "RECIST 1.1 only (PR/SD/PD/NE)",
            "do_not_substitute_recist_for_mpr": True,
        },
        "provenance": provenance,
    }
    write_json(tables / "summary.json", summary)

    honest = pd.DataFrame(
        [
            {"item": "GEO patients", "n": int(meta["patient"].nunique()), "note": "26 patients / 33 GSM"},
            {"item": "GEO samples", "n": int(len(meta)), "note": "not the DA unit"},
            {"item": "cells in processed matrix", "n": int(len(cells)), "note": "do not cite as n"},
            {"item": "normal-tissue samples excluded from graph", "n": int(dropped_normal["orig.ident"].nunique()), "note": "Normal LN / Normal Brain / Normal Lung"},
            {"item": "patients in graph", "n": int(len(sample_levels)), "note": "non-normal tissues; biopsies pooled"},
            {"item": "patients with malignant CLDN4 score (>=10 mal. cells)", "n": n_with_score, "note": "primary Spearman DA n"},
            {"item": "RECIST R (PR)", "n": n_r, "note": "secondary; not MPR"},
            {"item": "RECIST NR (SD+PD)", "n": n_nr, "note": "secondary; not NMPR"},
            {"item": "RECIST NE", "n": n_ne, "note": "kept in continuous CLDN4 DA; dropped from R vs NR"},
            {"item": "NSCLC ADC+SQ with CLDN4 score", "n": n_nsclc, "note": "sensitivity; drops SCLC/NUT"},
            {"item": "neighbourhoods", "n": int(len(members)), "note": "overlapping; not independent"},
            {"item": "testable nhoods vs malignant CLDN4", "n": int(da_cldn4["testable"].sum()), "note": ">=5 patients present"},
            {"item": "disjoint nhood subset", "n": int(indep.size), "note": "greedy zero-overlap"},
            {"item": "MPR/NMPR labelled patients", "n": 0, "note": "not a GEO or identity-table field"},
        ]
    )
    honest.to_csv(tables / "honest_n.tsv", sep="\t", index=False)

    plot_volcano(
        da_cldn4.rename(columns={"spearman_rho": "effect"}),
        "effect",
        "GSE205335 neighbourhood DA vs malignant CLDN4",
        figures / "fig_da_cldn4_volcano.png",
        "Spearman ρ (nhood abundance vs patient malignant CLDN4)",
    )
    plot_volcano(
        da_recist,
        "logFC_B_minus_A",
        "GSE205335 neighbourhood DA: RECIST NR vs R (not MPR)",
        figures / "fig_da_recist_volcano.png",
        "log2 FC (NR − R patient proportion)",
    )
    plot_scatter(
        comp.loc[iface, "cldn4_malig_mean"],
        comp.loc[iface, "frac_tnk"],
        figures / "fig_interface_cldn4_vs_tnk.png",
        "Neighbourhood malignant CLDN4 (mean log1p CP10k)",
        "Neighbourhood T/NK fraction",
        f"Interface nhoods (n={int(iface.sum())}); transcriptional, not spatial",
    )
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    ok = paired["malignant_cldn4"].notna() & paired["sample_frac_tnk"].notna()
    colors = {"R": "#4c72b0", "NR": "#c44e52", "NE": "#999999"}
    for lab, sub in paired.loc[ok].groupby("response"):
        ax.scatter(
            sub["malignant_cldn4"],
            sub["sample_frac_tnk"],
            s=36,
            c=colors.get(str(lab), "#4c72b0"),
            label=f"{lab} n={len(sub)}",
            edgecolors="white",
            linewidths=0.4,
        )
    ax.set_xlabel("Patient malignant CLDN4 (mean log1p CP10k)")
    ax.set_ylabel("T/NK fraction among (malignant + T/NK)")
    ax.set_title(f"Patient-level (n={int(ok.sum())}); unit = patient")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figures / "fig_patient_cldn4_vs_tnk.png", dpi=160)
    plt.close(fig)

    if paired["frac_tnk_high_nhood"].notna().any():
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        for _, r in paired.iterrows():
            if np.isfinite(r["frac_tnk_high_nhood"]) and np.isfinite(r["frac_tnk_low_nhood"]):
                c = colors.get(str(r["response"]), "#4c72b0")
                ax.plot(
                    [0, 1],
                    [r["frac_tnk_low_nhood"], r["frac_tnk_high_nhood"]],
                    "-o",
                    c=c,
                    alpha=0.75,
                    ms=5,
                )
        ax.set_xticks([0, 1], ["CLDN4-low nhoods", "CLDN4-high nhoods"])
        ax.set_ylabel("Patient T/NK fraction among cells in those nhoods")
        ax.set_title("GSE205335 patient-paired (unit = patient)")
        fig.tight_layout()
        fig.savefig(figures / "fig_patient_paired_tnk.png", dpi=160)
        plt.close(fig)

    print(json.dumps({
        "n_nhoods": summary["graph"]["n_nhoods"],
        "patients": summary["patients"],
        "da_cldn4": summary["da_malignant_cldn4"],
        "da_recist": summary["da_recist_r_vs_nr"],
        "iface": summary["composition"]["interface_nhoods"],
        "paired": summary["composition"]["sample_paired_TNK_in_CLDN4high_vs_low_nhoods"],
        "patient_sp": patient_sp,
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
