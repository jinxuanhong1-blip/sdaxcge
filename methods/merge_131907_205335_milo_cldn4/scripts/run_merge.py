#!/usr/bin/env python3
"""Merged GSE131907 + GSE205335 Milo-style DA vs malignant CLDN4.

Additive CLDN4-only. PR #320 T/NK ρ is given and is not re-audited.
GSE207422 is not run. miloR / edgeR are not used.

Primary graphs are per-dataset (sample/patient PCA centering).
Harmony is an extra joint graph on a documented per-unit cell cap.
"""

from __future__ import annotations

import argparse
import gc
import gzip
import json
import shutil
import sys
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from knn_nhood import (  # noqa: E402
    EPS,
    attach_fdr,
    count_matrix,
    greedy_independent,
    knn_from_embedding,
    make_nhoods,
    nhood_composition,
    nhood_membership,
    pca_knn,
    refine_indices,
    sample_paired_tnk_by_gene,
    sample_sizes,
    select_hvg,
    self_test,
    spearman_da,
    spearman_safe,
    welch_da,
    wilcoxon_paired,
    write_json,
)

MALIGNANT_131907 = {"Malignant cells", "tS1", "tS2", "tS3"}
KEEP_131907 = {
    "T lymphocytes",
    "NK cells",
    "B lymphocytes",
    "Myeloid cells",
    "MAST cells",
    "Epithelial cells",
}
TNK_131907 = {"T lymphocytes", "NK cells"}
TUMOR_NO_BRAIN = {"tLung", "tL/B", "mLN"}
MARKERS = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD3D", "CD3E", "NKG7", "KRT8", "KRT19"]
ALWAYS = ["CLDN4", "EPCAM", "PTPRC", "CD3D", "NKG7"]
RESPONSE_MAP = {"PR": "R", "CR": "R", "SD": "NR", "PD": "NR", "NE": "NE"}
NORMAL_TISSUE_PREFIX = "Normal"
GIVEN_PR320 = {
    "source": "PR #320 methods/cldn4_malig_q4_tnk (not re-audited)",
    "combo": "author %pos vs T/NK: GSE131907 + GSE205335",
    "k": 2,
    "N": 43,
    "rho": -0.4787109080847739,
    "p": 0.0015191063422345722,
    "I2": 0.0,
}


def fmt(x, nd=3, sci=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if sci:
        return f"{x:.2e}"
    return f"{x:.{nd}f}"


def fdr_counts(df: pd.DataFrame) -> dict:
    t = df[df["testable"]] if len(df) and "testable" in df.columns else df.iloc[0:0]
    rec = {
        "n_nhoods_total": int(len(df)),
        "n_testable": int(len(t)),
        "n_p_lt_0.05": int((t["p"] < 0.05).sum()) if len(t) else 0,
        "n_BH_FDR_lt_0.1": int((t["BH_FDR"] < 0.1).sum()) if len(t) else 0,
        "n_BH_FDR_lt_0.05": int((t["BH_FDR"] < 0.05).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.2": int((t["SpatialFDR"] < 0.2).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.1": int((t["SpatialFDR"] < 0.1).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.05": int((t["SpatialFDR"] < 0.05).sum()) if len(t) else 0,
        "min_p": float(t["p"].min()) if len(t) else None,
        "min_SpatialFDR": float(t["SpatialFDR"].min()) if len(t) else None,
        "min_BH_FDR": float(t["BH_FDR"].min()) if len(t) else None,
    }
    return rec


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
    ax.set_ylabel(r"−log10 p (sample/patient-level test)")
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


def gunzip_until_rds(src: Path, dest: Path) -> Path:
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
    print(f"matrix {matrix.shape} nnz={matrix.nnz} load={time.time() - t0:.1f}s", flush=True)
    return matrix, genes, barcodes


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


def gene_mean_var(matrix: sparse.spmatrix) -> tuple[np.ndarray, np.ndarray]:
    n = matrix.shape[1]
    mean = np.asarray(matrix.mean(axis=1)).ravel()
    mean_sq = np.asarray(matrix.multiply(matrix).mean(axis=1)).ravel()
    var = np.maximum(mean_sq - mean * mean, 0.0) * (n / max(n - 1, 1))
    return mean, var


def stream_pass1(matrix_path: Path, markers: set[str], cohort_masks: dict[str, np.ndarray]):
    t0 = time.time()
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        names: list[str] = []
        stats = {c: {"mean": [], "var": []} for c in cohort_masks}
        found: dict[str, np.ndarray] = {}
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} values, expected {n}")
            n_umi += arr
            names.append(gene)
            for c, mask in cohort_masks.items():
                sub = arr[mask]
                stats[c]["mean"].append(float(sub.mean()) if sub.size else 0.0)
                stats[c]["var"].append(float(sub.var()) if sub.size else 0.0)
            if gene in markers:
                found[gene] = arr
            if i % 2000 == 0:
                print(f"  pass1 {i} genes, markers={len(found)} [{time.time() - t0:.0f}s]", flush=True)
    out_stats = {c: (np.asarray(stats[c]["mean"]), np.asarray(stats[c]["var"])) for c in stats}
    return cell_ids, n_umi, names, out_stats, found


def stream_pass2(matrix_path: Path, keep: set[str], keep_idx: np.ndarray, n_cells: int) -> dict[str, np.ndarray]:
    t0 = time.time()
    out: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        handle.readline()
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            if gene not in keep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_cells}")
            out[gene] = arr[keep_idx].astype(np.float32, copy=False)
            if len(out) % 200 == 0:
                print(f"  pass2 stored {len(out)}/{len(keep)} [{time.time() - t0:.0f}s]", flush=True)
    missing = keep - set(out)
    if missing:
        print(f"  missing {sorted(missing)[:12]} ({len(missing)} genes)", flush=True)
    return out


def cap_cells(sample: np.ndarray, max_per_unit: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    keep = np.zeros(sample.size, dtype=bool)
    for s in np.unique(sample):
        idx = np.flatnonzero(sample == s)
        if idx.size <= max_per_unit:
            keep[idx] = True
        else:
            keep[rng.choice(idx, size=max_per_unit, replace=False)] = True
    return keep


def run_graph(
    name: str,
    X: np.ndarray,
    sample: np.ndarray,
    lineage: np.ndarray,
    is_malig: np.ndarray,
    is_tnk: np.ndarray,
    cldn: np.ndarray,
    extra_meta: pd.DataFrame | None,
    outdir: Path,
    figdir: Path,
    k: int,
    d: int,
    prop: float,
    n_hvg: int,
    min_interface: int,
    seed: int,
    min_malig: int,
    batch_note: str,
    center_by_sample: bool = True,
    embedding: np.ndarray | None = None,
) -> dict:
    print(f"{name}: cells={X.shape[0]} HVG={X.shape[1]} units={pd.unique(sample).size}", flush=True)
    if embedding is None:
        pcs, knn_idx, knn_dist = pca_knn(
            X, sample, n_pcs=d, k=k, random_state=seed, center_by_sample=center_by_sample
        )
    else:
        pcs = embedding
        knn_idx, knn_dist = knn_from_embedding(pcs, k=k)
    indices = refine_indices(pcs, knn_idx, prop=prop, random_state=seed)
    nhoods = make_nhoods(knn_idx, knn_dist, indices)
    members = nhood_membership(nhoods)
    print(f"{name}: nhoods={len(members)} median_size={np.median([m.size for m in members]):.0f}", flush=True)

    sample_levels = sorted(map(str, pd.unique(sample)))
    counts = count_matrix(members, sample.astype(str), sample_levels)
    sizes = sample_sizes(sample.astype(str), sample_levels)

    mal_score = []
    tnk_frac = []
    n_mal_s = []
    n_tnk_s = []
    for s in sample_levels:
        m = sample.astype(str) == s
        n_mal = int((m & is_malig).sum())
        n_mal_s.append(n_mal)
        n_tnk_s.append(int((m & is_tnk).sum()))
        mal_score.append(float(cldn[m & is_malig].mean()) if n_mal >= min_malig else np.nan)
        tnk_frac.append(float(is_tnk[m].mean()))
    mal_score = np.asarray(mal_score, dtype=float)
    tnk_frac = np.asarray(tnk_frac, dtype=float)

    da_cldn = attach_fdr(spearman_da(counts, sizes, mal_score, min_samples=5), nhoods.k_distance)
    scored = np.isfinite(mal_score)
    group = np.full(len(sample_levels), "unscored", dtype=object)
    if scored.sum() >= 4:
        med = float(np.median(mal_score[scored]))
        group[scored & (mal_score > med)] = "high"
        group[scored & (mal_score < med)] = "low"
    n_high = int((group == "high").sum())
    n_low = int((group == "low").sum())
    if n_high >= 2 and n_low >= 2:
        da_split = attach_fdr(
            welch_da(counts, sizes, group, group_a="low", group_b="high", min_per_group=2, min_samples=4),
            nhoods.k_distance,
        )
    else:
        da_split = pd.DataFrame(
            {
                "nhood": np.arange(len(members)),
                "p": np.nan,
                "testable": False,
                "logFC_B_minus_A": np.nan,
                "BH_FDR": np.nan,
                "SpatialFDR": np.nan,
            }
        )

    comp = nhood_composition(members, is_tnk, is_malig, cldn, lineage, gene_prefix="cldn4")
    merged = comp.merge(
        da_cldn.rename(
            columns={
                "spearman_rho": "cldn4_rho",
                "p": "cldn4_p",
                "testable": "cldn4_testable",
                "BH_FDR": "cldn4_BH_FDR",
                "SpatialFDR": "cldn4_SpatialFDR",
                "n_cells": "n_cells_da",
            }
        ),
        on="nhood",
        how="left",
    )
    merged.insert(0, "graph", name)
    dest = outdir / name
    dest.mkdir(parents=True, exist_ok=True)
    merged.to_csv(dest / "nhoods.tsv", sep="\t", index=False)
    da_cldn.to_csv(dest / "da_malignant_cldn4.tsv", sep="\t", index=False)
    da_split.to_csv(dest / "da_cldn4_median_split.tsv", sep="\t", index=False)

    iface = (comp["n_malig"] >= min_interface) & (comp["n_tnk"] >= min_interface)
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
    med_n = float(np.nanmedian(comp.loc[has_mal, "cldn4_malig_mean"])) if has_mal.any() else np.nan
    high = has_mal.to_numpy() & (comp["cldn4_malig_mean"].to_numpy() >= med_n)
    low = has_mal.to_numpy() & (comp["cldn4_malig_mean"].to_numpy() < med_n)
    paired = sample_paired_tnk_by_gene(members, sample.astype(str), is_tnk, high, low)
    sample_tab = pd.DataFrame(
        {
            "unit": sample_levels,
            "graph": name,
            "n_cells": sizes.astype(int),
            "n_malignant": n_mal_s,
            "n_tnk": n_tnk_s,
            "malignant_cldn4": mal_score,
            "frac_tnk": tnk_frac,
            "cldn4_arm": group,
        }
    )
    if extra_meta is not None and len(extra_meta):
        sample_tab = sample_tab.merge(extra_meta, left_on="unit", right_on="unit", how="left")
    paired = paired.rename(columns={"sample": "unit"}).merge(sample_tab, on="unit")
    paired.to_csv(dest / "sample_paired_tnk_by_nhood_cldn4.tsv", sep="\t", index=False)
    sample_tab.to_csv(dest / "sample_scores.tsv", sep="\t", index=False)
    paired_test = wilcoxon_paired(paired["frac_tnk_high_nhood"], paired["frac_tnk_low_nhood"])
    sample_level = spearman_safe(mal_score, tnk_frac)

    plot_volcano(
        da_cldn.rename(columns={"spearman_rho": "stat"}),
        "stat",
        f"{name}: nhood abundance vs malignant CLDN4",
        figdir / f"fig_da_cldn4_volcano_{name}.png",
        "Spearman ρ (nhood proportion vs unit malignant CLDN4)",
    )
    if da_split["testable"].any():
        plot_volcano(
            da_split,
            "logFC_B_minus_A",
            f"{name}: nhood DA CLDN4-high vs low",
            figdir / f"fig_da_cldn4_split_volcano_{name}.png",
            "log2 FC (high − low unit proportion)",
        )
    plot_scatter(
        comp.loc[iface, "cldn4_malig_mean"],
        comp.loc[iface, "frac_tnk"],
        figdir / f"fig_interface_cldn4_vs_tnk_{name}.png",
        "Neighbourhood malignant CLDN4 (mean log1p CP10k)",
        "Neighbourhood T/NK fraction",
        f"{name} interface nhoods (n={int(iface.sum())}); transcriptional, not PR #320",
    )
    if paired["frac_tnk_high_nhood"].notna().any():
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        for _, r in paired.iterrows():
            if np.isfinite(r["frac_tnk_high_nhood"]) and np.isfinite(r["frac_tnk_low_nhood"]):
                ax.plot(
                    [0, 1],
                    [r["frac_tnk_low_nhood"], r["frac_tnk_high_nhood"]],
                    "-o",
                    c="#4c72b0",
                    alpha=0.75,
                    ms=5,
                )
        ax.set_xticks([0, 1], ["CLDN4-low nhoods", "CLDN4-high nhoods"])
        ax.set_ylabel("Unit T/NK fraction among cells in those nhoods")
        ax.set_title(f"{name} paired (unit = sample/patient; not PR #320)")
        fig.tight_layout()
        fig.savefig(figdir / f"fig_sample_paired_tnk_{name}.png", dpi=160)
        plt.close(fig)

    dropped = sample_tab.loc[~np.isfinite(sample_tab["malignant_cldn4"])]
    drop_bits = [f"{r.unit}={int(r.n_malignant)}" for r in dropped.itertuples()]
    summary = {
        "graph": name,
        "n_cells_in_graph": int(X.shape[0]),
        "n_malignant": int(is_malig.sum()),
        "n_tnk": int(is_tnk.sum()),
        "lineage_counts": {k: int(v) for k, v in pd.Series(lineage).value_counts().items()},
        "graph_params": {
            "k": k,
            "d": min(d, X.shape[1], int(X.shape[0]) - 1) if embedding is None else int(pcs.shape[1]),
            "prop": prop,
            "n_hvg": int(X.shape[1]),
            "n_hvg_requested": n_hvg,
            "n_nhoods": int(len(members)),
            "nhood_size": k + 1,
            "median_nhood_size": float(np.median([m.size for m in members])),
            "batch": batch_note,
        },
        "units": {
            "n": int(len(sample_levels)),
            "n_with_malignant_ge10": int(np.isfinite(mal_score).sum()),
            "n_cldn4_high": n_high,
            "n_cldn4_low": n_low,
            "dropped_lt10_malignant": drop_bits,
            "median_split_cut": float(np.median(mal_score[scored])) if scored.any() else None,
        },
        "unit_level_malignant_CLDN4_vs_TNK_fraction": sample_level,
        "note_tnk": "descriptive only; PR #320 T/NK ρ is given and is not re-audited",
        "da_malignant_cldn4": fdr_counts(da_cldn),
        "da_cldn4_median_split": fdr_counts(da_split),
        "composition": {
            "note": "transcriptional kNN != spatial niche; not the PR #320 T/NK audit",
            "all_nhoods_maligCLDN4_vs_fracTNK": all_sp,
            "interface_nhoods": {
                "definition": f">={min_interface} malignant AND >={min_interface} T/NK cells",
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
                "median_cut": med_n,
                "n_high_nhoods": int(high.sum()),
                "n_low_nhoods": int(low.sum()),
                "wilcoxon_high_vs_low": paired_test,
            },
        },
        "honest_n": {
            "independent_unit_DA": f"{name} units (n={len(sample_levels)}; scored={int(np.isfinite(mal_score).sum())})",
            "nhoods_are_overlapping": True,
            "SpatialFDR_is_overlap_aware": True,
            "disjoint_nhood_subset_n": int(indep.size),
            "do_not_cite_n_cells_as_n": int(X.shape[0]),
        },
    }
    write_json(dest / "summary.json", summary)
    return {"summary": summary, "nhoods": merged, "scores": sample_tab, "da": da_cldn}


def load_gse205335(datadir: Path, n_hvg: int) -> dict:
    ident_path = datadir / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = datadir / "GSE205335_family.soft.gz"
    matrix_gz = datadir / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    for p in (ident_path, soft_path, matrix_gz):
        if not p.exists():
            raise SystemExit(f"missing {p}; run scripts/download.py")
    rds_path = gunzip_until_rds(matrix_gz, datadir / "GSE205335_Lung_IO_UMI_matrix.rds")
    meta = parse_geo_soft(soft_path)
    ident = pd.read_csv(ident_path, sep="\t")
    matrix, genes, barcodes = load_dgcmatrix(rds_path)
    if len(barcodes) != len(ident) or not np.array_equal(barcodes, ident["barcode"].to_numpy(dtype=str)):
        order = pd.Index(ident["barcode"]).get_indexer(barcodes)
        if (order < 0).any():
            raise ValueError("identity barcodes do not match matrix barcodes")
        ident = ident.iloc[order].reset_index(drop=True)
    cells = ident.merge(
        meta[["orig.ident", "gsm", "patient", "tissue", "recist", "platform", "cancer_subtype", "tumor_stage"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        missing = cells.loc[cells["patient"].isna(), "orig.ident"].unique().tolist()
        raise ValueError(f"identity samples missing GEO metadata: {missing}")
    cells["is_normal_tissue"] = cells["tissue"].astype(str).str.startswith(NORMAL_TISSUE_PREFIX)
    cells["is_malig"] = cells["lineage.sub"].eq("Malignant cells")
    cells["is_tnk"] = cells["lineage.total"].eq("T/NK cells")
    gene_idx = {g: i for i, g in enumerate(genes)}
    for g in ALWAYS:
        if g not in gene_idx:
            raise ValueError(f"{g} missing from GSE205335 UMI matrix")
    keep_cells = ~cells["is_normal_tissue"].to_numpy()
    print(f"GSE205335 non-normal {int(keep_cells.sum())} / {len(cells)}", flush=True)
    matrix = matrix[:, keep_cells]
    cells = cells.loc[keep_cells].reset_index(drop=True)
    gene_mean, gene_var = gene_mean_var(matrix)
    hvg_idx = select_hvg(gene_mean, gene_var, n_hvg=n_hvg)
    hvg_names = [str(genes[i]) for i in hvg_idx]
    keep_genes = sorted(set(hvg_idx.tolist()) | {gene_idx[g] for g in ALWAYS if g in gene_idx})
    lib = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    scale = np.where(lib > 0, 1e4 / lib, 0.0)
    sub = matrix[keep_genes, :].astype(np.float32)
    del matrix
    dense = sub.toarray()
    del sub
    log_all = np.log1p(dense * scale).astype(np.float32)
    del dense
    name_of = {keep_genes[i]: i for i in range(len(keep_genes))}
    cldn4 = log_all[name_of[gene_idx["CLDN4"]]]
    hvg_rows = [name_of[i] for i in hvg_idx if i in name_of]
    X = log_all[hvg_rows].T.copy()
    gene_log = {str(genes[i]): log_all[name_of[i]] for i in keep_genes}
    del log_all
    gc.collect()
    return {
        "dataset": "GSE205335",
        "X": X,
        "hvg": hvg_names,
        "gene_log": gene_log,
        "sample": cells["patient"].to_numpy(dtype=object),
        "lineage": cells["lineage.total"].to_numpy(dtype=object),
        "is_malig": cells["is_malig"].to_numpy(),
        "is_tnk": cells["is_tnk"].to_numpy(),
        "cldn4": cldn4.astype(np.float32),
        "n_cells_matrix": 96505,
        "n_genes_matrix": int(len(genes)),
        "matrix_bytes": int(matrix_gz.stat().st_size),
    }


def load_gse131907(datadir: Path, n_hvg: int) -> dict:
    ann_path = datadir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = datadir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    for p in (ann_path, umi_path):
        if not p.exists():
            raise SystemExit(f"missing {p}; run scripts/download.py")
    with gzip.open(umi_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    cell_ids = header[1:]
    print(f"GSE131907 matrix header cells={len(cell_ids)}", flush=True)
    ann = pd.read_csv(ann_path, sep="\t", dtype=str)
    if "Index" not in ann.columns:
        raise SystemExit(f"annotation missing Index: {list(ann.columns)}")
    ann = ann.set_index("Index").reindex(cell_ids).reset_index()
    if ann["Sample"].isna().any():
        raise SystemExit(f"{int(ann['Sample'].isna().sum())} matrix cell IDs missing from annotation")
    cell_type = ann["Cell_type"].to_numpy()
    subtype = ann["Cell_subtype"].fillna("").to_numpy()
    origin = ann["Sample_Origin"].fillna("").to_numpy()
    sample = ann["Sample"].to_numpy()
    keep_type = np.isin(cell_type, list(KEEP_131907))
    is_malig_all = np.isin(subtype, list(MALIGNANT_131907))
    is_tnk_all = np.isin(cell_type, list(TNK_131907))
    masks = {
        "tLung": (origin == "tLung") & keep_type,
        "tumor_no_brain": np.isin(origin, list(TUMOR_NO_BRAIN)) & keep_type,
    }
    for c, m in masks.items():
        print(
            f"  {c}: graph_cells={int(m.sum())} malignant={int((m & is_malig_all).sum())} "
            f"samples={pd.unique(sample[m]).size}",
            flush=True,
        )
    print("GSE131907 pass1: gene stats + markers", flush=True)
    cell_ids2, n_umi, gene_names, cohort_stats, marker_expr = stream_pass1(umi_path, set(MARKERS), masks)
    if cell_ids2 != cell_ids:
        raise SystemExit("matrix header changed between reads")
    if "CLDN4" not in marker_expr:
        raise SystemExit("CLDN4 is absent from the GSE131907 UMI matrix")
    hvg_names = {}
    keep_genes = set(MARKERS)
    for c in masks:
        mu, va = cohort_stats[c]
        idx = select_hvg(mu, va, n_hvg=n_hvg)
        hvg_names[c] = [gene_names[i] for i in idx]
        keep_genes.update(hvg_names[c])
        print(f"  {c} HVG={len(hvg_names[c])}", flush=True)
    keep_any = masks["tumor_no_brain"]
    keep_idx = np.flatnonzero(keep_any)
    print(f"GSE131907 pass2: {len(keep_genes)} genes × {keep_idx.size} kept cells", flush=True)
    expr = stream_pass2(umi_path, keep_genes, keep_idx, len(cell_ids))
    for g, arr in marker_expr.items():
        expr.setdefault(g, arr[keep_idx].astype(np.float32, copy=False))
    n_umi_k = n_umi[keep_idx]
    scale_k = np.where(n_umi_k > 0, 1e4 / n_umi_k, 0.0)
    gene_log = {g: np.log1p(arr * scale_k).astype(np.float32) for g, arr in expr.items()}
    cldn_k = gene_log["CLDN4"]
    return {
        "dataset": "GSE131907",
        "expr_log": gene_log,
        "hvg": hvg_names,
        "sample": sample[keep_idx],
        "origin": origin[keep_idx],
        "lineage": cell_type[keep_idx],
        "is_malig": is_malig_all[keep_idx],
        "is_tnk": is_tnk_all[keep_idx],
        "cldn4": cldn_k,
        "local_tLung": origin[keep_idx] == "tLung",
        "n_cells_matrix": int(len(cell_ids)),
        "n_genes_matrix": int(len(gene_names)),
        "matrix_bytes": int(umi_path.stat().st_size),
    }


def matrix_from_hvg(gene_log: dict[str, np.ndarray], hvg: list[str], local: np.ndarray | None = None) -> np.ndarray:
    present = [g for g in hvg if g in gene_log]
    mats = []
    for g in present:
        arr = gene_log[g]
        mats.append(arr[local] if local is not None else arr)
    return np.vstack(mats).T.astype(np.float32)


def try_harmony(
    packs: list[dict],
    outdir: Path,
    figdir: Path,
    args,
) -> dict | None:
    """Joint Harmony graph on a per-unit cell cap. Returns None if it cannot run."""
    try:
        import harmonypy as hm
    except ImportError:
        print("harmonypy missing; skip joint graph", flush=True)
        return {
            "ran": False,
            "reason": "harmonypy not importable",
        }
    shared = set(packs[0]["hvg_used"])
    for p in packs[1:]:
        shared &= set(p["hvg_used"])
    shared = sorted(g for g in shared if all(g in p["gene_log"] for p in packs))
    print(f"Harmony shared HVG={len(shared)}", flush=True)
    if len(shared) < 500:
        return {"ran": False, "reason": f"shared HVG only {len(shared)} (<500)"}
    Xs = []
    samples = []
    lineages = []
    maligs = []
    tnks = []
    cldns = []
    datasets = []
    n_before = []
    n_after = []
    for p in packs:
        local = p.get("local")
        sample = p["sample"] if local is None else p["sample"][local]
        cap = cap_cells(sample.astype(str), args.harmony_max_per_unit, args.seed)
        n_before.append(int(sample.size))
        n_after.append(int(cap.sum()))
        X = matrix_from_hvg(p["gene_log"], shared, None if local is None else local)[cap]
        Xs.append(X)
        samples.append(np.array([f"{p['dataset']}:{s}" for s in sample.astype(str)[cap]], dtype=object))
        lin = p["lineage"] if local is None else p["lineage"][local]
        lineages.append(lin[cap])
        mal = p["is_malig"] if local is None else p["is_malig"][local]
        maligs.append(mal[cap])
        tnk = p["is_tnk"] if local is None else p["is_tnk"][local]
        tnks.append(tnk[cap])
        cld = p["cldn4"] if local is None else p["cldn4"][local]
        cldns.append(cld[cap])
        datasets.append(np.array([p["dataset"]] * int(cap.sum()), dtype=object))
    X = np.vstack(Xs)
    sample = np.concatenate(samples)
    lineage = np.concatenate(lineages)
    is_malig = np.concatenate(maligs)
    is_tnk = np.concatenate(tnks)
    cldn = np.concatenate(cldns)
    dataset = np.concatenate(datasets)
    print(f"Harmony cells={X.shape[0]} genes={X.shape[1]} (capped {n_after} from {n_before})", flush=True)
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler

        z = StandardScaler().fit_transform(X)
        n_pcs = min(args.d, z.shape[0] - 1, z.shape[1])
        pcs = PCA(n_components=n_pcs, svd_solver="randomized", random_state=args.seed).fit_transform(z)
        meta = pd.DataFrame({"dataset": dataset, "unit": sample})
        ho = hm.run_harmony(pcs, meta, "dataset", max_iter_harmony=20)
        Z = np.asarray(ho.Z_corr).T.astype(np.float32)
        fig, ax = plt.subplots(figsize=(5.6, 4.6))
        for ds, col in (("GSE131907", "#4c72b0"), ("GSE205335", "#c44e52")):
            m = dataset == ds
            ax.scatter(Z[m, 0], Z[m, 1], s=3, c=col, alpha=0.25, label=ds)
        ax.set_xlabel("Harmony PC1")
        ax.set_ylabel("Harmony PC2")
        ax.set_title(f"Harmony-aligned (cap {args.harmony_max_per_unit}/unit)")
        ax.legend(frameon=False, markerscale=3)
        fig.tight_layout()
        fig.savefig(figdir / "fig_harmony_pc12.png", dpi=160)
        plt.close(fig)
        result = run_graph(
            name="harmony_tumor_no_brain_plus_GSE205335",
            X=X,
            sample=sample,
            lineage=lineage,
            is_malig=is_malig,
            is_tnk=is_tnk,
            cldn=cldn,
            extra_meta=None,
            outdir=outdir,
            figdir=figdir,
            k=args.k,
            d=args.d,
            prop=args.prop,
            n_hvg=len(shared),
            min_interface=args.min_interface,
            seed=args.seed,
            min_malig=args.min_malignant,
            batch_note=f"Harmony batch=dataset; cap {args.harmony_max_per_unit} cells/unit; shared HVG={len(shared)}",
            center_by_sample=False,
            embedding=Z,
        )
        result["summary"]["harmony"] = {
            "ran": True,
            "shared_hvg": len(shared),
            "n_cells_before_cap": n_before,
            "n_cells_after_cap": n_after,
            "max_per_unit": args.harmony_max_per_unit,
            "n_units": int(pd.unique(sample).size),
        }
        return result
    except Exception as exc:  # noqa: BLE001 — Harmony is optional
        print(f"Harmony failed: {exc!r}", flush=True)
        return {"ran": False, "reason": repr(exc)}


def write_finding(overall: dict, out_path: Path) -> None:
    graphs = overall["graphs"]

    def block(tag: str, s: dict) -> str:
        da = s["da_malignant_cldn4"]
        sp = s["da_cldn4_median_split"]
        c = s["composition"]
        iface = c["interface_nhoods"]["spearman"]
        paired = c["sample_paired_TNK_in_CLDN4high_vs_low_nhoods"]
        w = paired["wilcoxon_high_vs_low"]
        sl = s["unit_level_malignant_CLDN4_vs_TNK_fraction"]
        drop = ", ".join(s["units"]["dropped_lt10_malignant"]) or "none"
        return f"""### {tag}

| Contrast | Arm | n units | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | units with ≥10 author-malignant cells | **{s['units']['n_with_malignant_ge10']}** of {s['units']['n']} | Mean log1p-CP10k CLDN4 in author malignant cells. Dropped if <10 ({drop}). |
| Median split (Welch, secondary) | high / low | **{s['units']['n_cldn4_high']} vs {s['units']['n_cldn4_low']}** | Median of the scored units. Ties at the median unlabeled. |

Do not cite {s['n_cells_in_graph']:,} graph cells or {s['graph_params']['n_nhoods']:,} neighbourhoods as *n*. Neighbourhoods overlap. A disjoint subset has n={s['honest_n']['disjoint_nhood_subset_n']}.

Graph: {s['n_cells_in_graph']:,} epithelium+immune cells; {s['n_malignant']:,} author-malignant; {s['n_tnk']:,} T/NK; {s['graph_params']['n_nhoods']:,} neighbourhoods of size {s['graph_params']['nhood_size']}; k={s['graph_params']['k']}, d={s['graph_params']['d']}, {s['graph_params']['n_hvg']} HVG. Batch: {s['graph_params']['batch']}.

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | {da['n_testable']} | {da['n_p_lt_0.05']} | {fmt(da.get('min_p'), sci=True)} | {fmt(da.get('min_SpatialFDR'))} | **{da['n_SpatialFDR_lt_0.1']}** | {da['n_SpatialFDR_lt_0.05']} | {da['n_BH_FDR_lt_0.1']} |
| Median split high vs low | {sp['n_testable']} | {sp['n_p_lt_0.05']} | {fmt(sp.get('min_p'), sci=True)} | {fmt(sp.get('min_SpatialFDR'))} | **{sp['n_SpatialFDR_lt_0.1']}** | {sp['n_SpatialFDR_lt_0.05']} | {sp['n_BH_FDR_lt_0.1']} |

Unit-level malignant CLDN4 vs T/NK fraction (descriptive; **not** the PR #320 audit): ρ={fmt(sl.get('rho'))}, p={fmt(sl.get('p'), sci=True)}, n={sl.get('n')}.

Interface neighbourhoods (≥3 malignant and ≥3 T/NK): n={c['interface_nhoods']['n']}; Spearman neighbourhood malignant CLDN4 vs T/NK fraction ρ={fmt(iface.get('rho'))}, p={fmt(iface.get('p'), sci=True)}, n={iface.get('n')}. Transcriptional kNN, not histology.

Disjoint interface subset: n={c['disjoint_interface']['n_disjoint_interface']} of {c['disjoint_interface']['n_disjoint_nhoods']} disjoint neighbourhoods; ρ={fmt(c['disjoint_interface']['spearman'].get('rho'))}, p={fmt(c['disjoint_interface']['spearman'].get('p'), sci=True)}.

Sample/patient-paired T/NK in CLDN4-high vs CLDN4-low neighbourhoods (unit = sample/patient; extra, not PR #320): n={w.get('n')}; median T/NK high={fmt(w.get('median_a'))}, low={fmt(w.get('median_b'))}; Wilcoxon p={fmt(w.get('p'), sci=True)}.
"""

    rows = []
    for key, s in graphs.items():
        da = s["da_malignant_cldn4"]
        rows.append(
            f"| {key} | {s['honest_n']['independent_unit_DA']} | **{s['units']['n_with_malignant_ge10']}** | "
            f"{da['n_testable']} | **{da['n_SpatialFDR_lt_0.1']}** | {fmt(da.get('min_SpatialFDR'))} |"
        )
    table = "\n".join(rows)
    harmony = overall.get("harmony") or {}
    if isinstance(harmony, dict) and harmony.get("ran"):
        h_note = (
            f"Harmony joint graph ran: shared HVG={harmony.get('shared_hvg')}, "
            f"cap={harmony.get('max_per_unit')} cells/unit, "
            f"cells after cap={harmony.get('n_cells_after_cap')}."
        )
    else:
        reason = harmony.get("reason") if isinstance(harmony, dict) else "not attempted"
        h_note = f"Harmony joint graph was **not** used as the primary result ({reason}). Per-dataset graphs are the DA."

    blocks = []
    order = [
        "GSE131907_tLung",
        "GSE131907_tumor_no_brain",
        "GSE205335",
        "harmony_tumor_no_brain_plus_GSE205335",
    ]
    for key in order:
        if key in graphs:
            blocks.append(block(key, graphs[key]))
    for key, s in graphs.items():
        if key not in order:
            blocks.append(block(key, s))

    text = f"""# FINDING — merged GSE131907 + GSE205335 Milo vs malignant CLDN4

Additive **CLDN4-only**. Public UMIs from Kim et al. 2020 (GSE131907; PMID
32385277) and Ahn / Lee (GSE205335; eLife 98366). Neighbourhood abundance
is tested against **sample/patient malignant CLDN4**.

**PR #320 T/NK ρ is given and is not re-audited.** Author %pos vs T/NK on
GSE131907+GSE205335: k=2, N=43, ρ=−0.479, p=0.00152, I²=0%. This folder
does not recompute or re-rank that combinatorial pool.

No dual-high TACSTD2×CLDN4 score. **GSE207422 is not run** (that Milo was
SpatialFDR-null because only n=7 samples had ≥10 malignant-like cells).

The independent unit is the **sample** (GSE131907) or the **patient**
(GSE205335). Do not cite cell count as *n*. Per-dataset SpatialFDR tables
are not pooled into one *n*.

**Not miloR.** DA = sample/patient Spearman (primary) or Welch t-test
(median split) on neighbourhood proportions. SpatialFDR = miloR
`graphSpatialFDR` k-distance weights (Dann et al. 2022 / cydar).

{h_note}

Public matrices: GSE131907 {overall['gse131907']['n_cells_matrix']:,} cells,
{overall['gse131907']['n_genes_matrix']:,} genes ({overall['gse131907']['matrix_bytes']:,} bytes gzip);
GSE205335 {overall['gse205335']['n_cells_matrix']:,} cells,
{overall['gse205335']['n_genes_matrix']:,} genes ({overall['gse205335']['matrix_bytes']:,} bytes gzip).
The 3 GB GSE131907 log2TPM text and EGA FASTQ were not used.

## One-row SpatialFDR

| Graph | Unit | Honest n (scored) | testable nhoods | SpatialFDR<0.1 | min SpatialFDR |
|---|---|---:|---:|---:|---:|
{table}

## n and SpatialFDR by graph

{"".join(blocks)}

## What this does not say

- It does not re-audit PR #320 T/NK ρ.
- It does not invent an ICI / MPR contrast on GSE131907 (treatment-naive).
- GSE205335 MPR is not labelled; RECIST is not substituted as the primary DA.
- It does not treat 200k cells as *n*.
- Neighbourhood “next to” is kNN co-membership, not histology.
- miloR / edgeR QLF numbers are not claimed.
- GSE207422 was not merged in (n=7 SpatialFDR-null).
- Unrestricted neighbourhood CLDN4 vs T/NK is partly lineage geometry.

See `tables/nhoods.tsv`, `tables/honest_n.tsv`, `tables/one_row.tsv`,
and per-graph folders under `tables/`.
"""
    out_path.write_text(text)


def extra_summary_figures(graphs: dict, figdir: Path, tables: Path) -> None:
    names = list(graphs)
    n_scored = [graphs[n]["units"]["n_with_malignant_ge10"] for n in names]
    n_hit = [graphs[n]["da_malignant_cldn4"]["n_SpatialFDR_lt_0.1"] for n in names]
    n_test = [graphs[n]["da_malignant_cldn4"]["n_testable"] for n in names]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    y = np.arange(len(names))
    ax.barh(y - 0.18, n_test, height=0.35, color="#4c72b0", label="testable nhoods")
    ax.barh(y + 0.18, n_hit, height=0.35, color="#c44e52", label="SpatialFDR<0.1")
    ax.set_yticks(y, names)
    ax.set_xlabel("Neighbourhoods")
    ax.set_title("Merged CLDN4 Milo: testable vs SpatialFDR<0.1 (honest n in FINDING)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figdir / "fig_spatialfdr_counts.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.barh(np.arange(len(names)), n_scored, color="#55a868")
    ax.set_yticks(np.arange(len(names)), names)
    ax.set_xlabel("Scored samples/patients (≥10 malignant cells)")
    ax.set_title("Honest n (unit = sample or patient)")
    fig.tight_layout()
    fig.savefig(figdir / "fig_honest_n.png", dpi=160)
    plt.close(fig)

    rows = []
    for n, s in graphs.items():
        da = s["da_malignant_cldn4"]
        rows.append(
            {
                "graph": n,
                "n_units": s["units"]["n"],
                "n_scored": s["units"]["n_with_malignant_ge10"],
                "n_cells_not_n": s["n_cells_in_graph"],
                "n_nhoods": da["n_nhoods_total"],
                "n_testable": da["n_testable"],
                "n_p_lt_0.05": da["n_p_lt_0.05"],
                "n_SpatialFDR_lt_0.1": da["n_SpatialFDR_lt_0.1"],
                "n_SpatialFDR_lt_0.05": da["n_SpatialFDR_lt_0.05"],
                "min_SpatialFDR": da["min_SpatialFDR"],
                "min_p": da["min_p"],
                "n_BH_FDR_lt_0.1": da["n_BH_FDR_lt_0.1"],
            }
        )
    pd.DataFrame(rows).to_csv(tables / "honest_n.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(tables / "one_row.tsv", sep="\t", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir-131907", type=Path, default=Path("/tmp/GSE131907"))
    ap.add_argument("--datadir-205335", type=Path, default=Path("/tmp/GSE205335"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/merge_131907_205335_milo_cldn4"))
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--d", type=int, default=30)
    ap.add_argument("--prop", type=float, default=0.1)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-interface", type=int, default=3)
    ap.add_argument("--min-malignant", type=int, default=10)
    ap.add_argument("--harmony-max-per-unit", type=int, default=1200)
    ap.add_argument("--skip-harmony", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    tables = args.outdir / "tables"
    figures = args.outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    self_test()

    print("=== GSE205335 ===", flush=True)
    d205 = load_gse205335(args.datadir_205335, args.n_hvg)
    r205 = run_graph(
        name="GSE205335",
        X=d205["X"],
        sample=d205["sample"],
        lineage=d205["lineage"],
        is_malig=d205["is_malig"],
        is_tnk=d205["is_tnk"],
        cldn=d205["cldn4"],
        extra_meta=None,
        outdir=tables,
        figdir=figures,
        k=args.k,
        d=args.d,
        prop=args.prop,
        n_hvg=args.n_hvg,
        min_interface=args.min_interface,
        seed=args.seed,
        min_malig=args.min_malignant,
        batch_note="PCA per-patient mean centering (not Harmony)",
    )
    pack205 = {
        "dataset": "GSE205335",
        "hvg_used": d205["hvg"],
        "gene_log": d205["gene_log"],
        "sample": d205["sample"],
        "lineage": d205["lineage"],
        "is_malig": d205["is_malig"],
        "is_tnk": d205["is_tnk"],
        "cldn4": d205["cldn4"],
        "local": None,
    }
    d205_meta = {
        "n_cells_matrix": d205["n_cells_matrix"],
        "n_genes_matrix": d205["n_genes_matrix"],
        "matrix_bytes": d205["matrix_bytes"],
    }
    del d205
    gc.collect()

    print("=== GSE131907 ===", flush=True)
    d131 = load_gse131907(args.datadir_131907, args.n_hvg)
    graphs = {"GSE205335": r205["summary"]}
    nhood_frames = [r205["nhoods"]]
    score_frames = [r205["scores"]]

    for graph_name, local in (
        ("GSE131907_tLung", d131["local_tLung"]),
        ("GSE131907_tumor_no_brain", np.ones(d131["sample"].size, dtype=bool)),
    ):
        key = "tLung" if graph_name.endswith("tLung") else "tumor_no_brain"
        X = matrix_from_hvg(d131["expr_log"], d131["hvg"][key], local)
        r = run_graph(
            name=graph_name,
            X=X,
            sample=d131["sample"][local],
            lineage=d131["lineage"][local],
            is_malig=d131["is_malig"][local],
            is_tnk=d131["is_tnk"][local],
            cldn=d131["cldn4"][local],
            extra_meta=None,
            outdir=tables,
            figdir=figures,
            k=args.k,
            d=args.d,
            prop=args.prop,
            n_hvg=args.n_hvg,
            min_interface=args.min_interface,
            seed=args.seed,
            min_malig=args.min_malignant,
            batch_note="PCA per-sample mean centering (not Harmony)"
            + ("" if key == "tLung" else "; mixed sites tLung+tL/B+mLN, not one tissue"),
        )
        graphs[graph_name] = r["summary"]
        nhood_frames.append(r["nhoods"])
        score_frames.append(r["scores"])
        del X, r
        gc.collect()

    harmony_info: dict = {"ran": False, "reason": "skipped by --skip-harmony"}
    if not args.skip_harmony:
        pack131 = {
            "dataset": "GSE131907",
            "hvg_used": d131["hvg"]["tumor_no_brain"],
            "gene_log": d131["expr_log"],
            "sample": d131["sample"],
            "lineage": d131["lineage"],
            "is_malig": d131["is_malig"],
            "is_tnk": d131["is_tnk"],
            "cldn4": d131["cldn4"],
            "local": None,
        }
        print("=== Harmony extra ===", flush=True)
        h = try_harmony([pack131, pack205], tables, figures, args)
        if isinstance(h, dict) and h.get("summary"):
            graphs[h["summary"]["graph"]] = h["summary"]
            nhood_frames.append(h["nhoods"])
            score_frames.append(h["scores"])
            harmony_info = h["summary"].get("harmony", {"ran": True})
        elif isinstance(h, dict):
            harmony_info = h

    nhoods_all = pd.concat(nhood_frames, ignore_index=True)
    nhoods_all.to_csv(tables / "nhoods.tsv", sep="\t", index=False)
    pd.concat(score_frames, ignore_index=True).to_csv(tables / "sample_scores.tsv", sep="\t", index=False)
    extra_summary_figures(graphs, figures, tables)

    overall = {
        "datasets": ["GSE131907", "GSE205335"],
        "gene": "CLDN4",
        "dual_high": False,
        "gse207422_run": False,
        "gse207422_reason": "SpatialFDR-null at n=7; out of scope",
        "pr320_tnk_rho_given": GIVEN_PR320,
        "miloR": False,
        "da_model": "sample/patient Spearman or Welch on nhood proportions; not edgeR QLF",
        "spatial_fdr": "miloR graphSpatialFDR k-distance reimplementation (Dann 2022 / cydar)",
        "gse131907": {
            "n_cells_matrix": d131["n_cells_matrix"],
            "n_genes_matrix": d131["n_genes_matrix"],
            "matrix_bytes": d131["matrix_bytes"],
            "citation": "Kim et al. Nat Commun 2020 PMID 32385277",
        },
        "gse205335": {
            "n_cells_matrix": int(d205_meta["n_cells_matrix"]),
            "n_genes_matrix": int(d205_meta["n_genes_matrix"]),
            "matrix_bytes": int(d205_meta["matrix_bytes"]),
            "citation": "Ahn / Lee GSE205335 eLife 98366",
        },
        "harmony": harmony_info,
        "graphs": graphs,
        "skipped": [
            "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz",
            "EGAD00001005054",
            "EGAD00001008703",
            "GSE207422",
        ],
    }
    write_json(tables / "summary.json", overall)
    write_finding(overall, args.outdir / "FINDING.md")
    print(
        json.dumps(
            {
                "finding": str(args.outdir / "FINDING.md"),
                "nhoods": str(tables / "nhoods.tsv"),
                "graphs": {k: v["da_malignant_cldn4"] for k, v in graphs.items()},
                "harmony": harmony_info,
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
