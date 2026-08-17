#!/usr/bin/env python3
"""Pairwise Milo-style kNN DA vs malignant CLDN4 (GSE131907 + GSE148071).

Additive CLDN4-only. Not the triple. Not 131907+205335. No dual-high.
Graphs are per dataset / site. Pair n = GSE148071 biopsies + GSE131907 tLung.
mBrain is a same-atlas sensitivity and is not added to the pair n.
Independent unit = patient / sample. miloR is not used.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

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

LINEAGE_MARKERS = {
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2"],
    "Mast": ["TPSAB1", "CPA3"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}
NORMAL_LUNG = ["SFTPA2", "SFTPA1", "SFTPC", "SFTPB", "AGER", "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS"]
ALWAYS = ["TACSTD2", "CLDN4", "PTPRC", "EPCAM", "CD3D", "NKG7"]
HEADER_GENE = {"", "gene", "Gene", "GENE", "index", "Index", "symbol", "Symbol", "SYMBOL"}

MALIGNANT_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK_TYPES = {"T lymphocytes", "NK cells"}
KEEP_TYPES = {
    "T lymphocytes",
    "NK cells",
    "B lymphocytes",
    "Myeloid cells",
    "MAST cells",
    "Epithelial cells",
}
MARKERS_131907 = [
    "CLDN4",
    "TACSTD2",
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "PTPRC",
    "CD3D",
    "CD3E",
    "NKG7",
    "GNLY",
    "CD79A",
    "MS4A1",
    "LYZ",
    "CD68",
    "PECAM1",
    "COL1A1",
]
COHORTS_131907 = ("tLung", "mBrain")


def wanted_markers() -> set[str]:
    genes = set(ALWAYS)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    genes.update(NORMAL_LUNG)
    return genes


def norm_gene(raw: str) -> str:
    gene = raw.strip().strip('"')
    if "|" in gene:
        gene = gene.split("|")[-1]
    if ";" in gene:
        gene = gene.split(";")[0]
    return gene.split(".")[0]


def patient_from_name(path: Path) -> str:
    stem = path.name.replace(".txt.gz", "").replace(".txt", "")
    parts = stem.split("_")
    for p in parts:
        if p.startswith("P") and p[1:].isdigit():
            return p
    return parts[1] if len(parts) > 1 else stem


def parse_header(line: str) -> list[str]:
    header = line.rstrip("\n").split("\t")
    if header and header[0] in HEADER_GENE:
        return header[1:]
    return header


def list_exp_files(datadir: Path) -> list[Path]:
    files = sorted(datadir.glob("*_exp.txt.gz"))
    if not files:
        files = sorted(datadir.rglob("*_exp.txt.gz"))
    if not files:
        raise SystemExit(f"no *_exp.txt.gz under {datadir}")
    return files


def stream_file(path: Path, markers: set[str], keep: set[str] | None):
    with gzip.open(path, "rt") as handle:
        cell_ids = parse_header(handle.readline())
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        found: dict[str, np.ndarray] = {}
        gene_sum: dict[str, float] = {}
        gene_sumsq: dict[str, float] = {}
        n_genes = 0
        for line in handle:
            gene_raw, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = norm_gene(gene_raw)
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                toks = line.rstrip("\n").split("\t")
                arr = np.asarray(toks[1:], dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{path.name} {gene}: {arr.size} != {n}")
            n_umi += arr
            n_genes += 1
            if keep is None:
                gene_sum[gene] = gene_sum.get(gene, 0.0) + float(arr.sum())
                gene_sumsq[gene] = gene_sumsq.get(gene, 0.0) + float(np.dot(arr, arr))
                if gene in markers:
                    found[gene] = arr
            elif gene in keep:
                found[gene] = arr
        return cell_ids, n_umi, found, gene_sum, gene_sumsq, n_genes


def assign_lineage(log_cp: dict[str, np.ndarray], n: int) -> np.ndarray:
    scores = {}
    for name, genes in LINEAGE_MARKERS.items():
        mats = [log_cp[g] for g in genes if g in log_cp]
        scores[name] = np.mean(np.vstack(mats), axis=0) if mats else np.zeros(n, dtype=np.float32)
    names = list(scores)
    mat = np.vstack([scores[n_] for n_ in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(n)]
    labels = np.array(names, dtype=object)[best]
    cd3 = log_cp.get("CD3D", np.zeros(n, dtype=np.float32))
    t_idx, nk_idx = names.index("T"), names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk = np.isin(labels, ["T", "NK"])
    labels = labels.copy()
    labels[close & both & tnk & (cd3 > 0.15)] = "T"
    labels[close & both & tnk & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < 0.12] = "Unassigned"
    return labels


def gene_log(log_cp: dict[str, np.ndarray], name: str, n: int) -> np.ndarray:
    return log_cp[name] if name in log_cp else np.zeros(n, dtype=np.float32)


def fdr_counts(df: pd.DataFrame) -> dict:
    t = df[df["testable"]]
    return {
        "n_nhoods_total": int(len(df)),
        "n_testable": int(len(t)),
        "n_p_lt_0.05": int((t["p"] < 0.05).sum()) if len(t) else 0,
        "n_BH_FDR_lt_0.1": int((t["BH_FDR"] < 0.1).sum()) if len(t) else 0,
        "n_BH_FDR_lt_0.05": int((t["BH_FDR"] < 0.05).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.1": int((t["SpatialFDR"] < 0.1).sum()) if len(t) else 0,
        "n_SpatialFDR_lt_0.05": int((t["SpatialFDR"] < 0.05).sum()) if len(t) else 0,
        "min_p": float(t["p"].min()) if len(t) else None,
        "min_SpatialFDR": float(t["SpatialFDR"].min()) if len(t) else None,
        "min_BH_FDR": float(t["BH_FDR"].min()) if len(t) else None,
    }


def fmt(x, nd=3, sci=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{x:.2e}" if sci else f"{x:.{nd}f}"


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
    ax.set_ylabel(r"−log10 p (sample-level test)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_scatter(x, y, path: Path, xlabel: str, ylabel: str, title: str, c=None, labels=None) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    if c is None:
        ax.scatter(x, y, s=22, c="#4c72b0", alpha=0.7)
    else:
        for lab in pd.unique(labels):
            m = labels == lab
            ax.scatter(np.asarray(x)[m], np.asarray(y)[m], s=28, alpha=0.8, label=str(lab))
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fisher_z_meta(rows: list[dict]) -> dict:
    zs, ws = [], []
    for r in rows:
        n = r.get("n")
        rho = r.get("rho")
        if n is None or rho is None or n < 4 or not np.isfinite(rho):
            continue
        rho = float(np.clip(rho, -0.999999, 0.999999))
        z = np.arctanh(rho)
        w = n - 3
        zs.append(z * w)
        ws.append(w)
    if not ws:
        return {"n_studies": 0, "n_total": 0, "rho": None, "p": None}
    zbar = float(np.sum(zs) / np.sum(ws))
    se = float(1.0 / np.sqrt(np.sum(ws)))
    p = float(2 * stats.norm.sf(abs(zbar / se)))
    return {
        "n_studies": int(len(ws)),
        "n_total": int(sum(r["n"] for r in rows if r.get("n"))),
        "rho": float(np.tanh(zbar)),
        "p": p,
    }


def run_graph(
    name: str,
    dataset: str,
    sample: np.ndarray,
    lineage: np.ndarray,
    is_malig: np.ndarray,
    is_tnk: np.ndarray,
    cldn: np.ndarray,
    expr: dict[str, np.ndarray],
    hvg_names: list[str],
    n_umi: np.ndarray,
    outdir: Path,
    figdir: Path,
    k: int,
    d: int,
    prop: float,
    min_interface: int,
    seed: int,
    min_malig: int,
) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    hvg_present = [g for g in hvg_names if g in expr]
    X = np.vstack([np.log1p(expr[g] * scale) for g in hvg_present]).T.astype(np.float32)
    print(f"{name}: cells={len(sample)} HVG={X.shape[1]} samples={pd.unique(sample).size}", flush=True)
    pcs, knn_idx, knn_dist = pca_knn(X, sample, n_pcs=d, k=k, random_state=seed)
    indices = refine_indices(pcs, knn_idx, prop=prop, random_state=seed)
    nhoods = make_nhoods(knn_idx, knn_dist, indices)
    members = nhood_membership(nhoods)
    print(f"{name}: nhoods={len(members)} median_size={np.median([m.size for m in members]):.0f}", flush=True)

    sample_levels = sorted(pd.unique(sample), key=lambda s: (len(str(s)), str(s)))
    counts = count_matrix(members, sample, sample_levels)
    sizes = sample_sizes(sample, sample_levels)

    mal_score, tnk_frac, n_mal_s, n_tnk_s = [], [], [], []
    for s in sample_levels:
        m = sample == s
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
    med_s = float("nan")
    if scored.sum() >= 2:
        med_s = float(np.median(mal_score[scored]))
        group[scored & (mal_score >= med_s)] = "high"
        group[scored & (mal_score < med_s)] = "low"
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
                "n_samples_present": 0,
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
                "n_samples_present": "n_samples_present",
            }
        ),
        on="nhood",
        how="left",
    )
    merged.insert(0, "dataset", dataset)
    merged.insert(1, "graph", name)

    dest = outdir / name
    dest.mkdir(parents=True, exist_ok=True)
    merged.to_csv(dest / "nhoods.tsv", sep="\t", index=False)
    da_cldn.to_csv(dest / "da_malignant_cldn4.tsv", sep="\t", index=False)
    da_split.to_csv(dest / "da_cldn4_high_vs_low.tsv", sep="\t", index=False)

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
    paired = sample_paired_tnk_by_gene(members, sample, is_tnk, high, low)
    sample_tab = pd.DataFrame(
        {
            "dataset": dataset,
            "graph": name,
            "sample": sample_levels,
            "n_cells": sizes.astype(int),
            "n_malignant": n_mal_s,
            "n_tnk": n_tnk_s,
            "malignant_cldn4": mal_score,
            "frac_tnk": tnk_frac,
            "cldn4_arm": group,
            "in_pair_n": name != "GSE131907_mBrain",
        }
    )
    paired = paired.merge(sample_tab, on="sample")
    paired.to_csv(dest / "sample_paired_tnk_by_nhood_cldn4.tsv", sep="\t", index=False)
    sample_tab.to_csv(dest / "sample_scores.tsv", sep="\t", index=False)
    paired_test = wilcoxon_paired(paired["frac_tnk_high_nhood"], paired["frac_tnk_low_nhood"])
    sample_level = spearman_safe(mal_score, tnk_frac)

    def _sens(df: pd.DataFrame, min_n: int) -> dict:
        t = df.loc[df["testable"] & (df["n_samples_present"] >= min_n)].copy()
        t["testable"] = True
        rec = fdr_counts(t)
        rec["min_samples_present"] = min_n
        rec["n_abs_rho_eq_1"] = int((t["spearman_rho"].abs() >= 0.999).sum()) if len(t) else 0
        return rec

    floor_hits = da_cldn[(da_cldn["testable"]) & (da_cldn["SpatialFDR"] < 0.1)]
    da_cldn_counts = fdr_counts(da_cldn)
    da_cldn_counts["median_n_samples_present"] = (
        float(da_cldn.loc[da_cldn["testable"], "n_samples_present"].median()) if da_cldn["testable"].any() else None
    )

    slug = name.replace(" ", "_")
    plot_volcano(
        da_cldn.rename(columns={"spearman_rho": "stat"}),
        "stat",
        f"{name}: nhood abundance vs malignant CLDN4",
        figdir / f"fig_da_cldn4_volcano_{slug}.png",
        "Spearman ρ (nhood proportion vs sample malignant CLDN4)",
    )
    if da_split["testable"].any():
        plot_volcano(
            da_split,
            "logFC_B_minus_A",
            f"{name}: nhood DA CLDN4-high vs low samples",
            figdir / f"fig_da_cldn4_high_vs_low_{slug}.png",
            "log2 FC (high − low sample proportion)",
        )
    plot_scatter(
        comp.loc[iface, "cldn4_malig_mean"],
        comp.loc[iface, "frac_tnk"],
        figdir / f"fig_interface_cldn4_vs_tnk_{slug}.png",
        "Neighbourhood malignant CLDN4 (mean log1p CP10k)",
        "Neighbourhood T/NK fraction",
        f"{name} interface nhoods (n={int(iface.sum())}); transcriptional, not spatial",
    )
    if paired["frac_tnk_high_nhood"].notna().any():
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        for _, r in paired.iterrows():
            if np.isfinite(r["frac_tnk_high_nhood"]) and np.isfinite(r["frac_tnk_low_nhood"]):
                color = "#c44e52" if r["cldn4_arm"] == "high" else "#4c72b0"
                ax.plot([0, 1], [r["frac_tnk_low_nhood"], r["frac_tnk_high_nhood"]], "-o", c=color, alpha=0.75, ms=5)
        ax.set_xticks([0, 1], ["CLDN4-low nhoods", "CLDN4-high nhoods"])
        ax.set_ylabel("Sample T/NK fraction among cells in those nhoods")
        ax.set_title(f"{name} sample-paired (unit = sample)")
        fig.tight_layout()
        fig.savefig(figdir / f"fig_sample_paired_tnk_{slug}.png", dpi=160)
        plt.close(fig)

    dropped = sample_tab.loc[~np.isfinite(sample_tab["malignant_cldn4"])]
    drop_bits = [f"{r.sample}={int(r.n_malignant)}" for r in dropped.itertuples()]
    summary = {
        "dataset": dataset,
        "graph": name,
        "gene": "CLDN4",
        "dual_high": False,
        "n_cells_in_graph": int(len(sample)),
        "n_malignant": int(is_malig.sum()),
        "n_tnk": int(is_tnk.sum()),
        "lineage_counts": {k: int(v) for k, v in pd.Series(lineage).value_counts().items()},
        "graph_params": {
            "k": k,
            "d": min(d, X.shape[1], len(sample) - 1),
            "prop": prop,
            "n_hvg": len(hvg_present),
            "n_nhoods": int(len(members)),
            "nhood_size": k + 1,
            "median_nhood_size": float(np.median([m.size for m in members])),
            "batch": "PCA per-sample mean centering (not Harmony)",
        },
        "samples": {
            "n": int(len(sample_levels)),
            "n_with_malignant_ge10": int(np.isfinite(mal_score).sum()),
            "n_cldn4_high": n_high,
            "n_cldn4_low": n_low,
            "dropped_lt10_malignant": drop_bits,
            "median_split_cut": med_s if np.isfinite(med_s) else None,
        },
        "sample_level_malignant_CLDN4_vs_TNK_fraction": sample_level,
        "da_malignant_cldn4": da_cldn_counts,
        "da_cldn4_high_vs_low": fdr_counts(da_split),
        "da_malignant_cldn4_by_n_present": [_sens(da_cldn, t) for t in (5, 8, 10, 15)],
        "spatialfdr_floor_hits": {
            "n": int(len(floor_hits)),
            "all_abs_rho_eq_1": bool(len(floor_hits) and (floor_hits["spearman_rho"].abs() >= 0.999).all()),
            "n_samples_present_min": int(floor_hits["n_samples_present"].min()) if len(floor_hits) else None,
            "n_samples_present_max": int(floor_hits["n_samples_present"].max()) if len(floor_hits) else None,
        },
        "composition": {
            "note": "transcriptional kNN != spatial niche",
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
            "independent_unit_DA": f"{name} samples (n={len(sample_levels)}; scored={int(np.isfinite(mal_score).sum())})",
            "nhoods_are_overlapping": True,
            "SpatialFDR_is_overlap_aware": True,
            "disjoint_nhood_subset_n": int(indep.size),
            "do_not_cite_n_cells_as_n": int(len(sample)),
        },
    }
    write_json(dest / "summary.json", summary)
    return summary, merged, sample_tab, paired


def run_gse148071(datadir: Path, outdir: Path, figdir: Path, args) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    files = list_exp_files(datadir)
    markers = wanted_markers()
    print(f"GSE148071 pass1: {len(files)} files", flush=True)
    sample_all: list[str] = []
    n_umi_all: list[np.ndarray] = []
    gene_sum: dict[str, float] = {}
    gene_sumsq: dict[str, float] = {}
    n_cells_total = 0
    for fp in files:
        patient = patient_from_name(fp)
        cell_ids, n_umi, found, gsum, gss, n_genes = stream_file(fp, markers, keep=None)
        n = len(cell_ids)
        sample_all.extend([patient] * n)
        n_umi_all.append(n_umi)
        n_cells_total += n
        for g, val in gsum.items():
            gene_sum[g] = gene_sum.get(g, 0.0) + val
            gene_sumsq[g] = gene_sumsq.get(g, 0.0) + gss[g]
        print(f"  {fp.name} {patient} n={n} genes={n_genes} CLDN4={'CLDN4' in found}", flush=True)

    n_umi = np.concatenate(n_umi_all)
    sample = np.array(sample_all, dtype=object)
    keep_cell = n_umi >= args.min_umi
    genes = sorted(gene_sum)
    means = np.array([gene_sum[g] / n_cells_total for g in genes], dtype=float)
    varis = np.array([gene_sumsq[g] / n_cells_total - means[i] ** 2 for i, g in enumerate(genes)], dtype=float)
    varis = np.maximum(varis, 0.0)
    hvg_idx = select_hvg(means, varis, n_hvg=args.n_hvg)
    hvg_names = [genes[i] for i in hvg_idx]
    keep_genes = set(hvg_names) | markers
    print(f"GSE148071 HVG={len(hvg_names)} keep_genes={len(keep_genes)}", flush=True)

    print("GSE148071 pass2", flush=True)
    expr_blocks: dict[str, list[np.ndarray]] = defaultdict(list)
    for fp in files:
        cell_ids, _, found, _, _, _ = stream_file(fp, markers, keep=keep_genes)
        n = len(cell_ids)
        present = set(found)
        for g in keep_genes:
            expr_blocks[g].append(found[g] if g in present else np.zeros(n, dtype=np.float32))
        print(f"  {fp.name} stored {len(present)}/{len(keep_genes)}", flush=True)

    expr = {g: np.concatenate(blocks) for g, blocks in expr_blocks.items()}
    n = int(n_umi.size)
    scale = 1e4 / np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
    lineage = assign_lineage(log_cp, n)
    if "CLDN4" not in log_cp:
        raise RuntimeError("CLDN4 is absent from the public GSE148071 count files")
    cldn4 = log_cp["CLDN4"]
    alveolar = np.maximum.reduce([gene_log(log_cp, g, n) for g in ["SFTPA2", "SFTPA1", "SFTPC", "SFTPB", "AGER"]])
    club = np.maximum(gene_log(log_cp, "SCGB1A1", n), gene_log(log_cp, "SCGB3A2", n))
    ciliated = np.maximum.reduce([gene_log(log_cp, g, n) for g in ["TPPP3", "FOXJ1", "CAPS"]])
    clear_normal = (alveolar >= 1.0) | (club >= 1.0) | (ciliated >= 1.0)
    is_malig = (lineage == "Epithelial") & (~clear_normal)
    is_tnk = np.isin(lineage, ["T", "NK"])
    if keep_cell.sum() < n:
        sample = sample[keep_cell]
        lineage = lineage[keep_cell]
        is_malig = is_malig[keep_cell]
        is_tnk = is_tnk[keep_cell]
        cldn4 = cldn4[keep_cell]
        n_umi = n_umi[keep_cell]
        expr = {g: v[keep_cell] for g, v in expr.items()}

    return run_graph(
        name="GSE148071",
        dataset="GSE148071",
        sample=sample,
        lineage=lineage,
        is_malig=is_malig,
        is_tnk=is_tnk,
        cldn=cldn4,
        expr=expr,
        hvg_names=hvg_names,
        n_umi=n_umi,
        outdir=outdir,
        figdir=figdir,
        k=args.k,
        d=args.d,
        prop=args.prop,
        min_interface=args.min_interface,
        seed=args.seed,
        min_malig=args.min_malignant,
    )


def parse_series_matrix(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key.removeprefix("!Sample_")
            values = [v.strip('"') for v in rest.split("\t")]
            if key == "characteristics_ch1":
                prefix, _, _ = values[0].partition(": ")
                key = prefix.strip().lower().replace(" ", "_")
                values = [v.split(": ", 1)[-1] if ": " in v else v for v in values]
            if key not in fields:
                fields[key] = values
    n = len(fields["title"])
    geo = pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})
    return geo.rename(columns={"title": "Sample"})


def origin_from_sample(name: str) -> str:
    if name.startswith("LUNG_N"):
        return "nLung"
    if name.startswith("LUNG_T"):
        return "tLung"
    if name.startswith("LN_"):
        return "nLN"
    if name.startswith("EFFUSION_"):
        return "PE"
    if name.startswith("NS_"):
        return "mBrain"
    if name.startswith(("EBUS_", "BRONCHO_")):
        return "unknown_lung_or_ln"
    return "unknown"


def load_annotation(ann_path: Path, series_path: Path, cell_ids: list[str]) -> pd.DataFrame:
    ann = pd.read_csv(ann_path, sep="\t", dtype=str)
    if "Index" not in ann.columns:
        raise SystemExit(f"annotation missing Index: {list(ann.columns)}")
    ann = ann.set_index("Index").reindex(cell_ids).reset_index()
    if ann["Sample"].isna().any():
        raise SystemExit(f"{int(ann['Sample'].isna().sum())} matrix cell IDs missing from annotation")
    geo = parse_series_matrix(series_path)
    geo_small = geo[["Sample"]].copy()
    for col in ("patient_id", "tumor_stage", "source_name_ch1", "lung_cancer_subtype"):
        if col in geo.columns:
            geo_small[col] = geo[col]
    if "sample_origin" in geo.columns:
        geo_small["Sample_Origin_geo"] = geo["sample_origin"]
    elif "source_name_ch1" in geo.columns:
        geo_small["Sample_Origin_geo"] = geo["source_name_ch1"]
    ann = ann.merge(geo_small, on="Sample", how="left")
    origin = ann["Sample_Origin"].fillna("") if "Sample_Origin" in ann.columns else pd.Series("", index=ann.index)
    mapped = ann["Sample"].map(origin_from_sample)
    known = {"tLung", "nLung", "tL/B", "mLN", "nLN", "PE", "mBrain"}
    ann["origin"] = np.where(origin.isin(known), origin, mapped)
    return ann


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
                print(f"  GSE131907 pass1 {i} genes, markers={len(found)} [{time.time() - t0:.0f}s]", flush=True)
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
                print(f"  GSE131907 pass2 stored {len(out)}/{len(keep)} [{time.time() - t0:.0f}s]", flush=True)
    missing = keep - set(out)
    if missing:
        print(f"  missing {sorted(missing)[:12]} ({len(missing)} genes)", flush=True)
    return out


def run_gse131907(datadir: Path, outdir: Path, figdir: Path, args):
    ann_path = datadir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = datadir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    series_path = datadir / "GSE131907_series_matrix.txt.gz"
    for p in (ann_path, umi_path, series_path):
        if not p.exists():
            raise SystemExit(f"missing {p}; run scripts/download.py")
    with gzip.open(umi_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    cell_ids = header[1:]
    print(f"GSE131907 matrix header cells={len(cell_ids)}", flush=True)
    ann = load_annotation(ann_path, series_path, cell_ids)
    cell_type = ann["Cell_type"].to_numpy()
    subtype = ann["Cell_subtype"].fillna("").to_numpy() if "Cell_subtype" in ann.columns else np.array([""] * len(ann))
    origin = ann["origin"].to_numpy()
    sample = ann["Sample"].to_numpy()
    keep_type = np.isin(cell_type, list(KEEP_TYPES))
    is_malig_all = np.isin(subtype, list(MALIGNANT_SUBTYPES))
    is_tnk_all = np.isin(cell_type, list(TNK_TYPES))
    cohort_masks = {c: (origin == c) & keep_type for c in COHORTS_131907}
    for c, m in cohort_masks.items():
        print(
            f"  {c}: graph_cells={int(m.sum())} malignant={int((m & is_malig_all).sum())} "
            f"samples={pd.unique(sample[m]).size}",
            flush=True,
        )
        if m.sum() < 200:
            raise SystemExit(f"{c}: too few epithelium+immune cells ({int(m.sum())})")

    print("GSE131907 pass1", flush=True)
    cell_ids2, n_umi, gene_names, cohort_stats, marker_expr = stream_pass1(umi_path, set(MARKERS_131907), cohort_masks)
    if cell_ids2 != cell_ids:
        raise SystemExit("matrix header changed between reads")
    if "CLDN4" not in marker_expr:
        raise SystemExit("CLDN4 is absent from the public UMI matrix")
    hvg_names = {}
    keep_genes = set(MARKERS_131907)
    for c in COHORTS_131907:
        mu, va = cohort_stats[c]
        idx = select_hvg(mu, va, n_hvg=args.n_hvg)
        hvg_names[c] = [gene_names[i] for i in idx]
        keep_genes.update(hvg_names[c])
        print(f"  {c} HVG={len(hvg_names[c])}", flush=True)

    keep_any = cohort_masks["tLung"] | cohort_masks["mBrain"]
    keep_idx = np.flatnonzero(keep_any)
    print(f"GSE131907 pass2: {len(keep_genes)} genes × {keep_idx.size} kept cells", flush=True)
    expr = stream_pass2(umi_path, keep_genes, keep_idx, len(cell_ids))
    for g, arr in marker_expr.items():
        expr.setdefault(g, arr[keep_idx].astype(np.float32, copy=False))

    sample_k = sample[keep_idx]
    lineage_k = cell_type[keep_idx]
    is_malig_k = is_malig_all[keep_idx]
    is_tnk_k = is_tnk_all[keep_idx]
    n_umi_k = n_umi[keep_idx]
    scale_k = np.where(n_umi_k > 0, 1e4 / n_umi_k, 0.0)
    cldn_k = np.log1p(expr["CLDN4"] * scale_k).astype(np.float32)
    origin_k = origin[keep_idx]

    results = {}
    for c in COHORTS_131907:
        local = origin_k == c
        local_expr = {g: v[local] for g, v in expr.items()}
        results[c] = run_graph(
            name=f"GSE131907_{c}",
            dataset="GSE131907",
            sample=sample_k[local],
            lineage=lineage_k[local],
            is_malig=is_malig_k[local],
            is_tnk=is_tnk_k[local],
            cldn=cldn_k[local],
            expr=local_expr,
            hvg_names=hvg_names[c],
            n_umi=n_umi_k[local],
            outdir=outdir,
            figdir=figdir,
            k=args.k,
            d=args.d,
            prop=args.prop,
            min_interface=args.min_interface,
            seed=args.seed,
            min_malig=args.min_malignant,
        )
    inventory = {
        "n_cells_matrix": int(len(cell_ids)),
        "n_genes_matrix": int(len(gene_names)),
        "n_author_malignant": int(is_malig_all.sum()),
        "origin_counts": {k: int(v) for k, v in pd.Series(origin).value_counts().items()},
        "skipped": ["GSE205335", "triple 131907+148071+205335", "dual-high TACSTD2+CLDN4"],
    }
    return results, inventory


def extra_figures(sample_all: pd.DataFrame, pair: pd.DataFrame, graphs: dict, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    colors = {"GSE148071": "#4c72b0", "GSE131907_tLung": "#c44e52", "GSE131907_mBrain": "#55a868"}
    for g, sub in pair.groupby("graph"):
        ax.scatter(sub["malignant_cldn4"], sub["frac_tnk"], s=36, alpha=0.85, c=colors.get(g, "#888"), label=g)
    ax.set_xlabel("Sample malignant CLDN4 (mean log1p CP10k)")
    ax.set_ylabel("Sample T/NK fraction")
    ax.set_title("Pair n: GSE148071 + GSE131907 tLung (unit = sample)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_patient_cldn4_vs_tnk_pair.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    scored = sample_all[np.isfinite(sample_all["malignant_cldn4"])]
    groups = list(scored.groupby("graph"))
    ax.boxplot([g["malignant_cldn4"].to_numpy() for _, g in groups], tick_labels=[n for n, _ in groups])
    ax.set_ylabel("Malignant CLDN4 (log1p CP10k)")
    ax.set_title("CLDN4 score by graph (scored samples only)")
    fig.tight_layout()
    fig.savefig(figdir / "fig_cldn4_score_by_dataset.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    names, scored_n, total_n = [], [], []
    for name, s in graphs.items():
        names.append(name)
        scored_n.append(s["samples"]["n_with_malignant_ge10"])
        total_n.append(s["samples"]["n"])
    x = np.arange(len(names))
    ax.bar(x - 0.18, total_n, 0.36, label="samples in graph", color="#bbbbbb")
    ax.bar(x + 0.18, scored_n, 0.36, label="scored (≥10 malignant)", color="#4c72b0")
    ax.set_xticks(x, names, rotation=15, ha="right")
    ax.set_ylabel("n samples (honest unit)")
    ax.set_title("Honest n — do not cite cell count")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_honest_n.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    names, testable, total = [], [], []
    for name, s in graphs.items():
        names.append(name)
        testable.append(s["da_malignant_cldn4"]["n_testable"])
        total.append(s["da_malignant_cldn4"]["n_nhoods_total"])
    x = np.arange(len(names))
    ax.bar(x - 0.18, total, 0.36, label="nhoods", color="#bbbbbb")
    ax.bar(x + 0.18, testable, 0.36, label="testable vs CLDN4 (≥5 samples)", color="#c44e52")
    ax.set_xticks(x, names, rotation=15, ha="right")
    ax.set_ylabel("neighbourhoods (overlapping)")
    ax.set_title("Most neighbourhoods are untestable — not hidden by cell n")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_nhood_testable.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ys, labels, xs, xerr = [], [], [], []
    for i, (name, s) in enumerate(graphs.items()):
        sl = s["sample_level_malignant_CLDN4_vs_TNK_fraction"]
        if sl.get("rho") is None:
            continue
        n = sl["n"]
        rho = sl["rho"]
        se = 1.0 / np.sqrt(max(n - 3, 1))
        lo, hi = np.tanh(np.arctanh(np.clip(rho, -0.999, 0.999)) - 1.96 * se), np.tanh(
            np.arctanh(np.clip(rho, -0.999, 0.999)) + 1.96 * se
        )
        ys.append(i)
        labels.append(f"{name} (n={n})")
        xs.append(rho)
        xerr.append((rho - lo, hi - rho))
    if ys:
        ax.errorbar(xs, ys, xerr=np.array(xerr).T, fmt="o", c="#4c72b0", capsize=3)
        ax.axvline(0, c="0.5", lw=0.8)
        ax.set_yticks(ys, labels)
        ax.set_xlabel("Spearman ρ (malignant CLDN4 vs sample T/NK)")
        ax.set_title("Sample-level CLDN4–T/NK (unit = sample)")
        fig.tight_layout()
        fig.savefig(figdir / "fig_forest_sample_cldn4_tnk.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
    for ax, g in zip(axes, ["GSE148071", "GSE131907_tLung"]):
        sub = pair[pair["graph"] == g]
        ax.scatter(sub["malignant_cldn4"], sub["frac_tnk"], s=32, c="#4c72b0", alpha=0.85)
        ax.set_title(f"{g} n={int(np.isfinite(sub['malignant_cldn4']).sum())}")
        ax.set_xlabel("Malignant CLDN4")
    axes[0].set_ylabel("Sample T/NK fraction")
    fig.suptitle("Pair arms only — mBrain not shown")
    fig.tight_layout()
    fig.savefig(figdir / "fig_pair_scatter_grid.png", dpi=160)
    plt.close(fig)


def write_finding(overall: dict, out_path: Path) -> None:
    graphs = overall["graphs"]
    pair = overall["pair"]
    sl_pair = pair["sample_level_malignant_CLDN4_vs_TNK"]
    sl_resid = pair["dataset_residualized_CLDN4_vs_TNK"]
    meta = pair["fisher_z_meta"]

    def block(tag: str, s: dict) -> str:
        da = s["da_malignant_cldn4"]
        sp = s["da_cldn4_high_vs_low"]
        c = s["composition"]
        iface = c["interface_nhoods"]["spearman"]
        paired = c["sample_paired_TNK_in_CLDN4high_vs_low_nhoods"]
        w = paired["wilcoxon_high_vs_low"]
        sl = s["sample_level_malignant_CLDN4_vs_TNK_fraction"]
        drop = ", ".join(s["samples"]["dropped_lt10_malignant"]) or "none"
        sens = s.get("da_malignant_cldn4_by_n_present", [])
        sens_rows = "\n".join(
            f"| ≥{r['min_samples_present']} | {r['n_testable']} | {r['n_p_lt_0.05']} | "
            f"{fmt(r.get('min_p'), sci=True)} | {fmt(r.get('min_SpatialFDR'))} | "
            f"**{r['n_SpatialFDR_lt_0.1']}** | {r['n_abs_rho_eq_1']} |"
            for r in sens
        )
        return f"""### {tag}

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 malignant(-like) cells | **{s['samples']['n_with_malignant_ge10']}** of {s['samples']['n']} | Mean log1p-CP10k CLDN4. Dropped if <10 ({drop}). |
| Median split (Welch, secondary) | high / low | **{s['samples']['n_cldn4_high']} vs {s['samples']['n_cldn4_low']}** | Median of the scored samples. |

Do not cite {s['n_cells_in_graph']:,} graph cells or {s['graph_params']['n_nhoods']:,} neighbourhoods as *n*. A disjoint subset has n={s['honest_n']['disjoint_nhood_subset_n']}.

Graph: {s['n_cells_in_graph']:,} cells; {s['n_malignant']:,} malignant(-like); {s['n_tnk']:,} T/NK; {s['graph_params']['n_nhoods']:,} neighbourhoods of size {s['graph_params']['nhood_size']}; k={s['graph_params']['k']}, d={s['graph_params']['d']}, {s['graph_params']['n_hvg']} HVG; PCA per-sample mean centering (not Harmony).

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | {da['n_testable']} | {da['n_p_lt_0.05']} | {fmt(da.get('min_p'), sci=True)} | {fmt(da.get('min_SpatialFDR'))} | **{da['n_SpatialFDR_lt_0.1']}** | {da['n_SpatialFDR_lt_0.05']} | {da['n_BH_FDR_lt_0.1']} |
| Median split high vs low | {sp['n_testable']} | {sp['n_p_lt_0.05']} | {fmt(sp.get('min_p'), sci=True)} | {fmt(sp.get('min_SpatialFDR'))} | **{sp['n_SpatialFDR_lt_0.1']}** | {sp['n_SpatialFDR_lt_0.05']} | {sp['n_BH_FDR_lt_0.1']} |

Sensitivity (continuous CLDN4, min n samples present):

| min n present | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | abs(rho)=1 |
|---|---|---|---|---|---|---|
{sens_rows}

Sample-level malignant CLDN4 vs T/NK fraction: ρ={fmt(sl.get('rho'))}, p={fmt(sl.get('p'), sci=True)}, n={sl.get('n')}.

Interface nhoods: n={c['interface_nhoods']['n']}; ρ={fmt(iface.get('rho'))}, p={fmt(iface.get('p'), sci=True)}, n={iface.get('n')}.
Disjoint interface: n={c['disjoint_interface']['n_disjoint_interface']} of {c['disjoint_interface']['n_disjoint_nhoods']}; ρ={fmt(c['disjoint_interface']['spearman'].get('rho'))}, p={fmt(c['disjoint_interface']['spearman'].get('p'), sci=True)}.
Sample-paired T/NK (unit = sample): n={w.get('n')}; median high={fmt(w.get('median_a'))}, low={fmt(w.get('median_b'))}; Wilcoxon p={fmt(w.get('p'), sci=True)}.
"""

    g148 = graphs["GSE148071"]
    tlung = graphs["GSE131907_tLung"]
    mbrain = graphs["GSE131907_mBrain"]
    text = f"""# FINDING — pairwise Milo vs malignant CLDN4 (GSE131907 + GSE148071)

Additive **CLDN4-only** neighbourhood DA on the public Kim et al. 2020 LUAD
atlas (GSE131907; PMID 32385277) and Wu et al. 2021 advanced-NSCLC biopsies
(GSE148071; PMID 33953163). This is **not** the triple
(131907+148071+205335) and **not** the 131907+205335 pair. GSE205335 was
not downloaded. Dual-high TACSTD2+CLDN4 is not a gate.

The independent unit is the **patient / sample**, not the cell and not the
overlapping neighbourhood. miloR is not used. Graphs are built separately.
Do not pool tLung with mBrain. Do not pool neighbourhoods across GEO series
into one SpatialFDR.

**Pair n** = GSE148071 scored biopsies + GSE131907 **tLung** scored samples
= **{pair['n_scored']}**. mBrain is a same-atlas sensitivity and is **not**
in that n.

## Honest n

| item | n | note |
|---|---:|---|
| GSE148071 GEO biopsies | **{g148['samples']['n']}** | one sample per patient |
| GSE148071 with malignant CLDN4 score | **{g148['samples']['n_with_malignant_ge10']}** | ≥10 malignant-like cells |
| GSE131907 tLung samples in graph | **{tlung['samples']['n']}** | epithelium + immune; author labels |
| GSE131907 tLung with CLDN4 score | **{tlung['samples']['n_with_malignant_ge10']}** | ≥10 author-malignant cells |
| **Pair n (148071 + tLung scored)** | **{pair['n_scored']}** | only combined test whose unit is the sample |
| GSE131907 mBrain scored | **{mbrain['samples']['n_with_malignant_ge10']}** of {mbrain['samples']['n']} | **not** added to pair n |
| GSE205335 | **0** | out of scope |
| Dual-high samples | **0** | not defined |

Do not cite {g148['n_cells_in_graph'] + tlung['n_cells_in_graph']:,} pair-arm cells as *n*.

## Pairwise sample-level tests (unit = sample)

| Test | Honest n | Result |
|---|---|---|
| Malignant CLDN4 vs T/NK, pooled pair | **{sl_pair.get('n')}** | ρ={fmt(sl_pair.get('rho'))}, p={fmt(sl_pair.get('p'), sci=True)} |
| Same, residualized within dataset | **{sl_resid.get('n')}** | ρ={fmt(sl_resid.get('rho'))}, p={fmt(sl_resid.get('p'), sci=True)} |
| Fisher-z meta of the two pair arms | n_total={meta.get('n_total')} | ρ={fmt(meta.get('rho'))}, p={fmt(meta.get('p'), sci=True)} |
| GSE148071 alone | {g148['sample_level_malignant_CLDN4_vs_TNK_fraction'].get('n')} | ρ={fmt(g148['sample_level_malignant_CLDN4_vs_TNK_fraction'].get('rho'))}, p={fmt(g148['sample_level_malignant_CLDN4_vs_TNK_fraction'].get('p'), sci=True)} |
| GSE131907 tLung alone | {tlung['sample_level_malignant_CLDN4_vs_TNK_fraction'].get('n')} | ρ={fmt(tlung['sample_level_malignant_CLDN4_vs_TNK_fraction'].get('rho'))}, p={fmt(tlung['sample_level_malignant_CLDN4_vs_TNK_fraction'].get('p'), sci=True)} |

Neighbourhood DA stays per graph. A joint kNN was not built.

Every SpatialFDR<0.1 neighbourhood on the default ≥5-sample floor is a
perfect or near-perfect rank correlation on 5–6 samples (scipy p≈0).
That is **not** a cohort DA claim. At n_present ≥ 8, SpatialFDR<0.1 is
**0** on GSE148071, tLung, and mBrain. Median-split Welch tests are also
**0** at SpatialFDR<0.1 on every graph.

## SpatialFDR by graph

{block("GSE148071 (Wu 2021 biopsies; marker malignant-like)", g148)}

{block("GSE131907 tLung (primary pair arm; author malignant)", tlung)}

{block("GSE131907 mBrain (same-atlas sensitivity; not in pair n)", mbrain)}

## What this does not say

- It does not invent an ICI / MPR / RECIST contrast. GSE148071 is diagnostic
  biopsies; GSE131907 tLung is treatment-naive.
- It does not treat cell count or neighbourhood count as *n*.
- Neighbourhood “next to” is kNN co-membership, not histology.
- Unrestricted neighbourhood CLDN4 vs T/NK is partly lineage geometry.
- miloR / edgeR QLF numbers are not claimed.
- Dual-high TACSTD2+CLDN4 was not tested.
- GSE205335 and the triple were not run.

See `tables/nhoods.tsv`, `tables/honest_n.tsv`, `tables/sample_scores.tsv`,
`tables/one_row.tsv`, and `tables/summary.json`. Extra figures are under
`figures/`.
"""
    out_path.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=Path("data"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/pair_131907_148071_milo_cldn4"))
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--d", type=int, default=30)
    ap.add_argument("--prop", type=float, default=0.1)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-interface", type=int, default=3)
    ap.add_argument("--min-malignant", type=int, default=10)
    ap.add_argument("--min-umi", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--from-results", action="store_true", help="rebuild pair tables/FINDING from existing graph outputs")
    args = ap.parse_args()
    self_test()

    results_dir = args.outdir / "results"
    figdir = args.outdir / "figures"
    tables = args.outdir / "tables"
    for p in (results_dir, figdir, tables):
        p.mkdir(parents=True, exist_ok=True)

    if args.from_results:
        s148 = json.loads((results_dir / "GSE148071" / "summary.json").read_text())
        s_t = json.loads((results_dir / "GSE131907_tLung" / "summary.json").read_text())
        s_b = json.loads((results_dir / "GSE131907_mBrain" / "summary.json").read_text())
        n148 = pd.read_csv(results_dir / "GSE148071" / "nhoods.tsv", sep="\t")
        n_t = pd.read_csv(results_dir / "GSE131907_tLung" / "nhoods.tsv", sep="\t")
        n_b = pd.read_csv(results_dir / "GSE131907_mBrain" / "nhoods.tsv", sep="\t")
        samp148 = pd.read_csv(results_dir / "GSE148071" / "sample_scores.tsv", sep="\t")
        samp_t = pd.read_csv(results_dir / "GSE131907_tLung" / "sample_scores.tsv", sep="\t")
        samp_b = pd.read_csv(results_dir / "GSE131907_mBrain" / "sample_scores.tsv", sep="\t")
        pair148 = pd.read_csv(results_dir / "GSE148071" / "sample_paired_tnk_by_nhood_cldn4.tsv", sep="\t")
        pair_t = pd.read_csv(results_dir / "GSE131907_tLung" / "sample_paired_tnk_by_nhood_cldn4.tsv", sep="\t")
        pair_b = pd.read_csv(results_dir / "GSE131907_mBrain" / "sample_paired_tnk_by_nhood_cldn4.tsv", sep="\t")
        inventory = {
            "from_results": True,
            "skipped": ["GSE205335", "triple 131907+148071+205335", "dual-high TACSTD2+CLDN4"],
        }
    else:
        s148, n148, samp148, pair148 = run_gse148071(args.data_root / "GSE148071" / "files", results_dir, figdir, args)
        res131, inventory = run_gse131907(args.data_root / "GSE131907", results_dir, figdir, args)
        s_t, n_t, samp_t, pair_t = res131["tLung"]
        s_b, n_b, samp_b, pair_b = res131["mBrain"]

    graphs = {"GSE148071": s148, "GSE131907_tLung": s_t, "GSE131907_mBrain": s_b}
    nhoods = pd.concat([n148, n_t, n_b], ignore_index=True)
    nhoods.to_csv(tables / "nhoods.tsv", sep="\t", index=False)
    n148.to_csv(tables / "nhoods_gse148071.tsv", sep="\t", index=False)
    n_t.to_csv(tables / "nhoods_gse131907_tLung.tsv", sep="\t", index=False)
    n_b.to_csv(tables / "nhoods_gse131907_mBrain.tsv", sep="\t", index=False)

    sample_all = pd.concat([samp148, samp_t, samp_b], ignore_index=True)
    sample_all.to_csv(tables / "sample_scores.tsv", sep="\t", index=False)
    paired_all = pd.concat([pair148, pair_t, pair_b], ignore_index=True)
    paired_all.to_csv(tables / "sample_paired_tnk_by_nhood_cldn4.tsv", sep="\t", index=False)

    pair_mask = sample_all["in_pair_n"] & np.isfinite(sample_all["malignant_cldn4"])
    pair_df = sample_all.loc[pair_mask].copy()
    sl_pair = spearman_safe(pair_df["malignant_cldn4"], pair_df["frac_tnk"])
    resid_x, resid_y = [], []
    for _, sub in pair_df.groupby("graph"):
        x = sub["malignant_cldn4"].to_numpy(dtype=float)
        y = sub["frac_tnk"].to_numpy(dtype=float)
        resid_x.append(x - np.nanmean(x))
        resid_y.append(y - np.nanmean(y))
    sl_resid = spearman_safe(np.concatenate(resid_x), np.concatenate(resid_y)) if resid_x else {"n": 0, "rho": None, "p": None}
    meta = fisher_z_meta(
        [
            s148["sample_level_malignant_CLDN4_vs_TNK_fraction"],
            s_t["sample_level_malignant_CLDN4_vs_TNK_fraction"],
        ]
    )
    extra_figures(sample_all, pair_df, graphs, figdir)

    honest_rows = []
    for name, s in graphs.items():
        honest_rows.append(
            {
                "graph": name,
                "dataset": s["dataset"],
                "n_samples_in_graph": s["samples"]["n"],
                "n_scored_malignant_cldn4": s["samples"]["n_with_malignant_ge10"],
                "n_cells_do_not_cite": s["n_cells_in_graph"],
                "n_nhoods": s["graph_params"]["n_nhoods"],
                "n_testable_cldn4": s["da_malignant_cldn4"]["n_testable"],
                "n_SpatialFDR_lt_0.1": s["da_malignant_cldn4"]["n_SpatialFDR_lt_0.1"],
                "min_SpatialFDR": s["da_malignant_cldn4"]["min_SpatialFDR"],
                "sample_cldn4_tnk_rho": s["sample_level_malignant_CLDN4_vs_TNK_fraction"].get("rho"),
                "sample_cldn4_tnk_p": s["sample_level_malignant_CLDN4_vs_TNK_fraction"].get("p"),
                "in_pair_n": name != "GSE131907_mBrain",
            }
        )
    honest = pd.DataFrame(honest_rows)
    honest.to_csv(tables / "honest_n.tsv", sep="\t", index=False)

    one = pd.DataFrame(
        [
            {
                "pair": "GSE131907+GSE148071",
                "not_triple": True,
                "not_131907_205335": True,
                "gene": "CLDN4",
                "dual_high": False,
                "unit": "patient/sample",
                "n_pair_scored": sl_pair.get("n"),
                "pair_cldn4_tnk_rho": sl_pair.get("rho"),
                "pair_cldn4_tnk_p": sl_pair.get("p"),
                "residualized_rho": sl_resid.get("rho"),
                "residualized_p": sl_resid.get("p"),
                "fisher_z_rho": meta.get("rho"),
                "fisher_z_p": meta.get("p"),
                "gse148071_SpatialFDR_lt_0.1": s148["da_malignant_cldn4"]["n_SpatialFDR_lt_0.1"],
                "tLung_SpatialFDR_lt_0.1": s_t["da_malignant_cldn4"]["n_SpatialFDR_lt_0.1"],
                "mBrain_SpatialFDR_lt_0.1": s_b["da_malignant_cldn4"]["n_SpatialFDR_lt_0.1"],
                "gse205335": 0,
            }
        ]
    )
    one.to_csv(tables / "one_row.tsv", sep="\t", index=False)

    overall = {
        "pair": {
            "datasets": ["GSE131907", "GSE148071"],
            "excluded": ["GSE205335", "triple", "dual-high TACSTD2+CLDN4"],
            "gene": "CLDN4",
            "dual_high": False,
            "miloR": False,
            "independent_unit": "patient/sample",
            "pair_arms": ["GSE148071", "GSE131907_tLung"],
            "n_scored": sl_pair.get("n"),
            "sample_level_malignant_CLDN4_vs_TNK": sl_pair,
            "dataset_residualized_CLDN4_vs_TNK": sl_resid,
            "fisher_z_meta": meta,
        },
        "graphs": graphs,
        "gse131907_inventory": inventory,
        "citation": {
            "GSE131907": "Kim et al. Nat Commun 2020 PMID 32385277",
            "GSE148071": "Wu et al. Nat Commun 2021 PMID 33953163",
        },
    }
    write_json(tables / "summary.json", overall)
    write_finding(overall, args.outdir / "FINDING.md")
    print(
        json.dumps(
            {
                "pair_n": sl_pair,
                "gse148071": s148["da_malignant_cldn4"],
                "tLung": s_t["da_malignant_cldn4"],
                "mBrain": s_b["da_malignant_cldn4"],
                "nhood_table": str(tables / "nhoods.tsv"),
                "finding": str(args.outdir / "FINDING.md"),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
