#!/usr/bin/env python3
"""Triple-merge Milo-style kNN neighbourhood DA vs malignant CLDN4.

Additive CLDN4-only. Public processed files from GSE131907 + GSE148071 +
GSE205335. No dual-high gate. No GSE207422. Patient/sample is the unit.

Graphs are built per dataset (per-sample PCA centering). Joint Harmony of
the three public matrices is not run: a concatenated log-normalized HVG
matrix of ~370k cells does not fit in 15 GB. miloR / edgeR are not used.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import time
import warnings
from collections import defaultdict
from dataclasses import dataclass
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

# --- lineage / marker constants ------------------------------------------------

GSE131907_MALIGNANT = {"Malignant cells", "tS1", "tS2", "tS3"}
GSE131907_TNK = {"T lymphocytes", "NK cells"}
GSE131907_KEEP = {
    "T lymphocytes",
    "NK cells",
    "B lymphocytes",
    "Myeloid cells",
    "MAST cells",
    "Epithelial cells",
}
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
ALWAYS = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD3D", "NKG7"]
HEADER_GENE = {"", "gene", "Gene", "GENE", "index", "Index", "symbol", "Symbol", "SYMBOL"}
RESPONSE_MAP = {"PR": "R", "CR": "R", "SD": "NR", "PD": "NR", "NE": "NE"}
DS_COLORS = {"GSE131907": "#4c72b0", "GSE148071": "#dd8452", "GSE205335": "#55a868"}


@dataclass
class Packed:
    dataset: str
    X: np.ndarray
    sample: np.ndarray
    lineage: np.ndarray
    is_malig: np.ndarray
    is_tnk: np.ndarray
    cldn4: np.ndarray
    notes: dict


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


def fdr_counts(df: pd.DataFrame) -> dict:
    t = df[df["testable"]] if "testable" in df.columns else df.iloc[0:0]
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
        "median_n_samples_present": float(t["n_samples_present"].median()) if len(t) and "n_samples_present" in t else None,
        "n_abs_rho_eq_1": int((t["spearman_rho"].abs() >= 0.999).sum()) if len(t) and "spearman_rho" in t else 0,
    }


def _fmt(x, nd=3, sci=False):
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
        alpha=0.35,
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


def plot_scatter(x, y, path: Path, xlabel: str, ylabel: str, title: str, c="#4c72b0") -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.scatter(x, y, s=10, c=c, alpha=0.45)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


# --- shared DA -----------------------------------------------------------------


def run_one_graph(pack: Packed, outdir: Path, args) -> dict:
    dest = outdir / pack.dataset
    dest.mkdir(parents=True, exist_ok=True)
    X, sample = pack.X, pack.sample
    print(
        f"{pack.dataset}: cells={X.shape[0]} HVG={X.shape[1]} samples={pd.unique(sample).size}",
        flush=True,
    )
    pcs, knn_idx, knn_dist = pca_knn(X, sample, n_pcs=args.d, k=args.k, random_state=args.seed)
    indices = refine_indices(pcs, knn_idx, prop=args.prop, random_state=args.seed)
    nhoods = make_nhoods(knn_idx, knn_dist, indices)
    members = nhood_membership(nhoods, X.shape[0])
    print(
        f"{pack.dataset}: nhoods={len(members)} median_size={np.median([m.size for m in members]):.0f}",
        flush=True,
    )

    sample_levels = sorted(pd.unique(sample), key=lambda s: (str(s),))
    counts = count_matrix(members, sample, sample_levels)
    sizes = sample_sizes(sample, sample_levels)

    mal_score, tnk_frac, n_mal_s, n_tnk_s = [], [], [], []
    for s in sample_levels:
        m = sample == s
        n_mal = int((m & pack.is_malig).sum())
        n_tnk = int((m & pack.is_tnk).sum())
        n_mal_s.append(n_mal)
        n_tnk_s.append(n_tnk)
        mal_score.append(float(pack.cldn4[m & pack.is_malig].mean()) if n_mal >= args.min_malignant else np.nan)
        denom = n_mal + n_tnk
        tnk_frac.append(float(n_tnk / denom) if denom else float(pack.is_tnk[m].mean()))
    mal_score = np.asarray(mal_score, dtype=float)
    tnk_frac = np.asarray(tnk_frac, dtype=float)

    da_cldn = attach_fdr(spearman_da(counts, sizes, mal_score, min_samples=5), nhoods.k_distance)
    scored = np.isfinite(mal_score)
    group = np.full(len(sample_levels), "unscored", dtype=object)
    med_s = float("nan")
    if scored.sum() >= 4:
        med_s = float(np.median(mal_score[scored]))
        group[scored & (mal_score > med_s)] = "high"
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

    comp = nhood_composition(members, pack.is_tnk, pack.is_malig, pack.cldn4, pack.lineage)
    comp = comp.rename(columns={"tacstd2_all_mean": "cldn4_all_mean", "tacstd2_malig_mean": "cldn4_malig_mean"})
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
    merged.insert(0, "dataset", pack.dataset)
    merged.to_csv(dest / "nhoods.tsv", sep="\t", index=False)
    da_cldn.to_csv(dest / "da_malignant_cldn4.tsv", sep="\t", index=False)
    da_split.to_csv(dest / "da_cldn4_median_split.tsv", sep="\t", index=False)

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
    med_n = float(np.nanmedian(comp.loc[has_mal, "cldn4_malig_mean"])) if has_mal.any() else np.nan
    high = has_mal.to_numpy() & (comp["cldn4_malig_mean"].to_numpy() >= med_n)
    low = has_mal.to_numpy() & (comp["cldn4_malig_mean"].to_numpy() < med_n)
    paired = sample_paired_tnk_by_tacstd2(members, sample, pack.is_tnk, comp["cldn4_malig_mean"].to_numpy(), high, low)
    sample_tab = pd.DataFrame(
        {
            "dataset": pack.dataset,
            "sample": sample_levels,
            "n_cells": sizes.astype(int),
            "n_malignant": n_mal_s,
            "n_tnk": n_tnk_s,
            "malignant_cldn4": mal_score,
            "frac_tnk": tnk_frac,
            "cldn4_arm": group,
        }
    )
    paired = paired.merge(sample_tab, on="sample")
    paired.to_csv(dest / "sample_paired_tnk_by_nhood_cldn4.tsv", sep="\t", index=False)
    sample_tab.to_csv(dest / "sample_scores.tsv", sep="\t", index=False)
    paired_test = wilcoxon_paired(paired["frac_tnk_high_nhood"], paired["frac_tnk_low_nhood"])
    sample_level = spearman_safe(mal_score, tnk_frac)

    def _sens(min_n: int) -> dict:
        t = da_cldn.loc[da_cldn["testable"] & (da_cldn["n_samples_present"] >= min_n)].copy()
        t["testable"] = True
        rec = fdr_counts(t)
        rec["min_samples_present"] = min_n
        return rec

    floor_hits = da_cldn[(da_cldn["testable"]) & (da_cldn["SpatialFDR"] < 0.1)]
    da_counts = fdr_counts(da_cldn)

    plot_volcano(
        da_cldn.rename(columns={"spearman_rho": "stat"}),
        "stat",
        f"{pack.dataset}: nhood abundance vs malignant CLDN4",
        dest / "fig_da_cldn4_volcano.png",
        "Spearman ρ (nhood proportion vs sample malignant CLDN4)",
    )
    if da_split["testable"].any():
        plot_volcano(
            da_split,
            "logFC_B_minus_A",
            f"{pack.dataset}: nhood DA CLDN4-high vs low samples",
            dest / "fig_da_cldn4_split_volcano.png",
            "log2 FC (high − low sample proportion)",
        )
    plot_scatter(
        comp.loc[iface, "cldn4_malig_mean"],
        comp.loc[iface, "frac_tnk"],
        dest / "fig_interface_cldn4_vs_tnk.png",
        "Neighbourhood malignant CLDN4 (mean log1p CP10k)",
        "Neighbourhood T/NK fraction",
        f"{pack.dataset} interface nhoods (n={int(iface.sum())}); transcriptional, not spatial",
        c=DS_COLORS.get(pack.dataset, "#4c72b0"),
    )
    if paired["frac_tnk_high_nhood"].notna().any():
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        for _, r in paired.iterrows():
            if np.isfinite(r["frac_tnk_high_nhood"]) and np.isfinite(r["frac_tnk_low_nhood"]):
                ax.plot(
                    [0, 1],
                    [r["frac_tnk_low_nhood"], r["frac_tnk_high_nhood"]],
                    "-o",
                    c=DS_COLORS.get(pack.dataset, "#4c72b0"),
                    alpha=0.7,
                    ms=5,
                )
        ax.set_xticks([0, 1], ["CLDN4-low nhoods", "CLDN4-high nhoods"])
        ax.set_ylabel("Sample T/NK fraction among cells in those nhoods")
        ax.set_title(f"{pack.dataset} sample-paired (unit = sample)")
        fig.tight_layout()
        fig.savefig(dest / "fig_sample_paired_tnk.png", dpi=160)
        plt.close(fig)

    dropped = sample_tab.loc[~np.isfinite(sample_tab["malignant_cldn4"])]
    drop_bits = [f"{r.sample}={int(r.n_malignant)}" for r in dropped.itertuples()]
    summary = {
        "dataset": pack.dataset,
        "notes": pack.notes,
        "n_cells_in_graph": int(X.shape[0]),
        "n_malignant": int(pack.is_malig.sum()),
        "n_tnk": int(pack.is_tnk.sum()),
        "lineage_counts": {k: int(v) for k, v in pd.Series(pack.lineage).value_counts().items()},
        "graph": {
            "k": args.k,
            "d": min(args.d, X.shape[1], X.shape[0] - 1),
            "prop": args.prop,
            "n_hvg": int(X.shape[1]),
            "n_nhoods": int(len(members)),
            "nhood_size": args.k + 1,
            "median_nhood_size": float(np.median([m.size for m in members])),
            "batch": "PCA per-sample mean centering (not Harmony); graph per dataset",
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
        "da_malignant_cldn4": da_counts,
        "da_cldn4_median_split": fdr_counts(da_split),
        "da_malignant_cldn4_by_n_present": [_sens(t) for t in (5, 8, 10, 15)],
        "spatialfdr_floor_hits": {
            "n": int(len(floor_hits)),
            "all_abs_rho_eq_1": bool(len(floor_hits) and (floor_hits["spearman_rho"].abs() >= 0.999).all()),
            "n_samples_present_min": int(floor_hits["n_samples_present"].min()) if len(floor_hits) else None,
            "n_samples_present_max": int(floor_hits["n_samples_present"].max()) if len(floor_hits) else None,
            "n_tnk_max": int(comp.set_index("nhood").loc[floor_hits["nhood"], "n_tnk"].max()) if len(floor_hits) else None,
        },
        "composition": {
            "note": "transcriptional kNN != spatial niche",
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
                "median_cut": med_n,
                "n_high_nhoods": int(high.sum()),
                "n_low_nhoods": int(low.sum()),
                "wilcoxon_high_vs_low": paired_test,
            },
        },
        "honest_n": {
            "independent_unit_DA": f"{pack.dataset} samples (n={len(sample_levels)}; scored={int(np.isfinite(mal_score).sum())})",
            "nhoods_are_overlapping": True,
            "SpatialFDR_is_overlap_aware": True,
            "disjoint_nhood_subset_n": int(indep.size),
            "do_not_cite_n_cells_as_n": int(X.shape[0]),
        },
    }
    write_json(dest / "summary.json", summary)
    return {
        "summary": summary,
        "nhoods": merged,
        "da_cldn": da_cldn,
        "da_split": da_split,
        "sample_tab": sample_tab,
        "paired": paired,
        "comp": comp,
        "iface": iface,
    }


# --- GSE131907 ----------------------------------------------------------------


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
    return "unknown"


def load_gse131907_annotation(ann_path: Path, cell_ids: list[str]) -> pd.DataFrame:
    ann = pd.read_csv(ann_path, sep="\t", dtype=str)
    if "Index" not in ann.columns:
        raise SystemExit(f"GSE131907 annotation missing Index: {list(ann.columns)}")
    ann = ann.set_index("Index").reindex(cell_ids).reset_index()
    if ann["Sample"].isna().any():
        raise SystemExit(f"{int(ann['Sample'].isna().sum())} matrix cell IDs missing from annotation")
    if "Sample_Origin" in ann.columns:
        origin = ann["Sample_Origin"].fillna("")
    else:
        origin = pd.Series("", index=ann.index)
    mapped = ann["Sample"].map(origin_from_sample)
    known = {"tLung", "nLung", "tL/B", "mLN", "nLN", "PE", "mBrain"}
    use_author = origin.isin(known)
    ann["origin"] = np.where(use_author, origin, mapped)
    return ann


def load_gse131907(datadir: Path, n_hvg: int) -> Packed:
    ann_path = datadir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi_path = datadir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    for p in (ann_path, umi_path):
        if not p.exists():
            raise SystemExit(f"missing {p}; run scripts/download.py")
    t0 = time.time()
    with gzip.open(umi_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    cell_ids = header[1:]
    print(f"GSE131907 header cells={len(cell_ids)}", flush=True)
    ann = load_gse131907_annotation(ann_path, cell_ids)
    cell_type = ann["Cell_type"].to_numpy()
    subtype = ann["Cell_subtype"].fillna("").to_numpy() if "Cell_subtype" in ann.columns else np.array([""] * len(ann))
    origin = ann["origin"].to_numpy()
    sample_all = ann["Sample"].to_numpy()
    keep_type = np.isin(cell_type, list(GSE131907_KEEP))
    is_tlung = origin == "tLung"
    keep = keep_type & is_tlung
    keep_idx = np.flatnonzero(keep)
    print(
        f"GSE131907 tLung epithelium+immune {int(keep.sum())} / {len(cell_ids)} "
        f"(tLung={int(is_tlung.sum())})",
        flush=True,
    )
    markers = wanted_markers()
    n_all = len(cell_ids)
    n_umi = np.zeros(n_all, dtype=np.float64)
    names: list[str] = []
    means: list[float] = []
    varis: list[float] = []
    found: dict[str, np.ndarray] = {}
    with gzip.open(umi_path, "rt") as handle:
        handle.readline()
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_all:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_all}")
            n_umi += arr
            names.append(gene)
            sub = arr[keep]
            means.append(float(sub.mean()) if sub.size else 0.0)
            varis.append(float(sub.var()) if sub.size else 0.0)
            if gene in markers:
                found[gene] = arr
            if i % 2000 == 0:
                print(f"  GSE131907 pass1 {i} genes markers={len(found)} [{time.time() - t0:.0f}s]", flush=True)
    hvg_idx = select_hvg(np.asarray(means), np.asarray(varis), n_hvg=n_hvg)
    hvg_names = [names[i] for i in hvg_idx]
    keep_genes = set(hvg_names) | markers
    print(f"GSE131907 HVG={len(hvg_names)} keep_genes={len(keep_genes)}", flush=True)
    expr: dict[str, np.ndarray] = {}
    with gzip.open(umi_path, "rt") as handle:
        handle.readline()
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            if gene not in keep_genes:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            expr[gene] = arr[keep]
            if i % 4000 == 0:
                print(f"  GSE131907 pass2 stored {len(expr)} [{time.time() - t0:.0f}s]", flush=True)
    n_umi_k = n_umi[keep]
    scale = np.where(n_umi_k > 0, 1e4 / n_umi_k, 0.0)
    hvg_present = [g for g in hvg_names if g in expr]
    X = np.vstack([np.log1p(expr[g] * scale) for g in hvg_present]).T.astype(np.float32)
    if "CLDN4" not in expr:
        raise RuntimeError("CLDN4 missing from GSE131907 UMI matrix")
    cldn4 = np.log1p(expr["CLDN4"] * scale).astype(np.float32)
    lineage = cell_type[keep]
    is_malig = np.isin(subtype[keep], list(GSE131907_MALIGNANT))
    is_tnk = np.isin(lineage, list(GSE131907_TNK))
    sample = sample_all[keep]
    return Packed(
        dataset="GSE131907",
        X=X,
        sample=sample,
        lineage=lineage,
        is_malig=is_malig,
        is_tnk=is_tnk,
        cldn4=cldn4,
        notes={
            "citation": "Kim et al. Nat Commun 2020 PMID 32385277",
            "graph_cells": "tLung epithelium + immune; author Cell_type",
            "malignant_definition": "author Cell_subtype in {Malignant cells, tS1, tS2, tS3}",
            "author_barcode_labels": True,
            "ici_labels": False,
            "n_cells_matrix": n_all,
            "n_tlung": int(is_tlung.sum()),
            "sites_excluded": "nLung, nLN, mLN, PE, mBrain, tL/B (not mixed into the tLung graph)",
        },
    )


# --- GSE148071 ----------------------------------------------------------------


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


def stream_gse148071_file(path: Path, markers: set[str], keep: set[str] | None):
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


def load_gse148071(datadir: Path, n_hvg: int, min_umi: float = 1.0) -> Packed:
    files = sorted((datadir / "files").glob("*_exp.txt.gz"))
    if not files:
        files = sorted((datadir / "files").rglob("*_exp.txt.gz"))
    if not files:
        raise SystemExit(f"no *_exp.txt.gz under {datadir}/files; run scripts/download.py")
    markers = wanted_markers()
    print(f"GSE148071 pass1: {len(files)} files", flush=True)
    sample_all: list[str] = []
    n_umi_all: list[np.ndarray] = []
    gene_sum: dict[str, float] = {}
    gene_sumsq: dict[str, float] = {}
    n_cells_total = 0
    for fp in files:
        patient = patient_from_name(fp)
        cell_ids, n_umi, found, gsum, gss, n_genes = stream_gse148071_file(fp, markers, keep=None)
        n = len(cell_ids)
        sample_all.extend([patient] * n)
        n_umi_all.append(n_umi)
        n_cells_total += n
        for g, val in gsum.items():
            gene_sum[g] = gene_sum.get(g, 0.0) + val
            gene_sumsq[g] = gene_sumsq.get(g, 0.0) + gss[g]
        print(f"  {fp.name} {patient} n={n} genes={n_genes} markers={len(found)}", flush=True)
    n_umi = np.concatenate(n_umi_all)
    sample = np.array(sample_all, dtype=object)
    keep_cell = n_umi >= min_umi
    genes = sorted(gene_sum)
    means = np.array([gene_sum[g] / n_cells_total for g in genes], dtype=float)
    varis = np.array([gene_sumsq[g] / n_cells_total - means[i] ** 2 for i, g in enumerate(genes)], dtype=float)
    varis = np.maximum(varis, 0.0)
    hvg_idx = select_hvg(means, varis, n_hvg=n_hvg)
    hvg_names = [genes[i] for i in hvg_idx]
    keep_genes = set(hvg_names) | markers
    print(f"GSE148071 HVG={len(hvg_names)} cells={n_cells_total}", flush=True)
    expr_blocks: dict[str, list[np.ndarray]] = defaultdict(list)
    for fp in files:
        _, _, found, _, _, _ = stream_gse148071_file(fp, markers, keep=keep_genes)
        n = next(iter(found.values())).size if found else 0
        # n from file header if found empty
        if n == 0:
            with gzip.open(fp, "rt") as handle:
                n = len(parse_header(handle.readline()))
        present = set(found)
        for g in keep_genes:
            expr_blocks[g].append(found[g] if g in present else np.zeros(n, dtype=np.float32))
        print(f"  pass2 {fp.name} stored {len(present)}/{len(keep_genes)}", flush=True)
    expr = {g: np.concatenate(blocks) for g, blocks in expr_blocks.items()}
    scale = 1e4 / np.maximum(n_umi, 1.0)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
    n = int(n_umi.size)
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
        scale = scale[keep_cell]
        expr = {g: v[keep_cell] for g, v in expr.items()}
        n = int(keep_cell.sum())
    hvg_present = [g for g in hvg_names if g in expr]
    X = np.vstack([np.log1p(expr[g] * scale) for g in hvg_present]).T.astype(np.float32)
    return Packed(
        dataset="GSE148071",
        X=X,
        sample=sample,
        lineage=lineage,
        is_malig=is_malig,
        is_tnk=is_tnk,
        cldn4=cldn4,
        notes={
            "citation": "Wu et al. Nat Commun 2021 PMID 33953163",
            "author_barcode_labels": False,
            "lineage_source": "marker_argmax_canonical",
            "malignant_definition": "epithelial AND NOT (alveolar/club/ciliated log1p-CP10k >= 1); not inferCNV",
            "n_cells_before_umi_filter": n_cells_total,
            "n_geo_files": len(files),
        },
    )


# --- GSE205335 ----------------------------------------------------------------


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


def parse_geo_soft(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with gzip.open(path, "rt", errors="replace") as handle:
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
    return metadata.rename(columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"})


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
    print(f"GSE205335 matrix {matrix.shape} nnz={matrix.nnz} load={time.time() - t0:.1f}s", flush=True)
    return matrix, genes, barcodes


def gene_mean_var(matrix: sparse.spmatrix) -> tuple[np.ndarray, np.ndarray]:
    n = matrix.shape[1]
    mean = np.asarray(matrix.mean(axis=1)).ravel()
    mean_sq = np.asarray(matrix.multiply(matrix).mean(axis=1)).ravel()
    var = np.maximum(mean_sq - mean * mean, 0.0) * (n / max(n - 1, 1))
    return mean, var


def load_gse205335(datadir: Path, n_hvg: int) -> Packed:
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
    cells["is_normal_tissue"] = cells["tissue"].astype(str).str.startswith("Normal")
    cells["is_malig"] = cells["lineage.sub"].eq("Malignant cells")
    cells["is_tnk"] = cells["lineage.total"].eq("T/NK cells")
    gene_idx = {g: i for i, g in enumerate(genes)}
    if "CLDN4" not in gene_idx:
        raise RuntimeError("CLDN4 missing from GSE205335 UMI matrix")
    gene_mean, gene_var = gene_mean_var(matrix)
    hvg_idx = select_hvg(gene_mean, gene_var, n_hvg=n_hvg)
    keep_genes = sorted(set(hvg_idx.tolist()) | {gene_idx[g] for g in ALWAYS if g in gene_idx})
    lib = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    scale = np.where(lib > 0, 1e4 / lib, 0.0)
    sub = matrix[keep_genes, :].astype(np.float32)
    del matrix
    dense = sub.toarray()
    del sub
    log_all = np.log1p(dense * scale).astype(np.float32)
    name_of = {keep_genes[i]: i for i in range(len(keep_genes))}
    cldn4 = log_all[name_of[gene_idx["CLDN4"]]]
    hvg_rows = [name_of[i] for i in hvg_idx if i in name_of]
    X_all = log_all[hvg_rows].T.copy()
    del dense, log_all
    keep = ~cells["is_normal_tissue"].to_numpy()
    idx = np.flatnonzero(keep)
    print(f"GSE205335 non-normal {int(keep.sum())} / {len(cells)}", flush=True)
    return Packed(
        dataset="GSE205335",
        X=X_all[idx],
        sample=cells.loc[keep, "patient"].to_numpy(dtype=object),
        lineage=cells.loc[keep, "lineage.total"].to_numpy(dtype=object),
        is_malig=cells.loc[keep, "is_malig"].to_numpy(),
        is_tnk=cells.loc[keep, "is_tnk"].to_numpy(),
        cldn4=cldn4[idx],
        notes={
            "citation": "Ahn / Lee et al. eLife 2024 (GEO GSE205335)",
            "author_barcode_labels": True,
            "malignant_definition": "author lineage.sub == 'Malignant cells'; CNV not re-inferred",
            "tnk_definition": "author lineage.total == 'T/NK cells'",
            "independent_unit": "patient (biopsies from the same patient pooled)",
            "normal_cells_excluded": int((~keep).sum()),
            "n_cells_matrix": int(len(cells)),
            "mpr_labeled": False,
        },
    )


# --- combine / figures / FINDING ----------------------------------------------


def extra_figures(results: dict[str, dict], outdir: Path, stacked: pd.DataFrame, patients: pd.DataFrame) -> None:
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    # 3-panel volcano
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0), sharey=True)
    for ax, ds in zip(axes, ["GSE131907", "GSE148071", "GSE205335"]):
        da = results[ds]["da_cldn"]
        testable = da["testable"].to_numpy()
        x = da["spearman_rho"]
        y = -np.log10(np.clip(da["p"], EPS, 1))
        ax.scatter(x[~testable], y[~testable], s=5, c="#bbbbbb", alpha=0.3)
        t = da[testable]
        sig = t["SpatialFDR"] < 0.1
        ax.scatter(t.loc[~sig, "spearman_rho"], -np.log10(np.clip(t.loc[~sig, "p"], EPS, 1)), s=7, c=DS_COLORS[ds], alpha=0.55)
        ax.scatter(t.loc[sig, "spearman_rho"], -np.log10(np.clip(t.loc[sig, "p"], EPS, 1)), s=12, c="#c44e52", alpha=0.9)
        ax.axhline(-np.log10(0.05), ls="--", c="0.5", lw=0.7)
        ax.set_title(ds)
        ax.set_xlabel("Spearman ρ")
    axes[0].set_ylabel("−log10 p")
    fig.suptitle("Neighbourhood abundance vs malignant CLDN4 (per-dataset graph)", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "fig_triple_volcano.png", dpi=160)
    plt.close(fig)

    # patient CLDN4 vs T/NK
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for ds, sub in patients.groupby("dataset"):
        ok = sub["malignant_cldn4"].notna() & sub["frac_tnk"].notna()
        ax.scatter(
            sub.loc[ok, "malignant_cldn4"],
            sub.loc[ok, "frac_tnk"],
            s=36,
            c=DS_COLORS.get(ds, "#333"),
            label=f"{ds} n={int(ok.sum())}",
            edgecolors="white",
            linewidths=0.4,
        )
    ax.set_xlabel("Sample malignant CLDN4 (mean log1p CP10k)")
    ax.set_ylabel("T/NK fraction (malignant + T/NK denom when available)")
    ax.set_title("Patient/sample-level; unit = sample (not cell)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_patient_cldn4_vs_tnk.png", dpi=160)
    plt.close(fig)

    # honest n bars
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    labels, scored, total = [], [], []
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        s = results[ds]["summary"]["samples"]
        labels.append(ds)
        scored.append(s["n_with_malignant_ge10"])
        total.append(s["n"])
    x = np.arange(len(labels))
    ax.bar(x - 0.18, total, 0.36, color="#bbbbbb", label="samples in graph")
    ax.bar(x + 0.18, scored, 0.36, color="#4c72b0", label="scored (≥10 malignant)")
    ax.set_xticks(x, labels)
    ax.set_ylabel("n samples (the DA unit)")
    ax.set_title("Honest n — do not cite cell count")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_honest_n.png", dpi=160)
    plt.close(fig)

    # testable fraction
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    nhoods, testable = [], []
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        d = results[ds]["summary"]["da_malignant_cldn4"]
        nhoods.append(d["n_nhoods_total"])
        testable.append(d["n_testable"])
    ax.bar(x - 0.18, nhoods, 0.36, color="#bbbbbb", label="all nhoods")
    ax.bar(x + 0.18, testable, 0.36, color="#dd8452", label="testable (≥5 samples)")
    ax.set_xticks(x, labels)
    ax.set_ylabel("neighbourhoods (overlapping)")
    ax.set_title("Most neighbourhoods are sample-private")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_testable_nhoods.png", dpi=160)
    plt.close(fig)

    # SpatialFDR summary
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    mins = []
    hits = []
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        d = results[ds]["summary"]["da_malignant_cldn4"]
        mins.append(d["min_SpatialFDR"] if d["min_SpatialFDR"] is not None else 1.0)
        hits.append(d["n_SpatialFDR_lt_0.1"])
    ax.bar(x, mins, color=[DS_COLORS[d] for d in labels])
    ax.axhline(0.1, ls="--", c="0.3", lw=0.8, label="SpatialFDR 0.1")
    ax.set_xticks(x, labels)
    ax.set_ylabel("min SpatialFDR")
    ax.set_title("DA vs malignant CLDN4 — min SpatialFDR per graph")
    ax.set_ylim(0, 1.05)
    for i, h in enumerate(hits):
        ax.text(i, min(mins[i] + 0.03, 1.0), f"{h} hits <0.1", ha="center", fontsize=8)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_min_spatialfdr.png", dpi=160)
    plt.close(fig)

    # stacked interface
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        comp = results[ds]["comp"]
        iface = results[ds]["iface"]
        ax.scatter(
            comp.loc[iface, "cldn4_malig_mean"],
            comp.loc[iface, "frac_tnk"],
            s=10,
            c=DS_COLORS[ds],
            alpha=0.45,
            label=f"{ds} n={int(iface.sum())}",
        )
    ax.set_xlabel("Neighbourhood malignant CLDN4 (mean log1p CP10k)")
    ax.set_ylabel("Neighbourhood T/NK fraction")
    ax.set_title("Interface nhoods (overlapping); transcriptional, not spatial")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_interface_cldn4_vs_tnk.png", dpi=160)
    plt.close(fig)

    # combined paired T/NK
    fig, ax = plt.subplots(figsize=(5.8, 4.6))
    for ds, sub in patients.groupby("dataset"):
        for _, r in results[ds]["paired"].iterrows():
            if np.isfinite(r["frac_tnk_high_nhood"]) and np.isfinite(r["frac_tnk_low_nhood"]):
                ax.plot(
                    [0, 1],
                    [r["frac_tnk_low_nhood"], r["frac_tnk_high_nhood"]],
                    "-o",
                    c=DS_COLORS[ds],
                    alpha=0.35,
                    ms=4,
                )
    ax.set_xticks([0, 1], ["CLDN4-low nhoods", "CLDN4-high nhoods"])
    ax.set_ylabel("Sample T/NK fraction among cells in those nhoods")
    ax.set_title("Sample-paired (geometry, not a spatial niche)")
    fig.tight_layout()
    fig.savefig(figdir / "fig_sample_paired_tnk.png", dpi=160)
    plt.close(fig)

    # rho distribution of testable nhoods
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        da = results[ds]["da_cldn"]
        t = da.loc[da["testable"], "spearman_rho"].dropna()
        if len(t):
            ax.hist(t, bins=30, histtype="step", density=True, color=DS_COLORS[ds], label=ds, lw=1.4)
    ax.axvline(0, c="0.4", lw=0.7)
    ax.set_xlabel("Spearman ρ (testable nhoods)")
    ax.set_ylabel("density")
    ax.set_title("Testable neighbourhood ρ vs malignant CLDN4")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_rho_hist_testable.png", dpi=160)
    plt.close(fig)


def write_finding(overall: dict, out_path: Path) -> None:
    rows = []
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        s = overall["per_dataset"][ds]
        da = s["da_malignant_cldn4"]
        sl = s["sample_level_malignant_CLDN4_vs_TNK_fraction"]
        c = s["composition"]
        rows.append(
            f"| {ds} | **{s['samples']['n_with_malignant_ge10']}** / {s['samples']['n']} | "
            f"{da['n_testable']} / {da['n_nhoods_total']} | {da['n_p_lt_0.05']} | "
            f"{_fmt(da.get('min_SpatialFDR'))} | **{da['n_SpatialFDR_lt_0.1']}** | "
            f"{_fmt(sl.get('rho'))} ({_fmt(sl.get('p'), sci=True)}) |"
        )
    stacked = overall["stacked_patient_CLDN4_vs_TNK"]
    n_scored = overall["honest_n"]["n_scored_samples"]
    n_cells = overall["honest_n"]["do_not_cite_n_cells"]
    floor_notes = []
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        fl = overall["per_dataset"][ds]["spatialfdr_floor_hits"]
        if fl["n"]:
            floor_notes.append(
                f"{ds}: {fl['n']} SpatialFDR<0.1 nhoods"
                + (
                    f" (all |ρ|=1 on {fl['n_samples_present_min']}–{fl['n_samples_present_max']} samples; T/NK max={fl['n_tnk_max']})"
                    if fl.get("all_abs_rho_eq_1")
                    else ""
                )
            )
    floor_txt = " ".join(floor_notes) if floor_notes else "No SpatialFDR<0.1 neighbourhoods in any graph."
    sens_lines = []
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        for rec in overall["per_dataset"][ds]["da_malignant_cldn4_by_n_present"]:
            sens_lines.append(
                f"| {ds} | ≥{rec['min_samples_present']} | {rec['n_testable']} | {rec['n_p_lt_0.05']} | "
                f"{_fmt(rec.get('min_SpatialFDR'))} | **{rec['n_SpatialFDR_lt_0.1']}** | {rec['n_abs_rho_eq_1']} |"
            )
    iface_lines = []
    for ds in ["GSE131907", "GSE148071", "GSE205335"]:
        c = overall["per_dataset"][ds]["composition"]
        iface = c["interface_nhoods"]["spearman"]
        di = c["disjoint_interface"]
        w = c["sample_paired_TNK_in_CLDN4high_vs_low_nhoods"]["wilcoxon_high_vs_low"]
        iface_lines.append(
            f"| {ds} | {c['interface_nhoods']['n']} | {_fmt(iface.get('rho'))} ({_fmt(iface.get('p'), sci=True)}) | "
            f"{di['n_disjoint_interface']} / {di['n_disjoint_nhoods']} | "
            f"{_fmt(di['spearman'].get('rho'))} ({_fmt(di['spearman'].get('p'), sci=True)}) | "
            f"{w.get('n')}; med {_fmt(w.get('median_a'))} vs {_fmt(w.get('median_b'))}; p={_fmt(w.get('p'), sci=True)} |"
        )
    text = f"""# FINDING — triple-merge Milo vs malignant CLDN4

**Additive CLDN4-only.** Neighbourhood abundance versus **sample-level malignant CLDN4** on the public processed files from GSE131907 (Kim 2020; tLung epithelium+immune), GSE148071 (Wu 2021; 42 advanced-NSCLC biopsies) and GSE205335 (Ahn/Lee; non-normal ICI biopsies, patient-pooled). **No dual-high TACSTD2+CLDN4 gate. No GSE207422.**

The independent unit is the **patient / sample**, not the cell and not the overlapping neighbourhood. Do not write n = {n_cells:,}.

**Not miloR.** Graph + SpatialFDR follow Dann et al. 2022; DA is a sample-level Spearman (primary) or Welch t-test (median split). Joint Harmony of the three matrices was **not** run (15 GB RAM; ~370k cells). Each dataset has its own kNN graph with per-sample PCA centering.

## Verdict

| Dataset | Honest n (scored / in graph) | testable / nhoods | P<0.05 | min SpatialFDR | SpatialFDR<0.1 | sample CLDN4–T/NK ρ (p) |
|---|---|---|---:|---|---:|---|
{chr(10).join(rows)}
| **stacked patients** | **{n_scored}** | — | — | — | — | {_fmt(stacked.get('rho'))} ({_fmt(stacked.get('p'), sci=True)}) |

{floor_txt}

A SpatialFDR<0.1 hit that is |ρ|=1 on the 5-sample testability floor is **not** a cohort DA claim. At n_present ≥ 8 the three graphs are the numbers to quote.

## Honest n

| item | n | note |
|---|---:|---|
| Datasets | **3** | GSE131907 tLung + GSE148071 + GSE205335; GSE207422 excluded |
| Samples in the three graphs | **{overall['honest_n']['n_samples_in_graphs']}** | the DA unit |
| Samples with malignant CLDN4 (≥10 malignant cells) | **{n_scored}** | primary Spearman *n* |
| Cells in the three graphs | {n_cells:,} | **do not cite as n** |
| Neighbourhoods (stacked) | {overall['honest_n']['n_nhoods']} | overlapping; not independent |
| Testable nhoods (stacked, still per-graph FDR) | {overall['honest_n']['n_testable']} | ≥5 samples present **inside that dataset** |

GSE131907 uses **tLung only** so brain mets / PE / mLN are not mixed into the lung-tumor graph. GSE205335 drops normal LN/brain/lung and pools biopsies by patient. GSE148071 is one biopsy per patient; malignant-like is a marker proxy (author CNA IDs are not public).

## SpatialFDR sensitivity (min samples present)

| Dataset | floor | testable | P<0.05 | min SpatialFDR | SpatialFDR<0.1 | \\|ρ\\|=1 |
|---|---|---:|---:|---|---:|---:|
{chr(10).join(sens_lines)}

## Composition (transcriptional kNN, not histology)

Epithelium clusters with epithelium. Unrestricted nhood CLDN4 vs T/NK is lineage geometry. The sample-paired Wilcoxon is the same geometry at patient resolution. Disjoint-interface and sample-level Spearman are the honest composition tests.

| Dataset | interface n | iface ρ (p) | disjoint iface | disjoint ρ (p) | paired T/NK high vs low |
|---|---:|---|---|---|---|
{chr(10).join(iface_lines)}

Stacked sample-level malignant CLDN4 vs T/NK: n={stacked.get('n')}, ρ={_fmt(stacked.get('rho'))}, p={_fmt(stacked.get('p'), sci=True)}. CLDN4 scales are not Harmony-aligned across studies; this is a mixed-study rank correlation of independently computed scores.

## What this does not say

1. It does not invent an ICI / MPR contrast on GSE131907 or GSE148071 (GEO has none). GSE205335 RECIST is not re-tested here.
2. It does not treat {n_cells:,} cells as *n*.
3. Neighbourhood “next to” is kNN co-membership, not histology.
4. miloR / edgeR QLF numbers are not claimed.
5. Dual-high TACSTD2+CLDN4 was not used. TACSTD2 is not a gate.
6. GSE207422 was not downloaded or analysed.
7. Joint Harmony was not feasible from the public processed files at this memory budget; per-dataset graphs are the designed fallback.
8. mRNA ≠ protein.

This is extra weight on whether **neighbourhood DA vs malignant CLDN4** appears when the three public lung scRNA matrices are analysed additively. It is not a GSE207422 MPR replicate and it is not a spatial-niche paper.

## Files

- `tables/nhoods.tsv` (stacked; dataset column)
- `tables/sample_scores.tsv`, `tables/honest_n.tsv`, `tables/one_row.tsv`, `tables/summary.json`
- per-dataset folders under `tables/` and `figures/`
- extra figures: `figures/fig_triple_volcano.png`, `fig_patient_cldn4_vs_tnk.png`, `fig_honest_n.png`, `fig_testable_nhoods.png`, `fig_min_spatialfdr.png`, `fig_interface_cldn4_vs_tnk.png`, `fig_sample_paired_tnk.png`, `fig_rho_hist_testable.png`

```bash
python3 methods/triple_scrna_milo_cldn4/scripts/download.py
python3 methods/triple_scrna_milo_cldn4/scripts/analyze.py
```
"""
    out_path.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/triple_scrna_milo_cldn4"))
    ap.add_argument("--finding", type=Path, default=Path("methods/triple_scrna_milo_cldn4/FINDING.md"))
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--d", type=int, default=30)
    ap.add_argument("--prop", type=float, default=0.1)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-interface", type=int, default=3)
    ap.add_argument("--min-malignant", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["GSE131907", "GSE148071", "GSE205335"],
    )
    args = ap.parse_args()
    tables = args.outdir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    (args.outdir / "figures").mkdir(parents=True, exist_ok=True)
    self_test()

    loaders = {
        "GSE131907": lambda: load_gse131907(args.datadir / "GSE131907", args.n_hvg),
        "GSE148071": lambda: load_gse148071(args.datadir / "GSE148071", args.n_hvg),
        "GSE205335": lambda: load_gse205335(args.datadir / "GSE205335", args.n_hvg),
    }
    results: dict[str, dict] = {}
    nhood_frames = []
    sample_frames = []
    for ds in args.datasets:
        print(f"==== {ds} ====", flush=True)
        pack = loaders[ds]()
        rec = run_one_graph(pack, tables, args)
        results[ds] = rec
        nhood_frames.append(rec["nhoods"])
        sample_frames.append(rec["sample_tab"])
        # free the HVG matrix
        del pack

    stacked_nhoods = pd.concat(nhood_frames, ignore_index=True)
    stacked_nhoods.to_csv(tables / "nhoods.tsv", sep="\t", index=False)
    patients = pd.concat(sample_frames, ignore_index=True)
    patients.to_csv(tables / "sample_scores.tsv", sep="\t", index=False)
    stacked_sp = spearman_safe(patients["malignant_cldn4"], patients["frac_tnk"])

    extra_figures(results, args.outdir, stacked_nhoods, patients)

    per = {ds: results[ds]["summary"] for ds in args.datasets}
    n_scored = int(patients["malignant_cldn4"].notna().sum())
    n_samples = int(len(patients))
    n_cells = int(patients["n_cells"].sum())
    n_nhoods = int(len(stacked_nhoods))
    n_testable = int(stacked_nhoods["cldn4_testable"].fillna(False).sum()) if "cldn4_testable" in stacked_nhoods else 0
    overall = {
        "merge": "GSE131907+GSE148071+GSE205335",
        "gene": "CLDN4",
        "dual_high": False,
        "gse207422": False,
        "miloR": False,
        "harmony": False,
        "harmony_reason": "joint Harmony of the three public processed matrices is not feasible at 15 GB RAM; graphs are per dataset with per-sample PCA centering",
        "da_model": "sample-level Spearman of nhood proportion vs malignant CLDN4; SpatialFDR = miloR graphSpatialFDR k-distance",
        "independent_unit": "patient/sample, not cell",
        "per_dataset": per,
        "stacked_patient_CLDN4_vs_TNK": stacked_sp,
        "honest_n": {
            "n_datasets": len(args.datasets),
            "n_samples_in_graphs": n_samples,
            "n_scored_samples": n_scored,
            "n_nhoods": n_nhoods,
            "n_testable": n_testable,
            "do_not_cite_n_cells": n_cells,
        },
    }
    write_json(tables / "summary.json", overall)

    honest = pd.DataFrame(
        [
            {"item": "datasets", "n": len(args.datasets), "note": "GSE131907 tLung + GSE148071 + GSE205335; no GSE207422"},
            {"item": "samples in graphs", "n": n_samples, "note": "DA unit"},
            {"item": "samples with malignant CLDN4 (>=10 mal. cells)", "n": n_scored, "note": "primary Spearman n"},
            {"item": "cells in graphs", "n": n_cells, "note": "do not cite as n"},
            {"item": "neighbourhoods stacked", "n": n_nhoods, "note": "overlapping; FDR is per-graph"},
            {"item": "testable nhoods stacked", "n": n_testable, "note": ">=5 samples present inside that dataset"},
        ]
    )
    for ds in args.datasets:
        s = per[ds]["samples"]
        d = per[ds]["da_malignant_cldn4"]
        honest = pd.concat(
            [
                honest,
                pd.DataFrame(
                    [
                        {
                            "item": f"{ds} scored samples",
                            "n": s["n_with_malignant_ge10"],
                            "note": f"{s['n']} in graph; SpatialFDR<0.1 = {d['n_SpatialFDR_lt_0.1']}",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
    honest.to_csv(tables / "honest_n.tsv", sep="\t", index=False)

    one = pd.DataFrame(
        [
            {
                "merge": "GSE131907+GSE148071+GSE205335",
                "unit": "patient/sample",
                "n_scored": n_scored,
                "n_samples": n_samples,
                "n_nhoods": n_nhoods,
                "n_testable": n_testable,
                "stacked_CLDN4_TNK_rho": stacked_sp.get("rho"),
                "stacked_CLDN4_TNK_p": stacked_sp.get("p"),
                "harmony": False,
                "miloR": False,
                "dual_high": False,
                "gse207422": False,
            }
        ]
    )
    one.to_csv(tables / "one_row.tsv", sep="\t", index=False)
    write_finding(overall, args.finding)
    print(
        json.dumps(
            {
                "honest_n": overall["honest_n"],
                "stacked": stacked_sp,
                "per_dataset_da": {ds: per[ds]["da_malignant_cldn4"] for ds in args.datasets},
                "finding": str(args.finding),
                "nhoods": str(tables / "nhoods.tsv"),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
