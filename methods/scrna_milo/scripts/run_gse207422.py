#!/usr/bin/env python3
"""Milo-style kNN neighbourhood DA on GSE207422 processed UMI.

Public only: GEO supplementary UMI matrix + sample xlsx. Author barcode
labels are not deposited; lineages are reconstructed from Hu et al. 2023
canonical markers (same scheme as prior A3 work on this repo).
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


def wanted_markers() -> set[str]:
    genes = set(ALWAYS)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    genes.update(NORMAL_LUNG)
    return genes


def stream_pass1(matrix_path: Path, markers: set[str]):
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        names: list[str] = []
        means: list[float] = []
        varis: list[float] = []
        found: dict[str, np.ndarray] = {}
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} values, expected {n}")
            n_umi += arr
            names.append(gene)
            mu = float(arr.mean())
            means.append(mu)
            varis.append(float(arr.var()) if n else 0.0)
            if gene in markers:
                found[gene] = arr
            if i % 4000 == 0:
                print(f"  pass1 {i} genes, markers={len(found)}", flush=True)
    return cell_ids, n_umi, names, np.asarray(means), np.asarray(varis), found


def stream_pass2(matrix_path: Path, keep: set[str], n_cells: int) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        handle.readline()
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if gene not in keep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_cells}")
            out[gene] = arr
            if len(out) % 200 == 0:
                print(f"  pass2 stored {len(out)}/{len(keep)}", flush=True)
    missing = keep - set(out)
    if missing:
        print(f"  missing {len(missing)} requested genes (ok if absent)", flush=True)
    return out


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


def plot_volcano(df: pd.DataFrame, x: str, title: str, path: Path, xlabel: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    testable = df["testable"].to_numpy()
    ax.scatter(df.loc[~testable, x], -np.log10(np.clip(df.loc[~testable, "p"], EPS, 1)), s=6, c="#bbbbbb", alpha=0.4, label="not testable")
    t = df[testable]
    sig = t["SpatialFDR"] < 0.1
    ax.scatter(t.loc[~sig, x], -np.log10(np.clip(t.loc[~sig, "p"], EPS, 1)), s=8, c="#4c72b0", alpha=0.55, label="testable SpatialFDR≥0.1")
    ax.scatter(t.loc[sig, x], -np.log10(np.clip(t.loc[sig, "p"], EPS, 1)), s=14, c="#c44e52", alpha=0.85, label="SpatialFDR<0.1")
    ax.axhline(-np.log10(0.05), ls="--", c="0.5", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"−log10 p (sample-level test)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_scatter(x, y, path: Path, xlabel: str, ylabel: str, title: str, hue=None) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    if hue is None:
        ax.scatter(x, y, s=8, c="#4c72b0", alpha=0.45)
    else:
        ax.scatter(x, y, s=8, c=hue, cmap="viridis", alpha=0.55)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"))
    ap.add_argument("--metadata", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/scrna_milo/results/GSE207422"))
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--d", type=int, default=30)
    ap.add_argument("--prop", type=float, default=0.1)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-interface", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    self_test()

    meta = pd.read_excel(args.metadata).dropna(subset=["Sample", "Patient"]).copy()
    meta["response_paper"] = meta["Pathologic Response"].replace({"pCR": "MPR"})
    meta["is_post"] = meta["Resource"].astype(str).str.contains("Post", case=False)
    meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    markers = wanted_markers()
    print("pass1: gene stats + markers", flush=True)
    cell_ids, n_umi, gene_names, gene_mean, gene_var, marker_expr = stream_pass1(args.matrix, markers)
    n = len(cell_ids)
    print(f"cells={n} genes={len(gene_names)} median_nUMI={np.median(n_umi):.0f}", flush=True)

    hvg_idx = select_hvg(gene_mean, gene_var, n_hvg=args.n_hvg)
    hvg_names = [gene_names[i] for i in hvg_idx]
    keep = set(hvg_names) | markers
    print(f"pass2: {len(hvg_names)} HVG + {len(markers)} markers", flush=True)
    expr = stream_pass2(args.matrix, keep, n)
    # reuse pass1 marker arrays if pass2 missed nothing
    for g, arr in marker_expr.items():
        expr.setdefault(g, arr)

    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
    lineage = assign_lineage(log_cp, n)
    alveolar = np.maximum.reduce([gene_log(log_cp, g, n) for g in ["SFTPA2", "SFTPA1", "SFTPC", "SFTPB", "AGER"]])
    club = np.maximum(gene_log(log_cp, "SCGB1A1", n), gene_log(log_cp, "SCGB3A2", n))
    ciliated = np.maximum.reduce([gene_log(log_cp, g, n) for g in ["TPPP3", "FOXJ1", "CAPS"]])
    clear_normal = (alveolar >= 1.0) | (club >= 1.0) | (ciliated >= 1.0)
    is_epi = lineage == "Epithelial"
    is_malig = is_epi & (~clear_normal)
    is_tnk = np.isin(lineage, ["T", "NK"])
    tacstd2 = log_cp.get("TACSTD2", np.zeros(n, dtype=np.float32))

    sample = np.array([c.rsplit("_", 1)[0] for c in cell_ids])
    smap = meta.set_index("Sample")
    response = np.array([smap.loc[s, "response_paper"] if s in smap.index else "NA" for s in sample], dtype=object)
    is_post = np.array([bool(smap.loc[s, "is_post"]) if s in smap.index else False for s in sample])
    patient = np.array([smap.loc[s, "Patient"] if s in smap.index else s for s in sample], dtype=object)

    # Primary Milo graph: post-treatment surgery cells only (paper Fig. 1).
    keep_cells = is_post
    print(f"post-treatment cells: {keep_cells.sum()} / {n}", flush=True)
    idx_map = np.flatnonzero(keep_cells)
    sample_k = sample[keep_cells]
    lineage_k = lineage[keep_cells]
    is_malig_k = is_malig[keep_cells]
    is_tnk_k = is_tnk[keep_cells]
    tac_k = tacstd2[keep_cells]
    resp_k = response[keep_cells]
    n_k = int(keep_cells.sum())

    hvg_present = [g for g in hvg_names if g in expr]
    X = np.vstack([np.log1p(expr[g][keep_cells] * scale[keep_cells]) for g in hvg_present]).T.astype(np.float32)
    print(f"HVG matrix {X.shape}; PCA+kNN k={args.k} d={args.d}", flush=True)
    pcs, knn_idx, knn_dist = pca_knn(X, sample_k, n_pcs=args.d, k=args.k, random_state=args.seed)
    indices = refine_indices(pcs, knn_idx, prop=args.prop, random_state=args.seed)
    nhoods = make_nhoods(knn_idx, knn_dist, indices)
    members = nhood_membership(nhoods, n_k)
    print(f"neighbourhoods={len(members)} median_size={np.median([m.size for m in members]):.0f}", flush=True)

    sample_levels = sorted(pd.unique(sample_k))
    counts = count_matrix(members, sample_k, sample_levels)
    sizes = sample_sizes(sample_k, sample_levels)
    group = np.array([str(smap.loc[s, "response_paper"]) for s in sample_levels], dtype=object)
    # sample-level malignant TACSTD2 (post samples; NaN if <10 malignant cells)
    mal_score = []
    tnk_frac_sample = []
    n_mal_s = []
    for s in sample_levels:
        m = sample_k == s
        n_mal = int((m & is_malig_k).sum())
        n_mal_s.append(n_mal)
        mal_score.append(float(tac_k[m & is_malig_k].mean()) if n_mal >= 10 else np.nan)
        tnk_frac_sample.append(float(is_tnk_k[m].mean()))
    mal_score = np.asarray(mal_score, dtype=float)
    tnk_frac_sample = np.asarray(tnk_frac_sample, dtype=float)

    da_mpr = welch_da(counts, sizes, group, group_a="MPR", group_b="NMPR", min_per_group=2, min_samples=4)
    da_mpr = attach_fdr(da_mpr, nhoods.k_distance)
    da_tac = spearman_da(counts, sizes, mal_score, min_samples=5)
    da_tac = attach_fdr(da_tac, nhoods.k_distance)

    comp = nhood_composition(members, is_tnk_k, is_malig_k, tac_k, lineage_k)
    merged = comp.merge(da_mpr, on="nhood", suffixes=("", "_mpr")).merge(
        da_tac[["nhood", "spearman_rho", "p", "testable", "BH_FDR", "SpatialFDR"]].rename(
            columns={
                "spearman_rho": "tac_rho",
                "p": "tac_p",
                "testable": "tac_testable",
                "BH_FDR": "tac_BH_FDR",
                "SpatialFDR": "tac_SpatialFDR",
            }
        ),
        on="nhood",
    )
    merged.to_csv(args.outdir / "nhoods.tsv", sep="\t", index=False)
    da_mpr.to_csv(args.outdir / "da_mpr_vs_nmpr.tsv", sep="\t", index=False)
    da_tac.to_csv(args.outdir / "da_malignant_tacstd2.tsv", sep="\t", index=False)

    # Composition tests
    iface = (comp["n_malig"] >= args.min_interface) & (comp["n_tnk"] >= args.min_interface)
    all_sp = spearman_safe(comp["tacstd2_malig_mean"], comp["frac_tnk"])
    iface_sp = spearman_safe(comp.loc[iface, "tacstd2_malig_mean"], comp.loc[iface, "frac_tnk"])
    indep = greedy_independent(members, max_shared=0)
    indep_mask = np.zeros(len(members), dtype=bool)
    indep_mask[indep] = True
    indep_sp = spearman_safe(
        comp.loc[indep_mask & iface, "tacstd2_malig_mean"],
        comp.loc[indep_mask & iface, "frac_tnk"],
    )
    # TACSTD2-high / low among nhoods with ≥5 malignant cells
    has_mal = comp["n_malig"] >= 5
    med = float(np.nanmedian(comp.loc[has_mal, "tacstd2_malig_mean"])) if has_mal.any() else np.nan
    high = has_mal.to_numpy() & (comp["tacstd2_malig_mean"].to_numpy() >= med)
    low = has_mal.to_numpy() & (comp["tacstd2_malig_mean"].to_numpy() < med)
    paired = sample_paired_tnk_by_tacstd2(members, sample_k, is_tnk_k, comp["tacstd2_malig_mean"].to_numpy(), high, low)
    paired = paired.merge(
        pd.DataFrame(
            {
                "sample": sample_levels,
                "response": group,
                "malignant_tacstd2": mal_score,
                "sample_frac_tnk": tnk_frac_sample,
                "n_malignant": n_mal_s,
                "n_cells": sizes.astype(int),
            }
        ),
        on="sample",
    )
    paired.to_csv(args.outdir / "sample_paired_tnk_by_nhood_tacstd2.tsv", sep="\t", index=False)
    paired_test = wilcoxon_paired(paired["frac_tnk_high_nhood"], paired["frac_tnk_low_nhood"])

    # Among DA-vs-TACSTD2 nhoods: do T/NK-rich nhoods shrink as TACSTD2 rises?
    tac_ok = merged["tac_testable"].fillna(False).to_numpy()
    da_vs_comp = spearman_safe(merged.loc[tac_ok, "tac_rho"], merged.loc[tac_ok, "frac_tnk"])

    sample_tab = pd.DataFrame(
        {
            "sample": sample_levels,
            "patient": [str(smap.loc[s, "Patient"]) for s in sample_levels],
            "response": group,
            "n_cells": sizes.astype(int),
            "n_malignant": n_mal_s,
            "malignant_tacstd2": mal_score,
            "frac_tnk": tnk_frac_sample,
        }
    )
    sample_tab.to_csv(args.outdir / "sample_scores.tsv", sep="\t", index=False)

    def fdr_counts(df: pd.DataFrame, fdr_col: str, p_col: str = "p") -> dict:
        t = df[df["testable"]]
        return {
            "n_nhoods_total": int(len(df)),
            "n_testable": int(len(t)),
            "n_p_lt_0.05": int((t[p_col] < 0.05).sum()),
            "n_BH_FDR_lt_0.1": int((t["BH_FDR"] < 0.1).sum()),
            "n_BH_FDR_lt_0.05": int((t["BH_FDR"] < 0.05).sum()),
            "n_SpatialFDR_lt_0.1": int((t[fdr_col] < 0.1).sum()) if fdr_col in t else int((t["SpatialFDR"] < 0.1).sum()),
            "n_SpatialFDR_lt_0.05": int((t["SpatialFDR"] < 0.05).sum()),
        }

    summary = {
        "dataset": "GSE207422",
        "citation": "Hu et al. Genome Medicine 2023 PMID 36869384",
        "public_only": True,
        "author_barcode_labels": False,
        "lineage_source": "marker_argmax_Hu_canonical",
        "malignant_definition": "epithelial AND NOT (alveolar/club/ciliated log1p-CP10k >= 1); not CopyKAT",
        "miloR": False,
        "da_model": "sample-level Welch t-test (MPR vs NMPR) or Spearman (malignant TACSTD2); not edgeR QLF",
        "spatial_fdr": "miloR graphSpatialFDR k-distance reimplementation (Dann 2022 / cydar)",
        "graph": {
            "cells_all": int(n),
            "cells_post": int(n_k),
            "k": args.k,
            "d": min(args.d, X.shape[1], n_k - 1),
            "prop": args.prop,
            "n_hvg": len(hvg_present),
            "n_nhoods": int(len(members)),
            "nhood_size": args.k + 1,
            "batch": "PCA per-sample mean centering (not Harmony)",
        },
        "samples_post": {
            "n": int(len(sample_levels)),
            "MPR": int((group == "MPR").sum()),
            "NMPR": int((group == "NMPR").sum()),
            "n_with_malignant_ge10": int(np.isfinite(mal_score).sum()),
        },
        "lineage_counts_post": {k: int(v) for k, v in pd.Series(lineage_k).value_counts().items()},
        "n_malignant_post": int(is_malig_k.sum()),
        "n_tnk_post": int(is_tnk_k.sum()),
        "da_mpr_vs_nmpr": fdr_counts(da_mpr, "SpatialFDR"),
        "da_malignant_tacstd2": fdr_counts(da_tac, "SpatialFDR"),
        "composition": {
            "note": "transcriptional kNN != spatial niche; unrestricted TACSTD2 vs T/NK is partly lineage geometry",
            "all_nhoods_maligTACSTD2_vs_fracTNK": all_sp,
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
            "sample_paired_TNK_in_TACSTD2high_vs_low_nhoods": {
                "high_low_cut": "median malignant TACSTD2 among nhoods with >=5 malignant cells",
                "median_cut": med,
                "n_high_nhoods": int(high.sum()),
                "n_low_nhoods": int(low.sum()),
                "wilcoxon_high_vs_low": paired_test,
                "direction_expected": "frac_tnk_high < frac_tnk_low (immune-poor TROP2-high nhoods)",
            },
            "nhood_tacstd2_DA_rho_vs_fracTNK": da_vs_comp,
        },
        "honest_n": {
            "independent_unit_DA": "post-treatment samples (n=12; MPR=4 including pCR P06, NMPR=8)",
            "nhoods_are_overlapping": True,
            "SpatialFDR_is_overlap_aware": True,
            "disjoint_nhood_subset_n": int(indep.size),
            "do_not_cite_n_cells_as_n": int(n_k),
        },
    }
    write_json(args.outdir / "summary.json", summary)

    plot_volcano(
        da_mpr,
        "logFC_B_minus_A",
        "GSE207422 neighbourhood DA: NMPR vs MPR",
        args.outdir / "fig_da_mpr_volcano.png",
        "log2 FC (NMPR − MPR sample proportion)",
    )
    plot_volcano(
        da_tac.rename(columns={"spearman_rho": "logFC_B_minus_A", "p": "p", "testable": "testable"}),
        "logFC_B_minus_A",
        "GSE207422 neighbourhood DA vs malignant TACSTD2",
        args.outdir / "fig_da_tacstd2_volcano.png",
        "Spearman ρ (nhood abundance vs sample malignant TACSTD2)",
    )
    plot_scatter(
        comp.loc[iface, "tacstd2_malig_mean"],
        comp.loc[iface, "frac_tnk"],
        args.outdir / "fig_interface_tacstd2_vs_tnk.png",
        "Neighbourhood malignant TACSTD2 (mean log1p CP10k)",
        "Neighbourhood T/NK fraction",
        f"Interface nhoods (n={int(iface.sum())}); transcriptional, not spatial",
    )
    if paired["frac_tnk_high_nhood"].notna().any():
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        for _, r in paired.iterrows():
            if np.isfinite(r["frac_tnk_high_nhood"]) and np.isfinite(r["frac_tnk_low_nhood"]):
                c = "#c44e52" if r["response"] == "NMPR" else "#4c72b0"
                ax.plot([0, 1], [r["frac_tnk_low_nhood"], r["frac_tnk_high_nhood"]], "-o", c=c, alpha=0.75, ms=5)
        ax.set_xticks([0, 1], ["TACSTD2-low nhoods", "TACSTD2-high nhoods"])
        ax.set_ylabel("Sample T/NK fraction among cells in those nhoods")
        ax.set_title("GSE207422 sample-paired (unit = patient)")
        fig.tight_layout()
        fig.savefig(args.outdir / "fig_sample_paired_tnk.png", dpi=160)
        plt.close(fig)

    print(json_preview(summary), flush=True)


def json_preview(d: dict) -> str:
    import json

    return json.dumps(
        {
            "n_nhoods": d["graph"]["n_nhoods"],
            "da_mpr": d["da_mpr_vs_nmpr"],
            "da_tac": d["da_malignant_tacstd2"],
            "iface": d["composition"]["interface_nhoods"],
            "paired": d["composition"]["sample_paired_TNK_in_TACSTD2high_vs_low_nhoods"],
        },
        indent=2,
    )


if __name__ == "__main__":
    main()
