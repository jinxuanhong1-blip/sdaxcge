#!/usr/bin/env python3
"""Milo-style kNN neighbourhood DA on GSE241934 public processed MTX.

Cohorts
-------
- IIT (NEOTIDE/CTONG2104, 11 EGFR-mutant resected tumours): primary add-on.
- RWC (34 real-world neoadjuvant IO tumours): run if the MTX is present and
  memory allows. pCR is grouped with MPR.

Author major_cell_type is public (unlike GSE207422). Epithelial cells are
used as the malignant-like compartment; this is not a CNV call.
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
from scipy.io import mmread

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


def load_mtx(matrix: Path, barcodes: Path, features: Path, meta: Path) -> tuple[object, pd.DataFrame, np.ndarray]:
    print(f"reading {matrix}", flush=True)
    X = mmread(gzip.open(matrix, "rb")).tocsr()
    bc = pd.read_csv(barcodes, sep="\t", header=None)[0].astype(str)
    feat = pd.read_csv(features, sep="\t", header=None)
    genes = feat[0].astype(str)
    if X.shape[0] == len(bc) and X.shape[1] == len(genes):
        X = X.T.tocsr()
    if X.shape != (len(genes), len(bc)):
        raise ValueError(f"MTX {X.shape} vs genes={len(genes)} barcodes={len(bc)}")
    obs = pd.read_csv(meta, sep="\t")
    # align meta to MTX barcodes
    if "cellID" in obs.columns:
        obs = obs.set_index("cellID")
        obs = obs.reindex(bc)
    else:
        obs.index = bc
    if obs["sampleID"].isna().mean() > 0.05:
        raise RuntimeError("barcode/meta alignment failed")
    return X, obs, genes.to_numpy()


def gene_index(genes: np.ndarray, name: str) -> int | None:
    hits = np.flatnonzero(genes == name)
    return int(hits[0]) if hits.size else None


def extract_log1p_cp10k(X, genes: np.ndarray, n_umi: np.ndarray, name: str) -> np.ndarray:
    j = gene_index(genes, name)
    if j is None:
        return np.zeros(X.shape[1], dtype=np.float32)
    raw = np.asarray(X[j].todense()).ravel().astype(np.float32)
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    return np.log1p(raw * scale).astype(np.float32)


def plot_volcano(df: pd.DataFrame, x: str, title: str, path: Path, xlabel: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    testable = df["testable"].to_numpy()
    ax.scatter(df.loc[~testable, x], -np.log10(np.clip(df.loc[~testable, "p"], EPS, 1)), s=6, c="#bbbbbb", alpha=0.4)
    t = df[testable]
    sig = t["SpatialFDR"] < 0.1
    ax.scatter(t.loc[~sig, x], -np.log10(np.clip(t.loc[~sig, "p"], EPS, 1)), s=8, c="#4c72b0", alpha=0.55, label="SpatialFDR≥0.1")
    ax.scatter(t.loc[sig, x], -np.log10(np.clip(t.loc[sig, "p"], EPS, 1)), s=14, c="#c44e52", alpha=0.85, label="SpatialFDR<0.1")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"−log10 p")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def run_cohort(
    tag: str,
    X,
    obs: pd.DataFrame,
    genes: np.ndarray,
    outdir: Path,
    k: int,
    d: int,
    prop: float,
    n_hvg: int,
    min_interface: int,
    seed: int,
) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    n_umi = np.asarray(X.sum(axis=0)).ravel()
    # HVG on raw means/vars
    gmean = np.asarray(X.mean(axis=1)).ravel()
    # E[X^2] - mean^2
    gvar = np.asarray(X.multiply(X).mean(axis=1)).ravel() - gmean**2
    hvg_idx = select_hvg(gmean, gvar, n_hvg=n_hvg)
    print(f"{tag}: cells={X.shape[1]} HVG={len(hvg_idx)}", flush=True)
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    # dense HVG log1p CP10k (cells × HVG)
    hvg = X[hvg_idx].astype(np.float32).toarray()
    log_hvg = np.log1p(hvg * scale).T.astype(np.float32)
    del hvg

    sample = obs["sampleID"].astype(str).to_numpy()
    resp_raw = obs["Pathological Response"].astype(str).to_numpy()
    resp = np.where(np.isin(resp_raw, ["MPR", "pCR"]), "MPR", np.where(resp_raw == "non-MPR", "NMPR", resp_raw))
    lineage = obs["major_cell_type"].fillna(obs.get("major.cell.type", "NA")).astype(str).to_numpy()
    is_malig = np.isin(lineage, ["Epi", "Epithelial", "epithelial"])
    is_tnk = np.isin(lineage, ["T", "NK"])
    tac = extract_log1p_cp10k(X, genes, n_umi, "TACSTD2")

    pcs, knn_idx, knn_dist = pca_knn(log_hvg, sample, n_pcs=d, k=k, random_state=seed)
    del log_hvg
    indices = refine_indices(pcs, knn_idx, prop=prop, random_state=seed)
    nhoods = make_nhoods(knn_idx, knn_dist, indices)
    members = nhood_membership(nhoods, X.shape[1])
    print(f"{tag}: nhoods={len(members)}", flush=True)

    sample_levels = sorted(pd.unique(sample))
    counts = count_matrix(members, sample, sample_levels)
    sizes = sample_sizes(sample, sample_levels)
    group = []
    mal_score = []
    tnk_frac = []
    n_mal_s = []
    for s in sample_levels:
        m = sample == s
        group.append(str(pd.Series(resp[m]).mode().iloc[0]))
        n_mal = int((m & is_malig).sum())
        n_mal_s.append(n_mal)
        mal_score.append(float(tac[m & is_malig].mean()) if n_mal >= 10 else np.nan)
        tnk_frac.append(float(is_tnk[m].mean()))
    group = np.array(group, dtype=object)
    mal_score = np.asarray(mal_score, dtype=float)
    tnk_frac = np.asarray(tnk_frac, dtype=float)

    da_mpr = attach_fdr(welch_da(counts, sizes, group, "MPR", "NMPR", min_per_group=2, min_samples=4), nhoods.k_distance)
    da_tac = attach_fdr(spearman_da(counts, sizes, mal_score, min_samples=5), nhoods.k_distance)
    comp = nhood_composition(members, is_tnk, is_malig, tac, lineage)
    merged = comp.merge(da_mpr, on="nhood").merge(
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
    merged.to_csv(outdir / "nhoods.tsv", sep="\t", index=False)
    da_mpr.to_csv(outdir / "da_mpr_vs_nmpr.tsv", sep="\t", index=False)
    da_tac.to_csv(outdir / "da_malignant_tacstd2.tsv", sep="\t", index=False)

    iface = (comp["n_malig"] >= min_interface) & (comp["n_tnk"] >= min_interface)
    all_sp = spearman_safe(comp["tacstd2_malig_mean"], comp["frac_tnk"])
    iface_sp = spearman_safe(comp.loc[iface, "tacstd2_malig_mean"], comp.loc[iface, "frac_tnk"])
    indep = greedy_independent(members, max_shared=0)
    indep_mask = np.zeros(len(members), dtype=bool)
    indep_mask[indep] = True
    indep_sp = spearman_safe(comp.loc[indep_mask & iface, "tacstd2_malig_mean"], comp.loc[indep_mask & iface, "frac_tnk"])
    has_mal = comp["n_malig"] >= 5
    med = float(np.nanmedian(comp.loc[has_mal, "tacstd2_malig_mean"])) if has_mal.any() else np.nan
    high = has_mal.to_numpy() & (comp["tacstd2_malig_mean"].to_numpy() >= med)
    low = has_mal.to_numpy() & (comp["tacstd2_malig_mean"].to_numpy() < med)
    paired = sample_paired_tnk_by_tacstd2(members, sample, is_tnk, comp["tacstd2_malig_mean"].to_numpy(), high, low)
    paired = paired.merge(
        pd.DataFrame(
            {
                "sample": sample_levels,
                "response": group,
                "malignant_tacstd2": mal_score,
                "sample_frac_tnk": tnk_frac,
                "n_malignant": n_mal_s,
                "n_cells": sizes.astype(int),
            }
        ),
        on="sample",
    )
    paired.to_csv(outdir / "sample_paired_tnk_by_nhood_tacstd2.tsv", sep="\t", index=False)
    paired_test = wilcoxon_paired(paired["frac_tnk_high_nhood"], paired["frac_tnk_low_nhood"])
    tac_ok = merged["tac_testable"].fillna(False).to_numpy()
    da_vs_comp = spearman_safe(merged.loc[tac_ok, "tac_rho"], merged.loc[tac_ok, "frac_tnk"])

    pd.DataFrame(
        {
            "sample": sample_levels,
            "response": group,
            "n_cells": sizes.astype(int),
            "n_malignant": n_mal_s,
            "malignant_tacstd2": mal_score,
            "frac_tnk": tnk_frac,
        }
    ).to_csv(outdir / "sample_scores.tsv", sep="\t", index=False)

    def fdr_counts(df: pd.DataFrame) -> dict:
        t = df[df["testable"]]
        return {
            "n_nhoods_total": int(len(df)),
            "n_testable": int(len(t)),
            "n_p_lt_0.05": int((t["p"] < 0.05).sum()),
            "n_BH_FDR_lt_0.1": int((t["BH_FDR"] < 0.1).sum()),
            "n_BH_FDR_lt_0.05": int((t["BH_FDR"] < 0.05).sum()),
            "n_SpatialFDR_lt_0.1": int((t["SpatialFDR"] < 0.1).sum()),
            "n_SpatialFDR_lt_0.05": int((t["SpatialFDR"] < 0.05).sum()),
        }

    summary = {
        "dataset": f"GSE241934_{tag}",
        "citation": "Zhao et al. Signal Transduct Target Ther 2024; GEO GSE241934 NEOTIDE/CTONG2104",
        "public_only": True,
        "author_barcode_labels": True,
        "lineage_source": "author major_cell_type",
        "malignant_definition": "author Epi (not CopyKAT); TACSTD2 scored on Epi cells",
        "miloR": False,
        "da_model": "sample-level Welch t-test / Spearman; not edgeR QLF",
        "spatial_fdr": "miloR graphSpatialFDR k-distance reimplementation",
        "graph": {
            "cells": int(X.shape[1]),
            "k": k,
            "d": min(d, log_hvg.shape[1] if False else 30),
            "prop": prop,
            "n_hvg": int(len(hvg_idx)),
            "n_nhoods": int(len(members)),
            "nhood_size": k + 1,
        },
        "samples": {
            "n": int(len(sample_levels)),
            "MPR": int((group == "MPR").sum()),
            "NMPR": int((group == "NMPR").sum()),
            "n_with_epithelial_ge10": int(np.isfinite(mal_score).sum()),
        },
        "lineage_counts": {k_: int(v) for k_, v in pd.Series(lineage).value_counts().items()},
        "n_epi": int(is_malig.sum()),
        "n_tnk": int(is_tnk.sum()),
        "da_mpr_vs_nmpr": fdr_counts(da_mpr),
        "da_malignant_tacstd2": fdr_counts(da_tac),
        "composition": {
            "all_nhoods_epiTACSTD2_vs_fracTNK": all_sp,
            "interface_nhoods": {"n": int(iface.sum()), "spearman": iface_sp},
            "disjoint_interface": {
                "n_disjoint_nhoods": int(indep.size),
                "n_disjoint_interface": int((indep_mask & iface.to_numpy()).sum()),
                "spearman": indep_sp,
            },
            "sample_paired_TNK_in_TACSTD2high_vs_low_nhoods": {
                "median_cut": med,
                "n_high_nhoods": int(high.sum()),
                "n_low_nhoods": int(low.sum()),
                "wilcoxon_high_vs_low": paired_test,
            },
            "nhood_tacstd2_DA_rho_vs_fracTNK": da_vs_comp,
        },
        "honest_n": {
            "independent_unit_DA": f"samples (n={len(sample_levels)}; MPR={(group=='MPR').sum()}, NMPR={(group=='NMPR').sum()})",
            "nhoods_are_overlapping": True,
            "disjoint_nhood_subset_n": int(indep.size),
        },
    }
    # fix d in graph
    summary["graph"]["d"] = int(pcs.shape[1])
    write_json(outdir / "summary.json", summary)

    plot_volcano(da_mpr, "logFC_B_minus_A", f"GSE241934 {tag} DA: NMPR vs MPR", outdir / "fig_da_mpr_volcano.png", "log2 FC (NMPR − MPR)")
    plot_volcano(
        da_tac.rename(columns={"spearman_rho": "logFC_B_minus_A"}),
        "logFC_B_minus_A",
        f"GSE241934 {tag} DA vs epithelial TACSTD2",
        outdir / "fig_da_tacstd2_volcano.png",
        "Spearman ρ (abundance vs sample Epi TACSTD2)",
    )
    if iface.any():
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
        ax.scatter(comp.loc[iface, "tacstd2_malig_mean"], comp.loc[iface, "frac_tnk"], s=8, alpha=0.5, c="#4c72b0")
        ax.set_xlabel("Nhood epithelial TACSTD2 (mean log1p CP10k)")
        ax.set_ylabel("Nhood T/NK fraction")
        ax.set_title(f"{tag} interface nhoods n={int(iface.sum())}")
        fig.tight_layout()
        fig.savefig(outdir / "fig_interface_tacstd2_vs_tnk.png", dpi=160)
        plt.close(fig)
    if paired["frac_tnk_high_nhood"].notna().any():
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        for _, r in paired.iterrows():
            if np.isfinite(r["frac_tnk_high_nhood"]) and np.isfinite(r["frac_tnk_low_nhood"]):
                c = "#c44e52" if r["response"] == "NMPR" else "#4c72b0"
                ax.plot([0, 1], [r["frac_tnk_low_nhood"], r["frac_tnk_high_nhood"]], "-o", c=c, alpha=0.75, ms=5)
        ax.set_xticks([0, 1], ["TACSTD2-low nhoods", "TACSTD2-high nhoods"])
        ax.set_ylabel("Sample T/NK fraction in those nhoods")
        ax.set_title(f"{tag} sample-paired")
        fig.tight_layout()
        fig.savefig(outdir / "fig_sample_paired_tnk.png", dpi=160)
        plt.close(fig)
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data/GSE241934"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/scrna_milo/results/GSE241934"))
    ap.add_argument("--cohorts", nargs="+", default=["IIT"], choices=["IIT", "RWC"])
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--d", type=int, default=30)
    ap.add_argument("--prop", type=float, default=0.1)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-interface", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    self_test()

    specs = {
        "IIT": {
            "matrix": args.datadir / "GSE241934_IIT_Matrix.mtx.gz",
            "barcodes": args.datadir / "GSE241934_IIT_barcodes.tsv.gz",
            "features": args.datadir / "GSE241934_IIT_features.tsv.gz",
            "meta": args.datadir / "GSE241934_IIT_Meta.txt.gz",
        },
        "RWC": {
            "matrix": args.datadir / "GSE241934_Real_Matrix.mtx.gz",
            "barcodes": args.datadir / "GSE241934_RWC_barcodes.tsv.gz",
            "features": args.datadir / "GSE241934_RWC_features.tsv.gz",
            "meta": args.datadir / "GSE241934_Real_Meta.txt.gz",
        },
    }
    all_summ = {}
    for tag in args.cohorts:
        spec = specs[tag]
        missing = [p for p in spec.values() if not p.exists()]
        if missing:
            all_summ[tag] = {"skipped": True, "missing": [str(p) for p in missing]}
            write_json(args.outdir / tag / "feasibility.json", all_summ[tag])
            continue
        X, obs, genes = load_mtx(spec["matrix"], spec["barcodes"], spec["features"], spec["meta"])
        all_summ[tag] = run_cohort(
            tag, X, obs, genes, args.outdir / tag, args.k, args.d, args.prop, args.n_hvg, args.min_interface, args.seed
        )
        del X
    write_json(args.outdir / "summary.json", all_summ)
    print({k: (v.get("da_mpr_vs_nmpr") if isinstance(v, dict) else v) for k, v in all_summ.items()}, flush=True)


if __name__ == "__main__":
    main()
