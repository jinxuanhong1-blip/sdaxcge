#!/usr/bin/env python3
"""EXTRA series: GSE291670 neoadjuvant anlotinib + camrelizumab NSCLC scRNA.

Malignant-restricted TACSTD2 / CLDN4 vs MPR and vs T/NK fraction.
Public 10x MTX only. Author per-cell labels are not on GEO.

This is a second public neoadjuvant ICI lung series, not a GSE207422 audit.
GSE243013 (likely the intended 'GSE253013') is CD45+ immune-only.
GSE266035 is circulating T cells from one patient. Neither has malignant epithelium.
"""
from __future__ import annotations

import argparse
import gzip
import json
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

SAMPLES = [
    ("GSM8839599_MPR-1", "MPR", "MPR-1"),
    ("GSM8839600_MPR-2", "MPR", "MPR-2"),
    ("GSM8839601_MPR-3", "MPR", "MPR-3"),
    ("GSM8839602_Non-MPR-1", "NMPR", "Non-MPR-1"),
    ("GSM8839603_Non-MPR-2", "NMPR", "Non-MPR-2"),
    ("GSM8839604_Non-MPR-3", "NMPR", "Non-MPR-3"),
]

LINEAGE_MARKERS = {
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA"],
    "Mast": ["TPSAB1", "CPA3"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}
NORMAL_LUNG = ["SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3", "FOXJ1"]
EXTRA = ["TACSTD2", "CLDN4", "PTPRC", "CD8A", "CD4", "MKI67"]


def wanted() -> set[str]:
    g = set(EXTRA) | set(NORMAL_LUNG)
    for vs in LINEAGE_MARKERS.values():
        g.update(vs)
    return g


def read_10x(prefix: Path):
    mtx = scipy.io.mmread(str(prefix) + "_matrix.mtx.gz").tocsr().astype(np.float32)
    with gzip.open(str(prefix) + "_barcodes.tsv.gz", "rt") as fh:
        barcodes = [line.strip().split("\t")[0] for line in fh]
    genes = []
    with gzip.open(str(prefix) + "_features.tsv.gz", "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    # 10x MTX is genes x cells
    if mtx.shape[0] != len(genes) and mtx.shape[1] == len(genes):
        mtx = mtx.T.tocsr()
    if mtx.shape[0] != len(genes) or mtx.shape[1] != len(barcodes):
        raise ValueError(f"{prefix}: mtx {mtx.shape} genes {len(genes)} barcodes {len(barcodes)}")
    return mtx, np.array(genes), np.array(barcodes)


def gene_index(genes: np.ndarray) -> dict[str, int]:
    idx: dict[str, int] = {}
    for i, g in enumerate(genes):
        if g not in idx:
            idx[g] = i
    return idx


def extract_dense(mtx, gidx: dict[str, int], names: list[str]) -> dict[str, np.ndarray]:
    out = {}
    for g in names:
        if g in gidx:
            out[g] = np.asarray(mtx[gidx[g], :].todense()).ravel()
    return out


def module(log_cp: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [log_cp[g] for g in genes if g in log_cp]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(scores: dict[str, np.ndarray], cd3: np.ndarray) -> np.ndarray:
    names = list(scores)
    mat = np.vstack([scores[n] for n in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    labels = np.array(names, dtype=object)[best].copy()
    t_idx, nk_idx = names.index("T"), names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk = np.isin(labels, ["T", "NK"])
    labels[close & both & tnk & (cd3 > 0.15)] = "T"
    labels[close & both & tnk & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < 0.12] = "Unassigned"
    return labels


def exact_wilcoxon(values: np.ndarray, is_a: np.ndarray) -> dict:
    values = np.asarray(values, float)
    is_a = np.asarray(is_a, bool)
    ok = np.isfinite(values)
    values, is_a = values[ok], is_a[ok]
    n, n_a = len(values), int(is_a.sum())
    n_b = n - n_a
    if n_a < 1 or n_b < 1:
        return {"n_a": n_a, "n_b": n_b, "U": None, "p_exact": None, "mean_a": None, "mean_b": None}
    obs_u = float(stats.mannwhitneyu(values[is_a], values[~is_a], alternative="two-sided").statistic)
    expected = n_a * n_b / 2.0
    obs_ext = abs(obs_u - expected)
    count = total = 0
    mean_a = float(values[is_a].mean())
    mean_b = float(values[~is_a].mean())
    for combo in combinations(range(n), n_a):
        mask = np.zeros(n, dtype=bool)
        mask[list(combo)] = True
        u = float(stats.mannwhitneyu(values[mask], values[~mask], alternative="two-sided").statistic)
        if abs(u - expected) >= obs_ext - 1e-12:
            count += 1
        total += 1
    return {
        "n_a": n_a,
        "n_b": n_b,
        "U": obs_u,
        "p_exact": count / total,
        "mean_a": mean_a,
        "mean_b": mean_b,
        "delta_a_minus_b": mean_a - mean_b,
        "n_perm": total,
    }


def spearman_safe(x, y) -> dict:
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 4:
        return {"n": int(x.size), "rho": None, "p": None}
    rho, p = stats.spearmanr(x, y)
    return {"n": int(x.size), "rho": float(rho), "p": float(p)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data/GSE291670"))
    ap.add_argument("--outdir", type=Path, default=Path("results/extra/GSE291670"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    keep_genes = wanted()
    cell_rows = []
    sample_rows = []
    lineage_tables = []

    for prefix, response, short in SAMPLES:
        print(f"reading {prefix}", flush=True)
        mtx, genes, barcodes = read_10x(args.datadir / prefix)
        gidx = gene_index(genes)
        n_umi = np.asarray(mtx.sum(axis=0)).ravel()
        n_genes = np.asarray((mtx > 0).sum(axis=0)).ravel()
        mt_idx = [i for i, g in enumerate(genes) if str(g).startswith("MT-")]
        mt_umi = np.asarray(mtx[mt_idx, :].sum(axis=0)).ravel() if mt_idx else np.zeros(mtx.shape[1])
        pct_mt = np.where(n_umi > 0, 100.0 * mt_umi / n_umi, 0.0)
        qc = (n_genes >= 200) & (pct_mt < 20) & (n_genes < 8000)
        mtx = mtx[:, qc]
        barcodes = barcodes[qc]
        n_umi, n_genes, pct_mt = n_umi[qc], n_genes[qc], pct_mt[qc]
        n = mtx.shape[1]
        expr = extract_dense(mtx, gidx, list(keep_genes))
        scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
        log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
        scores = {name: module(log_cp, gs, n) for name, gs in LINEAGE_MARKERS.items()}
        lineage = assign_lineage(scores, log_cp.get("CD3E", log_cp.get("CD3D", np.zeros(n))))
        normal = module(log_cp, NORMAL_LUNG, n)
        tac = expr.get("TACSTD2", np.zeros(n))
        cld = expr.get("CLDN4", np.zeros(n))
        cd3e = expr.get("CD3E", np.zeros(n))
        cd8a = expr.get("CD8A", np.zeros(n))
        nkg7 = expr.get("NKG7", np.zeros(n))
        is_epi = lineage == "Epithelial"
        if is_epi.sum() >= 20:
            nl_cut = float(np.quantile(normal[is_epi], 0.75))
        else:
            nl_cut = 0.3
        is_mal = is_epi & (normal < nl_cut)
        is_tnk = np.isin(lineage, ["T", "NK"])
        is_umi_tnk = (cd3e >= 1) | (cd8a >= 1) | (nkg7 >= 1)

        df = pd.DataFrame(
            {
                "barcode": barcodes,
                "sample": prefix,
                "short": short,
                "response": response,
                "lineage": lineage,
                "n_umi": n_umi,
                "n_genes": n_genes,
                "pct_mt": pct_mt,
                "TACSTD2": tac,
                "CLDN4": cld,
                "tacstd2_log1p_cp10k": np.log1p(tac * scale),
                "cldn4_log1p_cp10k": np.log1p(cld * scale),
                "is_epithelial": is_epi.astype(int),
                "is_malignant": is_mal.astype(int),
                "is_lineage_tnk": is_tnk.astype(int),
                "is_umi_tnk": is_umi_tnk.astype(int),
            }
        )
        cell_rows.append(df)
        lineage_tables.append(df.groupby("lineage").size().rename(prefix))

        def pack(mask, gene_raw, gene_log, tag):
            sub = df.loc[mask]
            if len(sub) == 0:
                return {
                    f"{tag}_n": 0,
                    f"{tag}_mean_log1p_cp10k": np.nan,
                    f"{tag}_pct_pos": np.nan,
                }
            return {
                f"{tag}_n": int(len(sub)),
                f"{tag}_mean_log1p_cp10k": float(sub[gene_log].mean()),
                f"{tag}_pct_pos": float((sub[gene_raw] >= 1).mean()),
            }

        rec = {
            "sample": prefix,
            "short": short,
            "response": response,
            "n_cells": int(n),
            "n_epithelial": int(is_epi.sum()),
            "n_malignant": int(is_mal.sum()),
            "n_lineage_tnk": int(is_tnk.sum()),
            "n_umi_tnk": int(is_umi_tnk.sum()),
            "frac_epithelial": float(is_epi.mean()),
            "frac_malignant": float(is_mal.mean()),
            "frac_lineage_tnk": float(is_tnk.mean()),
            "frac_umi_tnk": float(is_umi_tnk.mean()),
        }
        rec.update(pack(is_mal, "TACSTD2", "tacstd2_log1p_cp10k", "mal_TACSTD2"))
        rec.update(pack(is_mal, "CLDN4", "cldn4_log1p_cp10k", "mal_CLDN4"))
        rec.update(pack(is_epi, "TACSTD2", "tacstd2_log1p_cp10k", "epi_TACSTD2"))
        rec.update(pack(is_epi, "CLDN4", "cldn4_log1p_cp10k", "epi_CLDN4"))
        sample_rows.append(rec)
        print(
            f"  qc_cells={n} epi={is_epi.sum()} mal={is_mal.sum()} "
            f"TNK={is_tnk.sum()} umiTNK={is_umi_tnk.sum()}",
            flush=True,
        )

    cells = pd.concat(cell_rows, ignore_index=True)
    per = pd.DataFrame(sample_rows)
    lin = pd.concat(lineage_tables, axis=1).fillna(0).astype(int)
    cells.to_csv(args.outdir / "cell_calls.tsv.gz", sep="\t", index=False, compression="gzip")
    per.to_csv(args.outdir / "per_sample.tsv", sep="\t", index=False)
    lin.to_csv(args.outdir / "lineage_counts.tsv", sep="\t")

    tests = []
    is_nmpr = per["response"].to_numpy() == "NMPR"

    def add_mw(col, gene, compartment):
        w = exact_wilcoxon(per[col].to_numpy(), is_nmpr)
        tests.append(
            {
                "family": "NMPR_vs_MPR",
                "gene": gene,
                "compartment": compartment,
                "metric": col,
                "n": int(per[col].notna().sum()),
                "n_NMPR": int(is_nmpr.sum()),
                "n_MPR": int((~is_nmpr).sum()),
                "mean_NMPR": w["mean_a"],
                "mean_MPR": w["mean_b"],
                "delta_NMPR_minus_MPR": w["delta_a_minus_b"],
                "U": w["U"],
                "p": w["p_exact"],
                "rho": None,
                "tnk_def": "",
                "note": f"exact C({w['n_a']+w['n_b']},{w['n_a']})={w['n_perm']}; two-sided min p at 3v3 is 0.10",
            }
        )

    def add_rho(ycol, tnkcol, gene, compartment, tnk_name):
        s = spearman_safe(per[ycol], per[tnkcol])
        tests.append(
            {
                "family": "vs_TNK",
                "gene": gene,
                "compartment": compartment,
                "metric": ycol,
                "n": s["n"],
                "n_NMPR": int(is_nmpr.sum()),
                "n_MPR": int((~is_nmpr).sum()),
                "mean_NMPR": None,
                "mean_MPR": None,
                "delta_NMPR_minus_MPR": None,
                "U": None,
                "p": s["p"],
                "rho": s["rho"],
                "tnk_def": tnk_name,
                "note": "patient-level Spearman; n=6 is underpowered for |ρ|≈0.45",
            }
        )

    for gene, mean_col, pct_col, comp in [
        ("TACSTD2", "mal_TACSTD2_mean_log1p_cp10k", "mal_TACSTD2_pct_pos", "malignant"),
        ("CLDN4", "mal_CLDN4_mean_log1p_cp10k", "mal_CLDN4_pct_pos", "malignant"),
        ("TACSTD2", "epi_TACSTD2_mean_log1p_cp10k", "epi_TACSTD2_pct_pos", "all_epithelial"),
        ("CLDN4", "epi_CLDN4_mean_log1p_cp10k", "epi_CLDN4_pct_pos", "all_epithelial"),
    ]:
        add_mw(mean_col, gene, comp)
        add_mw(pct_col, gene, comp)
        add_rho(mean_col, "frac_lineage_tnk", gene, comp, "lineage_TNK")
        add_rho(mean_col, "frac_umi_tnk", gene, comp, "umi_CD3E_CD8A_NKG7")
        add_rho(pct_col, "frac_lineage_tnk", gene, comp, "lineage_TNK")

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(args.outdir / "stats.tsv", sep="\t", index=False)

    # EXTRA figure
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.4))
    colors = {"MPR": "#2ca02c", "NMPR": "#d62728"}
    ax = axes[0, 0]
    for resp, sub in per.groupby("response"):
        ax.scatter(np.full(len(sub), 0 if resp == "MPR" else 1), sub["mal_TACSTD2_mean_log1p_cp10k"],
                   c=colors[resp], s=70, zorder=3)
        for _, r in sub.iterrows():
            ax.annotate(r["short"], (0 if resp == "MPR" else 1, r["mal_TACSTD2_mean_log1p_cp10k"]),
                        textcoords="offset points", xytext=(6, 0), fontsize=8)
    ax.set_xticks([0, 1], ["MPR", "NMPR"])
    ax.set_ylabel("malignant TACSTD2  mean log1p(CP10k)")
    row = tests_df[(tests_df.family == "NMPR_vs_MPR") & (tests_df.metric == "mal_TACSTD2_mean_log1p_cp10k")].iloc[0]
    ax.set_title(f"n=3 vs 3   p={row['p']:.2f}   Δ={row['delta_NMPR_minus_MPR']:+.3f}")

    ax = axes[0, 1]
    for resp, sub in per.groupby("response"):
        ax.scatter(np.full(len(sub), 0 if resp == "MPR" else 1), sub["mal_CLDN4_mean_log1p_cp10k"],
                   c=colors[resp], s=70, zorder=3)
        for _, r in sub.iterrows():
            ax.annotate(r["short"], (0 if resp == "MPR" else 1, r["mal_CLDN4_mean_log1p_cp10k"]),
                        textcoords="offset points", xytext=(6, 0), fontsize=8)
    ax.set_xticks([0, 1], ["MPR", "NMPR"])
    ax.set_ylabel("malignant CLDN4  mean log1p(CP10k)")
    row = tests_df[(tests_df.family == "NMPR_vs_MPR") & (tests_df.metric == "mal_CLDN4_mean_log1p_cp10k")].iloc[0]
    ax.set_title(f"n=3 vs 3   p={row['p']:.2f}   Δ={row['delta_NMPR_minus_MPR']:+.3f}")

    ax = axes[1, 0]
    for resp, sub in per.groupby("response"):
        ax.scatter(sub["frac_lineage_tnk"], sub["mal_TACSTD2_mean_log1p_cp10k"],
                   c=colors[resp], s=70, label=resp, zorder=3)
        for _, r in sub.iterrows():
            ax.annotate(r["short"], (r["frac_lineage_tnk"], r["mal_TACSTD2_mean_log1p_cp10k"]),
                        textcoords="offset points", xytext=(5, 2), fontsize=8)
    s = tests_df[(tests_df.family == "vs_TNK") & (tests_df.metric == "mal_TACSTD2_mean_log1p_cp10k")
                 & (tests_df.tnk_def == "lineage_TNK")].iloc[0]
    ax.set_xlabel("T/NK fraction (lineage)")
    ax.set_ylabel("malignant TACSTD2")
    ax.legend(fontsize=8)
    ax.set_title(f"n=6   ρ={s['rho']:.2f}   p={s['p']:.2f}")

    ax = axes[1, 1]
    for resp, sub in per.groupby("response"):
        ax.scatter(sub["frac_lineage_tnk"], sub["mal_CLDN4_mean_log1p_cp10k"],
                   c=colors[resp], s=70, label=resp, zorder=3)
        for _, r in sub.iterrows():
            ax.annotate(r["short"], (r["frac_lineage_tnk"], r["mal_CLDN4_mean_log1p_cp10k"]),
                        textcoords="offset points", xytext=(5, 2), fontsize=8)
    s = tests_df[(tests_df.family == "vs_TNK") & (tests_df.metric == "mal_CLDN4_mean_log1p_cp10k")
                 & (tests_df.tnk_def == "lineage_TNK")].iloc[0]
    ax.set_xlabel("T/NK fraction (lineage)")
    ax.set_ylabel("malignant CLDN4")
    ax.legend(fontsize=8)
    ax.set_title(f"n=6   ρ={s['rho']:.2f}   p={s['p']:.2f}")

    fig.suptitle("EXTRA  GSE291670  neoadjuvant anlotinib + camrelizumab  (not GSE207422)", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.outdir / "fig_extra_malignant_tacstd2_cldn4.png", dpi=160)
    plt.close(fig)

    summary = {
        "series": "GSE291670",
        "role": "EXTRA figure — second public neoadjuvant ICI lung scRNA",
        "not": "GSE207422 audit; author CopyKAT numbers are taken as given",
        "regimen": "neoadjuvant anlotinib + camrelizumab (PD-1)",
        "citation": "Xia et al. J Transl Med 2025; PMID from GEO GSE291670",
        "n_patients": 6,
        "n_MPR": 3,
        "n_NMPR": 3,
        "n_cells_qc": int(len(cells)),
        "n_malignant": int(cells["is_malignant"].sum()),
        "n_epithelial": int(cells["is_epithelial"].sum()),
        "n_lineage_tnk": int(cells["is_lineage_tnk"].sum()),
        "author_cell_labels_public": False,
        "why_not_preferred_alternatives": {
            "GSE253013": "No such live GEO series used. Closest name is GSE243013 (234 NSCLC post chemo-IO): CD45+ immune-only, 0 epithelial/malignant cells, 6.6 GB MTX.",
            "GSE266035": "Circulating T cells from one NSCLC patient on PD-1/CTLA-4. No tumor epithelium.",
            "GSE131907_tLung": "Treatment-naive LUAD atlas. No ICI / MPR labels. Already scored as leftover epithelial TACSTD2 vs T/NK (n=11, ρ=+0.09).",
        },
        "primary_tests": tests_df[tests_df.compartment == "malignant"][
            ["family", "gene", "metric", "tnk_def", "n", "mean_NMPR", "mean_MPR", "delta_NMPR_minus_MPR", "rho", "p"]
        ].to_dict(orient="records"),
    }
    with open(args.outdir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print("wrote", args.outdir, "n_cells", len(cells), flush=True)


if __name__ == "__main__":
    main()
