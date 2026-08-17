#!/usr/bin/env python3
"""ADDITIVE GSE293914 — H1975 EGFR xenograft snRNA (anti-PD-1 ± anti-CCL20).

Public series matrix + 10x MTX only. One GEO sample (GSM8893263) is a
Cell Ranger 7.1 raw Flex multiplex matrix (4 dominant probe barcodes).
GEO does not map probe barcode → treatment. Honest n = 1 library.

Scores human CLDN4 / TACSTD2 in tumor (epithelial) vs T/NK if present.
"""
from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Flex v1 translated probe barcodes (Cell Ranger collapses each BC to the
# lexicographically first 8-mer). Mapping is from 10x Flex v1 docs / 4-plex
# public examples. GEO does not say which BC is which treatment.
PROBE_ID = {
    "ACTTTAGG": "BC001",
    "AACGGGAA": "BC002",
    "AGTAGGCT": "BC003",
    "ATGTTGAC": "BC004",
}
FOCAL_PROBES = ("BC001", "BC002", "BC003", "BC004")

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
EXTRA = [
    "TACSTD2",
    "CLDN4",
    "PTPRC",
    "CD8A",
    "CD4",
    "MKI67",
    "CCL20",
    "EGFR",
    "NKX2-1",
    "NAPSA",
    "CD274",
    "PDCD1",
    "CCR6",
    "SFTPC",
    "AGER",
    "SCGB1A1",
]


def wanted() -> list[str]:
    g = list(EXTRA)
    for vs in LINEAGE_MARKERS.values():
        g.extend(vs)
    # unique, stable
    seen = set()
    out = []
    for x in g:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def load_features(path: Path) -> tuple[np.ndarray, dict[str, int], np.ndarray]:
    genes = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    genes = np.array(genes, dtype=object)
    idx: dict[str, int] = {}
    for i, g in enumerate(genes):
        if g not in idx:
            idx[g] = i
    mt = np.array([str(g).startswith("MT-") for g in genes])
    return genes, idx, mt


def load_barcodes(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    barcodes = []
    probes = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            b = line.strip().split("\t")[0]
            barcodes.append(b)
            core = b[:-2] if b.endswith("-1") else b
            probes.append(core[-8:] if len(core) >= 8 else "NA")
    barcodes = np.array(barcodes, dtype=object)
    probes = np.array(probes, dtype=object)
    probe_id = np.array([PROBE_ID.get(p, f"other_{p}") for p in probes], dtype=object)
    return barcodes, probes, probe_id


def stream_mtx(
    path: Path,
    n_genes: int,
    n_cells: int,
    keep_rows: dict[int, int],
    mt_rows: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """One-pass MTX stream. keep_rows maps 0-based gene row → dense column."""
    n_keep = 1 + max(keep_rows.values()) if keep_rows else 0
    umi = np.zeros(n_cells, dtype=np.uint32)
    nfeat = np.zeros(n_cells, dtype=np.uint16)
    mt_umi = np.zeros(n_cells, dtype=np.uint32)
    dense = np.zeros((n_keep, n_cells), dtype=np.uint16)

    with gzip.open(path, "rt") as fh:
        header = fh.readline()
        if "MatrixMarket" not in header:
            raise ValueError(f"not MatrixMarket: {header!r}")
        for line in fh:
            if line.startswith("%"):
                continue
            dims = line.split()
            ng, nc, nnz = int(dims[0]), int(dims[1]), int(dims[2])
            if ng != n_genes or nc != n_cells:
                raise ValueError(f"MTX dims {ng}x{nc} != features/barcodes {n_genes}x{n_cells}")
            break
        step = max(nnz // 20, 1)
        for i, line in enumerate(fh, 1):
            a, b, c = line.split()
            r = int(a) - 1
            col = int(b) - 1
            v = int(c)
            umi[col] += v
            nfeat[col] += 1
            if mt_rows[r]:
                mt_umi[col] += v
            j = keep_rows.get(r)
            if j is not None:
                if v > 65535:
                    v = 65535
                dense[j, col] = v
            if i % step == 0:
                print(f"  MTX {100 * i / nnz:.0f}% ({i:,}/{nnz:,})", flush=True)
    return umi, nfeat, mt_umi, dense


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


def summarize_expr(x: np.ndarray) -> dict:
    x = np.asarray(x, float)
    if x.size == 0:
        return {"n": 0, "mean": None, "median": None, "pct_pos": None}
    return {
        "n": int(x.size),
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "pct_pos": float(np.mean(x > 0) * 100.0),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("data/GSE293914"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/gse293914_cldn4"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    feat = args.datadir / "GSM8893263_features.tsv.gz"
    barc = args.datadir / "GSM8893263_barcodes.tsv.gz"
    mtxp = args.datadir / "GSM8893263_matrix.mtx.gz"
    for p in (feat, barc, mtxp):
        if not p.exists():
            raise SystemExit(f"missing {p}")

    print("loading features/barcodes", flush=True)
    genes, gidx, mt_mask = load_features(feat)
    barcodes, probe_seq, probe_id = load_barcodes(barc)
    n_genes, n_cells = len(genes), len(barcodes)
    print(f"  genes={n_genes:,} barcodes={n_cells:,}", flush=True)
    print("  raw probe counts:", dict(Counter(probe_id).most_common()), flush=True)

    keep_names = wanted()
    keep_rows: dict[int, int] = {}
    keep_order: list[str] = []
    missing = []
    for g in keep_names:
        if g in gidx:
            keep_rows[gidx[g]] = len(keep_order)
            keep_order.append(g)
        else:
            missing.append(g)
    print(f"  extracted genes={len(keep_order)} missing={missing}", flush=True)

    print("streaming MTX", flush=True)
    umi, nfeat, mt_umi, dense = stream_mtx(mtxp, n_genes, n_cells, keep_rows, mt_mask)
    pct_mt = np.where(umi > 0, 100.0 * mt_umi / umi, 0.0)

    # Author-stated Cell Ranger filters (GEO data processing).
    qc = (nfeat >= 200) & (nfeat < 8000) & (pct_mt < 10)
    print(
        f"QC pass {int(qc.sum()):,} / {n_cells:,} "
        f"(genes>=200, genes<8000, mt<10%)",
        flush=True,
    )

    barcodes = barcodes[qc]
    probe_id = probe_id[qc]
    umi, nfeat, pct_mt = umi[qc], nfeat[qc], pct_mt[qc]
    dense = dense[:, qc]
    n = int(qc.sum())

    expr = {g: dense[j].astype(np.float32) for j, g in enumerate(keep_order)}
    scale = np.where(umi > 0, 1e4 / umi.astype(np.float64), 0.0).astype(np.float32)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
    scores = {name: module(log_cp, gs, n) for name, gs in LINEAGE_MARKERS.items()}
    lineage = assign_lineage(scores, log_cp.get("CD3E", log_cp.get("CD3D", np.zeros(n))))

    is_tumor = lineage == "Epithelial"
    is_tnk = np.isin(lineage, ["T", "NK"])
    is_t = lineage == "T"
    is_nk = lineage == "NK"
    umi_tnk = (expr.get("CD3E", np.zeros(n)) >= 1) | (expr.get("CD8A", np.zeros(n)) >= 1) | (
        expr.get("NKG7", np.zeros(n)) >= 1
    )

    cells = pd.DataFrame(
        {
            "barcode": barcodes,
            "probe_id": probe_id,
            "n_umi": umi,
            "n_genes": nfeat,
            "pct_mt": pct_mt.astype(np.float32),
            "lineage": lineage,
            "is_tumor": is_tumor,
            "is_tnk": is_tnk,
            "TACSTD2_umi": expr.get("TACSTD2", np.zeros(n)),
            "CLDN4_umi": expr.get("CLDN4", np.zeros(n)),
            "TACSTD2_log1p_cp10k": log_cp.get("TACSTD2", np.zeros(n)),
            "CLDN4_log1p_cp10k": log_cp.get("CLDN4", np.zeros(n)),
            "CCL20_log1p_cp10k": log_cp.get("CCL20", np.zeros(n)),
            "EGFR_log1p_cp10k": log_cp.get("EGFR", np.zeros(n)),
            "PTPRC_log1p_cp10k": log_cp.get("PTPRC", np.zeros(n)),
        }
    )

    def gene_block(mask: np.ndarray, gene: str) -> dict:
        um = cells.loc[mask, f"{gene}_umi"].to_numpy()
        lg = cells.loc[mask, f"{gene}_log1p_cp10k"].to_numpy()
        s = summarize_expr(lg)
        s["pct_umi_ge1"] = float(np.mean(um >= 1) * 100.0) if um.size else None
        return s

    lineage_counts = (
        cells.groupby(["probe_id", "lineage"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(["probe_id", "n_cells"], ascending=[True, False])
    )
    lineage_counts.to_csv(args.outdir / "lineage_counts.tsv", sep="\t", index=False)

    rows = []
    for pid, sub in cells.groupby("probe_id", observed=True):
        row = {
            "probe_id": pid,
            "n_qc": int(len(sub)),
            "n_tumor": int(sub["is_tumor"].sum()),
            "n_tnk": int(sub["is_tnk"].sum()),
            "n_T": int((sub["lineage"] == "T").sum()),
            "n_NK": int((sub["lineage"] == "NK").sum()),
            "frac_tumor": float(sub["is_tumor"].mean()),
            "frac_tnk": float(sub["is_tnk"].mean()),
        }
        for gene in ("TACSTD2", "CLDN4"):
            for comp, mask in (
                ("tumor", sub["is_tumor"].to_numpy()),
                ("tnk", sub["is_tnk"].to_numpy()),
                ("other", (~sub["is_tumor"] & ~sub["is_tnk"]).to_numpy()),
            ):
                um = sub.loc[mask, f"{gene}_umi"].to_numpy()
                lg = sub.loc[mask, f"{gene}_log1p_cp10k"].to_numpy()
                row[f"{gene}_{comp}_n"] = int(mask.sum())
                row[f"{gene}_{comp}_mean_log1p"] = float(np.mean(lg)) if mask.sum() else np.nan
                row[f"{gene}_{comp}_pct_pos"] = float(np.mean(um >= 1) * 100.0) if mask.sum() else np.nan
        row["CCL20_tumor_mean_log1p"] = float(sub.loc[sub["is_tumor"], "CCL20_log1p_cp10k"].mean()) if sub["is_tumor"].any() else np.nan
        row["EGFR_tumor_mean_log1p"] = float(sub.loc[sub["is_tumor"], "EGFR_log1p_cp10k"].mean()) if sub["is_tumor"].any() else np.nan
        rows.append(row)
    per_probe = pd.DataFrame(rows).sort_values("n_qc", ascending=False)
    per_probe.to_csv(args.outdir / "per_probe.tsv", sep="\t", index=False)

    overall = {
        "series": "GSE293914",
        "sample": "GSM8893263",
        "label": "XENOGRAFT",
        "model": "H1975 EGFR L858R/T790M humanized NSG + human PBMC",
        "treatments_on_geo_title": ["isotype", "anti-PD-1", "anti-CCL20", "combination"],
        "geo_characteristics_treatment": "isotype",
        "n_geo_samples": 1,
        "n_libraries": 1,
        "honest_n": 1,
        "probe_to_treatment_map_on_geo": False,
        "raw_barcodes": int(n_cells),
        "qc_cells": n,
        "qc": {"min_genes": 200, "max_genes": 8000, "max_pct_mt": 10},
        "missing_genes": missing,
        "lineage_counts": {k: int(v) for k, v in Counter(lineage).items()},
        "n_tumor": int(is_tumor.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_umi_tnk_proxy": int(umi_tnk.sum()),
        "focal_probe_qc": {
            pid: int((probe_id == pid).sum()) for pid in FOCAL_PROBES
        },
        "TACSTD2_tumor": gene_block(is_tumor, "TACSTD2"),
        "TACSTD2_tnk": gene_block(is_tnk, "TACSTD2"),
        "CLDN4_tumor": gene_block(is_tumor, "CLDN4"),
        "CLDN4_tnk": gene_block(is_tnk, "CLDN4"),
        "note": (
            "One multiplexed Flex library. Four dominant probe barcodes match a "
            "4-plex (isotype / aPD-1 / aCCL20 / combo) but GEO MTX/series matrix "
            "do not label which barcode is which arm. Do not treat probe groups "
            "as independent biological n. Cell-level p-values are not reported."
        ),
    }
    (args.outdir / "summary.json").write_text(json.dumps(overall, indent=2) + "\n")

    stats_rows = []
    for gene in ("TACSTD2", "CLDN4"):
        for comp, mask in (("tumor", is_tumor), ("tnk", is_tnk), ("T", is_t), ("NK", is_nk)):
            blk = gene_block(mask, gene)
            stats_rows.append({"gene": gene, "compartment": comp, "scope": "pooled_library", **blk})
        for pid in FOCAL_PROBES:
            m = probe_id == pid
            for comp, mask in (("tumor", is_tumor & m), ("tnk", is_tnk & m)):
                blk = gene_block(mask, gene)
                stats_rows.append(
                    {"gene": gene, "compartment": comp, "scope": pid, **blk}
                )
    pd.DataFrame(stats_rows).to_csv(args.outdir / "stats.tsv", sep="\t", index=False)

    # Keep a compact QC table, not 2M raw barcodes.
    cells.to_csv(args.outdir / "cell_calls.tsv.gz", sep="\t", index=False, compression="gzip")

    # Figure
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.4))
    fig.suptitle(
        "GSE293914 H1975 EGFR xenograft (humanized) — CLDN4 / TACSTD2\n"
        "n=1 multiplexed Flex library; probe≠treatment map not on GEO",
        fontsize=10,
    )

    def _violin(ax, tumor, tnk, ylab, title):
        data = [tumor, tnk]
        parts = ax.violinplot(data, showmeans=True, showextrema=False, widths=0.7)
        for i, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(["#d95f02", "#7570b3"][i])
            pc.set_alpha(0.7)
        parts["cmeans"].set_color("black")
        ax.scatter([1, 2], [np.mean(tumor) if tumor.size else 0, np.mean(tnk) if tnk.size else 0], c="k", zorder=3, s=12)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(
            [f"tumor\nn={tumor.size:,}", f"T/NK\nn={tnk.size:,}"]
        )
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    _violin(
        axes[0, 0],
        cells.loc[is_tumor, "TACSTD2_log1p_cp10k"].to_numpy(),
        cells.loc[is_tnk, "TACSTD2_log1p_cp10k"].to_numpy(),
        "log1p(CP10k)",
        "TACSTD2  tumor vs T/NK (pooled library)",
    )
    _violin(
        axes[0, 1],
        cells.loc[is_tumor, "CLDN4_log1p_cp10k"].to_numpy(),
        cells.loc[is_tnk, "CLDN4_log1p_cp10k"].to_numpy(),
        "log1p(CP10k)",
        "CLDN4  tumor vs T/NK (pooled library)",
    )

    ax = axes[1, 0]
    x = np.arange(len(FOCAL_PROBES))
    w = 0.35
    t_means = []
    k_means = []
    t_n = []
    k_n = []
    for pid in FOCAL_PROBES:
        m = probe_id == pid
        t_means.append(float(cells.loc[is_tumor & m, "TACSTD2_log1p_cp10k"].mean()) if (is_tumor & m).any() else 0)
        k_means.append(float(cells.loc[is_tnk & m, "TACSTD2_log1p_cp10k"].mean()) if (is_tnk & m).any() else 0)
        t_n.append(int((is_tumor & m).sum()))
        k_n.append(int((is_tnk & m).sum()))
    ax.bar(x - w / 2, t_means, w, color="#d95f02", label="tumor")
    ax.bar(x + w / 2, k_means, w, color="#7570b3", label="T/NK")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p}\ntum {n:,}" for p, n in zip(FOCAL_PROBES, t_n)], fontsize=7)
    ax.set_ylabel("TACSTD2 mean log1p(CP10k)")
    ax.set_title("Unlabeled Flex arms (no GEO treatment map)", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1, 1]
    t_means = []
    k_means = []
    for pid in FOCAL_PROBES:
        m = probe_id == pid
        t_means.append(float(cells.loc[is_tumor & m, "CLDN4_log1p_cp10k"].mean()) if (is_tumor & m).any() else 0)
        k_means.append(float(cells.loc[is_tnk & m, "CLDN4_log1p_cp10k"].mean()) if (is_tnk & m).any() else 0)
    ax.bar(x - w / 2, t_means, w, color="#d95f02", label="tumor")
    ax.bar(x + w / 2, k_means, w, color="#7570b3", label="T/NK")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p}\nT/NK {n:,}" for p, n in zip(FOCAL_PROBES, k_n)], fontsize=7)
    ax.set_ylabel("CLDN4 mean log1p(CP10k)")
    ax.set_title("Unlabeled Flex arms (no GEO treatment map)", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(args.outdir / "fig_xenograft_cldn4_tacstd2.png", dpi=160)
    plt.close(fig)

    print(json.dumps({k: overall[k] for k in ("qc_cells", "n_tumor", "n_tnk", "lineage_counts", "TACSTD2_tumor", "CLDN4_tumor", "TACSTD2_tnk", "CLDN4_tnk")}, indent=2))
    print("wrote", args.outdir)


if __name__ == "__main__":
    main()
