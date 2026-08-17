#!/usr/bin/env python3
"""GSE316655 public MTX/TSV: CLDN / LILRB / myeloid NF-κB-STAT / T-NK / MHC-I.

Species/model (do not drop this label):
  Human CD45+ FACS scRNA from SK-MEL-5 tumors (and PB) grown in NSG-SGM3 mice
  humanized with human cord-blood CD34+ cells. Cell Ranger GRCh38_and_mm10-2020-A.
  Liu et al. Sci Immunol 2026 (PMID 41931598). GEO GSE316655.

Honest n: 4 public 10x libraries, 1 donor label ("ND"), 1 library per
tissue × treatment cell. Cell-level tests are descriptive / pseudoreplication.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread

from gene_sets import (
    BARRIER,
    CLAUDINS,
    CYTOTOX,
    EXHAUST,
    LILR,
    LINEAGE,
    LR_PAIRS,
    MHCI,
    MHCII,
    NFKB,
    STAT,
    SUPPRESS_MYELOID,
    all_panel_genes,
)

LIBRARIES = [
    {
        "gsm": "GSM9457798",
        "stem": "GSM9457798_ND_anti_B2",
        "title": "ND tumor LILRB2 treatment",
        "tissue": "tumor",
        "treatment": "anti-LILRB2",
        "library": "ND_anti_B2",
        "short": "Tumor anti-LILRB2",
    },
    {
        "gsm": "GSM9457799",
        "stem": "GSM9457799_ND_CTR",
        "title": "ND tumor Isotype",
        "tissue": "tumor",
        "treatment": "isotype",
        "library": "ND_Ctrl",
        "short": "Tumor isotype",
    },
    {
        "gsm": "GSM9457800",
        "stem": "GSM9457800_PB_anti_B2",
        "title": "PB LILRB2 treatment",
        "tissue": "PB",
        "treatment": "anti-LILRB2",
        "library": "PB_anti_B2",
        "short": "PB anti-LILRB2",
    },
    {
        "gsm": "GSM9457801",
        "stem": "GSM9457801_PB_CTR",
        "title": "PB Isotype",
        "tissue": "PB",
        "treatment": "isotype",
        "library": "PB_Ctrl",
        "short": "PB isotype",
    },
]

MODEL = (
    "Human CD45+ FACS | SK-MEL-5 in NSG-SGM3 humanized mice | "
    "GRCh38+mm10 | GSE316655"
)

QC_MIN_HUMAN_UMI = 200
QC_MIN_HUMAN_GENES = 100
QC_MIN_HUMAN_FRAC = 0.80
QC_MAX_MITO_FRAC = 0.20

LINEAGE_MIN = 0.20
LINEAGE_MARGIN = 0.05

COLORS = {
    "Tumor anti-LILRB2": "#1b9e77",
    "Tumor isotype": "#d95f02",
    "PB anti-LILRB2": "#7570b3",
    "PB isotype": "#e7298a",
    "myeloid": "#e69f00",
    "T": "#56b4e9",
    "NK": "#009e73",
    "B": "#f0e442",
    "pDC": "#0072b2",
    "melanoma": "#d55e00",
    "epithelial": "#cc79a7",
    "other": "#999999",
}


def strip_name(raw: str) -> tuple[str, str]:
    if raw.startswith("GRCh38_"):
        return "human", raw[7:]
    if raw.startswith("mm10___"):
        return "mouse", raw[7:]
    if raw.startswith("mm10_"):
        return "mouse", raw[5:]
    return "unknown", raw


def read_features(path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            gid = parts[0]
            gname = parts[1] if len(parts) > 1 else parts[0]
            species, symbol = strip_name(gname)
            rows.append((i, gid, gname, species, symbol))
    return pd.DataFrame(rows, columns=["idx", "gene_id", "raw_name", "species", "symbol"])


def read_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [line.strip() for line in fh if line.strip()]


def mean_log1p_cp10k(umi: np.ndarray, lib: np.ndarray, genes: list[str], present: dict[str, np.ndarray]) -> np.ndarray:
    mats = []
    for g in genes:
        if g not in present:
            continue
        cp = np.where(lib > 0, present[g] / lib * 1e4, 0.0)
        mats.append(np.log1p(cp))
    if not mats:
        return np.full(len(umi) if umi is not None else lib.shape[0], np.nan, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0).astype(np.float32)


def assign_lineage(present: dict[str, np.ndarray], lib: np.ndarray) -> np.ndarray:
    names = list(LINEAGE)
    scores = []
    for key in names:
        scores.append(mean_log1p_cp10k(None, lib, LINEAGE[key], present))
    scores = np.vstack(scores)
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < LINEAGE_MIN) | ((best_val - second) < LINEAGE_MARGIN)] = "other"
    return assigned


def load_library(raw: Path, spec: dict) -> dict:
    feat = read_features(raw / f"{spec['stem']}_features.tsv.gz")
    barcodes = read_barcodes(raw / f"{spec['stem']}_barcodes.tsv.gz")
    mtx_path = raw / f"{spec['stem']}_matrix.mtx.gz"
    print(f"reading {mtx_path}", flush=True)
    mat = mmread(mtx_path).tocsr().astype(np.float32)
    if mat.shape[0] != len(feat) or mat.shape[1] != len(barcodes):
        raise SystemExit(
            f"{spec['gsm']} shape {mat.shape} != features {len(feat)} x barcodes {len(barcodes)}"
        )

    human_mask = (feat["species"] == "human").to_numpy()
    mouse_mask = (feat["species"] == "mouse").to_numpy()
    mito_mask = human_mask & feat["symbol"].str.startswith("MT-").to_numpy()

    human_umi = np.asarray(mat[human_mask].sum(axis=0)).ravel()
    mouse_umi = np.asarray(mat[mouse_mask].sum(axis=0)).ravel()
    mito_umi = np.asarray(mat[mito_mask].sum(axis=0)).ravel()
    total_umi = human_umi + mouse_umi
    human_genes = np.asarray((mat[human_mask] > 0).sum(axis=0)).ravel()
    mouse_genes = np.asarray((mat[mouse_mask] > 0).sum(axis=0)).ravel()
    human_frac = np.divide(human_umi, total_umi, out=np.zeros_like(human_umi), where=total_umi > 0)
    mito_frac = np.divide(mito_umi, human_umi, out=np.zeros_like(mito_umi), where=human_umi > 0)

    keep = (
        (human_umi >= QC_MIN_HUMAN_UMI)
        & (human_genes >= QC_MIN_HUMAN_GENES)
        & (human_frac >= QC_MIN_HUMAN_FRAC)
        & (mito_frac <= QC_MAX_MITO_FRAC)
    )

    panel = all_panel_genes()
    symbol_to_idx = {}
    for i, row in feat.iterrows():
        if row["species"] != "human":
            continue
        # first occurrence wins (unique symbols in this ref)
        symbol_to_idx.setdefault(row["symbol"], int(row["idx"]))

    present = {}
    missing = []
    for g in panel:
        if g in symbol_to_idx:
            present[g] = np.asarray(mat[symbol_to_idx[g]].todense()).ravel()
        else:
            missing.append(g)

    qc_lib = human_umi[keep]
    qc_present = {g: v[keep] for g, v in present.items()}
    lineage = assign_lineage(qc_present, qc_lib)

    scores = {
        "cldn4": mean_log1p_cp10k(None, qc_lib, ["CLDN4"], qc_present),
        "claudin_mean": mean_log1p_cp10k(None, qc_lib, CLAUDINS, qc_present),
        "nfkb": mean_log1p_cp10k(None, qc_lib, NFKB, qc_present),
        "stat": mean_log1p_cp10k(None, qc_lib, STAT, qc_present),
        "mhci": mean_log1p_cp10k(None, qc_lib, MHCI, qc_present),
        "mhcii": mean_log1p_cp10k(None, qc_lib, MHCII, qc_present),
        "cytotox": mean_log1p_cp10k(None, qc_lib, CYTOTOX, qc_present),
        "exhaust": mean_log1p_cp10k(None, qc_lib, EXHAUST, qc_present),
        "suppress_myeloid": mean_log1p_cp10k(None, qc_lib, SUPPRESS_MYELOID, qc_present),
        "lilrb2": mean_log1p_cp10k(None, qc_lib, ["LILRB2"], qc_present),
        "lilrb1": mean_log1p_cp10k(None, qc_lib, ["LILRB1"], qc_present),
        "lilrb5": mean_log1p_cp10k(None, qc_lib, ["LILRB5"], qc_present),
        "barrier": mean_log1p_cp10k(None, qc_lib, BARRIER, qc_present),
    }

    cells = pd.DataFrame(
        {
            "barcode": np.asarray(barcodes)[keep],
            "gsm": spec["gsm"],
            "library": spec["library"],
            "short": spec["short"],
            "tissue": spec["tissue"],
            "treatment": spec["treatment"],
            "human_umi": qc_lib,
            "mouse_umi": mouse_umi[keep],
            "human_genes": human_genes[keep],
            "human_frac": human_frac[keep],
            "mito_frac": mito_frac[keep],
            "lineage": lineage,
        }
    )
    for k, v in scores.items():
        cells[k] = v
    for g, v in qc_present.items():
        cells[f"umi_{g}"] = v.astype(np.float32)

    raw_qc = pd.DataFrame(
        {
            "barcode": barcodes,
            "gsm": spec["gsm"],
            "short": spec["short"],
            "human_umi": human_umi,
            "mouse_umi": mouse_umi,
            "total_umi": total_umi,
            "human_genes": human_genes,
            "mouse_genes": mouse_genes,
            "human_frac": human_frac,
            "mito_frac": mito_frac,
            "pass_qc": keep,
        }
    )
    return {
        "spec": spec,
        "cells": cells,
        "raw_qc": raw_qc,
        "n_barcodes": len(barcodes),
        "n_qc": int(keep.sum()),
        "n_features": int(len(feat)),
        "n_human_features": int(human_mask.sum()),
        "n_mouse_features": int(mouse_mask.sum()),
        "missing_genes": missing,
        "present_genes": sorted(present),
    }


def pct_pos(series: pd.Series) -> float:
    if len(series) == 0:
        return float("nan")
    return float((series > 0).mean() * 100.0)


def mean_or_nan(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.nanmean(x)) if x.size else float("nan")


def describe_group(df: pd.DataFrame, prefix: str) -> dict:
    out = {
        f"{prefix}_n": int(len(df)),
        f"{prefix}_frac": float("nan"),
    }
    for col in [
        "cldn4",
        "claudin_mean",
        "nfkb",
        "stat",
        "mhci",
        "mhcii",
        "cytotox",
        "exhaust",
        "suppress_myeloid",
        "lilrb2",
        "lilrb1",
        "lilrb5",
        "barrier",
    ]:
        out[f"{prefix}_{col}_mean"] = mean_or_nan(df[col].to_numpy()) if col in df else float("nan")
    for g in ["CLDN4", "CLDN18", "LILRB2", "LILRB1", "LILRB5", "HLA-A", "B2M", "CD3D", "NKG7", "LYZ"]:
        col = f"umi_{g}"
        if col in df:
            out[f"{prefix}_{g}_pct"] = pct_pos(df[col])
            lib = df["human_umi"].to_numpy()
            umi = df[col].to_numpy()
            cp = np.where(lib > 0, umi / lib * 1e4, 0.0)
            out[f"{prefix}_{g}_mean_log1p_cp10k"] = float(np.mean(np.log1p(cp)))
        else:
            out[f"{prefix}_{g}_pct"] = float("nan")
            out[f"{prefix}_{g}_mean_log1p_cp10k"] = float("nan")
    return out


def cheap_lr(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for short, sub in cells.groupby("short", sort=False):
        lib = sub["human_umi"].to_numpy()
        myeloid = sub[sub["lineage"] == "myeloid"]
        tnk = sub[sub["lineage"].isin(["T", "NK"])]
        for lig, rec, note in LR_PAIRS:
            lcol, rcol = f"umi_{lig}", f"umi_{rec}"
            if lcol not in sub or rcol not in sub:
                continue
            l_all = sub[lcol].to_numpy()
            r_all = sub[rcol].to_numpy()
            l_cp = np.where(lib > 0, l_all / lib * 1e4, 0.0)
            r_cp = np.where(lib > 0, r_all / lib * 1e4, 0.0)
            rec_my = myeloid[rcol].to_numpy() if len(myeloid) and rcol in myeloid else np.array([])
            rec_tnk = tnk[rcol].to_numpy() if len(tnk) and rcol in tnk else np.array([])
            rows.append(
                {
                    "library": short,
                    "tissue": sub["tissue"].iloc[0],
                    "treatment": sub["treatment"].iloc[0],
                    "ligand": lig,
                    "receptor": rec,
                    "note": note,
                    "n_cells": int(len(sub)),
                    "n_myeloid": int(len(myeloid)),
                    "n_TNK": int(len(tnk)),
                    "pct_ligand": pct_pos(sub[lcol]),
                    "pct_receptor": pct_pos(sub[rcol]),
                    "pct_product": pct_pos(sub[lcol]) * pct_pos(sub[rcol]) / 100.0,
                    "mean_ligand_log1p_cp10k": float(np.mean(np.log1p(l_cp))),
                    "mean_receptor_log1p_cp10k": float(np.mean(np.log1p(r_cp))),
                    "pct_receptor_myeloid": pct_pos(myeloid[rcol]) if len(myeloid) else float("nan"),
                    "pct_receptor_TNK": pct_pos(tnk[rcol]) if len(tnk) else float("nan"),
                    "mean_receptor_myeloid": float(np.mean(np.log1p(np.where(myeloid["human_umi"] > 0, rec_my / myeloid["human_umi"] * 1e4, 0.0))))
                    if len(myeloid)
                    else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def claudin_table(cells: pd.DataFrame, present: set[str]) -> pd.DataFrame:
    rows = []
    for g in CLAUDINS:
        col = f"umi_{g}"
        if g not in present or col not in cells:
            for short, sub in cells.groupby("short", sort=False):
                rows.append(
                    {
                        "gene": g,
                        "library": short,
                        "present": False,
                        "n": int(len(sub)),
                        "pct_pos": np.nan,
                        "mean_log1p_cp10k": np.nan,
                    }
                )
            continue
        for short, sub in cells.groupby("short", sort=False):
            lib = sub["human_umi"].to_numpy()
            umi = sub[col].to_numpy()
            cp = np.where(lib > 0, umi / lib * 1e4, 0.0)
            rows.append(
                {
                    "gene": g,
                    "library": short,
                    "present": True,
                    "n": int(len(sub)),
                    "pct_pos": pct_pos(sub[col]),
                    "mean_log1p_cp10k": float(np.mean(np.log1p(cp))),
                }
            )
    return pd.DataFrame(rows)


def lilr_table(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for g in LILR:
        col = f"umi_{g}"
        if col not in cells:
            continue
        for (short, lin), sub in cells.groupby(["short", "lineage"], sort=False):
            lib = sub["human_umi"].to_numpy()
            umi = sub[col].to_numpy()
            cp = np.where(lib > 0, umi / lib * 1e4, 0.0)
            rows.append(
                {
                    "gene": g,
                    "library": short,
                    "lineage": lin,
                    "n": int(len(sub)),
                    "pct_pos": pct_pos(sub[col]),
                    "mean_log1p_cp10k": float(np.mean(np.log1p(cp))),
                }
            )
    return pd.DataFrame(rows)


def exploratory_mwu(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 5 or len(b) < 5:
        return {
            "n_a": int(len(a)),
            "n_b": int(len(b)),
            "mean_a": mean_or_nan(a),
            "mean_b": mean_or_nan(b),
            "median_a": float(np.median(a)) if len(a) else float("nan"),
            "median_b": float(np.median(b)) if len(b) else float("nan"),
            "U": float("nan"),
            "p": float("nan"),
            "note": "too_few_cells",
        }
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "U": float(U),
        "p": float(p),
        "note": "cell_level_pseudoreplication_n_library=1_vs_1",
    }


def style_axes(ax, title: str) -> None:
    ax.set_title(title, fontsize=11)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def add_footer(fig) -> None:
    fig.text(0.01, 0.01, MODEL, fontsize=7, color="#444444")


def plot_qc(raw: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.0))
    order = [s["short"] for s in LIBRARIES]
    raw["short"] = pd.Categorical(raw["short"], order)
    # human fraction
    ax = axes[0]
    data = [raw.loc[raw["short"] == s, "human_frac"].to_numpy() for s in order]
    bp = ax.boxplot(data, tick_labels=[s.replace(" ", "\n") for s in order], patch_artist=True, showfliers=False)
    for patch, s in zip(bp["boxes"], order):
        patch.set_facecolor(COLORS[s])
        patch.set_alpha(0.7)
    ax.axhline(QC_MIN_HUMAN_FRAC, color="black", ls="--", lw=0.8)
    ax.set_ylabel("Human UMI fraction")
    style_axes(ax, "Species mix (all barcodes)")
    # human UMI
    ax = axes[1]
    data = [np.log10(np.clip(raw.loc[raw["short"] == s, "human_umi"].to_numpy(), 1, None)) for s in order]
    bp = ax.boxplot(data, tick_labels=[s.replace(" ", "\n") for s in order], patch_artist=True, showfliers=False)
    for patch, s in zip(bp["boxes"], order):
        patch.set_facecolor(COLORS[s])
        patch.set_alpha(0.7)
    ax.set_ylabel("log10 human UMI")
    style_axes(ax, "Library size (human genes)")
    # pass rates
    ax = axes[2]
    rates = raw.groupby("short", observed=False)["pass_qc"].mean().reindex(order) * 100
    ax.bar(range(len(order)), rates.to_numpy(), color=[COLORS[s] for s in order])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([s.replace(" ", "\n") for s in order])
    ax.set_ylabel("% barcodes passing QC")
    style_axes(ax, "QC pass (human-high CD45-like)")
    fig.suptitle("GSE316655 QC — dual-genome humanized-mouse scRNA", fontsize=12)
    add_footer(fig)
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(out / "fig1_qc_species.png", dpi=160)
    fig.savefig(out / "fig1_qc_species.pdf")
    plt.close(fig)


def plot_composition(cells: pd.DataFrame, out: Path) -> None:
    order = [s["short"] for s in LIBRARIES]
    lin_order = ["myeloid", "T", "NK", "B", "pDC", "melanoma", "epithelial", "other"]
    tab = (
        cells.groupby(["short", "lineage"], observed=False)
        .size()
        .unstack(fill_value=0)
        .reindex(index=order, columns=lin_order, fill_value=0)
    )
    frac = tab.div(tab.sum(axis=1), axis=0) * 100
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    bottom = np.zeros(len(order))
    x = np.arange(len(order))
    for lin in lin_order:
        vals = frac[lin].to_numpy() if lin in frac else np.zeros(len(order))
        ax.bar(x, vals, bottom=bottom, color=COLORS[lin], label=lin)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace(" ", "\n") for s in order])
    ax.set_ylabel("% of QC cells")
    ax.legend(frameon=False, ncol=4, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    style_axes(ax, "Marker lineage composition (not author labels)")
    fig.suptitle("GSE316655 lineage mix — FACS human CD45+ (SK-MEL-5 / NSG-SGM3)", fontsize=12)
    add_footer(fig)
    fig.tight_layout(rect=(0, 0.05, 1, 0.90))
    fig.savefig(out / "fig2_lineage_composition.png", dpi=160)
    fig.savefig(out / "fig2_lineage_composition.pdf")
    plt.close(fig)


def plot_claudins(cldn: pd.DataFrame, out: Path) -> None:
    order = [s["short"] for s in LIBRARIES]
    genes = [g for g in CLAUDINS if g in set(cldn["gene"])]
    mat = (
        cldn.pivot_table(index="gene", columns="library", values="pct_pos", aggfunc="mean")
        .reindex(index=genes, columns=order)
    )
    fig, ax = plt.subplots(figsize=(7.2, 8.4))
    im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="YlOrRd", vmin=0, vmax=max(5, float(np.nanmax(mat.to_numpy()))))
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([s.replace(" ", "\n") for s in order], fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="% cells UMI>0")
    style_axes(ax, "Claudin detection in FACS CD45+ cells")
    fig.suptitle("GSE316655 human claudins — residual / leak expected", fontsize=12)
    add_footer(fig)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(out / "fig3_claudin_detection.png", dpi=160)
    fig.savefig(out / "fig3_claudin_detection.pdf")
    plt.close(fig)


def plot_lilrb(cells: pd.DataFrame, out: Path) -> None:
    genes = [g for g in LILR if f"umi_{g}" in cells.columns]
    order = [s["short"] for s in LIBRARIES]
    myeloid = cells[cells["lineage"] == "myeloid"]
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.4))
    ax = axes[0]
    x = np.arange(len(genes))
    width = 0.2
    for i, short in enumerate(order):
        sub = myeloid[myeloid["short"] == short]
        vals = [pct_pos(sub[f"umi_{g}"]) if len(sub) else np.nan for g in genes]
        ax.bar(x + (i - 1.5) * width, vals, width, color=COLORS[short], label=short)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("% myeloid cells UMI>0")
    style_axes(ax, "LILR family in myeloid cells")
    ax.legend(frameon=False, fontsize=7)
    ax = axes[1]
    for short in order:
        sub = myeloid[myeloid["short"] == short]
        if sub.empty:
            continue
        ax.scatter(
            np.random.default_rng(0).normal(order.index(short), 0.08, size=len(sub)),
            sub["lilrb2"],
            s=6,
            alpha=0.25,
            color=COLORS[short],
            label=short,
        )
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([s.replace(" ", "\n") for s in order], fontsize=8)
    ax.set_ylabel("LILRB2 log1p(CP10k)")
    style_axes(ax, "LILRB2 in myeloid cells (1 library / arm)")
    fig.suptitle("GSE316655 LILRB — human myeloid checkpoint (humanized mouse)", fontsize=12)
    add_footer(fig)
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    fig.savefig(out / "fig4_lilrb_family.png", dpi=160)
    fig.savefig(out / "fig4_lilrb_family.pdf")
    plt.close(fig)


def plot_myeloid_scores(cells: pd.DataFrame, out: Path) -> None:
    order = [s["short"] for s in LIBRARIES]
    myeloid = cells[cells["lineage"] == "myeloid"]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.2))
    for ax, col, title in [
        (axes[0], "nfkb", "Myeloid NF-κB score"),
        (axes[1], "stat", "Myeloid STAT score"),
        (axes[2], "suppress_myeloid", "Myeloid suppressors (IDO1/IL10/PD-L1/…)"),
    ]:
        data = [myeloid.loc[myeloid["short"] == s, col].to_numpy() for s in order]
        bp = ax.boxplot(data, tick_labels=[s.replace(" ", "\n") for s in order], patch_artist=True, showfliers=False)
        for patch, s in zip(bp["boxes"], order):
            patch.set_facecolor(COLORS[s])
            patch.set_alpha(0.7)
        ax.set_ylabel("mean log1p(CP10k)")
        style_axes(ax, title)
    fig.suptitle("GSE316655 myeloid programs — n=1 library per arm (descriptive)", fontsize=12)
    add_footer(fig)
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    fig.savefig(out / "fig5_myeloid_nfkb_stat.png", dpi=160)
    fig.savefig(out / "fig5_myeloid_nfkb_stat.pdf")
    plt.close(fig)


def plot_tnk_mhci(cells: pd.DataFrame, out: Path) -> None:
    order = [s["short"] for s in LIBRARIES]
    tnk = cells[cells["lineage"].isin(["T", "NK"])]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.2))
    for ax, df, col, title in [
        (axes[0], tnk, "cytotox", "T/NK cytotoxicity"),
        (axes[1], tnk, "exhaust", "T/NK exhaustion"),
        (axes[2], cells, "mhci", "MHC-I cassette (all QC cells)"),
    ]:
        data = [df.loc[df["short"] == s, col].to_numpy() for s in order]
        bp = ax.boxplot(data, tick_labels=[s.replace(" ", "\n") for s in order], patch_artist=True, showfliers=False)
        for patch, s in zip(bp["boxes"], order):
            patch.set_facecolor(COLORS[s])
            patch.set_alpha(0.7)
        ax.set_ylabel("mean log1p(CP10k)")
        style_axes(ax, title)
    fig.suptitle("GSE316655 T/NK and MHC-I — humanized-mouse CD45+", fontsize=12)
    add_footer(fig)
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    fig.savefig(out / "fig6_tnk_mhci.png", dpi=160)
    fig.savefig(out / "fig6_tnk_mhci.pdf")
    plt.close(fig)


def plot_lr(lr: pd.DataFrame, out: Path) -> None:
    focus = lr[lr["note"].str.contains("CLDN|MHC-I–LILRB2|paper")].copy()
    if focus.empty:
        focus = lr.copy()
    pairs = (focus["ligand"] + "–" + focus["receptor"]).drop_duplicates().tolist()
    order = [s["short"] for s in LIBRARIES]
    mat = (
        focus.pivot_table(index=focus["ligand"] + "–" + focus["receptor"], columns="library", values="pct_product")
        .reindex(index=pairs, columns=order)
    )
    fig, ax = plt.subplots(figsize=(7.6, max(4.2, 0.32 * len(pairs) + 1.8)))
    im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([s.replace(" ", "\n") for s in order], fontsize=8)
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels(pairs, fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="%L × %R / 100  (naive)")
    style_axes(ax, "Cheap LR potential (not CellChat; no permutation)")
    fig.suptitle("GSE316655 ligand–receptor — CLDN/LILRB + MHC-I/LILRB", fontsize=12)
    add_footer(fig)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(out / "fig7_ligand_receptor.png", dpi=160)
    fig.savefig(out / "fig7_ligand_receptor.pdf")
    plt.close(fig)


def plot_score_heatmap(lib_sum: pd.DataFrame, out: Path) -> None:
    cols = [
        "all_cldn4_mean",
        "myeloid_lilrb2_mean",
        "myeloid_nfkb_mean",
        "myeloid_stat_mean",
        "myeloid_suppress_myeloid_mean",
        "myeloid_mhci_mean",
        "tnk_cytotox_mean",
        "tnk_exhaust_mean",
        "all_mhci_mean",
    ]
    labels = [
        "CLDN4 (all)",
        "LILRB2 (myeloid)",
        "NF-κB (myeloid)",
        "STAT (myeloid)",
        "suppressors (myeloid)",
        "MHC-I (myeloid)",
        "cytotox (T/NK)",
        "exhaust (T/NK)",
        "MHC-I (all)",
    ]
    use = [c for c in cols if c in lib_sum.columns]
    labs = [labels[cols.index(c)] for c in use]
    mat = lib_sum.set_index("short")[use].T
    # z-score rows for display only
    z = mat.sub(mat.mean(axis=1), axis=0)
    sd = mat.std(axis=1).replace(0, np.nan)
    z = z.div(sd, axis=0)
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-1.5, vmax=1.5)
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels([c.replace(" ", "\n") for c in mat.columns], fontsize=8)
    ax.set_yticks(range(len(labs)))
    ax.set_yticklabels(labs, fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="row z (4 libraries)")
    style_axes(ax, "Library-level scores (z across n=4 libraries)")
    fig.suptitle("GSE316655 score overview — do not treat z as statistics", fontsize=12)
    add_footer(fig)
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    fig.savefig(out / "fig8_score_overview.png", dpi=160)
    fig.savefig(out / "fig8_score_overview.pdf")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("/tmp/gse316655"))
    ap.add_argument("--outdir", type=Path, default=Path("."))
    args = ap.parse_args()
    out = args.outdir
    figdir = out / "figures"
    tabdir = out / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    loaded = []
    for spec in LIBRARIES:
        loaded.append(load_library(args.raw, spec))

    cells = pd.concat([x["cells"] for x in loaded], ignore_index=True)
    raw_qc = pd.concat([x["raw_qc"] for x in loaded], ignore_index=True)
    present = set()
    missing = set()
    for x in loaded:
        present.update(x["present_genes"])
        missing.update(x["missing_genes"])

    # library-level summary
    rows = []
    for spec in LIBRARIES:
        sub = cells[cells["short"] == spec["short"]]
        raw = raw_qc[raw_qc["short"] == spec["short"]]
        myeloid = sub[sub["lineage"] == "myeloid"]
        tnk = sub[sub["lineage"].isin(["T", "NK"])]
        row = {
            "gsm": spec["gsm"],
            "library": spec["library"],
            "short": spec["short"],
            "tissue": spec["tissue"],
            "treatment": spec["treatment"],
            "title": spec["title"],
            "n_barcodes": int((raw_qc["short"] == spec["short"]).sum()),
            "n_qc": int(len(sub)),
            "qc_pass_pct": float(raw["pass_qc"].mean() * 100.0),
            "median_human_umi_qc": float(sub["human_umi"].median()) if len(sub) else float("nan"),
            "median_human_frac_all": float(raw["human_frac"].median()),
            "species": "Homo sapiens CD45+ in NSG-SGM3 humanized mouse",
            "tumor_model": "SK-MEL-5 (human melanoma cell line)",
            "host": "NSG-SGM3 + human CB CD34+",
            "alignment": "GRCh38_and_mm10-2020-A",
        }
        row.update(describe_group(sub, "all"))
        row.update(describe_group(myeloid, "myeloid"))
        row.update(describe_group(tnk, "tnk"))
        row["myeloid_n"] = int(len(myeloid))
        row["tnk_n"] = int(len(tnk))
        row["myeloid_frac"] = float(len(myeloid) / len(sub) * 100.0) if len(sub) else float("nan")
        row["tnk_frac"] = float(len(tnk) / len(sub) * 100.0) if len(sub) else float("nan")
        for lin in ["myeloid", "T", "NK", "B", "pDC", "melanoma", "epithelial", "other"]:
            row[f"n_{lin}"] = int((sub["lineage"] == lin).sum())
        rows.append(row)
    lib_sum = pd.DataFrame(rows)

    cldn = claudin_table(cells, present)
    lilr = lilr_table(cells)
    lr = cheap_lr(cells)

    # exploratory cell-level tumor myeloid / TNK contrasts (1 vs 1 library)
    tumor_my_anti = cells[(cells["tissue"] == "tumor") & (cells["treatment"] == "anti-LILRB2") & (cells["lineage"] == "myeloid")]
    tumor_my_iso = cells[(cells["tissue"] == "tumor") & (cells["treatment"] == "isotype") & (cells["lineage"] == "myeloid")]
    tumor_tnk_anti = cells[(cells["tissue"] == "tumor") & (cells["treatment"] == "anti-LILRB2") & (cells["lineage"].isin(["T", "NK"]))]
    tumor_tnk_iso = cells[(cells["tissue"] == "tumor") & (cells["treatment"] == "isotype") & (cells["lineage"].isin(["T", "NK"]))]
    tests = []
    for name, a, b in [
        ("tumor_myeloid_nfkb_anti_vs_iso", tumor_my_anti["nfkb"], tumor_my_iso["nfkb"]),
        ("tumor_myeloid_stat_anti_vs_iso", tumor_my_anti["stat"], tumor_my_iso["stat"]),
        ("tumor_myeloid_lilrb2_anti_vs_iso", tumor_my_anti["lilrb2"], tumor_my_iso["lilrb2"]),
        ("tumor_myeloid_mhci_anti_vs_iso", tumor_my_anti["mhci"], tumor_my_iso["mhci"]),
        ("tumor_myeloid_suppress_anti_vs_iso", tumor_my_anti["suppress_myeloid"], tumor_my_iso["suppress_myeloid"]),
        ("tumor_tnk_cytotox_anti_vs_iso", tumor_tnk_anti["cytotox"], tumor_tnk_iso["cytotox"]),
        ("tumor_tnk_exhaust_anti_vs_iso", tumor_tnk_anti["exhaust"], tumor_tnk_iso["exhaust"]),
        ("tumor_all_cldn4_anti_vs_iso", cells[(cells.tissue == "tumor") & (cells.treatment == "anti-LILRB2")]["cldn4"], cells[(cells.tissue == "tumor") & (cells.treatment == "isotype")]["cldn4"]),
    ]:
        rec = exploratory_mwu(a.to_numpy(), b.to_numpy())
        rec["contrast"] = name
        tests.append(rec)
    tests_df = pd.DataFrame(tests)

    gene_cov = pd.DataFrame(
        {
            "gene": sorted(set(all_panel_genes())),
            "present_human": [g in present for g in sorted(set(all_panel_genes()))],
        }
    )

    n_table = pd.DataFrame(
        [
            {"item": "public_10x_libraries", "n": 4, "note": "GEO MTX/TSV; series also lists these 4 GSM"},
            {"item": "humanized_donors_labeled", "n": 1, "note": "GEO title prefix ND; not independent mice"},
            {"item": "tumor_anti-LILRB2_libraries", "n": 1, "note": "GSM9457798 SK-MEL-5 TIL CD45+"},
            {"item": "tumor_isotype_libraries", "n": 1, "note": "GSM9457799 SK-MEL-5 TIL CD45+"},
            {"item": "PB_anti-LILRB2_libraries", "n": 1, "note": "GSM9457800"},
            {"item": "PB_isotype_libraries", "n": 1, "note": "GSM9457801; missing from some GEO HTML views"},
            {"item": "barcodes_pre_QC", "n": int(len(raw_qc)), "note": "Cell Ranger filtered matrices as deposited"},
            {"item": "cells_post_QC", "n": int(len(cells)), "note": "human UMI≥200, genes≥100, human_frac≥0.80, mito≤0.20"},
            {"item": "tumor_QC_cells", "n": int((cells.tissue == "tumor").sum()), "note": "anti + isotype"},
            {"item": "PB_QC_cells", "n": int((cells.tissue == "PB").sum()), "note": "anti + isotype"},
            {"item": "myeloid_QC_cells", "n": int((cells.lineage == "myeloid").sum()), "note": "marker score; not author annot"},
            {"item": "TNK_QC_cells", "n": int(cells.lineage.isin(["T", "NK"]).sum()), "note": "marker score"},
            {"item": "biological_replicates_per_arm", "n": 1, "note": "cannot support sample-level p-values"},
        ]
    )

    lib_sum.to_csv(tabdir / "library_summary.tsv", sep="\t", index=False)
    cldn.to_csv(tabdir / "claudin_detection.tsv", sep="\t", index=False)
    lilr.to_csv(tabdir / "lilr_by_lineage.tsv", sep="\t", index=False)
    lr.to_csv(tabdir / "ligand_receptor_naive.tsv", sep="\t", index=False)
    tests_df.to_csv(tabdir / "exploratory_cell_mwu.tsv", sep="\t", index=False)
    n_table.to_csv(tabdir / "n_table.tsv", sep="\t", index=False)
    gene_cov.to_csv(tabdir / "gene_coverage.tsv", sep="\t", index=False)
    raw_qc.to_csv(tabdir / "barcode_qc.tsv.gz", sep="\t", index=False)
    # keep a slim cell table (no all umi_* to limit size)
    slim_cols = [
        c
        for c in cells.columns
        if not c.startswith("umi_")
        or c
        in {
            "umi_CLDN4",
            "umi_CLDN18",
            "umi_LILRB2",
            "umi_LILRB1",
            "umi_LILRB5",
            "umi_HLA-A",
            "umi_B2M",
            "umi_CD3D",
            "umi_NKG7",
            "umi_LYZ",
            "umi_CD274",
            "umi_PDCD1",
        }
    ]
    cells[slim_cols].to_csv(tabdir / "cells_qc_slim.tsv.gz", sep="\t", index=False)

    plot_qc(raw_qc, figdir)
    plot_composition(cells, figdir)
    plot_claudins(cldn, figdir)
    plot_lilrb(cells, figdir)
    plot_myeloid_scores(cells, figdir)
    plot_tnk_mhci(cells, figdir)
    plot_lr(lr, figdir)
    plot_score_heatmap(lib_sum, figdir)

    summary = {
        "dataset": "GSE316655",
        "pmid": "41931598",
        "doi": "10.1126/sciimmunol.adt7832",
        "paper": "Liu et al. Sci Immunol 2026",
        "species": "Homo sapiens (FACS CD45+) in Mus musculus NSG-SGM3 host",
        "model": "Humanized NSG-SGM3 + human cord-blood CD34+; SK-MEL-5 subcutaneous tumors",
        "alignment": "Cell Ranger 5.0.0 GRCh38_and_mm10-2020-A",
        "input": "public MTX/TSV only (no SRA re-align)",
        "n_libraries": 4,
        "n_donors_labeled": 1,
        "qc_gates": {
            "min_human_umi": QC_MIN_HUMAN_UMI,
            "min_human_genes": QC_MIN_HUMAN_GENES,
            "min_human_frac": QC_MIN_HUMAN_FRAC,
            "max_mito_frac": QC_MAX_MITO_FRAC,
        },
        "n_barcodes": int(len(raw_qc)),
        "n_qc": int(len(cells)),
        "lineage_counts": cells["lineage"].value_counts().to_dict(),
        "library_n_qc": lib_sum.set_index("short")["n_qc"].to_dict(),
        "library_myeloid_n": lib_sum.set_index("short")["myeloid_n"].to_dict(),
        "library_tnk_n": lib_sum.set_index("short")["tnk_n"].to_dict(),
        "present_panel_n": len(present),
        "missing_panel": sorted(missing),
        "exploratory_tests": tests_df.to_dict(orient="records"),
        "note": (
            "Unit of biological replication is the library (n=1 per tissue×treatment). "
            "Cell-level p-values are exploratory and overstate certainty."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ["n_barcodes", "n_qc", "library_n_qc", "lineage_counts"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
