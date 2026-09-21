#!/usr/bin/env python3
"""CLDN4-related NHEJ / STING / IFN / APM extract from public GSE316655.

Species/model: human CD45+ FACS scRNA from SK-MEL-5 tumors and peripheral
blood in NSG-SGM3 mice humanized with cord-blood CD34+ cells. Cell Ranger
5.0.0, GRCh38_and_mm10-2020-A. Liu et al. Sci Immunol 2026 (PMID 41931598).

The deposited contrast is anti-LILRB2 vs isotype (1 library each). It is not
a CLDN4 perturbation. CLDN4-positive cells are counted, and gene-level
differences vs CLDN4-negative cells in the same compartment are written as a
descriptive extract. Biological n per arm is 1 library.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread

from gene_sets import ANCHORS, LINEAGE, PROGRAMS, lookup_names, panel_symbols

LIBRARIES = [
    {
        "gsm": "GSM9457798",
        "stem": "GSM9457798_ND_anti_B2",
        "tissue": "tumor",
        "treatment": "anti-LILRB2",
        "short": "Tumor anti-LILRB2",
    },
    {
        "gsm": "GSM9457799",
        "stem": "GSM9457799_ND_CTR",
        "tissue": "tumor",
        "treatment": "isotype",
        "short": "Tumor isotype",
    },
    {
        "gsm": "GSM9457800",
        "stem": "GSM9457800_PB_anti_B2",
        "tissue": "PB",
        "treatment": "anti-LILRB2",
        "short": "PB anti-LILRB2",
    },
    {
        "gsm": "GSM9457801",
        "stem": "GSM9457801_PB_CTR",
        "tissue": "PB",
        "treatment": "isotype",
        "short": "PB isotype",
    },
]

MODEL = (
    "Human CD45+ FACS | SK-MEL-5 in NSG-SGM3 humanized mice | "
    "GRCh38+mm10 | GSE316655 | CLDN4 extract is descriptive"
)

QC_MIN_HUMAN_UMI = 200
QC_MIN_HUMAN_GENES = 100
QC_MIN_HUMAN_FRAC = 0.80
QC_MAX_MITO_FRAC = 0.20
LINEAGE_MIN = 0.20
LINEAGE_MARGIN = 0.05

SHORT_ORDER = [x["short"] for x in LIBRARIES]
PROGRAM_ORDER = ["NHEJ", "STING", "IFN", "APM"]
PROGRAM_COLORS = {
    "NHEJ": "#4C78A8",
    "STING": "#F58518",
    "IFN": "#54A24B",
    "APM": "#E45756",
}


def strip_name(raw: str) -> tuple[str, str]:
    if raw.startswith("GRCh38_"):
        return "human", raw[len("GRCh38_") :]
    if raw.startswith("mm10___"):
        return "mouse", raw[len("mm10___") :]
    if raw.startswith("mm10_"):
        return "mouse", raw[len("mm10_") :]
    return "unknown", raw


def read_features(path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            species, symbol = strip_name(parts[1] if len(parts) > 1 else parts[0])
            rows.append((i, species, symbol))
    return pd.DataFrame(rows, columns=["idx", "species", "symbol"])


def read_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [line.strip() for line in fh if line.strip()]


def first_index(feat: pd.DataFrame, species: str, names: list[str]) -> tuple[int | None, str | None]:
    sub = feat.loc[feat["species"] == species]
    for name in names:
        hit = sub.loc[sub["symbol"] == name]
        if len(hit):
            return int(hit.iloc[0]["idx"]), name
    return None, None


def resolve_human(feat: pd.DataFrame, symbol: str) -> tuple[int | None, str | None]:
    return first_index(feat, "human", lookup_names(symbol))


def mean_score(log_mat: np.ndarray, gene_idx: list[int]) -> np.ndarray:
    if not gene_idx:
        return np.full(log_mat.shape[1], np.nan, dtype=np.float32)
    return log_mat[gene_idx].mean(axis=0).astype(np.float32)


def assign_lineage(log_mat: np.ndarray, marker_idx: dict[str, list[int]]) -> np.ndarray:
    names = list(LINEAGE)
    scores = np.vstack([mean_score(log_mat, marker_idx[k]) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    weak = (best_val < LINEAGE_MIN) | ((best_val - second) < LINEAGE_MARGIN)
    assigned[weak] = "other"
    return assigned


def load_library(raw: Path, spec: dict) -> tuple[pd.DataFrame, list[dict]]:
    feat = read_features(raw / f"{spec['stem']}_features.tsv.gz")
    barcodes = read_barcodes(raw / f"{spec['stem']}_barcodes.tsv.gz")
    mtx_path = raw / f"{spec['stem']}_matrix.mtx.gz"
    print(f"reading {mtx_path.name}", flush=True)
    mat = mmread(mtx_path, spmatrix=True).tocsr().astype(np.float32)
    if mat.shape[0] != len(feat) or mat.shape[1] != len(barcodes):
        raise SystemExit(
            f"{spec['gsm']} shape {mat.shape} != features {len(feat)} x barcodes {len(barcodes)}"
        )

    human = (feat["species"] == "human").to_numpy()
    mouse = (feat["species"] == "mouse").to_numpy()
    mito = human & feat["symbol"].str.startswith("MT-").to_numpy()
    human_umi = np.asarray(mat[human].sum(axis=0)).ravel()
    mouse_umi = np.asarray(mat[mouse].sum(axis=0)).ravel()
    mito_umi = np.asarray(mat[mito].sum(axis=0)).ravel()
    total = human_umi + mouse_umi
    human_genes = np.asarray((mat[human] > 0).sum(axis=0)).ravel()
    human_frac = np.divide(human_umi, total, out=np.zeros_like(human_umi), where=total > 0)
    mito_frac = np.divide(mito_umi, human_umi, out=np.zeros_like(mito_umi), where=human_umi > 0)
    keep = (
        (human_umi >= QC_MIN_HUMAN_UMI)
        & (human_genes >= QC_MIN_HUMAN_GENES)
        & (human_frac >= QC_MIN_HUMAN_FRAC)
        & (mito_frac <= QC_MAX_MITO_FRAC)
    )

    wanted = panel_symbols()
    resolved: list[dict] = []
    row_ids: list[int] = []
    for symbol in wanted:
        idx, used = resolve_human(feat, symbol)
        resolved.append(
            {
                "symbol": symbol,
                "matrix_symbol": used if used else "",
                "present": idx is not None,
                "gsm": spec["gsm"],
            }
        )
        if idx is not None:
            row_ids.append(idx)

    # Mouse Cldn4 is recorded only as a contaminant count, not as the contrast.
    mouse_cldn4_idx, mouse_cldn4_name = first_index(feat, "mouse", ["Cldn4"])
    mouse_cldn4 = (
        np.asarray(mat[mouse_cldn4_idx].todense()).ravel()
        if mouse_cldn4_idx is not None
        else np.zeros(mat.shape[1], dtype=np.float32)
    )

    sub = np.asarray(mat[row_ids][:, keep].todense())
    lib = human_umi[keep]
    cp = np.where(lib > 0, sub / lib * 1e4, 0.0)
    log_mat = np.log1p(cp).astype(np.float32)

    symbol_row = {}
    cursor = 0
    for rec in resolved:
        if rec["present"]:
            symbol_row[rec["symbol"]] = cursor
            cursor += 1

    marker_idx = {
        name: [symbol_row[g] for g in genes if g in symbol_row] for name, genes in LINEAGE.items()
    }
    lineage = assign_lineage(log_mat, marker_idx)

    data = {
        "barcode": np.array(barcodes, dtype=object)[keep],
        "gsm": spec["gsm"],
        "library": spec["short"],
        "tissue": spec["tissue"],
        "treatment": spec["treatment"],
        "human_umi": lib.astype(np.int32),
        "human_genes": human_genes[keep].astype(np.int32),
        "human_frac": human_frac[keep].astype(np.float32),
        "mito_frac": mito_frac[keep].astype(np.float32),
        "mouse_cldn4_umi": mouse_cldn4[keep].astype(np.float32),
        "mouse_cldn4_symbol": mouse_cldn4_name or "",
        "lineage": lineage,
    }
    for symbol, row in symbol_row.items():
        data[f"{symbol}__umi"] = sub[row].astype(np.float32)
        data[f"{symbol}__log"] = log_mat[row]
    for program, genes in PROGRAMS.items():
        idxs = [symbol_row[g] for g in genes if g in symbol_row]
        data[f"score_{program}"] = mean_score(log_mat, idxs)
        data[f"n_genes_{program}"] = len(idxs)
    frame = pd.DataFrame(data)

    print(
        f"  {spec['short']}: barcodes {len(barcodes)} QC {int(keep.sum())} "
        f"CLDN4>0 {int((frame['CLDN4__umi'] > 0).sum()) if 'CLDN4__umi' in frame else 'NA'}",
        flush=True,
    )
    del mat
    return frame, resolved


def coverage_table(resolved_rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(resolved_rows)
    # Presence is a property of the reference; require agreement across libraries.
    agg = (
        df.groupby("symbol", sort=False)
        .agg(
            present_libraries=("present", "sum"),
            matrix_symbol=("matrix_symbol", lambda s: next((x for x in s if x), "")),
        )
        .reset_index()
    )
    program = {}
    for name, genes in PROGRAMS.items():
        for g in genes:
            program[g] = name
    for g in ANCHORS:
        program.setdefault(g, "anchor")
    for genes in LINEAGE.values():
        for g in genes:
            program.setdefault(g, "lineage")
    agg["role"] = agg["symbol"].map(program)
    agg["present"] = agg["present_libraries"] == agg["present_libraries"].max()
    return agg


def _side_stats(values: np.ndarray) -> dict:
    if len(values) == 0:
        return {
            "n": 0,
            "n_detected": 0,
            "pct_detected": np.nan,
            "mean": np.nan,
            "median": np.nan,
        }
    detected = values > 0
    return {
        "n": int(len(values)),
        "n_detected": int(detected.sum()),
        "pct_detected": float(100.0 * detected.mean()),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
    }


def contrast_block(df: pd.DataFrame, compartment: str, pos: np.ndarray, neg: np.ndarray) -> list[dict]:
    rows = []
    for program, genes in PROGRAMS.items():
        for gene in genes:
            col = f"{gene}__log"
            if col not in df.columns:
                rows.append(
                    {
                        "compartment": compartment,
                        "program": program,
                        "gene": gene,
                        "in_matrix": False,
                        "n_pos": int(pos.sum()),
                        "n_neg": int(neg.sum()),
                        "n_pos_detected": 0,
                        "n_neg_detected": 0,
                        "pct_pos": np.nan,
                        "pct_neg": np.nan,
                        "mean_pos": np.nan,
                        "mean_neg": np.nan,
                        "median_pos": np.nan,
                        "median_neg": np.nan,
                        "delta_mean": np.nan,
                        "delta_median": np.nan,
                    }
                )
                continue
            a = _side_stats(df.loc[pos, col].to_numpy())
            b = _side_stats(df.loc[neg, col].to_numpy())
            rows.append(
                {
                    "compartment": compartment,
                    "program": program,
                    "gene": gene,
                    "in_matrix": True,
                    "n_pos": a["n"],
                    "n_neg": b["n"],
                    "n_pos_detected": a["n_detected"],
                    "n_neg_detected": b["n_detected"],
                    "pct_pos": a["pct_detected"],
                    "pct_neg": b["pct_detected"],
                    "mean_pos": a["mean"],
                    "mean_neg": b["mean"],
                    "median_pos": a["median"],
                    "median_neg": b["median"],
                    "delta_mean": a["mean"] - b["mean"] if a["n"] and b["n"] else np.nan,
                    "delta_median": a["median"] - b["median"] if a["n"] and b["n"] else np.nan,
                }
            )
    return rows


def program_contrast(df: pd.DataFrame, compartment: str, pos: np.ndarray, neg: np.ndarray) -> list[dict]:
    rows = []
    for program in PROGRAM_ORDER:
        col = f"score_{program}"
        a = _side_stats(df.loc[pos, col].to_numpy()) if pos.any() else _side_stats(np.array([]))
        b = _side_stats(df.loc[neg, col].to_numpy()) if neg.any() else _side_stats(np.array([]))
        # Scores are means, so "detected" above is not gene detection. Recompute
        # a real detection rate: fraction of cells with score > 0.
        rows.append(
            {
                "compartment": compartment,
                "program": program,
                "n_genes": int(df[f"n_genes_{program}"].iloc[0]) if len(df) else 0,
                "n_pos": int(pos.sum()),
                "n_neg": int(neg.sum()),
                "mean_pos": a["mean"],
                "mean_neg": b["mean"],
                "median_pos": a["median"],
                "median_neg": b["median"],
                "delta_mean": a["mean"] - b["mean"] if a["n"] and b["n"] else np.nan,
                "delta_median": a["median"] - b["median"] if a["n"] and b["n"] else np.nan,
                "pct_score_gt0_pos": a["pct_detected"],
                "pct_score_gt0_neg": b["pct_detected"],
            }
        )
    return rows


def build_contrasts(cells: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pos = cells["CLDN4__umi"] > 0
    gene_rows: list[dict] = []
    prog_rows: list[dict] = []

    def add(name: str, mask: pd.Series) -> None:
        m = mask.to_numpy()
        p = pos.to_numpy() & m
        n = (~pos).to_numpy() & m
        if m.sum() == 0:
            return
        gene_rows.extend(contrast_block(cells, name, p, n))
        prog_rows.extend(program_contrast(cells, name, p, n))

    add("QC_all", pd.Series(True, index=cells.index))
    add("tumor_QC", cells["tissue"] == "tumor")
    add("PB_QC", cells["tissue"] == "PB")
    for tissue in ("tumor", "PB"):
        for lineage in list(LINEAGE) + ["other"]:
            add(
                f"{tissue}_{lineage}",
                (cells["tissue"] == tissue) & (cells["lineage"] == lineage),
            )
    for lineage in sorted(cells.loc[pos, "lineage"].unique()):
        add(f"lineage_{lineage}", cells["lineage"] == lineage)
    return pd.DataFrame(gene_rows), pd.DataFrame(prog_rows)


def treatment_program_table(cells: pd.DataFrame) -> pd.DataFrame:
    """Anti-LILRB2 vs isotype inside tumor. n_library = 1 vs 1."""
    rows = []
    tumor = cells[cells["tissue"] == "tumor"]
    groups = [("tumor_QC", tumor)]
    for lineage in list(LINEAGE) + ["other"]:
        groups.append((f"tumor_{lineage}", tumor[tumor["lineage"] == lineage]))
    for name, sub in groups:
        if sub.empty:
            continue
        anti = sub["treatment"] == "anti-LILRB2"
        iso = sub["treatment"] == "isotype"
        for program in PROGRAM_ORDER:
            col = f"score_{program}"
            a = sub.loc[anti, col].to_numpy()
            b = sub.loc[iso, col].to_numpy()
            rows.append(
                {
                    "compartment": name,
                    "program": program,
                    "n_anti": int(anti.sum()),
                    "n_iso": int(iso.sum()),
                    "n_libraries_per_arm": 1,
                    "mean_anti": float(np.mean(a)) if len(a) else np.nan,
                    "mean_iso": float(np.mean(b)) if len(b) else np.nan,
                    "median_anti": float(np.median(a)) if len(a) else np.nan,
                    "median_iso": float(np.median(b)) if len(b) else np.nan,
                    "delta_mean_anti_minus_iso": (
                        float(np.mean(a) - np.mean(b)) if len(a) and len(b) else np.nan
                    ),
                    "note": "one library per arm; cell-level contrast is not a replicate test",
                }
            )
    return pd.DataFrame(rows)


def cldn4_cell_table(cells: pd.DataFrame) -> pd.DataFrame:
    pos = cells[cells["CLDN4__umi"] > 0].copy()
    if pos.empty:
        return pos
    for program in PROGRAM_ORDER:
        col = f"score_{program}"
        pct = []
        for _, row in pos.iterrows():
            pool = cells[
                (cells["library"] == row["library"]) & (cells["lineage"] == row["lineage"])
            ][col].to_numpy()
            # Percent of the matched pool at or below this cell.
            pct.append(float(stats.percentileofscore(pool, row[col], kind="mean")))
        pos[f"pctile_{program}_in_library_lineage"] = pct
    keep = [
        "barcode",
        "gsm",
        "library",
        "tissue",
        "treatment",
        "lineage",
        "human_umi",
        "CLDN4__umi",
        "CLDN18__umi",
        "LILRB2__umi",
        "mouse_cldn4_umi",
    ]
    for program in PROGRAM_ORDER:
        keep.append(f"score_{program}")
        keep.append(f"pctile_{program}_in_library_lineage")
    return pos[keep].sort_values(["tissue", "treatment", "lineage", "barcode"])


def cldn4_cell_genes(cells: pd.DataFrame) -> pd.DataFrame:
    pos = cells[cells["CLDN4__umi"] > 0]
    rows = []
    for _, cell in pos.iterrows():
        for program, genes in PROGRAMS.items():
            for gene in genes:
                umi_col = f"{gene}__umi"
                log_col = f"{gene}__log"
                rows.append(
                    {
                        "barcode": cell["barcode"],
                        "gsm": cell["gsm"],
                        "library": cell["library"],
                        "lineage": cell["lineage"],
                        "program": program,
                        "gene": gene,
                        "umi": float(cell[umi_col]) if umi_col in cells.columns else np.nan,
                        "log1p_cp10k": float(cell[log_col]) if log_col in cells.columns else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def detection_by_library(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for spec in LIBRARIES:
        sub = cells[cells["library"] == spec["short"]]
        for program, genes in PROGRAMS.items():
            for gene in genes:
                col = f"{gene}__umi"
                if col not in sub.columns or sub.empty:
                    pct = np.nan
                    mean = np.nan
                    n_det = 0
                else:
                    umi = sub[col].to_numpy()
                    n_det = int((umi > 0).sum())
                    pct = float(100.0 * (umi > 0).mean())
                    mean = float(sub[f"{gene}__log"].mean())
                rows.append(
                    {
                        "library": spec["short"],
                        "tissue": spec["tissue"],
                        "treatment": spec["treatment"],
                        "n_cells": int(len(sub)),
                        "program": program,
                        "gene": gene,
                        "n_detected": n_det,
                        "pct_detected": pct,
                        "mean_log1p_cp10k": mean,
                    }
                )
    return pd.DataFrame(rows)


def library_summary(cells: pd.DataFrame, n_barcodes: dict[str, int]) -> pd.DataFrame:
    rows = []
    for spec in LIBRARIES:
        sub = cells[cells["library"] == spec["short"]]
        row = {
            "gsm": spec["gsm"],
            "library": spec["short"],
            "tissue": spec["tissue"],
            "treatment": spec["treatment"],
            "n_barcodes": n_barcodes[spec["gsm"]],
            "n_qc": int(len(sub)),
            "n_cldn4_pos": int((sub["CLDN4__umi"] > 0).sum()) if len(sub) else 0,
            "cldn4_umi_sum": float(sub["CLDN4__umi"].sum()) if len(sub) else 0,
            "n_cldn18_pos": int((sub["CLDN18__umi"] > 0).sum()) if len(sub) else 0,
            "n_mouse_cldn4_pos": int((sub["mouse_cldn4_umi"] > 0).sum()) if len(sub) else 0,
        }
        for program in PROGRAM_ORDER:
            row[f"mean_{program}"] = float(sub[f"score_{program}"].mean()) if len(sub) else np.nan
        counts = sub["lineage"].value_counts()
        for lineage in list(LINEAGE) + ["other"]:
            row[f"n_{lineage}"] = int(counts.get(lineage, 0))
        rows.append(row)
    return pd.DataFrame(rows)


def depth_check(cells: pd.DataFrame) -> pd.DataFrame:
    """CLDN4 UMI>0 cells are deeper on average; NHEJ detection tracks depth."""
    masks = {
        "QC_all": pd.Series(True, index=cells.index),
        "lineage_T": cells["lineage"] == "T",
        "lineage_B": cells["lineage"] == "B",
        "lineage_melanoma": cells["lineage"] == "melanoma",
        "tumor_melanoma": (cells["tissue"] == "tumor") & (cells["lineage"] == "melanoma"),
    }
    rows = []
    for name, mask in masks.items():
        sub = cells.loc[mask]
        pos = sub["CLDN4__umi"] > 0
        neg = ~pos
        row = {
            "compartment": name,
            "n_pos": int(pos.sum()),
            "n_neg": int(neg.sum()),
            "median_human_umi_pos": float(sub.loc[pos, "human_umi"].median()) if pos.any() else np.nan,
            "median_human_umi_neg": float(sub.loc[neg, "human_umi"].median()) if neg.any() else np.nan,
        }
        if neg.sum() >= 20:
            depth = np.log1p(sub.loc[neg, "human_umi"].to_numpy())
            for program in PROGRAM_ORDER:
                rho = stats.spearmanr(depth, sub.loc[neg, f"score_{program}"].to_numpy()).statistic
                row[f"spearman_logumi_vs_{program}_in_cldn4neg"] = float(rho)
        rows.append(row)
    return pd.DataFrame(rows)


def add_footer(fig) -> None:
    fig.text(0.005, 0.004, MODEL, fontsize=7, color="#333333", ha="left", va="bottom")


def plot_detection(det: pd.DataFrame, out: Path) -> None:
    genes = []
    for program in PROGRAM_ORDER:
        genes.extend(PROGRAMS[program])
    mat = np.full((len(genes), len(SHORT_ORDER)), np.nan)
    for i, gene in enumerate(genes):
        for j, lib in enumerate(SHORT_ORDER):
            hit = det[(det["gene"] == gene) & (det["library"] == lib)]
            if len(hit):
                mat[i, j] = hit.iloc[0]["pct_detected"]
    fig_h = max(8.5, 0.18 * len(genes) + 1.6)
    fig, ax = plt.subplots(figsize=(8.2, fig_h))
    im = ax.imshow(mat, aspect="auto", cmap="YlGnBu", vmin=0, vmax=max(5, np.nanpercentile(mat, 98)))
    ax.set_xticks(range(len(SHORT_ORDER)))
    ax.set_xticklabels(SHORT_ORDER, rotation=30, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=6)
    # Program brackets via a colored y-tick.
    for tick, gene in zip(ax.get_yticklabels(), genes):
        for program, members in PROGRAMS.items():
            if gene in members:
                tick.set_color(PROGRAM_COLORS[program])
    ax.set_title("Percent of QC cells with UMI > 0")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("% cells UMI>0")
    fig.tight_layout(rect=(0, 0.012, 1, 1))
    add_footer(fig)
    fig.savefig(out / "fig1_detection_heatmap.png", dpi=160)
    fig.savefig(out / "fig1_detection_heatmap.pdf")
    plt.close(fig)


def plot_modules(cells: pd.DataFrame, out: Path) -> None:
    tumor = cells[cells["tissue"] == "tumor"].copy()
    order = [x for x in list(LINEAGE) + ["other"] if (tumor["lineage"] == x).sum() >= 15]
    fig, axes = plt.subplots(1, 4, figsize=(12.4, 4.4), sharey=False)
    rng = np.random.default_rng(0)
    for ax, program in zip(axes, PROGRAM_ORDER):
        col = f"score_{program}"
        data = [tumor.loc[tumor["lineage"] == lin, col].to_numpy() for lin in order]
        bp = ax.boxplot(
            data,
            tick_labels=order,
            showfliers=False,
            patch_artist=True,
            medianprops={"color": "black"},
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(PROGRAM_COLORS[program])
            patch.set_alpha(0.35)
        pos = tumor[tumor["CLDN4__umi"] > 0]
        for _, cell in pos.iterrows():
            if cell["lineage"] not in order:
                continue
            x = order.index(cell["lineage"]) + 1
            ax.scatter(
                x + rng.normal(0, 0.04),
                cell[col],
                s=28,
                c="black",
                zorder=3,
            )
        ax.set_title(program, color=PROGRAM_COLORS[program])
        ax.tick_params(axis="x", labelrotation=40)
        ax.set_ylabel("mean log1p(CP10k)" if program == "NHEJ" else "")
    fig.suptitle("Tumor QC module scores; black points are CLDN4 UMI>0 cells", fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 0.92))
    add_footer(fig)
    fig.savefig(out / "fig2_tumor_module_scores.png", dpi=160)
    fig.savefig(out / "fig2_tumor_module_scores.pdf")
    plt.close(fig)


def plot_cldn4_cells(cell_tbl: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    if cell_tbl.empty:
        ax.text(0.5, 0.5, "No CLDN4 UMI>0 QC cells", ha="center")
    else:
        labels = [
            f"{r.library} | {r.lineage} | UMI {int(r.CLDN4__umi)}"
            for r in cell_tbl.itertuples()
        ]
        y = np.arange(len(cell_tbl))
        for k, program in enumerate(PROGRAM_ORDER):
            vals = cell_tbl[f"pctile_{program}_in_library_lineage"].to_numpy()
            ax.scatter(
                vals,
                y + (k - 1.5) * 0.12,
                s=36,
                color=PROGRAM_COLORS[program],
                label=program,
                zorder=3,
            )
        ax.axvline(50, color="#bbbbbb", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Percentile of module score inside the same library and lineage")
        ax.set_xlim(0, 100)
        ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.16))
    ax.set_title("Where the CLDN4-positive cells sit")
    fig.tight_layout(rect=(0, 0.05, 1, 0.90))
    add_footer(fig)
    fig.savefig(out / "fig3_cldn4_cell_percentiles.png", dpi=160)
    fig.savefig(out / "fig3_cldn4_cell_percentiles.pdf")
    plt.close(fig)


def plot_program_delta(prog: pd.DataFrame, out: Path) -> None:
    # Compartments that actually contain a CLDN4-positive cell, plus tumor melanoma.
    use = []
    for name in ["QC_all", "tumor_QC", "tumor_melanoma", "lineage_T", "lineage_B", "lineage_melanoma", "lineage_other"]:
        sub = prog[prog["compartment"] == name]
        if sub.empty:
            continue
        if int(sub["n_pos"].iloc[0]) == 0 and name != "tumor_melanoma":
            continue
        use.append(name)
    # Fall back to whatever has n_pos > 0.
    if not use:
        use = sorted(prog.loc[prog["n_pos"] > 0, "compartment"].unique())
    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    x = np.arange(len(use))
    width = 0.18
    for k, program in enumerate(PROGRAM_ORDER):
        deltas = []
        for name in use:
            hit = prog[(prog["compartment"] == name) & (prog["program"] == program)]
            deltas.append(float(hit.iloc[0]["delta_median"]) if len(hit) else np.nan)
        ax.bar(
            x + (k - 1.5) * width,
            deltas,
            width=width,
            color=PROGRAM_COLORS[program],
            label=program,
        )
    ax.axhline(0, color="black", lw=0.6)
    labels = []
    for name in use:
        n_pos = int(prog.loc[prog["compartment"] == name, "n_pos"].iloc[0])
        labels.append(f"{name}\nn_pos={n_pos}")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Median score, CLDN4 UMI>0 minus CLDN4 UMI=0")
    ax.set_title("Descriptive CLDN4 split (each positive cell has 1 UMI)")
    ax.legend(frameon=False, ncol=4, loc="upper right")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    add_footer(fig)
    fig.savefig(out / "fig4_cldn4_program_delta.png", dpi=160)
    fig.savefig(out / "fig4_cldn4_program_delta.pdf")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("data/raw"))
    ap.add_argument("--outdir", type=Path, default=Path("."))
    args = ap.parse_args()
    tables = args.outdir / "tables"
    figures = args.outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    frames = []
    resolved_all: list[dict] = []
    n_barcodes = {}
    for spec in LIBRARIES:
        frame, resolved = load_library(args.raw, spec)
        frames.append(frame)
        resolved_all.extend(resolved)
        n_barcodes[spec["gsm"]] = len(read_barcodes(args.raw / f"{spec['stem']}_barcodes.tsv.gz"))
    cells = pd.concat(frames, ignore_index=True)

    cov = coverage_table(resolved_all)
    det = detection_by_library(cells)
    gene_de, prog_de = build_contrasts(cells)
    treat = treatment_program_table(cells)
    cell_tbl = cldn4_cell_table(cells)
    cell_genes = cldn4_cell_genes(cells)
    lib_sum = library_summary(cells, n_barcodes)
    depth = depth_check(cells)

    cov.to_csv(tables / "gene_coverage.tsv", sep="\t", index=False)
    det.to_csv(tables / "detection_by_library.tsv", sep="\t", index=False)
    gene_de.to_csv(tables / "cldn4_gene_delta.tsv", sep="\t", index=False)
    prog_de.to_csv(tables / "cldn4_program_delta.tsv", sep="\t", index=False)
    treat.to_csv(tables / "treatment_program_tumor.tsv", sep="\t", index=False)
    cell_tbl.to_csv(tables / "cldn4_positive_cells.tsv", sep="\t", index=False)
    cell_genes.to_csv(tables / "cldn4_positive_cell_genes.tsv", sep="\t", index=False)
    lib_sum.to_csv(tables / "library_summary.tsv", sep="\t", index=False)
    depth.to_csv(tables / "depth_check.tsv", sep="\t", index=False)

    # Slim per-library lineage counts for the write-up.
    lineage_counts = (
        cells.groupby(["library", "tissue", "treatment", "lineage"], observed=True)
        .size()
        .reset_index(name="n")
    )
    lineage_counts.to_csv(tables / "lineage_counts.tsv", sep="\t", index=False)

    plot_detection(det, figures)
    plot_modules(cells, figures)
    plot_cldn4_cells(cell_tbl, figures)
    plot_program_delta(prog_de, figures)

    pos = cells[cells["CLDN4__umi"] > 0]
    summary = {
        "n_barcodes": int(sum(n_barcodes.values())),
        "n_qc": int(len(cells)),
        "n_cldn4_pos": int(len(pos)),
        "cldn4_umi_values": sorted(pos["CLDN4__umi"].astype(int).tolist()),
        "n_cldn18_pos": int((cells["CLDN18__umi"] > 0).sum()),
        "n_mouse_cldn4_pos_in_qc": int((cells["mouse_cldn4_umi"] > 0).sum()),
        "cldn4_pos_by_library": pos.groupby("library").size().astype(int).to_dict() if len(pos) else {},
        "cldn4_pos_by_lineage": pos.groupby("lineage").size().astype(int).to_dict() if len(pos) else {},
        "genes_missing": cov.loc[~cov["present"], "symbol"].tolist(),
        "matrix_symbol_overrides": cov.loc[
            cov["matrix_symbol"].ne("") & cov["matrix_symbol"].ne(cov["symbol"]),
            ["symbol", "matrix_symbol"],
        ].to_dict(orient="records"),
        "program_n_genes": {p: int(cells[f"n_genes_{p}"].iloc[0]) for p in PROGRAM_ORDER},
        "qc": {
            "min_human_umi": QC_MIN_HUMAN_UMI,
            "min_human_genes": QC_MIN_HUMAN_GENES,
            "min_human_frac": QC_MIN_HUMAN_FRAC,
            "max_mito_frac": QC_MAX_MITO_FRAC,
        },
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
