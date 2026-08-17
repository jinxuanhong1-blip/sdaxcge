#!/usr/bin/env python3
"""CellChat-style focused LR panel. Runs only if a leftover merge/series triggers.

Trigger: Spearman ρ<0 with p<0.05, or a clear (non-thin) Q4 vs Q1 drop.
Not CellChat R. Public TISCH values. Focused barrier / checkpoint / chemokine pairs
used in the GSE207422 CLDN4 CellChat slice.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import DATASETS, EPITHELIAL_NONMALIGNANT, MALIGNANT, TNK_LINEAGES
from lib_io import keep_tumor_unit, load_tisch_genes, pick_tissue_col, pick_unit_columns

FOCUSED_PAIRS = [
    ("NECTIN2", "TIGIT", "out", "barrier/inhibitory"),
    ("PVR", "TIGIT", "out", "inhibitory"),
    ("CD274", "PDCD1", "out", "checkpoint"),
    ("CDH1", "ITGAE", "out", "barrier"),
    ("F11R", "ITGAL", "out", "barrier (JAM1)"),
    ("HLA-E", "CD8A", "out", "inhibitory"),
    ("CXCL16", "CXCR6", "out", "recruit"),
    ("CXCL9", "CXCR3", "out", "recruit"),
    ("CXCL10", "CXCR3", "out", "recruit"),
    ("HLA-DRA", "CD4", "out", "MHC-II"),
    ("IFNG", "IFNGR1", "in", "attack into epi"),
    ("FASLG", "FAS", "in", "attack"),
    ("GZMB", "PGRMC1", "in", "cytotoxicity proxy"),
]


def _norm(x) -> str:
    return str(x).strip()


def run_cellchat_if_triggered(trigger_rows: list[dict], units: pd.DataFrame, data_dir: Path, out_dir: Path) -> dict:
    """Score focused LR pairs on the largest triggered leftover merge."""
    if not trigger_rows:
        rec = {"ran": False, "reason": "no trigger"}
        (out_dir / "cellchat_skip.json").write_text(json.dumps(rec, indent=2) + "\n")
        return rec

    pick = sorted(trigger_rows, key=lambda r: (-int(r.get("n") or 0), r.get("p_mean") or 1))[0]
    gses = str(pick["cohorts"]).split("+")
    unit_kind = pick.get("unit_kind") or "patient"
    genes = sorted({g for pair in FOCUSED_PAIRS for g in pair[:2]} | {"CLDN4"})

    ds_by_gse = {v["gse"]: k for k, v in DATASETS.items()}
    sender_high, sender_low, recv = [], [], []
    n_cells = {"high": 0, "low": 0, "tnk": 0}
    used = []
    for gse in gses:
        ds = ds_by_gse.get(gse)
        if not ds:
            continue
        h5 = data_dir / f"{ds}_expression.h5"
        meta_path = data_dir / f"{ds}_CellMetainfo_table.tsv"
        if not h5.exists() or not meta_path.exists():
            continue
        expr, present, _ = load_tisch_genes(h5, genes)
        meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
        cell_col = meta.columns[0]
        meta = meta.set_index(cell_col)
        meta.index = meta.index.astype(str)
        lineage_col = next((c for c in meta.columns if "major-lineage" in c.lower()), None)
        if lineage_col is None:
            continue
        meta["_lineage"] = meta[lineage_col].map(_norm)
        common = meta.index.intersection(expr.index)
        meta = meta.loc[common]
        expr = expr.loc[common]
        patient_col, sample_col = pick_unit_columns(meta)
        tissue_col = pick_tissue_col(meta)
        sample_for_filter = meta[sample_col] if sample_col else (
            meta[patient_col] if patient_col else pd.Series([""] * len(meta), index=meta.index)
        )
        tissue_for_filter = meta[tissue_col] if tissue_col else pd.Series([None] * len(meta), index=meta.index)
        tumor_mask = [keep_tumor_unit(s, t) for s, t in zip(sample_for_filter.astype(str), tissue_for_filter)]
        if sum(tumor_mask) >= 50:
            meta = meta.loc[pd.Series(tumor_mask, index=meta.index)]
            expr = expr.loc[meta.index]
        is_mal = meta["_lineage"].isin(MALIGNANT)
        is_epi = meta["_lineage"].isin(MALIGNANT | EPITHELIAL_NONMALIGNANT)
        is_tnk = meta["_lineage"].isin(TNK_LINEAGES)
        sender_mask = is_mal if int(is_mal.sum()) >= 20 else is_epi
        if "CLDN4" not in expr.columns:
            continue
        cldn = expr.loc[sender_mask, "CLDN4"]
        if len(cldn) < 20:
            continue
        cut = float(cldn.median())
        high_idx = cldn.index[cldn > cut]
        low_idx = cldn.index[cldn <= cut]
        tnk_idx = meta.index[is_tnk]
        sender_high.append(expr.loc[high_idx])
        sender_low.append(expr.loc[low_idx])
        recv.append(expr.loc[tnk_idx])
        n_cells["high"] += int(len(high_idx))
        n_cells["low"] += int(len(low_idx))
        n_cells["tnk"] += int(len(tnk_idx))
        used.append(gse)

    if not sender_high:
        rec = {"ran": False, "reason": "trigger fired but leftover h5 could not be scored"}
        (out_dir / "cellchat_skip.json").write_text(json.dumps(rec, indent=2) + "\n")
        return rec

    high = pd.concat(sender_high)
    low = pd.concat(sender_low)
    tnk = pd.concat(recv)

    rows = []
    for lig, recp, direction, note in FOCUSED_PAIRS:
        if direction == "out":
            if lig not in high.columns or recp not in tnk.columns:
                rows.append({"pair": f"{lig}-{recp}", "dir": direction, "delta": np.nan, "note": "gene missing", "family": note})
                continue
            p_high = float(high[lig].mean()) * float(tnk[recp].mean())
            p_low = float(low[lig].mean()) * float(tnk[recp].mean())
        else:
            if lig not in tnk.columns or recp not in high.columns:
                rows.append({"pair": f"{lig}-{recp}", "dir": direction, "delta": np.nan, "note": "gene missing", "family": note})
                continue
            p_high = float(tnk[lig].mean()) * float(high[recp].mean())
            p_low = float(tnk[lig].mean()) * float(low[recp].mean())
        rows.append(
            {
                "pair": f"{lig}-{recp}",
                "dir": direction,
                "score_high": p_high,
                "score_low": p_low,
                "delta": p_high - p_low,
                "note": note,
                "family": note,
            }
        )
    lr = pd.DataFrame(rows)
    tables = out_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    lr.to_csv(tables / "cellchat_focused_lr.tsv", sep="\t", index=False)

    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    plot = lr.dropna(subset=["delta"]).sort_values("delta")
    if len(plot):
        fig, ax = plt.subplots(figsize=(7.2, max(3.4, 0.32 * len(plot) + 1.2)))
        colors = ["#8b1e3f" if d > 0 else "#1f4e79" for d in plot["delta"]]
        y = np.arange(len(plot))
        ax.barh(y, plot["delta"], color=colors)
        ax.axvline(0, color="0.3", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(plot["pair"], fontsize=8)
        ax.set_xlabel("score(CLDN4-high) − score(CLDN4-low)")
        ax.set_title(
            f"CellChat-style focused LR · {pick['cohorts']}\n"
            f"high n={n_cells['high']:,}  low n={n_cells['low']:,}  T/NK n={n_cells['tnk']:,}"
        )
        fig.tight_layout()
        fig.savefig(fig_dir / "fig_cellchat_focused_lr.png", dpi=160)
        fig.savefig(fig_dir / "fig_cellchat_focused_lr.pdf")
        plt.close(fig)

    rec = {
        "ran": True,
        "reason": f"triggered by {pick['family']} {pick['cohorts']} n={pick['n']}",
        "cohorts": pick["cohorts"],
        "unit_kind": unit_kind,
        "n_units": int(pick["n"]),
        "n_cells": n_cells,
        "series_used": used,
        "n_pairs": int(lr["delta"].notna().sum()),
        "method": "focused LR product of TISCH means; not CellChat R; no permutation",
    }
    rec["finding_md"] = (
        f"CellChat-style **ran** on leftover merge `{pick['cohorts']}` "
        f"({unit_kind} n={pick['n']}; trigger ρ={pick.get('rho_mean')} p={pick.get('p_mean')}). "
        f"Focused LR panel (not CellChat R). Sender cells: CLDN4-high {n_cells['high']:,} vs "
        f"low {n_cells['low']:,}; T/NK receivers {n_cells['tnk']:,}. "
        f"Table: `tables/cellchat_focused_lr.tsv`. This is a mean-product score, not a "
        f"permutation p-value."
    )
    (out_dir / "cellchat_run.json").write_text(
        json.dumps({**rec, "generated_at": datetime.now(timezone.utc).isoformat()}, indent=2) + "\n"
    )
    return rec
