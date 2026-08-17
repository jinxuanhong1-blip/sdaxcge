#!/usr/bin/env python3
"""Patient-level malignant CLDN4 vs T/NK on leftover TISCH NSCLC objects.

CLDN4-only. No dual-high. Per series, then leftover merges with honest n.
Combo + high-end (Q4 vs Q1). CellChat-style only if ρ<0 p<0.05 or a clear
Q4 vs Q1 drop.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import (
    DATASETS,
    EPITHELIAL_NONMALIGNANT,
    EXCLUDED_SET,
    GENES,
    MALIGNANT,
    MIN_EPI_CELLS,
    MIN_N_MERGE,
    MIN_N_Q4,
    MIN_N_SPEARMAN,
    MIN_TNK_CELLS,
    TNK_LINEAGES,
)
from lib_io import keep_tumor_unit, load_tisch_genes, pick_tissue_col, pick_unit_columns
from lib_stats import q4_vs_q1, spearman


def _norm_label(x) -> str:
    return str(x).strip()


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (math.isnan(p) or math.isinf(p))):
        return "—"
    if p < 0.001:
        return f"{p:.3g}"
    return f"{p:.3f}"


def fmt_rho(r) -> str:
    if r is None or (isinstance(r, float) and math.isnan(r)):
        return "—"
    return f"{r:+.3f}"


def score_units(ds: str, data_dir: Path) -> tuple[dict, pd.DataFrame]:
    spec = DATASETS[ds]
    meta_path = data_dir / f"{ds}_CellMetainfo_table.tsv"
    h5_path = data_dir / f"{ds}_expression.h5"
    rec: dict = {
        "dataset": ds,
        "gse": spec["gse"],
        "pmid": spec.get("pmid"),
        "paper": spec.get("paper"),
        "ici_gallery": spec.get("ici_gallery"),
        "in_excluded_set": spec["gse"] in EXCLUDED_SET,
        "skip_reason": None,
        "usable": False,
    }
    empty = pd.DataFrame()
    if spec["gse"] in EXCLUDED_SET:
        rec["skip_reason"] = "in 131907/148071/205335/207422 set"
        return rec, empty
    if not meta_path.exists():
        rec["skip_reason"] = "CellMetainfo missing"
        return rec, empty

    meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
    cell_col = meta.columns[0]
    meta = meta.set_index(cell_col)
    meta.index = meta.index.astype(str)
    lineage_col = next((c for c in meta.columns if "major-lineage" in c.lower()), None)
    if lineage_col is None:
        rec["skip_reason"] = "no major-lineage column"
        return rec, empty
    meta["_lineage"] = meta[lineage_col].map(_norm_label)
    rec["n_cells_meta"] = int(len(meta))
    rec["lineages"] = meta["_lineage"].value_counts().to_dict()
    rec["n_malignant"] = int(meta["_lineage"].isin(MALIGNANT).sum())
    rec["n_epithelial_nonmalig"] = int(meta["_lineage"].isin(EPITHELIAL_NONMALIGNANT).sum())
    rec["n_tnk"] = int(meta["_lineage"].isin(TNK_LINEAGES).sum())
    patient_col, sample_col = pick_unit_columns(meta)
    rec["patient_col"] = patient_col
    rec["sample_col"] = sample_col
    tissue_col = pick_tissue_col(meta)
    rec["tissue_col"] = tissue_col
    if tissue_col:
        rec["tissue_values"] = meta[tissue_col].astype(str).value_counts().to_dict()

    has_cldn4_comp = rec["n_malignant"] + rec["n_epithelial_nonmalig"] > 0
    has_tnk = rec["n_tnk"] > 0
    rec["has_cldn4_compartment"] = has_cldn4_comp
    rec["has_tnk"] = has_tnk
    if not has_cldn4_comp and not has_tnk:
        rec["skip_reason"] = "both CLDN4 compartment and T/NK absent"
        return rec, empty
    if not has_cldn4_comp:
        rec["skip_reason"] = "no epithelial/malignant cells; CLDN4 cannot be scored"
        return rec, empty
    if not has_tnk:
        rec["skip_reason"] = "no T/NK cells"
        return rec, empty
    if patient_col is None and sample_col is None:
        rec["skip_reason"] = "no Patient/Sample column"
        return rec, empty
    if not h5_path.exists():
        rec["skip_reason"] = "expression.h5 not downloaded (size/gate)"
        return rec, empty

    expr, present, h5_audit = load_tisch_genes(h5_path, GENES)
    rec["h5"] = h5_audit
    rec["genes_present"] = present
    if "CLDN4" not in present:
        rec["skip_reason"] = "CLDN4 absent from TISCH h5"
        return rec, empty

    common = meta.index.intersection(expr.index)
    rec["n_cells_aligned"] = int(len(common))
    if len(common) < 50:
        rec["skip_reason"] = f"barcode alignment too small (n={len(common)})"
        return rec, empty
    meta = meta.loc[common]
    expr = expr.loc[common]

    sample_for_filter = (
        meta[sample_col] if sample_col else (meta[patient_col] if patient_col else pd.Series([""] * len(meta), index=meta.index))
    )
    tissue_for_filter = meta[tissue_col] if tissue_col else pd.Series([None] * len(meta), index=meta.index)
    tumor_mask = pd.Series(
        [keep_tumor_unit(s, t) for s, t in zip(sample_for_filter.astype(str), tissue_for_filter)],
        index=meta.index,
    )
    rec["n_tumor_like_cells"] = int(tumor_mask.sum())
    if tumor_mask.sum() >= 50:
        meta = meta.loc[tumor_mask]
        expr = expr.loc[meta.index]
    else:
        rec["tissue_filter_note"] = "tumor filter left <50 cells; all cells kept"

    meta["_patient"] = meta[patient_col].astype(str) if patient_col else (
        meta[sample_col].astype(str) if sample_col else "unknown"
    )
    meta["_sample"] = meta[sample_col].astype(str) if sample_col else meta["_patient"]
    meta["_is_malig"] = meta["_lineage"].isin(MALIGNANT)
    meta["_is_epi"] = meta["_lineage"].isin(MALIGNANT | EPITHELIAL_NONMALIGNANT)
    meta["_is_tnk"] = meta["_lineage"].isin(TNK_LINEAGES)

    def _rows(unit_kind: str, unit_col: str) -> list[dict]:
        rows = []
        for unit, sub in meta.groupby(unit_col, sort=True):
            if str(unit) in {"nan", "None", ""}:
                continue
            n_cells = int(len(sub))
            n_mal = int(sub["_is_malig"].sum())
            n_epi = int(sub["_is_epi"].sum())
            n_tnk = int(sub["_is_tnk"].sum())
            if n_mal >= MIN_EPI_CELLS:
                score_mask = sub["_is_malig"].to_numpy()
                epi_def = "Malignant"
            else:
                score_mask = sub["_is_epi"].to_numpy()
                epi_def = "epithelial_like"
            n_score = int(score_mask.sum())
            eligible = n_score >= MIN_EPI_CELLS and n_tnk >= MIN_TNK_CELLS
            row = {
                "dataset": ds,
                "gse": spec["gse"],
                "unit_kind": unit_kind,
                "unit_id": str(unit),
                "patient": str(sub["_patient"].iloc[0]),
                "n_cells": n_cells,
                "n_malignant": n_mal,
                "n_epithelial_like": n_epi,
                "n_tnk": n_tnk,
                "frac_tnk": n_tnk / n_cells if n_cells else np.nan,
                "epi_definition": epi_def,
                "n_scored_epi": n_score,
                "eligible": eligible,
                "malignant_only": epi_def == "Malignant",
            }
            idx = sub.index
            for g in present:
                v = expr.loc[idx, g].to_numpy()
                scored = v[score_mask]
                row[f"{g}_mean"] = float(np.mean(scored)) if n_score else np.nan
                row[f"{g}_pctpos"] = float(np.mean(scored > 0) * 100) if n_score else np.nan
            rows.append(row)
        return rows

    rows = _rows("patient", "_patient") + _rows("sample", "_sample")
    units = pd.DataFrame(rows)
    rec["n_patients_total"] = int((units["unit_kind"] == "patient").sum()) if len(units) else 0
    rec["n_patients_eligible"] = int(
        ((units["unit_kind"] == "patient") & units["eligible"]).sum()
    ) if len(units) else 0
    rec["n_patients_malignant_eligible"] = int(
        ((units["unit_kind"] == "patient") & units["eligible"] & units["malignant_only"]).sum()
    ) if len(units) else 0
    rec["n_samples_eligible"] = int(
        ((units["unit_kind"] == "sample") & units["eligible"]).sum()
    ) if len(units) else 0
    rec["usable"] = rec["n_patients_eligible"] > 0 or rec["n_samples_eligible"] > 0
    if not rec["usable"] and rec["skip_reason"] is None:
        rec["skip_reason"] = (
            f"no eligible patient/sample (need ≥{MIN_EPI_CELLS} scored epi/malig "
            f"and ≥{MIN_TNK_CELLS} T/NK)"
        )
    return rec, units


def contrast_stats(sub: pd.DataFrame, label: str, cohorts: str, unit_kind: str, family: str) -> dict:
    n = int(len(sub))
    rho, p, n_sp = spearman(sub["CLDN4_mean"], sub["frac_tnk"])
    rho_pct, p_pct, _ = spearman(sub["CLDN4_pctpos"], sub["frac_tnk"])
    q = q4_vs_q1(sub["CLDN4_mean"], sub["frac_tnk"]) if n >= MIN_N_Q4 else None
    rec = {
        "label": label,
        "family": family,
        "cohorts": cohorts,
        "unit_kind": unit_kind,
        "k": int(sub["gse"].nunique()) if "gse" in sub.columns and len(sub) else 0,
        "n": n,
        "n_malignant": int(sub["malignant_only"].sum()) if len(sub) else 0,
        "n_epithelial_like": int((~sub["malignant_only"]).sum()) if len(sub) else 0,
        "rho_mean": rho,
        "p_mean": p,
        "rho_pct": rho_pct,
        "p_pct": p_pct,
        "n_spearman": n_sp,
        "spearman_ok": n >= MIN_N_SPEARMAN and math.isfinite(rho),
        "merge_ok": n >= MIN_N_MERGE,
        "r_rb": q["r_rb"] if q else np.nan,
        "p_q4q1": q["p"] if q else np.nan,
        "n_q1": q["n_q1"] if q else 0,
        "n_q4": q["n_q4"] if q else 0,
        "n_compared": q["n_compared"] if q else 0,
        "delta_median": q["delta_median"] if q else np.nan,
        "thin_q4q1": True if q is None else q["thin"],
        "q4q1_ok": q is not None,
        "clear_q4_drop": bool(
            q is not None
            and (not q["thin"])
            and q["r_rb"] < 0
            and (q["p"] < 0.05 or q["delta_median"] < -0.05)
        ),
        "rho_neg_sig": bool(n >= MIN_N_SPEARMAN and math.isfinite(rho) and rho < 0 and p < 0.05),
        "cellchat_trigger": False,
    }
    rec["cellchat_trigger"] = bool(rec["rho_neg_sig"] or rec["clear_q4_drop"])
    rec["note"] = ""
    if n < MIN_N_SPEARMAN:
        rec["note"] = f"n={n}<{MIN_N_SPEARMAN}; Spearman not a claim"
    elif n < MIN_N_MERGE:
        rec["note"] = f"n={n}<{MIN_N_MERGE}; merge threshold not met"
    if q is None and n >= MIN_N_Q4:
        rec["note"] = (rec["note"] + "; Q4 vs Q1 collapsed").strip("; ")
    elif q is not None and q["thin"]:
        rec["note"] = (rec["note"] + "; Q4 vs Q1 thin").strip("; ")
    return rec


def enumerate_combos(units: pd.DataFrame, unit_kind: str, malignant_only: bool) -> pd.DataFrame:
    work = units.loc[(units["unit_kind"] == unit_kind) & units["eligible"]].copy()
    if malignant_only:
        work = work.loc[work["malignant_only"]].copy()
    series = sorted(work["gse"].unique())
    rows = []
    family = f"{unit_kind}_{'malignant' if malignant_only else 'mixed'}"
    for k in range(1, len(series) + 1):
        for combo in itertools.combinations(series, k):
            sub = work.loc[work["gse"].isin(combo)]
            if len(sub) == 0:
                continue
            rec = contrast_stats(
                sub,
                label="+".join(combo),
                cohorts="+".join(combo),
                unit_kind=unit_kind,
                family=family,
            )
            rec["malignant_only_filter"] = malignant_only
            rows.append(rec)
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def highlighted(combos: pd.DataFrame) -> pd.DataFrame:
    if combos is None or len(combos) == 0:
        return pd.DataFrame()
    hit = combos.loc[combos["cellchat_trigger"] | (combos["rho_neg_sig"]) | (combos["clear_q4_drop"])].copy()
    return hit


def extra_figures(units: pd.DataFrame, audits: list[dict], combos: pd.DataFrame, fig_dir: Path) -> list[str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    written = []

    # Inventory
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    names, n_pat, n_mal, n_tnk = [], [], [], []
    for a in audits:
        names.append(a["gse"])
        n_pat.append(a.get("n_patients_eligible") or 0)
        n_mal.append(1 if a.get("n_malignant", 0) > 0 else 0)
        n_tnk.append(1 if a.get("has_tnk") else 0)
    x = np.arange(len(names))
    ax.bar(x, n_pat, color="#1f4e79")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=25, ha="right")
    ax.set_ylabel("eligible patients")
    ax.set_title("Leftover TISCH hunt — eligible patients (honest n)")
    for i, a in enumerate(audits):
        if a.get("skip_reason"):
            ax.text(i, 0.05, "skip", ha="center", va="bottom", fontsize=7, color="#8b1e3f", rotation=90)
    fig.tight_layout()
    p = fig_dir / "fig_inventory_eligible_patients.png"
    fig.savefig(p, dpi=160)
    fig.savefig(fig_dir / "fig_inventory_eligible_patients.pdf")
    plt.close(fig)
    written.append(str(p))

    # Per-series lineage
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.2))
    axes = axes.ravel()
    for ax, a in zip(axes, audits):
        lin = a.get("lineages") or {}
        items = sorted(lin.items(), key=lambda kv: -kv[1])[:8]
        if not items:
            ax.set_title(a["gse"])
            ax.text(0.5, 0.5, a.get("skip_reason") or "no lineages", ha="center", va="center", wrap=True)
            ax.axis("off")
            continue
        labs, vals = zip(*items)
        ax.barh(range(len(labs))[::-1], vals[::-1] if False else vals, color="#4a7c59")
        ax.set_yticks(range(len(labs)))
        ax.set_yticklabels(labs, fontsize=7)
        ax.invert_yaxis()
        ax.set_title(f"{a['gse']}\n{a.get('n_cells_meta', 0):,} cells", fontsize=8)
    fig.suptitle("TISCH major-lineage (leftover hunt)", fontsize=11)
    fig.tight_layout()
    p = fig_dir / "fig_lineage_bars.png"
    fig.savefig(p, dpi=160)
    fig.savefig(fig_dir / "fig_lineage_bars.pdf")
    plt.close(fig)
    written.append(str(p))

    # Scatter leftover patients
    pat = units.loc[(units["unit_kind"] == "patient") & units["eligible"]].copy()
    if len(pat):
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        colors = {
            "GSE117570": "#1f4e79",
            "GSE146100": "#c47b17",
            "GSE149655": "#4a7c59",
            "GSE143423": "#8b1e3f",
        }
        for gse, sub in pat.groupby("gse"):
            ax.scatter(
                sub["CLDN4_mean"],
                sub["frac_tnk"],
                s=70,
                c=colors.get(gse, "0.4"),
                label=f"{gse} n={len(sub)} ({sub['epi_definition'].iloc[0]})",
                edgecolors="k",
                linewidths=0.4,
            )
        if len(pat) >= 3:
            rho, p, n = spearman(pat["CLDN4_mean"], pat["frac_tnk"])
            ax.set_title(f"Leftover patients · CLDN4 mean vs T/NK\nn={n}  ρ={fmt_rho(rho)}  p={fmt_p(p)}")
        else:
            ax.set_title(f"Leftover patients · CLDN4 mean vs T/NK (n={len(pat)})")
        ax.set_xlabel("Malignant/epithelial CLDN4 mean (TISCH log2(TPM/10+1))")
        ax.set_ylabel("T/NK fraction")
        ax.legend(fontsize=7, loc="best")
        fig.tight_layout()
        p = fig_dir / "fig_scatter_leftover_patients.png"
        fig.savefig(p, dpi=160)
        fig.savefig(fig_dir / "fig_scatter_leftover_patients.pdf")
        plt.close(fig)
        written.append(str(p))

    samp = units.loc[(units["unit_kind"] == "sample") & units["eligible"]].copy()
    if len(samp):
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        colors = {
            "GSE117570": "#1f4e79",
            "GSE146100": "#c47b17",
            "GSE149655": "#4a7c59",
            "GSE143423": "#8b1e3f",
        }
        for gse, sub in samp.groupby("gse"):
            ax.scatter(
                sub["CLDN4_mean"],
                sub["frac_tnk"],
                s=70,
                c=colors.get(gse, "0.4"),
                label=f"{gse} n={len(sub)}",
                edgecolors="k",
                linewidths=0.4,
            )
        if len(samp) >= 3:
            rho, p, n = spearman(samp["CLDN4_mean"], samp["frac_tnk"])
            ax.set_title(f"Leftover samples/nodules · CLDN4 vs T/NK\nn={n}  ρ={fmt_rho(rho)}  p={fmt_p(p)}")
        ax.set_xlabel("CLDN4 mean (TISCH)")
        ax.set_ylabel("T/NK fraction")
        ax.legend(fontsize=7)
        fig.tight_layout()
        p = fig_dir / "fig_scatter_leftover_samples.png"
        fig.savefig(p, dpi=160)
        fig.savefig(fig_dir / "fig_scatter_leftover_samples.pdf")
        plt.close(fig)
        written.append(str(p))

    # Combo bars
    if combos is not None and len(combos):
        show = combos.loc[combos["family"] == "patient_mixed"].copy()
        if len(show):
            show = show.sort_values("n")
            fig, ax = plt.subplots(figsize=(8.8, max(3.2, 0.38 * len(show) + 1.2)))
            y = np.arange(len(show))
            cols = ["#8b1e3f" if (r < 0 if math.isfinite(r) else False) else "#1f4e79" for r in show["rho_mean"]]
            ax.barh(y, show["rho_mean"].fillna(0), color=cols)
            ax.axvline(0, color="0.3", lw=0.8)
            ax.set_yticks(y)
            ax.set_yticklabels([f"{c} (n={n})" for c, n in zip(show["cohorts"], show["n"])], fontsize=7)
            ax.set_xlabel("Spearman ρ (CLDN4 mean vs T/NK)")
            ax.set_title("Leftover patient-level combo hunt (mixed malignant + epi-like)")
            fig.tight_layout()
            p = fig_dir / "fig_combo_patient_mixed_bars.png"
            fig.savefig(p, dpi=160)
            fig.savefig(fig_dir / "fig_combo_patient_mixed_bars.pdf")
            plt.close(fig)
            written.append(str(p))

        show = combos.loc[combos["family"] == "sample_mixed"].copy()
        if len(show):
            show = show.sort_values("n")
            fig, ax = plt.subplots(figsize=(8.8, max(3.2, 0.38 * len(show) + 1.2)))
            y = np.arange(len(show))
            cols = ["#8b1e3f" if (r < 0 if math.isfinite(r) else False) else "#1f4e79" for r in show["rho_mean"]]
            ax.barh(y, show["rho_mean"].fillna(0), color=cols)
            ax.axvline(0, color="0.3", lw=0.8)
            ax.set_yticks(y)
            ax.set_yticklabels([f"{c} (n={n})" for c, n in zip(show["cohorts"], show["n"])], fontsize=7)
            ax.set_xlabel("Spearman ρ (CLDN4 mean vs T/NK)")
            ax.set_title("Leftover sample-level combo hunt (nodules not independent patients)")
            fig.tight_layout()
            p = fig_dir / "fig_combo_sample_mixed_bars.png"
            fig.savefig(p, dpi=160)
            fig.savefig(fig_dir / "fig_combo_sample_mixed_bars.pdf")
            plt.close(fig)
            written.append(str(p))

    # Q4 box on largest merge with Q4
    if combos is not None and len(combos):
        qok = combos.loc[combos["q4q1_ok"]].sort_values("n", ascending=False)
        if len(qok):
            top = qok.iloc[0]
            kind = top["unit_kind"]
            gses = top["cohorts"].split("+")
            sub = units.loc[(units["unit_kind"] == kind) & units["eligible"] & units["gse"].isin(gses)]
            if top.get("family", "").endswith("malignant"):
                sub = sub.loc[sub["malignant_only"]]
            q = q4_vs_q1(sub["CLDN4_mean"], sub["frac_tnk"])
            if q is not None:
                fig, ax = plt.subplots(figsize=(4.6, 4.2))
                labs = np.array(q["q_labels"])
                imm = np.array(q["immune"])
                ax.boxplot(
                    [imm[labs == "Q1"], imm[labs == "Q4"]],
                    tick_labels=[f"Q1 n={q['n_q1']}", f"Q4 n={q['n_q4']}"],
                    widths=0.55,
                )
                ax.set_ylabel("T/NK fraction")
                ax.set_title(
                    f"Q4 vs Q1 · {top['cohorts']}\n{kind} n={q['n']}  r={fmt_rho(q['r_rb'])}  p={fmt_p(q['p'])}"
                    + ("  thin" if q["thin"] else "")
                )
                fig.tight_layout()
                p = fig_dir / "fig_q4q1_largest_merge.png"
                fig.savefig(p, dpi=160)
                fig.savefig(fig_dir / "fig_q4q1_largest_merge.pdf")
                plt.close(fig)
                written.append(str(p))

    # Per-patient table figure
    if len(pat):
        fig, ax = plt.subplots(figsize=(9.2, max(2.6, 0.38 * len(pat) + 1.4)))
        ax.axis("off")
        cols = ["gse", "unit_id", "epi_definition", "n_scored_epi", "n_tnk", "CLDN4_mean", "CLDN4_pctpos", "frac_tnk"]
        tab = pat[cols].copy()
        tab["CLDN4_mean"] = tab["CLDN4_mean"].map(lambda v: f"{v:.3f}")
        tab["CLDN4_pctpos"] = tab["CLDN4_pctpos"].map(lambda v: f"{v:.1f}")
        tab["frac_tnk"] = tab["frac_tnk"].map(lambda v: f"{v:.3f}")
        table = ax.table(
            cellText=tab.values,
            colLabels=["series", "patient", "def", "n_epi", "n_TNK", "CLDN4 mean", "CLDN4 %pos", "frac T/NK"],
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1.0, 1.25)
        ax.set_title("Eligible leftover patients (honest rows)", pad=12)
        fig.tight_layout()
        p = fig_dir / "fig_patient_table.png"
        fig.savefig(p, dpi=160)
        fig.savefig(fig_dir / "fig_patient_table.pdf")
        plt.close(fig)
        written.append(str(p))

    return written


def write_finding(
    audits: list[dict],
    units: pd.DataFrame,
    combos: pd.DataFrame,
    hits: pd.DataFrame,
    figures: list[str],
    cellchat: dict,
    out: Path,
) -> None:
    pat = units.loc[(units["unit_kind"] == "patient") & units["eligible"]].copy() if len(units) else pd.DataFrame()
    mal = pat.loc[pat["malignant_only"]].copy() if len(pat) else pd.DataFrame()
    samp = units.loc[(units["unit_kind"] == "sample") & units["eligible"]].copy() if len(units) else pd.DataFrame()

    def _md_series_table() -> str:
        lines = [
            "| Series | cells | Malignant | Epi-like | T/NK | eligible patients | eligible samples | CLDN4 in h5 | skip |",
            "|---|---:|---:|---:|---:|---:|---:|---|---|",
        ]
        for a in audits:
            lines.append(
                "| {gse} | {n_cells} | {n_mal} | {n_epi} | {n_tnk} | {n_pe} | {n_se} | {cldn} | {skip} |".format(
                    gse=a["gse"],
                    n_cells=a.get("n_cells_meta", 0),
                    n_mal=a.get("n_malignant", 0),
                    n_epi=a.get("n_epithelial_nonmalig", 0),
                    n_tnk=a.get("n_tnk", 0),
                    n_pe=a.get("n_patients_eligible", 0) or 0,
                    n_se=a.get("n_samples_eligible", 0) or 0,
                    cldn="yes" if "CLDN4" in (a.get("genes_present") or []) else ("n/a" if a.get("skip_reason") else "no"),
                    skip=(a.get("skip_reason") or "usable").replace("|", "/"),
                )
            )
        return "\n".join(lines)

    def _md_patient_rows() -> str:
        if not len(pat):
            return "_No eligible leftover patients._"
        lines = [
            "| series | patient | def | n_scored | n_TNK | CLDN4 mean | CLDN4 %pos | frac T/NK |",
            "|---|---|---|---:|---:|---:|---:|---:|",
        ]
        for _, r in pat.sort_values(["gse", "unit_id"]).iterrows():
            lines.append(
                f"| {r['gse']} | {r['unit_id']} | {r['epi_definition']} | {int(r['n_scored_epi'])} | "
                f"{int(r['n_tnk'])} | {r['CLDN4_mean']:.3f} | {r['CLDN4_pctpos']:.1f} | {r['frac_tnk']:.3f} |"
            )
        return "\n".join(lines)

    def _combo_block(family: str, title: str) -> str:
        if combos is None or not len(combos):
            return f"### {title}\n\n_No combos._\n"
        sub = combos.loc[combos["family"] == family].sort_values(["k", "n"], ascending=[True, False])
        if not len(sub):
            return f"### {title}\n\n_No combos in this family._\n"
        lines = [
            f"### {title}",
            "",
            "| k | N | cohorts | ρ mean (p) | ρ %pos (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | trigger | note |",
            "|---:|---:|---|---|---|---|---|---|",
        ]
        for _, r in sub.iterrows():
            qtxt = "—"
            if r["q4q1_ok"]:
                qtxt = (
                    f"{fmt_rho(r['r_rb'])} ({fmt_p(r['p_q4q1'])}; "
                    f"{int(r['n_q1'])}/{int(r['n_q4'])})"
                    + (" thin" if r["thin_q4q1"] else "")
                )
            trig = "yes" if r["cellchat_trigger"] else "no"
            lines.append(
                f"| {int(r['k'])} | {int(r['n'])} | {r['cohorts']} | "
                f"{fmt_rho(r['rho_mean'])} ({fmt_p(r['p_mean'])}) | "
                f"{fmt_rho(r['rho_pct'])} ({fmt_p(r['p_pct'])}) | {qtxt} | {trig} | {r['note']} |"
            )
        return "\n".join(lines) + "\n"

    n_pat = int(len(pat))
    n_mal = int(len(mal))
    n_samp = int(len(samp))
    merge_pat = contrast_stats(pat, "all_leftover_patients", "+".join(sorted(pat["gse"].unique())) if n_pat else "", "patient", "patient_mixed") if n_pat else None
    merge_mal = contrast_stats(mal, "malignant_leftover_patients", "+".join(sorted(mal["gse"].unique())) if n_mal else "", "patient", "patient_malignant") if n_mal else None
    merge_samp = contrast_stats(samp, "all_leftover_samples", "+".join(sorted(samp["gse"].unique())) if n_samp else "", "sample", "sample_mixed") if n_samp else None

    hit_n = int(len(hits)) if hits is not None else 0
    if hit_n:
        hit_lines = ["Highlighted leftover combos (ρ<0 p<0.05 or clear Q4 vs Q1 drop):", ""]
        for _, r in hits.iterrows():
            hit_lines.append(
                f"- {r['family']} · {r['cohorts']} · k={int(r['k'])} · N={int(r['n'])} · "
                f"ρ={fmt_rho(r['rho_mean'])} p={fmt_p(r['p_mean'])} · "
                f"Q4vsQ1 r={fmt_rho(r['r_rb'])} p={fmt_p(r['p_q4q1'])}"
            )
        hit_md = "\n".join(hit_lines)
    else:
        hit_md = (
            "Highlighted leftover combos: **none**. Empty hunt. No leftover series or merge "
            "reached ρ<0 with p<0.05, and no clear (non-thin) Q4 vs Q1 drop."
        )

    cc_md = cellchat.get("finding_md") or (
        "CellChat-style **not run**. Trigger is ρ<0 with p<0.05 or a clear Q4 vs Q1 drop "
        "on a leftover series/merge. That trigger did not fire."
    )

    fig_list = "\n".join(f"- `{Path(p).name}`" for p in figures) if figures else "- (none)"

    def _merge_line(rec, name):
        if rec is None:
            return f"- **{name}:** n=0"
        return (
            f"- **{name}:** n={rec['n']} (malignant {rec['n_malignant']}, epi-like {rec['n_epithelial_like']}) · "
            f"ρ={fmt_rho(rec['rho_mean'])} p={fmt_p(rec['p_mean'])} · "
            f"Q4 vs Q1 r={fmt_rho(rec['r_rb'])} p={fmt_p(rec['p_q4q1'])} "
            f"(n_Q1={rec['n_q1']}, n_Q4={rec['n_q4']}"
            + ("; thin" if rec['thin_q4q1'] else "")
            + f") · n≥8 merge={'yes' if rec['merge_ok'] else 'no'} · {rec['note']}"
        )

    md = f"""# FINDING — leftover TISCH NSCLC merge, CLDN4-only combo + high-end

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. Public TISCH2
`expression.h5` + `CellMetainfo` only. Files >2 GB skipped. Series with
neither a CLDN4 compartment (epithelial/malignant) nor T/NK skipped.

Hunt list (not in the 131907/148071/205335/207422 set):
GSE117570, GSE139555, GSE146100, GSE149655, GSE143423, GSE127471.

Patient is the unit. Sample/nodule rows are a sensitivity (GSE146100 is
1 patient / 3 nodules). Eligible unit: ≥{MIN_EPI_CELLS} scored
malignant cells (else epithelial-like) and ≥{MIN_TNK_CELLS} T/NK.
Tumor-like tissue only. TISCH values are MAESTRO `log2(TPM/10+1)`.

p-values are descriptive. Combo table:
[`tables/highlighted_combos.tsv`](tables/highlighted_combos.tsv)
(empty header-only is an honest empty hunt).

## Verdict

{hit_md}

Patient-level leftover merge does **not** reach n≥8
(malignant-only n={n_mal}; mixed malignant+epi-like n={n_pat}).
Sample/nodule merge n={n_samp} is labeled as samples, not patients.

{_merge_line(merge_mal, "malignant leftover patients")}
{_merge_line(merge_pat, "mixed leftover patients")}
{_merge_line(merge_samp, "mixed leftover samples/nodules")}

## Series audit

{_md_series_table()}

GSE139555 is T/immune-sorted (T/NK present, no epithelium → CLDN4 not
scorable). GSE127471 is PBMC-only. Both fail the “CLDN4 present” gate;
h5 was not downloaded. GSE146100 and GSE149655 have **no TISCH
Malignant call** — they enter the mixed merge as epithelial-like only.

## Eligible leftover patients

{_md_patient_rows()}

## Per-series Spearman (patient-level, CLDN4 mean vs T/NK)

Spearman only if n≥{MIN_N_SPEARMAN}. Smaller n is listed, not tested.

"""
    # per-series rows
    md += "| series | def | n | ρ (p) | Q4 vs Q1 | note |\n|---|---|---:|---|---|---|\n"
    if len(pat):
        for gse, sub in pat.groupby("gse"):
            rec = contrast_stats(sub, gse, gse, "patient", "single")
            qtxt = "—"
            if rec["q4q1_ok"]:
                qtxt = f"{fmt_rho(rec['r_rb'])} ({fmt_p(rec['p_q4q1'])})"
            md += (
                f"| {gse} | {sub['epi_definition'].iloc[0]} | {rec['n']} | "
                f"{fmt_rho(rec['rho_mean'])} ({fmt_p(rec['p_mean'])}) | {qtxt} | {rec['note']} |\n"
            )
    else:
        md += "| — | — | 0 | — | — | no eligible leftover patients |\n"

    md += "\n" + _combo_block("patient_malignant", "Malignant-only leftover combos (patient)")
    md += "\n" + _combo_block("patient_mixed", "Mixed leftover combos (patient; malignant + epi-like)")
    md += "\n" + _combo_block("sample_malignant", "Malignant-only leftover combos (sample/nodule)")
    md += "\n" + _combo_block("sample_mixed", "Mixed leftover combos (sample/nodule; not independent patients)")

    md += f"""
## CellChat-style

{cc_md}

## What was not done

- No dual-high TACSTD2×CLDN4 score.
- GSE131907 / GSE148071 / GSE205335 / GSE207422 were not re-scored.
- Other TISCH NSCLC objects already in `methods/tisch_nsclc_pool/`
  (GSE127465, GSE153935, GSE162498, GSE150660, EMTAB6149) were not
  added to this hunt list.
- Files >2 GB skipped. GSE139555 / GSE127471 h5 skipped (no CLDN4 compartment).
- Cell-level p-values are not primary evidence.

## Extra figures

{fig_list}

Reproduce:

```bash
python3 -m pip install -r methods/tisch_leftover_merge_cldn4/requirements.txt
python3 methods/tisch_leftover_merge_cldn4/download.py
python3 methods/tisch_leftover_merge_cldn4/analyze.py
```

Generated: {datetime.now(timezone.utc).isoformat()}
"""
    out.write_text(md)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/tisch_leftover_merge_cldn4")
    ap.add_argument("--out-dir", default="methods/tisch_leftover_merge_cldn4")
    args = ap.parse_args()
    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    tables = out_dir / "tables"
    figs = out_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    audits = []
    unit_frames = []
    for ds in DATASETS:
        print(f"== {ds}", flush=True)
        rec, units = score_units(ds, data_dir)
        audits.append(rec)
        if len(units):
            unit_frames.append(units)
        print(
            f"   skip={rec.get('skip_reason')} eligible_patients={rec.get('n_patients_eligible')} "
            f"usable={rec.get('usable')}",
            flush=True,
        )

    units = pd.concat(unit_frames, ignore_index=True) if unit_frames else pd.DataFrame()
    combo_parts = [
        enumerate_combos(units, "patient", True),
        enumerate_combos(units, "patient", False),
        enumerate_combos(units, "sample", True),
        enumerate_combos(units, "sample", False),
    ]
    combos = pd.concat([c for c in combo_parts if len(c)], ignore_index=True) if any(len(c) for c in combo_parts) else pd.DataFrame()
    hits = highlighted(combos)

    units.to_csv(tables / "per_unit_metrics.tsv", sep="\t", index=False)
    if len(combos):
        combos.to_csv(tables / "combo_table.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(
            columns=[
                "label", "family", "cohorts", "unit_kind", "k", "n", "rho_mean", "p_mean",
                "r_rb", "p_q4q1", "cellchat_trigger", "note",
            ]
        ).to_csv(tables / "combo_table.tsv", sep="\t", index=False)
    if len(hits):
        hits.to_csv(tables / "highlighted_combos.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(
            columns=[
                "label", "family", "cohorts", "unit_kind", "k", "n", "rho_mean", "p_mean",
                "r_rb", "p_q4q1", "n_q1", "n_q4", "cellchat_trigger", "note",
            ]
        ).to_csv(tables / "highlighted_combos.tsv", sep="\t", index=False)

    figures = extra_figures(units, audits, combos, figs)

    trigger_rows = hits.to_dict(orient="records") if len(hits) else []
    cellchat = {"ran": False, "reason": "no leftover series/merge with ρ<0 p<0.05 or clear Q4 vs Q1 drop"}
    if trigger_rows:
        from cellchat_style import run_cellchat_if_triggered

        cellchat = run_cellchat_if_triggered(
            trigger_rows=trigger_rows,
            units=units,
            data_dir=data_dir,
            out_dir=out_dir,
        )
    else:
        (out_dir / "cellchat_skip.json").write_text(
            json.dumps(
                {
                    "ran": False,
                    "reason": cellchat["reason"],
                    "generated_at": datetime.now(timezone.utc).isoformat(),
                },
                indent=2,
            )
            + "\n"
        )

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "excluded_set": sorted(EXCLUDED_SET),
        "n_series_hunted": len(DATASETS),
        "n_usable_series": sum(1 for a in audits if a.get("usable")),
        "n_eligible_patients": int(((units["unit_kind"] == "patient") & units["eligible"]).sum()) if len(units) else 0,
        "n_eligible_malignant_patients": int(
            ((units["unit_kind"] == "patient") & units["eligible"] & units["malignant_only"]).sum()
        ) if len(units) else 0,
        "n_highlighted_combos": int(len(hits)) if hits is not None else 0,
        "cellchat": {k: cellchat.get(k) for k in ("ran", "reason", "n_pairs") if k in cellchat or k in ("ran", "reason")},
        "datasets": audits,
    }
    (out_dir / "summary.json").write_text(json.dumps(payload, indent=2) + "\n")

    write_finding(audits, units, combos, hits, figures, cellchat, out_dir / "FINDING.md")
    print(f"wrote {tables} {figs} {out_dir / 'FINDING.md'}")


if __name__ == "__main__":
    main()
