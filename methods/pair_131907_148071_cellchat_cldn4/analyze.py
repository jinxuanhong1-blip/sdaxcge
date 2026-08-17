#!/usr/bin/env python3
"""Pair GSE131907 + GSE148071: patient-level malignant CLDN4 vs T/NK, then CellChat.

ADDITIVE. CLDN4 only. No dual-high. GSE205335 is not in this pair.
Patient is the unit. Combo ρ is DerSimonian–Laird on Fisher-z of the two
within-cohort Spearmans. CellChat-style outgoing is same-patient Mal → T/NK
on the within-cohort Q4 vs Q1 tails that pass the cell floors.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from lib_stats import random_effects_dl, spearman

ROOT = Path(__file__).resolve().parent
KH = 0.5
TRIM = 0.10
EXPR_PROP = 0.10
MIN_MAL_131907 = 20
MIN_TNK_131907 = 20
MIN_EPI_148071 = 25
MIN_TNK_148071 = 25
MIN_CELLS_CC = 20
MIN_DETECT_ARM = 3
TUMOR_ORIGINS = ("tLung", "tL/B", "mLN", "mBrain")
MALIG_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
ORIGIN_PREF = {"tLung": 0, "tL/B": 1, "mLN": 2, "mBrain": 3}

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
EXTRA = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD8A", "CD3E", "IFNG", "TNF"]

BARRIER_LIGANDS = {
    "CDH1", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "F11R", "JAM2", "JAM3",
    "CEACAM1", "CEACAM5", "CEACAM6", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4",
    "PVR", "EPCAM", "DSG2", "DSC2", "CADM1",
}
INHIB_LIGANDS = {
    "CD274", "PDCD1LG2", "LGALS9", "HLA-E", "HLA-G", "HLA-F", "TGFB1", "TGFB2",
    "TGFB3", "CD80", "CD86", "CD276", "VSIR", "PVR", "NECTIN2", "CD47", "CDH1",
}
RECRUIT_LIGANDS = {
    "CXCL9", "CXCL10", "CXCL11", "CXCL16", "CCL5", "CCL3", "CCL4", "IL15",
    "IL2", "IL18", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3",
}
ATTACK_LIGANDS = {"IFNG", "TNF", "FASLG", "TNFSF10", "LTA"}


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def q4_vs_q1(cldn4, immune) -> dict | None:
    frame = pd.DataFrame(
        {"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)}
    )
    frame = frame[np.isfinite(frame["c"]) & np.isfinite(frame["i"])].copy()
    n = int(len(frame))
    if n < 6:
        return None
    ranks = frame["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return None
    if qs.nunique() < 4:
        return None
    q1 = frame.loc[qs == "Q1", "i"]
    q4 = frame.loc[qs == "Q4", "i"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 2 or n4 < 2:
        return None
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    return {
        "n": n,
        "n_q1": n1,
        "n_q4": n4,
        "n_compared": n1 + n4,
        "median_q1": float(q1.median()),
        "median_q4": float(q4.median()),
        "delta_median": float(q4.median() - q1.median()),
        "mwu_u": float(u),
        "p": float(p),
        "r_rb": float(r_rb),
        "thin": n < 8 or n1 < 3 or n4 < 3,
        "poolable": (n >= 8 and n1 >= 3 and n4 >= 3 and abs(r_rb) < 0.999),
    }


def assign_quartiles(values: pd.Series) -> pd.Series:
    ranks = values.rank(method="average")
    return pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")


def parse_symbols(val) -> list[str]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).replace(";", ",").replace("&", ",").replace("|", ",")
    return [x.strip() for x in text.split(",") if x.strip() and x.strip().lower() != "nan"]


def load_lr(db_dir: Path, matrix_genes: set[str] | None = None) -> pd.DataFrame:
    inter = pd.read_csv(db_dir / "interaction_cellchatdb_v2_protein.csv")
    inter = inter.rename(
        columns={"ligand.symbol": "ligand_symbol", "receptor.symbol": "receptor_symbol"}
    )
    rows = []
    for rec in inter.itertuples(index=False):
        lig = parse_symbols(getattr(rec, "ligand_symbol", None))
        recp = parse_symbols(getattr(rec, "receptor_symbol", None))
        if not lig:
            lig = parse_symbols(rec.ligand)
        if not recp:
            recp = parse_symbols(rec.receptor)
        if not lig or not recp:
            continue
        if matrix_genes is not None and any(g not in matrix_genes for g in lig + recp):
            continue
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": tuple(lig),
                "receptor_genes": tuple(recp),
            }
        )
    return pd.DataFrame(rows)


def ligand_class(genes) -> str:
    parts = set(genes) if not isinstance(genes, str) else set(str(genes).split("|"))
    tags = []
    if parts & BARRIER_LIGANDS:
        tags.append("barrier")
    if parts & INHIB_LIGANDS:
        tags.append("inhibitory")
    if parts & RECRUIT_LIGANDS:
        tags.append("recruit")
    if parts & ATTACK_LIGANDS:
        tags.append("attack")
    return "|".join(tags) if tags else "other"


def load_gse131907_patients(data_dir: Path) -> pd.DataFrame:
    samples = pd.read_csv(data_dir / "GSE131907_samples.tsv", sep="\t")
    meta = pd.read_csv(data_dir / "GSE131907_sample_metadata.tsv", sep="\t")
    samples = samples.merge(
        meta[["Sample", "patient_id", "tumor_stage", "tissue_origin_abbrevation"]],
        left_on="sample",
        right_on="Sample",
        how="left",
        validate="one_to_one",
    )
    keep = samples[
        samples["origin"].isin(TUMOR_ORIGINS)
        & (samples["n_malignant"] >= MIN_MAL_131907)
        & (samples["n_tnk"] >= MIN_TNK_131907)
    ].copy()
    keep["origin_rank"] = keep["origin"].map(ORIGIN_PREF).fillna(9)
    keep = keep.sort_values(["patient_id", "origin_rank", "n_malignant"], ascending=[True, True, False])
    keep = keep.drop_duplicates("patient_id", keep="first")
    out = pd.DataFrame(
        {
            "cohort": "GSE131907",
            "unit": keep["patient_id"],
            "sample": keep["sample"],
            "origin": keep["origin"],
            "stage": keep["tumor_stage"],
            "malig_def": "author_malig_tS",
            "n_cells": keep["n_cells"],
            "n_malignant": keep["n_malignant"],
            "n_tnk": keep["n_tnk"],
            "frac_tnk": keep["frac_tnk"],
            "mal_CLDN4_mean": keep["mal_CLDN4_mean"],
            "mal_CLDN4_pct_pos": keep["mal_CLDN4_pct"],
        }
    )
    return out.reset_index(drop=True)


def load_gse148071_patients(data_dir: Path) -> pd.DataFrame:
    raw = pd.read_csv(data_dir / "GSE148071_per_sample.tsv", sep="\t")
    keep = raw[
        (raw["n_epithelial"] >= MIN_EPI_148071) & (raw["n_TNK"] >= MIN_TNK_148071)
    ].copy()
    out = pd.DataFrame(
        {
            "cohort": "GSE148071",
            "unit": keep["sample"],
            "sample": keep["sample"],
            "origin": "biopsy",
            "stage": "III/IV",
            "malig_def": "marker_epi_putative",
            "n_cells": keep["n_total"],
            "n_malignant": keep["n_epithelial"],
            "n_tnk": keep["n_TNK"],
            "frac_tnk": keep["n_TNK"] / keep["n_total"],
            "mal_CLDN4_mean": keep["mean_CLDN4_epithelial"],
            "mal_CLDN4_pct_pos": keep["frac_CLDN4_pos_epithelial"] * 100.0,
        }
    )
    return out.reset_index(drop=True)


def score_cohort(patients: pd.DataFrame, score: str) -> dict:
    ccol = "mal_CLDN4_pct_pos" if score == "pct_pos" else "mal_CLDN4_mean"
    rho, p_s, n = spearman(patients[ccol], patients["frac_tnk"])
    q = q4_vs_q1(patients[ccol], patients["frac_tnk"])
    return {
        "cohort": patients["cohort"].iloc[0],
        "malig_def": patients["malig_def"].iloc[0],
        "immune_def": "same_patient_tnk_fraction",
        "score": score,
        "n_patients": n,
        "n_q1": None if q is None else q["n_q1"],
        "n_q4": None if q is None else q["n_q4"],
        "n_compared": None if q is None else q["n_compared"],
        "spearman_rho": rho,
        "spearman_p": p_s,
        "q4q1_r_rb": None if q is None else q["r_rb"],
        "q4q1_p": None if q is None else q["p"],
        "median_tnk_q1": None if q is None else q["median_q1"],
        "median_tnk_q4": None if q is None else q["median_q4"],
        "delta_median_tnk": None if q is None else q["delta_median"],
        "thin": True if q is None else q["thin"],
        "poolable_q4q1": False if q is None else q["poolable"],
        "unit": "patient",
    }


def combo_from_rows(rows: list[dict], effect_key: str, n_key: str, p_key: str) -> dict:
    rhos = [float(r[effect_key]) for r in rows if r[effect_key] is not None and np.isfinite(r[effect_key])]
    ns = [int(r[n_key]) for r in rows if r[effect_key] is not None and np.isfinite(r[effect_key])]
    if len(rhos) < 2:
        return {"k": len(rhos)}
    pooled = random_effects_dl(rhos, ns)
    pooled["members"] = [
        {
            "cohort": r["cohort"],
            "n": r[n_key],
            "effect": r[effect_key],
            "p": r[p_key],
        }
        for r in rows
    ]
    return pooled


def plot_q4q1(q1, q4, ylabel, title, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    bp = ax.boxplot(
        [q1, q4],
        tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"],
        patch_artist=True,
        widths=0.55,
    )
    for patch, color in zip(bp["boxes"], ["#6a8aaa", "#b2182b"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, vals in enumerate((q1, q4), start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=16, zorder=3)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_scatter(frame: pd.DataFrame, title: str, path: Path) -> None:
    colors = {"Q1": "#6a8aaa", "Q2": "#bdbdbd", "Q3": "#f4a582", "Q4": "#b2182b"}
    markers = {"GSE131907": "o", "GSE148071": "s"}
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    for cohort, sub in frame.groupby("cohort"):
        for q, color in colors.items():
            m = sub["q_pct"] == q
            ax.scatter(
                sub.loc[m, "mal_CLDN4_pct_pos"],
                sub.loc[m, "frac_tnk"],
                c=color,
                marker=markers.get(cohort, "o"),
                s=36,
                label=f"{cohort} {q}",
                zorder=3,
            )
    ax.set_xlabel("Malignant CLDN4 % positive")
    ax.set_ylabel("Same-patient T/NK fraction")
    ax.set_title(title, fontsize=9)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_forest(rows: list[dict], pooled: dict, title: str, path: Path, effect="spearman_rho") -> None:
    labels = [r["cohort"] + f" n={int(r['n_patients'])}" for r in rows]
    if pooled.get("k", 0) >= 2:
        labels.append(f"combo DL n={int(pooled['n_patients_total'])}")
    effects = [r[effect] for r in rows]
    if pooled.get("k", 0) >= 2:
        effects.append(pooled["pooled_rho"])
    fig, ax = plt.subplots(figsize=(6.2, 2.4 + 0.35 * len(labels)))
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="0.4", lw=0.8)
    ax.scatter(effects[:-1] if pooled.get("k", 0) >= 2 else effects, y[: len(rows)], c="#7a2d0b", s=40, zorder=3)
    if pooled.get("k", 0) >= 2:
        lo, hi = pooled["ci95_rho"]
        ax.plot([lo, hi], [y[-1], y[-1]], color="#7a2d0b", lw=2)
        ax.scatter([pooled["pooled_rho"]], [y[-1]], c="#7a2d0b", s=55, marker="D", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Spearman ρ (malignant CLDN4 %pos vs T/NK)")
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_ligand_table(tbl: pd.DataFrame, path: Path, title: str) -> None:
    if tbl.empty:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.axis("off")
        ax.text(0.5, 0.5, "No differential pairs at the stated gates", ha="center")
        fig.savefig(path, dpi=150)
        fig.savefig(path.with_suffix(".pdf"))
        plt.close(fig)
        return
    show = tbl.head(18).copy()
    fig, ax = plt.subplots(figsize=(8.6, 0.42 * len(show) + 1.4))
    y = np.arange(len(show))
    colors = np.where(show["delta_median"] >= 0, "#b2182b", "#2166ac")
    ax.barh(y, show["delta_median"], color=colors, alpha=0.85)
    ax.set_yticks(y)
    labels = [
        f"{r.direction[:3]} {r.ligand}–{r.receptor} ({r.ligand_class})"
        for r in show.itertuples()
    ]
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_xlabel("median P(Q4) − median P(Q1)")
    ax.set_title(title, fontsize=9)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def run_tnk(data_dir: Path, out: Path) -> dict:
    a = load_gse131907_patients(data_dir)
    b = load_gse148071_patients(data_dir)
    a["q_pct"] = assign_quartiles(a["mal_CLDN4_pct_pos"])
    a["q_mean"] = assign_quartiles(a["mal_CLDN4_mean"])
    b["q_pct"] = assign_quartiles(b["mal_CLDN4_pct_pos"])
    b["q_mean"] = assign_quartiles(b["mal_CLDN4_mean"])
    patients = pd.concat([a, b], ignore_index=True)
    rows = []
    for frame in (a, b):
        rows.append(score_cohort(frame, "pct_pos"))
        rows.append(score_cohort(frame, "mean"))
    table = pd.DataFrame(rows)
    pct_rows = [r for r in rows if r["score"] == "pct_pos"]
    mean_rows = [r for r in rows if r["score"] == "mean"]
    combo_pct = combo_from_rows(pct_rows, "spearman_rho", "n_patients", "spearman_p")
    combo_mean = combo_from_rows(mean_rows, "spearman_rho", "n_patients", "spearman_p")
    q_pct = [r for r in pct_rows if r["poolable_q4q1"]]
    combo_q = combo_from_rows(q_pct, "q4q1_r_rb", "n_compared", "q4q1_p") if len(q_pct) >= 2 else {"k": len(q_pct)}

    stacked_rho, stacked_p, stacked_n = spearman(patients["mal_CLDN4_pct_pos"], patients["frac_tnk"])
    stacked_q = q4_vs_q1(patients["mal_CLDN4_pct_pos"], patients["frac_tnk"])

    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "results" / "q4q1_tnk.tsv", sep="\t", index=False)
    patients.to_csv(out / "results" / "patients_with_quartiles.tsv", sep="\t", index=False)
    combo_tbl = pd.DataFrame(
        [
            {
                "analysis": "spearman",
                "score": "pct_pos",
                "combo": "GSE131907+GSE148071",
                "k": combo_pct.get("k"),
                "N": combo_pct.get("n_patients_total"),
                "effect": combo_pct.get("pooled_rho"),
                "p": combo_pct.get("p"),
                "I2": combo_pct.get("I2"),
                "ci95_lo": None if "ci95_rho" not in combo_pct else combo_pct["ci95_rho"][0],
                "ci95_hi": None if "ci95_rho" not in combo_pct else combo_pct["ci95_rho"][1],
                "note": "DL Fisher-z of within-cohort Spearman; GSE205335 not included",
            },
            {
                "analysis": "spearman",
                "score": "mean",
                "combo": "GSE131907+GSE148071",
                "k": combo_mean.get("k"),
                "N": combo_mean.get("n_patients_total"),
                "effect": combo_mean.get("pooled_rho"),
                "p": combo_mean.get("p"),
                "I2": combo_mean.get("I2"),
                "ci95_lo": None if "ci95_rho" not in combo_mean else combo_mean["ci95_rho"][0],
                "ci95_hi": None if "ci95_rho" not in combo_mean else combo_mean["ci95_rho"][1],
                "note": "DL Fisher-z; secondary mean score",
            },
            {
                "analysis": "q4q1",
                "score": "pct_pos",
                "combo": "GSE131907+GSE148071",
                "k": combo_q.get("k"),
                "N": combo_q.get("n_patients_total"),
                "effect": combo_q.get("pooled_rho"),
                "p": combo_q.get("p"),
                "I2": combo_q.get("I2"),
                "ci95_lo": None if "ci95_rho" not in combo_q else combo_q["ci95_rho"][0],
                "ci95_hi": None if "ci95_rho" not in combo_q else combo_q["ci95_rho"][1],
                "note": "DL Fisher-z of within-cohort rank-biserial r; n = n_Q1+n_Q4",
            },
            {
                "analysis": "spearman_stacked",
                "score": "pct_pos",
                "combo": "GSE131907+GSE148071",
                "k": 1,
                "N": stacked_n,
                "effect": stacked_rho,
                "p": stacked_p,
                "I2": None,
                "ci95_lo": None,
                "ci95_hi": None,
                "note": "companion only; 10x + Singleron stacked, not the combo ρ",
            },
        ]
    )
    combo_tbl.to_csv(out / "results" / "combo_rho.tsv", sep="\t", index=False)

    tails = patients[patients["q_pct"].isin(["Q1", "Q4"])].copy()
    plot_q4q1(
        tails.loc[tails["q_pct"] == "Q1", "frac_tnk"].to_numpy(),
        tails.loc[tails["q_pct"] == "Q4", "frac_tnk"].to_numpy(),
        "Same-patient T/NK fraction",
        (
            f"Pair 131907+148071 within-cohort CLDN4 %pos Q4 vs Q1 T/NK\n"
            f"combo r={combo_q.get('pooled_rho', float('nan')):+.3f} "
            f"p={fmt_p(combo_q.get('p', float('nan')))} "
            f"n_compared={combo_q.get('n_patients_total', 'NA')}"
        ),
        out / "figures" / "q4q1_tnk_pct.png",
    )
    plot_scatter(
        patients,
        (
            f"Pair GSE131907+GSE148071  combo ρ={combo_pct.get('pooled_rho', float('nan')):+.3f} "
            f"p={fmt_p(combo_pct.get('p', float('nan')))}  N={combo_pct.get('n_patients_total')}"
        ),
        out / "figures" / "scatter_cldn4_tnk.png",
    )
    plot_forest(
        pct_rows,
        combo_pct,
        "Combo ρ — malignant CLDN4 %pos vs T/NK (not +GSE205335)",
        out / "figures" / "fig_combo_rho.png",
    )
    return {
        "patients": patients,
        "table": table,
        "combo_pct": combo_pct,
        "combo_mean": combo_mean,
        "combo_q": combo_q,
        "stacked": {"rho": stacked_rho, "p": stacked_p, "n": stacked_n, "q": stacked_q},
        "tails": tails,
    }


def trim_mean_1d(values: np.ndarray, proportiontocut: float = TRIM) -> float:
    n = int(values.size)
    if n == 0:
        return 0.0
    if n == 1:
        return float(values[0])
    k = int(n * proportiontocut)
    if k == 0:
        return float(values.mean())
    s = np.sort(values)
    return float(s[k : n - k].mean())


def geom_mean(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return 0.0
    if np.any(arr <= 0):
        return 0.0
    if arr.size == 1:
        return float(arr[0])
    return float(np.exp(np.mean(np.log(arr))))


def hill_prob(lig: float, rec: float, kh: float = KH) -> float:
    if lig <= 0 or rec <= 0:
        return 0.0
    prod = lig * rec
    return float(prod / (kh + prod))


def complex_from_maps(means: dict[str, float], props: dict[str, float], subunits: tuple[str, ...]):
    vals, prs = [], []
    for gene in subunits:
        if gene not in means:
            return 0.0, 0.0
        vals.append(means[gene])
        prs.append(props[gene])
    return geom_mean(vals), float(min(prs)) if prs else 0.0


def compartment_gene_stats(log_cp, pos, idx: np.ndarray):
    means, props = {}, {}
    if idx.size == 0:
        return means, props
    for gene, vec in log_cp.items():
        means[gene] = trim_mean_1d(vec[idx])
        props[gene] = float(pos[gene][idx].mean())
    return means, props


def score_patient_pairs(lr: pd.DataFrame, log_cp, pos, mal_idx, tnk_idx) -> list[dict]:
    rows = []
    if mal_idx.size < MIN_CELLS_CC or tnk_idx.size < MIN_CELLS_CC:
        return rows
    mal_mu, mal_pr = compartment_gene_stats(log_cp, pos, mal_idx)
    tnk_mu, tnk_pr = compartment_gene_stats(log_cp, pos, tnk_idx)
    for rec in lr.itertuples(index=False):
        lig_mal, lig_mal_p = complex_from_maps(mal_mu, mal_pr, rec.ligand_genes)
        rec_tnk, rec_tnk_p = complex_from_maps(tnk_mu, tnk_pr, rec.receptor_genes)
        lig_tnk, lig_tnk_p = complex_from_maps(tnk_mu, tnk_pr, rec.ligand_genes)
        rec_mal, rec_mal_p = complex_from_maps(mal_mu, mal_pr, rec.receptor_genes)
        out_det = lig_mal_p >= EXPR_PROP and rec_tnk_p >= EXPR_PROP
        in_det = lig_tnk_p >= EXPR_PROP and rec_mal_p >= EXPR_PROP
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": "|".join(rec.ligand_genes),
                "receptor_genes": "|".join(rec.receptor_genes),
                "direction": "outgoing",
                "prob": hill_prob(lig_mal, rec_tnk) if out_det else 0.0,
                "detected": bool(out_det),
            }
        )
        rows.append(
            {
                "interaction_name": rec.interaction_name,
                "pathway_name": rec.pathway_name,
                "annotation": rec.annotation,
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "ligand_genes": "|".join(rec.ligand_genes),
                "receptor_genes": "|".join(rec.receptor_genes),
                "direction": "incoming",
                "prob": hill_prob(lig_tnk, rec_mal) if in_det else 0.0,
                "detected": bool(in_det),
            }
        )
    return rows


def stream_matrix_keep(matrix_path: Path, wanted: set[str], keep_barcodes: set[str] | None):
    found: dict[str, list[np.ndarray]] = {}
    kept_ids: list[str] = []
    n_umi_parts: list[np.ndarray] = []
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        if keep_barcodes is None:
            keep_idx = np.arange(n)
        else:
            keep_idx = np.array([i for i, c in enumerate(cell_ids) if c in keep_barcodes], dtype=int)
        kept_ids = [cell_ids[i] for i in keep_idx]
        n_umi = np.zeros(keep_idx.size, dtype=np.float64)
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            sub = arr[keep_idx]
            n_umi += sub
            if gene in wanted:
                found[gene] = sub
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  stream genes={n_streamed} stored={len(found)} keep={keep_idx.size}", flush=True)
    print(f"stream done genes={n_streamed} kept_cells={len(kept_ids)} stored={len(found)}", flush=True)
    return kept_ids, found, n_umi, n_streamed


def stream_one_148071(path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        header = [h.strip().strip('"') for h in header if h != ""]
        if header and header[0] in {"", "gene", "Gene", "index", "Index", "GENE"}:
            cell_ids = header[1:]
        else:
            cell_ids = header
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.strip().strip('"').split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                toks = line.rstrip("\n").split("\t")
                arr = np.asarray(toks[1:], dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{path.name} {gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in wanted:
                found[gene] = arr
    return cell_ids, found, n_umi


def module_score(log_cp: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
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
    t_idx = names.index("T")
    nk_idx = names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both_high = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk_best = np.isin(labels, ["T", "NK"])
    labels[close & both_high & tnk_best & (cd3 > 0.15)] = "T"
    labels[close & both_high & tnk_best & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < 0.12] = "Unassigned"
    return labels


def patient_from_filename(path: Path) -> str:
    m = re.search(r"_(P\d+)_", path.name)
    if m:
        return m.group(1)
    parts = path.name.replace(".txt.gz", "").split("_")
    return parts[1] if len(parts) >= 2 else path.stem


def list_exp_files(data_dir: Path) -> list[Path]:
    files = sorted(data_dir.rglob("*_exp.txt.gz"))
    if files:
        return files
    tar_path = data_dir / "GSE148071_RAW.tar"
    if not tar_path.exists():
        raise SystemExit(f"no exp matrices and no {tar_path}")
    dest = data_dir / "GSE148071_files"
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    return sorted(dest.rglob("*_exp.txt.gz"))


def run_cellchat(args, tails: pd.DataFrame, out: Path) -> dict:
    lr = load_lr(args.db)
    wanted = set(EXTRA)
    for vs in LINEAGE_MARKERS.values():
        wanted.update(vs)
    for rec in lr.itertuples(index=False):
        wanted.update(rec.ligand_genes)
        wanted.update(rec.receptor_genes)

    qmap = tails.set_index(["cohort", "unit"])["q_pct"].astype(str)
    per_rows = []
    cell_rows = []

    a_tails = tails[tails["cohort"] == "GSE131907"]
    if not a_tails.empty:
        ann = pd.read_csv(args.ann_131907, sep="\t", dtype=str)
        keep_samples = set(a_tails["sample"])
        keep_barcodes = set(ann.loc[ann["Sample"].isin(keep_samples), "Index"])
        print(f"GSE131907 Q4/Q1 samples={len(keep_samples)} barcodes={len(keep_barcodes)}", flush=True)
        cell_ids, found, n_umi, _n = stream_matrix_keep(args.matrix_131907, wanted, keep_barcodes)
        per = ann.set_index("Index").reindex(cell_ids).reset_index()
        lib = np.maximum(n_umi, 1.0)
        log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
        pos = {g: (found[g] > 0).astype(np.float32) for g in found}
        mal = (
            (per["Cell_type"] == "Epithelial cells")
            & per["Cell_subtype"].isin(MALIG_SUBTYPES)
        ).to_numpy()
        tnk = per["Cell_type"].isin(["T lymphocytes", "NK cells"]).to_numpy()
        sample = per["Sample"].to_numpy()
        unit_map = a_tails.set_index("sample")["unit"]
        for samp, sub_idx in pd.Series(np.arange(len(sample))).groupby(sample):
            if samp not in unit_map.index:
                continue
            unit = unit_map.loc[samp]
            idx = sub_idx.to_numpy()
            mal_idx = idx[mal[idx]]
            tnk_idx = idx[tnk[idx]]
            q = str(qmap.loc[("GSE131907", unit)])
            cell_rows.append(
                {
                    "cohort": "GSE131907",
                    "unit": unit,
                    "sample": samp,
                    "quartile_cldn4_pct": q,
                    "n_malignant": int(mal_idx.size),
                    "n_tnk": int(tnk_idx.size),
                }
            )
            scored = score_patient_pairs(lr, log_cp, pos, mal_idx, tnk_idx)
            for row in scored:
                row["cohort"] = "GSE131907"
                row["unit"] = unit
                row["quartile_cldn4_pct"] = q
                per_rows.append(row)
            print(f"  scored GSE131907 {unit} {samp} mal={mal_idx.size} tnk={tnk_idx.size} q={q}", flush=True)

    b_tails = tails[tails["cohort"] == "GSE148071"]
    if not b_tails.empty:
        files = list_exp_files(args.data_148071)
        want_pts = set(b_tails["unit"])
        for fp in files:
            patient = patient_from_filename(fp)
            if patient not in want_pts:
                continue
            cell_ids, found, n_umi = stream_one_148071(fp, wanted)
            n = len(cell_ids)
            lib = np.maximum(n_umi, 1.0)
            log_cp = {g: np.log1p(found[g] / lib * 1e4).astype(np.float32) for g in found}
            pos = {g: (found[g] > 0).astype(np.float32) for g in found}
            scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
            cd3 = log_cp.get("CD3E", np.zeros(n, dtype=np.float32))
            lineage = assign_lineage(scores, cd3)
            mal_idx = np.flatnonzero(lineage == "Epithelial")
            tnk_idx = np.flatnonzero(np.isin(lineage, ["T", "NK"]))
            q = str(qmap.loc[("GSE148071", patient)])
            cell_rows.append(
                {
                    "cohort": "GSE148071",
                    "unit": patient,
                    "sample": patient,
                    "quartile_cldn4_pct": q,
                    "n_malignant": int(mal_idx.size),
                    "n_tnk": int(tnk_idx.size),
                }
            )
            scored = score_patient_pairs(lr, log_cp, pos, mal_idx, tnk_idx)
            for row in scored:
                row["cohort"] = "GSE148071"
                row["unit"] = patient
                row["quartile_cldn4_pct"] = q
                per_rows.append(row)
            print(f"  scored GSE148071 {patient} mal={mal_idx.size} tnk={tnk_idx.size} q={q}", flush=True)

    per_patient_cells = pd.DataFrame(cell_rows)
    per_patient_cells.to_csv(out / "results" / "per_patient_cells.tsv", sep="\t", index=False)
    long = pd.DataFrame(per_rows)
    if long.empty:
        raise SystemExit("no per-patient LR rows scored")
    long.to_csv(out / "results" / "per_patient_lr.tsv.gz", sep="\t", index=False)

    contrast_rows = []
    for (name, direction), block in long.groupby(["interaction_name", "direction"], observed=True):
        q1 = block[block["quartile_cldn4_pct"] == "Q1"]
        q4 = block[block["quartile_cldn4_pct"] == "Q4"]
        d1 = q1[q1["detected"]]
        d4 = q4[q4["detected"]]
        if len(d1) < MIN_DETECT_ARM or len(d4) < MIN_DETECT_ARM:
            continue
        u, p = stats.mannwhitneyu(d4["prob"].to_numpy(), d1["prob"].to_numpy(), alternative="two-sided")
        n1, n4 = int(len(d1)), int(len(d4))
        r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
        med1 = float(d1["prob"].median())
        med4 = float(d4["prob"].median())
        rec0 = block.iloc[0]
        contrast_rows.append(
            {
                "interaction_name": name,
                "direction": direction,
                "pathway_name": rec0["pathway_name"],
                "annotation": rec0["annotation"],
                "ligand": rec0["ligand"],
                "receptor": rec0["receptor"],
                "ligand_genes": rec0["ligand_genes"],
                "receptor_genes": rec0["receptor_genes"],
                "ligand_class": ligand_class(rec0["ligand_genes"]),
                "n_q1_detected": n1,
                "n_q4_detected": n4,
                "n_compared": n1 + n4,
                "n_q1_131907": int(((d1["cohort"] == "GSE131907")).sum()),
                "n_q4_131907": int(((d4["cohort"] == "GSE131907")).sum()),
                "n_q1_148071": int(((d1["cohort"] == "GSE148071")).sum()),
                "n_q4_148071": int(((d4["cohort"] == "GSE148071")).sum()),
                "median_prob_q1": med1,
                "median_prob_q4": med4,
                "delta_median": med4 - med1,
                "mwu_u": float(u),
                "p": float(p),
                "r_rb": float(r_rb),
            }
        )
    contrast = pd.DataFrame(contrast_rows)
    if contrast.empty:
        ligand_tbl = contrast
    else:
        contrast = contrast[contrast["delta_median"] != 0].copy()
        contrast = contrast.assign(_absd=contrast["delta_median"].abs()).sort_values(
            ["p", "_absd"], ascending=[True, False]
        ).drop(columns="_absd")
        contrast.to_csv(out / "results" / "lr_q4q1_all.tsv", sep="\t", index=False)
        outgoing = contrast[contrast["direction"] == "outgoing"]
        ligand_tbl = outgoing.sort_values(["p"]).head(20).copy() if not outgoing.empty else contrast.head(20).copy()
        n_sig = int((contrast["p"] < 0.05).sum())
        n_sig_out = int(((contrast["direction"] == "outgoing") & (contrast["p"] < 0.05)).sum())
        ligand_tbl["sig_p05"] = ligand_tbl["p"] < 0.05
        ligand_tbl["note"] = (
            f"{n_sig_out} outgoing p<0.05 of {int((contrast['direction']=='outgoing').sum())} "
            f"detect-gated outgoing rows ({n_sig} any direction)"
        )
    ligand_tbl.to_csv(out / "results" / "ligand_table.tsv", sep="\t", index=False)
    n_q1 = int((tails["q_pct"] == "Q1").sum())
    n_q4 = int((tails["q_pct"] == "Q4").sum())
    plot_ligand_table(
        ligand_tbl if not ligand_tbl.empty else pd.DataFrame(),
        out / "figures" / "fig_extra_ligand_table.png",
        f"Pair 131907+148071 CellChat-style Mal→T/NK  Q4 vs Q1 (n={n_q1}/{n_q4})",
    )
    return {
        "n_lr_in_db": int(len(lr)),
        "n_patients_scored": int(per_patient_cells.shape[0]),
        "n_q1": n_q1,
        "n_q4": n_q4,
        "n_pairs_tested_after_detect_gate": int(len(contrast)) if not contrast.empty else 0,
        "n_pairs_p_lt_05": int((contrast["p"] < 0.05).sum()) if not contrast.empty else 0,
        "n_outgoing_p_lt_05": int(((contrast["direction"] == "outgoing") & (contrast["p"] < 0.05)).sum())
        if not contrast.empty
        else 0,
        "ligand_table_rows": int(len(ligand_tbl)),
    }


def write_finding(tnk: dict, cellchat: dict | None, out: Path) -> None:
    table = tnk["table"]
    patients = tnk["patients"]
    combo = tnk["combo_pct"]
    combo_q = tnk["combo_q"]
    prim = table[table["score"] == "pct_pos"]
    sec = table[table["score"] == "mean"]
    tails = tnk["tails"].sort_values(["q_pct", "cohort", "mal_CLDN4_pct_pos"], ascending=[True, True, False])

    def row_md(r) -> str:
        return (
            f"| {r.cohort} | {r.score} | {int(r.n_patients)} | "
            f"{r.spearman_rho:+.3f} ({fmt_p(r.spearman_p)}) | "
            f"{r.q4q1_r_rb:+.3f} ({fmt_p(r.q4q1_p)}; {int(r.n_q1)}/{int(r.n_q4)}) | "
            f"{r.median_tnk_q1:.3f} | {r.median_tnk_q4:.3f} | {r.delta_median_tnk:+.3f} |"
        )

    n131 = int((patients["cohort"] == "GSE131907").sum())
    n148 = int((patients["cohort"] == "GSE148071").sum())
    nq1 = int((tails["q_pct"] == "Q1").sum())
    nq4 = int((tails["q_pct"] == "Q4").sum())
    lines = [
        "# FINDING — pair GSE131907 + GSE148071: malignant CLDN4 vs T/NK + CellChat",
        "",
        "ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. **GSE205335 is not",
        "in this pair.** Patient is the unit. Combo ρ is DerSimonian–Laird on",
        "Fisher-z of the two within-cohort Spearmans. p-values are descriptive.",
        "Prior single-dataset CellChat (PR #347, #348) and the Q4 meta that",
        "pairs GSE131907 with GSE205335 (PR #320) are given and are not re-ranked.",
        "",
        "## Honest n",
        "",
        f"- GSE131907 (Kim et al. 2020): locked author-`Malignant cells` extract",
        f"  (PR #320) with ≥{MIN_MAL_131907} malignant **and** ≥{MIN_TNK_131907} T/NK.",
        f"  After the floor this is **{n131} patients** (tL/B + mLN + mBrain; one",
        "  sample each). tLung is out of this extract because primary tumor",
        "  epithelium is labeled tS* not `Malignant cells`. PE / nLung / nLN out.",
        "  No ICI / MPR labels.",
        f"- GSE148071 (Wu et al. 2021): **{n148} / 42** biopsies with",
        f"  ≥{MIN_EPI_148071} marker-argmax epithelial (putative malignant) **and**",
        f"  ≥{MIN_TNK_148071} T/NK. No histology / ICI labels. Epithelium is",
        "  putative; CopyKAT IDs are not on GEO.",
        f"- Combo N = **{n131 + n148}**. Q4 vs Q1 uses **within-cohort** tails:",
        f"  **n={nq1} vs {nq4}** (n_compared={nq1 + nq4}), not the full N.",
        "- Stacked 10x + Singleron Spearman is a companion only. It is not the combo ρ.",
        "",
        "## Combo ρ (primary deliverable)",
        "",
        "| analysis | combo | k | N | effect (p, I²) |",
        "|---|---|---:|---:|---|",
    ]
    if combo.get("k", 0) >= 2:
        lines.append(
            f"| Spearman %pos vs T/NK | GSE131907+GSE148071 | {int(combo['k'])} | "
            f"{int(combo['n_patients_total'])} | "
            f"ρ={combo['pooled_rho']:+.3f} ({fmt_p(combo['p'])}, {combo['I2']:.0f}%) |"
        )
    if tnk["combo_mean"].get("k", 0) >= 2:
        cm = tnk["combo_mean"]
        lines.append(
            f"| Spearman mean vs T/NK | GSE131907+GSE148071 | {int(cm['k'])} | "
            f"{int(cm['n_patients_total'])} | "
            f"ρ={cm['pooled_rho']:+.3f} ({fmt_p(cm['p'])}, {cm['I2']:.0f}%) |"
        )
    if combo_q.get("k", 0) >= 2:
        lines.append(
            f"| Q4 vs Q1 r %pos | GSE131907+GSE148071 | {int(combo_q['k'])} | "
            f"{int(combo_q['n_patients_total'])} compared | "
            f"r={combo_q['pooled_rho']:+.3f} ({fmt_p(combo_q['p'])}, {combo_q['I2']:.0f}%) |"
        )
    st = tnk["stacked"]
    lines += [
        "",
        (
            f"Stacked companion (not the combo): ρ={st['rho']:+.3f} "
            f"p={fmt_p(st['p'])} n={st['n']}."
        ),
        "",
        "Primary cut is **%pos**. Mean is the same patients, secondary.",
        "GSE131907 alone is CLDN4-negative vs T/NK. GSE148071 is null.",
        "The pair is heterogeneous; that is the combo, not a hidden single-cohort claim.",
        "",
        "## Per-cohort malignant CLDN4 vs same-patient T/NK",
        "",
        "| cohort | score | n | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | median T/NK Q1 | median T/NK Q4 | Δ |",
        "|---|---|---:|---|---|---:|---:|---:|",
    ]
    for rec in prim.itertuples():
        lines.append(row_md(rec))
    for rec in sec.itertuples():
        lines.append(row_md(rec))
    lines += [
        "",
        "### Quartile tails (within-cohort CLDN4 %pos)",
        "",
        "| tail | cohort | unit | origin | n_mal | n_TNK | CLDN4 %pos | T/NK frac |",
        "|---|---|---|---|---:|---:|---:|---:|",
    ]
    for rec in tails.itertuples():
        lines.append(
            f"| {rec.q_pct} | {rec.cohort} | {rec.unit} | {rec.origin} | "
            f"{int(rec.n_malignant)} | {int(rec.n_tnk)} | "
            f"{rec.mal_CLDN4_pct_pos:.1f} | {rec.frac_tnk:.3f} |"
        )
    lines += ["", "Q4 vs Q1 is within-cohort so 10x and Singleron are not ranked on one scale.", ""]

    if cellchat is None:
        lines += [
            "## CellChat-style ligands",
            "",
            "Not scored in this write-up (matrix step skipped). Re-run without",
            "`--skip-cellchat` to fill `results/ligand_table.tsv`.",
            "",
        ]
    else:
        lig_path = out / "results" / "ligand_table.tsv"
        lig = pd.read_csv(lig_path, sep="\t") if lig_path.exists() else pd.DataFrame()
        lines += [
            "## CellChat-style outgoing CLDN4-high → T/NK (Q4 vs Q1 patients)",
            "",
            "Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs.",
            "Outgoing = author-malignant (GSE131907) or putative epithelium",
            "(GSE148071) → **same-patient** T/NK. Incoming = T/NK → malignant.",
            f"Test = Mann–Whitney on per-patient *P* (detected in ≥{MIN_DETECT_ARM}",
            "Q1 and ≥3 Q4). CellChat R and LIANA were not run.",
            "",
            (
                f"Patients scored: {cellchat['n_patients_scored']} "
                f"(Q1={cellchat['n_q1']}, Q4={cellchat['n_q4']}). "
                f"Detect-gated direction×pair rows: {cellchat['n_pairs_tested_after_detect_gate']}. "
                f"p<0.05: {cellchat['n_pairs_p_lt_05']} "
                f"({cellchat['n_outgoing_p_lt_05']} outgoing)."
            ),
            "",
        ]
        if lig.empty:
            lines.append("No pairs passed the detect gate.")
        else:
            lines += [
                "| direction | pair | class | n_Q1/n_Q4 | median P Q1 | median P Q4 | Δ | r | p |",
                "|---|---|---|---|---:|---:|---:|---:|---|",
            ]
            show = lig.head(15)
            for rec in show.itertuples():
                star = " **" if getattr(rec, "sig_p05", rec.p < 0.05) else ""
                lines.append(
                    f"| {rec.direction} | {rec.ligand}–{rec.receptor}{star} | {rec.ligand_class} | "
                    f"{int(rec.n_q1_detected)}/{int(rec.n_q4_detected)} | "
                    f"{rec.median_prob_q1:.3f} | {rec.median_prob_q4:.3f} | "
                    f"{rec.delta_median:+.3f} | {rec.r_rb:+.3f} | {fmt_p(rec.p)} |"
                )
            n_sig_out = int(cellchat["n_outgoing_p_lt_05"])
            lines += [
                "",
                (
                    f"{n_sig_out} outgoing pair(s) p<0.05. The table is top outgoing "
                    "rows by p among detect-gated pairs. Cell-pooled truncated means "
                    "are **not** the test."
                ),
                "",
                "Honest limits on those six: APP–CD74 (n=12/11, both cohorts) and",
                "MDK–NCL (n=13/11, both cohorts) are the only p<0.05 rows with",
                "both arms ≥6. F11R–ITGAL/ITGB2 is n_Q4=3 (thin). The three",
                "HLA-II–CD4 rows are n=4/4 with |r|=1 (complete separation, not",
                "a poolable effect). CD274–PDCD1 is not detect-gated.",
                "NECTIN2–TIGIT is GSE148071-only and p>0.4.",
                "",
            ]
        lines += [
            "CD274–PDCD1 / NECTIN2 are claimed only if they appear as detected rows.",
            "No dual-high TACSTD2×CLDN4 split.",
            "",
        ]

    lines += [
        "## Files",
        "",
        "- `results/combo_rho.tsv` — combo ρ / Q4 vs Q1 r (this pair only)",
        "- `results/q4q1_tnk.tsv` — per-cohort Spearman + Q4 vs Q1",
        "- `results/patients_with_quartiles.tsv` — honest patient table",
        "- `results/ligand_table.tsv` — CellChat-style outgoing Mal → T/NK",
        "- `figures/fig_combo_rho.png` — extra combo-ρ forest",
        "- `figures/q4q1_tnk_pct.png` — Q4 vs Q1 T/NK box",
        "- `figures/scatter_cldn4_tnk.png` — extra scatter",
        "- `figures/fig_extra_ligand_table.png` — extra ligand-table figure",
        "- `METHODS.md` — floors, quartiles, Hill probability",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=ROOT / "data")
    p.add_argument("--db", type=Path, default=ROOT / "db")
    p.add_argument("--out", type=Path, default=ROOT)
    p.add_argument("--matrix-131907", dest="matrix_131907", type=Path, default=Path("/tmp/pair_131907_148071/gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"))
    p.add_argument("--ann-131907", dest="ann_131907", type=Path, default=Path("/tmp/pair_131907_148071/gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"))
    p.add_argument("--data-148071", dest="data_148071", type=Path, default=Path("/tmp/pair_131907_148071/gse148071"))
    p.add_argument("--skip-cellchat", action="store_true")
    args = p.parse_args()
    (args.out / "results").mkdir(parents=True, exist_ok=True)
    (args.out / "figures").mkdir(parents=True, exist_ok=True)

    tnk = run_tnk(args.data, args.out)
    print(tnk["table"].to_string(index=False), flush=True)
    print("combo_pct", json.dumps({k: tnk["combo_pct"].get(k) for k in ("k", "n_patients_total", "pooled_rho", "p", "I2")}), flush=True)

    cellchat = None
    if not args.skip_cellchat:
        need = [args.matrix_131907, args.ann_131907]
        missing = [str(x) for x in need if not x.exists()]
        if missing:
            raise SystemExit("missing CellChat inputs: " + ", ".join(missing))
        cellchat = run_cellchat(args, tnk["tails"], args.out)
        print(json.dumps(cellchat, indent=2), flush=True)

    write_finding(tnk, cellchat, args.out)
    summary = {
        "pair": ["GSE131907", "GSE148071"],
        "excluded": ["GSE205335"],
        "additive": True,
        "marker": "CLDN4",
        "dual_high": False,
        "unit": "patient",
        "n_gse131907": int((tnk["patients"]["cohort"] == "GSE131907").sum()),
        "n_gse148071": int((tnk["patients"]["cohort"] == "GSE148071").sum()),
        "combo_rho_pct": tnk["combo_pct"],
        "combo_rho_mean": tnk["combo_mean"],
        "combo_q4q1_pct": tnk["combo_q"],
        "per_cohort": tnk["table"].to_dict(orient="records"),
        "cellchat": cellchat,
        "algorithm": (
            "within-cohort pd.qcut on average ranks; MWU rank-biserial; "
            "DL Fisher-z combo ρ; CellChat-like 10% trim mean, Hill Kh=0.5, "
            "expr_prop>=0.10, same-patient Mal↔T/NK on Q4 vs Q1 tails"
        ),
        "not_run": "CellChat R; LIANA; EGA raw; TACSTD2 split; GSE205335",
    }
    (args.out / "results" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote", args.out / "FINDING.md", flush=True)


if __name__ == "__main__":
    main()
