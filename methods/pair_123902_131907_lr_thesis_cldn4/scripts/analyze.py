#!/usr/bin/env python3
"""ADDITIVE CLDN4-only thesis-aligned ligand test: GSE123902 + GSE131907.

The pair %pos Spearman is taken as given from PR #459 and is not re-audited:
  n=34, ρ=−0.575, Q4 vs Q1 r=−0.700.

No dual-high. Do not add GSE148071.
Patient/donor is the unit (GSE123902 donor; GSE131907 locked sample).
Both pre-specified families are scored. Done when both family tables exist.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lib import (  # noqa: E402
    EPI_MARKERS,
    EXPR_PROP,
    FAMILIES,
    KH,
    MALIG_SUBTYPES,
    MIN_CELLS_ARM,
    MIN_CELLS_COMP,
    MIN_MAL_Q4,
    TNK_MARKERS,
    TRIM,
    TUMOR_ORIGINS,
    dersimonian_laird,
    fmt_num,
    fmt_p,
    highend_split,
    load_lr,
    score_outgoing,
    thesis_genes,
    unit_delta_stats,
)

LOCKED = ROOT / "data"
DB = ROOT / "db"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

# PR #459 given pair — not re-audited
GIVEN = {
    "combo": "GSE123902+GSE131907",
    "score": "pct",
    "n": 34,
    "rho": -0.575072194708154,
    "p": 0.0005276801167933881,
    "I2": 0.0,
    "q4_r": -0.7,
    "q4_p": 0.011655011655011656,
    "n_q1": 10,
    "n_q4": 8,
    "source": "PR459",
    "re_audited": False,
    "member_rhos": {"GSE123902": -0.659, "GSE131907": -0.522},
    "member_ns": {"GSE123902": 13, "GSE131907": 21},
}
GIVEN_SINGLES = [
    {"cohort": "GSE123902", "n": 13, "rho": -0.659, "p": 0.014, "score": "pct"},
    {"cohort": "GSE131907", "n": 21, "rho": -0.522, "p": 0.015, "score": "pct"},
    {"cohort": "combo (given)", "n": 34, "rho": -0.575, "p": 0.00053, "score": "pct"},
]


def to_log_pos(extracted: dict[str, np.ndarray], library: np.ndarray):
    lib = np.maximum(library, 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}
    return log_cp, pos


def load_locked() -> dict[str, pd.DataFrame]:
    out = {}
    d = pd.read_csv(LOCKED / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"].copy()
    el["unit"] = "GSE123902"
    el["patient_id"] = el["patient"].astype(str)
    el["cldn4_given"] = el["mal_CLDN4_pct"].astype(float)
    el["frac_tnk_given"] = el["frac_tnk"].astype(float)
    out["GSE123902"] = el

    s = pd.read_csv(LOCKED / "GSE131907_samples.tsv", sep="\t")
    el2 = s[(s["n_malignant"] > 0) & (s["n_tnk"] > 0)].copy()
    el2["unit"] = "GSE131907"
    el2["patient_id"] = el2["sample"].astype(str)
    el2["cldn4_given"] = el2["mal_CLDN4_pct"].astype(float)
    el2["frac_tnk_given"] = el2["frac_tnk"].astype(float)
    out["GSE131907"] = el2
    return out


def _marker_mal_tnk(extracted: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    def col(g: str) -> np.ndarray:
        if g not in extracted:
            return np.zeros(n, dtype=float)
        return extracted[g]

    epi = np.zeros(n, dtype=bool)
    for g in EPI_MARKERS:
        epi |= col(g) > 0
    ptprc = col("PTPRC")
    mal = epi & (ptprc == 0)
    tnk = np.zeros(n, dtype=bool)
    for g in TNK_MARKERS:
        tnk |= col(g) > 0
    tnk = tnk & (~mal)
    return mal, tnk


def _score_one_123902_csv(handle, wanted: set[str]):
    df = pd.read_csv(handle, index_col=0)
    df.columns = [str(c).upper() for c in df.columns]
    gene_set = set(df.columns)
    libs = df.sum(axis=1).to_numpy(dtype=np.float64)
    extracted = {g: df[g].to_numpy(dtype=np.float32) for g in wanted if g in df.columns}
    return extracted, libs, gene_set


def load_gse123902(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE123902 ====", flush=True)
    keep = set(locked["patient_id"].astype(str))
    tar_path = raw / "GSE123902" / "GSE123902_RAW.tar"
    chosen: dict[str, str] = {}
    members = []
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".csv.gz"):
                continue
            name = Path(m.name).name
            parts = name.split("_")
            donor = parts[2] if len(parts) > 2 else name
            tissue = (
                "NORMAL"
                if "NORMAL" in name
                else (
                    "METASTASIS"
                    if "METASTASIS" in name
                    else ("PRIMARY" if "PRIMARY" in name else "OTHER")
                )
            )
            if tissue not in {"PRIMARY", "METASTASIS"}:
                continue
            if donor not in keep:
                continue
            members.append((donor, tissue, m, name))
        members.sort(key=lambda x: (x[0], 0 if x[1] == "PRIMARY" else 1))
        chunks = []
        gene_union: set[str] = set()
        for donor, tissue, m, name in members:
            if donor in chosen:
                continue
            chosen[donor] = tissue
            raw_f = gzip.GzipFile(fileobj=tf.extractfile(m))
            extracted, libs, genes = _score_one_123902_csv(raw_f, wanted)
            gene_union |= genes
            n = int(libs.size)
            chunks.append((donor, extracted, libs))
            print(f"  {name}: cells={n} stored={len(extracted)} donor={donor} {tissue}", flush=True)
    if not chunks:
        raise SystemExit("GSE123902: no locked tumor matrices read")
    units = []
    for donor, extracted, libs in chunks:
        n = int(libs.size)
        log_cp, pos = to_log_pos(extracted, libs)
        mal, tnk = _marker_mal_tnk(extracted, n)
        cldn4 = log_cp.get("CLDN4", np.zeros(n, dtype=np.float32))
        units.append(
            {
                "cohort": "GSE123902",
                "patient_id": donor,
                "n_cells": n,
                "n_mal": int(mal.sum()),
                "n_tnk": int(tnk.sum()),
                "mal": mal,
                "tnk": tnk,
                "cldn4": cldn4,
                "log_cp": log_cp,
                "pos": pos,
                "genes": set(extracted),
            }
        )
        print(f"  unit {donor}: n={n} mal={int(mal.sum())} tnk={int(tnk.sum())}", flush=True)
    return units, gene_union


def stream_131907(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in wanted:
                found[gene] = arr
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  stream genes={n_streamed} stored={len(found)}", flush=True)
    print(f"stream done genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    return cell_ids, found, n_umi, n_streamed


def load_gse131907(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE131907 ====", flush=True)
    keep = set(locked["patient_id"].astype(str))
    ann_path = raw / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat_path = raw / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t")
    cell_ids, found, n_umi, n_genes = stream_131907(mat_path, wanted)
    cell_index = {c: i for i, c in enumerate(cell_ids)}
    # annotation Index usually matches matrix column
    key = "Index" if "Index" in ann.columns else ann.columns[0]
    ann = ann.copy()
    ann["_i"] = ann[key].map(cell_index)
    ann = ann[ann["_i"].notna()].copy()
    ann["_i"] = ann["_i"].astype(int)
    sample_col = "Sample" if "Sample" in ann.columns else "sample"
    origin_col = "Sample_Origin" if "Sample_Origin" in ann.columns else None
    type_col = "Cell_type" if "Cell_type" in ann.columns else None
    sub_col = "Cell_subtype" if "Cell_subtype" in ann.columns else None

    mal_mask = np.zeros(len(cell_ids), dtype=bool)
    tnk_mask = np.zeros(len(cell_ids), dtype=bool)
    sample_of = np.array([""] * len(cell_ids), dtype=object)
    for rec in ann.itertuples(index=False):
        i = int(getattr(rec, "_i"))
        sample_of[i] = str(getattr(rec, sample_col))
        subtype = str(getattr(rec, sub_col)) if sub_col else ""
        ctype = str(getattr(rec, type_col)) if type_col else ""
        origin = str(getattr(rec, origin_col)) if origin_col else ""
        if subtype in MALIG_SUBTYPES and (origin in TUMOR_ORIGINS or origin_col is None):
            mal_mask[i] = True
        if ctype in {"T lymphocytes", "NK cells"}:
            tnk_mask[i] = True

    log_cp, pos = to_log_pos(found, n_umi)
    cldn4 = log_cp.get("CLDN4", np.zeros(len(cell_ids), dtype=np.float32))
    units = []
    for sample in sorted(keep):
        idx = np.flatnonzero(sample_of == sample)
        if idx.size == 0:
            print(f"  MISSING sample {sample}", flush=True)
            continue
        mal = mal_mask[idx]
        tnk = tnk_mask[idx]
        # restrict arrays to this sample
        extracted = {g: found[g][idx] for g in found}
        libs = n_umi[idx]
        log_u, pos_u = to_log_pos(extracted, libs)
        units.append(
            {
                "cohort": "GSE131907",
                "patient_id": sample,
                "n_cells": int(idx.size),
                "n_mal": int(mal.sum()),
                "n_tnk": int(tnk.sum()),
                "mal": mal,
                "tnk": tnk,
                "cldn4": cldn4[idx],
                "log_cp": log_u,
                "pos": pos_u,
                "genes": set(found),
            }
        )
        print(
            f"  unit {sample}: n={int(idx.size)} mal={int(mal.sum())} tnk={int(tnk.sum())}",
            flush=True,
        )
    return units, set(found), n_genes


def score_unit(unit: dict, lr: pd.DataFrame, mode: str) -> dict | None:
    mal_idx = np.flatnonzero(unit["mal"])
    tnk_idx = np.flatnonzero(unit["tnk"])
    if tnk_idx.size < MIN_CELLS_ARM:
        return None
    split = highend_split(unit["cldn4"], mal_idx, mode)
    if split is None:
        return None
    hi, lo = split
    hi_rows = score_outgoing(lr, unit["log_cp"], unit["pos"], hi, tnk_idx)
    lo_rows = score_outgoing(lr, unit["log_cp"], unit["pos"], lo, tnk_idx)
    if not hi_rows or not lo_rows:
        return None
    hi_map = {r["interaction_name"]: r for r in hi_rows}
    lo_map = {r["interaction_name"]: r for r in lo_rows}
    pair_rows = []
    for name in hi_map:
        h = hi_map[name]
        l = lo_map[name]
        detected = bool(h["detected"] or l["detected"])
        pair_rows.append(
            {
                "cohort": unit["cohort"],
                "patient_id": unit["patient_id"],
                "split": mode,
                "interaction_name": name,
                "ligand": h["ligand"],
                "receptor": h["receptor"],
                "prob_high": h["prob"],
                "prob_low": l["prob"],
                "delta": h["prob"] - l["prob"],
                "detected": detected,
                "detected_high": bool(h["detected"]),
                "detected_low": bool(l["detected"]),
                "n_mal_high": int(hi.size),
                "n_mal_low": int(lo.size),
                "n_tnk": int(tnk_idx.size),
                "n_mal": int(mal_idx.size),
            }
        )
    return {
        "cohort": unit["cohort"],
        "patient_id": unit["patient_id"],
        "split": mode,
        "n_mal": int(mal_idx.size),
        "n_mal_high": int(hi.size),
        "n_mal_low": int(lo.size),
        "n_tnk": int(tnk_idx.size),
        "pairs": pair_rows,
    }


def meta_pair(df: pd.DataFrame) -> dict:
    """Patient-level ΔP meta. Honest n = number of units with a finite delta."""
    sub = df[df["detected"]].copy() if "detected" in df.columns else df.copy()
    if sub.empty:
        return {
            "n": 0,
            "n_123902": 0,
            "n_131907": 0,
            "mean": np.nan,
            "median": np.nan,
            "se": np.nan,
            "wilcoxon_p": np.nan,
            "meta_p": np.nan,
            "i2": np.nan,
            "n_pos": 0,
            "n_neg": 0,
        }
    stats = unit_delta_stats(sub["delta"])
    # per-unit variance is not known; use sample variance for DL (same n)
    if stats["n"] >= 2 and np.isfinite(stats["sd"]) and stats["sd"] > 0:
        v = np.full(stats["n"], stats["sd"] ** 2)
        dl = dersimonian_laird(sub["delta"].to_numpy(float), v)
        meta_p = dl["p"]
        i2 = dl["i2"]
        mean = dl["mean"]
        se = dl["se"]
    else:
        meta_p, i2, mean, se = stats["wilcoxon_p"], np.nan, stats["mean"], stats["se"]
    return {
        "n": stats["n"],
        "n_123902": int((sub["cohort"] == "GSE123902").sum()),
        "n_131907": int((sub["cohort"] == "GSE131907").sum()),
        "mean": mean,
        "median": stats["median"],
        "se": se,
        "wilcoxon_p": stats["wilcoxon_p"],
        "meta_p": meta_p,
        "i2": i2,
        "n_pos": int((sub["delta"] > 0).sum()),
        "n_neg": int((sub["delta"] < 0).sum()),
    }


def agrees(expect: str, mean_delta: float, n: int, p: float) -> str:
    if n < 3:
        return "thin"
    if not np.isfinite(mean_delta):
        return "undetected"
    if expect == "high>low":
        direction_ok = mean_delta > 0
    else:
        direction_ok = mean_delta < 0
    if direction_ok and np.isfinite(p) and p < 0.05:
        return "yes"
    if direction_ok:
        return "direction_only"
    if np.isfinite(p) and p < 0.05:
        return "opposite"
    return "no"


def family_table(long_df: pd.DataFrame, family_key: str, mode: str) -> pd.DataFrame:
    spec = FAMILIES[family_key]
    rows = []
    pair_long = []
    for rec in spec["pairs"]:
        name = rec["interaction_name"]
        sub = long_df[(long_df["interaction_name"] == name) & (long_df["split"] == mode)]
        m = meta_pair(sub)
        observed = (
            "high>low"
            if np.isfinite(m["mean"]) and m["mean"] > 0
            else ("low>high" if np.isfinite(m["mean"]) and m["mean"] < 0 else "flat/NA")
        )
        rows.append(
            {
                "family": family_key,
                "row": "pair",
                "interaction_name": name,
                "ligand": rec["ligand"],
                "receptor": rec["receptor"],
                "axis": rec["axis"],
                "thesis_expect": rec["thesis_expect"],
                "split": mode,
                "n_units": m["n"],
                "n_gse123902": m["n_123902"],
                "n_gse131907": m["n_131907"],
                "mean_delta": m["mean"],
                "median_delta": m["median"],
                "se": m["se"],
                "wilcoxon_p": m["wilcoxon_p"],
                "meta_p": m["meta_p"],
                "I2": m["i2"],
                "n_pos": m["n_pos"],
                "n_neg": m["n_neg"],
                "observed_direction": observed,
                "agrees_thesis": agrees(rec["thesis_expect"], m["mean"], m["n"], m["wilcoxon_p"]),
            }
        )
        if not sub.empty:
            pair_long.append(sub)

    # Family aggregate: per-unit mean ΔP of detected pairs, then patient test
    if pair_long:
        cat = pd.concat(pair_long, ignore_index=True)
        det = cat[cat["detected"]].copy()
        if not det.empty:
            fam_unit = (
                det.groupby(["cohort", "patient_id"], as_index=False)
                .agg(delta=("delta", "mean"), n_pairs=("delta", "size"))
            )
            fam_unit["detected"] = True
            m = meta_pair(fam_unit)
            observed = (
                "high>low"
                if np.isfinite(m["mean"]) and m["mean"] > 0
                else ("low>high" if np.isfinite(m["mean"]) and m["mean"] < 0 else "flat/NA")
            )
            rows.append(
                {
                    "family": family_key,
                    "row": "FAMILY_AGGREGATE",
                    "interaction_name": f"FAMILY:{family_key}",
                    "ligand": "|".join(sorted({r["ligand"] for r in spec["pairs"]})),
                    "receptor": "T/NK",
                    "axis": "family",
                    "thesis_expect": spec["thesis_expect"],
                    "split": mode,
                    "n_units": m["n"],
                    "n_gse123902": m["n_123902"],
                    "n_gse131907": m["n_131907"],
                    "mean_delta": m["mean"],
                    "median_delta": m["median"],
                    "se": m["se"],
                    "wilcoxon_p": m["wilcoxon_p"],
                    "meta_p": m["meta_p"],
                    "I2": m["i2"],
                    "n_pos": m["n_pos"],
                    "n_neg": m["n_neg"],
                    "observed_direction": observed,
                    "agrees_thesis": agrees(spec["thesis_expect"], m["mean"], m["n"], m["wilcoxon_p"]),
                }
            )
    return pd.DataFrame(rows)


def write_family_md_rows(df: pd.DataFrame) -> str:
    lines = [
        "| pair | axis | expect | n | mean ΔP | p_W | observed | agrees |",
        "|---|---|---|---:|---:|---|---|---|",
    ]
    for rec in df.itertuples(index=False):
        name = rec.interaction_name.replace("FAMILY:", "**FAMILY** ")
        lines.append(
            f"| {name} | {rec.axis} | {rec.thesis_expect} | {rec.n_units} | "
            f"{fmt_num(rec.mean_delta)} | {fmt_p(rec.wilcoxon_p)} | "
            f"{rec.observed_direction} | {rec.agrees_thesis} |"
        )
    return "\n".join(lines)


def fig_given(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    names = [r["cohort"] for r in GIVEN_SINGLES]
    rhos = [r["rho"] for r in GIVEN_SINGLES]
    ns = [r["n"] for r in GIVEN_SINGLES]
    colors = ["#4C78A8", "#4C78A8", "#F58518"]
    ax.barh(names[::-1], rhos[::-1], color=colors[::-1])
    ax.axvline(0, color="k", lw=0.8)
    for i, (n, r) in enumerate(zip(ns[::-1], rhos[::-1])):
        ax.text(r - 0.02 if r < 0 else r + 0.02, i, f"n={n}  ρ={r:.3f}", va="center", ha="right" if r < 0 else "left", fontsize=8)
    ax.set_xlabel("Spearman ρ (malignant CLDN4 %pos vs T/NK frac)")
    ax.set_title("PR #459 given pair — not re-audited")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_honest_n(cov: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    labels = ["given combo", "GSE123902 locked", "GSE131907 locked", "median LR", "tertile LR", "Q4 vs Q1 LR"]
    vals = [
        GIVEN["n"],
        int((cov["cohort"] == "GSE123902").sum()),
        int((cov["cohort"] == "GSE131907").sum()),
        int(cov["median_ok"].sum()),
        int(cov["tertile_ok"].sum()),
        int(cov["q4q1_ok"].sum()),
    ]
    ax.barh(labels[::-1], vals[::-1], color="#4C78A8")
    for i, v in enumerate(vals[::-1]):
        ax.text(v + 0.3, i, str(v), va="center", fontsize=8)
    ax.set_xlabel("honest n (patient/donor or locked sample)")
    ax.set_title("Given T/NK n is not the ligand n")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_n_per_unit(cov: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.2), sharey=False)
    for ax, cohort in zip(axes, ["GSE123902", "GSE131907"]):
        sub = cov[cov["cohort"] == cohort].sort_values("patient_id")
        x = np.arange(len(sub))
        ax.bar(x - 0.2, sub["n_mal"], width=0.4, label="malignant", color="#4C78A8")
        ax.bar(x + 0.2, sub["n_tnk"], width=0.4, label="T/NK", color="#F58518")
        ax.axhline(MIN_CELLS_ARM, color="k", ls="--", lw=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(sub["patient_id"], rotation=90, fontsize=6)
        ax.set_title(cohort)
        ax.set_ylabel("cells")
        ax.legend(fontsize=7)
    fig.suptitle("Per-unit malignant and T/NK floors (10-cell arm)")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_family_bars(fam: pd.DataFrame, title: str, path: Path) -> None:
    plot = fam[fam["row"] == "pair"].copy()
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    y = np.arange(len(plot))
    colors = ["#54A24B" if d > 0 else "#E45756" for d in plot["mean_delta"].fillna(0)]
    ax.barh(y, plot["mean_delta"].fillna(0), color=colors)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot["interaction_name"], fontsize=8)
    ax.set_xlabel("mean patient ΔP (CLDN4-high − CLDN4-low → T/NK)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_family_forest(long_df: pd.DataFrame, family_key: str, mode: str, path: Path) -> None:
    spec = FAMILIES[family_key]
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    ylabels = []
    means = []
    ses = []
    for rec in spec["pairs"]:
        sub = long_df[(long_df["interaction_name"] == rec["interaction_name"]) & (long_df["split"] == mode)]
        m = meta_pair(sub)
        ylabels.append(rec["interaction_name"])
        means.append(m["mean"] if np.isfinite(m["mean"]) else 0.0)
        ses.append(m["se"] if np.isfinite(m["se"]) else 0.0)
    y = np.arange(len(ylabels))
    ax.errorbar(means, y, xerr=ses, fmt="o", color="#4C78A8", capsize=3)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("mean ΔP ± SE (patient/donor unit)")
    ax.set_title(f"{spec['label']} — {mode}")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_family_patient(long_df: pd.DataFrame, family_key: str, mode: str, path: Path) -> None:
    spec = FAMILIES[family_key]
    names = {r["interaction_name"] for r in spec["pairs"]}
    sub = long_df[(long_df["interaction_name"].isin(names)) & (long_df["split"] == mode) & (long_df["detected"])]
    if sub.empty:
        return
    fam = sub.groupby(["cohort", "patient_id"], as_index=False).agg(delta=("delta", "mean"))
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0), sharex=True)
    for ax, cohort, color in zip(axes, ["GSE123902", "GSE131907"], ["#4C78A8", "#F58518"]):
        part = fam[fam["cohort"] == cohort].sort_values("delta")
        y = np.arange(len(part))
        ax.scatter(part["delta"], y, color=color, s=22)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(part["patient_id"], fontsize=6)
        ax.set_title(cohort)
        ax.set_xlabel("family mean ΔP")
    fig.suptitle(f"Per-unit family score — {family_key} / {mode}")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(
    cov: pd.DataFrame,
    fam_bar: pd.DataFrame,
    fam_ifn: pd.DataFrame,
    summary: dict,
) -> str:
    n_med = int(cov["median_ok"].sum())
    n_ter = int(cov["tertile_ok"].sum())
    n_q = int(cov["q4q1_ok"].sum())
    n123 = int((cov["cohort"] == "GSE123902").sum())
    n131 = int((cov["cohort"] == "GSE131907").sum())
    bar_agg = fam_bar[fam_bar["row"] == "FAMILY_AGGREGATE"]
    ifn_agg = fam_ifn[fam_ifn["row"] == "FAMILY_AGGREGATE"]
    bar_txt = "NA"
    ifn_txt = "NA"
    if not bar_agg.empty:
        r = bar_agg.iloc[0]
        bar_txt = f"n={int(r.n_units)} mean ΔP={fmt_num(r.mean_delta)} p_W={fmt_p(r.wilcoxon_p)} agrees={r.agrees_thesis}"
    if not ifn_agg.empty:
        r = ifn_agg.iloc[0]
        ifn_txt = f"n={int(r.n_units)} mean ΔP={fmt_num(r.mean_delta)} p_W={fmt_p(r.wilcoxon_p)} agrees={r.agrees_thesis}"

    md = f"""# FINDING — CLDN4-only thesis-aligned ligand test on GSE123902 + GSE131907

ADDITIVE. **CLDN4 only. No dual-high.** Patient/donor is the unit
(GSE123902 donor; GSE131907 locked sample from PR #459).
Do **not** add GSE148071. Do **not** re-audit the T/NK Spearman.

Thesis (already correct; this folder tests it on the given cut):

- CLDN4-high → more barrier/inhibitory outgoing to T/NK (F11R, NECTIN2–TIGIT, CDH1, LGALS9)
- CLDN4-low / KD-like → more IFN / T-recruit / MHC-I outgoing (CXCL9/10–CXCR3, CCL5, HLA–CD8)

The pair that **differs** is taken as given from PR #459:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE131907 | %pos | 34 | −0.575 (0.00053, 0%) | −0.700 (0.012; 10 vs 8) |

Singles (given, same PR): GSE123902 n=13 ρ=−0.659; GSE131907 n=21 ρ=−0.522.

This folder scores **both pre-specified families** as CellChat-style outgoing
Hill probabilities (Jin et al. 2021; 10% truncated mean, Kh={KH}, expr_prop≥{EXPR_PROP})
from CLDN4-high vs CLDN4-low malignant cells to same-unit T/NK.
CellChat R was not run. Pairs were not discovered.

## Honest n

Given combo n=34 is **not** the ligand n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (median: ≥{2 * MIN_CELLS_ARM}
malignant and ≥{MIN_CELLS_ARM}/arm plus ≥{MIN_CELLS_ARM} T/NK; Q4 vs Q1:
n_mal≥{MIN_MAL_Q4} and ≥{MIN_CELLS_ARM}/arm).

| gate | n | note |
|---|---:|---|
| Given combo (do not re-audit) | 34 | 13 donors + 21 samples |
| Locked GSE123902 donors | {n123} | marker-malignant; PRIMARY preferred over METASTASIS |
| Locked GSE131907 samples | {n131} | author Malignant cells / tS*; T lymphocytes + NK cells |
| Median-split LR | **{n_med}** | primary ligand n |
| Tertile extra | {n_ter} | exclusive arms |
| Q4 vs Q1 extra | {n_q} | n_mal≥40 |

GSE123902 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0).
GSE131907 malignant = author `Cell_subtype` in {{Malignant cells, tS1, tS2, tS3}}.
TACSTD2 is never a gate.

## Family 1 — barrier / inhibitory (expect CLDN4-high > low)

Primary table: `results/family_barrier_inhibitory.tsv`.

{write_family_md_rows(fam_bar)}

Family aggregate: {bar_txt}.

## Family 2 — IFN / T-recruit / MHC-I (expect CLDN4-low / KD-like > high)

Primary table: `results/family_ifn_recruit_mhci.tsv`.

{write_family_md_rows(fam_ifn)}

Family aggregate: {ifn_txt}.

## Extra figures

- `figures/fig_given_combo_rho.png` — given %pos ρ (not re-audited)
- `figures/fig_honest_n.png` — given n vs ligand-eligible n
- `figures/fig_n_per_unit.png` — per-unit malignant / T/NK floors
- `figures/fig_family_barrier_bars.png` / `fig_family_barrier_forest.png`
- `figures/fig_family_ifn_bars.png` / `fig_family_ifn_forest.png`
- `figures/fig_family_barrier_units.png` / `fig_family_ifn_units.png`

## What is not claimed

- The n=34 ρ=−0.575 / Q4 r=−0.700 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071 is not added.
- Cell-pooled permutations are not the test. Patient/donor is the unit.
- This is not a CellChat discovery screen. Only the two pre-specified families.

## Reproduce

```bash
python3 methods/pair_123902_131907_lr_thesis_cldn4/scripts/download.py
python3 methods/pair_123902_131907_lr_thesis_cldn4/scripts/analyze.py
```

Hill constants: trim={TRIM}, Kh={KH}, expr_prop={EXPR_PROP}.
"""
    return md


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("/tmp/pair_123902_131907_raw"))
    args = ap.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    locked = load_locked()
    print(
        f"locked GSE123902 n={len(locked['GSE123902'])} GSE131907 n={len(locked['GSE131907'])}",
        flush=True,
    )
    wanted = thesis_genes()
    names = {r["interaction_name"] for fam in FAMILIES.values() for r in fam["pairs"]}
    # load LR against union of wanted; matrix genes filled after load
    lr_all = load_lr(DB, names, matrix_genes=None)
    print(f"pre-specified pairs in CellChatDB: {len(lr_all)} / {len(names)}", flush=True)
    missing = names - set(lr_all["interaction_name"])
    if missing:
        print("WARNING missing interactions:", sorted(missing), flush=True)

    u123, genes123 = load_gse123902(args.raw, wanted, locked["GSE123902"])
    u131, genes131, n_streamed = load_gse131907(args.raw, wanted, locked["GSE131907"])
    units = u123 + u131
    gene_union = genes123 | genes131
    lr = load_lr(DB, names, matrix_genes=None)
    # drop pairs whose genes are absent from BOTH matrices
    keep_rows = []
    for rec in lr.itertuples(index=False):
        if all(g in gene_union for g in rec.ligand_genes + rec.receptor_genes):
            keep_rows.append(rec.Index if hasattr(rec, "Index") else None)
    # keep all DB pairs; undetected genes score as 0 via complex_from_maps
    print(f"units loaded: {len(units)}  genes stored 123902={len(genes123)} 131907={len(genes131)}", flush=True)

    coverage_rows = []
    long_rows = []
    for unit in units:
        rec = {
            "cohort": unit["cohort"],
            "patient_id": unit["patient_id"],
            "n_cells": unit["n_cells"],
            "n_mal": unit["n_mal"],
            "n_tnk": unit["n_tnk"],
            "median_ok": False,
            "tertile_ok": False,
            "q4q1_ok": False,
        }
        for mode in ("median", "tertile", "q4q1"):
            scored = score_unit(unit, lr, mode)
            if scored is None:
                continue
            rec[f"{mode}_ok"] = True
            rec[f"{mode}_n_high"] = scored["n_mal_high"]
            rec[f"{mode}_n_low"] = scored["n_mal_low"]
            long_rows.extend(scored["pairs"])
        coverage_rows.append(rec)

    cov = pd.DataFrame(coverage_rows)
    long_df = pd.DataFrame(long_rows)
    cov.to_csv(RESULTS / "patient_coverage.tsv", sep="\t", index=False)
    long_df.to_csv(RESULTS / "patient_lr_long.tsv", sep="\t", index=False)

    fam_bar = family_table(long_df, "barrier_inhibitory", "median")
    fam_ifn = family_table(long_df, "ifn_recruit_mhci", "median")
    fam_bar.to_csv(RESULTS / "family_barrier_inhibitory.tsv", sep="\t", index=False)
    fam_ifn.to_csv(RESULTS / "family_ifn_recruit_mhci.tsv", sep="\t", index=False)

    # extras: other splits
    extras = []
    for key in FAMILIES:
        for mode in ("median", "tertile", "q4q1"):
            extras.append(family_table(long_df, key, mode))
    pd.concat(extras, ignore_index=True).to_csv(RESULTS / "family_all_splits.tsv", sep="\t", index=False)

    fig_given(FIGURES / "fig_given_combo_rho")
    fig_honest_n(cov, FIGURES / "fig_honest_n")
    fig_n_per_unit(cov, FIGURES / "fig_n_per_unit")
    fig_family_bars(fam_bar, "Barrier/inhibitory outgoing ΔP (median)", FIGURES / "fig_family_barrier_bars")
    fig_family_bars(fam_ifn, "IFN/T-recruit/MHC-I outgoing ΔP (median)", FIGURES / "fig_family_ifn_bars")
    fig_family_forest(long_df, "barrier_inhibitory", "median", FIGURES / "fig_family_barrier_forest")
    fig_family_forest(long_df, "ifn_recruit_mhci", "median", FIGURES / "fig_family_ifn_forest")
    fig_family_patient(long_df, "barrier_inhibitory", "median", FIGURES / "fig_family_barrier_units")
    fig_family_patient(long_df, "ifn_recruit_mhci", "median", FIGURES / "fig_family_ifn_units")

    summary = {
        "given": GIVEN,
        "n_locked": {"GSE123902": int((cov.cohort == "GSE123902").sum()), "GSE131907": int((cov.cohort == "GSE131907").sum())},
        "n_lr": {
            "median": int(cov["median_ok"].sum()),
            "tertile": int(cov["tertile_ok"].sum()),
            "q4q1": int(cov["q4q1_ok"].sum()),
        },
        "family_tables": [
            "results/family_barrier_inhibitory.tsv",
            "results/family_ifn_recruit_mhci.tsv",
        ],
        "n_streamed_131907_genes": int(n_streamed),
        "re_audited_tnk_rho": False,
        "dual_high": False,
        "gse148071": False,
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    finding = write_finding(cov, fam_bar, fam_ifn, summary)
    (ROOT / "FINDING.md").write_text(finding)
    (RESULTS / "FINDING.md").write_text(finding)
    print("wrote", ROOT / "FINDING.md", flush=True)
    print("family tables:", RESULTS / "family_barrier_inhibitory.tsv", RESULTS / "family_ifn_recruit_mhci.tsv", flush=True)
    print(fam_bar.to_string(index=False), flush=True)
    print(fam_ifn.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
