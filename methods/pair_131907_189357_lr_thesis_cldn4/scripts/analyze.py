#!/usr/bin/env python3
"""ADDITIVE CLDN4-only thesis-aligned ligand test: GSE131907 + GSE189357.

The pair %pos Spearman is taken as given from PR #459 and is not re-audited:
  n=30, ρ=−0.542, Q4 vs Q1 r=−0.619.

No dual-high. Do not add GSE148071.
Patient is the unit (GSE131907 locked sample; GSE189357 patient TD1–TD9).
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
    "combo": "GSE131907+GSE189357",
    "score": "pct",
    "n": 30,
    "rho": -0.5424924995584873,
    "p": 0.002910690018538355,
    "I2": 0.0,
    "q4_r": -0.6190476190476191,
    "q4_p": 0.04178321678321679,
    "n_q1": 9,
    "n_q4": 7,
    "source": "PR459",
    "re_audited": False,
    "member_rhos": {"GSE131907": -0.522, "GSE189357": -0.600},
    "member_ns": {"GSE131907": 21, "GSE189357": 9},
}
GIVEN_SINGLES = [
    {"cohort": "GSE131907", "n": 21, "rho": -0.522, "p": 0.015, "score": "pct"},
    {"cohort": "GSE189357", "n": 9, "rho": -0.600, "p": 0.088, "score": "pct"},
    {"cohort": "combo (given)", "n": 30, "rho": -0.542, "p": 0.00291, "score": "pct"},
]


def to_log_pos(extracted: dict[str, np.ndarray], library: np.ndarray):
    lib = np.maximum(library, 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}
    return log_cp, pos


def load_locked() -> dict[str, pd.DataFrame]:
    out = {}
    s = pd.read_csv(LOCKED / "GSE131907_samples.tsv", sep="\t")
    el = s[
        s["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])
        & (s["n_malignant"] >= 20)
    ].copy()
    el["unit"] = "GSE131907"
    el["patient_id"] = el["sample"].astype(str)
    el["cldn4_given"] = el["mal_CLDN4_pct"].astype(float)
    el["frac_tnk_given"] = el["frac_tnk"].astype(float)
    out["GSE131907"] = el

    d = pd.read_csv(LOCKED / "GSE189357_marker_units.tsv", sep="\t")
    el2 = d[d["eligible"].astype(str).str.lower() == "true"].copy()
    el2["unit"] = "GSE189357"
    el2["patient_id"] = el2["patient"].astype(str)
    el2["cldn4_given"] = el2["mal_CLDN4_pct"].astype(float)
    el2["frac_tnk_given"] = el2["frac_tnk"].astype(float)
    out["GSE189357"] = el2
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


def stream_131907(matrix_path: Path, wanted: set[str], cache_dir: Path | None = None):
    cache_npz = None
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_npz = cache_dir / "gse131907_thesis_panel.npz"
        if cache_npz.exists():
            print(f"  loading cached panel {cache_npz}", flush=True)
            z = np.load(cache_npz, allow_pickle=True)
            cell_ids = z["cell_ids"].tolist()
            n_umi = z["n_umi"]
            n_streamed = int(z["n_streamed"])
            found = {k: z[k] for k in z.files if k not in {"cell_ids", "n_umi", "n_streamed"}}
            print(f"  cache genes={n_streamed} cells={len(cell_ids)} stored={len(found)}", flush=True)
            return cell_ids, found, n_umi, n_streamed
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
    if cache_npz is not None:
        np.savez_compressed(
            cache_npz,
            cell_ids=np.asarray(cell_ids, dtype=object),
            n_umi=n_umi,
            n_streamed=np.asarray(n_streamed),
            **found,
        )
        print(f"  wrote cache {cache_npz}", flush=True)
    return cell_ids, found, n_umi, n_streamed


def load_gse131907(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE131907 ====", flush=True)
    keep = set(locked["patient_id"].astype(str))
    ann_path = raw / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    mat_path = raw / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t")
    cell_ids, found, n_umi, n_genes = stream_131907(mat_path, wanted, cache_dir=raw / "GSE131907")
    cell_index = {c: i for i, c in enumerate(cell_ids)}
    key = "Index" if "Index" in ann.columns else ann.columns[0]
    ann = ann.copy()
    ann["cell_i"] = ann[key].map(cell_index)
    ann = ann[ann["cell_i"].notna()].copy()
    ann["cell_i"] = ann["cell_i"].astype(int)
    sample_col = "Sample" if "Sample" in ann.columns else "sample"
    origin_col = "Sample_Origin" if "Sample_Origin" in ann.columns else None
    type_col = "Cell_type" if "Cell_type" in ann.columns else None
    sub_col = "Cell_subtype" if "Cell_subtype" in ann.columns else None

    mal_mask = np.zeros(len(cell_ids), dtype=bool)
    tnk_mask = np.zeros(len(cell_ids), dtype=bool)
    sample_of = np.array([""] * len(cell_ids), dtype=object)
    idx_arr = ann["cell_i"].to_numpy()
    sample_arr = ann[sample_col].astype(str).to_numpy()
    subtype_arr = ann[sub_col].astype(str).to_numpy() if sub_col else np.array([""] * len(ann), dtype=object)
    ctype_arr = ann[type_col].astype(str).to_numpy() if type_col else np.array([""] * len(ann), dtype=object)
    origin_arr = ann[origin_col].astype(str).to_numpy() if origin_col else np.array([""] * len(ann), dtype=object)
    sample_of[idx_arr] = sample_arr
    mal_ok = np.isin(subtype_arr, list(MALIG_SUBTYPES)) & (
        np.isin(origin_arr, list(TUMOR_ORIGINS)) if origin_col else True
    )
    tnk_ok = np.isin(ctype_arr, ["T lymphocytes", "NK cells"])
    mal_mask[idx_arr[mal_ok]] = True
    tnk_mask[idx_arr[tnk_ok]] = True

    cldn4_all = None
    units = []
    for sample in sorted(keep):
        idx = np.flatnonzero(sample_of == sample)
        if idx.size == 0:
            print(f"  MISSING sample {sample}", flush=True)
            continue
        mal = mal_mask[idx]
        tnk = tnk_mask[idx]
        extracted = {g: found[g][idx] for g in found}
        libs = n_umi[idx]
        log_u, pos_u = to_log_pos(extracted, libs)
        cldn4 = log_u.get("CLDN4", np.zeros(int(idx.size), dtype=np.float32))
        if cldn4_all is None:
            cldn4_all = True
        units.append(
            {
                "cohort": "GSE131907",
                "patient_id": sample,
                "n_cells": int(idx.size),
                "n_mal": int(mal.sum()),
                "n_tnk": int(tnk.sum()),
                "mal": mal,
                "tnk": tnk,
                "cldn4": cldn4,
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


def _read_10x_features(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = next(n for n in members if f"_{sample}_features" in n or f"_{sample}_genes" in n)
    genes = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as f:
        for line in f:
            p = line.decode().strip().split("\t")
            genes.append((p[1] if len(p) > 1 else p[0]).split(".")[0].upper())
    return genes


def _read_10x_nbarcodes(tf: tarfile.TarFile, members: dict, sample: str) -> int:
    name = next(n for n in members if f"_{sample}_barcodes" in n)
    n = 0
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as f:
        for _ in f:
            n += 1
    return n


def stream_sample_mtx(tf, members, sample: str, genes: set[str]):
    symbols = _read_10x_features(tf, members, sample)
    n = _read_10x_nbarcodes(tf, members, sample)
    want_rows = {i for i, s in enumerate(symbols) if s in genes}
    keep: dict[str, np.ndarray] = {}
    for i in want_rows:
        keep.setdefault(symbols[i], np.zeros(n, dtype=np.float32))
    lib = np.zeros(n, dtype=np.float64)
    mtx_name = next(n_ for n_ in members if f"_{sample}_matrix.mtx" in n_)
    print(f"  {sample} cells={n} want_rows={len(want_rows)}", flush=True)
    with gzip.GzipFile(fileobj=tf.extractfile(members[mtx_name])) as f:
        for line in f:
            if line.decode().startswith("%"):
                continue
            break
        n_lines = 0
        for line in f:
            parts = line.decode().split()
            if len(parts) < 3:
                continue
            r = int(parts[0]) - 1
            c = int(parts[1]) - 1
            v = float(parts[2])
            lib[c] += v
            if r in want_rows:
                keep[symbols[r]][c] += v
            n_lines += 1
            if n_lines % 5_000_000 == 0:
                print(f"    mtx entries={n_lines}", flush=True)
    return keep, lib


def load_gse189357(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE189357 ====", flush=True)
    keep = set(locked["patient_id"].astype(str))
    tar_path = raw / "GSE189357" / "GSE189357_RAW.tar"
    units = []
    gene_union: set[str] = set()
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers()}
        samples = [f"TD{i}" for i in range(1, 10)]
        for sample in samples:
            if sample not in keep:
                print(f"  skip {sample} (not locked/eligible)", flush=True)
                continue
            extracted, libs = stream_sample_mtx(tf, members, sample, wanted)
            gene_union |= set(extracted)
            n = int(libs.size)
            log_cp, pos = to_log_pos(extracted, libs)
            mal, tnk = _marker_mal_tnk(extracted, n)
            cldn4 = log_cp.get("CLDN4", np.zeros(n, dtype=np.float32))
            units.append(
                {
                    "cohort": "GSE189357",
                    "patient_id": sample,
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
            print(f"  unit {sample}: n={n} mal={int(mal.sum())} tnk={int(tnk.sum())}", flush=True)
    if not units:
        raise SystemExit("GSE189357: no locked patient matrices read")
    return units, gene_union


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
            "n_131907": 0,
            "n_189357": 0,
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
        "n_131907": int((sub["cohort"] == "GSE131907").sum()),
        "n_189357": int((sub["cohort"] == "GSE189357").sum()),
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
                "n_gse131907": m["n_131907"],
                "n_gse189357": m["n_189357"],
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
                    "n_gse131907": m["n_131907"],
                    "n_gse189357": m["n_189357"],
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
        ax.text(
            r - 0.02 if r < 0 else r + 0.02,
            i,
            f"n={n}  ρ={r:.3f}",
            va="center",
            ha="right" if r < 0 else "left",
            fontsize=8,
        )
    ax.set_xlabel("Spearman ρ (malignant CLDN4 %pos vs T/NK frac)")
    ax.set_title("PR #459 given pair — not re-audited")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_honest_n(cov: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    labels = ["given combo", "GSE131907 locked", "GSE189357 locked", "median LR", "tertile LR", "Q4 vs Q1 LR"]
    vals = [
        GIVEN["n"],
        int((cov["cohort"] == "GSE131907").sum()),
        int((cov["cohort"] == "GSE189357").sum()),
        int(cov["median_ok"].sum()),
        int(cov["tertile_ok"].sum()),
        int(cov["q4q1_ok"].sum()),
    ]
    ax.barh(labels[::-1], vals[::-1], color="#4C78A8")
    for i, v in enumerate(vals[::-1]):
        ax.text(v + 0.3, i, str(v), va="center", fontsize=8)
    ax.set_xlabel("honest n (patient or locked sample)")
    ax.set_title("Given T/NK n is not the ligand n")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_n_per_unit(cov: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.2), sharey=False)
    for ax, cohort in zip(axes, ["GSE131907", "GSE189357"]):
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
    ax.set_xlabel("mean ΔP ± SE (patient/sample unit)")
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
    for ax, cohort, color in zip(axes, ["GSE131907", "GSE189357"], ["#4C78A8", "#F58518"]):
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
    n131 = int((cov["cohort"] == "GSE131907").sum())
    n189 = int((cov["cohort"] == "GSE189357").sum())
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

    md = f"""# FINDING — CLDN4-only thesis-aligned ligand test on GSE131907 + GSE189357

ADDITIVE. **CLDN4 only. No dual-high.** Patient is the unit
(GSE131907 locked sample from PR #459; GSE189357 patient TD1–TD9).
Do **not** add GSE148071. Do **not** re-audit the T/NK Spearman.

Thesis (already correct; this folder tests it on the given cut):

- CLDN4-high → more barrier/inhibitory outgoing to T/NK (F11R, NECTIN2–TIGIT, CDH1, LGALS9)
- CLDN4-low / KD-like → more IFN / T-recruit / MHC-I outgoing (CXCL9/10–CXCR3, CCL5, HLA–CD8)

The pair that **differs** is taken as given from PR #459:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE131907+GSE189357 | %pos | 30 | −0.542 (0.00291, 0%) | −0.619 (0.042; 9 vs 7) |

Singles (given, same PR): GSE131907 n=21 ρ=−0.522; GSE189357 n=9 ρ=−0.600.

This folder scores **both pre-specified families** as CellChat-style outgoing
Hill probabilities (Jin et al. 2021; 10% truncated mean, Kh={KH}, expr_prop≥{EXPR_PROP})
from CLDN4-high vs CLDN4-low malignant cells to same-unit T/NK.
CellChat R was not run. Pairs were not discovered.

## Honest n

Given combo n=30 is **not** the ligand n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (median: ≥{2 * MIN_CELLS_ARM}
malignant and ≥{MIN_CELLS_ARM}/arm plus ≥{MIN_CELLS_ARM} T/NK; Q4 vs Q1:
n_mal≥{MIN_MAL_Q4} and ≥{MIN_CELLS_ARM}/arm).

| gate | n | note |
|---|---:|---|
| Given combo (do not re-audit) | 30 | 21 samples + 9 patients |
| Locked GSE131907 samples | {n131} | author Malignant cells / tS*; T lymphocytes + NK cells |
| Locked GSE189357 patients | {n189} | marker-malignant; TD1–TD9 |
| Median-split LR | **{n_med}** | primary ligand n |
| Tertile extra | {n_ter} | exclusive arms |
| Q4 vs Q1 extra | {n_q} | n_mal≥40 |

GSE131907 malignant = author `Cell_subtype` in {{Malignant cells, tS1, tS2, tS3}}.
GSE189357 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0).
TACSTD2 is never a gate.

Median can drop units whose CLDN4 is almost all zero (high arm <10).
Ligand n is pair-specific after the expr_prop≥{EXPR_PROP} detection floor.

## Verdict

Family tables below. Barrier/inhibitory is tested as high>low; IFN/T-recruit/MHC-I
as low>high. Thin n and opposite calls are reported, not dropped.

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

- The n=30 ρ=−0.542 / Q4 r=−0.619 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071 is not added.
- Cell-pooled permutations are not the test. Patient is the unit.
- This is not a CellChat discovery screen. Only the two pre-specified families.

## Reproduce

```bash
python3 methods/pair_131907_189357_lr_thesis_cldn4/scripts/download.py
python3 methods/pair_131907_189357_lr_thesis_cldn4/scripts/analyze.py
```

Hill constants: trim={TRIM}, Kh={KH}, expr_prop={EXPR_PROP}.
"""
    return md


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("/tmp/pair_131907_189357_raw"))
    args = ap.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    locked = load_locked()
    print(
        f"locked GSE131907 n={len(locked['GSE131907'])} GSE189357 n={len(locked['GSE189357'])}",
        flush=True,
    )
    wanted = thesis_genes()
    names = {r["interaction_name"] for fam in FAMILIES.values() for r in fam["pairs"]}
    lr_all = load_lr(DB, names, matrix_genes=None)
    print(f"pre-specified pairs in CellChatDB: {len(lr_all)} / {len(names)}", flush=True)
    missing = names - set(lr_all["interaction_name"])
    if missing:
        print("WARNING missing interactions:", sorted(missing), flush=True)

    u131, genes131, n_streamed = load_gse131907(args.raw, wanted, locked["GSE131907"])
    u189, genes189 = load_gse189357(args.raw, wanted, locked["GSE189357"])
    units = u131 + u189
    gene_union = genes131 | genes189
    lr = load_lr(DB, names, matrix_genes=None)
    print(
        f"units loaded: {len(units)}  genes stored 131907={len(genes131)} 189357={len(genes189)} union={len(gene_union)}",
        flush=True,
    )

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
        "n_locked": {
            "GSE131907": int((cov.cohort == "GSE131907").sum()),
            "GSE189357": int((cov.cohort == "GSE189357").sum()),
        },
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
    print(
        "family tables:",
        RESULTS / "family_barrier_inhibitory.tsv",
        RESULTS / "family_ifn_recruit_mhci.tsv",
        flush=True,
    )
    print(fam_bar.to_string(index=False), flush=True)
    print(fam_ifn.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
