#!/usr/bin/env python3
"""ADDITIVE CLDN4-only thesis-aligned ligand test: GSE189357 + GSE205335.

The pair %pos Spearman is taken as given from PR #459 and is not re-audited:
  n=31, ρ=−0.478, Q4 vs Q1 r=−0.750.

No dual-high. Patient is the unit (GSE189357 TD1–TD9; GSE205335 locked patients).
Both pre-specified families are scored. Done when both family tables exist.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import tarfile
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lib import (  # noqa: E402
    EPI_MARKERS,
    EXPR_PROP,
    FAMILIES,
    KH,
    MIN_CELLS_ARM,
    MIN_MAL_Q4,
    TNK_MARKERS,
    TRIM,
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
    "combo": "GSE189357+GSE205335",
    "score": "pct",
    "n": 31,
    "rho": -0.4783764213621875,
    "p": 0.00920396129366116,
    "I2": 0.0,
    "q4_r": -0.75,
    "q4_p": 0.010411810411810413,
    "n_q1": 8,
    "n_q4": 8,
    "source": "PR459",
    "re_audited": False,
    "member_rhos": {"GSE189357": -0.600, "GSE205335": -0.435},
    "member_ns": {"GSE189357": 9, "GSE205335": 22},
}
GIVEN_SINGLES = [
    {"cohort": "GSE189357", "n": 9, "rho": -0.600, "p": 0.088, "score": "pct"},
    {"cohort": "GSE205335", "n": 22, "rho": -0.435, "p": 0.043, "score": "pct"},
    {"cohort": "combo (given)", "n": 31, "rho": -0.478, "p": 0.009, "score": "pct"},
]
COHORT_A = "GSE189357"
COHORT_B = "GSE205335"


def to_log_pos(extracted: dict[str, np.ndarray], library: np.ndarray):
    lib = np.maximum(library, 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}
    return log_cp, pos


def load_locked() -> dict[str, pd.DataFrame]:
    out = {}
    d = pd.read_csv(LOCKED / "GSE189357_marker_units.tsv", sep="\t")
    el = d[d["eligible"].astype(str).str.lower() == "true"].copy()
    el["unit"] = "GSE189357"
    el["patient_id"] = el["patient"].astype(str)
    el["cldn4_given"] = el["mal_CLDN4_pct"].astype(float)
    el["frac_tnk_given"] = el["frac_tnk"].astype(float)
    out["GSE189357"] = el

    s = pd.read_csv(LOCKED / "GSE205335_patients.tsv", sep="\t")
    el2 = s[(s["n_malignant"] > 0) & (s["n_tnk"] > 0)].copy()
    el2["unit"] = "GSE205335"
    el2["patient_id"] = el2["patient"].astype(str)
    el2["cldn4_given"] = el2["mal_CLDN4_pct_pos"].astype(float)
    el2["frac_tnk_given"] = el2["frac_tnk"].astype(float)
    out["GSE205335"] = el2
    return out


def _marker_mal_tnk(extracted: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    def col(g: str) -> np.ndarray:
        if g not in extracted:
            return np.zeros(n, dtype=float)
        return extracted[g]

    epi = np.zeros(n, dtype=bool)
    for g in EPI_MARKERS:
        epi |= col(g) > 0
    mal = epi & (col("PTPRC") == 0)
    tnk = np.zeros(n, dtype=bool)
    for g in TNK_MARKERS:
        tnk |= col(g) > 0
    tnk = tnk & (~mal)
    return mal, tnk


def _member(members: dict[str, tarfile.TarInfo], sample: str, kind: str) -> str:
    for name in members:
        if f"_{sample}_{kind}" in Path(name).name:
            return name
    raise KeyError(f"no {kind} for {sample}")


def _read_10x_features(tf: tarfile.TarFile, members: dict, sample: str) -> list[str]:
    name = _member(members, sample, "features")
    genes = []
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as handle:
        for line in handle:
            parts = line.decode().strip().split("\t")
            genes.append((parts[1] if len(parts) > 1 else parts[0]).upper())
    return genes


def _read_10x_nbarcodes(tf: tarfile.TarFile, members: dict, sample: str) -> int:
    name = _member(members, sample, "barcodes")
    n = 0
    with gzip.GzipFile(fileobj=tf.extractfile(members[name])) as handle:
        for _ in handle:
            n += 1
    return n


def _stream_10x_wanted(tf, members, sample, genes, wanted) -> tuple[dict[str, np.ndarray], np.ndarray]:
    n_cells = _read_10x_nbarcodes(tf, members, sample)
    keep_idx = {i for i, g in enumerate(genes) if g in wanted}
    extracted = {genes[i]: np.zeros(n_cells, dtype=np.float32) for i in keep_idx}
    library = np.zeros(n_cells, dtype=np.float64)
    mtx_name = _member(members, sample, "matrix.mtx")
    with gzip.GzipFile(fileobj=tf.extractfile(members[mtx_name])) as handle:
        for line in handle:
            if line.startswith(b"%"):
                continue
            break
        for line in handle:
            a, b, v = line.decode().split()
            gi = int(a) - 1
            ci = int(b) - 1
            val = float(v)
            library[ci] += val
            if gi in keep_idx:
                extracted[genes[gi]][ci] = val
    return extracted, library


def load_gse189357(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE189357 ====", flush=True)
    keep = set(locked["patient_id"].astype(str))
    tar_path = raw / "GSE189357" / "GSE189357_RAW.tar"
    units = []
    gene_union: set[str] = set()
    with tarfile.open(tar_path) as tf:
        members = {m.name: m for m in tf.getmembers() if m.isfile()}
        for sample in [f"TD{i}" for i in range(1, 10)]:
            if sample not in keep:
                print(f"  skip {sample} (not in locked extract)", flush=True)
                continue
            genes = _read_10x_features(tf, members, sample)
            gene_union |= set(genes)
            extracted, libs = _stream_10x_wanted(tf, members, sample, genes, wanted)
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
            print(
                f"  {sample}: n={n} mal={int(mal.sum())} tnk={int(tnk.sum())} stored={len(extracted)}",
                flush=True,
            )
    if not units:
        raise SystemExit("GSE189357: no locked tumor matrices read")
    return units, gene_union


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def load_rds_genes(path: Path, wanted: set[str]):
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            print(f"decompress {path.name}", flush=True)
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        print("read GSE205335 RDS", flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel()
    print(f"GSE205335 extracted {len(extracted)} / {len(wanted)} genes", flush=True)
    return extracted, library_umi, barcodes, set(genes.tolist())


def load_gse205335(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE205335 ====", flush=True)
    keep = set(locked["patient_id"].astype(str))
    ident_path = raw / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    soft_path = raw / "GSE205335" / "GSE205335_family.soft.gz"
    mat_path = raw / "GSE205335" / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    identities = pd.read_csv(ident_path, sep="\t")
    if identities["barcode"].duplicated().any():
        raise ValueError("GSE205335 identity barcodes are not unique")
    metadata = parse_geo_soft(soft_path)
    extracted, library_umi, barcodes, genes = load_rds_genes(mat_path, wanted)
    if "CLDN4" not in extracted:
        raise SystemExit("CLDN4 missing from GSE205335 UMI")
    indexed = identities.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise ValueError(
            f"GSE205335 matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra"
        )
    cells = indexed.loc[barcodes].reset_index()
    cells = cells.merge(
        metadata[["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"]],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        raise ValueError("GSE205335 identity samples did not match GEO metadata")
    log_cp, pos = to_log_pos(extracted, library_umi)
    mal_all = cells["lineage.sub"].eq("Malignant cells").to_numpy()
    tnk_all = cells["lineage.total"].eq("T/NK cells").to_numpy()
    patient = cells["patient"].astype(str).to_numpy()
    cldn4 = log_cp["CLDN4"]
    units = []
    for pid in sorted(keep):
        idx = np.flatnonzero(patient == pid)
        if idx.size == 0:
            print(f"  MISSING patient {pid}", flush=True)
            continue
        mal = mal_all[idx]
        tnk = tnk_all[idx]
        extracted_u = {g: extracted[g][idx] for g in extracted}
        libs = library_umi[idx]
        log_u, pos_u = to_log_pos(extracted_u, libs)
        units.append(
            {
                "cohort": "GSE205335",
                "patient_id": pid,
                "n_cells": int(idx.size),
                "n_mal": int(mal.sum()),
                "n_tnk": int(tnk.sum()),
                "mal": mal,
                "tnk": tnk,
                "cldn4": cldn4[idx],
                "log_cp": log_u,
                "pos": pos_u,
                "genes": set(extracted),
            }
        )
        print(
            f"  {pid}: n={int(idx.size)} mal={int(mal.sum())} tnk={int(tnk.sum())}",
            flush=True,
        )
    return units, genes


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
            "n_189357": 0,
            "n_205335": 0,
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
        "n_189357": int((sub["cohort"] == COHORT_A).sum()),
        "n_205335": int((sub["cohort"] == COHORT_B).sum()),
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
                "n_gse189357": m["n_189357"],
                "n_gse205335": m["n_205335"],
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
                    "n_gse189357": m["n_189357"],
                    "n_gse205335": m["n_205335"],
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
        "| pair | axis | expect | n | n_189357/n_205335 | mean ΔP | p_W | observed | agrees |",
        "|---|---|---|---:|---|---:|---|---|---|",
    ]
    for rec in df.itertuples(index=False):
        name = rec.interaction_name.replace("FAMILY:", "**FAMILY** ")
        lines.append(
            f"| {name} | {rec.axis} | {rec.thesis_expect} | {rec.n_units} | "
            f"{int(rec.n_gse189357)}/{int(rec.n_gse205335)} | "
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
    labels = [
        "given combo",
        "GSE189357 locked",
        "GSE205335 locked",
        "median LR",
        "tertile LR",
        "Q4 vs Q1 LR",
    ]
    vals = [
        GIVEN["n"],
        int((cov["cohort"] == COHORT_A).sum()),
        int((cov["cohort"] == COHORT_B).sum()),
        int(cov["median_ok"].sum()),
        int(cov["tertile_ok"].sum()),
        int(cov["q4q1_ok"].sum()),
    ]
    ax.barh(labels[::-1], vals[::-1], color="#4C78A8")
    for i, v in enumerate(vals[::-1]):
        ax.text(v + 0.3, i, str(v), va="center", fontsize=8)
    ax.set_xlabel("honest n (patient)")
    ax.set_title("Given T/NK n is not the ligand n")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_n_per_unit(cov: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.2), sharey=False)
    for ax, cohort in zip(axes, [COHORT_A, COHORT_B]):
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
    ax.set_xlabel("mean ΔP ± SE (patient unit)")
    ax.set_title(f"{spec['label']} — {mode}")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig_family_patient(long_df: pd.DataFrame, family_key: str, mode: str, path: Path) -> None:
    spec = FAMILIES[family_key]
    names = {r["interaction_name"] for r in spec["pairs"]}
    sub = long_df[
        (long_df["interaction_name"].isin(names))
        & (long_df["split"] == mode)
        & (long_df["detected"])
    ]
    if sub.empty:
        return
    fam = sub.groupby(["cohort", "patient_id"], as_index=False).agg(delta=("delta", "mean"))
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0), sharex=True)
    for ax, cohort, color in zip(axes, [COHORT_A, COHORT_B], ["#4C78A8", "#F58518"]):
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
    n189 = int((cov["cohort"] == COHORT_A).sum())
    n205 = int((cov["cohort"] == COHORT_B).sum())
    bar_agg = fam_bar[fam_bar["row"] == "FAMILY_AGGREGATE"]
    ifn_agg = fam_ifn[fam_ifn["row"] == "FAMILY_AGGREGATE"]
    bar_txt = "NA"
    ifn_txt = "NA"
    if not bar_agg.empty:
        r = bar_agg.iloc[0]
        bar_txt = (
            f"n={int(r.n_units)} mean ΔP={fmt_num(r.mean_delta)} "
            f"p_W={fmt_p(r.wilcoxon_p)} agrees={r.agrees_thesis}"
        )
    if not ifn_agg.empty:
        r = ifn_agg.iloc[0]
        ifn_txt = (
            f"n={int(r.n_units)} mean ΔP={fmt_num(r.mean_delta)} "
            f"p_W={fmt_p(r.wilcoxon_p)} agrees={r.agrees_thesis}"
        )

    md = f"""# FINDING — CLDN4-only thesis-aligned ligand test on GSE189357 + GSE205335

ADDITIVE. **CLDN4 only. No dual-high.** Patient is the unit
(GSE189357 TD1–TD9 marker-malignant; GSE205335 locked author-malignant patients).
Do **not** re-audit the T/NK Spearman.

Thesis (already correct; this folder tests it on the given cut):

- CLDN4-high → more barrier/inhibitory outgoing to T/NK (F11R, NECTIN2–TIGIT, CDH1, LGALS9)
- CLDN4-low / KD-like → more IFN / T-recruit / MHC-I outgoing (CXCL9/10–CXCR3, CCL5, HLA–CD8)

The pair that **differs** is taken as given from PR #459:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE189357+GSE205335 | %pos | 31 | −0.478 (0.009, 0%) | −0.750 (0.010; 8 vs 8) |

Singles (given, same PR): GSE189357 n=9 ρ=−0.600; GSE205335 n=22 ρ=−0.435.

This folder scores **both pre-specified families** as CellChat-style outgoing
Hill probabilities (Jin et al. 2021; 10% truncated mean, Kh={KH}, expr_prop≥{EXPR_PROP})
from CLDN4-high vs CLDN4-low malignant cells to same-patient T/NK.
CellChat R was not run. Pairs were not discovered.

## Honest n

Given combo n=31 is **not** the ligand n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (median: ≥{2 * MIN_CELLS_ARM}
malignant and ≥{MIN_CELLS_ARM}/arm plus ≥{MIN_CELLS_ARM} T/NK; Q4 vs Q1:
n_mal≥{MIN_MAL_Q4} and ≥{MIN_CELLS_ARM}/arm).

| gate | n | note |
|---|---:|---|
| Given combo (do not re-audit) | 31 | 9 + 22 patients |
| Locked GSE189357 patients | {n189} | marker-malignant; TD1–TD9 |
| Locked GSE205335 patients | {n205} | author Malignant cells; T/NK cells |
| Median-split LR | **{n_med}** | primary ligand n |
| Tertile extra | {n_ter} | exclusive arms |
| Q4 vs Q1 extra | {n_q} | n_mal≥40 |

GSE189357 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0).
GSE205335 malignant = author `lineage.sub` == `Malignant cells`.
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

- The n=31 ρ=−0.478 / Q4 r=−0.750 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- Cell-pooled permutations are not the test. Patient is the unit.
- This is not a CellChat discovery screen. Only the two pre-specified families.
- Malignant definitions differ (marker vs author) and are stated on every row.

## Reproduce

```bash
python3 methods/pair_189357_205335_lr_thesis_cldn4/scripts/download.py
python3 methods/pair_189357_205335_lr_thesis_cldn4/scripts/analyze.py
```

Hill constants: trim={TRIM}, Kh={KH}, expr_prop={EXPR_PROP}.
"""
    return md


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("/tmp/pair_189357_205335_raw"))
    args = ap.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    locked = load_locked()
    print(
        f"locked GSE189357 n={len(locked['GSE189357'])} GSE205335 n={len(locked['GSE205335'])}",
        flush=True,
    )
    wanted = thesis_genes()
    names = {r["interaction_name"] for fam in FAMILIES.values() for r in fam["pairs"]}
    lr = load_lr(DB, names, matrix_genes=None)
    print(f"pre-specified pairs in CellChatDB: {len(lr)} / {len(names)}", flush=True)
    missing = names - set(lr["interaction_name"])
    if missing:
        print("WARNING missing interactions:", sorted(missing), flush=True)

    u189, genes189 = load_gse189357(args.raw, wanted, locked["GSE189357"])
    u205, genes205 = load_gse205335(args.raw, wanted, locked["GSE205335"])
    units = u189 + u205
    print(
        f"units loaded: {len(units)}  genes stored 189357={len(genes189)} 205335={len(genes205)}",
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
            "GSE189357": int((cov.cohort == COHORT_A).sum()),
            "GSE205335": int((cov.cohort == COHORT_B).sum()),
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
        "re_audited_tnk_rho": False,
        "dual_high": False,
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
