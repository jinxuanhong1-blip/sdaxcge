#!/usr/bin/env python3
"""ADDITIVE CLDN4-only high-end CellChat + LIANA on GSE123902 + GSE131907.

The pair %pos Spearman is taken as given from PR #459 and is not re-audited:
  n=34, ρ=−0.575, Q4 vs Q1 r=−0.700.

No dual-high. Do not add GSE148071. Do not pile all seven.
Patient/donor is the unit (GSE131907 sample).
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
    EXPR_PROP,
    EXTRA,
    KH,
    LINEAGE_MARKERS,
    MIN_CELLS_ARM,
    MIN_CELLS_COMP,
    MIN_MAL_Q4,
    NORMAL_LUNG,
    TRIM,
    assign_patient_quartiles,
    dersimonian_laird,
    fmt_num,
    fmt_p,
    highend_split,
    ligand_class,
    load_lr,
    parse_symbols,
    score_outgoing,
    unit_delta_stats,
)

LOCKED = ROOT / "data"
DB = ROOT / "db"
UNITS = ["GSE123902", "GSE131907"]

# PR #459 given pair — not re-audited
GIVEN = {
    "combo": "GSE123902+GSE131907",
    "score": "pct",
    "n": 34,
    "rho": -0.575,
    "p": 0.001,
    "I2": 0.0,
    "q4_r": -0.700,
    "q4_p": 0.012,
    "n_q1": 10,
    "n_q4": 8,
    "source": "PR459",
    "re_audited": False,
}
GIVEN_SINGLES = [
    {"cohort": "GSE123902", "n": 13, "rho": -0.659, "p": 0.014, "score": "pct"},
    {"cohort": "GSE131907", "n": 21, "rho": -0.522, "p": 0.015, "score": "pct"},
    {"cohort": "combo (given)", "n": 34, "rho": -0.575, "p": 0.001, "score": "pct"},
]

EPI_MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19"]
TNK_MARKERS = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
TUMOR_ORIGINS_131907 = ("tLung", "tL/B", "mLN", "PE", "mBrain")
TNK_TYPES_131907 = {"T lymphocytes", "NK cells"}
SCORE_DROP = {
    "prob",
    "liana_score",
    "detected",
    "ligand_mean",
    "receptor_mean",
    "ligand_prop",
    "receptor_prop",
}


def to_log_pos(extracted: dict[str, np.ndarray], library: np.ndarray):
    lib = np.maximum(library, 1.0)
    log_cp = {g: np.log1p(extracted[g] / lib * 1e4).astype(np.float32) for g in extracted}
    pos = {g: (extracted[g] > 0).astype(np.float32) for g in extracted}
    return log_cp, pos


def wanted_genes() -> set[str]:
    inter = pd.read_csv(DB / "interaction_cellchatdb_v2_protein.csv")
    genes: set[str] = set(EXTRA) | set(NORMAL_LUNG) | set(EPI_MARKERS) | set(TNK_MARKERS) | {"PTPRC"}
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    for rec in inter.itertuples(index=False):
        genes.update(parse_symbols(getattr(rec, "ligand.symbol", None) or rec.ligand))
        genes.update(parse_symbols(getattr(rec, "receptor.symbol", None) or rec.receptor))
    return genes


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

    h = pd.read_csv(LOCKED / "GSE131907_samples.tsv", sep="\t")
    tumor = h[h["origin"].isin(TUMOR_ORIGINS_131907)].copy()
    mal = tumor[tumor["n_malignant"] >= 20].copy()
    mal["unit"] = "GSE131907"
    mal["patient_id"] = mal["sample"].astype(str)
    mal["cldn4_given"] = mal["mal_CLDN4_pct"].astype(float)
    mal["frac_tnk_given"] = mal["frac_tnk"].astype(float)
    out["GSE131907"] = mal
    return out


def _marker_mal_tnk(extracted: dict[str, np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    """Locked PR #459 gates: count>0 EPCAM|KRT* and PTPRC==0; T/NK markers and not mal."""

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


def _score_one_123902_csv(handle, wanted: set[str]) -> tuple[dict[str, np.ndarray], np.ndarray, set[str]]:
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
        chunks: list[tuple[str, dict[str, np.ndarray], np.ndarray]] = []
        gene_union: set[str] = set()
        seen_wanted: set[str] = set()
        for donor, tissue, m, name in members:
            if donor in chosen:
                continue
            chosen[donor] = tissue
            raw_f = gzip.GzipFile(fileobj=tf.extractfile(m))
            extracted, libs, genes = _score_one_123902_csv(raw_f, wanted)
            gene_union |= genes
            seen_wanted |= set(extracted)
            n = int(libs.size)
            chunks.append((donor, extracted, libs))
            print(f"  {name}: cells={n} stored={len(extracted)} donor={donor} {tissue}", flush=True)
    if not chunks:
        raise SystemExit("GSE123902: no locked tumor matrices read")
    libs_all = []
    patients = []
    extracted_all: dict[str, list[np.ndarray]] = {g: [] for g in seen_wanted}
    for donor, extracted, libs in chunks:
        n = int(libs.size)
        libs_all.append(libs)
        patients.append(np.full(n, donor, dtype=object))
        for g in seen_wanted:
            extracted_all[g].append(extracted.get(g, np.zeros(n, dtype=np.float32)))
    library = np.concatenate(libs_all)
    extracted = {g: np.concatenate(vs) for g, vs in extracted_all.items()}
    patient = np.concatenate(patients)
    n = int(patient.size)
    mal, tnk = _marker_mal_tnk(extracted, n)
    log_cp, pos = to_log_pos(extracted, library)
    print(
        f"  donors={len(chosen)} cells={n} mal={int(mal.sum())} tnk={int(tnk.sum())}",
        flush=True,
    )
    return {
        "patient": patient,
        "mal": mal,
        "tnk": tnk,
        "log_cp": log_cp,
        "pos": pos,
        "genes": gene_union,
        "keep": keep,
        "note": "Laughney 2020; marker-malignant (EPCAM|KRT8/18/19>0 & PTPRC==0); donor-level; normals dropped",
    }


def stream_131907(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = np.asarray(header[1:], dtype=str)
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
                print(f"  131907 stream genes={n_streamed} stored={len(found)}", flush=True)
    print(f"131907 stream done genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    return cell_ids, found, n_umi, n_streamed


def load_gse131907(raw: Path, wanted: set[str], locked: pd.DataFrame):
    print("==== GSE131907 ====", flush=True)
    keep = set(locked["patient_id"].astype(str))
    ann = pd.read_csv(raw / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    cell_ids, extracted, n_umi, n_genes = stream_131907(
        raw / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", wanted
    )
    if "CLDN4" not in extracted:
        raise SystemExit("CLDN4 not in GSE131907 UMI matrix")
    bc_to_i = {b: i for i, b in enumerate(cell_ids)}
    keep_ann = ann["Index"].map(bc_to_i)
    ann = ann.loc[keep_ann.notna()].copy()
    ann["idx"] = keep_ann.loc[ann.index].astype(int)
    print(f"  annotation matched to matrix: {len(ann)} / {len(cell_ids)}", flush=True)
    n = int(cell_ids.size)
    mal = np.zeros(n, dtype=bool)
    tnk = np.zeros(n, dtype=bool)
    patient = np.full(n, "", dtype=object)
    sample = ann["Sample"].astype(str).to_numpy()
    idx = ann["idx"].to_numpy()
    mal_mask = ann["Cell_subtype"].astype(str).eq("Malignant cells").to_numpy()
    tnk_mask = ann["Cell_type"].astype(str).isin(TNK_TYPES_131907).to_numpy()
    mal[idx] = mal_mask
    tnk[idx] = tnk_mask
    patient[idx] = sample
    log_cp, pos = to_log_pos(extracted, n_umi)
    print(
        f"  mal={int(mal.sum())} tnk={int(tnk.sum())} locked_samples={len(keep)} genes_streamed={n_genes}",
        flush=True,
    )
    return {
        "patient": patient,
        "mal": mal,
        "tnk": tnk,
        "log_cp": log_cp,
        "pos": pos,
        "genes": set(extracted) | {"CLDN4"},
        "keep": keep,
        "note": "author Cell_subtype == Malignant cells; Cell_type in {T lymphocytes, NK cells}; sample-level; n_mal>=20",
    }


def score_bundle(unit: str, bundle: dict, lr: pd.DataFrame, modes: list[str]):
    patient = np.asarray(bundle["patient"], dtype=object)
    mal = np.asarray(bundle["mal"], dtype=bool)
    tnk = np.asarray(bundle["tnk"], dtype=bool)
    log_cp = bundle["log_cp"]
    pos = bundle["pos"]
    keep = set(bundle["keep"])
    cldn4 = log_cp.get("CLDN4")
    if cldn4 is None:
        raise SystemExit(f"{unit}: CLDN4 missing")
    n = patient.size
    idx_all = np.arange(n)
    coverage_rows = []
    pair_rows = []
    for pid in sorted(keep):
        mask = patient == pid
        mal_idx = idx_all[mask & mal]
        tnk_idx = idx_all[mask & tnk]
        rec = {
            "unit": unit,
            "patient_id": pid,
            "n_mal": int(mal_idx.size),
            "n_tnk": int(tnk_idx.size),
            "mean_cldn4_mal": float(cldn4[mal_idx].mean()) if mal_idx.size else np.nan,
            "eligible_allmal": int(mal_idx.size >= MIN_CELLS_COMP and tnk_idx.size >= MIN_CELLS_COMP),
        }
        if rec["eligible_allmal"]:
            scored = score_outgoing(lr, log_cp, pos, mal_idx, tnk_idx)
            for row in scored:
                row.update({"unit": unit, "patient_id": pid, "split": "all_mal", "arm": "all"})
                pair_rows.append(row)
        for mode in modes:
            split = highend_split(cldn4, mal_idx, mode)
            rec[f"eligible_{mode}"] = int(split is not None and tnk_idx.size >= MIN_CELLS_ARM)
            if split is None or tnk_idx.size < MIN_CELLS_ARM:
                rec[f"n_high_{mode}"] = 0
                rec[f"n_low_{mode}"] = 0
                continue
            hi, lo = split
            rec[f"n_high_{mode}"] = int(hi.size)
            rec[f"n_low_{mode}"] = int(lo.size)
            hi_rows = score_outgoing(lr, log_cp, pos, hi, tnk_idx)
            lo_rows = score_outgoing(lr, log_cp, pos, lo, tnk_idx)
            lo_map = {r["interaction_name"]: r for r in lo_rows}
            for hr in hi_rows:
                lr_low = lo_map.get(hr["interaction_name"])
                if lr_low is None:
                    continue
                pair_rows.append(
                    {
                        **{k: hr[k] for k in hr if k not in SCORE_DROP},
                        "unit": unit,
                        "patient_id": pid,
                        "split": mode,
                        "arm": "delta",
                        "prob_high": hr["prob"],
                        "prob_low": lr_low["prob"],
                        "delta": hr["prob"] - lr_low["prob"],
                        "liana_high": hr["liana_score"],
                        "liana_low": lr_low["liana_score"],
                        "delta_liana": hr["liana_score"] - lr_low["liana_score"],
                        "detected_high": hr["detected"],
                        "detected_low": lr_low["detected"],
                        "detected_either": bool(hr["detected"] or lr_low["detected"]),
                    }
                )
        coverage_rows.append(rec)
        print(
            f"  {unit} {pid} mal={rec['n_mal']} tnk={rec['n_tnk']} "
            f"q4={rec.get('eligible_q4q1', 0)} med={rec.get('eligible_median', 0)}",
            flush=True,
        )
    return pd.DataFrame(coverage_rows), pd.DataFrame(pair_rows)


def meta_patient_deltas(pairs: pd.DataFrame, split: str, delta_col: str = "delta") -> pd.DataFrame:
    sub = pairs[(pairs["split"] == split) & (pairs["arm"] == "delta") & (pairs["detected_either"])].copy()
    if sub.empty:
        return pd.DataFrame()
    rows = []
    for name, g in sub.groupby("interaction_name", sort=False):
        unit_stats = []
        for unit, ug in g.groupby("unit"):
            st = unit_delta_stats(ug[delta_col].to_numpy())
            st["unit"] = unit
            if st["n"] >= 2 and np.isfinite(st["se"]) and st["se"] > 0:
                unit_stats.append(st)
        st_all = unit_delta_stats(g[delta_col].to_numpy())
        first = g.iloc[0]
        if not unit_stats:
            rows.append(
                {
                    "interaction_name": name,
                    "pathway_name": first["pathway_name"],
                    "ligand": first["ligand"],
                    "receptor": first["receptor"],
                    "lr_class": first["lr_class"],
                    "split": split,
                    "score": "cellchat_hill" if delta_col == "delta" else "liana_mean",
                    "k_units": 0,
                    "n_patients": st_all["n"],
                    "mean_delta": st_all["mean"],
                    "median_delta": st_all["median"],
                    "se": st_all["se"],
                    "p_meta": np.nan,
                    "p_wilcoxon_patients": st_all["wilcoxon_p"],
                    "i2": np.nan,
                    "units": "",
                }
            )
            continue
        meta = dersimonian_laird(
            np.array([s["mean"] for s in unit_stats]),
            np.array([s["se"] ** 2 for s in unit_stats]),
        )
        rows.append(
            {
                "interaction_name": name,
                "pathway_name": first["pathway_name"],
                "ligand": first["ligand"],
                "receptor": first["receptor"],
                "lr_class": first["lr_class"],
                "split": split,
                "score": "cellchat_hill" if delta_col == "delta" else "liana_mean",
                "k_units": meta["k"],
                "n_patients": st_all["n"],
                "mean_delta": meta["mean"],
                "median_delta": st_all["median"],
                "se": meta["se"],
                "p_meta": meta["p"],
                "p_wilcoxon_patients": st_all["wilcoxon_p"],
                "i2": meta["i2"],
                "units": "+".join(sorted(s["unit"] for s in unit_stats)),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["p_meta", "n_patients"], ascending=[True, False], na_position="last")


def meta_q4q1_between(pairs: pd.DataFrame, locked: dict[str, pd.DataFrame], score_col: str = "prob") -> pd.DataFrame:
    allmal = pairs[(pairs["split"] == "all_mal") & (pairs["detected"])].copy()
    if allmal.empty:
        return pd.DataFrame()
    qmap = {}
    for unit, df in locked.items():
        try:
            qs = assign_patient_quartiles(df["cldn4_given"])
        except ValueError:
            continue
        for pid, q in zip(df["patient_id"], qs.astype(str)):
            qmap[(unit, str(pid))] = q
    allmal["quartile"] = [qmap.get((u, str(p))) for u, p in zip(allmal["unit"], allmal["patient_id"])]
    rows = []
    for name, g in allmal.groupby("interaction_name", sort=False):
        unit_effects = []
        n_q1 = n_q4 = 0
        for unit, ug in g.groupby("unit"):
            a = ug.loc[ug["quartile"] == "Q1", score_col]
            b = ug.loc[ug["quartile"] == "Q4", score_col]
            if len(a) < 2 or len(b) < 2:
                continue
            n_q1 += int(len(a))
            n_q4 += int(len(b))
            delta = float(b.median() - a.median())
            if (a.var(ddof=1) + b.var(ddof=1)) > 0:
                se = float(np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)))
            else:
                se = 0.05
            if se <= 0:
                se = 0.05
            unit_effects.append({"unit": unit, "delta": delta, "se": se})
        if not unit_effects:
            continue
        meta = dersimonian_laird(
            np.array([e["delta"] for e in unit_effects]),
            np.array([e["se"] ** 2 for e in unit_effects]),
        )
        first = g.iloc[0]
        rows.append(
            {
                "interaction_name": name,
                "pathway_name": first["pathway_name"],
                "ligand": first["ligand"],
                "receptor": first["receptor"],
                "lr_class": first.get("lr_class", ligand_class(str(first.get("ligand_genes", "")))),
                "split": "between_q4q1_given_pct",
                "score": "cellchat_hill" if score_col == "prob" else "liana_mean",
                "k_units": meta["k"],
                "n_q1": n_q1,
                "n_q4": n_q4,
                "n_compared": n_q1 + n_q4,
                "mean_delta": meta["mean"],
                "se": meta["se"],
                "p_meta": meta["p"],
                "i2": meta["i2"],
                "units": "+".join(sorted(e["unit"] for e in unit_effects)),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["p_meta", "n_compared"], ascending=[True, False], na_position="last")


def plot_n(coverage: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    units = [u for u in UNITS if u in set(coverage["unit"])]
    x = np.arange(len(units))
    locked_n = coverage.groupby("unit")["patient_id"].nunique().reindex(units).fillna(0)
    elig = (
        coverage.groupby("unit")["eligible_q4q1"].sum().reindex(units).fillna(0)
        if "eligible_q4q1" in coverage
        else pd.Series(0, index=units)
    )
    elig_m = (
        coverage.groupby("unit")["eligible_median"].sum().reindex(units).fillna(0)
        if "eligible_median" in coverage
        else pd.Series(0, index=units)
    )
    ax.bar(x - 0.25, locked_n, 0.24, label="locked n (given %pos)", color="#9e9e9e")
    ax.bar(x, elig_m, 0.24, label="median-split eligible", color="#6a8aaa")
    ax.bar(x + 0.25, elig, 0.24, label="Q4 vs Q1 eligible", color="#b2182b")
    ax.set_xticks(x, units, rotation=15, ha="right")
    ax.set_ylabel("patients / donors / samples")
    ax.set_title("Honest n: given combo vs LR-eligible units")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_n_per_patient(coverage: pd.DataFrame, path: Path) -> None:
    df = coverage.sort_values(["unit", "patient_id"]).copy()
    fig, ax = plt.subplots(figsize=(max(8.0, 0.28 * len(df)), 3.8))
    x = np.arange(len(df))
    ax.bar(x - 0.18, df["n_mal"], 0.36, label="malignant", color="#8C6D31")
    ax.bar(x + 0.18, df["n_tnk"], 0.36, label="T/NK", color="#4C72B0")
    ax.axhline(MIN_MAL_Q4, color="#b2182b", ls="--", lw=0.8, label=f"Q4 floor n_mal≥{MIN_MAL_Q4}")
    ax.axhline(MIN_CELLS_ARM, color="#2166ac", ls=":", lw=0.8, label=f"T/NK floor ≥{MIN_CELLS_ARM}")
    ax.set_xticks(x)
    labels = [f"{r.unit.replace('GSE', '')}:{r.patient_id}" for r in df.itertuples()]
    ax.set_xticklabels(labels, rotation=80, ha="right", fontsize=6)
    ax.set_ylabel("cells")
    ax.set_title("Per-unit cell floors (patient / donor / sample is the unit)")
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_forest(meta: pd.DataFrame, path: Path, title: str, n=12) -> None:
    show = meta.head(n).copy() if meta is not None and not meta.empty else pd.DataFrame()
    fig, ax = plt.subplots(figsize=(8.6, 0.42 * max(len(show), 1) + 1.6))
    if show.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No meta rows", ha="center")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    y = np.arange(len(show))
    colors = np.where(show["mean_delta"] >= 0, "#b2182b", "#2166ac")
    ax.axvline(0, color="#444", lw=0.8)
    ax.barh(y, show["mean_delta"], color=colors, alpha=0.85)
    if "se" in show:
        ax.errorbar(show["mean_delta"], y, xerr=1.96 * show["se"].fillna(0), fmt="none", ecolor="#333", lw=0.8)
    labels = []
    for r in show.itertuples():
        n_lab = int(getattr(r, "n_patients", getattr(r, "n_compared", 0)))
        labels.append(f"{r.interaction_name}  n={n_lab}  p={fmt_p(r.p_meta)}")
    ax.set_yticks(y, labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("meta Δ (CLDN4-high − low)")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_patient_deltas(pairs: pd.DataFrame, interaction: str, split: str, path: Path, delta_col="delta") -> None:
    sub = pairs[(pairs["split"] == split) & (pairs["arm"] == "delta") & (pairs["interaction_name"] == interaction)]
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    if sub.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, f"No patient deltas for {interaction}", ha="center")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return
    units = [u for u in UNITS if u in set(sub["unit"])]
    data = [sub.loc[sub["unit"] == u, delta_col].to_numpy() for u in units]
    bp = ax.boxplot(data, tick_labels=units, patch_artist=True, widths=0.55)
    for patch in bp["boxes"]:
        patch.set_facecolor("#f4a582")
        patch.set_alpha(0.6)
    rng = np.random.default_rng(0)
    for i, vals in enumerate(data, start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=14, zorder=3)
    ax.axhline(0, color="#444", lw=0.8)
    ylab = "patient ΔP" if delta_col == "delta" else "patient ΔS (LIANA)"
    ax.set_ylabel(ylab)
    ax.set_title(f"{interaction}  within-unit {split}  (patient is the unit)")
    ax.tick_params(axis="x", rotation=15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_ligand_table(tbl: pd.DataFrame, path: Path, title: str) -> None:
    show = tbl.head(18).copy() if tbl is not None and not tbl.empty else pd.DataFrame()
    fig, ax = plt.subplots(figsize=(8.8, 0.42 * max(len(show), 1) + 1.4))
    if show.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No differential pairs at the stated gates", ha="center")
        fig.savefig(path, dpi=150)
        fig.savefig(path.with_suffix(".pdf"))
        plt.close(fig)
        return
    y = np.arange(len(show))
    colors = np.where(show["mean_delta"] >= 0, "#b2182b", "#2166ac")
    ax.barh(y, show["mean_delta"], color=colors, alpha=0.85)
    ax.axvline(0, color="#444", lw=0.7)
    ax.set_yticks(
        y,
        [
            f"{r.interaction_name}  [{r.lr_class}]  n={int(r.n_patients)}  p={fmt_p(r.p_meta)}"
            for r in show.itertuples()
        ],
        fontsize=8,
    )
    ax.invert_yaxis()
    ax.set_xlabel("meta mean patient Δ")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_given_rho(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 2.8))
    y = np.arange(len(GIVEN_SINGLES))
    for i, r in enumerate(GIVEN_SINGLES):
        ax.plot(r["rho"], i, "o", color="#b2182b", ms=8)
        ax.text(
            r["rho"] - 0.02,
            i + 0.18,
            f"n={r['n']} p={fmt_p(r['p'])}",
            fontsize=7,
            ha="right",
        )
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_yticks(y, [r["cohort"] for r in GIVEN_SINGLES], fontsize=8)
    ax.set_xlabel("Spearman ρ (given %pos, not re-audited)")
    ax.set_title("PR #459 GSE123902+GSE131907 %pos (given)")
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _honest_pairs(df: pd.DataFrame, min_n: int = 8) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    out = df[df["n_patients"] >= min_n].copy()
    return out.sort_values(["p_meta", "n_patients"], ascending=[True, False], na_position="last")


def _top_md(df: pd.DataFrame, n=12) -> str:
    df = _honest_pairs(df)
    if df is None or df.empty:
        return "_No pairs passed the detect gate with honest n≥8._\n"
    lines = [
        "| pair | class | k | n_patients | mean Δ | p_meta | p_Wilcoxon | I² | units |",
        "|---|---|---:|---:|---:|---|---|---:|---|",
    ]
    for r in df.head(n).itertuples():
        lines.append(
            f"| {r.interaction_name} | {r.lr_class} | {int(r.k_units)} | {int(r.n_patients)} | "
            f"{fmt_num(r.mean_delta)} | {fmt_p(r.p_meta)} | {fmt_p(r.p_wilcoxon_patients)} | "
            f"{fmt_num(r.i2, 0) if np.isfinite(r.i2) else 'NA'} | {r.units} |"
        )
    return "\n".join(lines) + "\n"


def _top_between(df: pd.DataFrame, n=10) -> str:
    if df is None or df.empty:
        return "_No between-patient Q4 vs Q1 pairs._\n"
    df = df[df["n_compared"] >= 8].copy() if "n_compared" in df.columns else df
    if df.empty:
        return "_No between-patient Q4 vs Q1 pairs with n_compared≥8._\n"
    lines = [
        "| pair | class | k | n_Q1/n_Q4 | mean Δ | p_meta | I² | units |",
        "|---|---|---:|---|---:|---|---:|---|",
    ]
    for r in df.head(n).itertuples():
        lines.append(
            f"| {r.interaction_name} | {r.lr_class} | {int(r.k_units)} | {int(r.n_q1)}/{int(r.n_q4)} | "
            f"{fmt_num(r.mean_delta)} | {fmt_p(r.p_meta)} | "
            f"{fmt_num(r.i2, 0) if np.isfinite(r.i2) else 'NA'} | {r.units} |"
        )
    return "\n".join(lines) + "\n"


def write_finding(coverage, meta_q4, meta_med, meta_between, meta_q4_liana, meta_between_liana, skip_notes, out: Path) -> None:
    n_q4 = int(coverage["eligible_q4q1"].sum()) if "eligible_q4q1" in coverage else 0
    n_med = int(coverage["eligible_median"].sum()) if "eligible_median" in coverage else 0
    n_all = int(coverage["eligible_allmal"].sum()) if "eligible_allmal" in coverage else 0
    units_q4 = (
        sorted(coverage.loc[coverage["eligible_q4q1"] == 1, "unit"].unique())
        if not coverage.empty and "eligible_q4q1" in coverage
        else []
    )
    cov_lines = [
        "| unit | locked n | all-mal LR | Q4-eligible | median-eligible | note |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for unit in UNITS:
        sub = coverage[coverage["unit"] == unit] if not coverage.empty else pd.DataFrame()
        locked_n = int(sub["patient_id"].nunique()) if not sub.empty else 0
        q4 = int(sub["eligible_q4q1"].sum()) if not sub.empty and "eligible_q4q1" in sub else 0
        med = int(sub["eligible_median"].sum()) if not sub.empty and "eligible_median" in sub else 0
        am = int(sub["eligible_allmal"].sum()) if not sub.empty and "eligible_allmal" in sub else 0
        note = skip_notes.get(unit, "")
        cov_lines.append(f"| {unit} | {locked_n} | {am} | {q4} | {med} | {note} |")

    text = f"""# FINDING — CLDN4-only high-end CellChat / LIANA on GSE123902 + GSE131907

ADDITIVE. **CLDN4 only. No dual-high.** Patient/donor is the unit
(GSE123902 donor; GSE131907 sample). Do **not** add GSE148071. Do
**not** pile the seven-cohort pool.

The pair %pos Spearman is **taken as given** from PR #459 and is **not
re-audited**:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE131907 | %pos | 34 | −0.575 (0.001, 0%) | −0.700 (0.012; 10 vs 8) |

Singles (given, same PR): GSE123902 n=13 ρ=−0.659; GSE131907 n=21 ρ=−0.522.

This folder adds **outgoing CLDN4-high malignant → same-unit T/NK**:

- CellChat-style: Jin et al. 2021 Hill probability on CellChatDB v2 protein
  pairs (10% truncated mean, Kh=0.5, expr_prop ≥ 0.10).
- LIANA-style: CellPhoneDB mean-of-means on the **same** complexes
  (Efremova 2020 / Garcia-Alonso 2022). Partner = complex mean; pair =
  0.5 × (L + R).

CellChat R and LIANA `mt.cellphonedb` permutations were not run.

Primary high-end: **within-unit** malignant CLDN4 Q4 vs Q1. Sensitivity:
median split. Extra: between-unit Q4 vs Q1 using the **given** %pos
scores (same vectors as the n=34 row).

## Honest n

Given combo n=34 is **not** the LR n. Both CLDN4-high and CLDN4-low
malignant arms plus T/NK must meet the cell floor (Q4 vs Q1: n_mal≥{MIN_MAL_Q4}
and ≥{MIN_CELLS_ARM}/arm; median: ≥{2 * MIN_CELLS_ARM} malignant and
≥{MIN_CELLS_ARM} T/NK; all-mal outgoing: n_mal≥{MIN_CELLS_COMP} and
n_tnk≥{MIN_CELLS_COMP}).

Within-unit Q4 vs Q1 eligible: **n={n_q4}**. Median-split eligible: **n={n_med}**.
All-malignant outgoing (between-unit extra): **n={n_all}**.
Units with ≥1 Q4-eligible patient: {', '.join(units_q4) if units_q4 else 'none'}.

{chr(10).join(cov_lines)}

GSE123902 malignant = marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and
PTPRC == 0). GSE131907 malignant = author `Malignant cells` (not tS*;
tLung primaries are almost unlabeled as malignant and are not in the
locked n=21). TACSTD2 is never a gate.

## Primary: CellChat-style within-unit Q4 vs Q1 ΔP (outgoing Mal → T/NK)

ΔP = P(CLDN4-high mal → same-unit T/NK) − P(CLDN4-low mal → T/NK).
Each patient/donor/sample is one delta. Unit means are random-effects
pooled (DerSimonian–Laird). p-values are descriptive.

{_top_md(meta_q4)}

Primary table (`results/lr_table.tsv`) keeps pairs with **n_patients ≥ 8**.
Tiny-n rows stay in `results/lr_meta_q4q1.tsv` and are not the claim.

## LIANA-style within-unit Q4 vs Q1 ΔS (outgoing Mal → T/NK)

Same patients, same detect gate, CellPhoneDB-style mean-of-means
S = 0.5 × (L + R). ΔS = S_high − S_low.

{_top_md(meta_q4_liana)}

LIANA table (`results/lr_table_liana.tsv`) keeps pairs with **n_patients ≥ 8**.
Full meta: `results/lr_meta_q4q1_liana.tsv`.

## Sensitivity: within-unit median split (CellChat)

{_top_md(meta_med)}

Full table: `results/lr_meta_median.tsv`.

## Extra: between-unit high-end (given %pos Q4 vs Q1)

Quartiles use the **given** malignant CLDN4 %pos scores from the locked
tables (the same vectors as the n=34 ρ=−0.575 row). Per-unit P (or S)
is scored from all malignant cells → that unit’s T/NK. This is not a
re-audit of the Spearman.

CellChat-style:

{_top_between(meta_between)}

LIANA-style:

{_top_between(meta_between_liana)}

Full tables: `results/lr_meta_between_q4q1.tsv`, `results/lr_meta_between_q4q1_liana.tsv`.

## Extra figures

- `figures/fig_given_combo_rho.png` — given %pos ρ (not re-audited)
- `figures/fig_honest_n.png` — locked vs LR-eligible n
- `figures/fig_n_per_patient.png` — per-unit malignant / T/NK floors
- `figures/fig_extra_ligand_table.png` — top within-unit Q4 ΔP pairs (CellChat)
- `figures/fig_liana_outgoing.png` — top within-unit Q4 ΔS pairs (LIANA)
- `figures/fig_forest_q4q1.png` — CellChat meta forest
- `figures/fig_forest_liana_q4q1.png` — LIANA meta forest
- `figures/fig_patient_delta_top.png` — patient ΔP by unit for the top CellChat pair
- `figures/fig_forest_between_q4q1.png` — between-unit high-end extra
- `figures/fig_forest_median.png` — median-split sensitivity

## What is not claimed

- The n=34 ρ=−0.575 / Q4 r=−0.700 row is given. It is not recomputed here.
- TACSTD2 is not used to define high/low. This is not dual-high.
- GSE148071 is not added (it dilutes). The seven-cohort pool is not piled.
- Cell-pooled permutations are not the test. Patient/donor/sample is the unit.
- CellChat R visualizations and LIANA `mt.cellphonedb` were not generated.
- GSE131907 tS1–tS3 tLung epithelium is not in the locked author-malignant n=21.

## Reproduce

```bash
python3 methods/pair_123902_131907_hiend_cldn4/scripts/download.py
python3 methods/pair_123902_131907_hiend_cldn4/scripts/analyze.py
```

Hill constants: trim={TRIM}, Kh={KH}, expr_prop={EXPR_PROP}.
LIANA-style: mean-of-means on the same CellChat complexes.
"""
    (out / "FINDING.md").write_text(text)
    (ROOT / "FINDING.md").write_text(text)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--raw", type=Path, default=Path("/tmp/pair_123902_131907_raw"))
    p.add_argument("--out", type=Path, default=ROOT)
    args = p.parse_args()
    out = args.out
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    locked = load_locked()
    print(
        "locked n",
        {u: len(locked[u]) for u in UNITS},
        "sum",
        sum(len(locked[u]) for u in UNITS),
        flush=True,
    )
    if sum(len(locked[u]) for u in UNITS) != GIVEN["n"]:
        print(
            f"WARNING: locked n={sum(len(locked[u]) for u in UNITS)} != given n={GIVEN['n']}",
            flush=True,
        )

    wanted = wanted_genes()
    loaders = {
        "GSE123902": load_gse123902,
        "GSE131907": load_gse131907,
    }
    skip_notes = {}
    coverage_all = []
    pairs_all = []
    for unit, loader in loaders.items():
        bundle = loader(args.raw, wanted, locked[unit])
        skip_notes[unit] = bundle.get("note", "")
        lr = load_lr(DB, bundle["genes"])
        print(f"  {unit} LR pairs in matrix: {len(lr)}", flush=True)
        cov, pairs = score_bundle(unit, bundle, lr, modes=["q4q1", "median", "tertile"])
        coverage_all.append(cov)
        pairs_all.append(pairs)
        del bundle

    coverage = pd.concat(coverage_all, ignore_index=True) if coverage_all else pd.DataFrame()
    pairs = pd.concat(pairs_all, ignore_index=True) if pairs_all else pd.DataFrame()
    coverage.to_csv(out / "results" / "patient_coverage.tsv", sep="\t", index=False)
    if not pairs.empty:
        pairs.to_csv(out / "results" / "patient_lr_long.tsv.gz", sep="\t", index=False, compression="gzip")

    meta_q4 = meta_patient_deltas(pairs, "q4q1", "delta") if not pairs.empty else pd.DataFrame()
    meta_med = meta_patient_deltas(pairs, "median", "delta") if not pairs.empty else pd.DataFrame()
    meta_ter = meta_patient_deltas(pairs, "tertile", "delta") if not pairs.empty else pd.DataFrame()
    meta_q4_liana = meta_patient_deltas(pairs, "q4q1", "delta_liana") if not pairs.empty else pd.DataFrame()
    meta_between = meta_q4q1_between(pairs, locked, "prob") if not pairs.empty else pd.DataFrame()
    meta_between_liana = meta_q4q1_between(pairs, locked, "liana_score") if not pairs.empty else pd.DataFrame()

    if not meta_q4.empty:
        meta_q4.to_csv(out / "results" / "lr_meta_q4q1.tsv", sep="\t", index=False)
        primary = _honest_pairs(meta_q4)
        primary.to_csv(out / "results" / "lr_table.tsv", sep="\t", index=False)
    if not meta_q4_liana.empty:
        meta_q4_liana.to_csv(out / "results" / "lr_meta_q4q1_liana.tsv", sep="\t", index=False)
        _honest_pairs(meta_q4_liana).to_csv(out / "results" / "lr_table_liana.tsv", sep="\t", index=False)
    if not meta_med.empty:
        meta_med.to_csv(out / "results" / "lr_meta_median.tsv", sep="\t", index=False)
    if not meta_ter.empty:
        meta_ter.to_csv(out / "results" / "lr_meta_tertile.tsv", sep="\t", index=False)
    if not meta_between.empty:
        meta_between.to_csv(out / "results" / "lr_meta_between_q4q1.tsv", sep="\t", index=False)
    if not meta_between_liana.empty:
        meta_between_liana.to_csv(out / "results" / "lr_meta_between_q4q1_liana.tsv", sep="\t", index=False)

    if not pairs.empty:
        unit_rows = []
        sub = pairs[(pairs["split"] == "q4q1") & (pairs["arm"] == "delta") & (pairs["detected_either"])]
        for (unit, name), g in sub.groupby(["unit", "interaction_name"]):
            st = unit_delta_stats(g["delta"])
            first = g.iloc[0]
            unit_rows.append(
                {
                    "unit": unit,
                    "interaction_name": name,
                    "pathway_name": first["pathway_name"],
                    "lr_class": first["lr_class"],
                    "score": "cellchat_hill",
                    **st,
                }
            )
            st_l = unit_delta_stats(g["delta_liana"])
            unit_rows.append(
                {
                    "unit": unit,
                    "interaction_name": name,
                    "pathway_name": first["pathway_name"],
                    "lr_class": first["lr_class"],
                    "score": "liana_mean",
                    **st_l,
                }
            )
        pd.DataFrame(unit_rows).to_csv(out / "results" / "lr_unit_q4q1.tsv", sep="\t", index=False)

    plot_given_rho(out / "figures" / "fig_given_combo_rho.png")
    plot_n(coverage, out / "figures" / "fig_honest_n.png")
    if not coverage.empty:
        plot_n_per_patient(coverage, out / "figures" / "fig_n_per_patient.png")
    primary = _honest_pairs(meta_q4)
    primary_liana = _honest_pairs(meta_q4_liana)
    plot_forest(primary, out / "figures" / "fig_forest_q4q1.png", "Within-unit Q4 vs Q1 outgoing ΔP (CellChat, RE meta, n≥8)")
    plot_ligand_table(primary, out / "figures" / "fig_extra_ligand_table.png", "Extra: top outgoing patient-ΔP pairs (CellChat Q4 vs Q1, n≥8)")
    plot_forest(primary_liana, out / "figures" / "fig_forest_liana_q4q1.png", "Within-unit Q4 vs Q1 outgoing ΔS (LIANA mean-of-means, n≥8)")
    plot_ligand_table(primary_liana, out / "figures" / "fig_liana_outgoing.png", "Extra: top outgoing patient-ΔS pairs (LIANA Q4 vs Q1, n≥8)")
    between_show = (
        meta_between[meta_between["n_compared"] >= 8]
        if not meta_between.empty and "n_compared" in meta_between.columns
        else meta_between
    )
    plot_forest(between_show, out / "figures" / "fig_forest_between_q4q1.png", "Between-unit given-%pos Q4 vs Q1 (extra, n≥8)")
    if not primary.empty:
        plot_patient_deltas(pairs, primary.iloc[0]["interaction_name"], "q4q1", out / "figures" / "fig_patient_delta_top.png")
    if not meta_med.empty:
        plot_forest(_honest_pairs(meta_med), out / "figures" / "fig_forest_median.png", "Within-unit median-split outgoing ΔP (sensitivity, n≥8)")

    summary = {
        "given_spearman": GIVEN,
        "locked_n": {u: int(len(locked[u])) for u in UNITS},
        "eligible_q4q1": int(coverage["eligible_q4q1"].sum()) if "eligible_q4q1" in coverage else 0,
        "eligible_median": int(coverage["eligible_median"].sum()) if "eligible_median" in coverage else 0,
        "eligible_allmal": int(coverage["eligible_allmal"].sum()) if "eligible_allmal" in coverage else 0,
        "n_meta_q4_pairs": int(len(meta_q4)),
        "n_meta_q4_liana_pairs": int(len(meta_q4_liana)),
        "n_meta_between_pairs": int(len(meta_between)),
        "skipped": skip_notes,
        "dual_high": False,
        "added_gse148071": False,
        "piled_seven": False,
    }
    (out / "results" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_finding(
        coverage,
        meta_q4,
        meta_med,
        meta_between,
        meta_q4_liana,
        meta_between_liana,
        skip_notes,
        out / "results",
    )
    print("done", json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
