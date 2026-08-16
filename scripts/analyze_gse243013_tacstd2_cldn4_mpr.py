#!/usr/bin/env python3
"""Extract TACSTD2/CLDN4 from GSE243013 processed immune MTX and test vs MPR.

Public processed matrix:
  GSE243013_NSCLC_immune_scRNA_counts.mtx.gz
  (cells x genes, Matrix Market coordinate, 1-based)

This is CD45+ immune scRNA only (T/NK, B, myeloid). It is not a tumor-epithelial
bulk matrix. TACSTD2 and CLDN4 are epithelial genes; low detection is expected
and is reported rather than hidden.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import statistics
import time
from collections import Counter, defaultdict
from pathlib import Path


GENE_COL = {"TACSTD2": 1057, "CLDN4": 11797}  # 1-based MTX column; verified from genes.csv
MPR_ANY = {"MPR", "pCR"}
KNOWN = {"MPR", "pCR", "non-MPR"}


def mannwhitney_u(x: list[float], y: list[float]) -> tuple[float, float, float]:
    """Two-sided Mann-Whitney U with normal approximation and tie correction.

    Returns (U_x, z, two_sided_p). p is nan if either group is empty or variance is 0.
    """
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return float("nan"), float("nan"), float("nan")
    ranks = _average_ranks(x + y)
    ux = sum(ranks[:nx]) - nx * (nx + 1) / 2.0
    uy = nx * ny - ux
    n = nx + ny
    ties = Counter(ranks)
    tie_term = sum(t * t * t - t for t in ties.values())
    var = nx * ny * (n + 1) / 12.0
    if n > 1:
        var -= nx * ny * tie_term / (12.0 * n * (n - 1))
    if var <= 0:
        return ux, float("nan"), float("nan")
    mu = nx * ny / 2.0
    # continuity correction toward the mean
    z = (ux - mu - math.copysign(0.5, ux - mu)) / math.sqrt(var)
    p = 2.0 * _norm_sf(abs(z))
    if p > 1.0:
        p = 1.0
    return ux, z, p


def _average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j + 2) / 2.0  # 1-based ranks
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def _norm_sf(z: float) -> float:
    """Standard normal survival function P(Z > z)."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def cliffs_delta(x: list[float], y: list[float]) -> float:
    """Cliff's delta: (P(x>y) - P(x<y)). Positive means x > y."""
    if not x or not y:
        return float("nan")
    gt = lt = 0
    for a in x:
        for b in y:
            if a > b:
                gt += 1
            elif a < b:
                lt += 1
    return (gt - lt) / (len(x) * len(y))


def auc_from_cliffs(delta: float) -> float:
    if math.isnan(delta):
        return float("nan")
    return 0.5 + 0.5 * delta


def median(xs: list[float]) -> float:
    if not xs:
        return float("nan")
    return statistics.median(xs)


def mean(xs: list[float]) -> float:
    if not xs:
        return float("nan")
    return statistics.fmean(xs)


def qtile(xs: list[float], p: float) -> float:
    if not xs:
        return float("nan")
    ys = sorted(xs)
    if len(ys) == 1:
        return ys[0]
    k = (len(ys) - 1) * p
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return ys[int(k)]
    return ys[f] * (c - k) + ys[c] * (k - f)


def load_patient_table(meta_path: Path) -> tuple[dict[str, dict], list[str]]:
    """Return (patient_record_by_sample, cell_sampleIDs aligned to MTX rows)."""
    patients: dict[str, dict] = {}
    cell_sample: list[str] = []
    cell_lib: list[float] = []
    with gzip.open(meta_path, "rt", newline="") as f:
        for row in csv.DictReader(f):
            sid = row["sampleID"]
            cell_sample.append(sid)
            try:
                lib = float(row["total_counts"])
            except ValueError:
                lib = float("nan")
            cell_lib.append(lib)
            if sid not in patients:
                patients[sid] = {
                    "sampleID": sid,
                    "n_cells": 0,
                    "pathological_response": row["pathological_response"],
                    "pathological_response_rate": row["pathological_response_rate"],
                    "radiological_response": row["radiological_response"],
                    "cancer_type": row["cancer_type"],
                    "gender": row["gender"],
                    "age": row["age"],
                    "smoking_history": row["smoking_history"],
                    "pre_treatment_staging": row["pre_treatment_staging"],
                    "anti_PD1_therapy": row["anti-PD1_therapy"],
                    "chemotherapy": row["chemotherapy"],
                    "targeted_therapy": row["targeted_therapy"],
                    "cycles": row["cycles"],
                    "sum_total_counts": 0.0,
                }
            patients[sid]["n_cells"] += 1
            if not math.isnan(lib):
                patients[sid]["sum_total_counts"] += lib
    return patients, cell_sample, cell_lib


def extract_gene_columns(mtx_path: Path, out_path: Path) -> dict:
    """Stream MTX and keep only TACSTD2/CLDN4 columns. Write TSV: cell_index gene count."""
    want = {str(v).encode(): name for name, v in GENE_COL.items()}
    n_lines = 0
    n_kept = 0
    header = None
    t0 = time.time()
    with gzip.open(mtx_path, "rb") as fin, open(out_path, "w", newline="") as fout:
        # skip comments / banner
        while True:
            line = fin.readline()
            if not line:
                raise RuntimeError("MTX ended before header")
            if line.startswith(b"%"):
                continue
            header = line.decode("ascii").strip().split()
            break
        if header is None or len(header) < 3:
            raise RuntimeError(f"bad MTX header: {header}")
        n_rows, n_cols, n_nz = (int(x) for x in header[:3])
        fout.write("cell_index_1based\tgene\tcount\n")
        for line in fin:
            n_lines += 1
            parts = line.split()
            if len(parts) < 3:
                continue
            gene_name = want.get(parts[1])
            if gene_name is None:
                continue
            n_kept += 1
            fout.write(f"{int(parts[0])}\t{gene_name}\t{float(parts[2])}\n")
            if n_lines % 50_000_000 == 0:
                elapsed = time.time() - t0
                print(
                    f"[extract] scanned {n_lines} nnz, kept {n_kept}, {elapsed:.1f}s",
                    flush=True,
                )
    return {
        "mtx_rows": n_rows,
        "mtx_cols": n_cols,
        "mtx_nnz_declared": n_nz,
        "nnz_scanned": n_lines,
        "nnz_kept": n_kept,
        "seconds": time.time() - t0,
    }


def aggregate(
    extract_path: Path,
    patients: dict[str, dict],
    cell_sample: list[str],
    cell_lib: list[float],
) -> None:
    for sid, rec in patients.items():
        rec["TACSTD2_sum"] = 0.0
        rec["CLDN4_sum"] = 0.0
        rec["TACSTD2_n_pos"] = 0
        rec["CLDN4_n_pos"] = 0
        rec["TACSTD2_cpm_sum"] = 0.0
        rec["CLDN4_cpm_sum"] = 0.0
    with open(extract_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            idx = int(row["cell_index_1based"]) - 1
            gene = row["gene"]
            count = float(row["count"])
            sid = cell_sample[idx]
            rec = patients[sid]
            rec[f"{gene}_sum"] += count
            if count > 0:
                rec[f"{gene}_n_pos"] += 1
            lib = cell_lib[idx]
            if lib and not math.isnan(lib) and lib > 0:
                rec[f"{gene}_cpm_sum"] += count / lib * 1e6
    for rec in patients.values():
        n = rec["n_cells"]
        rec["TACSTD2_frac_pos"] = rec["TACSTD2_n_pos"] / n if n else float("nan")
        rec["CLDN4_frac_pos"] = rec["CLDN4_n_pos"] / n if n else float("nan")
        rec["TACSTD2_mean_count"] = rec["TACSTD2_sum"] / n if n else float("nan")
        rec["CLDN4_mean_count"] = rec["CLDN4_sum"] / n if n else float("nan")
        rec["TACSTD2_mean_cpm"] = rec["TACSTD2_cpm_sum"] / n if n else float("nan")
        rec["CLDN4_mean_cpm"] = rec["CLDN4_cpm_sum"] / n if n else float("nan")
        rec["mpr_any"] = (
            "yes"
            if rec["pathological_response"] in MPR_ANY
            else ("no" if rec["pathological_response"] == "non-MPR" else "unknown")
        )


def group_values(patients: dict[str, dict], field: str, group: str) -> list[float]:
    out = []
    for rec in patients.values():
        if rec["pathological_response"] not in KNOWN:
            continue
        if group == "mpr_any" and rec["mpr_any"] != "yes":
            continue
        if group == "non-MPR" and rec["pathological_response"] != "non-MPR":
            continue
        if group in KNOWN and rec["pathological_response"] != group:
            continue
        out.append(float(rec[field]))
    return out


def summarize_field(patients: dict[str, dict], field: str) -> list[dict]:
    rows = []
    pairs = [
        ("mpr_any_vs_nonMPR", "mpr_any", "non-MPR"),
        ("pCR_vs_nonMPR", "pCR", "non-MPR"),
        ("MPR_excl_pCR_vs_nonMPR", "MPR", "non-MPR"),
        ("pCR_vs_MPR_excl_pCR", "pCR", "MPR"),
    ]
    for contrast, a, b in pairs:
        xa = group_values(patients, field, a)
        xb = group_values(patients, field, b)
        u, z, p = mannwhitney_u(xa, xb)
        delta = cliffs_delta(xa, xb)
        rows.append(
            {
                "metric": field,
                "contrast": contrast,
                "n_a": len(xa),
                "n_b": len(xb),
                "median_a": median(xa),
                "median_b": median(xb),
                "mean_a": mean(xa),
                "mean_b": mean(xb),
                "q25_a": qtile(xa, 0.25),
                "q75_a": qtile(xa, 0.75),
                "q25_b": qtile(xb, 0.25),
                "q75_b": qtile(xb, 0.75),
                "mannwhitney_U_a": u,
                "mannwhitney_z": z,
                "mannwhitney_p_twosided": p,
                "cliffs_delta_a_minus_b": delta,
                "auc_a_gt_b": auc_from_cliffs(delta),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if not rows:
        path.write_text("")
        return
    fieldnames = fieldnames or list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            out = {}
            for k in fieldnames:
                v = row.get(k, "")
                if isinstance(v, float):
                    out[k] = "" if math.isnan(v) else f"{v:.10g}"
                else:
                    out[k] = v
            w.writerow(out)


def fmt(x: float, digits: int = 4) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if abs(x) >= 0.001 or x == 0:
        return f"{x:.{digits}g}"
    return f"{x:.3e}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--mtx", required=True)
    ap.add_argument("--extract", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--skip-extract", action="store_true")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    patients, cell_sample, cell_lib = load_patient_table(Path(args.meta))
    extract_info = None
    extract_path = Path(args.extract)
    if args.skip_extract and extract_path.exists():
        extract_info = {"skipped": True, "path": str(extract_path)}
    else:
        extract_info = extract_gene_columns(Path(args.mtx), extract_path)
    aggregate(extract_path, patients, cell_sample, cell_lib)

    patient_rows = [patients[k] for k in sorted(patients)]
    write_tsv(
        outdir / "patient_level_tacstd2_cldn4.tsv",
        patient_rows,
        [
            "sampleID",
            "pathological_response",
            "mpr_any",
            "pathological_response_rate",
            "radiological_response",
            "cancer_type",
            "gender",
            "age",
            "smoking_history",
            "pre_treatment_staging",
            "anti_PD1_therapy",
            "chemotherapy",
            "targeted_therapy",
            "cycles",
            "n_cells",
            "sum_total_counts",
            "TACSTD2_n_pos",
            "TACSTD2_frac_pos",
            "TACSTD2_sum",
            "TACSTD2_mean_count",
            "TACSTD2_mean_cpm",
            "CLDN4_n_pos",
            "CLDN4_frac_pos",
            "CLDN4_sum",
            "CLDN4_mean_count",
            "CLDN4_mean_cpm",
        ],
    )

    metrics = [
        "TACSTD2_frac_pos",
        "TACSTD2_mean_count",
        "TACSTD2_mean_cpm",
        "CLDN4_frac_pos",
        "CLDN4_mean_count",
        "CLDN4_mean_cpm",
    ]
    stat_rows = []
    for m in metrics:
        stat_rows.extend(summarize_field(patients, m))
    write_tsv(outdir / "stats_tacstd2_cldn4_vs_mpr.tsv", stat_rows)

    known = [r for r in patient_rows if r["pathological_response"] in KNOWN]
    resp_n = Counter(r["pathological_response"] for r in patient_rows)
    global_pos = {
        "TACSTD2_cells_pos": sum(r["TACSTD2_n_pos"] for r in patient_rows),
        "CLDN4_cells_pos": sum(r["CLDN4_n_pos"] for r in patient_rows),
        "TACSTD2_patients_any_pos": sum(1 for r in known if r["TACSTD2_n_pos"] > 0),
        "CLDN4_patients_any_pos": sum(1 for r in known if r["CLDN4_n_pos"] > 0),
        "n_cells_total": sum(r["n_cells"] for r in patient_rows),
        "n_patients_total": len(patient_rows),
        "n_patients_known_response": len(known),
    }

    run_info = {
        "dataset": "GSE243013",
        "matrix": "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz",
        "matrix_orientation": "cells x genes (1254749 x 31831)",
        "compartment": "CD45+ immune cells only (T/NK, B, myeloid)",
        "genes": GENE_COL,
        "response_counts_all_metadata_patients": dict(resp_n),
        "paper_n_patients": 234,
        "geo_metadata_n_sampleIDs": len(patient_rows),
        "excluded_from_tests": [
            r["sampleID"]
            for r in patient_rows
            if r["pathological_response"] not in KNOWN
        ],
        "mpr_definition_used": (
            "GEO pathological_response MPR or pCR = MPR-any; "
            "paper defines MPR as residual viable tumor <=10% and pCR as a subset of MPR"
        ),
        "extract": extract_info,
        "detection": global_pos,
    }
    (outdir / "run_info.json").write_text(json.dumps(run_info, indent=2) + "\n")
    print(json.dumps(run_info, indent=2))


if __name__ == "__main__":
    main()
