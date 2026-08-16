#!/usr/bin/env python3
"""GSE241934: residual epithelial TACSTD2/CLDN4 vs T/NK (IIT and Real kept separate).

Streams only TACSTD2 and CLDN4 from the 10x MTX. Cell types come from the
author meta (`major.cell.type`: Epi / T / NK). Residual-tumor epithelium is
the malignant/epithelial compartment available on GEO (no separate CNV call).
"""
from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from pathlib import Path

from lib_stats import spearman

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/"
GENES = ("TACSTD2", "CLDN4")
MIN_EPI = 20
MIN_TNK = 20


def gene_rows(features_path: Path) -> dict[int, str]:
    wanted: dict[int, str] = {}
    with gzip.open(features_path, "rt") as handle:
        for i, line in enumerate(handle, start=1):
            gene = line.split("\t", 1)[0].strip()
            if gene in GENES:
                wanted[i] = gene
    missing = set(GENES) - set(wanted.values())
    if missing:
        raise SystemExit(f"missing genes in {features_path}: {missing}")
    return wanted


def load_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def stream_genes(mtx_gz: Path, wanted_rows: dict[int, str], n_barcodes: int) -> dict[str, list[float]]:
    counts = {gene: [0.0] * n_barcodes for gene in GENES}
    with gzip.open(mtx_gz, "rt") as handle:
        header_seen = False
        for line in handle:
            if line.startswith("%"):
                continue
            row_s, col_s, val_s = line.split()
            if not header_seen:
                header_seen = True
                continue
            row = int(row_s)
            if row not in wanted_rows:
                continue
            counts[wanted_rows[row]][int(col_s) - 1] += float(val_s)
    return counts


def load_meta(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def patient_rows(meta, barcodes, counts) -> list[dict]:
    import math

    idx = {b: i for i, b in enumerate(barcodes)}
    buckets: dict[str, dict] = {}
    for cell in meta:
        pid = cell["orig.ident"]
        b = buckets.get(pid)
        if b is None:
            b = {
                "patient": pid,
                "n_cells": 0,
                "n_epi": 0,
                "n_tnk": 0,
                "n_cd8": 0,
                "epi_lib": 0.0,
                "TACSTD2_sum": 0.0,
                "CLDN4_sum": 0.0,
                "TACSTD2_pos": 0,
                "CLDN4_pos": 0,
                "TACSTD2_log1p_cp10k_sum": 0.0,
                "CLDN4_log1p_cp10k_sum": 0.0,
                "pathology": cell.get("Pathological Response", ""),
                "PD1": cell.get("PD1", ""),
                "EGFR": cell.get("EGFR", ""),
            }
            buckets[pid] = b
        b["n_cells"] += 1
        major = cell["major.cell.type"]
        if major in {"T", "NK"}:
            b["n_tnk"] += 1
        ctype = cell.get("cell.type", "")
        if ctype.startswith("CD8"):
            b["n_cd8"] += 1
        if major != "Epi":
            continue
        i = idx[cell["cellID"]]
        t = counts["TACSTD2"][i]
        c = counts["CLDN4"][i]
        lib = float(cell["nCount_RNA"])
        b["n_epi"] += 1
        b["epi_lib"] += lib
        b["TACSTD2_sum"] += t
        b["CLDN4_sum"] += c
        b["TACSTD2_pos"] += int(t > 0)
        b["CLDN4_pos"] += int(c > 0)
        if lib > 0:
            b["TACSTD2_log1p_cp10k_sum"] += math.log1p(1e4 * t / lib)
            b["CLDN4_log1p_cp10k_sum"] += math.log1p(1e4 * c / lib)

    rows = []
    for pid, b in sorted(buckets.items()):
        n_epi = b["n_epi"]
        lib = b["epi_lib"]
        rows.append(
            {
                "patient": pid,
                "n_cells": b["n_cells"],
                "n_epi": n_epi,
                "n_tnk": b["n_tnk"],
                "n_cd8": b["n_cd8"],
                "frac_tnk": b["n_tnk"] / b["n_cells"] if b["n_cells"] else float("nan"),
                "frac_cd8": b["n_cd8"] / b["n_cells"] if b["n_cells"] else float("nan"),
                "TACSTD2_log1p_cp10k": (b["TACSTD2_log1p_cp10k_sum"] / n_epi) if n_epi else float("nan"),
                "CLDN4_log1p_cp10k": (b["CLDN4_log1p_cp10k_sum"] / n_epi) if n_epi else float("nan"),
                "TACSTD2_pct_pos": (b["TACSTD2_pos"] / n_epi) if n_epi else float("nan"),
                "CLDN4_pct_pos": (b["CLDN4_pos"] / n_epi) if n_epi else float("nan"),
                "TACSTD2_sum": b["TACSTD2_sum"],
                "CLDN4_sum": b["CLDN4_sum"],
                "epi_lib": lib,
                "pathology": b["pathology"],
                "PD1": b["PD1"],
                "EGFR": b["EGFR"],
                "eligible": int(n_epi >= MIN_EPI and b["n_tnk"] >= MIN_TNK),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def cohort_stats(rows: list[dict], cohort: str) -> list[dict]:
    elig = [r for r in rows if r["eligible"]]
    out = []
    for gene, expr in (
        ("TACSTD2", "TACSTD2_log1p_cp10k"),
        ("TACSTD2", "TACSTD2_pct_pos"),
        ("CLDN4", "CLDN4_log1p_cp10k"),
        ("CLDN4", "CLDN4_pct_pos"),
    ):
        rho, p, n = spearman([r[expr] for r in elig], [r["frac_tnk"] for r in elig])
        out.append(
            {
                "cohort": cohort,
                "gene": gene,
                "metric": expr,
                "immune": "frac_tnk",
                "n": n,
                "rho": rho,
                "p": p,
            }
        )
        rho, p, n = spearman([r[expr] for r in elig], [r["frac_cd8"] for r in elig])
        out.append(
            {
                "cohort": cohort,
                "gene": gene,
                "metric": expr,
                "immune": "frac_cd8",
                "n": n,
                "rho": rho,
                "p": p,
            }
        )
    return out


def run_one(data: Path, prefix: str, feat: str, barc: str, mtx: str, meta: str, cohort: str, outdir: Path):
    print(f"== {cohort} ==", flush=True)
    wanted = gene_rows(data / feat)
    barcodes = load_barcodes(data / barc)
    print(f"  barcodes={len(barcodes)} gene_rows={wanted}", flush=True)
    counts = stream_genes(data / mtx, wanted, len(barcodes))
    meta_rows = load_meta(data / meta)
    rows = patient_rows(meta_rows, barcodes, counts)
    for r in rows:
        r["cohort"] = cohort
    write_tsv(outdir / f"{cohort}_patients.tsv", rows)
    stats_rows = cohort_stats(rows, cohort)
    write_tsv(outdir / f"{cohort}_spearman.tsv", stats_rows)
    print("  patients", len(rows), "eligible", sum(r["eligible"] for r in rows), flush=True)
    for s in stats_rows:
        print(f"  {s['gene']} {s['metric']} vs {s['immune']}: n={s['n']} ρ={s['rho']:.3f} p={s['p']:.4g}", flush=True)
    return rows, stats_rows


def main() -> None:
    data = Path("/tmp/scrna_meta_tnk/GSE241934")
    outdir = Path(__file__).resolve().parent / "results" / "cohorts"
    outdir.mkdir(parents=True, exist_ok=True)
    run_one(
        data,
        "IIT",
        "GSE241934_IIT_features.tsv.gz",
        "GSE241934_IIT_barcodes.tsv.gz",
        "GSE241934_IIT_Matrix.mtx.gz",
        "GSE241934_IIT_Meta.txt.gz",
        "GSE241934_IIT",
        outdir,
    )
    run_one(
        data,
        "Real",
        "GSE241934_RWC_features.tsv.gz",
        "GSE241934_RWC_barcodes.tsv.gz",
        "GSE241934_Real_Matrix.mtx.gz",
        "GSE241934_Real_Meta.txt.gz",
        "GSE241934_Real",
        outdir,
    )


if __name__ == "__main__":
    main()
