#!/usr/bin/env python3
"""Patient-level TACSTD2/CLDN4 in GSE241934 residual epithelial cells.

Leftover 2024 GEO series: post-neoadjuvant PD-1 + chemotherapy residual
tumors. IIT (EGFR-mutant, uniform sintilimab) and Real/RWC (EGFR-WT, mixed
PD-1 agents) are analyzed separately. Residual-tumor bias is explicit:
MPR/pCR patients have few leftover epithelial cells by definition.

Only the Python standard library is required. The script streams two gene
rows from the 10x MTX instead of loading the full sparse matrix.
"""
from __future__ import annotations

import csv
import gzip
import math
import statistics
import urllib.request
from itertools import combinations
from pathlib import Path


DATA = Path("/tmp/geo2024b")
OUT = Path(__file__).resolve().parent / "data"
GENES = ("TACSTD2", "CLDN4")
MPR_LABELS = {"MPR", "pCR"}
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/"
FILES = (
    "GSE241934_IIT_features.tsv.gz",
    "GSE241934_IIT_barcodes.tsv.gz",
    "GSE241934_IIT_Matrix.mtx.gz",
    "GSE241934_IIT_Meta.txt.gz",
    "GSE241934_RWC_features.tsv.gz",
    "GSE241934_RWC_barcodes.tsv.gz",
    "GSE241934_Real_Matrix.mtx.gz",
    "GSE241934_Real_Meta.txt.gz",
)


def ensure_inputs() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        path = DATA / name
        if path.exists() and path.stat().st_size > 0:
            continue
        print(f"downloading {name}", flush=True)
        urllib.request.urlretrieve(FTP + name, path)


def gene_rows(features_path: Path) -> dict[int, str]:
    wanted: dict[int, str] = {}
    with gzip.open(features_path, "rt") as handle:
        for i, line in enumerate(handle, start=1):
            gene = line.split("\t", 1)[0]
            if gene in GENES:
                wanted[i] = gene
    if set(wanted.values()) != set(GENES):
        raise SystemExit(f"missing genes in {features_path}: {wanted}")
    return wanted


def load_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def stream_two_genes_indexed(
    mtx_gz: Path, wanted_rows: dict[int, str], barcodes: list[str]
) -> dict[str, list[float]]:
    n = len(barcodes)
    counts = {gene: [0.0] * n for gene in GENES}
    with gzip.open(mtx_gz, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            row_s, col_s, val_s = line.split()
            row = int(row_s)
            if row not in wanted_rows:
                # header line is nrow ncol nnz
                if int(col_s) == n:
                    continue
                continue
            counts[wanted_rows[row]][int(col_s) - 1] += float(val_s)
    return counts


def load_meta(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def patient_table(
    meta: list[dict[str, str]],
    barcodes: list[str],
    counts: dict[str, list[float]],
    min_epi: int,
) -> list[dict[str, object]]:
    barcode_index = {barcode: i for i, barcode in enumerate(barcodes)}
    buckets: dict[str, dict[str, object]] = {}
    for cell in meta:
        if cell["major.cell.type"] != "Epi":
            continue
        pid = cell["orig.ident"]
        idx = barcode_index[cell["cellID"]]
        bucket = buckets.get(pid)
        if bucket is None:
            bucket = {
                "patientID": pid,
                "n_epithelial": 0,
                "libsize": 0.0,
                "TACSTD2_count": 0.0,
                "CLDN4_count": 0.0,
                "TACSTD2_pos": 0,
                "CLDN4_pos": 0,
                "pathology": cell["Pathological Response"],
                "PD1": cell["PD1"],
                "EGFR": cell["EGFR"],
            }
            buckets[pid] = bucket
        t = counts["TACSTD2"][idx]
        c = counts["CLDN4"][idx]
        bucket["n_epithelial"] = int(bucket["n_epithelial"]) + 1
        bucket["libsize"] = float(bucket["libsize"]) + float(cell["nCount_RNA"])
        bucket["TACSTD2_count"] = float(bucket["TACSTD2_count"]) + t
        bucket["CLDN4_count"] = float(bucket["CLDN4_count"]) + c
        bucket["TACSTD2_pos"] = int(bucket["TACSTD2_pos"]) + int(t > 0)
        bucket["CLDN4_pos"] = int(bucket["CLDN4_pos"]) + int(c > 0)

    rows = []
    for pid in sorted(buckets):
        bucket = buckets[pid]
        n_epi = int(bucket["n_epithelial"])
        lib = float(bucket["libsize"])
        t = float(bucket["TACSTD2_count"])
        c = float(bucket["CLDN4_count"])
        rows.append(
            {
                "patientID": pid,
                "n_epithelial": n_epi,
                "libsize": f"{lib:.1f}",
                "TACSTD2_count": f"{t:.1f}",
                "CLDN4_count": f"{c:.1f}",
                "TACSTD2_log2cpm": f"{math.log2(1e6 * t / lib + 1):.6f}" if lib > 0 else "",
                "CLDN4_log2cpm": f"{math.log2(1e6 * c / lib + 1):.6f}" if lib > 0 else "",
                "TACSTD2_detect": f"{int(bucket['TACSTD2_pos']) / n_epi:.6f}" if n_epi else "",
                "CLDN4_detect": f"{int(bucket['CLDN4_pos']) / n_epi:.6f}" if n_epi else "",
                "pathology": bucket["pathology"],
                "PD1": bucket["PD1"],
                "EGFR": bucket["EGFR"],
                "keep_min20": int(n_epi >= 20),
                "keep_min50": int(n_epi >= 50),
            }
        )
        if n_epi < min_epi:
            continue
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def pearson(x: list[float], y: list[float]) -> float:
    mx, my = statistics.mean(x), statistics.mean(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    den = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
    return num / den if den else math.nan


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2
        for index in order[start:end]:
            result[index] = rank
        start = end
    return result


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(ranks(x), ranks(y))


def exact_perm(x: list[float], y: list[float]) -> tuple[float, float, float]:
    n1 = len(x)
    pooled = x + y
    n = len(pooled)
    obs = statistics.mean(x) - statistics.mean(y)
    extreme = 0
    total = 0
    for idx in combinations(range(n), n1):
        chosen = set(idx)
        a = [pooled[i] for i in idx]
        b = [pooled[i] for i in range(n) if i not in chosen]
        if abs(statistics.mean(a) - statistics.mean(b)) >= abs(obs) - 1e-12:
            extreme += 1
        total += 1
    return statistics.mean(x), statistics.mean(y), extreme / total


def summarize(rows: list[dict[str, object]], cohort: str, min_epi: int) -> list[dict[str, object]]:
    use = [row for row in rows if int(row["n_epithelial"]) >= min_epi]
    mpr = [row for row in use if row["pathology"] in MPR_LABELS]
    non = [row for row in use if row["pathology"] not in MPR_LABELS]
    out = []
    if len(use) >= 3:
        rho = spearman(
            [float(row["TACSTD2_log2cpm"]) for row in use],
            [float(row["CLDN4_log2cpm"]) for row in use],
        )
        out.append(
            {
                "cohort": cohort,
                "min_epithelial": min_epi,
                "n_patients": len(use),
                "n_mpr_pcr": len(mpr),
                "n_non_mpr": len(non),
                "metric": "spearman_TACSTD2_CLDN4",
                "value": f"{rho:.6f}",
                "p_exact": "",
                "mpr_mean": "",
                "non_mean": "",
            }
        )
    if len(mpr) >= 2 and len(non) >= 2:
        for metric in ("TACSTD2_log2cpm", "CLDN4_log2cpm", "TACSTD2_detect", "CLDN4_detect"):
            x = [float(row[metric]) for row in mpr]
            y = [float(row[metric]) for row in non]
            mx, my, p = exact_perm(x, y)
            out.append(
                {
                    "cohort": cohort,
                    "min_epithelial": min_epi,
                    "n_patients": len(use),
                    "n_mpr_pcr": len(mpr),
                    "n_non_mpr": len(non),
                    "metric": metric,
                    "value": f"{mx - my:.6f}",
                    "p_exact": f"{p:.6f}",
                    "mpr_mean": f"{mx:.6f}",
                    "non_mean": f"{my:.6f}",
                }
            )
    return out


def analyze_cohort(
    label: str,
    features: Path,
    barcodes_path: Path,
    mtx: Path,
    meta_path: Path,
) -> list[dict[str, object]]:
    print(f"streaming {label} MTX", flush=True)
    wanted = gene_rows(features)
    barcodes = load_barcodes(barcodes_path)
    counts = stream_two_genes_indexed(mtx, wanted, barcodes)
    meta = load_meta(meta_path)
    rows = patient_table(meta, barcodes, counts, min_epi=20)
    for row in rows:
        row["cohort"] = label
        # move cohort first
    ordered = []
    for row in rows:
        ordered.append({"cohort": label, **{k: v for k, v in row.items() if k != "cohort"}})
    return ordered


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ensure_inputs()
    iit = analyze_cohort(
        "IIT_EGFRmut",
        DATA / "GSE241934_IIT_features.tsv.gz",
        DATA / "GSE241934_IIT_barcodes.tsv.gz",
        DATA / "GSE241934_IIT_Matrix.mtx.gz",
        DATA / "GSE241934_IIT_Meta.txt.gz",
    )
    write_tsv(OUT / "GSE241934_IIT_patient_targets.tsv", iit)
    rwc = analyze_cohort(
        "Real_EGFRWT",
        DATA / "GSE241934_RWC_features.tsv.gz",
        DATA / "GSE241934_RWC_barcodes.tsv.gz",
        DATA / "GSE241934_Real_Matrix.mtx.gz",
        DATA / "GSE241934_Real_Meta.txt.gz",
    )
    write_tsv(OUT / "GSE241934_Real_patient_targets.tsv", rwc)
    stats = []
    stats += summarize(iit, "IIT_EGFRmut", 20)
    stats += summarize(rwc, "Real_EGFRWT", 20)
    stats += summarize(rwc, "Real_EGFRWT", 50)
    write_tsv(OUT / "GSE241934_statistics.tsv", stats)
    print("IIT patients", len(iit), "Real patients", len(rwc), flush=True)
    for row in stats:
        print(row, flush=True)


if __name__ == "__main__":
    main()
