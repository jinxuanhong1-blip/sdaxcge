#!/usr/bin/env python3
"""Quantify leftover 2024 lung IO series that were not in the first pass.

This script covers bulk/array leftovers that can be downloaded as small
processed tables. GSE241934 is handled separately because it is single-cell.
"""

from __future__ import annotations

import csv
import gzip
import math
import statistics
import tempfile
import urllib.request
from pathlib import Path


OUT = Path(__file__).resolve().parent / "data"
LOCAL_CACHE = Path("/tmp/geo2024b")
GENES = ("TACSTD2", "CLDN4")
ENSEMBL = {"ENSG00000184292": "TACSTD2", "ENSG00000189143": "CLDN4"}
URLS = {
    "GSE253564_fpkm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz",
    "GSE253564_meta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/matrix/GSE253564_series_matrix.txt.gz",
    "GSE255144_counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255144/suppl/GSE255144_raw_counts.txt.gz",
    "GSE260598_counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE260nnn/GSE260598/suppl/GSE260598_allData_ReadCounts.txt.gz",
    "GSE270711_fpkm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270711/suppl/GSE270711_matrix_mRNA.fpkm.txt.gz",
}
LOCAL_NAMES = {
    "GSE253564_fpkm": "GSE253564_FPKM.txt.gz",
    "GSE253564_meta": "GSE253564_series_matrix.txt.gz",
    "GSE255144_counts": "GSE255144_counts.txt.gz",
    "GSE260598_counts": "GSE260598_counts.txt.gz",
    "GSE270711_fpkm": "GSE270711_mRNA.fpkm.txt.gz",
}


def download(cache: Path, key: str) -> Path:
    local = LOCAL_CACHE / LOCAL_NAMES[key]
    if local.exists():
        return local
    path = cache / Path(URLS[key]).name
    if not path.exists():
        urllib.request.urlretrieve(URLS[key], path)
    return path


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


def series_chars(path: Path, prefix: str) -> dict[str, list[str]]:
    titles: list[str] = []
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            row = next(csv.reader([line], delimiter="\t"))
            if not row:
                continue
            if row[0] == "!Sample_title":
                titles = row[1:]
            elif row[0] == "!Sample_characteristics_ch1":
                key = row[1].split(":", 1)[0]
                fields[key] = [value.split(": ", 1)[-1] for value in row[1:]]
    fields["title"] = titles
    return fields


def summarize(rows: list[dict[str, object]], group: str, note: str) -> dict[str, object]:
    selected = rows if group == "all" else [row for row in rows if row["group"] == group]
    x = [float(row["TACSTD2"]) for row in selected]
    y = [float(row["CLDN4"]) for row in selected]
    return {
        "accession": selected[0]["accession"],
        "subset": group,
        "n": len(selected),
        "TACSTD2_median": statistics.median(x),
        "CLDN4_median": statistics.median(y),
        "spearman_rho": spearman(x, y),
        "note": note,
    }


def read_253564(path: Path, meta: Path) -> list[dict[str, object]]:
    fields = series_chars(meta, "GSE253564")
    arm = dict(zip(fields["title"], fields["treatment"]))
    found: dict[str, list[float]] = {}
    with gzip.open(path, "rt") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
        samples = header[2:]
        for row in csv.reader(handle, delimiter="\t"):
            if row and row[0] in GENES:
                found[row[0]] = [float(value) for value in row[2:]]
    return [
        {
            "accession": "GSE253564",
            "sample": sample,
            "group": arm[sample],
            "patient": sample,
            "TACSTD2": found["TACSTD2"][i],
            "CLDN4": found["CLDN4"][i],
        }
        for i, sample in enumerate(samples)
    ]


def read_255144(path: Path) -> list[dict[str, object]]:
    found: dict[str, list[float]] = {}
    with gzip.open(path, "rt") as handle:
        header = next(handle).rstrip("\n").split("\t")
        samples = header[1:]
        for line in handle:
            row = line.rstrip("\n").split("\t")
            gene = ENSEMBL.get(row[0])
            if gene:
                found[gene] = [float(value.replace(",", ".")) for value in row[1:]]
    return [
        {
            "accession": "GSE255144",
            "sample": sample,
            "group": "baseline" if sample.startswith("ADK17") else "hyperprogression",
            "patient": "single_HPD_patient",
            "TACSTD2": found["TACSTD2"][i],
            "CLDN4": found["CLDN4"][i],
        }
        for i, sample in enumerate(samples)
    ]


def read_260598(path: Path) -> list[dict[str, object]]:
    found: dict[str, list[float]] = {}
    with gzip.open(path, "rt") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
        samples = header[1:]
        for row in csv.reader(handle, delimiter="\t"):
            if row and row[0] in GENES:
                found[row[0]] = [float(value) for value in row[1:]]
    rows = []
    for i, sample in enumerate(samples):
        if not sample.startswith(("43728", "43881")):
            continue
        if "tumor" in sample and sample.endswith("neg"):
            group = "lung_tumor_nest_CD68neg"
        elif "stroma" in sample and sample.endswith("neg"):
            group = "lung_stroma_CD68neg"
        else:
            group = "lung_other"
        rows.append(
            {
                "accession": "GSE260598",
                "sample": sample,
                "group": group,
                "patient": sample.split("_")[0],
                "TACSTD2": found["TACSTD2"][i],
                "CLDN4": found["CLDN4"][i],
            }
        )
    return rows


def read_270711(path: Path) -> list[dict[str, object]]:
    found: dict[str, list[float]] = {}
    with gzip.open(path, "rt") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
        samples = header[1:]
        for row in csv.reader(handle, delimiter="\t"):
            if row and row[0] in GENES:
                found[row[0]] = [float(value) for value in row[1:]]
    return [
        {
            "accession": "GSE270711",
            "sample": sample,
            "group": "adjacent" if sample.startswith("C") else "tumor",
            "patient": sample[1:],
            "TACSTD2": found["TACSTD2"][i],
            "CLDN4": found["CLDN4"][i],
        }
        for i, sample in enumerate(samples)
    ]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="geo2024_leftover_") as tmp:
        cache = Path(tmp)
        rows_253 = read_253564(download(cache, "GSE253564_fpkm"), download(cache, "GSE253564_meta"))
        rows_255 = read_255144(download(cache, "GSE255144_counts"))
        rows_260 = read_260598(download(cache, "GSE260598_counts"))
        rows_270 = read_270711(download(cache, "GSE270711_fpkm"))

    all_rows = rows_253 + rows_255 + rows_260 + rows_270
    with (OUT / "leftover_series_target_values.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    summary = [
        summarize(rows_253, "all", "Pretreatment tumor FPKM; GEO has arm labels only, not MPR"),
        summarize(rows_253, "Arm1", "Durvalumab monotherapy arm"),
        summarize(rows_253, "Arm2", "Durvalumab plus SBRT arm"),
        summarize(rows_255, "baseline", "One-patient HPD model; three culture replicates"),
        summarize(rows_255, "hyperprogression", "Same patient after ICI hyperprogression"),
        summarize(rows_260, "lung_tumor_nest_CD68neg", "Leftover lung GeoMx; not an ICI-outcome cohort"),
        summarize(rows_260, "lung_stroma_CD68neg", "Leftover lung GeoMx stroma AOIs"),
        summarize(rows_270, "tumor", "Three surgical LUAD samples; not ICI-treated"),
        summarize(rows_270, "adjacent", "Matched adjacent tissue"),
    ]
    with (OUT / "leftover_series_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)


if __name__ == "__main__":
    main()
