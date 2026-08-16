#!/usr/bin/env python3
"""Reproduce the lightweight TACSTD2/CLDN4 checks outside GSE205335.

Only Python's standard library is required. Source files are downloaded to a
temporary cache and are not committed.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
import statistics
import tempfile
import urllib.request
from pathlib import Path


OUT = Path(__file__).resolve().parent / "data"
GENES = ("TACSTD2", "CLDN4")
URLS = {
    "GSE248249_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248249/matrix/GSE248249_series_matrix.txt.gz",
    "GSE225620_counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225620/suppl/GSE225620_pre_vs_post_featureCounts.txt.gz",
    "GSE260770_fpkm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE260nnn/GSE260770/suppl/GSE260770_mRNA_FPKM.txt.gz",
    "GSE260770_meta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE260nnn/GSE260770/matrix/GSE260770_series_matrix.txt.gz",
    "GSE265899_q3": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE265nnn/GSE265899/suppl/GSE265899_Q3Norm.csv.gz",
    "GSE265899_meta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE265nnn/GSE265899/matrix/GSE265899_series_matrix.txt.gz",
}
PROBES_248249 = {"TACSTD2": "TC0100014340.hg.1", "CLDN4": "TC0700007993.hg.1"}


def download(cache: Path, key: str) -> Path:
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


def series_metadata(path: Path) -> tuple[list[str], list[str]]:
    titles: list[str] = []
    groups: list[str] = []
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            row = next(csv.reader([line], delimiter="\t"))
            if row and row[0] == "!Sample_title":
                titles = row[1:]
            elif row and row[0] == "!Sample_characteristics_ch1":
                values = [value.removeprefix("group: ") for value in row[1:]]
                if any(value != values[0] for value in values[1:]):
                    groups = values
    return titles, groups


def read_248249(path: Path) -> list[dict[str, object]]:
    titles: list[str] = []
    samples: list[str] = []
    values: dict[str, list[float]] = {}
    in_table = False
    reverse = {probe: gene for gene, probe in PROBES_248249.items()}
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            row = next(csv.reader([line], delimiter="\t"))
            if not row:
                continue
            if row[0] == "!Sample_title":
                titles = row[1:]
            elif row[0] == "!Sample_geo_accession":
                samples = row[1:]
            elif row[0] == "!series_matrix_table_begin":
                in_table = True
            elif in_table and row[0] in reverse:
                values[reverse[row[0]]] = [float(value) for value in row[1:]]
    rows = []
    for i, (sample, title) in enumerate(zip(samples, titles)):
        patient = re.search(r"Patient (\d+)", title).group(1)
        timepoint = "pre" if "pre-immunotherapy" in title else "post"
        rows.append(
            {
                "accession": "GSE248249",
                "sample": sample,
                "group": timepoint,
                "patient": patient,
                "TACSTD2": values["TACSTD2"][i],
                "CLDN4": values["CLDN4"][i],
            }
        )
    return rows


def read_gene_table(
    path: Path,
    accession: str,
    delimiter: str,
    symbol_column: int,
    value_start: int,
    labels: list[str] | None = None,
    groups: list[str] | None = None,
) -> list[dict[str, object]]:
    header: list[str] = []
    found: dict[str, list[float]] = {}
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        for row in reader:
            if len(row) > symbol_column and row[symbol_column] in GENES:
                found[row[symbol_column]] = [float(value) for value in row[value_start:]]
    sample_names = header[value_start:] if labels is None else labels
    if groups is None:
        groups = [
            "pre" if name.startswith("pre") else "post" if name.startswith("post") else "normal"
            for name in sample_names
        ]
    return [
        {
            "accession": accession,
            "sample": sample,
            "group": groups[i],
            "patient": "",
            "TACSTD2": found["TACSTD2"][i],
            "CLDN4": found["CLDN4"][i],
        }
        for i, sample in enumerate(sample_names)
    ]


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
        "TACSTD2_nonzero_fraction": sum(value > 0 for value in x) / len(x),
        "CLDN4_nonzero_fraction": sum(value > 0 for value in y) / len(y),
        "spearman_rho": spearman(x, y),
        "note": note,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="geo2024_") as tmp:
        cache = Path(tmp)
        rows_248 = read_248249(download(cache, "GSE248249_matrix"))
        rows_225 = read_gene_table(
            download(cache, "GSE225620_counts"), "GSE225620", "\t", 0, 1
        )
        titles_260, groups_260 = series_metadata(download(cache, "GSE260770_meta"))
        rows_260 = read_gene_table(
            download(cache, "GSE260770_fpkm"),
            "GSE260770",
            "\t",
            1,
            2,
        )
        group_map_260 = dict(zip(titles_260, groups_260))
        for row in rows_260:
            row["group"] = group_map_260[row["sample"]]
        titles_265, groups_265 = series_metadata(download(cache, "GSE265899_meta"))
        rows_265 = read_gene_table(
            download(cache, "GSE265899_q3"),
            "GSE265899",
            ",",
            0,
            1,
        )
        # The informative spatial split is encoded in the sample title.
        for row in rows_265:
            row["group"] = "tumor" if str(row["sample"]).startswith("tumor") else "immune"

    all_rows = rows_248 + rows_225 + rows_260 + rows_265
    with (OUT / "other_series_target_values.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    summary = [
        summarize(rows_248, "all", "RMA signal; nonzero fraction is not a detection rate"),
        summarize(rows_248, "pre", "baseline tumors; 13 patients"),
        summarize(rows_248, "post", "acquired-resistance tumors; 29 patients"),
        summarize(rows_225, "all", "whole blood raw counts; epithelial-marker readout is not credible"),
        summarize(rows_260, "Responsed to sintilimab", "peripheral-blood FPKM"),
        summarize(rows_260, "Non-responsed to sintilimab", "peripheral-blood FPKM"),
        summarize(rows_265, "tumor", "context-only LUAD tumor AOIs; patients were not an ICI cohort"),
        summarize(rows_265, "immune", "context-only LUAD immune AOIs; patients were not an ICI cohort"),
    ]

    paired = {}
    for row in rows_248:
        paired.setdefault(row["patient"], {})[row["group"]] = row
    pairs = [pair for pair in paired.values() if set(pair) == {"pre", "post"}]
    for gene in GENES:
        changes = [float(pair["post"][gene]) - float(pair["pre"][gene]) for pair in pairs]
        positive = sum(change > 0 for change in changes)
        summary.append(
            {
                "accession": "GSE248249",
                "subset": f"paired_post_minus_pre_{gene}",
                "n": len(changes),
                "TACSTD2_median": statistics.median(changes) if gene == "TACSTD2" else "",
                "CLDN4_median": statistics.median(changes) if gene == "CLDN4" else "",
                "TACSTD2_nonzero_fraction": "",
                "CLDN4_nonzero_fraction": "",
                "spearman_rho": "",
                "note": f"{positive}/{len(changes)} changes were positive; descriptive only",
            }
        )

    with (OUT / "other_series_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)


if __name__ == "__main__":
    main()
