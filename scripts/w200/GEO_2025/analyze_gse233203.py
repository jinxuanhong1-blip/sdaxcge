#!/usr/bin/env python3
"""Reproduce the strict-2025 lung ICI TACSTD2/CLDN4 analysis.

GSE233203 is the only reviewed 2025 GEO series with all three requirements:
tumor-containing material, both genes on the assay, and a patient-level ICI
response label. The script downloads public GEO files, streams the 10x sparse
matrices, builds one pseudobulk observation per patient, and runs exact
two-sided permutation Mann-Whitney tests.

Only Python's standard library is required.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import itertools
import math
import statistics
import tarfile
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "w200" / "GEO_2025"
DATA = OUT / "data"
GSE = "GSE233203"
GENES = ("TACSTD2", "CLDN4")
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233203/matrix/"
    "GSE233203_series_matrix.txt.gz"
)
RAW_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233203/suppl/"
    "GSE233203_RAW.tar"
)


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "w200-geo-2025/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as out:
        while chunk := response.read(1024 * 1024):
            out.write(chunk)
    partial.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def parse_metadata(path: Path) -> dict[str, dict[str, str]]:
    fields: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                fields["gsm"] = next(csv.reader([line], delimiter="\t"))[1:]
            elif line.startswith("!Sample_title"):
                fields["title"] = next(csv.reader([line], delimiter="\t"))[1:]
            elif line.startswith("!Sample_characteristics"):
                characteristics.append(next(csv.reader([line], delimiter="\t"))[1:])
            elif line.startswith("!series_matrix_table_begin"):
                break
    if set(fields) != {"gsm", "title"}:
        raise ValueError("GEO matrix metadata is missing sample accessions or titles")
    rows = {
        gsm: {"gsm": gsm, "title": title}
        for gsm, title in zip(fields["gsm"], fields["title"], strict=True)
    }
    gsms = fields["gsm"]
    for characteristic_row in characteristics:
        for gsm, cell in zip(gsms, characteristic_row):
            key, separator, value = cell.partition(":")
            if separator:
                rows[gsm][key.strip().lower()] = value.strip()
    return rows


def target_rows(features_member: io.BufferedReader) -> dict[str, set[int]]:
    rows = {gene: set() for gene in GENES}
    with gzip.GzipFile(fileobj=features_member) as compressed:
        text = io.TextIOWrapper(compressed, encoding="utf-8", errors="replace")
        for row_number, line in enumerate(text, start=1):
            fields = line.rstrip("\n").split("\t")
            symbol = (fields[1] if len(fields) > 1 else fields[0]).upper()
            if symbol in rows:
                rows[symbol].add(row_number)
    missing = [gene for gene, indices in rows.items() if not indices]
    if missing:
        raise ValueError(f"Missing target feature(s): {', '.join(missing)}")
    return rows


def stream_matrix(
    matrix_member: io.BufferedReader, rows: dict[str, set[int]]
) -> tuple[int, float, dict[str, float], dict[str, int]]:
    inverse = {row: gene for gene, indices in rows.items() for row in indices}
    sums = {gene: 0.0 for gene in GENES}
    detected_cells = {gene: set() for gene in GENES}
    total_counts = 0.0
    dimensions_seen = False
    n_cells = 0
    with gzip.GzipFile(fileobj=matrix_member) as compressed:
        text = io.TextIOWrapper(compressed, encoding="ascii", errors="strict")
        for line in text:
            if line.startswith("%"):
                continue
            fields = line.split()
            if not dimensions_seen:
                if len(fields) != 3:
                    raise ValueError("Invalid MatrixMarket dimensions")
                _, n_cells, _ = map(int, fields)
                dimensions_seen = True
                continue
            row, column = int(fields[0]), int(fields[1])
            value = float(fields[2])
            total_counts += value
            gene = inverse.get(row)
            if gene is not None:
                sums[gene] += value
                detected_cells[gene].add(column)
    if not dimensions_seen or total_counts <= 0:
        raise ValueError("MatrixMarket file is empty or has no counts")
    return n_cells, total_counts, sums, {
        gene: len(columns) for gene, columns in detected_cells.items()
    }


def analyze_raw(raw_path: Path, metadata: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    with tarfile.open(raw_path) as archive:
        names = set(archive.getnames())
        prefixes = sorted(
            name.removesuffix("_matrix.mtx.gz")
            for name in names
            if name.endswith("_matrix.mtx.gz")
        )
        for prefix in prefixes:
            gsm = prefix.split("_", 1)[0]
            feature_name = f"{prefix}_features.tsv.gz"
            matrix_name = f"{prefix}_matrix.mtx.gz"
            if feature_name not in names:
                raise ValueError(f"Missing {feature_name}")
            feature_handle = archive.extractfile(feature_name)
            matrix_handle = archive.extractfile(matrix_name)
            if feature_handle is None or matrix_handle is None:
                raise ValueError(f"Could not read files for {prefix}")
            rows = target_rows(feature_handle)
            n_cells, total, sums, detected = stream_matrix(matrix_handle, rows)
            meta = metadata.get(gsm)
            if meta is None:
                raise ValueError(f"No GEO metadata found for {gsm}")
            row: dict[str, object] = {
                "gsm": gsm,
                "sample": meta["title"],
                "response": meta.get("therapeutic response", ""),
                "n_cells": n_cells,
                "total_counts": int(total),
            }
            for gene in GENES:
                cpm = sums[gene] / total * 1_000_000
                row[f"{gene}_counts"] = int(sums[gene])
                row[f"{gene}_log2_cpm1"] = math.log2(cpm + 1)
                row[f"{gene}_detected_cells"] = detected[gene]
                row[f"{gene}_detected_fraction"] = detected[gene] / n_cells
            results.append(row)
    return results


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position + 1
        while end < len(order) and values[order[end]] == values[order[position]]:
            end += 1
        rank = ((position + 1) + end) / 2
        for index in order[position:end]:
            ranks[index] = rank
        position = end
    return ranks


def exact_mann_whitney(group_a: list[float], group_b: list[float]) -> tuple[float, float]:
    """Return U for group A and an exact two-sided permutation p-value."""
    values = group_a + group_b
    ranks = average_ranks(values)
    n_a, n_b = len(group_a), len(group_b)
    observed = sum(ranks[:n_a]) - n_a * (n_a + 1) / 2
    center = n_a * n_b / 2
    extreme = 0
    total = 0
    tolerance = 1e-12
    for selected in itertools.combinations(range(len(values)), n_a):
        u_value = sum(ranks[index] for index in selected) - n_a * (n_a + 1) / 2
        extreme += abs(u_value - center) + tolerance >= abs(observed - center)
        total += 1
    return observed, extreme / total


def cliffs_delta(group_a: list[float], group_b: list[float]) -> float:
    greater = sum(a > b for a in group_a for b in group_b)
    less = sum(a < b for a in group_a for b in group_b)
    return (greater - less) / (len(group_a) * len(group_b))


def holm_adjust(p_values: list[float]) -> list[float]:
    adjusted = [0.0] * len(p_values)
    running = 0.0
    for rank, index in enumerate(sorted(range(len(p_values)), key=p_values.__getitem__)):
        value = min(1.0, (len(p_values) - rank) * p_values[index])
        running = max(running, value)
        adjusted[index] = running
    return adjusted


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_svg(path: Path, gene: str, rows: list[dict[str, object]]) -> None:
    groups = ("Response", "Non-response")
    values = {
        group: [
            float(row[f"{gene}_log2_cpm1"])
            for row in rows
            if row["response"] == group
        ]
        for group in groups
    }
    all_values = values[groups[0]] + values[groups[1]]
    lower, upper = min(all_values), max(all_values)
    padding = max(0.5, (upper - lower) * 0.15)
    lower, upper = lower - padding, upper + padding
    width, height = 560, 420
    left, right, top, bottom = 75, 25, 45, 65

    def y(value: float) -> float:
        return top + (upper - value) / (upper - lower) * (height - top - bottom)

    colors = ("#2f855a", "#c05621")
    x_positions = (210, 410)
    jitter = (-22, 0, 22, -11, 11)
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="24" text-anchor="middle" '
        f'font-family="sans-serif" font-size="16">GSE233203: {gene}</text>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" '
        'stroke="#333"/>',
        f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" '
        f'y2="{height-bottom}" stroke="#333"/>',
    ]
    for tick in range(math.floor(lower), math.ceil(upper) + 1):
        tick_y = y(tick)
        elements.extend(
            [
                f'<line x1="{left-5}" y1="{tick_y:.1f}" x2="{left}" '
                f'y2="{tick_y:.1f}" stroke="#333"/>',
                f'<text x="{left-10}" y="{tick_y+5:.1f}" text-anchor="end" '
                f'font-family="sans-serif" font-size="12">{tick}</text>',
            ]
        )
    for group_index, group in enumerate(groups):
        x_value = x_positions[group_index]
        median = statistics.median(values[group])
        for index, value in enumerate(values[group]):
            point_x = x_value + jitter[index]
            elements.append(
                f'<circle cx="{point_x}" cy="{y(value):.1f}" r="5" '
                f'fill="{colors[group_index]}" fill-opacity="0.8"/>'
            )
        elements.extend(
            [
                f'<line x1="{x_value-40}" y1="{y(median):.1f}" '
                f'x2="{x_value+40}" y2="{y(median):.1f}" '
                f'stroke="{colors[group_index]}" stroke-width="4"/>',
                f'<text x="{x_value}" y="{height-bottom+25}" text-anchor="middle" '
                f'font-family="sans-serif" font-size="13">{group} '
                f'(n={len(values[group])})</text>',
            ]
        )
    elements.append(
        f'<text transform="translate(18 {height/2}) rotate(-90)" '
        'text-anchor="middle" font-family="sans-serif" font-size="13">'
        'log2(pseudobulk CPM + 1)</text>'
    )
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    matrix_path = DATA / f"{GSE}_series_matrix.txt.gz"
    raw_path = DATA / f"{GSE}_RAW.tar"
    download(MATRIX_URL, matrix_path)
    download(RAW_URL, raw_path)

    metadata = parse_metadata(matrix_path)
    patient_rows = analyze_raw(raw_path, metadata)
    patient_fields = [
        "gsm", "sample", "response", "n_cells", "total_counts",
        "TACSTD2_counts", "TACSTD2_log2_cpm1", "TACSTD2_detected_cells",
        "TACSTD2_detected_fraction", "CLDN4_counts", "CLDN4_log2_cpm1",
        "CLDN4_detected_cells", "CLDN4_detected_fraction",
    ]
    write_tsv(OUT / "GSE233203_pseudobulk.tsv", patient_rows, patient_fields)

    tests: list[dict[str, object]] = []
    p_values: list[float] = []
    for gene in GENES:
        response = [
            float(row[f"{gene}_log2_cpm1"])
            for row in patient_rows
            if row["response"] == "Response"
        ]
        nonresponse = [
            float(row[f"{gene}_log2_cpm1"])
            for row in patient_rows
            if row["response"] == "Non-response"
        ]
        u_value, p_value = exact_mann_whitney(response, nonresponse)
        p_values.append(p_value)
        tests.append(
            {
                "dataset": GSE,
                "gene": gene,
                "comparison": "Response vs Non-response",
                "unit": "patient-level whole-sample pseudobulk",
                "n_response": len(response),
                "n_nonresponse": len(nonresponse),
                "median_response": statistics.median(response),
                "median_nonresponse": statistics.median(nonresponse),
                "mann_whitney_u": u_value,
                "exact_p": p_value,
                "holm_p_two_genes": "",
                "cliffs_delta": cliffs_delta(response, nonresponse),
            }
        )
        write_svg(OUT / f"GSE233203_{gene}_response.svg", gene, patient_rows)
    for row, adjusted in zip(tests, holm_adjust(p_values), strict=True):
        row["holm_p_two_genes"] = adjusted
    write_tsv(OUT / "marker_outcome_tests.tsv", tests, list(tests[0]))

    manifest = [
        {
            "accession": GSE,
            "file": matrix_path.name,
            "url": MATRIX_URL,
            "bytes": matrix_path.stat().st_size,
            "sha256": sha256(matrix_path),
        },
        {
            "accession": GSE,
            "file": raw_path.name,
            "url": RAW_URL,
            "bytes": raw_path.stat().st_size,
            "sha256": sha256(raw_path),
        },
    ]
    write_tsv(OUT / "input_manifest.tsv", manifest, list(manifest[0]))
    print(f"Wrote {len(patient_rows)} patients and {len(tests)} tests to {OUT}")


if __name__ == "__main__":
    main()
