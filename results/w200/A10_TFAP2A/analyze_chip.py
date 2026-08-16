#!/usr/bin/env python3
"""Public TFAP2A ChIP overlap at TACSTD2 and CLDN4.

No lung TFAP2A ChIP-seq was found. This script records that inventory and
checks the available non-lung public peak sets plus ChIP-Atlas target-gene
scores. Large assembled BED files are streamed and not stored in the repo.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import urllib.parse
import urllib.request
from pathlib import Path


OUT = Path(__file__).resolve().parent
CACHE = Path("/tmp/a10_chip")
CACHE.mkdir(parents=True, exist_ok=True)

GENES = {
    "TACSTD2": {
        "ensembl_id": "ENSG00000184292",
        "chrom": "chr1",
        "gene_start": 58575433,
        "gene_end": 58577252,
        "strand": "-",
        "canonical_tss": 58577252,
        "source": "Ensembl REST lookup, GRCh38, 2026-08-16",
    },
    "CLDN4": {
        "ensembl_id": "ENSG00000189143",
        "chrom": "chr7",
        "gene_start": 73799542,
        "gene_end": 73832690,
        "strand": "+",
        "canonical_tss": 73830996,
        "source": "Ensembl REST lookup, GRCh38, 2026-08-16; TSS from canonical ENST00000340958",
    },
}

REMAP_DATASETS = [
    {
        "source": "ReMap2022",
        "dataset_id": "GSE60270.TFAP2A.MCF-7",
        "cell": "MCF-7",
        "tissue_class": "Breast",
        "lung": "no",
        "url": "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/DATASET/GSE60270.TFAP2A.MCF-7.bed.gz",
    },
    {
        "source": "ReMap2022",
        "dataset_id": "GSE60270.TFAP2A.MCF-7_E2",
        "cell": "MCF-7_E2",
        "tissue_class": "Breast",
        "lung": "no",
        "url": "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/DATASET/GSE60270.TFAP2A.MCF-7_E2.bed.gz",
    },
    {
        "source": "ReMap2022",
        "dataset_id": "GSE105081.TFAP2A.WA09",
        "cell": "WA09",
        "tissue_class": "Pluripotent stem cell",
        "lung": "no",
        "url": "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/DATASET/GSE105081.TFAP2A.WA09.bed.gz",
    },
]

CHIPATLAS_QVALS = ("05", "10")
WINDOWS_KB = (1, 5, 10)


def get_json(url: str):
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def download(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as response, tmp.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    tmp.replace(dest)
    return dest


def write_tsv(name: str, fieldnames: list[str], rows: list[dict]) -> None:
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def windows_for(gene: dict) -> dict[str, tuple[int, int]]:
    tss = gene["canonical_tss"]
    out = {
        "gene_body": (gene["gene_start"], gene["gene_end"]),
    }
    for kb in WINDOWS_KB:
        out[f"tss_pm_{kb}kb"] = (tss - kb * 1000, tss + kb * 1000)
    return out


def overlaps(start: int, end: int, window: tuple[int, int]) -> bool:
    return start < window[1] and end > window[0]


def open_bed(path: Path):
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return path.open(encoding="utf-8")


def parse_chipatlas_name(field: str) -> dict[str, str]:
    decoded = urllib.parse.unquote(field)
    parts = {}
    for item in decoded.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parts[key] = value
    return {
        "experiment_id": parts.get("ID", ""),
        "cell": parts.get("Name", "").replace("TFAP2A (@ ", "").rstrip(")"),
        "cell_group": parts.get("Cell group", ""),
        "title": parts.get("Title", ""),
    }


def encode_chip_antigen_count(cl_class: str) -> int:
    query = urllib.parse.urlencode(
        {
            "genome": "hg38",
            "agClass": "TFs and others",
            "clClass": cl_class,
        }
    )
    data = get_json(f"https://chip-atlas.org/data/chip_antigen?{query}")
    hit = next((item for item in data if item.get("id") == "TFAP2A"), None)
    return int(hit["count"]) if hit else 0


def chipatlas_bed_url(qval: str) -> str:
    payload = json.dumps(
        {
            "condition": {
                "genome": "hg38",
                "agClass": "TFs and others",
                "agSubClass": "TFAP2A",
                "clClass": "All cell types",
                "qval": qval,
            }
        }
    ).encode()
    request = urllib.request.Request(
        "https://chip-atlas.org/download",
        data=payload,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)["url"]


def scan_bed(path: Path, source: str, dataset_id: str, extra: dict[str, str]) -> list[dict]:
    by_chrom: dict[str, list[tuple[str, int, int]]] = {}
    for symbol, gene in GENES.items():
        start, end = windows_for(gene)["tss_pm_10kb"]
        by_chrom.setdefault(gene["chrom"], []).append((symbol, start, end))

    hits: list[dict] = []
    with open_bed(path) as handle:
        for line in handle:
            if not line or line.startswith(("track", "browser", "#")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start_s, end_s = fields[0], fields[1], fields[2]
            if chrom not in by_chrom:
                continue
            start, end = int(start_s), int(end_s)
            for symbol, window_start, window_end in by_chrom[chrom]:
                if overlaps(start, end, (window_start, window_end)):
                    meta = {
                        "experiment_id": dataset_id,
                        "cell": extra.get("cell", ""),
                        "cell_group": extra.get("tissue_class", ""),
                        "title": extra.get("title", ""),
                    }
                    if source == "ChIP-Atlas" and len(fields) > 3:
                        meta.update(parse_chipatlas_name(fields[3]))
                    score = fields[4] if len(fields) > 4 else ""
                    gene = GENES[symbol]
                    row = {
                        "source": source,
                        "dataset_id": dataset_id if source != "ChIP-Atlas" else meta["experiment_id"],
                        "threshold": extra.get("threshold", "NA"),
                        "gene": symbol,
                        "chrom": chrom,
                        "peak_start": start,
                        "peak_end": end,
                        "score": score,
                        "cell": meta["cell"],
                        "cell_group": meta["cell_group"],
                        "lung": "yes" if "lung" in f"{meta['cell']} {meta['cell_group']}".lower() else "no",
                    }
                    for window_name, window in windows_for(gene).items():
                        row[window_name] = "yes" if overlaps(start, end, window) else "no"
                    hits.append(row)
    return hits


def target_gene_rows() -> list[dict]:
    rows: list[dict] = []
    for kb in WINDOWS_KB:
        url = f"https://chip-atlas.dbcls.jp/data/hg38/target/TFAP2A.{kb}.tsv"
        dest = CACHE / f"TFAP2A.{kb}.tsv"
        download(url, dest)
        with dest.open(encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            found: set[str] = set()
            n_experiments = len(header) - 2
            for line in handle:
                cols = line.rstrip("\n").split("\t")
                if cols[0] not in GENES:
                    continue
                found.add(cols[0])
                nonzero = []
                for name, value in zip(header[2:], cols[2:]):
                    try:
                        score = float(value)
                    except ValueError:
                        continue
                    if score > 0:
                        nonzero.append((name, score))
                nonzero.sort(key=lambda item: item[1], reverse=True)
                rows.append(
                    {
                        "source": "ChIP-Atlas target genes",
                        "distance_from_tss": f"{kb}kb",
                        "gene": cols[0],
                        "n_experiments": n_experiments,
                        "average_score": cols[1],
                        "n_nonzero_experiments": len(nonzero),
                        "top_experiments": ";".join(
                            f"{name}:{score:g}" for name, score in nonzero[:8]
                        ),
                        "lung_experiments_nonzero": 0,
                    }
                )
            for gene in GENES:
                if gene not in found:
                    rows.append(
                        {
                            "source": "ChIP-Atlas target genes",
                            "distance_from_tss": f"{kb}kb",
                            "gene": gene,
                            "n_experiments": n_experiments,
                            "average_score": "0",
                            "n_nonzero_experiments": 0,
                            "top_experiments": "",
                            "lung_experiments_nonzero": 0,
                        }
                    )
    rows.sort(key=lambda row: (row["gene"], row["distance_from_tss"]))
    return rows


def main() -> None:
    inventory = [
        {
            "source": "ENCODE",
            "dataset_id": "ENCSR000EVP",
            "cell": "HeLa-S3",
            "tissue_class": "Uterus/cervix",
            "lung": "no",
            "status": "revoked",
            "n_peaks": "",
            "note": "Only ENCODE TFAP2A TF ChIP-seq experiment; revoked, not used for overlap",
            "url": "https://www.encodeproject.org/experiments/ENCSR000EVP/",
        }
    ]

    hit_rows: list[dict] = []
    for dataset in REMAP_DATASETS:
        dest = CACHE / Path(dataset["url"]).name
        download(dataset["url"], dest)
        n_peaks = 0
        with gzip.open(dest, "rt", encoding="utf-8") as handle:
            n_peaks = sum(1 for line in handle if line and not line.startswith(("track", "#")))
        inventory.append(
            {
                "source": dataset["source"],
                "dataset_id": dataset["dataset_id"],
                "cell": dataset["cell"],
                "tissue_class": dataset["tissue_class"],
                "lung": dataset["lung"],
                "status": "released",
                "n_peaks": n_peaks,
                "note": "ReMap 2022 uniformly processed public ChIP-seq; no lung biotype",
                "url": dataset["url"],
            }
        )
        hit_rows.extend(
            scan_bed(
                dest,
                "ReMap2022",
                dataset["dataset_id"],
                {
                    "cell": dataset["cell"],
                    "tissue_class": dataset["tissue_class"],
                    "threshold": "NA",
                },
            )
        )

    classes = [
        "Lung",
        "Breast",
        "Digestive tract",
        "Epidermis",
        "Neural",
        "Pluripotent stem cell",
        "Uterus",
        "Muscle",
        "Others",
        "Placenta",
        "Kidney",
        "Liver",
        "Prostate",
        "Blood",
    ]
    atlas_counts = {cl: encode_chip_antigen_count(cl) for cl in classes}
    for qval in CHIPATLAS_QVALS:
        url = chipatlas_bed_url(qval)
        dest = CACHE / f"Oth.ALL.{qval}.TFAP2A.AllCell.bed"
        download(url, dest)
        inventory.append(
            {
                "source": "ChIP-Atlas",
                "dataset_id": f"Oth.ALL.{qval}.TFAP2A.AllCell",
                "cell": "All cell types",
                "tissue_class": "mixed; Lung=" + str(atlas_counts["Lung"]),
                "lung": "no",
                "status": "released",
                "n_peaks": "",
                "note": (
                    f"Assembled MACS2 peaks, q-value 1e-{int(qval)}. "
                    f"ChIP-Atlas TFAP2A experiment counts by class: "
                    + ", ".join(f"{k}={v}" for k, v in atlas_counts.items() if v or k == "Lung")
                ),
                "url": url,
            }
        )
        hit_rows.extend(
            scan_bed(
                dest,
                "ChIP-Atlas",
                f"q{qval}",
                {"cell": "", "tissue_class": f"q{qval}", "threshold": f"q{qval}"},
            )
        )

    window_rows = []
    for symbol, gene in GENES.items():
        for window_name, (start, end) in windows_for(gene).items():
            window_rows.append(
                {
                    "gene": symbol,
                    "ensembl_id": gene["ensembl_id"],
                    "chrom": gene["chrom"],
                    "strand": gene["strand"],
                    "canonical_tss": gene["canonical_tss"],
                    "window": window_name,
                    "window_start": start,
                    "window_end": end,
                    "coordinate_note": gene["source"],
                }
            )

    summary_rows = []
    keys = sorted(
        {
            (row["source"], row["dataset_id"], row["threshold"], row["gene"])
            for row in hit_rows
        }
    )
    window_names = ["gene_body", *[f"tss_pm_{kb}kb" for kb in WINDOWS_KB]]
    for source, dataset_id, threshold, gene in keys:
        subset = [
            row
            for row in hit_rows
            if row["source"] == source
            and row["dataset_id"] == dataset_id
            and row["threshold"] == threshold
            and row["gene"] == gene
        ]
        summary_rows.append(
            {
                "source": source,
                "dataset_id": dataset_id,
                "threshold": threshold,
                "gene": gene,
                "n_overlapping_peaks_in_tss_pm_10kb": len(subset),
                **{
                    f"n_{window}": sum(row[window] == "yes" for row in subset)
                    for window in window_names
                },
                "cells": ";".join(sorted({row["cell"] for row in subset if row["cell"]})),
                "lung_peak_n": sum(row["lung"] == "yes" for row in subset),
            }
        )

    write_tsv(
        "chip_inventory.tsv",
        [
            "source",
            "dataset_id",
            "cell",
            "tissue_class",
            "lung",
            "status",
            "n_peaks",
            "note",
            "url",
        ],
        inventory,
    )
    write_tsv(
        "chip_locus_windows.tsv",
        [
            "gene",
            "ensembl_id",
            "chrom",
            "strand",
            "canonical_tss",
            "window",
            "window_start",
            "window_end",
            "coordinate_note",
        ],
        window_rows,
    )
    write_tsv(
        "chip_peak_hits.tsv",
        [
            "source",
            "dataset_id",
            "threshold",
            "gene",
            "chrom",
            "peak_start",
            "peak_end",
            "score",
            "cell",
            "cell_group",
            "lung",
            "gene_body",
            "tss_pm_1kb",
            "tss_pm_5kb",
            "tss_pm_10kb",
        ],
        hit_rows,
    )
    write_tsv(
        "chip_hit_summary.tsv",
        [
            "source",
            "dataset_id",
            "threshold",
            "gene",
            "n_overlapping_peaks_in_tss_pm_10kb",
            "n_gene_body",
            "n_tss_pm_1kb",
            "n_tss_pm_5kb",
            "n_tss_pm_10kb",
            "cells",
            "lung_peak_n",
        ],
        summary_rows,
    )
    write_tsv(
        "chip_target_gene_scores.tsv",
        [
            "source",
            "distance_from_tss",
            "gene",
            "n_experiments",
            "average_score",
            "n_nonzero_experiments",
            "top_experiments",
            "lung_experiments_nonzero",
        ],
        target_gene_rows(),
    )


if __name__ == "__main__":
    main()
