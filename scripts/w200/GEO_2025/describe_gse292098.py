#!/usr/bin/env python3
"""Leftover descriptive pass for GSE292098 GeoMx WTA.

The Greek lung ICI spatial series measures TACSTD2 and CLDN4 but GEO does not
provide a per-ROI response or survival key. This script only reports
compartment-level leftover expression. It does not invent outcomes.
"""

from __future__ import annotations

import csv
import hashlib
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "w200" / "GEO_2025"
DATA = OUT / "data"
GENES = ("TACSTD2", "CLDN4")
NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
FILES = {
    "CK": "GSE292098_CK_0_filter_Greek.xlsx",
    "CD45": "GSE292098_CD45_0_filter_Greek.xlsx",
    "CD68": "GSE292098_CD68_0_filter_Greek.xlsx",
}


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


def cell_text(cell: ET.Element, shared: list[str]) -> str:
    value_el = cell.find("m:v", NS)
    if cell.attrib.get("t") == "s" and value_el is not None and value_el.text:
        return shared[int(value_el.text)]
    if cell.attrib.get("t") == "inlineStr":
        return "".join(node.text or "" for node in cell.iter("{%s}t" % NS["m"]))
    return value_el.text if value_el is not None and value_el.text else ""


def read_sheet(path: Path, sheet_name: str) -> list[list[str]]:
    with ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {item.attrib["Id"]: item.attrib["Target"] for item in rels}
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = [
                "".join(node.text or "" for node in item.iter("{%s}t" % NS["m"]))
                for item in root
            ]
        target = None
        for sheet in workbook.find("m:sheets", NS):
            if sheet.attrib["name"] == sheet_name:
                rel_id = sheet.attrib["{%s}id" % NS["r"]]
                target = targets[rel_id]
                break
        if target is None:
            raise ValueError(f"{path.name} is missing sheet {sheet_name}")
        sheet_path = "xl/" + target.lstrip("/") if not target.startswith("/xl/") else target.lstrip("/")
        root = ET.fromstring(archive.read(sheet_path))
        rows: list[list[str]] = []
        for row in root.findall(".//m:sheetData/m:row", NS):
            rows.append([cell_text(cell, shared) for cell in row.findall("m:c", NS)])
        return rows


def median(values: list[float]) -> float:
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    long_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    manifest: list[dict[str, object]] = []
    for compartment, filename in FILES.items():
        url = (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE292nnn/GSE292098/suppl/"
            + filename
        )
        path = DATA / filename
        download(url, path)
        manifest.append(
            {
                "accession": "GSE292098",
                "file": filename,
                "url": url,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
        matrix = read_sheet(path, "TargetCountMatrix")
        header = matrix[0]
        gene_index = {row[0].upper(): row for row in matrix[1:]}
        for gene in GENES:
            if gene not in gene_index:
                raise ValueError(f"{filename} is missing {gene}")
            values = [float(item) for item in gene_index[gene][1:] if item]
            for segment, raw in zip(header[1:], gene_index[gene][1:]):
                if not raw:
                    continue
                long_rows.append(
                    {
                        "dataset": "GSE292098",
                        "compartment": compartment,
                        "segment": segment,
                        "gene": gene,
                        "raw_count": float(raw),
                    }
                )
            summary_rows.append(
                {
                    "dataset": "GSE292098",
                    "compartment": compartment,
                    "gene": gene,
                    "n_segments": len(values),
                    "median_raw_count": median(values),
                    "min_raw_count": min(values),
                    "max_raw_count": max(values),
                    "outcome_test": "not_performed",
                    "honest_reason": "GEO has no per-ROI ICI outcome key",
                }
            )
    with (OUT / "GSE292098_compartment_long.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(long_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(long_rows)
    with (OUT / "GSE292098_compartment_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    existing = []
    manifest_path = OUT / "input_manifest.tsv"
    if manifest_path.exists():
        with manifest_path.open(encoding="utf-8") as handle:
            existing = list(csv.DictReader(handle, delimiter="\t"))
    existing = [row for row in existing if row.get("accession") != "GSE292098"]
    existing.extend(manifest)
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(existing[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(existing)
    print(f"Wrote {len(long_rows)} leftover GeoMx rows without outcome tests")


if __name__ == "__main__":
    main()
