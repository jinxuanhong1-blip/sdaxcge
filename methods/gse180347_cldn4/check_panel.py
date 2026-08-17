#!/usr/bin/env python3
"""Verify whether CLDN4 or TACSTD2 is on GSE180347 (GPL29738) and report honest n.

Stop condition: if both genes are absent from the NanoString panel / series matrix,
do not score vs hot/cold or PD-L1.
"""
from __future__ import annotations

import gzip
import json
import urllib.request
from collections import Counter
from pathlib import Path

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180347/"
    "matrix/GSE180347_series_matrix.txt.gz"
)
PLATFORM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL29nnn/GPL29738/"
    "soft/GPL29738_family.soft.gz"
)
TARGETS = ("CLDN4", "TACSTD2")
ALIASES = {
    "CLDN4": ("CLDN4", "NM_001305"),
    "TACSTD2": ("TACSTD2", "TROP2", "TROP-2", "EGP1", "EGP-1", "GA733-1", "M1S1", "NM_002353"),
}


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    urllib.request.urlretrieve(url, dest)
    return dest


def parse_platform(path: Path) -> list[dict]:
    rows: list[dict] = []
    in_table = False
    header: list[str] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = next(fh).rstrip("\n").split("\t")
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                parts = line.rstrip("\n").split("\t")
                rec = {header[i]: (parts[i] if i < len(parts) else "") for i in range(len(header))}
                rows.append(rec)
    return rows


def parse_matrix(path: Path) -> dict:
    titles: list[str] = []
    groups: list[str] = []
    design = ""
    summary = ""
    expr_genes: list[str] = []
    in_table = False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Series_overall_design"):
                design = line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Series_summary"):
                summary = line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and "group:" in line:
                groups = [
                    x.strip('"').replace("group: ", "")
                    for x in line.rstrip("\n").split("\t")[1:]
                ]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                gene_id = line.split("\t", 1)[0].strip().strip('"')
                if gene_id == "ID_REF" or not gene_id:
                    continue
                expr_genes.append(gene_id)
    return {
        "titles": titles,
        "groups": groups,
        "design": design,
        "summary": summary,
        "expr_genes": expr_genes,
    }


def search_rows(rows: list[dict], needles: tuple[str, ...]) -> list[str]:
    hits = []
    for rec in rows:
        blob = " ".join(rec.values()).upper()
        if any(n.upper() in blob for n in needles):
            hits.append(rec.get("ID", ""))
    return sorted(set(hits))


def main() -> None:
    out_dir = Path(__file__).resolve().parent
    cache = out_dir / "cache"
    plat_path = _download(PLATFORM_URL, cache / "GPL29738_family.soft.gz")
    mat_path = _download(MATRIX_URL, cache / "GSE180347_series_matrix.txt.gz")

    plat = parse_platform(plat_path)
    mat = parse_matrix(mat_path)
    codeclass = Counter(r.get("CodeClass", "") for r in plat)
    group_n = dict(Counter(mat["groups"]))
    title_n = dict(Counter(t.rsplit("_", 1)[-1] for t in mat["titles"]))

    target_hits = {}
    for gene, needles in ALIASES.items():
        plat_hits = search_rows(plat, needles)
        expr_hits = [g for g in mat["expr_genes"] if any(n.upper() in g.upper() for n in needles)]
        target_hits[gene] = {
            "on_platform": bool(plat_hits),
            "in_series_matrix": bool(expr_hits),
            "platform_ids": plat_hits,
            "matrix_ids": expr_hits,
        }

    present = [g for g, h in target_hits.items() if h["on_platform"] or h["in_series_matrix"]]
    verdict = "score" if present else "panel_missing_stop"

    report = {
        "accession": "GSE180347",
        "platform": "GPL29738",
        "panel": "NanoString nCounter PanCancer Immune Profiling Panel",
        "verdict": verdict,
        "targets_present": present,
        "targets": target_hits,
        "n_geo_arrays": len(mat["titles"]),
        "n_expr_genes": len(mat["expr_genes"]),
        "platform_codeclass": dict(codeclass),
        "group_n_from_geo_labels": group_n,
        "group_n_from_sample_titles": title_n,
        "design_text_n": {
            "PH": 23,
            "PC": 7,
            "NH": 53,
            "NC": 55,
            "sum": 138,
            "source": "GEO !Series_overall_design",
        },
        "paper_abstract_n": 138,
        "nearby_on_panel": {
            "EPCAM": "EPCAM" in mat["expr_genes"],
            "CDH1": "CDH1" in mat["expr_genes"],
            "CD274": "CD274" in mat["expr_genes"],
        },
        "note": (
            "NM_002354 on this platform is EPCAM (TACSTD1), not TACSTD2 (NM_002353). "
            "Do not treat EPCAM as a TACSTD2 surrogate."
        ),
    }

    (out_dir / "panel_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if verdict != "score":
        print("\nSTOP: CLDN4 and TACSTD2 are absent from the panel. No phenotype scores.")


if __name__ == "__main__":
    main()
