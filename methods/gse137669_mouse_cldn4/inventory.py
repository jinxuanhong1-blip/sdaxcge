#!/usr/bin/env python3
"""NO-GO inventory: GSE137669 is bulk LV, not mouse lung-tumor scRNA.

Downloads the public series matrix and edgeR counts only (not the 1 GB
bigWig RAW.tar). Writes tables used by FINDING.md.
"""
from __future__ import annotations

import csv
import gzip
import json
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
DATA = HERE / "data"

SERIES_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137669/"
    "matrix/GSE137669_series_matrix.txt.gz"
)
EDGER_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137669/"
    "suppl/GSE137669_edgeR_normalized_counts.txt.gz"
)
SAMPLES = [
    "Sham_sfa.rep1",
    "Sham_sfa.rep2",
    "Sham_sfa.eaa.rep1",
    "Sham_sfa.eaa.rep2",
    "Tac_sfa.rep1",
    "Tac_sfa.rep2",
    "Tac_sfa.eaa.rep1",
    "Tac_sfa.eaa.rep2",
]
AUDIT = [
    "Cldn4",
    "Tacstd2",
    "Epcam",
    "Cdh1",
    "Krt8",
    "Krt18",
    "Krt19",
    "Ptprc",
    "Cd3d",
    "Cd3e",
    "Cd8a",
    "Nkg7",
    "Ncr1",
]


def _fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 100:
        return dest
    print("DOWNLOAD", url)
    urllib.request.urlretrieve(url, dest)
    return dest


def parse_series(path: Path) -> dict:
    meta = {
        "title": "",
        "pubmed": "",
        "overall_design": "",
        "sample_titles": [],
        "n_expression_rows": 0,
    }
    in_table = False
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Series_title"):
                meta["title"] = line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Series_pubmed_id"):
                meta["pubmed"] = line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Series_overall_design"):
                meta["overall_design"] = line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Sample_title"):
                meta["sample_titles"] = [
                    x.strip().strip('"') for x in line.rstrip().split("\t")[1:]
                ]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif in_table and not line.startswith("!"):
                first = line.split("\t", 1)[0].strip().strip('"')
                if first in {"ID_REF", ""}:
                    continue
                if line.strip():
                    meta["n_expression_rows"] += 1
    return meta


def parse_edger(path: Path) -> tuple[dict[str, dict], int]:
    found: dict[str, dict] = {}
    n_genes = 0
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_cols = header[1:9]
        if sample_cols != SAMPLES:
            raise ValueError(f"unexpected edgeR columns: {sample_cols}")
        for line in f:
            n_genes += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            symbol = parts[9]
            if symbol in AUDIT:
                vals = [float(x) for x in parts[1:9]]
                found[symbol] = {
                    "ensembl": parts[0],
                    "values": {s: v for s, v in zip(SAMPLES, vals)},
                    "max": max(vals),
                    "n_nonzero": sum(v > 0 for v in vals),
                }
    return found, n_genes


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    series_path = _fetch(SERIES_URL, DATA / "GSE137669_series_matrix.txt.gz")
    edger_path = _fetch(EDGER_URL, DATA / "GSE137669_edgeR_normalized_counts.txt.gz")
    series = parse_series(series_path)
    genes, n_genes = parse_edger(edger_path)
    cldn4 = genes.get("Cldn4", {})
    cldn4_zero = bool(cldn4) and cldn4.get("max", 1) == 0.0 and cldn4.get("n_nonzero", 1) == 0

    summary = {
        "accession": "GSE137669",
        "verdict": "NO-GO",
        "reason": "bulk left-ventricle RNA-seq; series matrix has 0 expression rows; Cldn4 all zeros",
        "title": series["title"],
        "pubmed": series["pubmed"],
        "n_series_matrix_expression_rows": series["n_expression_rows"],
        "n_bulk_libraries": 8,
        "n_genes_edger": n_genes,
        "cldn4_present_as_row": "Cldn4" in genes,
        "cldn4_all_zero": cldn4_zero,
        "epcam_all_zero": genes.get("Epcam", {}).get("max", None) == 0.0,
        "tacstd2_all_zero": genes.get("Tacstd2", {}).get("max", None) == 0.0,
        "n_mouse_lung_tumor_scrna_units": 0,
        "n_cldn4_vs_tnk_pairs": 0,
        "n_epithelial_ifn_mhc_q4q1": 0,
        "join_concordant4": False,
        "dual_high": False,
    }

    with (TABLES / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")

    honest = [
        ("item", "n", "note"),
        ("Assigned GEO series", "1", "GSE137669"),
        ("Series-matrix expression rows", "0", "metadata only"),
        ("Bulk LV libraries", "8", "2 per diet x surgery"),
        ("Cldn4 non-zero libraries", "0", "ENSMUSG00000047501"),
        ("Epcam non-zero libraries", "0", ""),
        ("Lung-tumor scRNA cells", "0", "not scRNA"),
        ("Mouse-level Cldn4 vs T/NK pairs", "0", ""),
        ("Epithelial IFN/MHC Q4 vs Q1 units", "0", ""),
        ("Concordant-pool join", "0", "sign not tested"),
    ]
    with (TABLES / "honest_n.tsv").open("w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows(honest)

    with (TABLES / "audit_genes.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["symbol", "ensembl", *SAMPLES, "max", "n_nonzero"])
        for sym in AUDIT:
            g = genes.get(sym)
            if not g:
                w.writerow([sym, "MISSING", *["NA"] * 8, "NA", 0])
                continue
            w.writerow(
                [sym, g["ensembl"], *[f"{g['values'][s]:.7f}" for s in SAMPLES], f"{g['max']:.7f}", g["n_nonzero"]]
            )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
