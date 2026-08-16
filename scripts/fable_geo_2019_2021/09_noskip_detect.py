#!/usr/bin/env python3
"""Detect TACSTD2/CLDN4 (and Ensembl IDs) in leftover downloaded files.

Also extract !Sample_characteristics keys from series_matrix files so we can
honestly say which leftover series have a deposited ICI outcome.
"""
import csv
import gzip
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "noskip" / "GEO_2019_2021"
DL = OUT / "downloads"

GENE_RE = re.compile(
    rb"(TACSTD2|CLDN4|ENSG00000184292|ENSG00000189143)",
    re.I,
)


def open_any(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def detect_genes(path, max_hits=8):
    hits = []
    try:
        with open_any(path) as f:
            for i, line in enumerate(f):
                if GENE_RE.search(line):
                    snippet = line[:180].decode("utf-8", "replace").replace("\n", " ")
                    hits.append({"line": i + 1, "snippet": snippet})
                    if len(hits) >= max_hits:
                        break
    except Exception as e:  # noqa: BLE001
        return {"error": str(e), "hits": []}
    return {"error": None, "hits": hits}


def parse_matrix_keys(path):
    keys = {"n_samples": None, "char_keys": [], "source_names": [], "titles": []}
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("!Sample_geo_accession"):
                    keys["n_samples"] = len(line.rstrip().split("\t")) - 1
                elif line.startswith("!Sample_title") and not keys["titles"]:
                    keys["titles"] = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:6]]
                elif line.startswith("!Sample_source_name"):
                    vals = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
                    keys["source_names"] = sorted(set(vals))[:12]
                elif line.startswith("!Sample_characteristics"):
                    vals = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
                    field = None
                    for v in vals:
                        if ":" in v:
                            field = v.split(":", 1)[0].strip()
                            break
                    if field:
                        keys["char_keys"].append(field)
    except Exception as e:  # noqa: BLE001
        keys["error"] = str(e)
    return keys


def main():
    gene_rows = []
    matrix_rows = []
    for gse_dir in sorted(DL.iterdir()):
        if not gse_dir.is_dir():
            continue
        gse = gse_dir.name
        for path in sorted(gse_dir.iterdir()):
            if path.name.endswith("_series_matrix.txt.gz"):
                info = parse_matrix_keys(path)
                matrix_rows.append({"gse": gse, **{k: (";".join(v) if isinstance(v, list) else v)
                                                   for k, v in info.items()}})
                # series_matrix may also contain expression values
                det = detect_genes(path)
                gene_rows.append({
                    "gse": gse, "file": path.name, "kind": "series_matrix",
                    "has_TACSTD2_or_CLDN4": bool(det["hits"]),
                    "n_hits_shown": len(det["hits"]),
                    "error": det["error"] or "",
                    "first_hit": det["hits"][0]["snippet"] if det["hits"] else "",
                })
            else:
                det = detect_genes(path)
                gene_rows.append({
                    "gse": gse, "file": path.name, "kind": "suppl",
                    "has_TACSTD2_or_CLDN4": bool(det["hits"]),
                    "n_hits_shown": len(det["hits"]),
                    "error": det["error"] or "",
                    "first_hit": det["hits"][0]["snippet"] if det["hits"] else "",
                })
            print(f"{gse}/{path.name}: genes={bool(det['hits'])} n={len(det['hits'])} {det['error'] or ''}")

    with open(OUT / "gene_presence.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(gene_rows[0].keys()))
        w.writeheader()
        w.writerows(gene_rows)
    with open(OUT / "leftover_series_matrix_keys.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(matrix_rows[0].keys()) if matrix_rows else ["gse"])
        w.writeheader()
        w.writerows(matrix_rows)
    print("Wrote gene_presence.csv and leftover_series_matrix_keys.csv")


if __name__ == "__main__":
    main()
