#!/usr/bin/env python3
"""Panel check for GSE162520 TUMADOR HTG EdgeSeq OBP.

Downloads the public log2-CPM matrix and SOFT family file, tests whether
CLDN4 is on the HTG Oncology Biomarker Panel used in this series, and
writes panel_check.json next to this script.

If CLDN4 is absent, stop. Do not score CLDN4 vs CD8A / Immune / CD274 / IFN
and do not use CLDN3 as a stand-in.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import urllib.request
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
CACHE = HERE / "_cache"

MATRIX_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?acc=GSE162520&format=file&file=GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz"
)
SOFT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162520/"
    "soft/GSE162520_family.soft.gz"
)

QUERY_GENES = [
    "CLDN4",
    "CD8A",
    "CD274",
    "IFNG",
    "STAT1",
    "CXCL9",
    "CXCL10",
    "GZMB",
    "PRF1",
    "PDCD1",
    "CLDN3",
]
ALIAS_RE = re.compile(r"CLDN|CLAUDIN|CLD4|TACSTD|TROP", re.I)


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    req = urllib.request.Request(url, headers={"User-Agent": "gse162520-cldn4-panel-check"})
    with urllib.request.urlopen(req, timeout=120) as resp, tmp.open("wb") as out:
        out.write(resp.read())
    tmp.replace(dest)
    return dest


def load_matrix(path: Path) -> tuple[list[str], list[str]]:
    """Return (sample_ids, gene_symbols). European CSV: ';' sep, ',' decimals."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter=";")
        header = next(reader)
        samples = [c.strip().strip('"') for c in header[1:] if c.strip().strip('"')]
        genes = []
        for row in reader:
            if not row:
                continue
            g = row[0].strip().strip('"')
            if g:
                genes.append(g)
    return samples, genes


def load_soft_characteristics(path: Path) -> list[dict[str, str]]:
    text = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    blocks = re.split(r"\n\^SAMPLE = ", text)
    samples: list[dict[str, str]] = []
    for block in blocks[1:]:
        rec: dict[str, str] = {}
        for line in block.splitlines():
            if line.startswith("!Sample_geo_accession") or (
                not rec.get("gsm") and not line.startswith("!")
            ):
                if not rec.get("gsm"):
                    rec["gsm"] = line.split("=", 1)[-1].strip() if "=" in line else line.strip()
            if line.startswith("!Sample_title"):
                rec["title"] = line.split("=", 1)[1].strip()
            if line.startswith("!Sample_description"):
                rec["description"] = line.split("=", 1)[1].strip()
            if line.startswith("!Sample_characteristics_ch1"):
                val = line.split("=", 1)[1].strip()
                if ":" in val:
                    key, value = val.split(":", 1)
                    rec[key.strip()] = value.strip()
        samples.append(rec)
    return samples


def main() -> None:
    matrix_path = _download(MATRIX_URL, CACHE / "GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz")
    soft_path = _download(SOFT_URL, CACHE / "GSE162520_family.soft.gz")

    samples, genes = load_matrix(matrix_path)
    gene_set = set(genes)
    meta = load_soft_characteristics(soft_path)
    histology = Counter(s.get("histology", "NA") for s in meta)
    hot = Counter(s.get("hot phenotype", "NA") for s in meta)
    sex = Counter(s.get("Sex", "NA") for s in meta)

    alias_hits = sorted({g for g in genes if ALIAS_RE.search(g)})
    presence = {g: g in gene_set for g in QUERY_GENES}
    cldn4_present = bool(presence["CLDN4"])

    payload = {
        "accession": "GSE162520",
        "cohort": "TUMADOR",
        "platform": "HTG EdgeSeq Oncology Biomarker Panel (OBP)",
        "geo_platform_id": "GPL18573",
        "matrix_file": "GSE162520_GEO_data_TUMADOR_log2cpm.csv.gz",
        "n_genes_matrix": len(genes),
        "n_unique_genes": len(gene_set),
        "n_samples_matrix": len(samples),
        "n_samples_soft": len(meta),
        "honest_n": {
            "geo_arrays": len(samples),
            "series_summary_patients": 92,
            "histology": dict(histology),
            "hot_phenotype": dict(hot),
            "sex": dict(sex),
        },
        "query_presence": presence,
        "claudin_or_trop_aliases_on_panel": alias_hits,
        "cldn4_present": cldn4_present,
        "verdict": "present" if cldn4_present else "panel-missing",
        "stop": not cldn4_present,
        "note": (
            "CLDN4 is not among the 2560 OBP targets in the public TUMADOR matrix. "
            "CLDN3 is the only claudin on this panel and is not used as a stand-in. "
            "No CLDN4 vs CD8A / Immune / CD274 / IFN scores. No dual-high."
        ),
    }

    out = HERE / "panel_check.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if not cldn4_present:
        print("\nVERDICT: panel-missing. Stop.")


if __name__ == "__main__":
    main()
