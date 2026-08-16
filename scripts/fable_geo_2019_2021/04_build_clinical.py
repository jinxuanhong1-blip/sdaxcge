#!/usr/bin/env python3
"""Parse GEO series_matrix files into per-sample clinical tables.

Outputs one CSV per GSE under results/fable_geo_2019_2021/clinical/ plus a
combined verification summary. All values are extracted verbatim from the
official GEO series_matrix files (no fabricated fields).
"""
import gzip
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DL = ROOT / "results" / "fable_geo_2019_2021" / "downloads"
OUT = ROOT / "results" / "fable_geo_2019_2021" / "clinical"
OUT.mkdir(parents=True, exist_ok=True)

GSES = ["GSE126044", "GSE135222", "GSE136961", "GSE111414", "GSE182328"]


def parse_series_matrix(path):
    """Return dict of row-label -> list[str] for !Sample_* rows."""
    rows = {}
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0]
            vals = [p.strip().strip('"') for p in parts[1:]]
            rows.setdefault(key, []).append(vals)
    return rows


def build_table(gse):
    sm = DL / gse / f"{gse}_series_matrix.txt.gz"
    rows = parse_series_matrix(sm)
    titles = rows["!Sample_title"][0]
    gsms = rows["!Sample_geo_accession"][0]
    n = len(gsms)
    table = [{"gsm": gsms[i], "title": titles[i]} for i in range(n)]

    # source name
    if "!Sample_source_name_ch1" in rows:
        src = rows["!Sample_source_name_ch1"][0]
        for i in range(n):
            table[i]["source_name"] = src[i]

    # characteristics: each occurrence is a "key: value" list
    for occ in rows.get("!Sample_characteristics_ch1", []):
        # derive field name from first non-empty entry
        field = None
        for v in occ:
            if ":" in v:
                field = v.split(":", 1)[0].strip()
                break
        if field is None:
            continue
        col = re.sub(r"[^0-9a-zA-Z]+", "_", field).strip("_").lower()
        for i in range(n):
            v = occ[i]
            table[i][col] = v.split(":", 1)[1].strip() if ":" in v else v
    return table


def main():
    summary = []
    for gse in GSES:
        table = build_table(gse)
        cols = []
        for r in table:
            for k in r:
                if k not in cols:
                    cols.append(k)
        out = OUT / f"{gse}_clinical.csv"
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(table)
        print(f"{gse}: {len(table)} samples -> {out.name} (fields: {cols})")
        summary.append((gse, len(table), cols))
    print("\nDone. Clinical tables in", OUT)


if __name__ == "__main__":
    main()
