#!/usr/bin/env python3
"""Print the full GEO series-matrix header (title, design, summary, sample titles, characteristics).

Used for manual accession verification: `python3 inspect_series.py GSE161537 GSE162520 ...`
Re-uses the download cache created by 03_probe_metadata.py.
"""
from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

spec = importlib.util.spec_from_file_location("p3", Path(__file__).with_name("03_probe_metadata.py"))
p3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p3)


def header_of(acc: str) -> dict[str, list[str]]:
    base = f"{p3.FTP_BASE}/{p3.stub(acc)}/{acc}/matrix/"
    files = [f for f, _ in p3.list_dir(base) if f.endswith("series_matrix.txt.gz")]
    header: dict[str, list[str]] = {}
    for fname in files[:3]:
        cached = p3.CACHE / fname
        blob = cached.read_bytes() if cached.exists() else p3.fetch(base + fname, binary=True)
        if blob is None:
            continue
        if not cached.exists():
            cached.write_bytes(blob)
        for line in gzip.decompress(blob).decode("utf-8", "replace").splitlines():
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!"):
                key, _, val = line[1:].partition("\t")
                header.setdefault(key.strip(), []).append(val.strip())
    return header


def main(argv: list[str]) -> int:
    for acc in argv:
        h = header_of(acc)
        print("=" * 110)
        print(acc, "|", " ".join(h.get("Series_title", []))[:200])
        for key in (
            "Series_summary",
            "Series_overall_design",
            "Series_type",
            "Series_platform_id",
            "Series_pubmed_id",
            "Series_relation",
        ):
            if h.get(key):
                print(f"  {key}: {' '.join(h[key])[:900]}")
        for key in ("Sample_geo_accession", "Sample_title", "Sample_source_name_ch1"):
            if h.get(key):
                vals = h[key][0].replace('"', "").split("\t")
                print(f"  {key} (n={len(vals)}): {vals[:40]}")
        for line in h.get("Sample_characteristics_ch1", []):
            vals = [v.strip('"') for v in line.split("\t")]
            uniq = sorted(set(vals))
            print(f"  CHAR n={len(vals)} uniq={len(uniq)}: {uniq[:10]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
