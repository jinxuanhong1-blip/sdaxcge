#!/usr/bin/env python3
"""Fetch per-sample metadata (title, GSM, characteristics, description) from GEO
series-matrix headers so sample columns can be mapped to treatment groups.

Writes notes/fable_mouse_ici/sample_metadata.json and prints a summary.
"""
import gzip
import io
import json
import time
import urllib.request
from pathlib import Path

GEO_ACCS = [
    "GSE239485", "GSE222158", "GSE297630", "GSE330658", "GSE197260",
    "GSE133604", "GSE129297", "GSE297632",
]
NOTES = Path(__file__).resolve().parents[2] / "notes" / "fable_mouse_ici"
UA = {"User-Agent": "Mozilla/5.0 (fable-mouse-ici meta)"}

KEEP = (
    "!Sample_title", "!Sample_geo_accession", "!Sample_source_name_ch1",
    "!Sample_description",
)


def fetch(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** i)
    raise last


def get_matrix(acc):
    nnn = acc[:-3] + "nnn"
    url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{acc}/matrix/"
        f"{acc}_series_matrix.txt.gz"
    )
    raw = fetch(url)
    return gzip.decompress(raw).decode("utf-8", "replace")


def parse(acc, txt):
    rec = {"accession": acc, "fields": {}}
    for line in txt.splitlines():
        if line.startswith("!Sample_"):
            key = line.split("\t", 1)[0]
            vals = [v.strip('"') for v in line.split("\t")[1:]]
            if key in KEEP or key.startswith("!Sample_characteristics"):
                rec["fields"].setdefault(key, []).append(vals)
        if line.startswith("!series_matrix_table_begin"):
            break
    return rec


def main():
    out = []
    for acc in GEO_ACCS:
        print("[meta]", acc)
        try:
            txt = get_matrix(acc)
            out.append(parse(acc, txt))
        except Exception as e:  # noqa: BLE001
            out.append({"accession": acc, "error": str(e)})
    (NOTES / "sample_metadata.json").write_text(json.dumps(out, indent=2))

    for rec in out:
        print("=" * 70)
        print(rec["accession"], rec.get("error", ""))
        f = rec.get("fields", {})
        titles = f.get("!Sample_title", [[]])
        gsms = f.get("!Sample_geo_accession", [[]])
        titles = titles[0] if titles else []
        gsms = gsms[0] if gsms else []
        for i, t in enumerate(titles):
            g = gsms[i] if i < len(gsms) else ""
            print(f"   {g}  {t}")
        for k, v in f.items():
            if k.startswith("!Sample_characteristics"):
                print("   CHAR", k, "->", v[0][:8])
            if k == "!Sample_source_name_ch1":
                print("   SRC ->", v[0][:8])


if __name__ == "__main__":
    main()
