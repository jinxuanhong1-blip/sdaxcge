#!/usr/bin/env python3
"""Scan ALL candidate GSEs' series-matrix sample characteristics for ICI
response / outcome labels. Downloads only the small metadata matrices.
Writes results/fable_geo_2026/label_scan_all.tsv with, per GSE, whether a
response-like field exists and which characteristic keys were seen.
"""
import csv
import gzip
import io
import time
import urllib.request
from pathlib import Path

RES = Path(__file__).resolve().parents[2] / "results" / "fable_geo_2026"
META = RES / "series_matrix"
META.mkdir(exist_ok=True)

RESP_KEYS = ("respon", "recist", "pfs", "overall surviv", "os time", "os event",
             "dcb", "benefit", "pcr", " mpr", "pathologic", "progression",
             "sensitiv", "resistan", "outcome", "efficacy", "best response",
             "irrecist", "clinical benefit")


def fetch(gse: str):
    stub = gse[:-3] + "nnn"
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo"})
            with urllib.request.urlopen(req, timeout=90) as r:
                raw = r.read()
            return gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode(
                "utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None  # multi-platform superseries: no single matrix
            time.sleep(2 ** attempt)
        except Exception:  # noqa: BLE001
            time.sleep(2 ** attempt)
    return None


def char_keys(txt: str):
    keys = []
    for ln in txt.splitlines():
        if ln.startswith("!Sample_characteristics"):
            # take the key part before first ':' of the first quoted value
            parts = ln.split("\t")[1:]
            for p in parts:
                p = p.strip().strip('"')
                if ":" in p:
                    keys.append(p.split(":", 1)[0].strip().lower())
                    break
    return keys


def main() -> None:
    cands = [r["accession"] for r in csv.DictReader(
        (RES / "geo_search_candidates.tsv").open(), delimiter="\t")]
    rows = []
    for i, gse in enumerate(cands):
        txt = fetch(gse)
        if txt is None:
            rows.append((gse, "NO_SINGLE_MATRIX", "", ""))
            continue
        keys = char_keys(txt)
        low = "\n".join(keys)
        # also scan whole char block text for values like response fields
        char_block = "\n".join(
            ln.lower() for ln in txt.splitlines()
            if ln.startswith("!Sample_characteristics")
            or ln.startswith("!Sample_title")
            or ln.startswith("!Sample_description")
            or ln.startswith("!Sample_source_name"))
        hits = sorted({k.strip() for k in RESP_KEYS if k in char_block})
        n = txt.count('"GSM')
        rows.append((gse, "OK", ";".join(hits), ";".join(sorted(set(keys)))))
        if hits:
            print(f"{gse}: RESP {hits}")
        if (i + 1) % 25 == 0:
            print(f"  ...scanned {i+1}/{len(cands)}")
        time.sleep(0.15)
    with (RES / "label_scan_all.tsv").open("w") as fh:
        fh.write("gse\tstatus\tresp_hits\tchar_keys\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    got = [r for r in rows if r[2]]
    print(f"\n{len(got)} series with response-like labels:")
    for r in got:
        print(" ", r[0], r[2])


if __name__ == "__main__":
    main()
