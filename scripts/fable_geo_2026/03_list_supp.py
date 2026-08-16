#!/usr/bin/env python3
"""List supplementary files + sizes for a shortlist of GSEs via GEO FTP dir
listing (https). Records everything so nothing is invented and the <2GB /
open-processed constraints can be checked.
Writes results/fable_geo_2026/supp_listing.tsv.
"""
import json
import time
import urllib.request
from pathlib import Path

RES = Path(__file__).resolve().parents[2] / "results" / "fable_geo_2026"

SHORTLIST = [
    "GSE253564", "GSE309652", "GSE248249", "GSE285029", "GSE292421",
    "GSE266219", "GSE265899", "GSE280232", "GSE241934", "GSE243013",
    "GSE271689", "GSE292098", "GSE261345", "GSE261348", "GSE311200",
    "GSE300685", "GSE270711", "GSE285298", "GSE237087",
]


def ftp_dir(gse: str):
    stub = gse[:-3] + "nnn"
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/")
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-fetch"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace"), url
        except Exception as e:  # noqa: BLE001
            if attempt == 4:
                return f"ERROR {e}", url
            time.sleep(2 ** attempt)
    return "ERROR", url


def parse_apache(html: str):
    """Parse NCBI autoindex (JSON-ish or HTML). NCBI serves an HTML table."""
    import re
    files = []
    # NCBI lists rows like: <a href="file">file</a> ... size
    for m in re.finditer(r'<a href="([^"/][^"]*)">[^<]+</a>\s*([0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:]+)\s+([0-9.]+[KMGT]?|-)', html):
        name, _date, size = m.group(1), m.group(2), m.group(3)
        files.append((name, size))
    if not files:  # fallback: just hrefs
        for m in re.finditer(r'<a href="([^"/?][^"]*)">', html):
            if m.group(1) not in ("Parent Directory",):
                files.append((m.group(1), "?"))
    return files


def main() -> None:
    rows = []
    for gse in SHORTLIST:
        html, url = ftp_dir(gse)
        if html.startswith("ERROR"):
            rows.append((gse, html, "", url))
            print(gse, html)
            continue
        files = parse_apache(html)
        for name, size in files:
            rows.append((gse, name, size, url + name))
        print(f"{gse}: {len(files)} files")
        time.sleep(0.3)
    with (RES / "supp_listing.tsv").open("w") as fh:
        fh.write("gse\tfile\tsize\turl\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print("wrote", RES / "supp_listing.tsv")


if __name__ == "__main__":
    main()
