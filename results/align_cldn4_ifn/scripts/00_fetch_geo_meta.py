#!/usr/bin/env python3
"""Fetch GEO series + sample metadata and supplementary file listings.

Writes one SOFT-brief text file per series plus a parsed sample table so the
experimental design of every series can be checked before any analysis.
"""
import os
import re
import sys
import time
import urllib.request

OUT = os.path.join(os.path.dirname(__file__), "..", "meta")
OUT = os.path.abspath(OUT)

SERIES = [
    "GSE207704", "GSE22493", "GSE50927", "GSE274940",   # CLDN4 sets (claimed)
    "GSE334497", "GSE289287", "GSE245459",              # TACSTD2 sets (claimed)
]

ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ={targ}&form=text&view={view}"


def get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-meta/1.0"})
            with urllib.request.urlopen(req, timeout=120) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as exc:  # transient NCBI throttling is common
            if i == tries - 1:
                raise
            sys.stderr.write(f"retry {i+1} {url}: {exc}\n")
            time.sleep(4 * 2 ** i)


def ftp_dir(gse):
    stub = gse[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/"
    try:
        html = get(url, tries=2)
    except Exception:
        return url, []
    files = sorted(set(re.findall(r'href="([^"?/][^"]*)"', html)))
    return url, [f for f in files if not f.startswith("http")]


def main():
    os.makedirs(OUT, exist_ok=True)
    for gse in SERIES:
        soft = get(ACC.format(acc=gse, targ="gse", view="brief"))
        with open(os.path.join(OUT, f"{gse}.series.txt"), "w") as fh:
            fh.write(soft)

        allmeta = get(ACC.format(acc=gse, targ="all", view="brief"))
        with open(os.path.join(OUT, f"{gse}.all.txt"), "w") as fh:
            fh.write(allmeta)

        url, files = ftp_dir(gse)
        with open(os.path.join(OUT, f"{gse}.suppl.txt"), "w") as fh:
            fh.write(url + "\n" + "\n".join(files) + "\n")
        print(f"{gse}: soft={len(soft)}b all={len(allmeta)}b suppl={len(files)}")
        time.sleep(1)


if __name__ == "__main__":
    main()
