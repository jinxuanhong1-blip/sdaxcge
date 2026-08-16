#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
geo_query.py  -  find GEO Series for a gene perturbation, using only stdlib.

Queries NCBI E-utilities (esearch + esummary on db=gds) so you get REAL GEO
accessions (never invent them). Filters to Series (GSE), reports assay type,
organism and sample count, and can list the supplementary files + sizes so you
can pick a processed matrix that is small enough to download (<2 GB).

EXAMPLES
    # CLDN4 knockdown/knockout RNA-seq series
    python geo_query.py "CLDN4 knockdown" "CLDN4 knockout" "CLDN4 CRISPR"

    # TACSTD2 / TROP2 perturbations, list suppl files for the top hits
    python geo_query.py "TACSTD2 shRNA" "TROP2 knockdown" --suppl

Tips:
    * Set NCBI_API_KEY env var to raise the rate limit from 3 to 10 req/s.
    * Refine on the GEO website with fielded queries, e.g.
        CLDN4[Description] AND "expression profiling by high throughput
        sequencing"[DataSet Type] AND Homo sapiens[Organism]
"""
from __future__ import annotations
import argparse
import os
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{acc}/suppl/"
API_KEY = os.environ.get("NCBI_API_KEY", "")


def _get(url: str) -> str:
    if API_KEY:
        url += ("&" if "?" in url else "?") + "api_key=" + API_KEY
    req = urllib.request.Request(url, headers={"User-Agent": "geo_query/1.0"})
    return urllib.request.urlopen(req, timeout=60).read().decode("utf-8", "replace")


def esearch(term: str, retmax: int = 25) -> list[str]:
    q = urllib.parse.quote(term)
    url = f"{EUTILS}/esearch.fcgi?db=gds&term={q}&retmax={retmax}"
    return re.findall(r"<Id>(\d+)</Id>", _get(url))


def esummary(ids: list[str]) -> list[dict]:
    if not ids:
        return []
    url = f"{EUTILS}/esummary.fcgi?db=gds&id={','.join(ids)}"
    out = []
    for d in re.findall(r"<DocSum>.*?</DocSum>", _get(url), re.S):
        def g(name):
            m = re.search(r'<Item Name="%s"[^>]*>(.*?)</Item>' % name, d, re.S)
            return (m.group(1).strip() if m else "")
        out.append({
            "acc": g("Accession"), "title": g("title"), "taxon": g("taxon"),
            "type": g("gdsType"), "n": g("n_samples"), "summary": g("summary"),
        })
    return out


def list_suppl(acc: str) -> list[tuple[str, str]]:
    """Return [(filename, content-length)] for a GSE's suppl/ directory."""
    if not acc.startswith("GSE"):
        return []
    stub = acc[:-3] if len(acc) > 6 else acc[:3]  # GSE207704 -> GSE207
    stub = re.sub(r"\d{3}$", "", acc)
    url = FTP.format(stub=stub, acc=acc)
    try:
        html = _get(url)
    except Exception as e:
        return [("<listing failed: %r>" % e, "")]
    files = [f for f in re.findall(r'href="([^"]+)"', html)
             if not f.startswith("/") and "http" not in f and f != "../"]
    rows = []
    for f in files:
        try:
            req = urllib.request.Request(url + f, method="HEAD",
                                         headers={"User-Agent": "geo_query/1.0"})
            cl = urllib.request.urlopen(req, timeout=60).headers.get("Content-Length", "")
        except Exception:
            cl = ""
        rows.append((f, cl))
    return rows


def human(nbytes: str) -> str:
    if not nbytes.isdigit():
        return nbytes
    n = int(nbytes)
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f}{unit}"
        n /= 1024


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("terms", nargs="+")
    ap.add_argument("--retmax", type=int, default=25)
    ap.add_argument("--series-only", action="store_true", default=True)
    ap.add_argument("--suppl", action="store_true", help="list suppl files + sizes")
    args = ap.parse_args()

    seen = set()
    for term in args.terms:
        print(f"\n########## QUERY: {term} ##########")
        ids = esearch(term, args.retmax)
        time.sleep(0.4 if API_KEY else 0.75)
        for r in esummary(ids):
            if args.series_only and not r["acc"].startswith("GSE"):
                continue
            if r["acc"] in seen:
                continue
            seen.add(r["acc"])
            print(f"{r['acc']} | {r['taxon']} | n={r['n']} | {r['type'][:48]}")
            print(f"   {r['title'][:150]}")
            if args.suppl:
                for fn, cl in list_suppl(r["acc"]):
                    print(f"      - {fn}  ({human(cl)})")
                time.sleep(0.4)
        time.sleep(0.4 if API_KEY else 0.75)


if __name__ == "__main__":
    main()
