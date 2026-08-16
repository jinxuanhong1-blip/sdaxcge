#!/usr/bin/env python3
"""Fetch GEO/ArrayExpress/PRIDE metadata for candidate mouse lung ICI datasets.

Outputs raw metadata pages under notes/mouse/raw_meta/ so the catalog can be
built from real records (no invented accessions or values).
"""
import os
import sys
import time
import urllib.request
import urllib.error

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "notes", "mouse", "raw_meta")
OUT = os.path.abspath(OUT)
os.makedirs(OUT, exist_ok=True)

GEO_ACCS = [
    "GSE239485", "GSE133604", "GSE222158", "GSE129297", "GSE241978",
    "GSE330658", "GSE197260", "GSE297630", "GSE297632",
]
AE_ACCS = ["E-MTAB-13704"]
PXD_ACCS = ["PXD059688"]

UA = {"User-Agent": "Mozilla/5.0 (mouse-ici-catalog)"}


def fetch(url, timeout=60):
    req = urllib.request.Request(url, headers=UA)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:  # noqa: BLE001
            sys.stderr.write(f"  attempt {attempt+1} failed for {url}: {e}\n")
            time.sleep(2 ** attempt)
    return None


def geo_targzlist(acc):
    """Query the GEO supplementary FTP directory listing (via HTTPS)."""
    stub = acc[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{acc}/suppl/"
    return url, fetch(url)


def geo_acc_page(acc):
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
        f"{acc}&targ=self&form=text&view=brief"
    )
    return url, fetch(url)


def main():
    for acc in GEO_ACCS:
        print(f"== {acc} ==")
        u1, page = geo_acc_page(acc)
        if page:
            with open(os.path.join(OUT, f"{acc}.soft.txt"), "w") as f:
                f.write(f"# URL: {u1}\n{page}")
        u2, suppl = geo_targzlist(acc)
        if suppl:
            with open(os.path.join(OUT, f"{acc}.suppl.html"), "w") as f:
                f.write(f"<!-- URL: {u2} -->\n{suppl}")
        time.sleep(0.5)

    for acc in AE_ACCS:
        print(f"== {acc} ==")
        url = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{acc}"
        page = fetch(url)
        if page:
            with open(os.path.join(OUT, f"{acc}.json"), "w") as f:
                f.write(page)
        # files listing
        furl = f"https://www.ebi.ac.uk/biostudies/files/{acc}"
        fpage = fetch(furl)
        if fpage:
            with open(os.path.join(OUT, f"{acc}.files.json"), "w") as f:
                f.write(fpage)
        time.sleep(0.5)

    for acc in PXD_ACCS:
        print(f"== {acc} ==")
        url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}"
        page = fetch(url)
        if page:
            with open(os.path.join(OUT, f"{acc}.json"), "w") as f:
                f.write(page)
        furl = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files?pageSize=200"
        fpage = fetch(furl)
        if fpage:
            with open(os.path.join(OUT, f"{acc}.files.json"), "w") as f:
                f.write(fpage)
        time.sleep(0.5)

    print("done ->", OUT)


if __name__ == "__main__":
    main()
