#!/usr/bin/env python3
"""Fetch public GEO / ArrayExpress metadata for leftover 2024–2026 lung ICI scRNA candidates.

Public only. Does not download EGA/dbGaP payloads.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "metadata"
OUT.mkdir(parents=True, exist_ok=True)

GSE = [
    "GSE274584",
    "GSE274588",
    "GSE274595",
    "GSE302113",
    "GSE267108",
    "GSE176021",
    "GSE176022",
    "GSE186446",
    "GSE337519",
    "GSE308745",
]

EMTAB = ["E-MTAB-13526"]

UA = "Mozilla/5.0 leftover-scrna-2426/1.0 (public metadata fetch)"


def get(url: str, dest: Path, retries: int = 4) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    last_err = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=90) as r:
                dest.write_bytes(r.read())
            print(f"OK  {dest.name}  {dest.stat().st_size} bytes  <- {url}")
            return
        except Exception as e:
            last_err = e
            wait = 2 ** (i + 1)
            print(f"RETRY {i+1}/{retries} {url} :: {e}; sleep {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"failed {url}: {last_err}")


def main() -> None:
    for acc in GSE:
        n = acc.replace("GSE", "")
        get(
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={acc}[ACCN]&retmode=json",
            OUT / f"{acc}.esearch.json",
        )
        time.sleep(0.4)
        get(
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief",
            OUT / f"{acc}.soft.txt",
        )
        time.sleep(0.4)
        get(
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{n[:-3]}nnn/{acc}/miniml/{acc}_family.xml.tgz",
            OUT / f"{acc}_family.xml.tgz",
        )
        time.sleep(0.4)
        # file list if present
        try:
            get(
                f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{n[:-3]}nnn/{acc}/suppl/filelist.txt",
                OUT / f"{acc}_filelist.txt",
            )
        except Exception as e:
            (OUT / f"{acc}_filelist.txt").write_text(f"MISSING\n{e}\n")
            print(f"NO FILELIST {acc}: {e}")
        time.sleep(0.4)

    # ArrayExpress / BioStudies for Cvejic atlas
    for acc in EMTAB:
        get(
            f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{acc}.json",
            OUT / f"{acc}.biostudies.json",
        )
        time.sleep(0.4)
        try:
            get(
                f"https://www.ebi.ac.uk/biostudies/files/{acc}/{acc}.sdrf.txt",
                OUT / f"{acc}.sdrf.txt",
            )
        except Exception as e:
            print(f"SDRF fail {acc}: {e}")
        try:
            get(
                f"https://www.ebi.ac.uk/biostudies/files/{acc}/{acc}.idf.txt",
                OUT / f"{acc}.idf.txt",
            )
        except Exception as e:
            print(f"IDF fail {acc}: {e}")
        try:
            get(
                f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{acc}/files",
                OUT / f"{acc}.files.html",
            )
        except Exception as e:
            print(f"files html fail {acc}: {e}")

    print("done")


if __name__ == "__main__":
    main()
