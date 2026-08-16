#!/usr/bin/env python3
"""Fetch GEO sample tables + matrix URLs for C3-relevant accessions.

Writes omics/geo_samples/<GSE>.tsv and logs/geo_meta/<GSE>.json
"""

from __future__ import annotations

import json
import os
import time
import urllib.request
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "omics", "geo_samples")
RAW = os.path.join(ROOT, "logs", "geo_meta")
UA = "C3-PDL1-catalog/1.0"

ACCESSIONS = [
    "GSE50927", "GSE207704", "GSE22493", "GSE22421",
    "GSE274940", "GSE245459", "GSE334497",
    "GSE278664",  # SG + berzosertib ovarian
    "GSE222450",  # SN-38 + anti-PD1 HNSCC
    "GSE309617", "GSE309616",
    "GSE302284",
    "GSE312098", "GSE311016",
]


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=90) as fh:
        return fh.read()


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RAW, exist_ok=True)
    for acc in ACCESSIONS:
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=gds&term={acc}[Accession]&retmode=json"
        )
        ids = json.loads(get(url))["esearchresult"]["idlist"]
        time.sleep(0.34)
        # fetch GSE record via esummary + efetch for samples
        url2 = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=gds&term={acc}[Accession]+AND+gse[Entry+Type]&retmode=json"
        )
        gse_ids = json.loads(get(url2))["esearchresult"]["idlist"]
        time.sleep(0.34)
        samples = []
        # Use GEO SOFT brief via ftp series matrix header is later;
        # here pull GSM list from esearch of the series.
        url3 = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=gds&term={acc}[Accession]&retmax=200&retmode=json"
        )
        all_ids = json.loads(get(url3))["esearchresult"]["idlist"]
        time.sleep(0.34)
        if all_ids:
            for i in range(0, len(all_ids), 20):
                chunk = ",".join(all_ids[i:i + 20])
                url4 = (
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                    f"?db=gds&id={chunk}&retmode=json"
                )
                payload = json.loads(get(url4))
                time.sleep(0.34)
                res = payload.get("result", {})
                for uid in res.get("uids", []):
                    rec = res.get(uid, {})
                    samples.append({
                        "uid": uid,
                        "accession": rec.get("accession", ""),
                        "title": rec.get("title", ""),
                        "n_samples": rec.get("n_samples", ""),
                        "gdstype": rec.get("gdstype", ""),
                        "taxon": rec.get("taxon", ""),
                        "summary": (rec.get("summary") or "")[:300],
                    })
        with open(os.path.join(RAW, f"{acc}.json"), "w") as fh:
            json.dump({"accession": acc, "ids": all_ids, "records": samples}, fh, indent=1)
        with open(os.path.join(OUT, f"{acc}.tsv"), "w") as fh:
            fh.write("accession\tgdstype\ttaxon\tn_samples\ttitle\n")
            for s in samples:
                title = s["title"].replace("\t", " ")
                fh.write(f"{s['accession']}\t{s['gdstype']}\t{s['taxon']}\t{s['n_samples']}\t{title}\n")
        print(f"{acc}: {len(samples)} records")
        for s in samples[:12]:
            print(f"  {s['accession']:12} {s['gdstype'][:28]:28} {s['title'][:90]}")


if __name__ == "__main__":
    main()
