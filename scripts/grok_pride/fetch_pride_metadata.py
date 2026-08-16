#!/usr/bin/env python3
"""Fetch PRIDE / ProteomeXchange metadata and processed-file lists (no raw MS)."""

from __future__ import annotations

import json
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path("results/grok_pride/raw_api")
OUT.mkdir(parents=True, exist_ok=True)

UA = "grok-pride-slice/1.0 (research; no-raw-ms; +https://github.com/jinxuanhong1-blip/sdaxcge)"
CTX = ssl.create_default_context()

CORE = ["PXD042091", "PXD059688", "PXD019061"]
CANDIDATES = [
    "PXD019573",  # NSCLC PA200 / Durvalumab-associated
    "PXD020191",  # Lehtiö NSCLC proteogenomics (immune phenotypes; not ICI-treated)
    "PXD036226",
    "PXD033237",
    "PXD028354",
    "PXD016739",
    "PXD014882",
    "PXD029847",
    "PXD041312",
    "PXD045218",
    "PXD048912",
    "PXD051234",
]


def get(url: str, timeout: int = 60) -> tuple[int, bytes, str]:
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json,text/xml,*/*"})
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
            return r.status, r.read(), r.headers.get("Content-Type", "")
    except urllib.error.HTTPError as e:
        return e.code, e.read() if e.fp else b"", str(e)
    except Exception as e:
        return 0, str(e).encode(), "error"


def save_bytes(name: str, data: bytes) -> Path:
    p = OUT / name
    p.write_bytes(data)
    return p


def fetch_json(url: str) -> object | None:
    code, body, _ = get(url)
    if code != 200 or not body:
        return {"_http": code, "_url": url, "_preview": body[:400].decode("utf-8", "replace")}
    try:
        return json.loads(body.decode("utf-8", "replace"))
    except json.JSONDecodeError:
        return {"_http": code, "_url": url, "_text": body[:2000].decode("utf-8", "replace")}


def main() -> None:
    endpoints = {
        "project": "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}",
        "files": "https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={acc}",
        "files_page": "https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={acc}&pageSize=200&page=0",
        "px_json": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={acc}&outputMode=JSON",
        "px_xml": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={acc}&outputMode=XML",
    }

    searches = {
        "search_immunotherapy_nsclc": "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=NSCLC%20immunotherapy&pageSize=50&page=0",
        "search_lung_immunotherapy": "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=lung%20cancer%20immunotherapy&pageSize=50&page=0",
        "search_pembrolizumab": "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=pembrolizumab%20lung&pageSize=30&page=0",
        "search_anti_pd1_nsclc": "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=anti-PD1%20NSCLC&pageSize=30&page=0",
        "search_checkpoint_lung": "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=immune%20checkpoint%20lung&pageSize=30&page=0",
        "search_ici_lung": "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=ICI%20lung%20proteomics&pageSize=30&page=0",
    }

    for name, url in searches.items():
        print(f"SEARCH {name}")
        data = fetch_json(url)
        save_bytes(f"{name}.json", json.dumps(data, indent=2).encode())
        time.sleep(0.3)

    all_acc = CORE + CANDIDATES
    for acc in all_acc:
        print(f"==== {acc} ====")
        for key, tmpl in endpoints.items():
            url = tmpl.format(acc=acc)
            code, body, ctype = get(url)
            ext = "xml" if "xml" in key else "json"
            dest = OUT / f"{acc}_{key}.{ext}"
            dest.write_bytes(body)
            print(f"  {key}: HTTP {code} {len(body)} bytes {ctype}")
            time.sleep(0.2)

    # PRIDE HTTPS file index pages (HTML directory listings)
    ftp_guesses = {
        "PXD042091": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD042091/",
        "PXD059688": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD059688/",
        "PXD019573": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/10/PXD019573/",
    }
    for acc, url in ftp_guesses.items():
        code, body, _ = get(url)
        save_bytes(f"{acc}_ftp_index.html", body)
        print(f"FTP {acc}: HTTP {code} {len(body)} bytes")

    print("done")


if __name__ == "__main__":
    main()
