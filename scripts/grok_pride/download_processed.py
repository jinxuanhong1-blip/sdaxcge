#!/usr/bin/env python3
"""Download processed PRIDE/PX tables only. Never download raw MS (.raw/.wiff/.scan/.d)."""

from __future__ import annotations

import json
import ssl
import time
import urllib.error
import urllib.request
from pathlib import Path

UA = "grok-pride-slice/1.0 (processed-only; no-raw-ms)"
CTX = ssl.create_default_context()
OUT = Path("results/grok_pride/processed")
OUT.mkdir(parents=True, exist_ok=True)
API = Path("results/grok_pride/raw_api")

FILES = [
    # PXD042091 SEARCH / library (text tables, not vendor raw)
    (
        "PXD042091",
        "2_dat_cs_all.txt",
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD042091/2_dat_cs_all.txt",
    ),
    (
        "PXD042091",
        "6_data_tf_response.txt",
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/02/PXD042091/6_data_tf_response.txt",
    ),
    # PXD059688 mzTab processed identifications
    (
        "PXD059688",
        "20230904_Tumor_AA.mzTab.gz",
        "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD059688/20230904_Tumor_AA.mzTab.gz",
    ),
]

# extra accessions to metadata-fetch
EXTRA = [
    "PXD044740",
    "PXD019774",
    "PXD039141",
    "PXD028364",
    "PXD037365",
    "PXD035347",
]


def get(url: str, timeout: int = 90) -> tuple[int, bytes]:
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read() if e.fp else b""
    except Exception as e:
        return 0, str(e).encode()


def download(url: str, dest: Path, timeout: int = 300) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"skip exists {dest} ({dest.stat().st_size} bytes)")
        return
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    print(f"GET {url}")
    with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r, dest.open("wb") as f:
        while True:
            chunk = r.read(1024 * 256)
            if not chunk:
                break
            f.write(chunk)
    print(f"  -> {dest} {dest.stat().st_size} bytes")


def main() -> None:
    for acc, name, url in FILES:
        dest = OUT / acc / name
        try:
            download(url, dest)
        except Exception as e:
            print(f"FAIL {acc} {name}: {e}")

    for acc in EXTRA:
        for key, tmpl in {
            "project": "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}",
            "px_json": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={acc}&outputMode=JSON",
        }.items():
            code, body = get(tmpl.format(acc=acc))
            (API / f"{acc}_{key}.json").write_bytes(body)
            print(f"{acc} {key} HTTP {code} {len(body)}")
            time.sleep(0.2)

    # paper / repository landing pages
    pages = {
        "PXD019061_px_page.html": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD019061",
        "PXD019061_sci_rep.html": "https://www.nature.com/articles/s41598-020-66902-0",
        "PXD059688_frontiers.html": "https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1512605/full",
        "PXD042091_mcp.html": "https://www.mcponline.org/article/S1535-9476(24)00123-4/fulltext",
        "PXD019573_natcancer.html": "https://www.nature.com/articles/s43018-023-00557-4",
        "PXD039141_px.json": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD039141&outputMode=JSON",
        "PXD044740_px.json": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD044740&outputMode=JSON",
        "PXD019774_px.json": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD019774&outputMode=JSON",
    }
    for name, url in pages.items():
        code, body = get(url, timeout=60)
        (API / name).write_bytes(body)
        print(f"PAGE {name} HTTP {code} {len(body)}")
        time.sleep(0.3)

    print("done")


if __name__ == "__main__":
    main()
