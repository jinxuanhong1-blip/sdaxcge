#!/usr/bin/env python3
"""Download open protein tables. RAW spectra are not fetched.

ftp.pride.ebi.ac.uk TLS failed from this environment (unexpected EOF).
proteinGroups_Cldn4.txt is kept only when its SHA-1 matches the live PRIDE
file checksum. proteinGroups_PRISMA.txt (66,679,807 bytes) was not retrieved.
"""

from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
DATA.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)
UA = "cldn4-masswave/1.0 (protein-tables; no-raw)"

FILES = [
    {
        "name": "crc24_suppst1_bioid.docx",
        "url": "https://ndownloader.figshare.com/files/53292048",
        "md5": "70bcbba8b14311102973f11eda6cc41d",
    },
    {
        "name": "plos2015_s2_all_proteins.xlsx",
        "url": "https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0117074.s006",
    },
    {
        "name": "plos2015_s3_enriched.xlsx",
        "url": "https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0117074.s007",
    },
]


def sha1(path: Path) -> str:
    h = hashlib.sha1()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def try_download(url: str, dest: Path) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            dest.write_bytes(resp.read())
        return "downloaded"
    except Exception as exc:  # noqa: BLE001
        return f"failed: {exc}"


def main() -> None:
    log = []
    pride_url = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/10/PXD031094/proteinGroups_Cldn4.txt"
    dest = DATA / "proteinGroups_Cldn4.txt"
    status = try_download(pride_url, dest) if not dest.exists() else "already_present"
    digest = sha1(dest) if dest.exists() else ""
    log.append(
        {
            "name": dest.name,
            "url": pride_url,
            "status": status,
            "bytes": dest.stat().st_size if dest.exists() else 0,
            "sha1": digest,
            "pride_api_sha1": "4dbde6c0d827105412ef1f34cbb969b8dd8acb36",
            "sha1_matches_pride_api": digest == "4dbde6c0d827105412ef1f34cbb969b8dd8acb36",
        }
    )
    prisma_url = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/10/PXD031094/proteinGroups_PRISMA.txt"
    prisma_status = try_download(prisma_url, DATA / "proteinGroups_PRISMA.txt")
    log.append(
        {
            "name": "proteinGroups_PRISMA.txt",
            "url": prisma_url,
            "status": prisma_status,
            "pride_api_sha1": "770554ea6737df2a850b08fa2404d0f5be0f75b6",
            "pride_api_bytes": 66679807,
            "note": "C-terminal PRISMA peptide matrix. Not scanned when the download fails.",
        }
    )
    for spec in FILES:
        dest = DATA / spec["name"]
        status = "already_present" if dest.exists() and dest.stat().st_size > 0 else try_download(spec["url"], dest)
        log.append({"name": spec["name"], "url": spec["url"], "status": status, "bytes": dest.stat().st_size if dest.exists() else 0})
    (RESULTS / "download_log.json").write_text(json.dumps(log, indent=2) + "\n")
    print(json.dumps(log, indent=2))


if __name__ == "__main__":
    main()
