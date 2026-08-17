#!/usr/bin/env python3
"""Download public GSE166449 TPM matrix and GEO metadata (n=22)."""

from __future__ import annotations

import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(parents=True, exist_ok=True)

FILES = {
    "GSE166449_Raw_gene_TPM_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/"
        "suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz"
    ),
    "GSE166449_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/"
        "matrix/GSE166449_series_matrix.txt.gz"
    ),
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest} ({dest.stat().st_size} bytes)")
        return
    print(f"GET {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "gse166449-cldn4-gsea/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as out:
        out.write(resp.read())
    print(f"wrote {dest} ({dest.stat().st_size} bytes)")


def main() -> None:
    manifest = {
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "source": "NCBI GEO GSE166449 (public)",
        "pmid": "33857424",
        "files": {},
    }
    for name, url in FILES.items():
        dest = DATA / name
        fetch(url, dest)
        manifest["files"][name] = {
            "url": url,
            "bytes": dest.stat().st_size,
            "path": str(dest.relative_to(HERE)),
        }
    (DATA / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("done")


if __name__ == "__main__":
    main()
