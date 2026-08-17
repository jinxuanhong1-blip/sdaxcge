#!/usr/bin/env python3
"""Download public GSE285029 author WTS matrix (n=234)."""

from __future__ import annotations

import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(parents=True, exist_ok=True)

FILES = {
    "GSE285029_WTS_expr_count_235_032820.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285029/"
        "suppl/GSE285029_WTS_expr_count_235_032820.txt.gz"
    ),
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest} ({dest.stat().st_size} bytes)")
        return
    print(f"GET {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "gse285029-cldn4-gsea/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as out:
        out.write(resp.read())
    print(f"wrote {dest} ({dest.stat().st_size} bytes)")


def main() -> None:
    manifest = {
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "source": "NCBI GEO GSE285029 (public)",
        "pmid": "40050048",
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
