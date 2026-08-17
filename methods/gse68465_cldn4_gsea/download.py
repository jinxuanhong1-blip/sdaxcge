#!/usr/bin/env python3
"""Download public GSE68465 series matrix + official GPL96 annotation."""
from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(parents=True, exist_ok=True)

FILES = {
    "GSE68465_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68465/matrix/"
        "GSE68465_series_matrix.txt.gz"
    ),
    "GPL96.annot.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"
    ),
}


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    manifest = {}
    for name, url in FILES.items():
        dest = DATA / name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"have {dest} ({dest.stat().st_size} bytes)")
        else:
            print(f"fetch {url}")
            urllib.request.urlretrieve(url, dest)
            print(f"wrote {dest} ({dest.stat().st_size} bytes)")
        manifest[name] = {"url": url, "bytes": dest.stat().st_size, "md5": md5(dest)}
    (DATA / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("done")


if __name__ == "__main__":
    main()
