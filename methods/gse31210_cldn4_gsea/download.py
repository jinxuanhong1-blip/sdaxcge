#!/usr/bin/env python3
"""Download public GSE31210 series matrix + GPL570 annotation."""

from __future__ import annotations

import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(parents=True, exist_ok=True)

FILES = {
    "GSE31210_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/"
        "matrix/GSE31210_series_matrix.txt.gz"
    ),
    "GPL570.annot.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
    ),
}

UA = (
    "sdaxcge-gse31210-cldn4-gsea/1.0 "
    "(+https://github.com/jinxuanhong1-blip/sdaxcge)"
)


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 10_000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)")
        return
    print(f"GET {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)")


def main() -> None:
    manifest = {
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "source": "NCBI GEO GSE31210 (public) + GPL570.annot",
        "pmid": "22261853",
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
