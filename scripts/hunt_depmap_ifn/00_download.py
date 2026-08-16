#!/usr/bin/env python3
"""Download DepMap 24Q4 public files used by hunt_depmap_ifn.

Source: DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+.
https://doi.org/10.25452/figshare.plus.27993248.v1
Hallmark GMT: MSigDB 2024.1 Hs (Broad / UC San Diego).
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path

FILES = {
    "Model.csv": {
        "url": "https://ndownloader.figshare.com/files/51065297",
        "md5": "675210d17675f3517b0ce39a3c274f16",
        "min_bytes": 100_000,
    },
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": {
        "url": "https://ndownloader.figshare.com/files/51065489",
        "md5": "71794802b750ce77c422dad0720a40af",
        "min_bytes": 400_000_000,
    },
    "OmicsSomaticMutationsMatrixHotspot.csv": {
        "url": "https://ndownloader.figshare.com/files/51065750",
        "md5": "a5aeb1deef897ead3c955c372148d840",
        "min_bytes": 1_000_000,
    },
    "OmicsSomaticMutationsMatrixDamaging.csv": {
        "url": "https://ndownloader.figshare.com/files/51065747",
        "md5": "cb20fdbe1cf3b9b0d8ed4f53e1f399b6",
        "min_bytes": 50_000_000,
    },
    "h.all.v2024.1.Hs.symbols.gmt": {
        "url": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt",
        "md5": None,
        "min_bytes": 10_000,
    },
}


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"GET {url} -> {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "hunt_depmap_ifn/1.0"})
    with urllib.request.urlopen(req, timeout=600) as r, tmp.open("wb") as out:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/hunt_depmap_ifn")
    ap.add_argument("--skip-md5", action="store_true")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for name, meta in FILES.items():
        dest = outdir / name
        if dest.exists() and dest.stat().st_size >= meta["min_bytes"]:
            if meta["md5"] and not args.skip_md5:
                got = md5sum(dest)
                if got != meta["md5"]:
                    print(f"MD5 mismatch {name}: {got} != {meta['md5']}; re-downloading", flush=True)
                    dest.unlink()
                else:
                    print(f"OK exists {name} ({dest.stat().st_size} bytes)", flush=True)
                    continue
            else:
                print(f"OK exists {name} ({dest.stat().st_size} bytes)", flush=True)
                continue
        download(meta["url"], dest)
        if dest.stat().st_size < meta["min_bytes"]:
            print(f"ERROR {name} too small: {dest.stat().st_size}", file=sys.stderr)
            return 1
        if meta["md5"] and not args.skip_md5:
            got = md5sum(dest)
            if got != meta["md5"]:
                print(f"ERROR MD5 {name}: {got} != {meta['md5']}", file=sys.stderr)
                return 1
        print(f"OK {name} ({dest.stat().st_size} bytes)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
