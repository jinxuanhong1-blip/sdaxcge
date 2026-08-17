#!/usr/bin/env python3
"""Download public TCGA-LUSC RNA and official ESTIMATE scores.

Sources (same URLs as prior LUSC slices in this repo):
  - UCSC Xena TCGA-LUSC HiSeqV2 log2(RSEM normalized_count + 1)
  - MD Anderson official ESTIMATE RNAseqV2 ImmuneScore
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from urllib.request import Request, urlopen

URLS = {
    "expr": {
        "url": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap%2FHiSeqV2.gz",
        "file": "TCGA.LUSC.HiSeqV2.gz",
        "min_bytes": 1_000_000,
        "desc": "UCSC Xena TCGA-LUSC HiSeqV2, log2(RSEM normalized_count + 1)",
    },
    "estimate": {
        "url": "https://ibl.mdanderson.org/estimate/tables/lung_squamous_cell_carcinoma_RNAseqV2.txt",
        "file": "ESTIMATE_LUSC_RNAseqV2.txt",
        "min_bytes": 5_000,
        "desc": "Official ESTIMATE RNAseqV2 scores, TCGA-LUSC (Yoshihara 2013 / MD Anderson)",
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, min_bytes: int) -> None:
    if dest.exists() and dest.stat().st_size >= min_bytes:
        print(f"  cached {dest.name} ({dest.stat().st_size} bytes)")
        return
    last = None
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(5):
        try:
            req = Request(url, headers={"User-Agent": "tcga-lusc-cldn4-q4/1.0"})
            with urlopen(req, timeout=180) as r:
                data = r.read()
            if len(data) < min_bytes:
                raise RuntimeError(f"too small: {len(data)} bytes from {url}")
            tmp.write_bytes(data)
            tmp.replace(dest)
            print(f"  downloaded {dest.name} ({len(data)} bytes)")
            return
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2**attempt)
    raise RuntimeError(f"failed to download {url}: {last}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/tcga_lusc_cldn4_q4")
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    provenance = {}
    for key, meta in URLS.items():
        dest = out / meta["file"]
        download(meta["url"], dest, meta["min_bytes"])
        provenance[key] = {
            "url": meta["url"],
            "file": meta["file"],
            "desc": meta["desc"],
            "bytes": dest.stat().st_size,
            "sha256": sha256(dest),
            "downloaded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({k: {"bytes": v["bytes"], "sha256": v["sha256"]} for k, v in provenance.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
