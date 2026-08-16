#!/usr/bin/env python3
"""Download the exact public Nusinow 2020 / Gygi CCLE proteomics tables.

Primary table (Protein Quantitation, TSV format of manuscript Table S2):
  https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz

Sample information (manuscript Table S1):
  https://gygi.hms.harvard.edu/data/ccle/Table_S1_Sample_Information.xlsx

Lineage labels for NSCLC vs SCLC (DepMap Public 24Q4 Model.csv):
  https://ndownloader.figshare.com/files/51065297
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

FILES = {
    "protein_quant_current_normalized.csv.gz": {
        "url": "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz",
        "role": "Nusinow 2020 / Gygi normalized protein quantitation (Table S2 TSV)",
        "citation": "Nusinow et al. Cell 2020;180:387-402.e16",
    },
    "Table_S1_Sample_Information.xlsx": {
        "url": "https://gygi.hms.harvard.edu/data/ccle/Table_S1_Sample_Information.xlsx",
        "role": "Nusinow 2020 Table S1 sample information (Tissue of Origin)",
        "citation": "Nusinow et al. Cell 2020;180:387-402.e16",
    },
    "Model.csv": {
        "url": "https://ndownloader.figshare.com/files/51065297",
        "role": "DepMap Public 24Q4 Model.csv (Oncotree NSCLC / SCLC labels)",
        "citation": "DepMap Public 24Q4, Figshare+ 10.25452/figshare.plus.27993248.v1 file 51065297",
    },
}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "B2_n118-recompute/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r, dest.open("wb") as out:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="data/B2_n118")
    args = p.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = {}
    for name, meta in FILES.items():
        dest = outdir / name
        print(f"GET {meta['url']} -> {dest}", flush=True)
        if dest.exists() and dest.stat().st_size > 0:
            print(f"  exists ({dest.stat().st_size} bytes), skip", flush=True)
        else:
            download(meta["url"], dest)
        rec = {
            **meta,
            "path": name,
            "bytes": dest.stat().st_size,
            "sha256": sha256_of(dest),
        }
        manifest[name] = rec
        print(f"  {rec['bytes']} bytes  sha256={rec['sha256']}", flush=True)

    (outdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("wrote", outdir / "download_manifest.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
