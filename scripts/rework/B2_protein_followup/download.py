#!/usr/bin/env python3
"""Download public CCLE/DepMap protein matrices for the B2 follow-up.

Sources (portal is bot-gated; these are direct public files):
  Gygi/Nusinow CCLE MS (Cell 2020)
    https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz
    https://gygi.hms.harvard.edu/data/ccle/Table_S1_Sample_Information.xlsx
  CCLE RPPA 20180123
    https://data.broadinstitute.org/ccle/CCLE_RPPA_20180123.csv
    https://data.broadinstitute.org/ccle/CCLE_RPPA_Ab_info_20180123.csv
  ProCan-DepMapSanger 949-line MS (Gonçalves et al., Cancer Cell 2022)
    figshare 19345397 averaged 8498-protein matrix + mapping
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

UA = "sdaxcge-rework-B2-protein/1.0"
FILES = {
    "protein_quant_current_normalized.csv.gz": "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz",
    "Table_S1_Sample_Information.xlsx": "https://gygi.hms.harvard.edu/data/ccle/Table_S1_Sample_Information.xlsx",
    "CCLE_RPPA_20180123.csv": "https://data.broadinstitute.org/ccle/CCLE_RPPA_20180123.csv",
    "CCLE_RPPA_Ab_info_20180123.csv": "https://data.broadinstitute.org/ccle/CCLE_RPPA_Ab_info_20180123.csv",
    "ProCan-DepMapSanger_protein_matrix_8498_averaged.txt": "https://ndownloader.figshare.com/files/34411172",
    "ProCan-DepMapSanger_mapping_file_averaged.txt": "https://ndownloader.figshare.com/files/34411133",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"GET {url} -> {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=600) as resp, dest.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/ccle_protein")
    p.add_argument("--out-dir", default="results/rework/B2_protein_followup")
    args = p.parse_args()
    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    manifest = {"files": {}}
    for name, url in FILES.items():
        dest = cache / name
        if not dest.exists() or dest.stat().st_size < 1000:
            download(url, dest)
        manifest["files"][name] = {
            "url": url,
            "bytes": dest.stat().st_size,
            "sha256": sha256_file(dest),
        }
    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
