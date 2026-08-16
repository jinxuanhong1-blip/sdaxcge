#!/usr/bin/env python3
"""Download the open CPTAC LSCC freeze files needed for the xCell/purity slice.

Protein tumor + phenotype + meta only. RNA is not used here.
Filenames match LinkedOmics CPTAC-pancan-LSCC / PR23.
"""

from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "w200" / "CPTAC_LSCC_xcell"
NOTES = ROOT / "results" / "w200" / "CPTAC_LSCC_xcell"

BASE = "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC"
FILES = [
    "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "LSCC_phenotype.txt",
    "LSCC_meta.txt",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(name: str) -> dict:
    url = f"{BASE}/{name}"
    dest = DATA / name
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "w200-cptac-lscc-xcell"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        code = resp.status
        length = int(resp.headers.get("Content-Length") or 0)
    if code != 200:
        raise RuntimeError(f"HEAD {code} {url}")
    if dest.exists() and dest.stat().st_size == length:
        return {"file": name, "url": url, "bytes": length, "cached": True, "sha256": sha256(dest)}
    print(f"GET {name} ({length} bytes)", flush=True)
    urllib.request.urlretrieve(url, dest)
    return {"file": name, "url": url, "bytes": dest.stat().st_size, "cached": False, "sha256": sha256(dest)}


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    NOTES.mkdir(parents=True, exist_ok=True)
    records = [download(name) for name in FILES]
    for rec in records:
        print(json.dumps(rec), flush=True)
    manifest = {
        "freeze": "data_freeze_v1.2_reorganized",
        "bucket": "cptac-pancancer-data",
        "region": "us-west-2",
        "cohort": "LSCC",
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": records,
    }
    (NOTES / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
