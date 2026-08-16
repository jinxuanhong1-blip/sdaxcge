#!/usr/bin/env python3
"""Download CPTAC LSCC open S3 freeze v1.2 matrices (protein, RNA, phenotype).

Source (LinkedOmics CPTAC pan-cancer LSCC download page, verified HTTP 200):
https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/
"""

from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "grok_cptac_lscc"
NOTES = ROOT / "notes" / "grok_cptac_lscc"

BASE = "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC"
FILES = [
    "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt",
    "LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt",
    "LSCC_phenotype.txt",
    "LSCC_survival.txt",
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
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "grok-cptac-lscc"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        code = resp.status
        length = int(resp.headers.get("Content-Length") or 0)
    if code != 200:
        raise RuntimeError(f"HEAD {code} {url}")
    if length > 2 * 1024**3:
        return {"file": name, "url": url, "skipped": True, "reason": f">2GB ({length})", "bytes": length}
    if dest.exists() and dest.stat().st_size == length:
        rec = {"file": name, "url": url, "bytes": length, "cached": True, "sha256": sha256(dest)}
        return rec
    print(f"GET {name} ({length} bytes)", flush=True)
    urllib.request.urlretrieve(url, dest)
    rec = {"file": name, "url": url, "bytes": dest.stat().st_size, "cached": False, "sha256": sha256(dest)}
    return rec


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    NOTES.mkdir(parents=True, exist_ok=True)
    records = []
    for name in FILES:
        rec = download(name)
        records.append(rec)
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
