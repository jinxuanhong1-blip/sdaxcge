#!/usr/bin/env python3
"""Download the open CPTAC LUAD freeze files needed for this rework.

Protein tumor + RNA tumor + phenotype + meta.
Standalone LUAD_{xcell,cibersort,mcpcounter,estimate}.txt files are HTTP 403.
CIBERSORT and ESTIMATE live in LUAD_phenotype.txt (HTTP 200).
MCP-counter is not a freeze column; analyze.py computes it from public RNA.
"""

from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "rework" / "CPTAC_LUAD_defs"
NOTES = ROOT / "results" / "rework" / "CPTAC_LUAD_defs"

BASE = "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD"
FILES = [
    "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "LUAD_phenotype.txt",
    "LUAD_meta.txt",
]
UA = "rework-cptac-luad-defs"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(name: str) -> dict:
    url = f"{BASE}/{name}"
    dest = DATA / name
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as resp:
        code = resp.status
        length = int(resp.headers.get("Content-Length") or 0)
    if code != 200:
        raise RuntimeError(f"HEAD {code} {url}")
    if dest.exists() and dest.stat().st_size == length:
        return {
            "file": name,
            "url": url,
            "bytes": length,
            "cached": True,
            "http": 200,
            "sha256": sha256(dest),
        }
    print(f"GET {name} ({length} bytes)", flush=True)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    return {
        "file": name,
        "url": url,
        "bytes": dest.stat().st_size,
        "cached": False,
        "http": 200,
        "sha256": sha256(dest),
    }


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
        "cohort": "LUAD",
        "linkedomics_index": "https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/",
        "citation": "Gillette et al. Cell 2020 PMID 32649874; Li et al. Cell Syst 2023 freeze",
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": records,
        "notes": [
            "Standalone LUAD_xcell.txt / LUAD_cibersort.txt / LUAD_mcpcounter.txt / LUAD_estimate.txt HEAD as HTTP 403.",
            "CIBERSORT and ESTIMATE columns are used from LUAD_phenotype.txt.",
            "MCP-counter is computed in analyze.py from the public RNA matrix + Becht 2016 genes.",
        ],
    }
    (NOTES / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (DATA / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
