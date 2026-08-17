#!/usr/bin/env python3
"""Download public CPTAC-LSCC (LUSC) processed tables.

Open freeze v1.2 on AWS S3 (LinkedOmics CPTAC-pancan-LSCC filenames).
Protein tumor + matched RNA tumor + phenotype only. No LUAD files.
No dbGaP / CDS tokens.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

FREEZE = (
    "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/"
    "data_freeze_v1.2_reorganized/LSCC"
)
LINKEDOMICS = "https://www.linkedomics.org/data_download/CPTAC-pancan-LSCC/"

FILES = [
    "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "LSCC_phenotype.txt",
]

MAX_BYTES = 2 * 1024 * 1024 * 1024
UA = "cptac-lusc-cldn4-protein"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def head_length(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as resp:
        if resp.status != 200:
            raise RuntimeError(f"HEAD {resp.status} {url}")
        cl = resp.headers.get("Content-Length")
        return int(cl) if cl else None


def fetch(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    rec = {"url": url, "name": dest.name}
    nbytes = head_length(url)
    rec["content_length"] = nbytes
    if nbytes is not None and nbytes > MAX_BYTES:
        rec["status"] = "skipped_gt_2gb"
        return rec
    if dest.exists() and nbytes is not None and dest.stat().st_size == nbytes:
        rec["status"] = "cached"
    else:
        print(f"GET {dest.name}", file=sys.stderr)
        tmp = dest.with_suffix(dest.suffix + ".part")
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=300) as resp, tmp.open("wb") as out:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
        tmp.replace(dest)
        rec["status"] = "downloaded"
    rec["local_bytes"] = dest.stat().st_size
    rec["sha256"] = sha256(dest)
    return rec


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="data/cptac_lusc_cldn4_protein")
    args = p.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "cohort": "LSCC",
        "synonym": "LUSC",
        "freeze": "data_freeze_v1.2_reorganized",
        "freeze_base": FREEZE,
        "linkedomics_index": LINKEDOMICS,
        "luad_excluded": (
            "CPTAC LUAD TJ-15 vs ImmuneScore (ρ=−0.30, n=110; PR #245) "
            "is taken as given and is not re-downloaded."
        ),
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": [],
    }
    for name in FILES:
        rec = fetch(f"{FREEZE}/{name}", outdir / name)
        rec["path"] = str(outdir / name)
        manifest["files"].append(rec)
        print(json.dumps({k: rec[k] for k in rec if k != "sha256"}), flush=True)

    man_path = outdir / "manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(man_path)
    failed = [f for f in manifest["files"] if f.get("status") not in {"downloaded", "cached"}]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
