#!/usr/bin/env python3
"""Download CPTAC LUAD public TMT + RNA + phenotype (freeze v1.2).

LUAD extra only. LSCC Q4 is already reported in methods/cptac_cldn4_q4.
Open S3. Filenames from the LinkedOmics CPTAC-pancan index. No dbGaP tokens.
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
    "data_freeze_v1.2_reorganized"
)
MAX_BYTES = 2 * 1024 * 1024 * 1024

NEEDED = [
    "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "LUAD_phenotype.txt",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def head_length(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as resp:
        cl = resp.headers.get("Content-Length")
        return int(cl) if cl else None


def fetch(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    rec = {"url": url, "name": dest.name}
    try:
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
            urllib.request.urlretrieve(url, tmp)
            tmp.replace(dest)
            rec["status"] = "downloaded"
        rec["local_bytes"] = dest.stat().st_size
        rec["sha256"] = sha256(dest)
    except Exception as exc:
        rec["status"] = "failed"
        rec["error"] = str(exc)
        print(f"FAIL {dest.name}: {exc}", file=sys.stderr)
    return rec


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="data/cptac_luad_cldn4_q4")
    args = p.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "freeze": "data_freeze_v1.2_reorganized",
        "freeze_base": FREEZE,
        "linkedomics_luad": "https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/",
        "scope": "LUAD extra: CLDN4 protein Q4 vs ImmuneScore / CD274. LSCC Q4 already reported.",
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": [],
    }
    for name in NEEDED:
        url = f"{FREEZE}/LUAD/{name}"
        dest = outdir / "LUAD" / name
        rec = fetch(url, dest)
        rec["cohort"] = "LUAD"
        rec["path"] = str(dest)
        manifest["files"].append(rec)

    man_path = outdir / "manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(man_path)
    failed = [f for f in manifest["files"] if f.get("status") == "failed"]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
