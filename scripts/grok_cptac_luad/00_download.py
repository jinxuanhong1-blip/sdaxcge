#!/usr/bin/env python3
"""Download CPTAC LUAD open S3 freeze v1.2 protein, RNA, phenotype, clinical files.

Source (verified HTTP 200, 2026-08-16):
https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/

Filenames were taken from the LinkedOmics CPTAC-pancan-LUAD download page
(https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/), which points
at this same S3 prefix. Do not invent accessions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE = (
    "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/"
    "data_freeze_v1.2_reorganized/LUAD"
)

# Protein + RNA + immune/clinical only. Skip isoform/phospho/methylation/CNV/MAF.
FILES = [
    "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt",
    "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt",
    "LUAD_phenotype.txt",
    "LUAD_meta.txt",
    "LUAD_survival.txt",
]

MAX_BYTES = 2 * 1024 * 1024 * 1024


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


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--outdir",
        default="data/grok_cptac_luad",
        help="Local cache (not committed).",
    )
    args = p.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "source_base": BASE,
        "linkedomics_index": "https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/",
        "freeze": "data_freeze_v1.2_reorganized",
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": [],
    }

    for name in FILES:
        url = f"{BASE}/{name}"
        dest = outdir / name
        rec = {"name": name, "url": url}
        try:
            nbytes = head_length(url)
            rec["content_length"] = nbytes
            if nbytes is not None and nbytes > MAX_BYTES:
                rec["status"] = "skipped_gt_2gb"
                print(f"SKIP >2GB {name} ({nbytes} bytes)", file=sys.stderr)
                manifest["files"].append(rec)
                continue
            if dest.exists() and nbytes is not None and dest.stat().st_size == nbytes:
                rec["status"] = "cached"
            else:
                print(f"GET {name}", file=sys.stderr)
                download(url, dest)
                rec["status"] = "downloaded"
            rec["local_bytes"] = dest.stat().st_size
            rec["sha256"] = sha256(dest)
        except Exception as exc:
            rec["status"] = "failed"
            rec["error"] = str(exc)
            print(f"FAIL {name}: {exc}", file=sys.stderr)
        manifest["files"].append(rec)

    man_path = outdir / "manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(man_path)
    failed = [f for f in manifest["files"] if f["status"] == "failed"]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
