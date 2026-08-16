#!/usr/bin/env python3
"""Download CPTAC LUAD + LSCC tumor RNA, protein, and phenotype (freeze v1.2).

Source (HEAD HTTP 200, 2026-08-16):
https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/{LUAD,LSCC}/

Filenames match LinkedOmics CPTAC-pancan-{LUAD,LSCC} download pages, which
point at this S3 prefix. Tumor-only files — this slice does not use NAT.
Do not invent accessions. Cache is not committed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

S3_ROOT = (
    "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/"
    "data_freeze_v1.2_reorganized"
)
LINKEDOMICS = {
    "LUAD": "https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/",
    "LSCC": "https://www.linkedomics.org/data_download/CPTAC-pancan-LSCC/",
}
COHORTS = ("LUAD", "LSCC")
FILE_TEMPLATES = (
    "{c}_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "{c}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "{c}_phenotype.txt",
)
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
        default="data/hunt_cptac_partial",
        help="Local cache (not committed).",
    )
    args = p.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "source_root": S3_ROOT,
        "linkedomics_index": LINKEDOMICS,
        "freeze": "data_freeze_v1.2_reorganized",
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": [],
    }

    for cohort in COHORTS:
        for tmpl in FILE_TEMPLATES:
            name = tmpl.format(c=cohort)
            url = f"{S3_ROOT}/{cohort}/{name}"
            dest = outdir / cohort / name
            rec = {"cohort": cohort, "name": name, "url": url}
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
