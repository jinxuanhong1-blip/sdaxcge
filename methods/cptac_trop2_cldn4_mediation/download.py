#!/usr/bin/env python3
"""Download CPTAC LUAD and LSCC public TMT tumor protein (freeze v1.2).

Open S3 only. Same filenames as the LinkedOmics CPTAC-pancan index and as
the earlier MHC/IFN protein page. No dbGaP token. Protein matrices are not
committed; manifest.json records the hash.
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
UA = "cptac-trop2-cldn4-mediation"

# Hashes recorded in data/cptac_ifn_mhc_protein/manifest.json on the
# MHC/IFN protein page. A mismatch means this is not that freeze file.
EXPECTED_SHA256 = {
    "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt":
        "746900f54072198481207ea13f8a7ba14dcabd52345bace0bef7438120c03845",
    "LUAD_phenotype.txt":
        "9f777402f7d7d53856a5adfc73c8291136430ddfad85641584426921609b349a",
    "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt":
        "fd9f9cc0aed40942fd33d750d1d207bcdf9ea231fcddeff3cc8ac8fdf3e1384d",
    "LSCC_phenotype.txt":
        "34f17a653dd5ebce0cdf1257c4498fe6579e314d18d435d333e4bcc4db22ef8d",
}

NEEDED = {
    "LUAD": [
        "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "LUAD_phenotype.txt",
        "LUAD_meta.txt",
    ],
    "LSCC": [
        "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "LSCC_phenotype.txt",
        "LSCC_meta.txt",
    ],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def head_length(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
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
        expected = EXPECTED_SHA256.get(dest.name)
        if expected is not None:
            rec["sha256_matches_prior_manifest"] = rec["sha256"] == expected
    except Exception as exc:
        rec["status"] = "failed"
        rec["error"] = str(exc)
        print(f"FAIL {dest.name}: {exc}", file=sys.stderr)
    return rec


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="data/cptac_trop2_cldn4_mediation")
    args = p.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "freeze": "data_freeze_v1.2_reorganized",
        "freeze_base": FREEZE,
        "linkedomics_luad": "https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/",
        "linkedomics_lscc": "https://www.linkedomics.org/data_download/CPTAC-pancan-LSCC/",
        "scope": (
            "TROP2 (TACSTD2) protein vs CD8A and MHC-I protein; "
            "partial correlation given CLDN4; mediation TACSTD2→CLDN4→CD8A/MHC-I; "
            "EPCAM through the same paths. LUAD and LSCC kept separate."
        ),
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": [],
    }
    for cohort, names in NEEDED.items():
        for name in names:
            url = f"{FREEZE}/{cohort}/{name}"
            dest = outdir / cohort / name
            rec = fetch(url, dest)
            rec["cohort"] = cohort
            rec["path"] = str(dest)
            manifest["files"].append(rec)
    man_path = outdir / "manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(man_path)
    bad = [
        f for f in manifest["files"]
        if f.get("status") == "failed" or f.get("sha256_matches_prior_manifest") is False
    ]
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
