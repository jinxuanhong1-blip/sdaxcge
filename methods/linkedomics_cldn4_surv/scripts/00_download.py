#!/usr/bin/env python3
"""Download open LinkedOmics CPTAC freeze v1.2 files needed for CLDN4 OS/PFS.

Clinical tables (survival + meta) plus tumor RNA and tumor protein only.
Does not download phenotype / ImmuneScore matrices (already reported elsewhere).

Source (LinkedOmics download pages point at this S3 prefix; HEAD-checked):
https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/
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
    "data_freeze_v1.2_reorganized"
)

LINKEDOMICS = {
    "LUAD": "https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/",
    "LSCC": "https://www.linkedomics.org/data_download/CPTAC-pancan-LSCC/",
}

# Survival + meta (clinical) and CLDN4-bearing tumor matrices only.
FILES = {
    "LUAD": [
        "LUAD_survival.txt",
        "LUAD_meta.txt",
        "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
        "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    ],
    "LSCC": [
        "LSCC_survival.txt",
        "LSCC_meta.txt",
        "LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
        "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    ],
}

MAX_BYTES = 2 * 1024 * 1024 * 1024


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def head(url: str) -> tuple[int, int | None]:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "linkedomics-cldn4-surv"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        cl = resp.headers.get("Content-Length")
        return int(resp.status), (int(cl) if cl else None)


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--outdir",
        default=str(Path(__file__).resolve().parents[1] / "data"),
        help="Local cache (not committed).",
    )
    args = p.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "source_base": BASE,
        "linkedomics_index": LINKEDOMICS,
        "freeze": "data_freeze_v1.2_reorganized",
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "CLDN4 RNA/protein + OS/PFS clinical tables only; no ImmuneScore",
        "files": [],
    }

    for cohort, names in FILES.items():
        for name in names:
            url = f"{BASE}/{cohort}/{name}"
            dest = outdir / name
            rec = {"cohort": cohort, "name": name, "url": url}
            try:
                code, nbytes = head(url)
                rec["http_status"] = code
                rec["content_length"] = nbytes
                if code != 200:
                    rec["status"] = "not_open"
                    print(f"NOT OPEN HTTP {code} {name}", file=sys.stderr)
                    manifest["files"].append(rec)
                    continue
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
    failed = [f for f in manifest["files"] if f["status"] in {"failed", "not_open"}]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
