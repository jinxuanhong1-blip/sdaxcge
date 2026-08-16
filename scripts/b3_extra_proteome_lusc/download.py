#!/usr/bin/env python3
"""Download public tables for the B3 CPTAC-LUAD protein + remaining proteome + TCGA-LUSC extra.

All URLs are open HTTP 200 sources. No dbGaP / CDS tokens.
LSCC protein is not downloaded (already reported; not re-audited).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from genes import REMAINING_PROTEOME  # noqa: E402

FREEZE = (
    "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/"
    "data_freeze_v1.2_reorganized"
)

XENA_LUSC = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "TCGA.LUSC.sampleMap%2FHiSeqV2.gz"
)
ESTIMATE_LUSC = (
    "https://ibl.mdanderson.org/estimate/tables/"
    "lung_squamous_cell_carcinoma_RNAseqV2.txt"
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


def luad_jobs(outdir: Path) -> list[tuple[str, Path]]:
    names = [
        "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
        "LUAD_phenotype.txt",
    ]
    return [(f"{FREEZE}/LUAD/{n}", outdir / "cptac" / "LUAD" / n) for n in names]


def remaining_jobs(outdir: Path) -> list[tuple[str, Path]]:
    jobs = []
    for c in REMAINING_PROTEOME:
        prot = f"{c}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
        pheno = f"{c}_phenotype.txt"
        jobs.append((f"{FREEZE}/{c}/{prot}", outdir / "cptac" / c / prot))
        jobs.append((f"{FREEZE}/{c}/{pheno}", outdir / "cptac" / c / pheno))
    return jobs


def lusc_jobs(outdir: Path) -> list[tuple[str, Path]]:
    return [
        (XENA_LUSC, outdir / "tcga" / "TCGA.LUSC.HiSeqV2.gz"),
        (ESTIMATE_LUSC, outdir / "tcga" / "ESTIMATE_LUSC_RNAseqV2.txt"),
    ]


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="data/b3_extra_proteome_lusc")
    p.add_argument(
        "--skip-remaining",
        action="store_true",
        help="Download only CPTAC LUAD + TCGA-LUSC (skip other freeze proteomes).",
    )
    args = p.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    jobs = luad_jobs(outdir) + lusc_jobs(outdir)
    if not args.skip_remaining:
        jobs += remaining_jobs(outdir)

    manifest = {
        "freeze": "data_freeze_v1.2_reorganized",
        "freeze_base": FREEZE,
        "linkedomics_index": "https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/",
        "xena_lusc": XENA_LUSC,
        "estimate_lusc": ESTIMATE_LUSC,
        "lscc_excluded": (
            "CPTAC LSCC CLDN4 protein vs ImmuneScore is already public; "
            "not re-downloaded as an audit."
        ),
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": [],
    }
    for url, dest in jobs:
        rec = fetch(url, dest)
        rec["path"] = str(dest)
        manifest["files"].append(rec)

    man_path = outdir / "manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(man_path)
    failed = [f for f in manifest["files"] if f.get("status") == "failed"]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
