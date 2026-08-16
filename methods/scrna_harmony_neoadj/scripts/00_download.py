#!/usr/bin/env python3
"""Download public processed UMI/author matrices only. No FASTQ.

Hard cap: refuse any single file > 2.0 GB and any cumulative download > 8 GB.
GSE241934 Real MTX is ~1.2 GB gzipped (allowed). FASTQ / EGA / GSA are not used.
"""
from __future__ import annotations

import hashlib
import json
import sys
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA.mkdir(parents=True, exist_ok=True)

# Single-file cap (compressed). User rule: do not pull >15 GB FASTQ.
MAX_FILE_GB = 2.0
MAX_TOTAL_GB = 8.0

FILES = [
    # GSE207422 (Hu et al. Genome Med 2023) — User A3 source, taken as given
    {
        "cohort": "GSE207422",
        "name": "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "expected_mb": 176,
    },
    {
        "cohort": "GSE207422",
        "name": "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "expected_mb": 0.02,
    },
    # GSE241934 (NEOTIDE + real-world) — author Epi/T/NK + MPR
    {
        "cohort": "GSE241934",
        "name": "GSE241934_IIT_features.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_features.tsv.gz",
        "expected_mb": 0.2,
    },
    {
        "cohort": "GSE241934",
        "name": "GSE241934_IIT_barcodes.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_barcodes.tsv.gz",
        "expected_mb": 0.4,
    },
    {
        "cohort": "GSE241934",
        "name": "GSE241934_IIT_Meta.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Meta.txt.gz",
        "expected_mb": 1.9,
    },
    {
        "cohort": "GSE241934",
        "name": "GSE241934_IIT_Matrix.mtx.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx.gz",
        "expected_mb": 444,
    },
    {
        "cohort": "GSE241934",
        "name": "GSE241934_RWC_features.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_features.tsv.gz",
        "expected_mb": 0.2,
    },
    {
        "cohort": "GSE241934",
        "name": "GSE241934_RWC_barcodes.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_barcodes.tsv.gz",
        "expected_mb": 1.0,
    },
    {
        "cohort": "GSE241934",
        "name": "GSE241934_Real_Meta.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Meta.txt.gz",
        "expected_mb": 5.4,
    },
    {
        "cohort": "GSE241934",
        "name": "GSE241934_Real_Matrix.mtx.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Matrix.mtx.gz",
        "expected_mb": 1200,
    },
    # GSE291670 (Xia et al. J Transl Med 2025)
    {
        "cohort": "GSE291670",
        "name": "GSE291670_RAW.tar",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291670/suppl/GSE291670_RAW.tar",
        "expected_mb": 122,
    },
    # GSE205335 identity only (palliative ICI; decide inclusion from labels)
    {
        "cohort": "GSE205335",
        "name": "GSE205335_Lung_IO_CellIdentity.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
        "expected_mb": 0.7,
    },
]


def md5(path: Path, nbytes: int = 1_048_576) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            b = f.read(nbytes)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"GET {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "scrna-harmony-neoadj/1.0"})
    with urllib.request.urlopen(req, timeout=600) as r, tmp.open("wb") as out:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)


def main() -> int:
    manifest = []
    total = 0
    for spec in FILES:
        dest = DATA / spec["cohort"] / spec["name"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and dest.stat().st_size > 0:
            print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        else:
            if spec["expected_mb"] / 1024 > MAX_FILE_GB:
                raise SystemExit(f"refusing {spec['name']}: expected {spec['expected_mb']} MB > {MAX_FILE_GB} GB")
            if (total / 1e9 + spec["expected_mb"] / 1024) > MAX_TOTAL_GB:
                raise SystemExit(f"refusing cumulative download > {MAX_TOTAL_GB} GB")
            download(spec["url"], dest)
        size = dest.stat().st_size
        total += size
        rec = {
            "cohort": spec["cohort"],
            "name": spec["name"],
            "url": spec["url"],
            "bytes": size,
            "md5": md5(dest),
        }
        manifest.append(rec)
        print(f"ok {dest.name} {size} md5={rec['md5']}", flush=True)

    tar = DATA / "GSE291670" / "GSE291670_RAW.tar"
    extract_dir = DATA / "GSE291670" / "raw"
    if tar.exists() and not (extract_dir / "_extracted.ok").exists():
        extract_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar, "r") as tf:
            tf.extractall(extract_dir)
        (extract_dir / "_extracted.ok").write_text("ok\n")
        print(f"extracted {tar} -> {extract_dir}", flush=True)

    (DATA / "download_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"total_bytes={total}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
