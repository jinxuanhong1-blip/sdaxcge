#!/usr/bin/env python3
"""Download the GSE180963 processed 10x matrices from GEO and unpack them.

GSE180963: scRNA-seq of lung tumour nodules from KrasG12D/+ (K) and
KrasG12D/+;Lkb1fl/fl (KL) GEMMs (Bai et al., Southern Medical University).

Only GEO-hosted *processed* count matrices are used; no raw FASTQ / SRA
reprocessing is performed and no private data are touched.
"""

from __future__ import annotations

import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path

GEO_SUPPL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180963/suppl/GSE180963_RAW.tar"
)
DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "w200" / "GSE180963"

# GEO sample -> genotype label used throughout the analysis.
SAMPLES = {
    "GSM5481386": ("K", "KrasG12D/+"),
    "GSM5481387": ("KL", "KrasG12D/+;Lkb1fl/fl"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    tar_path = DATA_DIR / "GSE180963_RAW.tar"

    if not tar_path.exists():
        print(f"downloading {GEO_SUPPL}")
        urllib.request.urlretrieve(GEO_SUPPL, tar_path)
    print(f"{tar_path.name}: {tar_path.stat().st_size} bytes sha256={sha256(tar_path)}")

    with tarfile.open(tar_path) as tf:
        tf.extractall(DATA_DIR, filter="data")

    provenance = {"series": "GSE180963", "source": GEO_SUPPL, "samples": {}}
    for gsm, (label, genotype) in SAMPLES.items():
        inner = DATA_DIR / f"{gsm}_{label}.tar.gz"
        with tarfile.open(inner) as tf:
            tf.extractall(DATA_DIR, filter="data")
        provenance["samples"][gsm] = {
            "label": label,
            "genotype": genotype,
            "archive": inner.name,
            "archive_sha256": sha256(inner),
            "files": sorted(p.name for p in (DATA_DIR / label).iterdir()),
        }
        print(f"{gsm} ({label}) -> {DATA_DIR / label}")

    out = DATA_DIR / "provenance.json"
    out.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
