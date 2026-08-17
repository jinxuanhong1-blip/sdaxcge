#!/usr/bin/env python3
"""Download public processed GSE123902 (Laughney 2020).

GSE123902_RAW.tar is 90.4 MB. The author annotated H5 (36.5 GB) is skipped.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

FILES = {
    "GSE123902_RAW.tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar"
    ),
    "GSE123902_GEO_README.rtf": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_GEO_README.rtf"
    ),
}

SKIPPED = {
    "PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5": {
        "url": (
            "https://s3.amazonaws.com/dp-lab-data-public/"
            "lung-development-cancer-progression/PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5"
        ),
        "bytes": 36489164157,
        "reason": "author annotated H5 is 36.5 GB (>2 GB public-processed gate)",
    }
}


def curl_download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    cmd = ["curl", "-fL", "--retry", "4", "--retry-delay", "4", "-o", str(tmp), url]
    print(f"GET {url}", flush=True)
    subprocess.run(cmd, check=True)
    tmp.replace(dest)
    print(f"wrote {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("/tmp/gse123902_slingshot"))
    args = p.parse_args()
    for name, url in FILES.items():
        curl_download(url, args.out / name)
    print("skipped:")
    for name, rec in SKIPPED.items():
        print(f"  {name}: {rec['reason']}")


if __name__ == "__main__":
    main()
