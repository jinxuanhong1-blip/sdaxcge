#!/usr/bin/env python3
"""Download the public OncoSG LUAD expression and sample metadata."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import urllib.request
from pathlib import Path


DATA_DIR = Path(os.environ.get("ONCOSG_DATA_DIR", "/tmp/hunt_oncosg_data"))
DATASETS = {
    "expression_zscores.tsv": {
        "url": (
            "https://media.githubusercontent.com/media/cBioPortal/datahub/master/"
            "public/luad_oncosg_2020/"
            "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
        ),
        "sha256": "44093dbc6633e5280d0fcd18fb07f5c6ee5c2cff175cd76b48b3665c3d4ff54f",
    },
    "clinical_sample.tsv": {
        "url": (
            "https://media.githubusercontent.com/media/cBioPortal/datahub/master/"
            "public/luad_oncosg_2020/data_clinical_sample.txt"
        ),
        "sha256": "55739b69bf4f2e8e624aab51a0b5903ea797751905219b7bf6581afb4c96f368",
    },
}
API_SNAPSHOTS = {
    "study.json": "https://www.cbioportal.org/api/studies/luad_oncosg_2020",
    "molecular_profiles.json": (
        "https://www.cbioportal.org/api/studies/luad_oncosg_2020/"
        "molecular-profiles?projection=DETAILED"
    ),
    "clinical_attributes.json": (
        "https://www.cbioportal.org/api/studies/luad_oncosg_2020/"
        "clinical-attributes"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, destination: Path) -> None:
    request = urllib.request.Request(
        url, headers={"User-Agent": "hunt-oncosg-repro/1.0"}
    )
    temporary = destination.with_suffix(destination.suffix + ".part")
    with urllib.request.urlopen(request, timeout=120) as response:
        with temporary.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
    temporary.replace(destination)


def main() -> int:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for filename, metadata in DATASETS.items():
        path = DATA_DIR / filename
        if not path.exists() or sha256(path) != metadata["sha256"]:
            fetch(metadata["url"], path)
        observed = sha256(path)
        if observed != metadata["sha256"]:
            print(
                f"Checksum mismatch for {path}: {observed} != {metadata['sha256']}",
                file=sys.stderr,
            )
            return 1
        print(f"{filename}\t{path.stat().st_size}\t{observed}")

    for filename, url in API_SNAPSHOTS.items():
        path = DATA_DIR / filename
        fetch(url, path)
        # Fail early if a proxy returned HTML or malformed JSON.
        with path.open(encoding="utf-8") as handle:
            json.load(handle)
        print(f"{filename}\t{path.stat().st_size}\t{sha256(path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
