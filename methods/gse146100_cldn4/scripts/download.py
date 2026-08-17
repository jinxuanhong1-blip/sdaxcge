#!/usr/bin/env python3
"""Download public GSE146100 objects (GEO NormData + TISCH2 h5).

Both files are <2 GB. Nothing is committed.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

FILES = {
    "GSE146100_NormData.txt.gz": {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146100/"
            "suppl/GSE146100_NormData.txt.gz"
        ),
        "expected_bytes": 229_465_901,
    },
    "NSCLC_GSE146100_expression.h5": {
        "url": (
            "https://tisch.compbio.cn/static/data/NSCLC_GSE146100/"
            "NSCLC_GSE146100_expression.h5"
        ),
        "expected_bytes": 200_745_710,
        "max_bytes": 2_000_000_000,
    },
    "NSCLC_GSE146100_CellMetainfo_table.tsv": {
        "url": (
            "https://tisch.compbio.cn/static/data/NSCLC_GSE146100/"
            "NSCLC_GSE146100_CellMetainfo_table.tsv"
        ),
        "expected_bytes": None,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def curl(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    cmd = ["curl", "-fL", "--retry", "3", "--retry-delay", "4", "-o", str(tmp), url]
    subprocess.run(cmd, check=True)
    tmp.replace(dest)


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    provenance: dict = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "files": {},
    }
    for name, spec in FILES.items():
        dest = DATA / name
        url = spec["url"]
        expected = spec.get("expected_bytes")
        if dest.exists() and dest.stat().st_size > 0:
            if expected is None or dest.stat().st_size == expected:
                print(f"already present: {dest} ({dest.stat().st_size} bytes)")
                provenance["files"][name] = {
                    "url": url,
                    "path": str(dest),
                    "bytes": dest.stat().st_size,
                    "sha256": sha256(dest),
                    "skipped": True,
                }
                continue
        print(f"downloading {url}")
        curl(url, dest)
        size = dest.stat().st_size
        print(f"wrote {dest} ({size} bytes)")
        if spec.get("max_bytes") and size > spec["max_bytes"]:
            dest.unlink()
            raise SystemExit(f"{name} is {size} bytes (>2 GB); not used")
        if expected is not None and size != expected:
            raise SystemExit(f"unexpected size {size} != {expected} for {name}")
        provenance["files"][name] = {
            "url": url,
            "path": str(dest),
            "bytes": size,
            "sha256": sha256(dest),
            "skipped": False,
        }
    (DATA / "download_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
