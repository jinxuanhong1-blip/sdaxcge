#!/usr/bin/env python3
"""Download processed inputs for B6 radius rework. Raw FASTQ/SRA skipped.

Writes a SHA-256 manifest. Local cache is NOT committed (default: /tmp/b6_radius_data).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from urllib.request import Request, urlopen

EMTAB_FTP = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/530/E-MTAB-13530/Files"
TUMOR_SECTIONS = [
    "P10_T1", "P10_T2", "P10_T3", "P10_T4",
    "P11_T1", "P11_T2", "P11_T3", "P11_T4",
    "P15_T1", "P15_T2",
    "P16_T1", "P16_T2",
    "P17_T1", "P17_T2",
    "P19_T1", "P19_T2",
    "P24_T1", "P24_T2",
    "P25_T1", "P25_T2",
]
MAX_BYTES = 2_000_000_000

FILES = [
    (
        "gse271689/GSE271689_series_matrix.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/matrix/GSE271689_series_matrix.txt.gz",
    ),
    (
        "gse271689/GSE271689_RAW.tar",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/suppl/GSE271689_RAW.tar",
    ),
    (
        "geomx_pkc/Hs_R_NGS_WTA_v1.0.pkc",
        "https://zenodo.org/records/12752405/files/Hs_R_NGS_WTA_v1.0.pkc?download=1",
    ),
    (
        "gse271689/41588_2025_2351_MOESM10_ESM.xlsx",
        "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-025-02351-7/MediaObjects/41588_2025_2351_MOESM10_ESM.xlsx",
    ),
    (
        "emtab13530/E-MTAB-13530.sdrf.txt",
        f"{EMTAB_FTP}/E-MTAB-13530.sdrf.txt",
    ),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def head_length(url: str) -> int | None:
    try:
        with urlopen(Request(url, method="HEAD"), timeout=60) as r:
            cl = r.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:
        return None


def curl(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return {
            "url": url,
            "path": str(dest),
            "bytes": dest.stat().st_size,
            "sha256": sha256(dest),
            "decision": "accepted",
            "note": "already present",
        }
    size = head_length(url)
    if size is not None and size > MAX_BYTES:
        return {
            "url": url,
            "path": str(dest),
            "bytes": size,
            "sha256": "",
            "decision": "refused",
            "refusal_reason": f">2GB processed-data budget ({size} bytes)",
        }
    cmd = ["curl", "-fL", "--retry", "4", "--retry-delay", "4", "-o", str(dest), url]
    print("GET", url, flush=True)
    subprocess.check_call(cmd)
    return {
        "url": url,
        "path": str(dest),
        "bytes": dest.stat().st_size,
        "sha256": sha256(dest),
        "decision": "accepted",
        "note": "downloaded",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/tmp/b6_radius_data")
    ap.add_argument("--manifest", default=None)
    args = ap.parse_args()
    data = Path(args.data_dir)
    data.mkdir(parents=True, exist_ok=True)

    rows = []
    for rel, url in FILES:
        rows.append({"file_name": rel, **curl(url, data / rel)})

    for sec in TUMOR_SECTIONS:
        for suffix in ("-filtered_feature_bc_matrix.h5", "-spatial.tar"):
            name = f"{sec}{suffix}"
            rel = f"emtab13530/{name}"
            rows.append({"file_name": rel, **curl(f"{EMTAB_FTP}/{name}", data / rel)})

    here = Path(__file__).resolve().parent
    manifest_path = Path(args.manifest) if args.manifest else here / "download_manifest.json"
    payload = {
        "task": "B6_radius",
        "note": "Processed files only. Tumor Visium sections + GeoMx WTA. FASTQ/SRA not downloaded.",
        "data_dir": str(data),
        "n_files": len(rows),
        "n_accepted": sum(1 for r in rows if r["decision"] == "accepted"),
        "n_refused": sum(1 for r in rows if r["decision"] == "refused"),
        "total_accepted_bytes": sum(r["bytes"] for r in rows if r["decision"] == "accepted"),
        "files": rows,
    }
    manifest_path.write_text(json.dumps(payload, indent=2))
    print(f"wrote {manifest_path} accepted={payload['n_accepted']} bytes={payload['total_accepted_bytes']}", flush=True)
    return 0 if payload["n_refused"] == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
