#!/usr/bin/env python3
"""Download GSE302284 Cell Ranger filtered H5 matrices from GEO.

Public data only. Writes files + a provenance JSON (URL, bytes, md5).
Does not commit the H5 blobs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

SERIES = "GSE302284"
RAW_TAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE302nnn/GSE302284/suppl/"
    "GSE302284_RAW.tar"
)
SOFT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE302nnn/GSE302284/soft/"
    "GSE302284_family.soft.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE302nnn/GSE302284/matrix/"
    "GSE302284_series_matrix.txt.gz"
)

SAMPLES = [
    {
        "gsm": "GSM9101264",
        "title": "Lymph Node Sample",
        "library": "N357",
        "model": "patient_357",
        "tissue": "lymph_node",
        "treatment": "neoadjuvant_EGFR_TKI_residual",
        "h5": "GSM9101264_N357_filtered_feature_bc_matrix.h5",
    },
    {
        "gsm": "GSM9101265",
        "title": "Tumor Site Sample",
        "library": "T357",
        "model": "patient_357",
        "tissue": "tumor",
        "treatment": "neoadjuvant_EGFR_TKI_residual",
        "h5": "GSM9101265_T357_filtered_feature_bc_matrix.h5",
    },
    {
        "gsm": "GSM9128162",
        "title": "DFCI282_Ositreated",
        "library": "DFCI282_Ositreated",
        "model": "DFCI282",
        "tissue": "cell_line",
        "treatment": "osimertinib",
        "h5": "GSM9128162_DF282Osi_filtered_feature_bc_matrix.h5",
    },
    {
        "gsm": "GSM9128163",
        "title": "DFCI282_Vehicle",
        "library": "DFCI282_Vehicle",
        "model": "DFCI282",
        "tissue": "cell_line",
        "treatment": "vehicle",
        "h5": "GSM9128163_DF282Veh_filtered_feature_bc_matrix.h5",
    },
    {
        "gsm": "GSM9128164",
        "title": "PC9_Osi",
        "library": "PC9_Osi",
        "model": "PC9",
        "tissue": "cell_line",
        "treatment": "osimertinib",
        "h5": "GSM9128164_PC9Osi_filtered_feature_bc_matrix.h5",
    },
    {
        "gsm": "GSM9128165",
        "title": "PC9_Veh",
        "library": "PC9_Veh",
        "model": "PC9",
        "tissue": "cell_line",
        "treatment": "vehicle",
        "note": (
            "GEO source_name/cell_line fields say PC10; library name is PC9_Veh. "
            "Treated as the paired vehicle for PC9_Osi."
        ),
        "h5": "GSM9128165_PC9Veh_filtered_feature_bc_matrix.h5",
    },
]


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, retries: int = 5) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    last_err: Exception | None = None
    for i in range(retries):
        try:
            print(f"download {url} -> {dest}", flush=True)
            urllib.request.urlretrieve(url, tmp)
            tmp.replace(dest)
            return dest
        except Exception as e:  # noqa: BLE001
            last_err = e
            wait = 4 * (2**i)
            print(f"retry {i + 1}/{retries} after {wait}s: {e}", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"failed to download {url}: {last_err}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cache-dir",
        default="/tmp/gse302284",
        help="Directory for RAW tar + extracted H5 (not committed).",
    )
    p.add_argument(
        "--out-dir",
        default=None,
        help="Where to write provenance JSON (default: results/C_GSE302284).",
    )
    args = p.parse_args()

    cache = Path(args.cache_dir)
    raw = cache / "raw"
    raw.mkdir(parents=True, exist_ok=True)

    root = Path(__file__).resolve().parents[2]
    out = Path(args.out_dir) if args.out_dir else root / "results" / "C_GSE302284"
    out.mkdir(parents=True, exist_ok=True)

    tar_path = download(RAW_TAR_URL, raw / "GSE302284_RAW.tar")
    soft_path = download(SOFT_URL, cache / "GSE302284_family.soft.gz")
    mtx_path = download(MATRIX_URL, cache / "GSE302284_series_matrix.txt.gz")

    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(raw)

    files = []
    for path in (tar_path, soft_path, mtx_path):
        files.append(
            {
                "path": str(path),
                "url": {
                    tar_path: RAW_TAR_URL,
                    soft_path: SOFT_URL,
                    mtx_path: MATRIX_URL,
                }[path],
                "bytes": path.stat().st_size,
                "md5": md5sum(path),
            }
        )
    h5_files = []
    for rec in SAMPLES:
        h5 = raw / rec["h5"]
        if not h5.exists():
            raise FileNotFoundError(h5)
        h5_files.append(
            {
                **rec,
                "path": str(h5),
                "bytes": h5.stat().st_size,
                "md5": md5sum(h5),
            }
        )

    provenance = {
        "series": SERIES,
        "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302284",
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "note": (
            "No sacituzumab / SKB264 / TROP2-ADC treated libraries are deposited. "
            "Analyzable treatment contrast is osimertinib vs vehicle (DFCI282, PC9)."
        ),
        "files": files,
        "samples": h5_files,
    }
    dest = out / "provenance.json"
    dest.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"wrote {dest}", flush=True)


if __name__ == "__main__":
    main()
