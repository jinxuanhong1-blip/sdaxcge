#!/usr/bin/env python3
"""Download IMvigor210CoreBiologies 1.0.0 and dump counts + clinical tables.

Official public source (Mariathasan et al., Nature 2018; CC BY 3.0):
  http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz

Requires R (used only to extract the S4 CountDataSet). Processed matrices
are written under data/IMvigor210/ and are gitignored.

Usage:
  python3 scripts/w200/IMvigor_both/download.py
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = ROOT / "data" / "IMvigor210"
SCRIPT_DIR = Path(__file__).resolve().parent

PKG_URL = (
    "http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/"
    "IMvigor210CoreBiologies_1.0.0.tar.gz"
)
PKG_NAME = "IMvigor210CoreBiologies_1.0.0.tar.gz"
# Published package size from research-pub.gene.com (2019-02-12).
EXPECTED_SIZE = 122_127_298


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_existing_tarball() -> Path | None:
    candidates = [
        DATA_DIR / PKG_NAME,
        Path("/tmp") / PKG_NAME,
        Path("/tmp/IMvigor210CoreBiologies_1.0.0.tar.gz"),
    ]
    for p in candidates:
        if p.is_file() and p.stat().st_size == EXPECTED_SIZE:
            return p
    return None


def download_tarball(dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    existing = find_existing_tarball()
    if existing is not None:
        if existing.resolve() != dest.resolve():
            shutil.copy2(existing, dest)
        return dest
    print(f"Downloading {PKG_URL}")
    urllib.request.urlretrieve(PKG_URL, dest)
    if dest.stat().st_size != EXPECTED_SIZE:
        raise RuntimeError(
            f"Unexpected download size {dest.stat().st_size} "
            f"(expected {EXPECTED_SIZE})"
        )
    return dest


def extract_cds_rdata(tarball: Path, dest_rdata: Path) -> Path:
    dest_rdata.parent.mkdir(parents=True, exist_ok=True)
    if dest_rdata.is_file() and dest_rdata.stat().st_size > 0:
        return dest_rdata
    member = "IMvigor210CoreBiologies/data/cds.RData"
    with tarfile.open(tarball, "r:gz") as tf:
        src = tf.extractfile(member)
        if src is None:
            raise RuntimeError(f"{member} missing from tarball")
        dest_rdata.write_bytes(src.read())
    return dest_rdata


def dump_tables(rdata: Path, out_dir: Path) -> None:
    needed = [out_dir / "pData.csv", out_dir / "fData.csv", out_dir / "counts.csv.gz"]
    if all(p.is_file() and p.stat().st_size > 0 for p in needed):
        print("Extracted tables already present")
        return
    rscript = SCRIPT_DIR / "extract_cds.R"
    cmd = ["Rscript", str(rscript), str(rdata), str(out_dir)]
    print("Running", " ".join(cmd))
    subprocess.check_call(cmd)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    tarball = download_tarball(DATA_DIR / PKG_NAME)
    digest = sha256(tarball)
    (DATA_DIR / "tarball.sha256").write_text(f"{digest}  {PKG_NAME}\n")
    print(f"tarball sha256 {digest}")
    rdata = extract_cds_rdata(tarball, DATA_DIR / "cds.RData")
    dump_tables(rdata, DATA_DIR)
    print(f"Ready: {DATA_DIR}")


if __name__ == "__main__":
    main()
