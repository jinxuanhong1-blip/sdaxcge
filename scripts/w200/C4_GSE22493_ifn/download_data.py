#!/usr/bin/env python3
"""Download public GSE22493 files for the C4 IFN/MHC-I/APM slice.

All files are open GEO FTP objects. Nothing here is controlled-access.
Override the cache directory with W200_C4_GSE22493_DATA.
"""

from __future__ import annotations

import hashlib
import os
import tarfile
import urllib.request

DATA_DIR = os.environ.get("W200_C4_GSE22493_DATA", "/tmp/w200_c4_gse22493_ifn")
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493"

FILES = {
    "GSE22493_series_matrix.txt.gz": f"{BASE}/matrix/GSE22493_series_matrix.txt.gz",
    "GSE22493_family.soft.gz": f"{BASE}/soft/GSE22493_family.soft.gz",
    "GSE22493_RAW.tar": f"{BASE}/suppl/GSE22493_RAW.tar",
}


def md5_of(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: str) -> None:
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"exists  {dest}  md5={md5_of(dest)}  bytes={os.path.getsize(dest)}")
        return
    print(f"fetch   {url}")
    tmp = dest + ".part"
    urllib.request.urlretrieve(url, tmp)
    os.replace(tmp, dest)
    print(f"wrote   {dest}  md5={md5_of(dest)}  bytes={os.path.getsize(dest)}")


def extract_raw(tar_path: str, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    with tarfile.open(tar_path, "r") as tf:
        for member in tf.getmembers():
            name = os.path.basename(member.name)
            if not name.endswith(".txt.gz"):
                continue
            dest = os.path.join(out_dir, name)
            if os.path.exists(dest) and os.path.getsize(dest) > 0:
                print(f"exists  {dest}  md5={md5_of(dest)}")
                continue
            src = tf.extractfile(member)
            if src is None:
                continue
            with open(dest, "wb") as fh:
                fh.write(src.read())
            print(f"extract {dest}  md5={md5_of(dest)}")


def main() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"DATA_DIR={DATA_DIR}")
    for name, url in FILES.items():
        fetch(url, os.path.join(DATA_DIR, name))
    extract_raw(
        os.path.join(DATA_DIR, "GSE22493_RAW.tar"),
        os.path.join(DATA_DIR, "scanarray"),
    )
    print("done")


if __name__ == "__main__":
    main()
