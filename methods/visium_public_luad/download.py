#!/usr/bin/env python3
"""Download public Visium LUAD matrices (GSE277206 + Zenodo 13337961).

Raw archives stay under data/visium_public_luad/ (gitignored).
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "visium_public_luad"

GSE_TAR = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277206/suppl/GSE277206_RAW.tar"
)
ZENODO_FILES = {
    "lepidic-filtered_feature_bc_matrix.rar": (
        "https://zenodo.org/api/records/13337961/files/lepidic-filtered_feature_bc_matrix.rar/content"
    ),
    "solid-filtered_feature_bc_matrix.rar": (
        "https://zenodo.org/api/records/13337961/files/solid-filtered_feature_bc_matrix.rar/content"
    ),
}


def _fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print("have", dest, flush=True)
        return dest
    print("get", url, "->", dest, flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    return dest


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _unrar(archive: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    marker = dest / "filtered_feature_bc_matrix" / "matrix.mtx.gz"
    if marker.exists():
        return
    unrar = shutil.which("unrar")
    seven = shutil.which("7z")
    if unrar:
        subprocess.check_call([unrar, "x", "-o+", str(archive), str(dest) + "/"])
    elif seven:
        subprocess.check_call([seven, "x", f"-o{dest}", "-y", str(archive)])
    else:
        raise RuntimeError("need unrar or 7z to extract Zenodo .rar archives")


def download_gse277206() -> dict:
    d = DATA / "GSE277206"
    tar_path = _fetch(GSE_TAR, d / "GSE277206_RAW.tar")
    with tarfile.open(tar_path) as t:
        t.extractall(d)
    h5s = sorted(d.glob("GSM*_filtered_feature_bc_matrix.h5"))
    rec = {
        "source": "GSE277206",
        "tar_sha256": sha256(tar_path),
        "files": [{"name": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)} for p in h5s],
    }
    return rec


def download_zenodo() -> dict:
    d = DATA / "zenodo_13337961"
    rec = {"source": "zenodo_13337961", "files": []}
    for name, url in ZENODO_FILES.items():
        rar = _fetch(url, d / name)
        label = "lepidic" if name.startswith("lepidic") else "solid"
        _unrar(rar, d / label)
        rec["files"].append(
            {
                "name": name,
                "bytes": rar.stat().st_size,
                "sha256": sha256(rar),
                "extracted": str((d / label / "filtered_feature_bc_matrix").relative_to(DATA)),
            }
        )
    return rec


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    gse = download_gse277206()
    zen = download_zenodo()
    print("GSE277206", gse)
    print("Zenodo", zen)


if __name__ == "__main__":
    main()
