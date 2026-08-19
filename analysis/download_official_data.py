#!/usr/bin/env python3
"""Download official CosMx NSCLC flat files and GSE307534 invasive LUAD Visium."""

from __future__ import annotations

import argparse
import io
import json
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

COSMX_SAMPLES = [
    "Lung5_Rep1",
    "Lung5_Rep2",
    "Lung5_Rep3",
    "Lung6",
    "Lung9_Rep1",
    "Lung9_Rep2",
    "Lung12",
    "Lung13",
]

# Author-labeled invasive LUAD only (exclude AAH / AIS / MIA / Normal).
# Accessions and library names from GEO GSE307534 sample titles.
VISIUM_LUAD = [
    ("GSM9226169", "P1_LUAD"),
    ("GSM9226171", "P2_LUAD"),
    ("GSM9226173", "P3_LUAD"),
    ("GSM9226177", "P4_LUAD"),
    ("GSM9226179", "P5_LUAD"),
    ("GSM9226181", "P6_LUAD"),
    ("GSM9226182", "P7_LUAD"),
    ("GSM9226183", "P7_LUAD-1"),
    ("GSM9226185", "P8_LUAD"),
    ("GSM9226188", "P9_LUAD"),
    ("GSM9226190", "P10_LUAD"),
    ("GSM9226192", "P11_LUAD"),
    ("GSM9226194", "P12_LUAD"),
    ("GSM9226196", "P13_LUAD"),
    ("GSM9226198", "P14_LUAD"),
    ("GSM9226200", "P15_LUAD"),
    ("GSM9226202", "P16_LUAD"),
    ("GSM9226204", "P17_LUAD"),
    ("GSM9226206", "P18_LUAD"),
    ("GSM9226208", "P19_LUAD"),
    ("GSM9226210", "P20_LUAD"),
    ("GSM9226213", "P21_LUAD"),
    ("GSM9226216", "P22_LUAD"),
    ("GSM9226219", "P23_LUAD"),
    ("GSM9226221", "P24_LUAD"),
    ("GSM9226223", "P25_LUAD"),
]

KEEP_COSMX = ("_metadata_file.csv", "_exprMat_file.csv", "_fov_positions_file.csv")
KEEP_VISIUM = (
    "filtered_feature_bc_matrix.h5",
    "filtered_feature_bc_matrix/",
    "matrix.mtx",
    "barcodes.tsv",
    "features.tsv",
    "genes.tsv",
    "tissue_positions",
    "scalefactors_json.json",
    "tissue_hires_image.png",
    "tissue_lowres_image.png",
)


def _want_cosmx(name: str) -> bool:
    lower = name.lower()
    if lower.endswith("/"):
        return False
    return any(name.endswith(s) for s in KEEP_COSMX)


def _want_visium(name: str) -> bool:
    lower = name.lower()
    if lower.endswith("/"):
        return False
    if any(tok in lower for tok in KEEP_VISIUM):
        # skip huge hi-res images; keep lowres only if needed later
        if "tissue_hires_image" in lower:
            return False
        if lower.endswith((".tif", ".tiff", ".jpg", ".jpeg", ".cloupe", ".html")):
            return False
        return True
    return False


def _extract_filtered(tar: tarfile.TarFile, dest: Path, keep_fn) -> list[str]:
    dest.mkdir(parents=True, exist_ok=True)
    written = []
    for member in tar:
        if not member.isfile():
            continue
        if not keep_fn(member.name):
            continue
        # flatten to dest / basename, but keep unique names
        out = dest / Path(member.name).name
        if out.exists() and out.stat().st_size > 0:
            written.append(str(out))
            continue
        src = tar.extractfile(member)
        if src is None:
            continue
        out.write_bytes(src.read())
        written.append(str(out))
    return written


def download_cosmx(out_root: Path, samples: list[str] | None = None) -> dict:
    samples = samples or COSMX_SAMPLES
    manifest = {}
    for sample in samples:
        dest = out_root / "cosmx" / sample
        dest.mkdir(parents=True, exist_ok=True)
        marker = dest / "_DONE.json"
        if marker.exists() and (dest / f"{sample}_metadata_file.csv").exists():
            print(f"[cosmx] skip {sample} (already present)")
            manifest[sample] = json.loads(marker.read_text())
            continue
        key = f"SMI-Compressed/{sample}/{sample} SMI Flat data.tar.gz"
        url = "https://nanostring-public-share.s3.us-west-2.amazonaws.com/" + urllib.request.quote(
            key, safe="/"
        )
        print(f"[cosmx] streaming {sample} from official S3 …")
        req = urllib.request.Request(url, headers={"User-Agent": "cldn4-tumor-geography/1.0"})
        with urllib.request.urlopen(req, timeout=600) as resp:
            # tarfile needs a seekable stream for some encodings; buffer to a temp file
            tmp = dest / f"{sample}.tar.gz"
            with tmp.open("wb") as fh:
                while True:
                    chunk = resp.read(1024 * 1024 * 8)
                    if not chunk:
                        break
                    fh.write(chunk)
        print(f"[cosmx] extracted archive {tmp.stat().st_size/1e6:.1f} MB; pulling CSVs")
        with tarfile.open(tmp, "r:gz") as tar:
            written = _extract_filtered(tar, dest, _want_cosmx)
        tmp.unlink(missing_ok=True)
        info = {"sample": sample, "url": url, "files": written}
        marker.write_text(json.dumps(info, indent=2))
        manifest[sample] = info
        print(f"[cosmx] {sample}: {len(written)} files")
    return manifest


def _visium_urls(gsm: str, lib: str) -> list[str]:
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/"
    return [
        base + f"{gsm}_{lib}.tar.gz",
        base + f"{gsm}_{lib}.tar",
    ]


def download_visium(out_root: Path, samples: list[tuple[str, str]] | None = None) -> dict:
    samples = samples or VISIUM_LUAD
    manifest = {}
    for gsm, lib in samples:
        dest = out_root / "visium" / f"{gsm}_{lib}"
        dest.mkdir(parents=True, exist_ok=True)
        marker = dest / "_DONE.json"
        if marker.exists() and any(dest.iterdir()):
            print(f"[visium] skip {gsm}_{lib} (already present)")
            manifest[f"{gsm}_{lib}"] = json.loads(marker.read_text())
            continue
        last_err = None
        written = []
        used_url = None
        for url in _visium_urls(gsm, lib):
            print(f"[visium] trying {url}")
            req = urllib.request.Request(url, headers={"User-Agent": "cldn4-tumor-geography/1.0"})
            try:
                with urllib.request.urlopen(req, timeout=600) as resp:
                    data = resp.read()
                used_url = url
                bio = io.BytesIO(data)
                with tarfile.open(fileobj=bio, mode="r:*") as tar:
                    written = _extract_filtered(tar, dest, _want_visium)
                    # keep nested spatial folders if extract flattened names collide
                    if not written:
                        for member in tar:
                            if member.isfile() and _want_visium(member.name):
                                out = dest / Path(member.name).name
                                src = tar.extractfile(member)
                                if src is None:
                                    continue
                                out.write_bytes(src.read())
                                written.append(str(out))
                break
            except urllib.error.HTTPError as e:
                last_err = e
                print(f"  HTTP {e.code}")
                continue
        if not written:
            print(f"[visium] FAILED {gsm}_{lib}: {last_err}")
            continue
        info = {"gsm": gsm, "library": lib, "url": used_url, "files": written}
        marker.write_text(json.dumps(info, indent=2))
        manifest[f"{gsm}_{lib}"] = info
        print(f"[visium] {gsm}_{lib}: {len(written)} files")
    return manifest


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="/tmp/spatial_data")
    p.add_argument("--cosmx-only", action="store_true")
    p.add_argument("--visium-only", action="store_true")
    p.add_argument("--cosmx-samples", nargs="*", default=None)
    p.add_argument("--max-visium", type=int, default=None)
    args = p.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if not args.visium_only:
        download_cosmx(out, args.cosmx_samples)
    if not args.cosmx_only:
        vis = VISIUM_LUAD[: args.max_visium] if args.max_visium else VISIUM_LUAD
        download_visium(out, vis)


if __name__ == "__main__":
    main()
