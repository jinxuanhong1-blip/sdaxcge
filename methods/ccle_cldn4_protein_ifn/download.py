#!/usr/bin/env python3
"""Download public CCLE/DepMap protein matrices + 24Q4 Model.csv.

Additive extra to methods/depmap_cldn4_ifn (RNA ρ already reported).
This download is protein only: Gygi/Nusinow CCLE TMT-MS, CCLE RPPA,
ProCan-DepMapSanger (inventory), DepMap 24Q4 Model.csv, MSigDB Hallmark GMT.

Portal is bot-gated; URLs are the same public files used in
scripts/rework/B2_protein_followup and scripts/rework/DepMap_IFN_treated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

UA = "sdaxcge-ccle-cldn4-protein-ifn/1.0 (public extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"

FILES = {
    "Model.csv": {
        "url": "https://ndownloader.figshare.com/files/51065297",
        "min_bytes": 100_000,
        "note": "DepMap Public 24Q4 Model.csv (Figshare+ 10.25452/figshare.plus.27993248.v1)",
    },
    "protein_quant_current_normalized.csv.gz": {
        "url": "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz",
        "min_bytes": 5_000_000,
        "note": "Nusinow et al. Cell 2020 CCLE TMT-MS (Gygi lab current normalized)",
        "fallbacks": [
            "https://gygi.med.harvard.edu/sites/gygi.med.harvard.edu/files/documents/protein_quant_current_normalized.csv.gz",
        ],
    },
    "CCLE_RPPA_20180123.csv": {
        "url": "https://data.broadinstitute.org/ccle/CCLE_RPPA_20180123.csv",
        "min_bytes": 100_000,
        "note": "CCLE RPPA 20180123",
    },
    "CCLE_RPPA_Ab_info_20180123.csv": {
        "url": "https://data.broadinstitute.org/ccle/CCLE_RPPA_Ab_info_20180123.csv",
        "min_bytes": 1_000,
        "note": "CCLE RPPA antibody map 20180123",
    },
    "ProCan-DepMapSanger_protein_matrix_8498_averaged.txt": {
        "url": "https://ndownloader.figshare.com/files/34411172",
        "min_bytes": 1_000_000,
        "note": "ProCan-DepMapSanger 8498-protein averaged matrix (Gonçalves 2022)",
    },
    "ProCan-DepMapSanger_mapping_file_averaged.txt": {
        "url": "https://ndownloader.figshare.com/files/34411133",
        "min_bytes": 1_000,
        "note": "ProCan-DepMapSanger sample mapping",
    },
    "h.all.v2024.1.Hs.symbols.gmt": {
        "url": (
            "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
            "h.all.v2024.1.Hs.symbols.gmt"
        ),
        "min_bytes": 5_000,
        "note": "MSigDB 2024.1 Hs Hallmark GMT",
    },
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"GET {url} -> {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)


def fetch(spec: dict, dest: Path) -> str:
    urls = [spec["url"]] + spec.get("fallbacks", [])
    last_err: Exception | None = None
    for url in urls:
        try:
            download(url, dest)
            if dest.stat().st_size < spec["min_bytes"]:
                raise RuntimeError(f"{dest} too small ({dest.stat().st_size} bytes)")
            return url
        except Exception as e:
            last_err = e
            print(f"WARN {url}: {e}", flush=True)
            if dest.exists():
                dest.unlink()
    raise SystemExit(f"failed to download {dest.name}: {last_err}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/ccle_cldn4_protein_ifn")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = {"files": {}}
    for name, spec in FILES.items():
        dest = outdir / name
        used_url = spec["url"]
        if dest.exists() and dest.stat().st_size >= spec["min_bytes"]:
            print(f"OK exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        else:
            used_url = fetch(spec, dest)
        manifest["files"][name] = {
            "url": used_url,
            "note": spec["note"],
            "bytes": dest.stat().st_size,
            "sha256": sha256_file(dest),
        }

    (outdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: v["bytes"] for k, v in manifest["files"].items()}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
