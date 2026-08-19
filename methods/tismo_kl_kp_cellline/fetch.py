#!/usr/bin/env python3
"""Download public TISMO + GEO cell-line matrices (no private 8-KL, no FASTQ)."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

ALIYUN_TOKEN = "https://bj21400.api.aliyunfile.com/v2/share_link/get_share_token"
ALIYUN_LIST = "https://bj21400.api.aliyunfile.com/v2/file/list"

TISMO_SHARES = {
    "cell_lines_meta": "YzQB2DQonQE",
    "vitro_sample_meta": "MXMDDtkQLqQ",
    "vivo_sample_meta": "voEA1DXBEFo",
    "vitro_expression": "AuDf6PnXFtV",
    "vivo_expression": "AQKzyRCi1Jb",
}

GEO_FILES = {
    "GSE274352_IFNB.tsv.gz": (
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE274352"
        "&format=file&file=GSE274352_normalizedcounts_genes_IFNB_vs_empty.tsv.gz"
    ),
    "GSE274352_STING.tsv.gz": (
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE274352"
        "&format=file&file=GSE274352_normalizedcounts_genes_STING_vs_emtpy.tsv.gz"
    ),
    "GSE274352_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274352/matrix/"
        "GSE274352_series_matrix.txt.gz"
    ),
    "GSE295685_TPM.txt.gz": (
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE295685"
        "&format=file&file=GSE295685_TPM.txt.gz"
    ),
    "GSE295685_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295685/matrix/"
        "GSE295685_series_matrix.txt.gz"
    ),
    "GSE167381_bRNA-Matrices-Initial-LKB1-Resoration.xlsx": (
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE167381"
        "&format=file&file=GSE167381_bRNA-Matrices-Initial-LKB1-Resoration.xlsx"
    ),
}


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _post_json(url: str, body: dict, headers: dict | None = None) -> dict:
    req = urllib.request.Request(
        url,
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", **(headers or {})},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.loads(resp.read().decode())


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "tismo-kl-kp-cellline/1.0"})
    with urllib.request.urlopen(req, timeout=600) as resp, dest.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)


def fetch_tismo(out: Path) -> list[dict]:
    rows = []
    for key, share_id in TISMO_SHARES.items():
        tok = _post_json(ALIYUN_TOKEN, {"share_id": share_id, "ignoreError": True})
        token = tok.get("share_token")
        if not token:
            raise RuntimeError(f"TISMO share {share_id} ({key}): no token {tok}")
        listing = _post_json(
            ALIYUN_LIST,
            {
                "limit": 100,
                "marker": "",
                "share_id": share_id,
                "parent_file_id": "root",
                "fields": "user_name,dir_size,url,content_type,upload_id,crc64_hash,revision_id,description",
                "url_expire_sec": 7200,
            },
            headers={"x-share-token": token},
        )
        item = (listing.get("items") or [None])[0]
        if not item:
            raise RuntimeError(f"empty share {share_id}")
        dest = out / item["name"]
        if not dest.exists() or dest.stat().st_size != item.get("size"):
            print(f"TISMO {key}: {item['name']} ({item.get('size')} bytes)", flush=True)
            _download(item["download_url"], dest)
        rows.append(
            {
                "key": key,
                "share_id": share_id,
                "name": item["name"],
                "path": str(dest),
                "bytes": dest.stat().st_size,
                "md5": _md5(dest),
                "source": "https://tismo.pku-genomics.org/ (Data Download Aliyun share)",
            }
        )
    return rows


def fetch_geo(out: Path) -> list[dict]:
    rows = []
    for name, url in GEO_FILES.items():
        dest = out / name
        if not dest.exists() or dest.stat().st_size < 1000:
            print(f"GEO {name}", flush=True)
            _download(url, dest)
        rows.append(
            {
                "name": name,
                "path": str(dest),
                "bytes": dest.stat().st_size,
                "md5": _md5(dest),
                "source": url,
            }
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tismo-dir", type=Path, default=Path("/tmp/tismo_dl"))
    ap.add_argument("--geo-dir", type=Path, default=Path("/tmp/geo"))
    ap.add_argument("--manifest", type=Path, default=None)
    args = ap.parse_args()
    manifest = {
        "tismo": fetch_tismo(args.tismo_dir),
        "geo": fetch_geo(args.geo_dir),
        "notes": [
            "TISMO Zeng et al. NAR 2022 PMID 34534350",
            "LLC is labeled LLC, not KL.",
            "No private 8-KL matrices. No GSE76628.",
        ],
    }
    dest = args.manifest or (Path("results/tismo_kl_kp_cellline/tables/download_manifest.json"))
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(manifest, indent=2))
    print("wrote", dest)


if __name__ == "__main__":
    main()
