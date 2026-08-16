#!/usr/bin/env python3
"""Minimal TISMO HTTP client (official PKU mirror).

Gene ICB CSVs come from POST /rtismo/gene/downVivoExprn.
Design metadata comes from POST /tismo/metaData/*.
Bulk matrices live on Aliyun shares (tokens expire; re-resolve each run).
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

BASE = "https://tismo.pku-genomics.org"
ALIYUN = "https://bj21400.api.aliyunfile.com"

ALIYUN_SHARES = {
    "cellline": "YzQB2DQonQE",
    "vitro_ann": "MXMDDtkQLqQ",
    "vivo_ann": "voEA1DXBEFo",
    "vitro_expr": "AuDf6PnXFtV",
    "vivo_expr": "AQKzyRCi1Jb",
    "immune": "kCmAULWGuJh",
}

UA = "Mozilla/5.0 (TISMO-lung-hunt; research reuse)"


def _post_form(url: str, data: dict[str, Any], timeout: int = 180) -> bytes:
    body = urllib.parse.urlencode(data).encode()
    req = urllib.request.Request(url, data=body, method="POST")
    req.add_header("Content-Type", "application/x-www-form-urlencoded")
    req.add_header("User-Agent", UA)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _post_json(url: str, data: dict[str, Any], headers: dict[str, str] | None = None, timeout: int = 60) -> Any:
    req = urllib.request.Request(url, data=json.dumps(data).encode(), method="POST")
    req.add_header("Content-Type", "application/json")
    req.add_header("User-Agent", UA)
    if headers:
        for k, v in headers.items():
            req.add_header(k, v)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def fetch_json(path: str) -> Any:
    raw = _post_form(BASE + path, {})
    return json.loads(raw.decode())


def download_vivo_gene_csv(gene: str, dest: Path) -> Path:
    """Download the official ICB-comparison CSV for one gene.

    Always request tumorList=['All']. Filtering to LLC on the server
    silently drops SRX8918393 (Setdb1-KO ICB, Tacstd2=5.41).
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    raw = _post_form(
        BASE + "/rtismo/gene/downVivoExprn",
        {
            "filename": "genetreatment_vivo.csv",
            "type": "3",
            "gene": gene,
            "icbList": '["All"]',
            "tumorList": '["All"]',
        },
    )
    dest.write_bytes(raw)
    return dest


def aliyun_download(share_key: str, dest_dir: Path) -> Path:
    sid = ALIYUN_SHARES[share_key]
    tok = _post_json(
        f"{ALIYUN}/v2/share_link/get_share_token",
        {"share_id": sid, "ignoreError": True},
    )["share_token"]
    listing = _post_json(
        f"{ALIYUN}/v2/file/list",
        {"share_id": sid, "parent_file_id": "root", "limit": 20},
        headers={"x-share-token": tok},
    )
    item = listing["items"][0]
    info = _post_json(
        f"{ALIYUN}/v2/file/get_download_url",
        {"share_id": sid, "file_id": item["file_id"], "file_name": item["name"]},
        headers={"x-share-token": tok},
    )
    url = info["url"]
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / item["name"]
    req = urllib.request.Request(url, method="GET")
    req.add_header("User-Agent", UA)
    with urllib.request.urlopen(req, timeout=600) as resp:
        dest.write_bytes(resp.read())
    return dest
