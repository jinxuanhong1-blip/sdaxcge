#!/usr/bin/env python3
"""Fetch public TISMO processed tables for lung / LLC Cldn4 work.

Sources (no private 8 KL, no FASTQ):
  - https://tismo.pku-genomics.org/  (Zeng et al., NAR 2022, PMID 34534350)
  - Gene-module per-sample CSV: /rtismo/gene/downVivoExprn (type=3)
  - Data-download Aliyun shares: vivo sample annotations + immune RDS
"""

from __future__ import annotations

import argparse
import json
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

PORT = "https://tismo.pku-genomics.org/tismo"
RPORT = "https://tismo.pku-genomics.org/rtismo"
ALIYUN_TOKEN = "https://bj21400.api.aliyunfile.com/v2/share_link/get_share_token"
ALIYUN_LIST = "https://bj21400.api.aliyunfile.com/v2/file/list"

ICB_TREATMENTS = [
    "antiCTLA4",
    "antiCTLA4&antiPD1",
    "antiCTLA4&antiPDL1",
    "antiPD1",
    "antiPDL1",
    "antiPDL2",
]

# Cldn4-only claim. T/NK and IFN/MHC are the correlates, not Tacstd2.
GENES = [
    "Cldn4",
    "Cd8a",
    "Cd3e",
    "Cd3d",
    "Cd2",
    "Nkg7",
    "Gzmb",
    "Prf1",
    "Ifng",
    "Ncr1",
    "Stat1",
    "Cxcl9",
    "Cxcl10",
    "Ido1",
    "H2-Aa",
    "H2-K1",
    "H2-D1",
    "H2-Q4",
    "B2m",
    "Tap1",
    "Tap2",
    "Cd274",
    "Actb",
    "Epcam",
    "Ptprc",
]

ALIYUN_SHARES = {
    "vivo_meta": "voEA1DXBEFo",
    "immune": "kCmAULWGuJh",
    "cell_lines": "YzQB2DQonQE",
}

CTX = ssl.create_default_context()


def _request(url: str, data: bytes | None, headers: dict[str, str], timeout: int) -> bytes:
    req = urllib.request.Request(url, data=data, headers=headers, method="POST" if data is not None else "GET")
    last_err: Exception | None = None
    for attempt in range(1, 5):
        try:
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError) as err:
            last_err = err
            time.sleep(2**attempt)
    raise RuntimeError(f"request {url} failed: {last_err}")


def post_form(url: str, data: dict[str, str], timeout: int = 180) -> bytes:
    payload = urllib.parse.urlencode(data).encode()
    return _request(url, payload, {"Content-Type": "application/x-www-form-urlencoded"}, timeout)


def post_json(url: str, body: dict, headers: dict[str, str] | None = None, timeout: int = 60) -> dict:
    hdrs = {"Content-Type": "application/json"}
    if headers:
        hdrs.update(headers)
    blob = _request(url, json.dumps(body).encode(), hdrs, timeout)
    return json.loads(blob.decode())


def fetch_vivo_meta(raw_dir: Path) -> None:
    blob = _request(
        f"{PORT}/metaData/vivoMeta",
        json.dumps({}).encode(),
        {"Content-Type": "application/json"},
        180,
    )
    payload = json.loads(blob)
    if payload.get("status") != 600200:
        raise RuntimeError(f"vivoMeta unexpected status: {payload.get('status')}")
    dest = raw_dir / "vivo_meta.json"
    dest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {dest}  n={len(payload['data'])}")


def fetch_gene_icb(raw_dir: Path, gene: str, tumor: str, dest_name: str) -> None:
    blob = post_form(
        f"{RPORT}/gene/downVivoExprn",
        {
            "filename": dest_name,
            "type": "3",
            "gene": gene,
            "icbList": json.dumps(ICB_TREATMENTS),
            "tumorList": json.dumps([tumor]),
        },
    )
    if blob[:1] not in (b"S", b"s") and b"Samples" not in blob[:120]:
        raise RuntimeError(f"{gene}/{tumor} download did not look like a CSV ({blob[:120]!r})")
    dest = raw_dir / dest_name
    dest.write_bytes(blob)
    print(f"wrote {dest}  rows≈{blob.count(b\"\\n\") - 1}")


def fetch_aliyun_share(raw_dir: Path, share_id: str) -> Path:
    tok = post_json(ALIYUN_TOKEN, {"share_id": share_id, "ignoreError": True})
    token = tok.get("share_token")
    if not token:
        raise RuntimeError(f"no share_token for {share_id}: {tok}")
    listing = post_json(
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
    items = listing.get("items") or []
    if not items:
        raise RuntimeError(f"empty Aliyun share {share_id}")
    item = items[0]
    url = item["download_url"]
    name = item["name"]
    dest = raw_dir / name
    print(f"downloading {name} ({item.get('size')} bytes) …")
    req = urllib.request.Request(url, method="GET")
    last_err: Exception | None = None
    for attempt in range(1, 5):
        try:
            with urllib.request.urlopen(req, timeout=300, context=CTX) as resp:
                dest.write_bytes(resp.read())
            print(f"wrote {dest}")
            return dest
        except (urllib.error.URLError, TimeoutError) as err:
            last_err = err
            time.sleep(2**attempt)
    raise RuntimeError(f"Aliyun download failed: {last_err}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("methods/tismo_llc_cldn4/raw"))
    args = parser.parse_args()
    raw_dir = args.out
    raw_dir.mkdir(parents=True, exist_ok=True)

    fetch_vivo_meta(raw_dir)
    for key, sid in ALIYUN_SHARES.items():
        try:
            fetch_aliyun_share(raw_dir, sid)
        except Exception as err:
            print(f"WARN {key}: {err}")

    # All-models Cldn4 ICB export: find every lung line that has Cldn4 + ICI.
    fetch_gene_icb(raw_dir, "Cldn4", "All", "cldn4_icb_all.csv")
    for gene in GENES:
        safe = gene.lower().replace("-", "")
        fetch_gene_icb(raw_dir, gene, "LLC", f"{safe}_llc.csv")
        try:
            fetch_gene_icb(raw_dir, gene, "CMT-167", f"{safe}_cmt167.csv")
        except Exception as err:
            print(f"WARN {gene} CMT-167: {err}")

    manifest = {
        "source": "https://tismo.pku-genomics.org/",
        "paper": "Zeng et al. Nucleic Acids Res 2022; PMID 34534350",
        "icb_treatments": ICB_TREATMENTS,
        "genes": GENES,
        "note": (
            "Public processed only. type=3 is TISMO gene-module per-sample CSV. "
            "tumorList=['All'] for Cldn4 keeps SRX8918393. No user-private 8 KL."
        ),
    }
    (raw_dir / "fetch_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("done")


if __name__ == "__main__":
    main()
