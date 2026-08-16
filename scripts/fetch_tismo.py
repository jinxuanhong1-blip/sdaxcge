#!/usr/bin/env python3
"""Download TISMO metadata and expression matrices.

TISMO (http://tismo.cistrome.org, now served from tismo.pku-genomics.org) hosts its
bulk files on an Aliyun share rather than plain HTTP. The share ids below are the ones
the Data Download page passes to its `getAiLiYunFile` helper.
"""
import json
import os
import sys
import urllib.request

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "tismo")

SHARES = {
    "cell_lines_meta": "YzQB2DQonQE",
    "sample_in_vitro_meta": "MXMDDtkQLqQ",
    "sample_in_vivo_meta": "voEA1DXBEFo",
    "expression_in_vitro": "AuDf6PnXFtV",
    "expression_in_vivo": "AQKzyRCi1Jb",
    "immune_cell": "kCmAULWGuJh",
}

TOKEN_URL = "https://bj21400.api.aliyunfile.com/v2/share_link/get_share_token"
LIST_URL = "https://bj21400.api.aliyunfile.com/v2/file/list"


def post_json(url, payload, headers=None):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=data, method="POST")
    req.add_header("Content-Type", "application/json")
    for k, v in (headers or {}).items():
        req.add_header(k, v)
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.loads(resp.read().decode())


def list_share(share_id):
    tok = post_json(TOKEN_URL, {"share_id": share_id, "ignoreError": True})
    token = tok.get("share_token")
    if not token:
        raise RuntimeError(f"no share token for {share_id}: {tok}")
    body = {
        "limit": 100,
        "marker": "",
        "share_id": share_id,
        "parent_file_id": "root",
        "fields": "user_name,dir_size,url,content_type,upload_id,crc64_hash,revision_id,description",
        "url_expire_sec": tok.get("expires_in", 7200),
    }
    return post_json(LIST_URL, body, {"x-share-token": token}).get("items", [])


def download(url, dest):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=1800) as resp, open(dest, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    wanted = sys.argv[1:] or list(SHARES)
    manifest = {}
    for name in wanted:
        share_id = SHARES[name]
        try:
            items = list_share(share_id)
        except Exception as exc:  # noqa: BLE001 - report and continue to next share
            print(f"[FAIL] {name}: {exc}")
            manifest[name] = {"share_id": share_id, "error": str(exc)}
            continue
        recs = []
        for item in items:
            fname = item.get("name")
            url = item.get("download_url") or item.get("url")
            print(f"[{name}] {fname} size={item.get('size')}")
            if not url:
                recs.append({"name": fname, "error": "no download url"})
                continue
            dest = os.path.join(OUT_DIR, fname)
            if os.path.exists(dest) and os.path.getsize(dest) == item.get("size"):
                print("   already present, skipping")
            else:
                download(url, dest)
            recs.append({"name": fname, "size": item.get("size"), "path": dest})
        manifest[name] = {"share_id": share_id, "files": recs}
    with open(os.path.join(OUT_DIR, "download_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps(manifest, indent=2)[:4000])


if __name__ == "__main__":
    main()
