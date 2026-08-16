#!/usr/bin/env python3
"""
Step 3: download every OPEN, <2GB processed series matrix found in step 2.

Files land in results/geo_sweep/matrices/ (git-ignored raw blobs). A manifest
of what was downloaded (path, bytes, sha-lite) is written to
notes/geo_sweep/download_manifest.json.
"""
import json
import time
import urllib.request
from pathlib import Path

NOTES = Path("notes/geo_sweep")
MATRIX_DIR = Path("results/geo_sweep/matrices")
MATRIX_DIR.mkdir(parents=True, exist_ok=True)
TWO_GB = 2 * 1024 ** 3


def download(url, dest, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-sweep/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
                f.write(r.read())
            return dest.stat().st_size
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"download failed {url}: {last}")


def main():
    probe = json.load(open(NOTES / "matrix_probe.json"))
    manifest = []
    for rec in probe:
        acc = rec["accession"]
        if not rec.get("open"):
            manifest.append({"accession": acc, "status": "not_open"})
            continue
        if rec["total_size"] >= TWO_GB:
            manifest.append({"accession": acc, "status": "skip_over_2gb",
                             "total_size": rec["total_size"]})
            continue
        files = []
        for f in rec["files"]:
            dest = MATRIX_DIR / f["name"]
            try:
                sz = download(f["url"], dest)
                files.append({"name": f["name"], "bytes": sz,
                              "path": str(dest)})
                print(f"{acc:12s} {f['name']:55s} {sz/1e6:8.3f}MB")
            except Exception as e:  # noqa: BLE001
                files.append({"name": f["name"], "error": str(e)})
                print(f"{acc:12s} {f['name']} ERROR {e}")
            time.sleep(0.15)
        manifest.append({"accession": acc, "status": "ok", "files": files})
    (NOTES / "download_manifest.json").write_text(json.dumps(manifest, indent=2))
    ok = sum(1 for m in manifest if m.get("status") == "ok")
    print(f"\nDownloaded matrices for {ok} series -> {MATRIX_DIR}")


if __name__ == "__main__":
    main()
