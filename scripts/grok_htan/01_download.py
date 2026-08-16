#!/usr/bin/env python3
"""Download a curated subset of eligible HTAN/HCA lung h5ad objects (<2GB each)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from lib_htan import ensure_dirs

# Curated for TACSTD2/CLDN4 vs immune in sc + spatial. URLs resolved live from CXG.
DOWNLOADS = [
    {"key": "hca_travaglini_10x", "dataset_id": "8c42cfd0-0b0a-46d5-910c-fc833d83c45e",
     "collection_id": "5d445965-6f1a-4b68-ba3a-b8f765155d3a", "kind": "scrna", "source": "HCA"},
    {"key": "hca_madissoon_lung", "dataset_id": "2672b679-8048-4f5e-9786-f1b196ccfd08",
     "collection_id": "4d74781b-8186-4c9a-b659-ff4dc4601d91", "kind": "scrna", "source": "HCA"},
    {"key": "htan_sclc_combined", "dataset_id": "576f193c-75d0-4a11-bd25-8676587e6dc2",
     "collection_id": "62e8f058-9c37-48bc-9200-e767f318a8ec", "kind": "scrna", "source": "HTAN"},
    {"key": "luad_histology_scrna", "dataset_id": "01ff5cf0-730f-4ddc-b1be-7b407211f544",
     "collection_id": "0bebef1a-4607-4584-9070-dacf89a0d635", "kind": "scrna", "source": "open_CXG_lung"},
    {"key": "htan_visium_Ctrl1_A1", "dataset_id": "b681d612-87ba-4234-8d37-5a7fd83526cb",
     "collection_id": "efd94500-1fdc-4e28-9e9f-a309d0154e21", "kind": "visium", "source": "HTAN"},
    {"key": "htan_visium_Ctrl1_B1", "dataset_id": "5e044266-bba1-4cc6-aa95-6b1a5b0811d7",
     "collection_id": "efd94500-1fdc-4e28-9e9f-a309d0154e21", "kind": "visium", "source": "HTAN"},
    {"key": "htan_visium_Ctrl2_C1", "dataset_id": "76df68b2-9608-45ce-9dd8-06b66905558f",
     "collection_id": "efd94500-1fdc-4e28-9e9f-a309d0154e21", "kind": "visium", "source": "HTAN"},
    {"key": "htan_visium_Ctrl2_D1", "dataset_id": "e48cdb8c-c436-418c-bfa9-7946044d408a",
     "collection_id": "efd94500-1fdc-4e28-9e9f-a309d0154e21", "kind": "visium", "source": "HTAN"},
    {"key": "htan_visium_DT1_A1", "dataset_id": "3dd69434-56fa-45da-8a23-579dcf3dc0e1",
     "collection_id": "efd94500-1fdc-4e28-9e9f-a309d0154e21", "kind": "visium", "source": "HTAN"},
    {"key": "htan_visium_DT1_B1", "dataset_id": "89b7ced8-7d2e-4c75-b0a8-45bb07d15c80",
     "collection_id": "efd94500-1fdc-4e28-9e9f-a309d0154e21", "kind": "visium", "source": "HTAN"},
    {"key": "htan_visium_DT2_C1", "dataset_id": "939794ef-3106-4126-a40c-541ee3c4b077",
     "collection_id": "efd94500-1fdc-4e28-9e9f-a309d0154e21", "kind": "visium", "source": "HTAN"},
    {"key": "htan_visium_DT2_D1", "dataset_id": "c9d294c7-857f-4ae5-8abe-6480d1584c6f",
     "collection_id": "efd94500-1fdc-4e28-9e9f-a309d0154e21", "kind": "visium", "source": "HTAN"},
    {"key": "hca_visium_LngSP10193347", "dataset_id": "708a7f62-260b-43f4-837a-d329c29f830b",
     "collection_id": "c1241244-b22d-483d-875b-75699efb9f3c", "kind": "visium", "source": "HCA"},
]


def curl_download(url: str, dest: Path, expected: int | None) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and expected and dest.stat().st_size >= 0.95 * expected:
        return {"key": dest.name, "status": "cached", "bytes": dest.stat().st_size, "path": str(dest)}
    tmp = dest.with_suffix(dest.suffix + ".part")
    cmd = [
        "curl", "-L", "--fail", "--retry", "5", "--retry-delay", "4",
        "--retry-all-errors", "-o", str(tmp), url,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not tmp.exists():
        err = (proc.stderr or proc.stdout or "")[-500:]
        if tmp.exists():
            tmp.unlink()
        return {"key": dest.name, "status": "fail", "error": err, "url": url}
    tmp.rename(dest)
    return {"key": dest.name, "status": "ok", "bytes": dest.stat().st_size, "path": str(dest)}


def resolve_asset(collection_id: str, dataset_id: str):
    import urllib.request

    req = urllib.request.Request(
        f"https://api.cellxgene.cziscience.com/curation/v1/collections/{collection_id}",
        headers={"Accept": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=90) as r:
        c = json.loads(r.read().decode())
    for d in c.get("datasets") or []:
        if d.get("dataset_id") == dataset_id:
            for a in d.get("assets") or []:
                if str(a.get("filetype", "")).upper() == "H5AD":
                    return a.get("url"), a.get("filesize")
    return None, None


def main():
    paths = ensure_dirs()
    manifest = []
    cache = {}
    for item in DOWNLOADS:
        cid = item["collection_id"]
        if cid not in cache:
            cache[cid] = True
        url, sz = resolve_asset(cid, item["dataset_id"])
        if not url:
            rec = {"status": "fail", "error": "no H5AD asset in collection API", **item}
            manifest.append(rec)
            print(f"FAIL {item['key']} no URL", flush=True)
            continue
        item["url"] = url
        item["bytes"] = sz
        dest = paths["data"] / f"{item['key']}.h5ad"
        print(f"GET {item['key']} {sz} <- {url}", flush=True)
        rec = curl_download(url, dest, sz)
        rec.update({k: item[k] for k in ("key", "dataset_id", "kind", "source", "url")})
        print(f"  {rec['status']} {rec.get('bytes')}", flush=True)
        manifest.append(rec)
    out = paths["notes"] / "download_manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    n_ok = sum(1 for m in manifest if m["status"] in ("ok", "cached"))
    print(f"done {n_ok}/{len(manifest)} -> {out}")
    if n_ok == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
