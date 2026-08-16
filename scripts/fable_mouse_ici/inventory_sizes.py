#!/usr/bin/env python3
"""Inventory supplementary-file sizes for each accession.

GEO: parse the HTTPS suppl directory listing (Apache autoindex) which reports a
size for every file, including the contents behind ``*_RAW.tar`` bundles are not
expanded but the tar size itself is reported. ArrayExpress/BioStudies: use the
file-list JSON endpoint.

Writes notes/fable_mouse_ici/file_inventory.json and prints a summary flagging
any file > 2 GB (to be skipped per task rules).
"""
import json
import re
import time
import urllib.request
from pathlib import Path

GEO_ACCS = [
    "GSE239485", "GSE133604", "GSE222158", "GSE129297", "GSE241978",
    "GSE330658", "GSE197260", "GSE297630", "GSE297632",
]
AE_ACCS = ["E-MTAB-13704"]
TWO_GB = 2 * 1024 ** 3

NOTES = Path(__file__).resolve().parents[2] / "notes" / "fable_mouse_ici"
UA = {"User-Agent": "Mozilla/5.0 (fable-mouse-ici inventory)"}


def fetch(url, tries=4, binary=False):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=90) as r:
                data = r.read()
                return data if binary else data.decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** i)
    raise last


def head_size(url, tries=3):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA, method="HEAD")
            with urllib.request.urlopen(req, timeout=60) as r:
                cl = r.headers.get("Content-Length")
                return int(cl) if cl is not None else None
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** i)
    return None


def geo_suppl_dir(acc):
    nnn = acc[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{acc}/suppl/"


def inventory_geo(acc):
    base = geo_suppl_dir(acc)
    rec = {"accession": acc, "source": "GEO", "suppl_url": base, "files": []}
    try:
        html = fetch(base)
    except Exception as e:  # noqa: BLE001
        rec["error"] = str(e)
        return rec
    # Apache autoindex rows: <a href="name">name</a> ... date  size
    for m in re.finditer(r'<a href="([^"?/][^"]*)">', html):
        name = m.group(1)
        if name in ("Parent Directory",) or name.endswith("/"):
            continue
        url = base + name
        size = head_size(url)
        rec["files"].append({"name": name, "url": url, "size_bytes": size})
        time.sleep(0.2)
    return rec


def inventory_ae(acc):
    rec = {"accession": acc, "source": "ArrayExpress/BioStudies", "files": []}
    url = f"https://www.ebi.ac.uk/biostudies/files/{acc}/dir?"  # fallback below
    # Reliable: the study JSON lists a filelist; use the files API.
    api = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{acc}"
    try:
        study = json.loads(fetch(api))
    except Exception as e:  # noqa: BLE001
        rec["error"] = str(e)
        return rec

    def walk(sec):
        files = sec.get("files", [])
        for f in files:
            entries = f if isinstance(f, list) else [f]
            for fe in entries:
                path = fe.get("path")
                size = fe.get("size")
                rec["files"].append({
                    "name": path,
                    "url": f"https://www.ebi.ac.uk/biostudies/files/{acc}/{path}",
                    "size_bytes": int(size) if size is not None else None,
                })
        for sub in sec.get("subsections", []) or []:
            if isinstance(sub, dict):
                walk(sub)
            elif isinstance(sub, list):
                for s in sub:
                    if isinstance(s, dict):
                        walk(s)

    walk(study.get("section", {}))
    return rec


def main():
    results = []
    for acc in GEO_ACCS:
        print("[inv]", acc)
        results.append(inventory_geo(acc))
    for acc in AE_ACCS:
        print("[inv]", acc)
        results.append(inventory_ae(acc))

    (NOTES / "file_inventory.json").write_text(json.dumps(results, indent=2))

    print("\n==== FILE INVENTORY (flag >2GB) ====")
    for r in results:
        print("=" * 70)
        print(r["accession"], r.get("error", ""))
        for f in r["files"]:
            sz = f["size_bytes"]
            human = f"{sz/1024**2:.1f} MB" if sz else "?"
            flag = "  <<< SKIP >2GB" if (sz and sz > TWO_GB) else ""
            print(f"   {human:>12}  {f['name']}{flag}")


if __name__ == "__main__":
    main()
