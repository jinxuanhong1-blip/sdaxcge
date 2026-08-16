"""Download OPEN PROCESSED supplementary/matrix files (< 2 GB) for leftover
human-lung immunotherapy series (and, separately, the mouse Pdl1 lung-SCC
superseries if it has a small processed matrix).

Files are written under results/w200/GEO_2010_2014/data/<ACC>/. Files >= 2 GB
are recorded as skipped. Nothing is fabricated.
"""
import json
import os
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "w200", "GEO_2010_2014")
DATA = os.path.join(OUT, "data")
os.makedirs(DATA, exist_ok=True)

TWO_GB = 2 * 1024**3


def ftp_prefix(acc):
    return acc[:-3] + "nnn"


def load_scope():
    recs = json.load(open(os.path.join(OUT, "verified_series.json")))
    # Download leftover human lung immunotherapy (even if not checkpoint ICI)
    # plus mouse lung Pdl1-expression series (biology leftover, not treatment).
    scope = set()
    for r in recs:
        if r["is_human"] and r["is_lung"] and r["is_immunotherapy"]:
            scope.add(r["accession"])
        if (not r["is_human"]) and r["is_lung"] and r.get("is_pdl1_gene_mention"):
            scope.add(r["accession"])
    return scope


def load_files():
    rows = []
    with open(os.path.join(OUT, "suppl_files.tsv")) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, parts)))
    return rows


def download(url, dest, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-2010-2014-leftover/1.0"})
            with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as out:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            return True
        except Exception as e:  # noqa: BLE001
            print(f"  retry {i+1}: {e}", flush=True)
            time.sleep(2 ** i)
    return False


def main():
    scope = load_scope()
    files = load_files()
    manifest = []
    for r in files:
        acc = r["accession"]
        if acc not in scope:
            continue
        size_b = int(r["size_bytes"])
        name = r["filename"]
        sub = r["location"]
        url = (
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/{ftp_prefix(acc)}/"
            f"{acc}/{sub}/{name}"
        )
        dest_dir = os.path.join(DATA, acc)
        os.makedirs(dest_dir, exist_ok=True)
        dest = os.path.join(dest_dir, name)
        entry = {
            "accession": acc,
            "location": sub,
            "filename": name,
            "listed_size": r["size"],
            "listed_bytes": size_b,
            "url": url,
        }
        if size_b >= TWO_GB:
            entry["status"] = "skipped_ge_2gb"
            entry["stored_bytes"] = 0
            print(f"SKIP >=2GB {acc}/{name} ({r['size']})", flush=True)
        elif os.path.exists(dest) and os.path.getsize(dest) == size_b:
            entry["status"] = "already_present"
            entry["stored_bytes"] = os.path.getsize(dest)
            print(f"HAVE {acc}/{name}", flush=True)
        else:
            print(f"GET  {acc}/{name} ({r['size']}) ...", flush=True)
            ok = download(url, dest)
            if ok:
                entry["status"] = "downloaded"
                entry["stored_bytes"] = os.path.getsize(dest)
            else:
                entry["status"] = "failed"
                entry["stored_bytes"] = os.path.getsize(dest) if os.path.exists(dest) else 0
        manifest.append(entry)

    cols = [
        "accession",
        "location",
        "filename",
        "listed_size",
        "listed_bytes",
        "stored_bytes",
        "status",
        "url",
    ]
    with open(os.path.join(OUT, "download_manifest.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for e in manifest:
            f.write("\t".join(str(e.get(c, "")) for c in cols) + "\n")
    with open(os.path.join(OUT, "download_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    n_ok = sum(1 for e in manifest if e["status"] in ("downloaded", "already_present"))
    n_skip = sum(1 for e in manifest if e["status"] == "skipped_ge_2gb")
    n_fail = sum(1 for e in manifest if e["status"] == "failed")
    print(f"download: ok={n_ok} skipped_ge_2gb={n_skip} failed={n_fail}")


if __name__ == "__main__":
    main()
