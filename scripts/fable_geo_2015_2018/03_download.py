"""Download all OPEN PROCESSED supplementary/matrix files (< 2 GB each) for the
in-scope verified series (human AND lung AND immune-checkpoint relevant).

Files are written under results/fable_geo_2015_2018/data/<ACC>/. Files >= 2 GB
are recorded in the manifest as skipped (per the "< 2 GB" rule); nothing is
fabricated. A download_manifest.tsv records the exact byte size actually stored.
"""
import json
import os
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "fable_geo_2015_2018")
DATA = os.path.join(OUT, "data")
os.makedirs(DATA, exist_ok=True)

TWO_GB = 2 * 1024**3


def ftp_prefix(acc):
    return acc[:-3] + "nnn"


def load_scope():
    recs = json.load(open(os.path.join(OUT, "verified_series.json")))
    scope = {r["accession"] for r in recs
             if r["is_human"] and r["is_lung"] and r["is_ici"]}
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
            req = urllib.request.Request(url, headers={"User-Agent": "fable-geo/1.0"})
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
        url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{ftp_prefix(acc)}/"
               f"{acc}/{sub}/{name}")
        dest_dir = os.path.join(DATA, acc)
        os.makedirs(dest_dir, exist_ok=True)
        dest = os.path.join(dest_dir, name)
        entry = {"accession": acc, "location": sub, "filename": name,
                 "listed_size": r["size"], "listed_bytes": size_b, "url": url}
        if size_b >= TWO_GB:
            entry["status"] = "skipped_ge_2gb"
            entry["stored_bytes"] = 0
            print(f"SKIP (>=2GB) {acc}/{name} ({r['size']})", flush=True)
            manifest.append(entry)
            continue
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            entry["status"] = "already_present"
            entry["stored_bytes"] = os.path.getsize(dest)
            manifest.append(entry)
            continue
        print(f"GET {acc}/{sub}/{name} ({r['size']})", flush=True)
        ok = download(url, dest)
        entry["status"] = "downloaded" if ok else "failed"
        entry["stored_bytes"] = os.path.getsize(dest) if ok and os.path.exists(dest) else 0
        manifest.append(entry)

    with open(os.path.join(OUT, "download_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    cols = ["accession", "location", "filename", "listed_size", "listed_bytes",
            "status", "stored_bytes", "url"]
    with open(os.path.join(OUT, "download_manifest.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for e in manifest:
            f.write("\t".join(str(e.get(c, "")) for c in cols) + "\n")
    dl = [e for e in manifest if e["status"] in ("downloaded", "already_present")]
    total = sum(e["stored_bytes"] for e in dl)
    print(f"\nstored {len(dl)} files, {total/1024**2:.1f} MB total; "
          f"skipped>=2GB={sum(1 for e in manifest if e['status']=='skipped_ge_2gb')}; "
          f"failed={sum(1 for e in manifest if e['status']=='failed')}")


if __name__ == "__main__":
    main()
