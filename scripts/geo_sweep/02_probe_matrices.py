#!/usr/bin/env python3
"""
Step 2: probe each candidate series' processed *series matrix* on the GEO FTP.

For every GSE from step 1 we list the /matrix/ directory over HTTPS, capture the
series_matrix.txt.gz file name(s) and byte size(s). This tells us:
  * whether the processed series matrix is publicly downloadable (OPEN), and
  * the total size, so we only download matrices < 2 GB (task constraint).

Writes notes/geo_sweep/matrix_probe.json.
"""
import json
import re
import time
import urllib.request
from pathlib import Path

NOTES = Path("notes/geo_sweep")
FTP_HTTPS = "https://ftp.ncbi.nlm.nih.gov/geo/series"
TWO_GB = 2 * 1024 ** 3


def ftp_dir_url(acc):
    # e.g. GSE50927 -> GSE50nnn ; GSE108417 -> GSE108nnn
    stub = acc[:-3] + "nnn" if len(acc) > 6 else "GSEnnn"
    return f"{FTP_HTTPS}/{stub}/{acc}/matrix/"


def http_get(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-sweep/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace"), r.status
        except Exception as e:  # noqa: BLE001
            last = e
            code = getattr(e, "code", None)
            if code == 404:
                return None, 404
            time.sleep(2 * (i + 1))
    return None, f"ERR:{last}"


# Apache autoindex rows look like:
# <a href="GSE50927_series_matrix.txt.gz">...</a>  DATE  SIZE
ROW = re.compile(
    r'href="(?P<name>[^"?][^"]*\.txt\.gz)".*?'
    r'(?P<date>\d{4}-\d{2}-\d{2} \d{2}:\d{2})\s+(?P<size>[\d.]+[KMG]?)',
)

UNIT = {"K": 1024, "M": 1024 ** 2, "G": 1024 ** 3, "": 1}


def parse_size(s):
    s = s.strip()
    m = re.match(r"([\d.]+)([KMG]?)", s)
    if not m:
        return None
    return int(float(m.group(1)) * UNIT[m.group(2)])


def probe(acc):
    url = ftp_dir_url(acc)
    html, status = http_get(url)
    rec = {"accession": acc, "matrix_dir": url, "dir_status": status,
           "files": [], "total_size": 0, "open": False}
    if html is None:
        return rec
    for m in ROW.finditer(html):
        name = m.group("name")
        size = parse_size(m.group("size"))
        rec["files"].append({"name": name, "size_bytes": size,
                             "url": url + name})
    rec["total_size"] = sum(f["size_bytes"] or 0 for f in rec["files"])
    rec["open"] = len(rec["files"]) > 0
    rec["under_2gb"] = rec["total_size"] < TWO_GB
    return rec


def main():
    d = json.load(open(NOTES / "search_raw.json"))
    accs = [r["Accession"] for r in d["records"]]
    out = []
    for acc in accs:
        rec = probe(acc)
        out.append(rec)
        sz = rec["total_size"]
        print(f"{acc:12s} open={rec['open']!s:5s} files={len(rec['files'])} "
              f"size={sz/1e6:8.2f}MB status={rec['dir_status']}")
        time.sleep(0.2)
    (NOTES / "matrix_probe.json").write_text(json.dumps(out, indent=2))
    print(f"\nWrote {NOTES/'matrix_probe.json'}")
    n_open = sum(1 for r in out if r["open"])
    n_over = sum(1 for r in out if r["open"] and not r["under_2gb"])
    print(f"open matrices: {n_open}/{len(out)}   over-2GB: {n_over}")


if __name__ == "__main__":
    main()
